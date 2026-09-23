# Ingest the kagome half of `kagome.superconductors` -- Hubbard on kagome
# cylinders, run by the METTS.jl driver
# `kagome.superconductors/code/metts_hubbard.jl`.
#
#   run path:  <data>/<latdir>/t.<t>.U.<U>.V.<V>/nparticles.<n>/T.<T>/
#              tau.<tau>.cutoff.<c>.maxdim.<m>/outfile.seed.<s>.samples.txt
#   states:    `<step>: [1, 2, 3, ...]`, 1-based *Electron* indices
#              (1 Emp, 2 Up, 3 Dn, 4 UpDn)
#   basis:     always X -- the driver calls collapse_with_qn!(psi, "X")
#              unconditionally, at line 104 and again at 196
#   energy:    companion outfile.seed.<s>.h5, dataset `energy`
#
# The lattice needs converting. The project keeps its own TOML schema
# (`coordinates`, `[hoppings]`, `[exchanges]`) with **1-based** site indices,
# against the library's `Coordinates` / `Interactions` and 0-based indices.
#
# Only `[hoppings]` survives the conversion. The driver reads
# `lattice["hoppings"]["sites"]` (metts_hubbard.jl:52) and nothing else, then
# builds fermi_hubbard_opsum(N, t, U, V, neighbors): -t hopping on those bonds,
# U on-site Nupdn, and V*Ntot*Ntot on *those same bonds*. The `[exchanges]`
# section, coupling "J", is never read -- and in all four lattice files its
# site list is byte-identical to the hoppings list, so dropping it loses no
# geometry, only a coupling no run ever gave a value to. Keeping it would make
# the lattice declare a J that `validate` then demands a parameter for.
#
# U and V are therefore not lattice couplings either; like U in
# Hubbard/superconductors they ride in `parameters`, which validate allows to
# name more than the lattice does. V acts on the lattice's own bond list.

using METTSLibrary, HDF5, TOML
const ML = METTSLibrary

const ROOT  = "/data/condmat/awietek/Research/Projects/kagome.superconductors/hubbard/data"
const LATD  = "/home/awietek/Research/Projects/kagome.superconductors/lattice-files"
const DRIVER = "kagome.superconductors/code/metts_hubbard.jl"

const RE_LAT = r"^kagome\.L\.(\d+)\.W\.(\d+)\.cylinder$"
const RE_PAR = r"^t\.(-?\d+\.\d+)\.U\.(-?\d+\.\d+)\.V\.(-?\d+\.\d+)$"
const RE_N   = r"^nparticles\.(\d+)$"
const RE_T   = r"^T\.(\d+\.\d+)$"
const RE_ALG = r"^tau\.([0-9.eE+-]+)\.cutoff\.([0-9.eE+-]+)\.maxdim\.(\d+)$"
const RE_SEED = r"^outfile\.seed\.(\d+)\.samples\.txt$"

"""
Convert one of the project's lattice TOMLs to the library's schema.

The project writes `coordinates` plus `[hoppings]`/`[exchanges]` tables whose
`sites` are 1-based pairs; the library wants `Coordinates` and
`Interactions = [[coupling, type, i, j], ...]` with 0-based sites. Only the
hoppings are carried over -- see the note at the top of this file.
"""
function kagome_lattice_toml(src::AbstractString)
    d = TOML.parsefile(src)
    coords = d["coordinates"]
    hop    = d["hoppings"]
    bonds  = hop["sites"]
    cname, ctype = hop["coupling"], hop["type"]

    n = length(coords)
    flat = [s for b in bonds for s in b]
    # The 1-based -> 0-based shift is the whole point of this function, so
    # check the source really is 1-based rather than trusting the schema.
    (minimum(flat) == 1 && maximum(flat) == n) ||
        error("'$src': sites span $(minimum(flat))..$(maximum(flat)), expected 1..$n")
    all(b -> length(b) == 2, bonds) || error("'$src': not every bond is a pair")

    io = IOBuffer()
    println(io, "Coordinates = [")
    for c in coords
        println(io, "  [", join(c, ", "), "],")
    end
    println(io, "]")
    println(io)
    println(io, "Interactions = [")
    for b in bonds
        println(io, "  ['", cname, "', '", ctype, "', ", b[1] - 1, ", ", b[2] - 1, "],")
    end
    println(io, "]")
    return String(take!(io))
end

function parse_run(p::AbstractString)
    parts = splitpath(p)
    length(parts) < 6 && return nothing
    mF = match(RE_SEED, parts[end]);   mF === nothing && return nothing
    mA = match(RE_ALG,  parts[end-1]); mA === nothing && return nothing
    mT = match(RE_T,    parts[end-2]); mT === nothing && return nothing
    mN = match(RE_N,    parts[end-3]); mN === nothing && return nothing
    mP = match(RE_PAR,  parts[end-4]); mP === nothing && return nothing
    mL = match(RE_LAT,  parts[end-5]); mL === nothing && return nothing

    L, W = parse(Int, mL[1]), parse(Int, mL[2])
    lat = joinpath(LATD, parts[end-5] * ".toml")
    isfile(lat) || error("no lattice file '$lat' for '$p'")
    return (; path = p, lattice = lat,
            lattice_name = "kagome.L$(L).W$(W).cyl",
            parameters = Dict("T" => parse(Float64, mP[1]), "U" => parse(Float64, mP[2]),
                              "V" => parse(Float64, mP[3])),
            sector      = Dict("n" => parse(Int, mN[1])),
            temperature = parse(Float64, mT[1]),
            algorithm = Dict{String,Any}(
                "tau"    => parse(Float64, mA[1]), "cutoff" => parse(Float64, mA[2]),
                "maxdim" => parse(Int, mA[3]),     "seed"   => parse(Int, mF[1]),
                "driver" => DRIVER, "time_evolution" => "TDVP tau 0.1 (applyexp, kkrylov 2, tau0 0.01)"))
end

"""
Per-step energies for a run, aligned to its `nsamp` states, or nothing.

Each loop iteration measures the energy, `dump!`s it, *then* collapses and
saves the sample (metts_hubbard.jl:131-202), so energy k belongs to sample k.
Two things break that correspondence:

  * a run killed between the energy dump and the sample save leaves one extra
    energy at the END of the array -- harmless, trim it;
  * `dump!` opens the file "cw" and appends to an extensible dataset, so a run
    that died there and was then RESUMED re-dumps that step's energy, leaving
    an extra INSIDE the array that shifts every entry after it.

The two cannot be told apart by value: on resume the driver re-runs sloppy
DMRG from the last saved sample (line 93-111, unconditional) rather than
continuing from the same state, so the re-measured energy differs from the one
already on disk -- an interior duplicate is not a repeated number.

What does separate them is which file was written last. The sample save is the
final write of a completed step, so a run that ran on after resuming leaves
samples.txt newest; a run whose last act was the energy dump leaves the .h5
newest. Measured over all 3,161 runs, that holds without a single exception:
of the 2,297 whose lengths already agree, **zero** have the .h5 written last.
So one extra energy is trimmed only when the .h5 is the newer file.

Everything else is dropped rather than misaligned. Over the archive that is
2,297 exact + 341 trimmed = 2,638 runs with energy, and 523 (16.5%) without:
251 off by one the wrong way round, 71 off by two, and 201 off by anywhere up
to 418 -- runs resumed so often that the energy array bears no relation to the
surviving samples.

Returns `(energies_or_nothing, status)`; the status goes into provenance so a
file without energy says why.
"""
function energies(h5path, samples_path, nsamp)
    isfile(h5path) || return nothing, "no_companion_h5"
    e = try
        h5open(h5path, "r") do h
            haskey(h, "energy") || return nothing
            return Float64.(vec(read(h["energy"])))
        end
    catch
        return nothing, "unreadable_h5"
    end
    e === nothing && return nothing, "no_energy_dataset"
    d = length(e) - nsamp
    d == 0 && return e, "exact"
    d == 1 && mtime(h5path) > mtime(samples_path) && return e[1:nsamp], "trimmed_one_trailing"
    return nothing, "dropped_unalignable_by_$d"
end

function convert_one(r, project)
    e = from_samples_txt(r.path;
        lattice = kagome_lattice_toml(r.lattice), lattice_name = r.lattice_name,
        model = "Hubbard", project = project, site_type = "Electron",
        temperature = r.temperature, basis = "X",
        parameters = r.parameters, sector = r.sector, algorithm = r.algorithm)
    # Align against the converted ensemble's own sample count, not a second
    # pass over the file: from_samples_txt drops empty bodies, so a line count
    # is not guaranteed to agree with it.
    en, status = energies(replace(r.path, ".samples.txt" => ".h5"), r.path, nsamples(e))
    obs = en === nothing ? Dict{String,Array{Float64}}() :
                           Dict{String,Array{Float64}}("energy" => en)
    return Ensemble(model = e.model, project = e.project, site_type = e.site_type,
                    local_states = e.local_states, lattice = e.lattice,
                    lattice_name = e.lattice_name, beta = e.beta, parameters = e.parameters,
                    sector = e.sector, algorithm = e.algorithm,
                    provenance = merge(e.provenance,
                                       Dict{String,Any}("energy_alignment" => status)),
                    collapse_bases = e.collapse_bases, states = e.states, basis = e.basis,
                    observables = obs)
end

function main()
    out     = get(ENV, "OUT", "")
    project = get(ENV, "PROJECT", "kagome.superconductors")
    maxf    = parse(Int, get(ENV, "MAXFILES", "0"))
    dry     = get(ENV, "DRYRUN", "0") == "1"
    resume  = get(ENV, "RESUME", "0") == "1"

    all_runs = String[]
    for (d, _, fs) in walkdir(ROOT)
        occursin(Base.Filesystem.path_separator * "prune", d) && continue
        for f in fs
            endswith(f, ".samples.txt") && push!(all_runs, joinpath(d, f))
        end
    end
    sort!(all_runs)
    println("samples.txt found: ", length(all_runs))

    rs = Any[]; nunparsed = 0
    for p in all_runs
        r = parse_run(p)
        r === nothing ? (nunparsed += 1; nunparsed <= 5 && println("  UNPARSED: ", p)) : push!(rs, r)
    end
    println("parsed: ", length(rs), "   unparsed: ", nunparsed)

    # One tree, one file per run: no backup subtree to deduplicate here. The
    # tag must still be unique inside each ensemble.
    byens = Dict{String,Vector{String}}()
    for r in rs
        k = string(r.lattice_name, "|", sort(collect(r.parameters)), "|",
                   sort(collect(r.sector)), "|", r.temperature)
        tag = string("basis=X_maxdim=", r.algorithm["maxdim"], "_tau=", r.algorithm["tau"],
                     "_cutoff=", r.algorithm["cutoff"], "_seed=", r.algorithm["seed"])
        push!(get!(byens, k, String[]), tag)
    end
    dups = [k for (k, v) in byens if length(unique(v)) < length(v)]
    println("ensembles: ", length(byens), "   with a duplicate tag: ", length(dups))
    for k in first(dups, 3); println("  DUP ", k); end
    isempty(dups) || error("duplicate tags; the tag field set is missing something")

    sort!(rs; by = r -> r.path)
    nsh = parse(Int, get(ENV, "NSHARDS", "1")); sh = parse(Int, get(ENV, "SHARD", "0"))
    0 <= sh < nsh || error("bad shard $sh of $nsh")
    nsh > 1 && (rs = rs[(sh+1):nsh:end]; println("shard $sh/$nsh: ", length(rs), " runs"))
    sel = maxf > 0 ? first(rs, maxf) : rs
    println("converting: ", length(sel), dry ? "  (DRY RUN, writing nothing)" : "  -> $out")

    nok = 0; nerr = 0; nskip = 0; tot = 0; noen = 0
    for (i, r) in enumerate(sel)
        try
            e = convert_one(r, project)
            rel = ML.relpath_for(e)
            tot += nsamples(e)
            haskey(e.observables, "energy") || (noen += 1)
            if resume && !dry && isfile(joinpath(out, rel))
                nskip += 1; continue
            end
            if dry
                # Electron: 0 = Emp, 3 = UpDn counts twice
                nelec = sum(s -> s == 0 ? 0 : (s == 3 ? 2 : 1), e.states) / nsamples(e)
                println("  ", rel)
                println("      ", nsamples(e), " samples, ", nsites(e), " sites, basis ",
                        e.collapse_bases, ", <n> per sample ", round(nelec, digits=2),
                        " (sector n=", e.sector["n"], "), energy ",
                        e.provenance["energy_alignment"])
            else
                write_ensemble(out, e)
            end
            nok += 1
        catch err
            nerr += 1
            println("  ERROR ", r.path, "\n        ", sprint(showerror, err))
        end
        i % 200 == 0 && !dry && println("  ... $i/$(length(sel))")
    end
    println("done: ok=", nok, " skipped=", nskip, " err=", nerr,
            " samples=", tot, " without-energy=", noen)
    nerr == 0 || exit(1)
end

main()
