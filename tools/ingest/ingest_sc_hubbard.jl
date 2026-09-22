# Ingest the Hubbard half of the `superconductors` project, into
# Hubbard/superconductors/ -- the same project name under a different model
# directory, which is exactly what the <model>/<project>/ layout is for.
#
#   dump path:  <tree>/<latdir>/t.<t>.tp.<tp>.U.<U>.holes.<h>/T.<T>/outfile.….dump.h5
#   lattice:    <latdir>.hubbard.lat, declaring T and Tp on HUBBARDHOP bonds
#   site type:  Electron -- states take values 0..3 (Emp, Up, Dn, UpDn),
#               against tJ's 0..2
#   sector:     n = nx*ny - holes
#
# U is NOT a lattice coupling. The .lat declares only T and Tp; U is the onsite
# repulsion, which the C++ code reads from the generated coupling file. It is
# carried in `parameters` anyway -- validate allows parameters the lattice does
# not name, and U is the knob these runs vary, so it belongs in the path.
#
# Time evolution: metts_hubbard takes no --tau. It reads a TEBD -> TDVP
# schedule from a generated timeevofile, identical in every surviving Hubbard
# script: TEBD 0->0.1 at tau 0.02, then TDVP 0.1->beta/2 at tau 0.5. The bulk
# step is therefore 0.5 with 0.02 as the initial stage -- the analogue of the
# t-J tau/init_tau pair, and not the 0.2/0.1 those runs used.

using METTSLibrary, HDF5
const ML = METTSLibrary

const TREES = ["/data/condmat/awietek/Research/Projects/superconductors/hubbard/metts/outfiles",
               "/data/condmat/awietek/Research/Projects/superconductors/hubbard/metts/outfiles.metts",
               "/data/condmat/awietek/Research/Projects/superconductors/hubbard_bkup/hubbard/metts/outfiles.metts"]
const LATDIR = "/home/awietek/Research/Projects/superconductors/hubbard/lattice-files"
const RUN_SCRIPT = "superconductors/hubbard/metts/run_metts_hubbard_auto.sh " *
                   "(timeevofile: TEBD 0->0.1 tau 0.02 cutoff 1e-12; TDVP 0.1->beta/2 tau 0.5)"

# U is written without decimals ("U.10"), t and tp with them ("t.1.00").
const RE_PAR  = r"^t\.(-?\d+\.\d+)\.tp\.(-?\d+\.\d+)\.U\.(-?\d+(?:\.\d+)?)\.holes\.(\d+)$"
const RE_T    = r"^T\.(\d+\.\d+)$"
const RE_LAT  = r"^square\.op\.nx\.(\d+)\.ny\.(\d+)$"
const RE_FILE = r"maxm\.(\d+)\.cutoff\.([0-9.eE+-]+?)\.updates\.([xz])\.nmetts\.(\d+)\.nwarm\.(\d+)\.seed\.(\d+)"

function parse_run(p::AbstractString)
    parts = splitpath(p)
    length(parts) < 4 && return nothing
    mT = match(RE_T,    parts[end-1]); mT === nothing && return nothing
    mP = match(RE_PAR,  parts[end-2]); mP === nothing && return nothing
    mL = match(RE_LAT,  parts[end-3]); mL === nothing && return nothing
    mF = match(RE_FILE, parts[end]);   mF === nothing && return nothing

    nx, ny = parse(Int, mL[1]), parse(Int, mL[2])
    lat = joinpath(LATDIR, parts[end-3] * ".hubbard.lat")
    isfile(lat) || error("no lattice file '$lat' for '$p'")
    return (; path = p, lattice = lat,
            parameters = Dict("T"  => parse(Float64, mP[1]), "Tp" => parse(Float64, mP[2]),
                              "U"  => parse(Float64, mP[3])),
            sector      = Dict("n" => nx * ny - parse(Int, mP[4])),
            temperature = parse(Float64, mT[1]),
            basis       = uppercase(mF[3]),
            algorithm = Dict{String,Any}(
                "maxdim" => parse(Int, mF[1]), "cutoff" => parse(Float64, mF[2]),
                "nmetts" => parse(Int, mF[4]), "nwarm"  => parse(Int, mF[5]),
                "seed"   => parse(Int, mF[6]),
                "tau" => 0.5, "init_tau" => 0.02,
                "time_evolution" => "TEBD 0->0.1 tau 0.02 cutoff 1e-12; TDVP 0.1->beta/2 tau 0.5"),
            # The dedup key must be exactly what determines the output path:
            # parsed values, never the raw directory or filename strings. The
            # two sibling trees write the same run with cosmetically different
            # spellings -- notably cutoff.1e-6 against cutoff.1.0e-6, identical
            # as Float64 and therefore identical in the tag -- so a
            # string-based key leaves them as "distinct" runs that then collide
            # on the same path. Keying on parsed values makes a collision
            # impossible by construction: same path <=> same key <=> deduped.
            key = (parts[end-3],
                   parse(Float64, mP[1]), parse(Float64, mP[2]), parse(Float64, mP[3]),
                   parse(Int, mP[4]), parse(Float64, mT[1]),
                   parse(Int, mF[1]), uppercase(mF[3]), parse(Int, mF[6]),
                   parse(Float64, mF[2])))
end

"Number of METTS steps in a dump, 0 if it holds none. Header read only."
function nsteps(p)
    try
        return h5open(p, "r") do h
            haskey(h, "ProductState") || return 0
            return size(dataspace(h["ProductState"]))[2]
        end
    catch
        return 0
    end
end

function convert_one(r, project)
    e = from_legacy_cpp(r.path; lattice = r.lattice, model = "Hubbard", project = project,
                        site_type = "Electron", temperature = r.temperature, basis = r.basis,
                        parameters = r.parameters, sector = r.sector, algorithm = r.algorithm)
    return Ensemble(model = e.model, project = e.project, site_type = e.site_type,
                    local_states = e.local_states, lattice = e.lattice,
                    lattice_name = e.lattice_name, beta = e.beta, parameters = e.parameters,
                    sector = e.sector, algorithm = e.algorithm,
                    provenance = merge(e.provenance,
                                       Dict{String,Any}("time_evolution_source" => RUN_SCRIPT)),
                    collapse_bases = e.collapse_bases, states = e.states, basis = e.basis,
                    step = e.step, observables = e.observables)
end

function main()
    out     = get(ENV, "OUT", "")
    project = get(ENV, "PROJECT", "superconductors")
    maxf    = parse(Int, get(ENV, "MAXFILES", "0"))
    dry     = get(ENV, "DRYRUN", "0") == "1"
    resume  = get(ENV, "RESUME", "0") == "1"

    all_dumps = String[]
    for t in TREES
        isdir(t) || (println("  (no such tree, skipping: ", t, ")"); continue)
        for (d, _, fs) in walkdir(t), f in fs
            endswith(f, ".h5") && push!(all_dumps, joinpath(d, f))
        end
    end
    sort!(all_dumps)
    println("dumps found: ", length(all_dumps))

    rs = Any[]; nunparsed = 0
    for p in all_dumps
        r = parse_run(p)
        r === nothing ? (nunparsed += 1; nunparsed <= 5 && println("  UNPARSED: ", p)) : push!(rs, r)
    end
    println("parsed: ", length(rs), "   unparsed: ", nunparsed)

    # Same deduplication as the t-J side: hubbard/ and hubbard_bkup/ overlap,
    # so keep, per distinct run, the copy with the most steps.
    best = Dict{Any,Any}(); steps = Dict{Any,Int}()
    for r in rs
        n = nsteps(r.path)
        n == 0 && continue
        if !haskey(best, r.key) || n > steps[r.key]
            best[r.key] = r; steps[r.key] = n
        end
    end
    kept = collect(values(best))
    println("distinct runs with states: ", length(kept),
            "   (dropped ", length(rs) - length(kept), " duplicate or stateless copies)")
    nb = count(r -> occursin("_bkup", r.path), kept)
    println("   from hubbard/: ", length(kept) - nb, "   from hubbard_bkup/: ", nb)

    byens = Dict{String,Vector{String}}()
    for r in kept
        k = string(r.lattice, "|", sort(collect(r.parameters)), "|",
                   sort(collect(r.sector)), "|", r.temperature)
        tag = string("basis=", r.basis, "_maxdim=", r.algorithm["maxdim"],
                     "_tau=", r.algorithm["tau"], "_cutoff=", r.algorithm["cutoff"],
                     "_seed=", r.algorithm["seed"])
        push!(get!(byens, k, String[]), tag)
    end
    dups = [k for (k, v) in byens if length(unique(v)) < length(v)]
    println("ensembles: ", length(byens), "   with a duplicate tag: ", length(dups))
    for k in first(dups, 3); println("  DUP ", k); end
    isempty(dups) || error("duplicate tags; the tag field set is missing something")

    sort!(kept; by = r -> r.path)
    nsh = parse(Int, get(ENV, "NSHARDS", "1")); sh = parse(Int, get(ENV, "SHARD", "0"))
    nsh > 1 && (kept = kept[(sh+1):nsh:end]; println("shard $sh/$nsh: ", length(kept), " runs"))
    sel = maxf > 0 ? first(kept, maxf) : kept
    println("converting: ", length(sel), dry ? "  (DRY RUN, writing nothing)" : "  -> $out")

    nok = 0; nerr = 0; nskip = 0; tot = 0
    for (i, r) in enumerate(sel)
        try
            e = convert_one(r, project)
            rel = ML.relpath_for(e)
            tot += nsamples(e)
            if resume && !dry && isfile(joinpath(out, rel))
                nskip += 1; continue
            end
            if dry
                # Electron: 0 = Emp, so occupancy counts singly AND doubly occupied
                nelec = sum(s -> s == 0 ? 0 : (s == 3 ? 2 : 1), e.states) / nsamples(e)
                println("  ", rel)
                println("      ", nsamples(e), " samples, ", nsites(e), " sites, basis ",
                        e.collapse_bases, ", <n> per sample ", round(nelec, digits=2),
                        " (sector n=", e.sector["n"], ")  <- ",
                        occursin("_bkup", r.path) ? "hubbard_bkup" : "hubbard")
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
    println("done: ok=", nok, " skipped=", nskip, " err=", nerr, " samples=", tot)
    nerr == 0 || exit(1)
end

main()
