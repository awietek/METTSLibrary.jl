# Ingest the t-J half of the `superconductors` project.
#
#   dump path:  <tree>/<latdir>/t.<t>.tp.<tp>.J.<J>.holes.<h>/T.<T>/outfile.….dump.h5
#   lattice:    <latdir>.ttpJ.lat, declaring T Tp J
#   sector:     n = nx*ny - holes
#   basis:      updates.x|z in the filename
#   algorithm:  maxm/cutoff/nmetts/nwarm/seed from the filename, plus
#               tau=0.2 init_tau=0.1 k=3 from run_metts_tj_auto.sh, which the
#               dump does not carry (see provenance.time_evolution_source).
#
# Deduplication. The runs live in two trees, tj/ and tj_bkup/, which overlap:
# of 3787 distinct runs, 2000 are only in tj/, 592 only in tj_bkup/ (all of the
# 32x6 12-hole family), and 1195 in both. Where a run is in both the copies are
# usually identical, but 198 are longer in tj/ and 2 are longer in tj_bkup/ --
# so neither tree wins outright and picking by tree would throw away samples.
# Each run is therefore taken from whichever copy has more METTS steps, read
# from the ProductState dataspace (a header read, no data transfer).
#
# Env: OUT, PROJECT, MAXFILES (0 = all), SHARD/NSHARDS, DRYRUN, RESUME, LATMATCH.

using METTSLibrary, HDF5
const ML = METTSLibrary

const TREES = ["/data/condmat/awietek/Research/Projects/superconductors/tj/metts/outfiles",
               "/data/condmat/awietek/Research/Projects/superconductors/tj/metts/outfiles.metts",
               "/data/condmat/awietek/Research/Projects/superconductors/tj_bkup/metts/outfiles",
               "/data/condmat/awietek/Research/Projects/superconductors/tj_bkup/metts/outfiles.metts"]
const LATDIR = "/home/awietek/Research/Projects/superconductors/tj/lattice-files"
const RUN_SCRIPT = "superconductors/tj/metts/run_metts_tj_auto.sh (tau=0.2, init_tau=0.1, k=3)"

# The 4x4 runs are a method benchmark, not production physics: 16 sites at
# maxdim 2000 is effectively exact, and the same system was run with TPQ
# (tj/tpq/) and DMRG (tj/dmrg_julia/) at the same couplings. They are also the
# largest single item in the project -- 80 ensembles, 7.1M samples, 310 MB,
# nearly half the t-J total -- because a 16-site chain runs the full 10000
# steps at every temperature. Excluded deliberately; they live only in
# metts/outfiles.metts, which holds nothing else.
const SKIP_LATTICES = ["square.op.nx.4.ny.4"]

const RE_PAR  = r"^t\.(-?\d+\.\d+)\.tp\.(-?\d+\.\d+)\.J\.(-?\d+\.\d+)\.holes\.(\d+)$"
const RE_T    = r"^T\.(\d+\.\d+)$"
const RE_LAT  = r"^square\.op\.nx\.(\d+)\.ny\.(\d+)$"
const RE_FILE = r"maxm\.(\d+)\.cutoff\.([0-9.eE+-]+?)\.updates\.([xz])\.nmetts\.(\d+)\.nwarm\.(\d+)\.seed\.(\d+)"

"Everything needed to convert one dump, parsed from its path."
function parse_run(p::AbstractString)
    parts = splitpath(p)
    length(parts) < 4 && return nothing
    mT = match(RE_T,    parts[end-1]); mT === nothing && return nothing
    mP = match(RE_PAR,  parts[end-2]); mP === nothing && return nothing
    mL = match(RE_LAT,  parts[end-3]); mL === nothing && return nothing
    parts[end-3] in SKIP_LATTICES && return :skip
    mF = match(RE_FILE, parts[end]);   mF === nothing && return nothing

    nx, ny = parse(Int, mL[1]), parse(Int, mL[2])
    lat = joinpath(LATDIR, parts[end-3] * ".ttpJ.lat")
    isfile(lat) || error("no lattice file '$lat' for '$p'")
    return (; path = p, lattice = lat,
            # the lattice declares T Tp J; the directory spells them t tp J
            parameters = Dict("T"  => parse(Float64, mP[1]), "Tp" => parse(Float64, mP[2]),
                              "J"  => parse(Float64, mP[3])),
            sector      = Dict("n" => nx * ny - parse(Int, mP[4])),
            temperature = parse(Float64, mT[1]),
            basis       = uppercase(mF[3]),
            algorithm = Dict{String,Any}(
                "maxdim" => parse(Int, mF[1]), "cutoff" => parse(Float64, mF[2]),
                "nmetts" => parse(Int, mF[4]), "nwarm"  => parse(Int, mF[5]),
                "seed"   => parse(Int, mF[6]),
                "tau" => 0.2, "init_tau" => 0.1, "k" => 3),
            # identifies the same physical run across the two trees
            key = (parts[end-3], parts[end-2], parts[end-1],
                   mF[1], mF[3], mF[6], mF[2]))
end

"Number of METTS steps in a dump, or 0 if it holds no states. Header read only."
function nsteps(p)
    try
        return h5open(p, "r") do h
            haskey(h, "ProductState") || return 0
            # Julia sees (nsites, nsteps)
            return size(dataspace(h["ProductState"]))[2]
        end
    catch
        return 0
    end
end

function convert_one(r, project)
    e = from_legacy_cpp(r.path; lattice = r.lattice, model = "tJ", project = project,
                        site_type = "tJ", temperature = r.temperature, basis = r.basis,
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
    latm    = get(ENV, "LATMATCH", "")

    all_dumps = String[]
    for t in TREES
        isdir(t) || (println("  (no such tree, skipping: ", t, ")"); continue)
        for (d, _, fs) in walkdir(t), f in fs
            endswith(f, ".h5") && push!(all_dumps, joinpath(d, f))
        end
    end
    sort!(all_dumps)
    println("dumps found: ", length(all_dumps))

    rs = Any[]; nskiplat = 0; nunparsed = 0
    for p in all_dumps
        r = parse_run(p)
        r === :skip     ? (nskiplat += 1) :
        r === nothing   ? (nunparsed += 1; println("  UNPARSED: ", p)) : push!(rs, r)
    end
    println("parsed: ", length(rs), "   unparsed: ", nunparsed,
            "   skipped by lattice (", join(SKIP_LATTICES, ", "), "): ", nskiplat)

    # Deduplicate: one entry per distinct run, the copy with the most steps.
    best = Dict{Any,Any}(); steps = Dict{Any,Int}()
    for r in rs
        n = nsteps(r.path)
        n == 0 && continue                      # stateless dump, nothing to convert
        if !haskey(best, r.key) || n > steps[r.key]
            best[r.key] = r; steps[r.key] = n
        end
    end
    kept = collect(values(best))
    println("distinct runs with states: ", length(kept),
            "   (dropped ", length(rs) - length(kept), " duplicate or stateless copies)")
    nb = count(r -> occursin("tj_bkup", r.path), kept)
    println("   taken from tj/: ", length(kept) - nb, "   from tj_bkup/: ", nb)

    isempty(latm) || (kept = [r for r in kept if occursin(latm, r.lattice)];
                      println("LATMATCH '", latm, "': ", length(kept), " runs"))

    # A tag must be unique inside its ensemble; check before writing anything.
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
    # Refuse a shard index outside the shard count rather than silently
    # processing everything: SHARD=18 with NSHARDS=1 means the caller believes
    # it is running one twentieth of the work and it is not.
    0 <= sh < nsh || error("SHARD=$sh is not in 0:$(nsh-1) (NSHARDS=$nsh). " *
                           "Requeueing one task must still pass the full array width.")
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
                occ = count(!=(0), e.states) / length(e.states)
                println("  ", rel)
                println("      ", nsamples(e), " samples, ", nsites(e), " sites, basis ",
                        e.collapse_bases, ", occupancy ", round(occ, digits=4),
                        " (sector ", e.sector["n"]/nsites(e), ")  <- ",
                        occursin("tj_bkup", r.path) ? "tj_bkup" : "tj")
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
