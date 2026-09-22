# Ingest the optical-lattice t-J runs (cluster project "hubbard.optical.lattice",
# which is t-J despite the name) into the library.
#
#   dump path:  outfiles/<latdir>/tx.<>.ty.<>.Jx.<>.Jy.<>.holes.<>/T.<>/outfile.….dump.h5
#   lattice:    <latdir>.ttpJ.mixedd.lat, declaring TX TY JX JY
#   sector:     n = nx*ny - holes
#   basis:      updates.x|z in the filename
#   algorithm:  maxm/cutoff/nmetts/nwarm/seed from the filename, and
#               tau=0.2 init_tau=0.1 k=3 from the run script, which the dump
#               does not carry (recorded in provenance.time_evolution_source).
#
# Env: OUT (staging root), PROJECT, MAXFILES (0 = all), WORKERS, DRYRUN=1.

using METTSLibrary
using HDF5
const ML = METTSLibrary          # relpath_for is not exported

# Runs that stopped before writing any state leave a dump with no /ProductState.
# In this project that is 36 of 1360, nearly all on the 32x6 lattice at the
# lowest temperature -- the longest, most expensive runs. They are not failures
# of the ingest, so classify them up front instead of letting them surface as
# errors. Many have a chkpt.txt holding their final state, recoverable separately.
function has_states(p)
    try
        return h5open(h -> haskey(h, "ProductState"), p, "r")
    catch
        return false
    end
end

const SRC = ["/data/condmat/awietek/Research/Projects/hubbard.optical.lattice/tj/metts/outfiles",
             "/data/condmat/awietek/flatiron/ceph/Research/Projects/hubbard.optical.lattice/tj/metts/outfiles"]
const LATDIR = "/home/awietek/Research/Projects/hubbard.optical.lattice/tj/lattice-files"
const RUN_SCRIPT = "hubbard.optical.lattice/tj/metts/run_metts_tj_auto.sh (tau=0.2, init_tau=0.1, k=3)"

const RE_COUP = r"tx\.(-?\d+\.\d+)\.ty\.(-?\d+\.\d+)\.Jx\.(-?\d+\.\d+)\.Jy\.(-?\d+\.\d+)\.holes\.(\d+)"
const RE_T    = r"^T\.(\d+\.\d+)$"
const RE_LAT  = r"^square\.op\.nx\.(\d+)\.ny\.(\d+)$"
const RE_FILE = r"maxm\.(\d+)\.cutoff\.([0-9.eE+-]+?)\.updates\.([xz])\.nmetts\.(\d+)\.nwarm\.(\d+)\.seed\.(\d+)"

"Everything needed to convert one dump, parsed from its path."
function parse_run(p::AbstractString)
    parts = splitpath(p)
    f = parts[end]
    mT = match(RE_T, parts[end-1]);      mT   === nothing && return nothing
    mC = match(RE_COUP, parts[end-2]);   mC   === nothing && return nothing
    mL = match(RE_LAT, parts[end-3]);    mL   === nothing && return nothing
    mF = match(RE_FILE, f);              mF   === nothing && return nothing

    nx, ny = parse(Int, mL[1]), parse(Int, mL[2])
    holes  = parse(Int, mC[5])
    lat    = joinpath(LATDIR, parts[end-3] * ".ttpJ.mixedd.lat")
    isfile(lat) || error("no lattice file '$lat' for '$p'")

    # The cluster file is called ...ttpJ.mixedd.lat, but it declares only
    # TX TY JX JY -- there is no t'. The library uses the honest name; the
    # cluster file is left alone. Without this, lattice_name would default to
    # the .lat basename and reintroduce "ttpJ".
    return (; path = p, lattice = lat,
            lattice_name = replace(parts[end-3] * ".ttpJ.mixedd",
                                   ".ttpJ.mixedd" => ".tJ.mixedd"),
            temperature = parse(Float64, mT[1]),
            # the lattice declares TX TY JX JY; the directory spells them tx ty Jx Jy
            parameters = Dict("TX" => parse(Float64, mC[1]), "TY" => parse(Float64, mC[2]),
                              "JX" => parse(Float64, mC[3]), "JY" => parse(Float64, mC[4])),
            sector = Dict("n" => nx * ny - holes),
            basis  = uppercase(mF[3]),
            algorithm = Dict{String,Any}(
                "maxdim" => parse(Int, mF[1]), "cutoff" => parse(Float64, mF[2]),
                "nmetts" => parse(Int, mF[4]), "nwarm"  => parse(Int, mF[5]),
                "seed"   => parse(Int, mF[6]),
                "tau" => 0.2, "init_tau" => 0.1, "k" => 3))
end

function dumps()
    out = String[]
    for root in SRC
        isdir(root) || (println("  (no such source tree, skipping: ", root, ")"); continue)
        for (d, _, files) in walkdir(root), f in files
            endswith(f, ".h5") && push!(out, joinpath(d, f))
        end
    end
    return sort!(out)
end

function convert_one(r, project)
    e = from_legacy_cpp(r.path; lattice = r.lattice, lattice_name = r.lattice_name,
                        model = "tJ", project = project,
                        site_type = "tJ", temperature = r.temperature, basis = r.basis,
                        parameters = r.parameters, sector = r.sector, algorithm = r.algorithm)
    return Ensemble(model = e.model, project = e.project, site_type = e.site_type,
                    local_states = e.local_states, lattice = e.lattice,
                    lattice_name = e.lattice_name, beta = e.beta, parameters = e.parameters,
                    sector = e.sector, algorithm = e.algorithm,
                    provenance = merge(e.provenance,
                                       Dict{String,Any}("time_evolution_source" => RUN_SCRIPT)),
                    collapse_bases = e.collapse_bases, states = e.states, basis = e.basis,
                    observables = e.observables)
end

function main()
    out     = get(ENV, "OUT", "")
    # The cluster tree calls this "hubbard.optical.lattice", which is wrong twice
    # over: the runs are t-J, not Hubbard. In the library it is named for what it
    # is -- t-J with hopping and exchange tunable independently along x and y.
    # The cluster directories are NOT renamed; provenance records where it came from.
    project = get(ENV, "PROJECT", "tj.mixed.dimension")
    maxf    = parse(Int, get(ENV, "MAXFILES", "0"))
    dry     = get(ENV, "DRYRUN", "0") == "1"

    all_dumps = dumps()
    println("dumps found: ", length(all_dumps))
    rs = Any[]
    for p in all_dumps
        r = parse_run(p)
        r === nothing ? println("  UNPARSED: ", p) : push!(rs, r)
    end
    println("parsed: ", length(rs), "   unparsed: ", length(all_dumps) - length(rs))

    # Drop the stateless runs first: checking tag uniqueness over files that are
    # never written could report a collision between a real run and a skipped one.
    nostates = [r for r in rs if !has_states(r.path)]
    rs       = [r for r in rs if has_states(r.path)]
    println("with /ProductState: ", length(rs), "   stateless (skipped): ", length(nostates))

    # A tag must be unique inside its ensemble; check before writing anything.
    byens = Dict{String,Vector{String}}()
    for r in rs
        key = string(r.lattice, "|", sort(collect(r.parameters)), "|",
                     sort(collect(r.sector)), "|", r.temperature)
        tag = string("basis=", r.basis, "_maxdim=", r.algorithm["maxdim"],
                     "_tau=", r.algorithm["tau"], "_cutoff=", r.algorithm["cutoff"],
                     "_seed=", r.algorithm["seed"])
        push!(get!(byens, key, String[]), tag)
    end
    dups = [k for (k, v) in byens if length(unique(v)) < length(v)]
    println("ensembles: ", length(byens), "   with a duplicate tag: ", length(dups))
    for k in first(dups, 3); println("  DUP ", k); end
    isempty(dups) || error("duplicate tags; the tag field set is missing something")

    # LATMATCH restricts to runs whose lattice file matches a substring. Used to
    # pre-create each lattice .toml serially: write_ensemble writes it only if
    # absent, so parallel tasks starting from nothing would race on the same file.
    latm = get(ENV, "LATMATCH", "")
    isempty(latm) || (rs = [r for r in rs if occursin(latm, r.lattice)];
                      println("LATMATCH '", latm, "': ", length(rs), " runs"))

    # Shard for parallel conversion: each worker takes every NSHARDS-th run.
    # Striding rather than blocking spreads the big files evenly, since run
    # length grows with temperature and the list is sorted by path.
    nsh = parse(Int, get(ENV, "NSHARDS", "1"))
    sh  = parse(Int, get(ENV, "SHARD", "0"))
    (0 <= sh < nsh) || error("SHARD must be in 0:$(nsh-1)")
    nsh > 1 && (rs = rs[(sh+1):nsh:end])
    nsh > 1 && println("shard $sh/$nsh: ", length(rs), " runs")

    sel = maxf > 0 ? first(rs, maxf) : rs
    println("converting: ", length(sel), dry ? "  (DRY RUN, writing nothing)" : "  -> $out")

    # Resume: skip runs already written. An array job that is requeued or
    # times out must be restartable, and write_ensemble's append-only check
    # would otherwise abort the task on the first already-converted run.
    resume = get(ENV, "RESUME", "0") == "1"

    nok = 0; nerr = 0; tot = 0; nskip = 0
    for (i, r) in enumerate(sel)
        try
            e = convert_one(r, project)
            rel = ML.relpath_for(e)
            tot += nsamples(e)
            if resume && !dry && isfile(joinpath(out, rel))
                nskip += 1
                continue
            end
            if dry
                occ = sum(e.states .!= 0) / length(e.states)   # tJ: 0 = Emp
                println("  ", rel)
                println("      ", nsamples(e), " samples, ", nsites(e), " sites, basis ",
                        e.collapse_bases, ", occupancy ", round(occ, digits=4),
                        " (sector n=", e.sector["n"], " -> ", round(e.sector["n"]/nsites(e), digits=4), ")")
            else
                write_ensemble(out, e)
            end
            nok += 1
        catch err
            nerr += 1
            println("  ERROR ", r.path, "\n        ", sprint(showerror, err))
        end
        i % 100 == 0 && !dry && println("  ... $i/$(length(sel))")
    end
    println("done: ok=", nok, " skipped=", nskip, " err=", nerr, " samples=", tot)
    nerr == 0 || exit(1)
end

main()
