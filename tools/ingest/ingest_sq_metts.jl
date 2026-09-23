# Ingest `hubbard.square.metts` into Hubbard/hubbard.square.metts.
#
# Hubbard on square cylinders with a diagonal t', from the Flatiron project of
# the same name. Much the easiest ingest so far, because the two things that
# cost most in the triangular project are both already present here:
#
#   * **The lattices survive.** `square.op.nx.{16,32}.ny.4.hubbard.lat` are in
#     `superconductors/hubbard/lattice-files/` and declare T and Tp on
#     HUBBARDHOP bonds, exactly matching the runs' t and tp. Nothing is
#     regenerated and nothing is inferred.
#   * **The collapse basis is in the filename** as `updates.x|z`, so the
#     Sz-variance test is not needed. This is the first project in the library
#     with genuine Z-basis runs (228 of them), which also means `validate`
#     actually checks nup/ndn/sz2 on those rather than only the charge.
#
# The 4x4 lattice (400 dumps) is EXCLUDED: it is a benchmark, as in
# `superconductors`. Only nx=16 and nx=32 at ny=4 are ingested.
#
# Time evolution: metts_hubbard reads a TEBD -> TDVP schedule from a generated
# timeevofile. No run script survives for this project, so the schedule is the
# one every surviving sibling Hubbard script writes, and it is flagged as
# inferred in provenance.
#
# Not every dump holds states -- 100% of the nx=16 runs do but only ~70% of
# nx=32 -- so runs are counted after reading headers, never from the filename.

using METTSLibrary, HDF5
const ML = METTSLibrary

const ROOT   = "/data/condmat/awietek/flatiron/ceph/Research/Projects/" *
               "hubbard.square.metts/metts"
const LATDIR = "/home/awietek/Research/Projects/superconductors/hubbard/lattice-files"
const SKIP_LATTICES = ["square.op.nx.4.ny.4"]      # 4x4 benchmark
const TIME_EVOLUTION = "TEBD 0->0.1 tau 0.02 cutoff 1e-12; TDVP 0.1->beta/2 tau 0.5"
const TIME_EVOLUTION_SOURCE =
    "INFERRED: no run script survives for hubbard.square.metts; schedule taken " *
    "from the sibling Hubbard scripts (superconductors/hubbard, hubbard.polaron), " *
    "which are identical to each other"

# One uniform filename carries every parameter, so the directory layout is not
# consulted. Non-greedy groups anchored on the next key parse values that
# contain dots or a leading minus (tp.-0.25, cutoff.1e-6, T.0.0125).
const RE_FILE = r"""^outfile\.
    (square\.\w+\.nx\.\d+\.ny\.\d+)\.
    t\.(.+?)\.tp\.(.+?)\.U\.(.+?)\.holes\.(\d+)\.T\.(.+?)\.
    maxm\.(\d+)\.cutoff\.(.+?)\.updates\.([xz])\.
    nmetts\.(\d+)\.nwarm\.(\d+)\.seed\.(\d+)\.dump\.h5$"""x

const RE_LAT = r"^square\.\w+\.nx\.(\d+)\.ny\.(\d+)$"

function parse_run(p::AbstractString)
    m = match(RE_FILE, basename(p))
    m === nothing && return nothing
    name = m[1]
    name in SKIP_LATTICES && return nothing
    ml = match(RE_LAT, name); ml === nothing && return nothing
    nx, ny = parse(Int, ml[1]), parse(Int, ml[2])
    lat = joinpath(LATDIR, name * ".hubbard.lat")
    isfile(lat) || error("no lattice file '$lat' for '$p'")
    return (; path = p, lattice = lat, lattice_name = name,
            parameters = Dict("T"  => parse(Float64, m[2]),
                              "Tp" => parse(Float64, m[3]),
                              "U"  => parse(Float64, m[4])),
            sector      = Dict("n" => nx * ny - parse(Int, m[5])),
            temperature = parse(Float64, m[6]),
            basis       = uppercase(m[9]),
            algorithm = Dict{String,Any}(
                "maxdim" => parse(Int, m[7]),  "cutoff" => parse(Float64, m[8]),
                "nmetts" => parse(Int, m[10]), "nwarm"  => parse(Int, m[11]),
                "seed"   => parse(Int, m[12]),
                "tau" => 0.5, "init_tau" => 0.02,
                "time_evolution" => TIME_EVOLUTION),
            # Dedup on parsed values, never path strings: the same run is spelled
            # differently in different places (U.10 against U.10.00, cutoff.1e-6
            # against 1.0e-6) and those are equal as Float64, so a string key
            # leaves duplicates looking distinct and they then collide on one path.
            key = (name, parse(Float64, m[2]), parse(Float64, m[3]),
                   parse(Float64, m[4]), parse(Int, m[5]), parse(Float64, m[6]),
                   parse(Int, m[7]), parse(Float64, m[8]), uppercase(m[9]),
                   parse(Int, m[12])))
end

"Steps in a dump, 0 if it holds none. Header read only."
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

_nelec(col) = sum(s -> s == 0 ? 0 : (s == 3 ? 2 : 1), col)

"""
Drop samples whose particle number does not match the declared sector.

Any collapse conserves charge, so an off-sector sample is not a sample: in the
triangular project these were rows of the extensible `ProductState` dataset
allocated but never written, reading back as all-empty. Cheap to check, and
`validate` would only catch it on Z-basis runs.
"""
function drop_offsector(e)
    want = get(e.sector, "n", nothing)
    want === nothing && return e, 0
    keep = [j for j in 1:nsamples(e) if _nelec(@view e.states[:, j]) == want]
    length(keep) == nsamples(e) && return e, 0
    isempty(keep) && error("every sample off-sector in $(e.lattice_name)")
    obs = Dict{String,Array{Float64}}(k => v[keep] for (k, v) in e.observables)
    return Ensemble(model = e.model, project = e.project, site_type = e.site_type,
                    local_states = e.local_states, lattice = e.lattice,
                    lattice_name = e.lattice_name, beta = e.beta, parameters = e.parameters,
                    sector = e.sector, algorithm = e.algorithm, provenance = e.provenance,
                    collapse_bases = e.collapse_bases, states = e.states[:, keep],
                    basis = e.basis[keep], observables = obs), nsamples(e) - length(keep)
end

function convert_one(r, project)
    e0 = from_legacy_cpp(r.path; lattice = r.lattice, model = "Hubbard", project = project,
                         site_type = "Electron", temperature = r.temperature, basis = r.basis,
                         parameters = r.parameters, sector = r.sector, algorithm = r.algorithm)
    e, ndrop = drop_offsector(e0)
    ndrop > 0 && println("      dropped ", ndrop, " off-sector sample(s) from ", r.path)
    return Ensemble(model = e.model, project = e.project, site_type = e.site_type,
                    local_states = e.local_states, lattice = e.lattice,
                    lattice_name = e.lattice_name, beta = e.beta, parameters = e.parameters,
                    sector = e.sector, algorithm = e.algorithm,
                    provenance = merge(e.provenance, Dict{String,Any}(
                        "time_evolution_source" => TIME_EVOLUTION_SOURCE,
                        (ndrop > 0 ? ("dropped_offsector_samples" => ndrop,) : ())...)),
                    collapse_bases = e.collapse_bases, states = e.states, basis = e.basis,
                    observables = e.observables)
end

function main()
    out     = get(ENV, "OUT", "")
    project = get(ENV, "PROJECT", "hubbard.square.metts")
    maxf    = parse(Int, get(ENV, "MAXFILES", "0"))
    dry     = get(ENV, "DRYRUN", "0") == "1"
    resume  = get(ENV, "RESUME", "0") == "1"

    all_dumps = String[]
    for (d, _, fs) in walkdir(ROOT), f in fs
        endswith(f, ".h5") && push!(all_dumps, joinpath(d, f))
    end
    sort!(all_dumps)
    println("dumps found: ", length(all_dumps))
    # Smoke tests only: the header pass opens every dump, which is slow over
    # NFS. It changes dedup, so never set it for a real ingest.
    scanmax = parse(Int, get(ENV, "SCANMAX", "0"))
    if scanmax > 0
        all_dumps = all_dumps[1:min(end, scanmax)]
        println("SCANMAX: truncated to ", length(all_dumps), " dumps (SMOKE TEST ONLY)")
    end

    rs = Any[]; nunparsed = 0; nskipped = 0
    for p in all_dumps
        m = match(RE_FILE, basename(p))
        if m !== nothing && m[1] in SKIP_LATTICES
            nskipped += 1; continue
        end
        r = parse_run(p)
        r === nothing ? (nunparsed += 1; nunparsed <= 5 && println("  UNPARSED: ", basename(p))) :
                        push!(rs, r)
    end
    println("parsed: ", length(rs), "   unparsed: ", nunparsed,
            "   skipped (4x4 benchmark): ", nskipped)

    best = Dict{Any,Any}(); steps = Dict{Any,Int}(); nostate = 0
    for (i, r) in enumerate(rs)
        n = nsteps(r.path)
        n == 0 && (nostate += 1; continue)
        if !haskey(best, r.key) || n > steps[r.key]
            best[r.key] = r; steps[r.key] = n
        end
        i % 1000 == 0 && println("  ... scanned $i/$(length(rs))")
    end
    kept = collect(values(best))
    println("no /ProductState: ", nostate, " of ", length(rs), " parsed")
    println("distinct runs with states: ", length(kept),
            "   (dropped ", length(rs) - nostate - length(kept), " duplicates)")
    nz = count(r -> r.basis == "Z", kept)
    println("basis: X ", length(kept) - nz, "   Z ", nz)

    byens = Dict{String,Vector{String}}()
    for r in kept
        k = string(r.lattice_name, "|", sort(collect(r.parameters)), "|",
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
    0 <= sh < nsh || error("bad shard $sh of $nsh")
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
                nelec = sum(_nelec(@view e.states[:, j]) for j in 1:nsamples(e)) / nsamples(e)
                println("  ", rel)
                println("      ", nsamples(e), " samples, ", nsites(e), " sites, basis ",
                        e.collapse_bases, ", <n> per sample ", round(nelec, digits = 3),
                        " (sector n=", e.sector["n"], ")")
            else
                write_ensemble(out, e)
            end
            nok += 1
        catch err
            nerr += 1
            nerr <= 20 && println("  ERROR ", r.path, "\n        ", sprint(showerror, err))
        end
        i % 500 == 0 && !dry && println("  ... $i/$(length(sel))")
    end
    println("done: ok=", nok, " skipped=", nskip, " err=", nerr, " samples=", tot)
    nerr == 0 || exit(1)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
