# Ingest `hubbard.triangular.metts.v2` into Hubbard/hubbard.triangular.metts.
#
# The cluster name carries a `.v2` that means nothing outside the author's own
# history, so the library drops it; provenance.source_path keeps the origin.
#
# Anisotropic triangular Hubbard cylinders: T on the two square directions,
# Tp on the diagonal, so tp=0 is the square lattice and tp=1 the isotropic
# triangular one. 36,151 dumps, by far the largest project in the archive.
#
# THREE THINGS DIFFER FROM THE OTHER DRIVERS.
#
# 1. **Parse the filename, not the path.** The four METTS subtrees use six
#    different directory shapes -- `<lat>.t.1/tp.1.U.3.holes.0/T.N`,
#    `<lat>/t.1.tp.1.U.12.holes.0/T.N`, `<lat>/t.N.tp.N.holes.N.T.N/U.N`,
#    `<lat>/t.N.tp.N.holes.N.U.N/T.N`, and two more under
#    `outfiles.metts.fixedt.save/outfiles.metts.fixedu/`. The FILENAME is
#    uniform across every one of them and carries the full parameter set, so
#    the directory layout is simply not consulted.
#
# 2. **The lattice files had to be regenerated** -- they lived in the Flatiron
#    home tree (`/mnt/home/...`, see run_dmrg_hubbard.sh) and only ceph was
#    copied here. They are produced by `tools/ingest/triangular/make_lattices.py`
#    from the original `latt.py`/`create_triangular.py`, which reproduces all
#    34 surviving hubbard.finitet triangular lattices byte-for-byte, and are
#    committed under `tools/ingest/triangular/lattices/`. `rect` is the XC
#    cylinder, plain is YC. The T/Tp split was recovered from the stored
#    energies, not assumed; see CLAUDE.md in the data repo for the numbers.
#
# 3. **The collapse basis is not in the filename.** There is no `updates.x|z`
#    here. X-basis collapse does not conserve total Sz and Z-basis does, so the
#    basis is read off the variance of the stored `TotalSz`, and runs too short
#    to judge inherit it from their parameter family.
#
# Time evolution: metts_hubbard reads a TEBD -> TDVP schedule from a generated
# timeevofile. No run script for THIS project survives, so the schedule is the
# one every surviving sibling Hubbard script writes -- TEBD 0->0.1 at tau 0.02,
# then TDVP 0.1->beta/2 at tau 0.5 -- and it is recorded as inferred.

using METTSLibrary, HDF5
const ML = METTSLibrary

# ONLY the `thermo` subtree is ingested. The project's four METTS subtrees come
# from two different run campaigns, distinguishable by the number formatting in
# their paths (`tp.1.00`/`T.0.26250` against `tp.1`/`T.0.1`), and only the later
# one dumped collapsed states at all: `outfiles.metts` and
# `outfiles.metts.fixedt.save` are 0% states across every file, and
# `outfiles.metts.fixedt` only 35% with a median chain of 27. `thermo` is 4,613
# usable runs of 5,095 with a median of 245 samples and a maximum of 9,263 --
# the best-sampled data in the archive -- so it is taken alone rather than
# diluted with short chains.
#
# Of the 482 thermo files without states, 471 contain ONLY an `InitState`
# dataset (no H, no measurements, typically 1,472 bytes): runs that died before
# completing a single METTS step. They are all ny=6, and concentrated at large
# bond dimension -- 100% of maxm=6000, 53% of maxm=4000 -- i.e. too expensive to
# reach their first sample. The other 11 are unreadable (OSError). Both are
# skipped by the `nsteps == 0` test, no special case needed.
const ROOT = "/data/condmat/awietek/flatiron/ceph/Research/Projects/" *
             "hubbard.triangular.metts.v2/metts/outfiles.metts.thermo"
const LATDIR = joinpath(@__DIR__, "triangular", "lattices")
const TIME_EVOLUTION = "TEBD 0->0.1 tau 0.02 cutoff 1e-12; TDVP 0.1->beta/2 tau 0.5"
const TIME_EVOLUTION_SOURCE =
    "INFERRED: no run script survives for hubbard.triangular.metts.v2; " *
    "schedule taken from the sibling Hubbard scripts (superconductors/hubbard, " *
    "hubbard.polaron), which are identical to each other"

# The keys are distinctive, so non-greedy groups anchored on the next key parse
# values that themselves contain dots (U.10.2, T.0.450, cutoff.1e-6).
const RE_FILE = r"""^outfile\.
    (triangular\.aniso(?:\.rect)?\.(?:op|oo|pp|po))\.
    nx\.(\d+)\.ny\.(\d+)\.
    t\.(.+?)\.tp\.(.+?)\.U\.(.+?)\.holes\.(\d+)\.T\.(.+?)\.
    maxm\.(\d+)\.cutoff\.(.+?)\.nmetts\.(\d+)\.nwarm\.(\d+)\.seed\.(\d+)
    \.dump\.h5$"""x

# The cluster's boundary codes do not say what the cylinders are. `op` and
# `rect.op` are the standard YC and XC cylinders -- established by matching bulk
# bond offsets against the surviving triangular.{XC4,YC3} lattices in
# triangular.heisenberg.dynamics, at ny=4 and ny=3 independently -- so the
# library names them YC and XC, which is the name a reader can act on. YC/XC
# already imply the cylinder, so `op` is dropped. `oo` (both open) and `pp`
# (torus) are not cylinders and keep their codes. This is a pure renaming: the
# cluster tree is untouched and provenance.source_path records the origin.
const LIBRARY_LATTICE = Dict(
    "triangular.aniso.op"      => "triangular.aniso.YC",
    "triangular.aniso.rect.op" => "triangular.aniso.XC",
    "triangular.aniso.oo"      => "triangular.aniso.oo",
    "triangular.aniso.pp"      => "triangular.aniso.pp",
    "triangular.aniso.po"      => "triangular.aniso.po",
)

function parse_run(p::AbstractString)
    m = match(RE_FILE, basename(p))
    m === nothing && return nothing
    lat, nx, ny = m[1], parse(Int, m[2]), parse(Int, m[3])
    haskey(LIBRARY_LATTICE, lat) || error("unmapped lattice type '$lat' in '$p'")
    name = "$(LIBRARY_LATTICE[lat]).nx.$(nx).ny.$(ny)"
    latfile = joinpath(LATDIR, name * ".toml")
    isfile(latfile) || error("no generated lattice '$latfile' for '$p'")
    t  = parse(Float64, m[4]);  tp = parse(Float64, m[5])
    U  = parse(Float64, m[6]);  holes = parse(Int, m[7])
    T  = parse(Float64, m[8])
    return (; path = p, lattice = latfile, lattice_name = name,
            parameters = Dict("T" => t, "Tp" => tp, "U" => U),
            sector      = Dict("n" => nx * ny - holes),
            temperature = T,
            algorithm = Dict{String,Any}(
                "maxdim" => parse(Int, m[9]),  "cutoff" => parse(Float64, m[10]),
                "nmetts" => parse(Int, m[11]), "nwarm"  => parse(Int, m[12]),
                "seed"   => parse(Int, m[13]),
                "tau" => 0.5, "init_tau" => 0.02,
                "time_evolution" => TIME_EVOLUTION),
            # Dedup key: exactly what determines the output path, on parsed
            # values. The subtrees overlap heavily (outfiles.metts.fixedt.save
            # duplicates much of outfiles.metts.fixedt) and spell the same
            # number differently -- U.10 against U.10.00, cutoff.1e-6 against
            # 1.0e-6 -- so a string key would leave duplicates looking distinct
            # and they would then collide on one path.
            key = (name, t, tp, U, holes, T,
                   parse(Int, m[9]), parse(Float64, m[10]), parse(Int, m[13])))
end

"""
Steps in a dump and the collapse basis, from one header-only pass.

Returns `(nsteps, varies)` where `varies` is true when the run's total Sz
changes across the chain.

**The evidence is asymmetric, and treating it as symmetric is wrong.** Z-basis
collapse conserves total Sz, so Sz varying *proves* the run is X. Sz staying
constant proves nothing — it is absence of evidence, not evidence of Z. A quiet
X run looks exactly like a Z run.

That is not hypothetical here. At the coldest temperature the state is
essentially a spin singlet and Sz barely moves: in one 316-run family at
T=0.0125, seeds 2 and 3 deviate from Sz=0 in a single sample out of ~27, and
seed 1 drew Sz=0 all 29 times — which at that rate happens about a third of the
time. Reading "constant" as Z tagged it Z against 315 X siblings.

So a run votes X only when it shows variation, and a family is called Z only
when nothing in it ever varies.
"""
function steps_and_basis(p)
    try
        return h5open(p, "r") do h
            haskey(h, "ProductState") || return (0, false)
            n = size(dataspace(h["ProductState"]))[2]
            n == 0 && return (0, false)
            varies = false
            if haskey(h, "TotalSz")
                sz = vec(read(h["TotalSz"]))
                varies = length(sz) >= 2 && (maximum(sz) - minimum(sz)) > 1e-6
            end
            return (n, varies)
        end
    catch
        return (0, false)
    end
end

"Particle number of one sample. Electron: 0 Emp, 1 Up, 2 Dn, 3 UpDn."
_nelec(col) = sum(s -> s == 0 ? 0 : (s == 3 ? 2 : 1), col)

"""
Drop samples whose particle number does not match the declared sector.

A collapse in ANY basis conserves charge -- an X collapse rotates spin, it does
not move electrons -- so a sample off the sector is not a sample. In practice
they are rows of zeros: the dumps write `ProductState` as an extensible
dataset, and a row that was allocated but never filled reads back as all-Emp,
i.e. n=0, which at half filling carries exactly zero weight.

Rare but real: 8 samples in 3,812,921 across 6 of 4,613 runs. Keeping them
would bias every average taken over the ensemble, and `validate` would not
catch them, because its sector check only inspects Z-basis samples.
"""
function drop_offsector(e)
    want = get(e.sector, "n", nothing)
    want === nothing && return e, 0
    keep = [j for j in 1:nsamples(e) if _nelec(@view e.states[:, j]) == want]
    length(keep) == nsamples(e) && return e, 0
    ndrop = nsamples(e) - length(keep)
    isempty(keep) && error("every sample off-sector in $(e.lattice_name)")
    obs = Dict{String,Array{Float64}}(k => v[keep] for (k, v) in e.observables)
    return Ensemble(model = e.model, project = e.project, site_type = e.site_type,
                    local_states = e.local_states, lattice = e.lattice,
                    lattice_name = e.lattice_name, beta = e.beta, parameters = e.parameters,
                    sector = e.sector, algorithm = e.algorithm, provenance = e.provenance,
                    collapse_bases = e.collapse_bases, states = e.states[:, keep],
                    basis = e.basis[keep], observables = obs), ndrop
end

function convert_one(r, project, basis)
    e0 = from_legacy_cpp(r.path; lattice = r.lattice, model = "Hubbard", project = project,
                        site_type = "Electron", temperature = r.temperature, basis = basis,
                        parameters = r.parameters, sector = r.sector, algorithm = r.algorithm)
    e, ndrop = drop_offsector(e0)
    ndrop > 0 && println("      dropped ", ndrop, " off-sector sample(s) from ", r.path)
    return Ensemble(model = e.model, project = e.project, site_type = e.site_type,
                    local_states = e.local_states, lattice = e.lattice,
                    lattice_name = e.lattice_name, beta = e.beta, parameters = e.parameters,
                    sector = e.sector, algorithm = e.algorithm,
                    provenance = merge(e.provenance, Dict{String,Any}(
                        "time_evolution_source" => TIME_EVOLUTION_SOURCE,
                        "lattice_source" => "regenerated by tools/ingest/triangular/" *
                                            "make_lattices.py; originals not on this cluster",
                        # only present when something was dropped, so its absence
                        # means the run was clean rather than unchecked
                        (ndrop > 0 ? ("dropped_offsector_samples" => ndrop,) : ())...)),
                    collapse_bases = e.collapse_bases, states = e.states, basis = e.basis,
                    observables = e.observables)
end

"Family key for propagating the basis: everything but the seed and bond dimension."
family(r) = (r.lattice_name, r.parameters["T"], r.parameters["Tp"],
             r.parameters["U"], r.sector["n"])

function main()
    out     = get(ENV, "OUT", "")
    project = get(ENV, "PROJECT", "hubbard.triangular.metts")
    maxf    = parse(Int, get(ENV, "MAXFILES", "0"))
    dry     = get(ENV, "DRYRUN", "0") == "1"
    resume  = get(ENV, "RESUME", "0") == "1"

    all_dumps = String[]
    for (d, _, fs) in walkdir(ROOT), f in fs
        endswith(f, ".h5") && push!(all_dumps, joinpath(d, f))
    end
    sort!(all_dumps)
    println("dumps found: ", length(all_dumps))
    # SCANMAX truncates the enumeration for smoke tests only. The header pass
    # opens every dump, which over 36k multi-GB files on NFS is hours, so a dry
    # run that only needs to exercise the parser and one conversion sets this.
    # It changes dedup and basis voting, so never set it for a real ingest.
    scanmax = parse(Int, get(ENV, "SCANMAX", "0"))
    if scanmax > 0
        all_dumps = all_dumps[1:min(end, scanmax)]
        println("SCANMAX: truncated to ", length(all_dumps), " dumps (SMOKE TEST ONLY)")
    end

    rs = Any[]; nunparsed = 0
    for p in all_dumps
        r = parse_run(p)
        r === nothing ? (nunparsed += 1; nunparsed <= 5 && println("  UNPARSED: ", basename(p))) :
                        push!(rs, r)
    end
    println("parsed: ", length(rs), "   unparsed: ", nunparsed)

    # Header pass: step counts, basis votes, and deduplication in one sweep.
    best = Dict{Any,Any}(); steps = Dict{Any,Int}(); bas = Dict{Any,String}()
    votes = Dict{Any,Dict{String,Int}}()
    nostate = 0
    for (i, r) in enumerate(rs)
        n, varies = steps_and_basis(r.path)
        # Most dumps in this project store measurements but NOT the collapsed
        # states -- counting runs by filename overstates it badly. Report that
        # separately from deduplication, they are different facts.
        n == 0 && (nostate += 1; continue)
        varies && (bas[r.key] = "X")
        v = get!(votes, family(r), Dict{String,Int}())
        v[varies ? "X" : "quiet"] = get(v, varies ? "X" : "quiet", 0) + 1
        if !haskey(best, r.key) || n > steps[r.key]
            best[r.key] = r; steps[r.key] = n
        end
        i % 2000 == 0 && println("  ... scanned $i/$(length(rs))")
    end
    kept = collect(values(best))
    println("no /ProductState (measurements only): ", nostate,
            "   of ", length(rs), " parsed")
    println("distinct runs with states: ", length(kept),
            "   (dropped ", length(rs) - nostate - length(kept), " duplicates)")

    # A family is X if ANY of its runs shows Sz varying, because that is proof;
    # it is Z only when every run in it is quiet. A quiet run in an X family is
    # an X run that happened not to fluctuate, not a Z run.
    famX = Set(k for (k, v) in votes if get(v, "X", 0) > 0)
    println("families: ", length(votes), "   X (some run varies): ", length(famX),
            "   all-quiet -> Z: ", length(votes) - length(famX))
    for k in first([k for (k, v) in votes if !(k in famX)], 3)
        println("  ALL-QUIET (called Z) ", k, " -> ", votes[k])
    end

    resolved = Dict{Any,String}(); nproved = 0; nfam = 0
    for r in kept
        if haskey(bas, r.key)
            resolved[r.key] = "X"; nproved += 1           # Sz varied: proof
        else
            resolved[r.key] = family(r) in famX ? "X" : "Z"; nfam += 1
        end
    end
    println("basis: proved X ", nproved, ", inherited from family ", nfam)
    nx = count(r -> get(resolved, r.key, "") == "X", kept)
    println("       X: ", nx, "   Z: ", length(resolved) - nx)

    byens = Dict{String,Vector{String}}()
    for r in kept
        haskey(resolved, r.key) || continue
        k = string(r.lattice_name, "|", sort(collect(r.parameters)), "|",
                   sort(collect(r.sector)), "|", r.temperature)
        tag = string("basis=", resolved[r.key], "_maxdim=", r.algorithm["maxdim"],
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
            b = get(resolved, r.key, nothing)
            b === nothing && (nskip += 1; continue)
            e = convert_one(r, project, b)
            rel = ML.relpath_for(e)
            tot += nsamples(e)
            if resume && !dry && isfile(joinpath(out, rel))
                nskip += 1; continue
            end
            if dry
                nelec = sum(s -> s == 0 ? 0 : (s == 3 ? 2 : 1), e.states) / nsamples(e)
                println("  ", rel)
                println("      ", nsamples(e), " samples, ", nsites(e), " sites, basis ",
                        e.collapse_bases, ", <n> per sample ", round(nelec, digits=2),
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

# Run only when this file IS the program, so the driver can also be `include`d
# to reuse parse_run/convert_one on a single dump -- regenerating one file
# should not cost a 5,095-file header scan.
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
