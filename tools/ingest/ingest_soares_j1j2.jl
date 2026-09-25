# Ingest Rafael Soares' J1-J2 Heisenberg METTS runs on triangular cylinders.
#
#   cluster tree: /data/condmat/soares/metts_data_final/J1_J2/YC6/12_6/raw/
#   run path:     triangularJ1J2.L.<L>.W.<W>.J2.<J2>/Beta_Collapse.<beta>/
#                 tau.<tau>.cutoff.<c>.maxdim.<m>/outfile_therm.seed.<s>.hdf5
#   driver:       /home/soares/metts_share/J1_J2/YC6/L.12.W.6/main_thermalize.jl
#                 with the Hamiltonian in j1j2_triangular.jl
#
# THE MODEL is the spin-1/2 J1-J2 Heisenberg antiferromagnet on a triangular
# cylinder, open in x and periodic in y:
#
#   H = J1 sum_<ij> S_i . S_j  +  J2 sum_<<ij>> S_i . S_j
#
# with J1 = 1 fixed in the code and J2 passed per run. This is the library's
# first non-fermionic project: site type "S=1/2", local states ["Up","Dn"].
#
# THE SOURCE FORMAT IS THE CLEANEST IN THE ARCHIVE. One HDF5 per run, written
# by Dumper.jl, carrying everything needed:
#
#   /product_state  {nsteps, nsites}   1-based indices into /local_states
#   /energy         {nsteps, 1}        real(inner(psi',H,psi)), unshifted
#   /local_states   {2, 1}             ["Up","Dn"] -- the file says what its
#                                      own integers mean, so nothing is guessed
#   /beta_collapse  SCALAR             the beta of the ensemble
#   /tau_therm      SCALAR
#
# `e0_shift = -L/2` is applied to the time evolution but NOT to the stored
# energy, which is computed from H directly, so the values need no correction.
#
# THE LATTICE arrives as a Latlib.jl TOML already in this library's own format
# -- Coordinates plus Interactions with 'J1'/'J2' and type 'SdotS'. Only one
# thing must change: **its site indices are 1-based** (1..N) and the library's
# are 0-based, so every bond is shifted by one on the way in. The generated
# header calls its periodicity vectors "Torus vectors", which is misleading:
# the bond counts prove a cylinder. For L=12, W=6, N=72:
#     J1 = (L-1)W + LW + (L-1)W = 66 + 72 + 66 = 204   (a torus would give 216)
#     J2 over the three next-nearest vectors = 66 + 60 + 66 = 192
# `verify_lattice` re-derives those counts and aborts if a file disagrees.
#
# BASIS is X: collapse_with_qn(psi, "X") produces the initial state and every
# sample. Note that despite the name, this does NOT fix the magnetization --
# measured across runs, the number of "Up" sites varies from sample to sample
# (35..37 in short chains, 26..43 in long ones, about N/2). So these ensembles
# have **no conserved sector** and are stored with an empty one.
#
# ENERGIES ARE SHIFTED BY ONE STEP, as in every METTS driver in this archive.
# The loop evolves the previous sample's product state, measures the energy of
# that evolved state, and only then collapses it into the new sample. So
# energy[k] belongs to the sample saved one step earlier; the ingest pairs
# sample[j] with energy[j+1] and drops the last sample. The first energy
# belongs to the initial DMRG collapse, which is never stored.
#
# THESE ARE THERMALIZATION RUNS. main_thermalize.jl performs the thermalization
# phase only (`nmetts` is passed but unused), starting from a DMRG ground state
# collapsed in X rather than from a random product state. Alexander's call is
# that they are usable as samples regardless. The `_pruning.toml` files beside
# the data record Rafael's own thermalization windows; they are deliberately
# NOT applied here.
#
# 134 of 1,112 runs wrote a header and no steps -- `/beta_collapse`,
# `/local_states`, `/tau_therm` and nothing else. They are skipped. They
# cluster at low temperature but are not confined to it: 41 of 106 in the 18x9
# (all at beta 10, 50, 100) and 77 of 517 in the 12x6 J2=0, the latter spread
# across every beta.

using METTSLibrary, HDF5, Printf
const ML = METTSLibrary

const ROOT = "/data/condmat/soares/metts_data_final/J1_J2/YC6/12_6/raw"
const LATDIR = "/home/soares/metts_share/J1_J2/YC6/L.12.W.6"
const DRIVER = "metts_share/J1_J2/YC6/L.12.W.6/main_thermalize.jl"

const RE_SYS  = r"^triangularJ1J2\.L\.(\d+)\.W\.(\d+)\.J2\.([0-9.]+)$"
const RE_BETA = r"^Beta_Collapse\.([0-9.]+)$"
const RE_ALG  = r"^tau\.([0-9.]+)\.cutoff\.([0-9.eE+-]+)\.maxdim\.(\d+)$"
const RE_SEED = r"^outfile_therm\.seed\.(\d+)\.hdf5$"

# --- the lattice -----------------------------------------------------------

"""
    lattice_toml(L, W) -> String

Rafael's Latlib file with its site indices shifted from 1-based to 0-based.

Everything else is carried through untouched: the coordinates, the coupling
names J1/J2 and the bond type SdotS are already exactly what this library
expects, which is why this is a shift and not a conversion.
"""
function lattice_toml(L::Int, W::Int)
    src = joinpath(LATDIR, "triangular.lattice.J1.J2.L.$(L).W.$(W)_YC.toml")
    isfile(src) || error("no lattice file '$src'")
    txt = try
        read(src, String)
    catch err
        error("lattice '$src' exists but cannot be read ($(sprint(showerror, err))) " *
              "-- ask for `chmod g+r` on it")
    end
    lat = parse_lattice(txt)
    N = L * W
    nsites(lat) == N || error("lattice $(basename(src)) has $(nsites(lat)) sites, expected $N")

    flat = [s for (_, _, ss) in lat.interactions for s in ss]
    lo, hi = minimum(flat), maximum(flat)
    # The 1-based -> 0-based shift is the whole point, so verify the source
    # really is 1-based rather than trusting the generator.
    (lo == 1 && hi == N) ||
        error("$(basename(src)): site indices span $lo..$hi, expected 1..$N (1-based)")

    verify_lattice(L, W, lat)

    coords = [[lat.coordinates[d, s] for d in 1:size(lat.coordinates, 1)] for s in 1:N]
    inter = Tuple{String,String,Vector{Int}}[]
    for (cpl, typ, ss) in lat.interactions
        push!(inter, (cpl, typ, [s - 1 for s in ss]))
    end
    comments = ["# Generated by tools/ingest/ingest_soares_j1j2.jl",
                "# Source: $(basename(src)) (Latlib.jl), site indices shifted 1-based -> 0-based",
                "# Triangular cylinder L=$L W=$W, open in x, periodic in y",
                "# J1 nearest neighbours, J2 next-nearest, both type SdotS"]
    return ML._lattice_toml(comments, coords, inter)
end

"Bond counts must match a triangular cylinder, open in x and periodic in y."
function verify_lattice(L::Int, W::Int, lat)
    n1 = count(x -> x[1] == "J1", lat.interactions)
    n2 = count(x -> x[1] == "J2", lat.interactions)
    want1 = (L - 1) * W + L * W + (L - 1) * W
    want2 = (L - 1) * W + (L - 2) * W + (L - 1) * W
    n1 == want1 || error("L=$L W=$W: $n1 J1 bonds, expected $want1 " *
                         "(a torus would give $(3L*W) -- is this really a cylinder?)")
    n2 == want2 || error("L=$L W=$W: $n2 J2 bonds, expected $want2")
    types = unique(x[2] for x in lat.interactions)
    types == ["SdotS"] || error("L=$L W=$W: unexpected bond types $types")
    println("  lattice triangular.L$(L).W$(W).YC.J1J2 verified: $n1 J1 + $n2 J2 bonds on $(L*W) sites")
    return n1, n2
end

const LATCACHE = Dict{Tuple{Int,Int},String}()
lattice_for(L, W) = get!(() -> lattice_toml(L, W), LATCACHE, (L, W))

# --- one run ---------------------------------------------------------------

"""
Read a run: states as (nsites, nsteps) 0-based, energies, beta, tau.

Returns `nothing` for the header-only files that carry no steps.
"""
function read_run(path::AbstractString, N::Int)
    return h5open(path, "r") do h
        haskey(h, "product_state") || return nothing
        ls = String.(vec(read(h["local_states"])))
        want = ML.LOCAL_STATES["S=1/2"]
        ls == want || error("local_states $ls != $want in '$path'")

        ps = read(h["product_state"])
        # HDF5.jl reverses the file's dimension order; orient to (sites, steps).
        dims = size(ps)
        si = findfirst(==(N), dims)
        si === nothing && error("no dimension equals nsites=$N (dims=$dims) in '$path'")
        states_raw = si == 1 ? ps : permutedims(ps)
        nsteps = size(states_raw, 2)
        nsteps >= 1 || return nothing

        en = haskey(h, "energy") ? Float64.(vec(read(h["energy"]))) : Float64[]
        beta = Float64(read(h["beta_collapse"]))
        tau  = haskey(h, "tau_therm") ? Float64(read(h["tau_therm"])) : NaN

        for v in states_raw
            1 <= v <= length(ls) || error("state index $v out of range in '$path'")
        end
        states = UInt8.(states_raw .- 1)          # 1-based -> 0-based
        return (; states, energies = en, beta, tau, nsteps)
    end
end

function convert_one(path, project; shift::Bool)
    q = splitpath(path)
    mF = match(RE_SEED, q[end]);   mF === nothing && return nothing
    mA = match(RE_ALG,  q[end-1]); mA === nothing && return nothing
    mB = match(RE_BETA, q[end-2]); mB === nothing && return nothing
    mS = match(RE_SYS,  q[end-3]); mS === nothing && return nothing

    L, W = parse(Int, mS[1]), parse(Int, mS[2])
    J2   = parse(Float64, mS[3])
    N    = L * W

    r = read_run(path, N)
    r === nothing && return (; skip = "no_steps")
    ns = r.nsteps
    isempty(r.energies) && return (; skip = "no_energy")

    if shift
        nkeep = min(ns, length(r.energies) - 1)
        nkeep >= 1 || return (; skip = "too_few_energies")
        obs = Dict{String,Array{Float64}}("energy" => r.energies[2:(nkeep+1)])
        status = "shifted_one_dropped_$(ns - nkeep)"
    else
        nkeep = min(ns, length(r.energies))
        nkeep >= 1 || return (; skip = "too_few_energies")
        obs = Dict{String,Array{Float64}}("energy" => r.energies[1:nkeep])
        status = "unshifted_dropped_$(ns - nkeep)"
    end

    alg = Dict{String,Any}("maxdim" => parse(Int, mA[3]), "cutoff" => parse(Float64, mA[2]),
                           "tau" => r.tau, "seed" => parse(Int, mF[1]),
                           "driver" => DRIVER, "phase" => "thermalization",
                           "time_evolution" =>
                               "TDVP timeevo_tdvp_extend (applyexp, kkrylov 3, tau0 0.02, nsubdiv 2)")
    prov = Dict{String,Any}("source_path" => path, "source_format" => "soares_dumper_hdf5",
                            "code" => "METTS.jl", "energy_alignment" => status,
                            "samples_dropped" => ns - nkeep,
                            "initial_state" => "collapse of a DMRG ground state, X basis",
                            "lattice_source" => "Latlib.jl, shifted 1-based -> 0-based")
    e = Ensemble(model = "Heisenberg", project = project, site_type = "S=1/2",
                 local_states = ML.LOCAL_STATES["S=1/2"],
                 lattice = lattice_for(L, W),
                 lattice_name = "triangular.L$(L).W$(W).YC.J1J2",
                 beta = r.beta,
                 parameters = Dict("J1" => 1.0, "J2" => J2),
                 sector = Dict{String,Int}(),        # X collapse conserves no Sz here
                 algorithm = alg, provenance = prov,
                 collapse_bases = ["X"],
                 states = r.states[:, 1:nkeep], basis = zeros(UInt8, nkeep),
                 observables = obs)
    return (; ensemble = e)
end

# --- main ------------------------------------------------------------------

function find_runs(root)
    out = String[]
    kids(p) = isdir(p) ? sort(readdir(p; join=true)) : String[]
    for s in kids(root)
        match(RE_SYS, basename(s)) === nothing && continue
        for b in kids(s), a in kids(b)
            isdir(a) || continue
            for f in kids(a)
                match(RE_SEED, basename(f)) === nothing || push!(out, f)
            end
        end
    end
    return out
end

function main()
    out     = get(ENV, "OUT", "")
    project = get(ENV, "PROJECT", "soares.j1j2.triangular")
    maxf    = parse(Int, get(ENV, "MAXFILES", "0"))
    dry     = get(ENV, "DRYRUN", "0") == "1"
    resume  = get(ENV, "RESUME", "0") == "1"
    shift   = get(ENV, "ENERGY_SHIFT", "1") == "1"

    runs = find_runs(ROOT)
    println("run files: ", length(runs))
    isempty(runs) && error("found nothing under $ROOT")

    # Build and verify every lattice once, serially, before anything is
    # written: two workers racing to create the same .toml is a real hazard,
    # and an unreadable lattice should stop us now rather than half way.
    sizes = sort(unique((parse(Int, m[1]), parse(Int, m[2]))
                        for m in (match(RE_SYS, splitpath(p)[end-3]) for p in runs)
                        if m !== nothing))
    ok = Set{Tuple{Int,Int}}()
    for (L, W) in sizes
        try
            lattice_for(L, W); push!(ok, (L, W))
        catch err
            println("  LATTICE UNAVAILABLE for $(L)x$(W): ", sprint(showerror, err))
        end
    end
    runs = [p for p in runs if (m = match(RE_SYS, splitpath(p)[end-3]);
                                m !== nothing && (parse(Int, m[1]), parse(Int, m[2])) in ok)]
    println("runs with a usable lattice: ", length(runs))

    sort!(runs)
    nsh = parse(Int, get(ENV, "NSHARDS", "1")); sh = parse(Int, get(ENV, "SHARD", "0"))
    0 <= sh < nsh || error("bad shard $sh of $nsh")
    nsh > 1 && (runs = runs[(sh+1):nsh:end]; println("shard $sh/$nsh: ", length(runs)))
    sel = maxf > 0 ? first(runs, maxf) : runs
    println("converting: ", length(sel), dry ? "  (DRY RUN)" : "  -> $out")
    println("energy pairing: ", shift ? "sample[j] <-> energy[j+1]" : "unshifted")

    # Tag uniqueness inside each ensemble, before a single file is written.
    byens = Dict{String,Vector{String}}()
    for p in sel
        q = splitpath(p)
        mS = match(RE_SYS, q[end-3]); mB = match(RE_BETA, q[end-2])
        mA = match(RE_ALG, q[end-1]); mF = match(RE_SEED, q[end])
        (mS === nothing || mB === nothing || mA === nothing || mF === nothing) && continue
        k = string(mS[1], "x", mS[2], "|J2=", mS[3], "|beta=", mB[1])
        push!(get!(byens, k, String[]),
              string("maxdim=", mA[3], "_tau=", mA[1], "_cutoff=", mA[2], "_seed=", mF[1]))
    end
    dups = [k for (k, v) in byens if length(unique(v)) < length(v)]
    println("ensembles: ", length(byens), "   with a duplicate tag: ", length(dups))
    for k in first(dups, 3); println("  DUP ", k); end
    isempty(dups) || error("duplicate tags; the tag field set is missing something")

    nok = 0; nerr = 0; nskip = 0; tot = 0; ndrop = 0
    why = Dict{String,Int}()
    for (i, p) in enumerate(sel)
        try
            r = convert_one(p, project; shift = shift)
            if r === nothing || haskey(r, :skip)
                nskip += 1
                k = r === nothing ? "unparsed" : r.skip
                why[k] = get(why, k, 0) + 1
                continue
            end
            e = r.ensemble
            rel = ML.relpath_for(e)
            tot += nsamples(e); ndrop += e.provenance["samples_dropped"]
            if resume && !dry && isfile(joinpath(out, rel))
                nskip += 1; why["already_written"] = get(why, "already_written", 0) + 1
                continue
            end
            if dry
                up = count(==(0), e.states) / nsamples(e)
                println("  ", rel)
                @printf("      %d samples, %d sites, <n_up>=%.1f of %d, J2=%g beta=%g tau=%g maxdim=%d\n",
                        nsamples(e), nsites(e), up, nsites(e), e.parameters["J2"],
                        e.beta, e.algorithm["tau"], e.algorithm["maxdim"])
            else
                write_ensemble(out, e)
            end
            nok += 1
        catch err
            nerr += 1
            nerr <= 20 && println("  ERROR ", p, "\n        ", sprint(showerror, err))
        end
        i % 200 == 0 && !dry && println("  ... $i/$(length(sel))")
    end
    println("done: ok=", nok, " skipped=", nskip, " err=", nerr,
            " samples=", tot, " dropped-for-alignment=", ndrop)
    for (k, v) in sort(collect(why); by = x -> -x[2])
        println("   skip: ", rpad(k, 22), v)
    end
    nerr == 0 || exit(1)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
