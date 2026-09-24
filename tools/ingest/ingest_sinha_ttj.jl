# Ingest Aritra Sinha's t-t'-J METTS runs on square cylinders.
#
#   cluster tree: /data/condmat/asinha/Research/Projects/metts.cylinder.tj.tjp/
#                 new.metts.tj.tjp/outfiles.metts/
#   run path:     L.<L>.W.<W>/J.<J>/t.<t>/t_prime.<tp>/filling.<f>/T.<T>/D.<D>/
#                 tau.0.2.cutoff.<c>.seed.<s>/{samples.txt, outfile.h5}
#   states:       `<meas>: [1, 2, 3, ...]`, 1-based tJ indices (1 Emp, 2 Up, 3 Dn)
#   drivers:      scripts/compute_metts_t_tp_J{,_2,_3,_4,_5}.jl
#
# Everything below was read out of those five drivers, which agree on every
# point that touches stored data.
#
# BASIS is X. `collapse_with_qn!(psi, "X")` produces every stored sample. The
# only "Z" collapses are the post-DMRG one and, in v1, the warmup branch --
# neither is saved. Confirmed independently from the states: total 2Sz takes
# 16-18 distinct values across a long chain, which only an X collapse does.
#
# THE PATH LIES ABOUT tau. The output filename hardcodes `tau.0.2` whatever tau
# the run was given -- v3 even says so: "kept your legacy directory tag tau.0.2
# exactly (even if tau != 0.2)". The real value is column 8 of the params_*.txt
# file that launched it, and across the runs on disk it is 0.2 (1818), 0.3
# (952), 0.1 (115), 0.4 (89), 0.5 (13). This driver ingests ONLY the runs a
# params line proves ran at tau = 0.2; the rest are left for a later decision.
# A directory named by two params lines with different tau is dropped outright:
# the driver resumes from the last saved sample, so such a chain was continued
# at a different time step and is not one algorithm.
#
# ENERGIES ARE OFFSET BY ONE STEP, by construction. Each iteration evolves the
# product state of the PREVIOUS sample, measures and dump!s the observables,
# and only then collapses and saves the new sample -- both under index `meas`.
# So energy[j] is the energy of the METTS state grown from samples[j-1], and
# the pairing that makes the stored scalar belong to its own sample is
#
#     sample[j]  <->  energy[j+1]
#
# which costs the last sample of every run (its METTS state was never built)
# and discards energy[1] (its state was never saved). Measured over the 1,592
# selected runs that have an h5: energy rows equal the sample count in 173 and
# exceed it in 1,419, and in NOT ONE run are there fewer -- exactly what
# dump-then-save predicts, with the surplus being crashes between the two.
# Set ENERGY_SHIFT=0 to store energy[1:nsamples] unshifted instead.
#
# THE outfile.h5 ARE ENORMOUS AND IT DOES NOT MATTER. The four-fermion
# correlators write one HDF5 object per index quadruple per step, so these
# files run to 832 GB -- 112 TB across the selected runs. `h5ls` on one needs
# 16-23 GB of RSS, because it enumerates the whole group hierarchy; sixteen of
# those in parallel put the login node into swap. Reading two datasets BY NAME
# touches none of that: measured at the 1st, 25th, 50th, 75th, 90th, 97th and
# 100th percentile of file size, `h5open` + `read(h["energy"])` + `read(h["svn"])`
# costs 0.4-4.1 s and a flat 0.42 GB maxrss -- Julia's own baseline -- for files
# from 40 MB to 832 GB. So there is no size guard by default. LARGE_H5_GB is
# kept only as an escape hatch; set it to a number of GB to skip bigger files.
#
# THE LATTICE is generated, not read: the drivers build the Hamiltonian in code
# from `square_lattice(L, W; yperiodic=true)` plus a hand-written t' term.
# `verify_bonds` checks the generated bond set against the drivers' own index
# arithmetic before anything is written, because a wrong geometry is the one
# error `validate` cannot catch -- site count and coupling names would both
# still agree.

using METTSLibrary, HDF5, Printf
const ML = METTSLibrary

const ROOT = "/data/condmat/asinha/Research/Projects/metts.cylinder.tj.tjp/" *
             "new.metts.tj.tjp/outfiles.metts"
const PARAMDIR = "/home/asinha/Research/Projects/ttJ/METTS"
const DRIVER = "ttJ/METTS/scripts/compute_metts_t_tp_J_{1..5}.jl"

# Dataset name in the run's h5 -> the library's observable name. The driver
# writes `svn` for entropy_von_neumann(psi, N/2); the library calls that
# `entropy` everywhere else, and a project spelling it differently would be
# invisible to `ensembles(...)` queries and to the catalogue's observable list.
const OBS_NAMES = ("energy" => "energy", "svn" => "entropy")

const RE_LW   = r"^L\.(\d+)\.W\.(\d+)$"
const RE_J    = r"^J\.([0-9.]+)$"
const RE_T1   = r"^t\.([0-9.]+)$"
const RE_TP   = r"^t_prime\.(-?[0-9.]+)$"
const RE_FIL  = r"^filling\.([0-9.]+)$"
const RE_TEMP = r"^T\.([0-9.]+)$"
const RE_D    = r"^D\.(\d+)$"
const RE_ALG  = r"^tau\.([0-9.]+)\.cutoff\.([0-9.eE+-]+)\.seed\.(\d+)$"

# --- the params table ------------------------------------------------------

"""
Rebuild the output directory of one params line exactly as the drivers do.

All five use the same @sprintf template, tau.0.2 included, so this is how a
params line is matched to a run on disk. Any change here silently unmatches
everything, which shows up as "runs with no params line".
"""
function params_dirname(L, W, J, fil, t, tp, T, D, cutoff, seed)
    # Argument order is the drivers' own: L W J t t_prime filling T D cutoff seed.
    return @sprintf("%s/L.%d.W.%d/J.%.4f/t.%.0f/t_prime.%.4f/filling.%.5f/T.%.5f/D.%.0f/tau.0.2.cutoff.%.1e.seed.%d",
                    ROOT, L, W, J, t, tp, fil, T, D, cutoff, seed)
end

"""
    load_params(dir) -> Dict{String,NamedTuple}

Map each output directory named by a params file to its (tau, nwarm, nmetts,
source). Directories named with conflicting tau map to `nothing`, so the caller
can drop them explicitly rather than silently picking one.

Scans `dir` and one level below it. The subdirectory level is load-bearing:
`test/params.txt` alone names 80 of the 16x6 runs at tau = 0.2, and reading
only the top level silently demotes them to "no params line".

Files with `hubbard` in the name are skipped: they drive
compute_metts_hubbard_tprime0.jl into a different tree, and their column order
was never verified here.
"""
function load_params(dir::AbstractString)
    seen = Dict{String,Set{NTuple{3,Any}}}()
    src  = Dict{String,Set{String}}()
    nlines = 0
    cands = String[]
    for p in sort(readdir(dir; join=true))
        if isdir(p)
            try
                append!(cands, sort(readdir(p; join=true)))
            catch
            end
        else
            push!(cands, p)
        end
    end
    for f in cands
        base = basename(f)
        (startswith(base, "params") && endswith(base, ".txt")) || continue
        occursin("hubbard", base) && continue
        isfile(f) || continue
        for ln in eachline(f)
            c = split(strip(ln))
            length(c) == 13 || continue
            v = tryparse.(Float64, c)
            any(isnothing, v) && continue
            L, W = Int(v[1]), Int(v[2])
            J, fil, t, tp, T, tau = v[3], v[4], v[5], v[6], v[7], v[8]
            D, cutoff, seed = Int(v[9]), v[10], Int(v[11])
            nmetts, nwarm = Int(v[12]), Int(v[13])
            d = params_dirname(L, W, J, fil, t, tp, T, D, cutoff, seed)
            push!(get!(seen, d, Set{NTuple{3,Any}}()), (tau, nwarm, nmetts))
            push!(get!(src, d, Set{String}()), base)
            nlines += 1
        end
    end
    out = Dict{String,Any}()
    for (d, s) in seen
        taus = unique(x[1] for x in s)
        if length(taus) > 1
            out[d] = nothing                   # conflicting tau -- caller drops
        else
            v = first(s)
            out[d] = (tau = v[1], nwarm = v[2], nmetts = v[3],
                      params_file = join(sort(collect(src[d])), ","))
        end
    end
    println("params lines parsed: $nlines   directories named: $(length(out))")
    return out
end

# --- enumerating runs ------------------------------------------------------

"Descend the eight fixed levels, so `eigdata/` is never walked."
function find_runs(root::AbstractString)
    out = String[]
    _kids(p) = isdir(p) ? sort(readdir(p; join=true)) : String[]
    for a in _kids(root)
        match(RE_LW, basename(a)) === nothing && continue
        for b in _kids(a), c in _kids(b), d in _kids(c), e in _kids(d),
            f in _kids(e), g in _kids(f), h in _kids(g)
            match(RE_ALG, basename(h)) === nothing || push!(out, h)
        end
    end
    return out
end

function parse_run(p::AbstractString)
    q = splitpath(p)
    length(q) < 8 && return nothing
    mA = match(RE_ALG,  q[end]);   mA === nothing && return nothing
    mD = match(RE_D,    q[end-1]); mD === nothing && return nothing
    mT = match(RE_TEMP, q[end-2]); mT === nothing && return nothing
    mF = match(RE_FIL,  q[end-3]); mF === nothing && return nothing
    mP = match(RE_TP,   q[end-4]); mP === nothing && return nothing
    mt = match(RE_T1,   q[end-5]); mt === nothing && return nothing
    mJ = match(RE_J,    q[end-6]); mJ === nothing && return nothing
    mL = match(RE_LW,   q[end-7]); mL === nothing && return nothing

    L, W = parse(Int, mL[1]), parse(Int, mL[2])
    fil  = parse(Float64, mF[1])
    n    = round(Int, fil * L * W)
    return (; path = p, L = L, W = W,
            lattice_name = "square.L$(L).W$(W).cyl.ttpJ",
            parameters = Dict("J" => parse(Float64, mJ[1]),
                              "t" => parse(Float64, mt[1]),
                              "t_prime" => parse(Float64, mP[1])),
            sector = Dict("n" => n),
            temperature = parse(Float64, mT[1]),
            maxdim = parse(Int, mD[1]),
            cutoff = parse(Float64, mA[2]),
            seed   = parse(Int, mA[3]))
end

# --- the lattice -----------------------------------------------------------

"""
The bond set the drivers actually build, from their own arithmetic.

nn comes from ITensorMPS `square_lattice(L, W; yperiodic=true)`: (n, n+W) for
x < L, (n, n+1) for y < W, and the wrap (n, n+W-1) at y == 1 when W > 2. t'
comes from the hand-written loop, `s1 = (i-1)W + j` against
`s2_r = iW + mod(j,W) + 1` and `s2_l = iW + mod(j-2,W) + 1`, i running only to
L-1 so x never wraps. Returned 0-based and unordered within a bond.
"""
function driver_bonds(L::Int, W::Int)
    nn = Set{Tuple{Int,Int}}()
    for n in 1:(L * W)
        x = (n - 1) ÷ W + 1
        y = mod(n - 1, W) + 1
        x < L && push!(nn, minmax(n - 1, n + W - 1))
        if W > 1
            y < W && push!(nn, minmax(n - 1, n))
            (W > 2 && y == 1) && push!(nn, minmax(n - 1, n + W - 2))
        end
    end
    nnn = Set{Tuple{Int,Int}}()
    for i in 1:(L - 1), j in 1:W
        s1 = (i - 1) * W + j
        push!(nnn, minmax(s1 - 1, i * W + mod(j, W)))
        push!(nnn, minmax(s1 - 1, i * W + mod(j - 2, W)))
    end
    return nn, nnn
end

"Abort unless the generated lattice reproduces the drivers' bonds exactly."
function verify_bonds(L::Int, W::Int, toml::AbstractString)
    lat = parse_lattice(toml)
    gen = Dict("t" => Set{Tuple{Int,Int}}(), "J" => Set{Tuple{Int,Int}}(),
               "t_prime" => Set{Tuple{Int,Int}}())
    for (cpl, _, sites) in lat.interactions
        push!(gen[cpl], minmax(sites[1], sites[2]))
    end
    nn, nnn = driver_bonds(L, W)
    gen["t"] == nn || error("L=$L W=$W: generated t bonds differ from the driver's " *
                            "($(length(gen["t"])) vs $(length(nn)))")
    gen["J"] == nn || error("L=$L W=$W: generated J bonds differ from the driver's")
    gen["t_prime"] == nnn ||
        error("L=$L W=$W: generated t' bonds differ from the driver's " *
              "($(length(gen["t_prime"])) vs $(length(nnn)))")
    return length(nn), length(nnn)
end

const LATCACHE = Dict{Tuple{Int,Int},String}()

function lattice_for(L::Int, W::Int)
    get!(LATCACHE, (L, W)) do
        toml = square_lattice_toml(L, W; yperiodic = true,
                   bonds = [("t", "HOP", :nn), ("J", "HB", :nn),
                            ("t_prime", "HOP", :nnn)])
        n_nn, n_nnn = verify_bonds(L, W, toml)
        println("  lattice square.L$(L).W$(W).cyl.ttpJ verified: ",
                n_nn, " nn bonds, ", n_nnn, " t' bonds")
        toml
    end
end

# --- observables -----------------------------------------------------------

"""
    observables(h5path, nsamp; shift, maxgb) -> (obs, nkeep, status)

`energy` and `entropy` (the run's `svn`) for one run, aligned to its states.
Returns how many samples survive: with the one-step shift the last sample loses
its energy and is dropped, so `nkeep` can be `nsamp - 1`.

Only the two scalar datasets are read, BY NAME -- never the correlators, and
never anything that enumerates the group hierarchy. That is what keeps the cost
flat at ~0.4 GB and a few seconds on files up to 832 GB; see the header.
"""
function observables(h5path::AbstractString, nsamp::Int; shift::Bool, maxgb::Float64)
    isfile(h5path) || return Dict{String,Array{Float64}}(), nsamp, "no_outfile_h5"
    gb = filesize(h5path) / 2^30
    gb > maxgb && return Dict{String,Array{Float64}}(), nsamp,
                         "deferred_large_h5_$(round(Int, gb))GB"
    raw = try
        h5open(h5path, "r") do h
            d = Dict{String,Vector{Float64}}()
            for (src, name) in OBS_NAMES
                haskey(h, src) && (d[name] = Float64.(vec(read(h[src]))))
            end
            d
        end
    catch err
        return Dict{String,Array{Float64}}(), nsamp, "unreadable_h5"
    end
    haskey(raw, "energy") || return Dict{String,Array{Float64}}(), nsamp, "no_energy_dataset"
    nr = length(raw["energy"])

    if shift
        # sample j <-> row j+1; needs rows 2 .. nkeep+1
        nkeep = min(nsamp, nr - 1)
        nkeep >= 1 || return Dict{String,Array{Float64}}(), nsamp, "too_few_rows_$(nr)_vs_$(nsamp)"
        obs = Dict{String,Array{Float64}}()
        for (k, v) in raw
            length(v) >= nkeep + 1 && (obs[k] = v[2:(nkeep + 1)])
        end
        status = nkeep == nsamp ? "shifted_one" : "shifted_one_dropped_$(nsamp - nkeep)"
        return obs, nkeep, status
    else
        nkeep = min(nsamp, nr)
        obs = Dict{String,Array{Float64}}()
        for (k, v) in raw
            length(v) >= nkeep && (obs[k] = v[1:nkeep])
        end
        status = nkeep == nsamp ? "unshifted" : "unshifted_dropped_$(nsamp - nkeep)"
        return obs, nkeep, status
    end
end

# --- conversion ------------------------------------------------------------

function convert_one(r, pinfo, project; shift::Bool, maxgb::Float64)
    alg = Dict{String,Any}("maxdim" => r.maxdim, "cutoff" => r.cutoff, "seed" => r.seed,
                           "tau" => pinfo.tau, "nwarm" => pinfo.nwarm,
                           "nmetts" => pinfo.nmetts, "driver" => DRIVER,
                           "time_evolution" =>
                               "TDVP timeevo_tdvp_extend (applyexp, kkrylov 2, tau0 0.03, nsubdiv 2)")
    e = from_samples_txt(joinpath(r.path, "samples.txt");
        lattice = lattice_for(r.L, r.W), lattice_name = r.lattice_name,
        model = "tJ", project = project, site_type = "tJ",
        temperature = r.temperature, basis = "X",
        parameters = r.parameters, sector = r.sector, algorithm = alg)

    ns = nsamples(e)
    obs, nkeep, status = observables(joinpath(r.path, "outfile.h5"), ns;
                                     shift = shift, maxgb = maxgb)
    states = nkeep == ns ? e.states : e.states[:, 1:nkeep]
    basis  = nkeep == ns ? e.basis  : e.basis[1:nkeep]

    prov = merge(e.provenance, Dict{String,Any}(
        "source_path"      => r.path,
        "energy_alignment" => status,
        "params_file"      => pinfo.params_file,
        "samples_dropped"  => ns - nkeep))
    return Ensemble(model = e.model, project = e.project, site_type = e.site_type,
                    local_states = e.local_states, lattice = e.lattice,
                    lattice_name = e.lattice_name, beta = e.beta,
                    parameters = e.parameters, sector = e.sector, algorithm = e.algorithm,
                    provenance = prov, collapse_bases = e.collapse_bases,
                    states = states, basis = basis, observables = obs)
end

# --- main ------------------------------------------------------------------

function main()
    out     = get(ENV, "OUT", "")
    project = get(ENV, "PROJECT", "sinha.ttprime.cylinder")
    maxf    = parse(Int, get(ENV, "MAXFILES", "0"))
    dry     = get(ENV, "DRYRUN", "0") == "1"
    resume  = get(ENV, "RESUME", "0") == "1"
    shift   = get(ENV, "ENERGY_SHIFT", "1") == "1"
    # No size limit by default: a named-dataset read is O(1) in the file size.
    maxgb   = parse(Float64, get(ENV, "LARGE_H5_GB", "Inf"))
    # A run with no energy is not written; REQUIRE_ENERGY=0 to keep it anyway.
    reqE    = get(ENV, "REQUIRE_ENERGY", "1") == "1"
    wanted  = parse(Float64, get(ENV, "WANT_TAU", "0.2"))

    pars = load_params(PARAMDIR)

    runs = find_runs(ROOT)
    println("run directories on disk: ", length(runs))

    rs = Any[]; nunparsed = 0; nnoparams = 0; nconflict = 0; nskip44 = 0
    nwrongtau = 0; nempty = 0
    for p in runs
        r = parse_run(p)
        if r === nothing
            nunparsed += 1; nunparsed <= 5 && println("  UNPARSED: ", p); continue
        end
        (r.L, r.W) == (4, 4) && (nskip44 += 1; continue)
        s = joinpath(p, "samples.txt")
        (isfile(s) && filesize(s) > 0) || (nempty += 1; continue)
        if !haskey(pars, p)
            nnoparams += 1; continue
        end
        pinfo = pars[p]
        pinfo === nothing && (nconflict += 1; continue)
        pinfo.tau == wanted || (nwrongtau += 1; continue)
        push!(rs, (r, pinfo))
    end
    println("selected: ", length(rs))
    println("  skipped: unparsed=", nunparsed, " 4x4=", nskip44, " empty=", nempty,
            " no-params=", nnoparams, " tau-conflict=", nconflict,
            " tau!=", wanted, "=", nwrongtau)
    isempty(rs) && error("nothing selected")

    # Tag uniqueness inside each ensemble, before a single file is written: a
    # duplicate would make write_ensemble refuse part way through an ingest.
    byens = Dict{String,Vector{String}}()
    for (r, pinfo) in rs
        k = string(r.lattice_name, "|", sort(collect(r.parameters)), "|",
                   sort(collect(r.sector)), "|", r.temperature)
        tag = string("basis=X_maxdim=", r.maxdim, "_tau=", pinfo.tau,
                     "_cutoff=", r.cutoff, "_seed=", r.seed)
        push!(get!(byens, k, String[]), tag)
    end
    dups = [k for (k, v) in byens if length(unique(v)) < length(v)]
    println("ensembles: ", length(byens), "   with a duplicate tag: ", length(dups))
    for k in first(dups, 3); println("  DUP ", k); end
    isempty(dups) || error("duplicate tags; the tag field set is missing something")

    # Generate and verify every lattice once, serially, before sharding: two
    # array tasks racing to create the same .toml is a real hazard.
    for lw in sort(unique((r.L, r.W) for (r, _) in rs))
        lattice_for(lw[1], lw[2])
    end

    sort!(rs; by = x -> x[1].path)
    nsh = parse(Int, get(ENV, "NSHARDS", "1")); sh = parse(Int, get(ENV, "SHARD", "0"))
    0 <= sh < nsh || error("bad shard $sh of $nsh")
    nsh > 1 && (rs = rs[(sh+1):nsh:end]; println("shard $sh/$nsh: ", length(rs), " runs"))
    sel = maxf > 0 ? first(rs, maxf) : rs
    println("converting: ", length(sel), dry ? "  (DRY RUN, writing nothing)" : "  -> $out")
    println("energy pairing: ", shift ? "sample[j] <-> energy[j+1]" : "unshifted",
            "   large-h5 cutoff: ", maxgb, " GB")

    nok = 0; nerr = 0; nskipped = 0; tot = 0; ndrop = 0; nnoen = 0
    stat = Dict{String,Int}()
    for (i, (r, pinfo)) in enumerate(sel)
        try
            e = convert_one(r, pinfo, project; shift = shift, maxgb = maxgb)
            rel = ML.relpath_for(e)
            st = e.provenance["energy_alignment"]
            key = replace(String(st), r"_\d+$" => "_N")
            stat[key] = get(stat, key, 0) + 1
            # An ensemble with no energy is states without thermodynamics, and
            # the library's rule is that energy is always present. Such runs are
            # not written at all rather than written and later pruned.
            if reqE && !haskey(e.observables, "energy")
                nnoen += 1
                continue
            end
            tot += nsamples(e)
            ndrop += e.provenance["samples_dropped"]
            if resume && !dry && isfile(joinpath(out, rel))
                nskipped += 1; continue
            end
            if dry
                nholes = count(==(0), e.states) / nsamples(e)
                println("  ", rel)
                println("      ", nsamples(e), " samples, ", nsites(e), " sites, ",
                        "holes/sample ", round(nholes, digits=3),
                        " (sector n=", e.sector["n"], " => ",
                        nsites(e) - e.sector["n"], "), tau=", pinfo.tau,
                        " nwarm=", pinfo.nwarm, ", energy ", st,
                        haskey(e.observables, "entropy") ? " +entropy" : "")
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
    println("done: ok=", nok, " skipped=", nskipped, " err=", nerr,
            " no-energy-not-written=", nnoen,
            " samples=", tot, " samples-dropped-for-alignment=", ndrop)
    println("energy alignment:")
    for (k, v) in sort(collect(stat); by = x -> -x[2])
        println("   ", rpad(k, 34), v)
    end
    nerr == 0 || exit(1)
end

# Guarded, like the other ingest drivers: `include`ing this file to inspect it
# must not launch an ingest. Without the guard, an include with OUT unset
# writes ensembles relative to the current directory.
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
