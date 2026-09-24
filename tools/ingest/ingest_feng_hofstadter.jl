# Ingest Chunhan Feng's Hofstadter-Hubbard METTS runs on triangular cylinders.
#
#   cluster tree: /data/condmat/chhfeng/julia_package/Task/MBD2000/
#                 Lx<Lx>Ly<Ly>/diffT/outfiles.metts/
#                 L.<L>.W.<W>/U.<U>/filling.<f>/T.<T>/D.<D>/tau.<tau>.seed.<s>/
#   per run:      samples.txt, checkpoint.bin, outfile.h5
#   driver:       /home/chhfeng/julia_package/Hubbard_Hofstadter_triangle_metts.jl
#
# THE MODEL is Hubbard on a triangular cylinder in a FIXED pi/2 flux:
#
#   H = -sum_{<ij>,s} t_ij (c+_is c_js + h.c.) + U sum_i n_iup n_idn
#
# with the Peierls factor t_ij taken from three bond lists the driver reads at
# /home/chhfeng/julia_package/Lx<Lx>Ly<Ly>_OPBC_flux_pi_2_bond_{arrow,squiggly,usual}.txt:
#
#   arrow     -> -t * exp(-i pi/2)   on the ordered pair, h.c. back
#   squiggly  -> -t * (-1)
#   usual     -> -t
#
# and t hardcoded to 1.0 inside the MPO builder (line 147), overriding the
# command-line t. `/meta/t` is the argument, not necessarily what was used;
# every production run passed 1.0 anyway.
#
# THE THREE FILES ARE PHASE CLASSES, NOT BOND DIRECTIONS. Decoded on Lx16Ly4
# with site = (x-1)*Ly + y: `arrow` is all 60 x-bonds; `usual` is 32 diagonals
# from odd columns plus 32 rungs on even columns; `squiggly` is 32 rungs on odd
# columns plus 28 diagonals from even columns. 60+64+60 = 184 = (L-1)W + LW +
# (L-1)W, i.e. a triangular cylinder open in x and periodic in y, with a
# Landau-type gauge alternating by column parity. So one geometric direction
# appears in two different files; do not read the names as directions.
#
# COMPLEX HOPPING AND THE LATTICE FORMAT. `parameters` is Dict{String,Float64},
# so exp(-i pi/2) cannot be a coupling value. The arrow bonds are therefore
# written with the bond type "HOPIM", defined as
#
#     HOPIM, coupling g, ordered pair (i,j):   -g * ( -i c+_i c_j + i c+_j c_i )
#
# which is exactly -t exp(-i pi/2) c+_i c_j + h.c. The squiggly bonds need no
# new type: phase pi is real, so they are ordinary "HOP" with coupling -1.
# This keeps couplings real and the lattice file self-describing, at the price
# of a type whose meaning is pinned to pi/2 flux. The alternative -- complex
# coupling values -- is a schema change and would mean rewriting every file
# already in the library.
#
# ENERGIES LIVE ONLY IN checkpoint.bin, as `:energies`, a Julia Serialization
# blob. The .h5 holds only /average (mean and stderr) and /meta; the per-step
# HDF5 writer exists in the driver at line 213 but its call site (line 547) is
# commented out. Measured over 205 checkpoints spanning the parameter space:
# length(energies) == Nsample == (samples.txt lines - 1) in ALL of them.
#
# ENERGIES ARE SHIFTED BY ONE STEP, as in every METTS driver in this archive.
# Each iteration evolves the product state of the PREVIOUS sample, measures the
# energy of that evolved state, and only then collapses it into samples[step].
# So energies[k] belongs to the sample saved one step earlier, and the first
# energy belongs to the last warm-up collapse, which is never saved. Pair
# sample[j] with energy[j+1]; that costs exactly one sample per run here.
#
# samples.txt CARRIES AN EXTRA LINE labelled 0: the post-DMRG Z-basis collapse
# (driver line 338). It is ground-state-directed, not thermal, and a different
# basis from everything else, so it is dropped. The remaining labels run
# Nwarm+1 .. nmetts contiguously.
#
# BASIS is X for every stored sample (collapse_with_qn!(psi,"X"), lines 435 and
# 554); the only "Z" is that unsaved label-0 collapse. Confirmed independently
# from the states: total 2Sz takes 15-20 distinct values in the long chains,
# and siteinds are conserve_qns=true, so a Z collapse would have fixed it.
#
# SAMPLING IS BIMODAL IN TEMPERATURE. nmetts/Nwarm are 15/3 at T <= 1.0 and
# 160/10 above it, so runs hold at most 12 samples below T = 1 and at most 150
# above. Of 1,309,367 samples in this tree, only 59,098 are at T <= 1.

using METTSLibrary, HDF5, Serialization, Printf
const ML = METTSLibrary

const ROOT = "/data/condmat/chhfeng/julia_package/Task/MBD2000"
const BONDDIR = "/home/chhfeng/julia_package"
const DRIVER = "julia_package/Hubbard_Hofstadter_triangle_metts.jl"

const RE_LW   = r"^L\.(\d+)\.W\.(\d+)$"
const RE_U    = r"^U\.([0-9.]+)$"
const RE_FIL  = r"^filling\.([0-9.]+)$"
const RE_T    = r"^T\.([0-9.]+)$"
const RE_D    = r"^D\.(\d+)$"
const RE_ALG  = r"^tau\.([0-9.]+)\.seed\.(\d+)$"

# --- the lattice -----------------------------------------------------------

"Read one of the driver's bond lists as 1-based ordered pairs."
function read_bonds(path::AbstractString)
    out = Tuple{Int,Int}[]
    for ln in eachline(path)
        f = split(strip(ln))
        length(f) == 2 || continue
        push!(out, (parse(Int, f[1]), parse(Int, f[2])))
    end
    return out
end

"""
    lattice_toml(Lx, Ly) -> String

Build the library lattice for one system size directly from the driver's own
bond files, so no geometry is inferred for the Hamiltonian -- only the
coordinates are, and those are checked.

Coordinates are the triangular lattice `pos(x,y) = (x-1)*a1 + (y-1)*a2` with
`a1 = (1,0)`, `a2 = (1/2, sqrt3/2)`, laid out for site `(x-1)*Ly + y`. Every
bond in the three files must then be a nearest-neighbour vector (`a1`, `a2` or
`a1-a2`), allowing the periodic wrap in y; `verify_geometry` aborts if not.
"""
function lattice_toml(Lx::Int, Ly::Int)
    N = Lx * Ly
    pre = joinpath(BONDDIR, "Lx$(Lx)Ly$(Ly)_OPBC_flux_pi_2_bond_")
    arrow    = read_bonds(pre * "arrow.txt")
    squiggly = read_bonds(pre * "squiggly.txt")
    usual    = read_bonds(pre * "usual.txt")
    nb = length(arrow) + length(squiggly) + length(usual)
    nb == 3N - 2Ly || @warn("bond count $nb is not 3N-2Ly = $(3N - 2Ly) for $(Lx)x$(Ly)")

    for (nm, bs) in (("arrow", arrow), ("squiggly", squiggly), ("usual", usual))
        for (i, j) in bs
            (1 <= i <= N && 1 <= j <= N) ||
                error("$(Lx)x$(Ly) $nm: site index out of range in ($i,$j), N=$N")
            i == j && error("$(Lx)x$(Ly) $nm: self-bond at $i")
        end
    end
    all = vcat(arrow, squiggly, usual)
    length(unique(minmax.(first.(all), last.(all)))) == nb ||
        error("$(Lx)x$(Ly): the three bond files overlap")

    verify_geometry(Lx, Ly, all)

    s3 = sqrt(3) / 2
    coords = [[Float64((n - 1) ÷ Ly) + 0.5 * ((n - 1) % Ly), s3 * ((n - 1) % Ly)]
              for n in 1:N]
    inter = Tuple{String,String,Vector{Int}}[]
    for (i, j) in arrow;    push!(inter, ("t_arrow",    "HOPIM", [i - 1, j - 1])); end
    for (i, j) in squiggly; push!(inter, ("t_squiggly", "HOP",   [i - 1, j - 1])); end
    for (i, j) in usual;    push!(inter, ("t_usual",    "HOP",   [i - 1, j - 1])); end

    comments = ["# Generated by tools/ingest/ingest_feng_hofstadter.jl",
                "# Triangular cylinder Lx=$Lx Ly=$Ly, open in x, periodic in y",
                "# Site ordering: n = (x-1)*Ly + y, 1-based in the source, 0-based here",
                "# Bonds copied verbatim from the simulation's own phase-class lists",
                "#   Lx$(Lx)Ly$(Ly)_OPBC_flux_pi_2_bond_{arrow,squiggly,usual}.txt",
                "# t_arrow  uses type HOPIM: -g(-i c+_i c_j + i c+_j c_i), i.e. flux pi/2",
                "# t_squiggly is ordinary HOP with a negative coupling (phase pi)",
                "# t_usual    is ordinary HOP (phase 0)"]
    return ML._lattice_toml(comments, coords, inter)
end

"Every bond must join nearest neighbours of the triangular lattice."
function verify_geometry(Lx::Int, Ly::Int, bonds)
    xy(n) = ((n - 1) ÷ Ly, (n - 1) % Ly)          # 0-based (x, y)
    ok = 0
    for (i, j) in bonds
        (xi, yi) = xy(i); (xj, yj) = xy(j)
        dx = xj - xi
        dy = mod(yj - yi + Ly ÷ 2, Ly) - Ly ÷ 2   # shortest wrap in y
        (dx, dy) in ((1, 0), (-1, 0), (0, 1), (0, -1), (1, -1), (-1, 1)) ||
            error("$(Lx)x$(Ly): bond ($i,$j) = ($xi,$yi)->($xj,$yj) is not nearest neighbour " *
                  "(dx=$dx dy=$dy)")
        ok += 1
    end
    return ok
end

const LATCACHE = Dict{Tuple{Int,Int},String}()
function lattice_for(Lx::Int, Ly::Int)
    get!(LATCACHE, (Lx, Ly)) do
        t = lattice_toml(Lx, Ly)
        println("  lattice triangular.Lx$(Lx).Ly$(Ly).op.fluxpi2 verified: ",
                count(==('\n'), t), " lines")
        t
    end
end

# --- one run ---------------------------------------------------------------

"Metadata straight from the run's own /meta group; the path is only a locator."
function read_meta(h5path::AbstractString)
    isfile(h5path) || return nothing
    return try
        h5open(h5path, "r") do f
            haskey(f, "meta") || return nothing
            g = f["meta"]
            d = Dict{String,Any}()
            for k in keys(g)
                d[k] = read(g[k])
            end
            d
        end
    catch
        nothing
    end
end

"""
Samples as (labels, states), with the label-0 post-DMRG Z collapse removed.

The remaining labels must be contiguous; the driver writes them by sorting a
Dict keyed by step, so a gap would mean lost steps rather than a different
convention, and silently renumbering would scramble the chain.
"""
function read_samples(path::AbstractString, nlabels::Int)
    lab = Int[]; vecs = Vector{Vector{Int}}()
    for ln in eachline(path)
        c = findfirst(':', ln)
        c === nothing && continue
        k = tryparse(Int, strip(ln[1:c-1]))
        k === nothing && continue
        body = strip(ln[c+1:end], ['[', ']', ' ', '\r'])
        isempty(body) && continue
        v = [parse(Int, strip(t)) for t in split(body, ',')]
        k == 0 && continue                      # post-DMRG Z sample: not thermal
        push!(lab, k); push!(vecs, v)
    end
    isempty(lab) && return nothing
    p = sortperm(lab); lab = lab[p]; vecs = vecs[p]
    lab == collect(lab[1]:lab[end]) ||
        error("step labels are not contiguous in '$path' (first $(lab[1]), last $(lab[end]), n=$(length(lab)))")
    N = length(vecs[1])
    all(length(v) == N for v in vecs) || error("ragged samples in '$path'")
    for v in vecs, x in v
        1 <= x <= nlabels || error("state index $x out of range in '$path'")
    end
    return lab, vecs
end

"The per-sample energy series, or nothing."
function read_energies(ckpt::AbstractString)
    isfile(ckpt) || return nothing
    d = try
        open(deserialize, ckpt)
    catch
        return nothing
    end
    d isa AbstractDict || return nothing
    e = get(d, :energies, nothing)
    e === nothing && return nothing
    return Float64.(e)
end

function convert_one(dir, project; shift::Bool)
    q = splitpath(dir)
    mA = match(RE_ALG, q[end]);   mA === nothing && return nothing
    mD = match(RE_D,   q[end-1]); mD === nothing && return nothing
    mT = match(RE_T,   q[end-2]); mT === nothing && return nothing
    mF = match(RE_FIL, q[end-3]); mF === nothing && return nothing
    mU = match(RE_U,   q[end-4]); mU === nothing && return nothing
    mL = match(RE_LW,  q[end-5]); mL === nothing && return nothing

    meta = read_meta(joinpath(dir, "outfile.h5"))
    # /meta is authoritative where it exists; the path is the fallback.
    g(k, alt) = meta === nothing ? alt : (haskey(meta, k) ? meta[k] : alt)
    Lx = Int(g("Lx", parse(Int, mL[1])))
    Ly = Int(g("Ly", parse(Int, mL[2])))
    U  = Float64(g("U", parse(Float64, mU[1])))
    fil = Float64(g("filling", parse(Float64, mF[1])))
    T  = Float64(g("T", parse(Float64, mT[1])))
    tau = Float64(g("tau", parse(Float64, mA[1])))
    D  = Int(g("maxD", parse(Int, mD[1])))
    seed = Int(g("seed", parse(Int, mA[2])))
    nwarm = Int(g("Nwarm", -1))
    N = Lx * Ly

    sm = read_samples(joinpath(dir, "samples.txt"), 4)
    sm === nothing && return nothing
    lab, vecs = sm
    ns = length(vecs)

    en = read_energies(joinpath(dir, "checkpoint.bin"))
    if en === nothing
        return (; skip = "no_energies")
    end
    if shift
        nkeep = min(ns, length(en) - 1)
        nkeep >= 1 || return (; skip = "too_few_energies")
        obs = Dict{String,Array{Float64}}("energy" => en[2:(nkeep+1)])
        status = "shifted_one_dropped_$(ns - nkeep)"
    else
        nkeep = min(ns, length(en))
        nkeep >= 1 || return (; skip = "too_few_energies")
        obs = Dict{String,Array{Float64}}("energy" => en[1:nkeep])
        status = "unshifted_dropped_$(ns - nkeep)"
    end

    states = Matrix{UInt8}(undef, N, nkeep)
    for j in 1:nkeep, i in 1:N
        states[i, j] = UInt8(vecs[j][i] - 1)      # 1-based Electron -> 0-based
    end

    alg = Dict{String,Any}("maxdim" => D, "tau" => tau, "seed" => seed,
                           "nwarm" => nwarm, "driver" => DRIVER,
                           "time_evolution" =>
                               "TDVP timeevo_tdvp_extend (applyexp, kkrylov 2, tau0 0.02, nsubdiv 2)")
    prov = Dict{String,Any}("source_path" => dir, "source_format" => "feng_metts_samples_txt",
                            "code" => "METTS.jl", "energy_alignment" => status,
                            "samples_dropped" => ns - nkeep,
                            "first_step_label" => lab[1], "last_step_label" => lab[end],
                            "flux_per_plaquette" => "pi/2")
    e = Ensemble(model = "Hubbard", project = project, site_type = "Electron",
                 local_states = ML.LOCAL_STATES["Electron"],
                 lattice = lattice_for(Lx, Ly),
                 lattice_name = "triangular.Lx$(Lx).Ly$(Ly).op.fluxpi2",
                 beta = 1 / T,
                 parameters = Dict("t_usual" => 1.0, "t_squiggly" => -1.0,
                                   "t_arrow" => 1.0, "U" => U),
                 sector = Dict("n" => round(Int, fil * N)),
                 algorithm = alg, provenance = prov,
                 collapse_bases = ["X"], states = states,
                 basis = zeros(UInt8, nkeep), observables = obs)
    return (; ensemble = e)
end

# --- main ------------------------------------------------------------------

find_runs(root) = begin
    out = String[]
    kids(p) = isdir(p) ? sort(readdir(p; join=true)) : String[]
    for a in kids(root)                      # Lx..Ly..
        for b in kids(a)                     # diffT
            for c in kids(b)                 # outfiles.metts | slurm_files
                endswith(c, "outfiles.metts") || continue
                for d in kids(c), e in kids(d), f in kids(e), g in kids(f), h in kids(g),
                    i in kids(h)
                    match(RE_ALG, basename(i)) === nothing || push!(out, i)
                end
            end
        end
    end
    return out
end

function main()
    out     = get(ENV, "OUT", "")
    project = get(ENV, "PROJECT", "feng.hofstadter.triangular")
    maxf    = parse(Int, get(ENV, "MAXFILES", "0"))
    dry     = get(ENV, "DRYRUN", "0") == "1"
    resume  = get(ENV, "RESUME", "0") == "1"
    shift   = get(ENV, "ENERGY_SHIFT", "1") == "1"
    tmax    = parse(Float64, get(ENV, "TMAX", "Inf"))   # e.g. 1.0 for the cold end only

    runs = find_runs(ROOT)
    println("run directories: ", length(runs))
    isempty(runs) && error("found nothing under $ROOT")

    if isfinite(tmax)
        runs = filter(p -> (m = match(r"/T\.([0-9.]+)/", p)) !== nothing &&
                           parse(Float64, m[1]) <= tmax, runs)
        println("after TMAX=$tmax: ", length(runs))
    end

    sort!(runs)
    nsh = parse(Int, get(ENV, "NSHARDS", "1")); sh = parse(Int, get(ENV, "SHARD", "0"))
    0 <= sh < nsh || error("bad shard $sh of $nsh")
    nsh > 1 && (runs = runs[(sh+1):nsh:end]; println("shard $sh/$nsh: ", length(runs)))
    sel = maxf > 0 ? first(runs, maxf) : runs
    println("converting: ", length(sel), dry ? "  (DRY RUN)" : "  -> $out")
    println("energy pairing: ", shift ? "sample[j] <-> energy[j+1]" : "unshifted")

    nok = 0; nerr = 0; nskip = 0; tot = 0; ndrop = 0
    why = Dict{String,Int}()
    for (i, d) in enumerate(sel)
        try
            r = convert_one(d, project; shift = shift)
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
                nel = sum(s -> s == 0 ? 0 : (s == 3 ? 2 : 1), e.states) / nsamples(e)
                println("  ", rel)
                @printf("      %d samples, %d sites, <n>=%.2f (sector n=%d), U=%g T=%g tau=%g nwarm=%s\n",
                        nsamples(e), nsites(e), nel, e.sector["n"],
                        e.parameters["U"], 1/e.beta, e.algorithm["tau"], e.algorithm["nwarm"])
            else
                write_ensemble(out, e)
            end
            nok += 1
        catch err
            nerr += 1
            nerr <= 20 && println("  ERROR ", d, "\n        ", sprint(showerror, err))
        end
        i % 500 == 0 && !dry && println("  ... $i/$(length(sel))")
    end
    println("done: ok=", nok, " skipped=", nskip, " err=", nerr,
            " samples=", tot, " dropped-for-alignment=", ndrop)
    for (k, v) in sort(collect(why); by = x -> -x[2])
        println("   skip: ", rpad(k, 24), v)
    end
    nerr == 0 || exit(1)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
