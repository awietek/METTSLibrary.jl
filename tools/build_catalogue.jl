#  Build an ensemble-level catalogue of every ingestable METTS run.
#  One row per run file; then grouped into ensembles
#  (model, lattice, parameters, sector, temperature).
using Printf

const W = get(ENV, "METTS_SCAN_DIR", joinpath(pwd(), "scan"))

mutable struct Run
    family::String; project::String; fmt::String
    model::String; site_type::String; lattice::String
    params::Vector{Pair{String,Float64}}
    sector::Vector{Pair{String,Int}}
    T::Float64; seed::Int; basis::String
    algo::Vector{Pair{String,Any}}
    nsites::Int; nsteps::Int; nscalars::Int
    path::String; bytes::Int
end

num(s) = parse(Float64, s)

# helper: pull "name.value" pairs like t.1.00.tp.0.20.J.0.40.holes.12 out of a dir
function kvdir(s::AbstractString)
    out = Pair{String,Float64}[]
    toks = split(s, '.')
    i = 1
    while i <= length(toks)
        name = toks[i]
        occursin(r"^[A-Za-z][A-Za-z_]*$", name) || (i += 1; continue)
        j = i + 1; digits = String[]
        while j <= length(toks) && occursin(r"^-?\d+$", toks[j])
            push!(digits, toks[j]); j += 1
            length(digits) == 2 && break          # "1.00" -> two tokens
        end
        isempty(digits) && (i += 1; continue)
        push!(out, name => num(join(digits, ".")))
        i = j
    end
    return out
end

seed_of(f)  = (m = match(r"seed\.(\d+)", f);      m === nothing ? -1 : parse(Int, m[1]))
maxm_of(f)  = (m = match(r"maxm\.(\d+)", f);      m === nothing ? 0  : parse(Int, m[1]))
cut_of(f)   = (m = match(r"cutoff\.([0-9.eE+-]+?)(?:\.|$)", f); m === nothing ? 0.0 : num(m[1]))
basis_of(f) = (m = match(r"updates\.([xz])", f);  m === nothing ? "?" : uppercase(m[1]))
nmetts_of(f)= (m = match(r"nmetts\.(\d+)", f);    m === nothing ? 0  : parse(Int, m[1]))
nwarm_of(f) = (m = match(r"nwarm\.(\d+)", f);     m === nothing ? 0  : parse(Int, m[1]))

# Locate the lattice directory by name rather than by offset: the trees differ
# in depth (kagome has an extra "hubbard/" level; superconductors has both
# depth 8 and depth 9), and a fixed offset silently grabs the wrong component.
# A lattice directory always carries a size (nx.32, L.12). Without that the
# project directories themselves match -- "kagome.superconductors" and
# "triangular.heisenberg.dynamics" both start with a lattice-looking prefix.
const LAT_RE = r"^(square|triangular|ladder|kagome|hubbard\.kanamori)\..*?(nx\.\d+|L\.\d+)"
function lat_index(p)
    for i in 2:length(p)-1          # never the project dir, never the filename
        occursin(LAT_RE, p[i]) && return i
    end
    return 0
end

# parse one path below /Projects/ ; returns (family, lattice, params, sector, T) or nothing
function parse_path(rel::String, nsites::Int)
    p = split(rel, '/'); f = p[end]
    proj = p[1]
    T = (m = match(r"/T\.([0-9.]+)/", "/" * rel); m === nothing ? NaN : num(m[1]))
    li = lat_index(p)
    li == 0 && return nothing

    if proj == "superconductors" || proj == "hubbard.square.metts" ||
       proj == "hubbard.polaron" || proj == "hubbard.optical.lattice"
        # <lattice>/<couplings+holes>/T.x/outfile...
        lat = String(p[li]); kv = kvdir(p[li+1])
        holes = 0; par = Pair{String,Float64}[]
        for (k, v) in kv
            k == "holes" ? (holes = Int(v)) : push!(par, k => v)
        end
        # The project name does not give the model: "superconductors" holds both
        # tj/ and hubbard/ subtrees, and "hubbard.optical.lattice" is entirely
        # t-J. The couplings decide -- a Hubbard U means Electron, an exchange
        # J (or Jx/Jy) with no U means t-J. Verified against the dumps:
        # the Hubbard ones carry DoubleOcc/Nupdn, the t-J ones PairingCorrelation.
        hasU = any(x -> first(x) == "U", par)
        hasJ = any(x -> startswith(first(x), "J"), par)
        model = hasU ? "Hubbard" : (hasJ ? "tJ" : "Hubbard")
        st    = model == "tJ" ? "tJ" : "Electron"
        return (proj, model, st, lat, par, ["n" => nsites - holes], T)

    elseif proj == "hubbard.triangular.metts.v2"
        # <lattice>/t.1.00.tp.1.00.holes.0.T.0.10000/U.10.60/outfile...
        # Two path shapes in this family:
        #   <lat>/t.1.00.tp.1.00.holes.0.T.0.10000/U.10.60/...   T in the couplings dir
        #   <lat>/t.1.00.tp.1.00.U.10.60.holes.0/T.0.10000/...   T in its own dir
        # so "T" can turn up in either component. It is never a coupling here;
        # filter it wherever it appears rather than only in the first one.
        lat = String(p[li])
        comps = String[String(p[li+1])]
        li + 2 <= length(p) - 1 && push!(comps, String(p[li+2]))
        holes = 0; par = Pair{String,Float64}[]; Tloc = T
        for c in comps, (k, v) in kvdir(c)
            if k == "holes"
                holes = Int(v)
            elseif k == "T"
                Tloc = v
            else
                push!(par, k => v)
            end
        end
        return (proj, "Hubbard", "Electron", lat, par, ["n" => nsites - holes], Tloc)

    elseif proj == "triangular.heisenberg.dynamics"
        lat = String(p[li])
        return (proj, "Heisenberg", "S=1/2", lat, ["J" => 1.0], Pair{String,Int}[], T)

    elseif proj == "hubbard.bfield"
        # <lattice>/t.1.00.U.8.00/holes.4/T.0.65/eta.0.005.B.0.0/outfile...
        lat = String(p[li]); par = kvdir(p[li+1])
        holes = (m = match(r"holes\.(\d+)", p[li+2]); m === nothing ? 0 : parse(Int, m[1]))
        append!(par, kvdir(p[li+4]))
        Tloc = (m = match(r"T\.([0-9.]+)", p[li+3]); m === nothing ? T : num(m[1]))
        return (proj, "Hubbard", "Electron", lat, par, ["n" => nsites - holes], Tloc)

    elseif proj == "kagome.superconductors"
        # .../<lattice>/<params>/nparticles.N/T.x/tau...maxdim.D/outfile.seed.S.samples.txt
        # NB two different physical systems live under this project directory:
        #   hubbard/data/kagome.*        genuine kagome cylinders (t, U, V)
        #   data/hubbard.kanamori.*      SQUARE cylinders, m orbitals (t, U, Uprime, J)
        # The latter is written here by hubbard.kanamori/run.sh; it is misfiled.
        latdir = String(p[li]); par = kvdir(p[li+1])
        np  = (mm = match(r"nparticles\.(\d+)", p[li+2]); mm === nothing ? 0 : parse(Int, mm[1]))
        Tloc = (mm = match(r"T\.([0-9.]+)", p[li+3]); mm === nothing ? T : num(mm[1]))
        if startswith(latdir, "hubbard.kanamori")
            mm = match(r"L\.(\d+)\.W\.(\d+)", latdir)
            L = mm === nothing ? 0 : parse(Int, mm[1]); Wd = mm === nothing ? 0 : parse(Int, mm[2])
            morb = 0
            par2 = Pair{String,Float64}[]
            for (k, v) in par
                k == "m" ? (morb = Int(v)) : push!(par2, k => v)
            end
            lat = "square.L$(L).W$(Wd).m$(morb).cyl"
            return ("hubbard.kanamori", "HubbardKanamori", "Electron", lat, par2, ["n" => np], Tloc)
        else
            mm = match(r"L\.(\d+)\.W\.(\d+)", latdir)
            L = mm === nothing ? 0 : parse(Int, mm[1]); Wd = mm === nothing ? 0 : parse(Int, mm[2])
            lat = "kagome.L$(L).W$(Wd).cyl"
            return ("kagome.superconductors", "Hubbard", "Electron", lat, par, ["n" => np], Tloc)
        end
    end
    return nothing
end

# Collapse basis inferred from the states themselves: X-basis collapse does not
# conserve total Sz, Z-basis does. Calibrated against 113 files whose filename
# records the basis: 80/80 updates.x showed Sz varying, 33/33 updates.z fixed.
# Sz fixed on fewer than 20 samples stays ambiguous -- it can happen by chance.
const SZINFER = Dict{String,String}()
let p = joinpath(W, "sz.tsv")
    if isfile(p)
        for l in eachline(p)
            f = split(l, '\t'); length(f) < 10 && continue
            f[3] == "OK" || continue
            k = parse(Int, f[6]); ndist = parse(Int, f[7])
            SZINFER[String(f[1])] = ndist > 1 ? "X" : (k >= 20 ? "Z" : "?")
        end
    end
end
basis_final(path, f) = (b = basis_of(f); b != "?" ? b : get(SZINFER, path, "?"))

# Provenance: the directory tree a run came from, between the project and the
# lattice directory. This is what separates live from backup subtrees, tj/ from
# hubbard/, and one outfiles.* variant from another -- distinctions the physics
# key deliberately throws away.
function prov_of(path::AbstractString)
    i = findfirst("/Projects/", path)
    i === nothing && return "?"
    p = split(String(path[last(i)+1:end]), '/')
    li = lat_index(p)
    li <= 1 && return String(p[1])
    return join(p[1:li-1], "/")
end

runs = Run[]
unparsed = String[]

# ---- HDF5 runs with /ProductState -----------------------------------------
for line in eachline(joinpath(W, "classified.tsv"))
    fs = split(line, '\t'); length(fs) < 7 && continue
    fs[2] == "STATES" || continue
    path = String(fs[1]); bytes = parse(Int, fs[3])
    nsites = parse(Int, fs[4]); nsteps = parse(Int, fs[5]); nsc = parse(Int, fs[6])
    i = findfirst("/Projects/", path); i === nothing && (push!(unparsed, path); continue)
    rel = String(path[last(i)+1:end])
    r = parse_path(rel, nsites)
    r === nothing && (push!(unparsed, path); continue)
    (proj, model, st, lat, par, sec, T) = r
    f = basename(path)
    push!(runs, Run(proj, proj, "legacy_cpp", model, st, lat, par, sec, T,
        seed_of(f), basis_final(path, f),
        ["maxdim" => maxm_of(f), "cutoff" => cut_of(f), "nmetts" => nmetts_of(f), "nwarm" => nwarm_of(f)],
        nsites, nsteps, nsc, path, bytes))
end

# ---- kagome samples.txt ----------------------------------------------------
for line in eachline(joinpath(W, "mates2.tsv"))
    fs = split(line, '\t'); length(fs) < 2 && continue
    fs[1] == "LABELLED" || continue
    path = String(fs[2])
    nsteps = 0; nsites = 0
    try
        for l in eachline(path)
            m = match(r"\[(.*)\]", l); m === nothing && continue
            nsteps += 1
            nsites == 0 && (nsites = count(==(','), m[1]) + 1)
        end
    catch; end
    nsteps == 0 && (push!(unparsed, path); continue)
    i = findfirst("/Projects/", path); i === nothing && (push!(unparsed, path); continue)
    rel = String(path[last(i)+1:end])
    r = parse_path(rel, nsites)
    r === nothing && (push!(unparsed, path); continue)
    (proj, model, st, lat, par, sec, T) = r
    f = basename(path)
    d = dirname(path)
    # non-greedy and anchored: "tau.0.1.cutoff.1e-6.maxdim.1000" must not
    # capture the trailing dot ("0.1." fails to parse as Float64)
    tau = (mm = match(r"tau\.([0-9.]+?)\.cutoff", d); mm === nothing ? 0.0 : num(mm[1]))
    md  = (mm = match(r"maxdim\.(\d+)", d);           mm === nothing ? 0   : parse(Int, mm[1]))
    cf  = (mm = match(r"cutoff\.([0-9.eE+-]+?)\.maxdim", d); mm === nothing ? 0.0 : num(mm[1]))
    push!(runs, Run(proj, proj, "metts_jl_samples_txt", model, st, lat, par, sec, T,
        seed_of(f), "X", ["tau" => tau, "maxdim" => md, "cutoff" => cf],
        nsites, nsteps, 1, path, filesize(path)))
end

# ---- checkpoint states (one final, fully thermalised state per run) --------
# Only those whose chain is not already stored: the exact redundancy test lives
# in chkpt_unique.paths (written as the *dump* path of the matching seed).
for line in eachline(joinpath(W, "names.paths"))
    path = String(strip(line)); isempty(path) && continue
    isfile(path) || (push!(unparsed, path); continue)
    state = ""
    try
        for l in eachline(path)
            if occursin(r"^(Up|Dn|Emp|UpDn)([ \t]+(Up|Dn|Emp|UpDn)){7,}", l); state = l; break; end
        end
    catch; end
    isempty(state) && (push!(unparsed, path); continue)
    nsites = length(split(state))
    i = findfirst("/Projects/", path); i === nothing && (push!(unparsed, path); continue)
    rel = String(path[last(i)+1:end])
    r = parse_path(rel, nsites)
    r === nothing && (push!(unparsed, path); continue)
    (proj, model, st, lat, par, sec, T) = r
    f = basename(path)
    push!(runs, Run(proj, proj, "legacy_cpp_chkpt", model, st, lat, par, sec, T,
        seed_of(f), basis_of(f),
        ["maxdim" => maxm_of(f), "cutoff" => cut_of(f), "nmetts" => nmetts_of(f), "nwarm" => nwarm_of(f)],
        nsites, 1, 0, path, filesize(path)))
end

# ---- propagate the basis within each parameter family ----------------------
# The collapse basis is a property of the run script, and a temperature sweep is
# one script invoked repeatedly, so every temperature in a family was collapsed
# the same way. Where a family determines a single basis anywhere, adopt it for
# the runs whose own samples were too few to decide.
const BASIS_SRC = Dict{String,String}()
for r in runs
    r.basis == "?" && continue
    BASIS_SRC[r.path] = basis_of(basename(r.path)) != "?" ? "file" : "sz"
end
famkey(r) = join([r.project, r.lattice,
    join(["$k=$v" for (k, v) in sort(r.params, by=first)], "_"),
    join(["$k=$v" for (k, v) in sort(r.sector, by=first)], "_")], "|")
function propagate_basis!(runs, src)
    fam = Dict{String,Set{String}}()
    for r in runs
        r.basis == "?" && continue
        push!(get!(fam, famkey(r), Set{String}()), r.basis)
    end
    nprop = 0; nconf = 0; norphan = 0
    for r in runs
        r.basis == "?" || continue
        s = get(fam, famkey(r), Set{String}())
        if length(s) == 1
            r.basis = first(s); src[r.path] = "family"; nprop += 1
        elseif length(s) > 1
            nconf += 1
        else
            norphan += 1
        end
    end
    return (nprop, nconf, norphan)
end
let (a, b, c) = propagate_basis!(runs, BASIS_SRC)
    @printf("basis propagated within families: %d runs resolved, %d in mixed-basis families, %d with no determined sibling\n", a, b, c)
end

# ---- drop checkpoints whose chain is already stored -------------------------
# Matched on parsed values, not filename text: the dump and its own checkpoint
# format numbers differently (t.1.00.tp.0.00...T.0.01250 vs t.1.tp.0...T.0.0125),
# so a string substitution silently fails to pair them.
maxd_of(x) = (for (kk, vv) in x.algo; kk == "maxdim" && return round(Int, vv); end; 0)
runkey(x) = join([x.model, x.lattice,
    join(["$k=$v" for (k, v) in sort(x.params, by=first)], "_"),
    join(["$k=$v" for (k, v) in sort(x.sector, by=first)], "_"),
    @sprintf("%.5f", x.T), string(x.seed), string(maxd_of(x)), x.basis], "|")
function drop_redundant_checkpoints!(rs)
    have = Set(runkey(x) for x in rs if x.fmt != "legacy_cpp_chkpt")
    before = length(rs)
    filter!(x -> x.fmt != "legacy_cpp_chkpt" || !(runkey(x) in have), rs)
    return before - length(rs)
end
let d = drop_redundant_checkpoints!(runs)
    @printf("checkpoints dropped as redundant: %d;  checkpoints kept: %d\n",
            d, count(x -> x.fmt == "legacy_cpp_chkpt", runs))
end

@printf("parsed runs: %d   unparsed: %d\n", length(runs), length(unparsed))
if !isempty(unparsed)
    println("first unparsed:")
    for u in first(unparsed, 5); println("   ", u); end
end

# ---- group into ensembles --------------------------------------------------
key(r) = join([r.model, r.lattice,
    join(["$k=$v" for (k, v) in sort(r.params, by=first)], "_"),
    join(["$k=$v" for (k, v) in sort(r.sector, by=first)], "_"),
    @sprintf("T=%.5f", r.T)], "|")

ens = Dict{String,Vector{Run}}()
for r in runs; push!(get!(ens, key(r), Run[]), r); end
@printf("ensembles: %d\n", length(ens))

const COLLIDE = open(joinpath(W, "collisions.tsv"), "w")
open(joinpath(W, "catalogue.tsv"), "w") do io
    println(io, join(["key","project","format","model","site_type","lattice","params","sector",
                      "T","beta","nruns","nsamples","nsites","seeds","bases","scalars",
                      "raw_bytes","est_bytes","flags","provenance","maxdims","example_tag"], "\t"))
    for (k, rs) in sort(collect(ens), by=first)
        r0 = rs[1]
        nsamp = sum(x.nsteps for x in rs)
        raw   = sum(x.bytes  for x in rs)
        est   = sum(0.242*x.nsites*x.nsteps + 4.4*x.nscalars*x.nsteps + 4*x.nsteps + 5000 for x in rs)
        bases = join(sort(unique(x.basis for x in rs)), ",")
        seeds = join(sort(unique(x.seed for x in rs)), ",")
        flags = String[]
        occursin("?", bases) && push!(flags, "basis_ambiguous")
        # where each run's basis came from: the filename, its own Sz, or its family
        any(x -> get(BASIS_SRC, x.path, "") == "sz", rs)     && push!(flags, "basis_inferred")
        any(x -> get(BASIS_SRC, x.path, "") == "family", rs) && push!(flags, "basis_from_family")
        occursin("Z", bases) && occursin("X", bases) && push!(flags, "mixed_basis")
        r0.project == "hubbard.kanamori" && push!(flags, "lattice_built_in_code")
        r0.project == "hubbard.kanamori" && push!(flags, "misfiled_under_kagome")
        r0.project == "kagome.superconductors" && push!(flags, "lattice_needs_conversion")
        length(unique(x.nsites for x in rs)) > 1 && push!(flags, "mixed_nsites")
        isempty(r0.sector) && push!(flags, "no_sector")
        any(occursin("_bkup/", x.path) for x in rs) && push!(flags, "backup_subtree")
        # The ensemble key is physics (model, lattice, couplings, sector, T), so
        # identical systems run under different project directories land in the
        # same ensemble -- and would collide on <seed>.h5 in the library.
        length(unique(x.project for x in rs)) > 1 && push!(flags, "multi_project")
        # Runs in one ensemble often differ only in bond dimension, so a seed
        # number alone does not identify a run: seed.1 at maxm 2000 and at 4000
        # are different runs. The library tag has to carry both.
        # algo holds Any values and mixes Int with Float, so maxdim arrives as
        # 2000.0; force it back to an integer for names and sorting.
        maxd(x) = (for (kk, vv) in x.algo; kk == "maxdim" && return round(Int, vv); end; 0)
        maxdims = sort(unique(maxd(x) for x in rs))
        length(maxdims) > 1 && push!(flags, "multi_maxdim")
        # Library tag: seed<n>_maxm<m>_<basis>_<provenance>. Each component is
        # forced by a real collision: the seed repeats across bond dimensions,
        # (seed, maxdim) repeats across provenance trees, and the X/Z comparison
        # runs in hubbard.square.metts repeat everything but the collapse basis.
        tag(x) = string("seed", x.seed, "_maxm", maxd(x), "_", x.basis, "_",
                        replace(prov_of(x.path), '/' => '.'))
        tags = [tag(x) for x in rs]
        if length(unique(tags)) < length(rs)
            # Genuinely distinct runs that every parameter in the name agrees on
            # (sibling directories T.0.1000 / T.0.10 with different lengths).
            # Suffix deterministically by source path so the mapping is stable.
            push!(flags, "tag_disambiguated")
            counts = Dict{String,Int}()
            for i in sortperm([x.path for x in rs])
                c = get(counts, tags[i], 0) + 1
                counts[tags[i]] = c
                c > 1 && (tags[i] = string(tags[i], "_v", c))
            end
            for (t, ps) in Dict(t => [rs[j].path for j in eachindex(rs) if tag(rs[j]) == t]
                                for t in unique(tag(x) for x in rs))
                length(ps) > 1 && println(COLLIDE, k, "\t", t, "\t", join(ps, "\t"))
            end
        end
        length(unique(tags)) < length(rs) && push!(flags, "tag_collision")
        length(unique((x.seed, maxd(x)) for x in rs)) < length(rs) &&
            push!(flags, "duplicate_seeds")
        all(x.fmt == "legacy_cpp_chkpt" for x in rs) && push!(flags, "checkpoint_only")
        any(x.fmt == "legacy_cpp_chkpt" for x in rs) &&
            !all(x.fmt == "legacy_cpp_chkpt" for x in rs) && push!(flags, "has_checkpoints")
        nsamp < 100 && push!(flags, "few_samples")
        println(io, join([k, join(sort(unique(x.project for x in rs)), "+"), r0.fmt,
            r0.model, r0.site_type, r0.lattice,
            join(["$a=$b" for (a,b) in sort(r0.params, by=first)], " "),
            join(["$a=$b" for (a,b) in sort(r0.sector, by=first)], " "),
            @sprintf("%.5f", r0.T), @sprintf("%.4f", 1/r0.T), length(rs), nsamp, r0.nsites,
            seeds, bases, r0.nscalars, raw, round(Int, est), join(flags, ","),
            join(sort(unique(prov_of(x.path) for x in rs)), "+"),
            join(maxdims, ","), first(sort(tags))], "\t"))
    end
end
println("wrote catalogue.tsv")

# ---- summary ---------------------------------------------------------------
println("\n=== ensembles per project")
byproj = Dict{String,Vector{Int}}()
for (k, rs) in ens
    v = get!(byproj, rs[1].project, [0,0,0])
    v[1] += 1; v[2] += length(rs); v[3] += sum(x.nsteps for x in rs)
end
for (p, v) in sort(collect(byproj), by=x->-x[2][1])
    @printf("%-34s %6d ensembles %7d runs %10d samples\n", p, v[1], v[2], v[3])
end
