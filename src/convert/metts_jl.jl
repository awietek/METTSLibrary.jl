# ---------------------------------------------------------------------------
# Converters for METTS.jl output.
#
# The `examples/ttJ_metts.jl` driver writes collapsed states to `samples.txt`,
# one line per measurement:   `<meas>: [1, 2, 3, ...]`
# with 1-based ITensors indices for the "tJ" site type (1 Emp, 2 Up, 3 Dn).
# All stored samples are collapsed in the X basis; the run parameters are
# encoded in the output directory path.
# ---------------------------------------------------------------------------

"""
    from_samples_txt(path; lattice, model, site_type="tJ", beta=nothing, temperature=nothing,
                     basis="X", parameters=Dict(), sector=Dict(), algorithm=Dict(),
                     chain=1, lattice_name=nothing, observables=Dict()) -> Ensemble

Convert a METTS.jl `samples.txt`. Integers in the file are 1-based ITensors
indices for `site_type` and are stored 0-based. The measurement index becomes
`step`. The lattice (path or TOML text) and all metadata must be supplied;
`lattice_name` defaults to the lattice file's basename. See `from_ttj_run`
for the path-based convenience wrapper.
"""
function from_samples_txt(path::AbstractString;
                          lattice::AbstractString, model::AbstractString, site_type::AbstractString="tJ",
                          beta=nothing, temperature=nothing, basis::AbstractString="X",
                          parameters=Dict{String,Float64}(), sector=Dict{String,Int}(),
                          algorithm=Dict{String,Any}(), chain::Integer=1, lattice_name=nothing,
                          observables=Dict{String,Array{Float64}}())
    nlabels = length(LOCAL_STATES[site_type])
    samples = Dict{Int,Vector{Int}}()
    for line in eachline(path)
        parts = split(strip(line), ":"; limit=2)
        length(parts) == 2 || continue
        body = strip(parts[2], ['[', ']', ' '])
        isempty(body) || (samples[parse(Int, strip(parts[1]))] = parse.(Int, split(body, ",")))
    end
    isempty(samples) && error("no samples found in '$path'")
    steps = sort!(collect(keys(samples)))
    N = length(samples[steps[1]])
    all(length(samples[s]) == N for s in steps) || error("inconsistent sample lengths in '$path'")
    states = Matrix{UInt8}(undef, N, length(steps))
    for (j, s) in enumerate(steps), i in 1:N
        v = samples[s][i]
        1 <= v <= nlabels || error("state index $v out of range for site type '$site_type' (line $s)")
        states[i, j] = v - 1
    end
    return _converted(states, steps;
        source=path, source_format="metts_jl_samples_txt", code="METTS.jl",
        lattice, lattice_name, model, site_type, beta, temperature, basis,
        parameters, sector, algorithm, chain, observables)
end

const _TTJ_PATH_RE = r"L\.(\d+)\.W\.(\d+)/J\.([0-9.]+)/t\.([0-9.]+)/t_prime\.(-?[0-9.]+)/filling\.([0-9.]+)/T\.([0-9.]+)/D\.(\d+)/tau\.([0-9.]+)\.cutoff\.([0-9.eE+-]+)\.seed\.(\d+)"

"""
    parse_ttj_path(path) -> NamedTuple

Extract run parameters from an output path of `examples/ttJ_metts.jl`, which
has the form
`.../L.<L>.W.<W>/J.<J>/t.<t>/t_prime.<tp>/filling.<n>/T.<T>/D.<D>/tau.<tau>.cutoff.<c>.seed.<s>/...`.
"""
function parse_ttj_path(path::AbstractString)
    m = match(_TTJ_PATH_RE, path)
    m === nothing && error("path does not match the ttJ_metts.jl layout: '$path'")
    return (L=parse(Int, m[1]), W=parse(Int, m[2]), J=parse(Float64, m[3]), t=parse(Float64, m[4]),
            t_prime=parse(Float64, m[5]), filling=parse(Float64, m[6]), T=parse(Float64, m[7]),
            D=parse(Int, m[8]), tau=parse(Float64, m[9]), cutoff=parse(Float64, m[10]), seed=parse(Int, m[11]))
end

"""
    from_ttj_run(dir; lattice=nothing, lattice_name=nothing, basis="X", chain=nothing) -> Ensemble

Convert one output directory of `examples/ttJ_metts.jl` (containing
`samples.txt`) to an Ensemble, taking all parameters from the path. The
driver builds the t-t'-J Hamiltonian in code on `square_lattice(L, W;
yperiodic=true)`, so unless a `lattice` is given one is generated with
`square_lattice_toml` using couplings `t` (HOP, nearest neighbours), `J`
(HB, nearest neighbours) and `t_prime` (HOP, diagonals), named
`square.L<L>.W<W>.cyl`. The chain id defaults to the seed.
"""
function from_ttj_run(dir::AbstractString; lattice=nothing, lattice_name=nothing,
                      basis::AbstractString="X", chain=nothing)
    p = parse_ttj_path(abspath(dir))
    if lattice === nothing
        lattice = square_lattice_toml(p.L, p.W; yperiodic=true,
                                      bonds=[("t", "HOP", :nn), ("J", "HB", :nn), ("t_prime", "HOP", :nnn)])
        lattice_name = something(lattice_name, "square.L$(p.L).W$(p.W).cyl")
    end
    return from_samples_txt(joinpath(dir, "samples.txt");
        lattice, lattice_name, model="tJ", site_type="tJ", temperature=p.T, basis,
        parameters=Dict("t" => p.t, "t_prime" => p.t_prime, "J" => p.J),
        sector=Dict("n" => round(Int, p.filling * p.L * p.W)),
        algorithm=Dict{String,Any}("driver" => "examples/ttJ_metts.jl", "tau" => p.tau,
                                   "maxdim" => p.D, "cutoff" => p.cutoff, "seed" => p.seed),
        chain=something(chain, p.seed))
end
