# Shared by the converters.

# A lattice argument may be a path to a .toml/.lat file or TOML text itself.
_is_path(x::AbstractString) = !occursin('\n', x) && length(x) < 1024 && isfile(x)
_lattice_text(x::AbstractString) = _is_path(x) ? read_lattice(x) : (parse_lattice(x); String(x))

# lattice_name: explicit, or the lattice file's basename without extension
function _lattice_name(lattice::AbstractString, name)
    name === nothing || return String(name)
    _is_path(lattice) || throw(ArgumentError("lattice_name is required when the lattice is given as text"))
    return first(splitext(basename(lattice)))
end

_now() = Dates.format(now(UTC), dateformat"yyyy-mm-ddTHH:MM:SS\Z")

# Assemble and validate an Ensemble from converted states and user-supplied metadata.
function _converted(states::Matrix{UInt8}, steps::AbstractVector{<:Integer};
                    source::AbstractString, source_format::String, code::String,
                    lattice::AbstractString, lattice_name, model, site_type,
                    beta, temperature, basis, parameters, sector, algorithm, observables)
    (beta === nothing) == (temperature === nothing) &&
        throw(ArgumentError("give exactly one of beta or temperature"))
    M = length(steps)
    prov = Dict{String,Any}("created" => _now(), "source_format" => source_format,
                            "source_path" => abspath(source), "source_sha256" => sha256_file(source))
    _is_path(lattice) && (prov["source_lattice"] = abspath(lattice))
    return validate(Ensemble(
        model=String(model), site_type=String(site_type),
        lattice=_lattice_text(lattice), lattice_name=_lattice_name(lattice, lattice_name),
        beta=Float64(beta === nothing ? 1 / temperature : beta),
        parameters=Dict{String,Float64}(String(k) => Float64(v) for (k, v) in parameters),
        sector=Dict{String,Int}(String(k) => Int(v) for (k, v) in sector),
        algorithm=merge(Dict{String,Any}("code" => code), Dict{String,Any}(algorithm)),
        provenance=prov,
        collapse_bases=[uppercase(String(basis))],
        states=states, basis=zeros(UInt8, M), step=Int32.(steps),
        observables=Dict{String,Array{Float64}}(String(k) => v for (k, v) in observables),
    ))
end
