# Layout of the library tree:
#
#   <model>/<lattice_name>/<lattice_name>.toml
#   <model>/<lattice_name>/<parameters and sector>/beta_<beta>_<tag>.h5

_fmt(x::Real) = x isa Integer ? string(x) : string(Float64(x))

function _parameter_dir(e::Ensemble)
    parts = ["$k$(_fmt(v))" for (k, v) in sort(collect(e.parameters))]
    append!(parts, ["$k$(v)" for (k, v) in sort(collect(e.sector))])
    return isempty(parts) ? "default" : join(parts, "_")
end

"""
    relpath_for(e::Ensemble; tag) -> String

Location of an ensemble inside the library, relative to its root. `tag`
distinguishes files with otherwise identical metadata (chains, batches).
"""
function relpath_for(e::Ensemble; tag::AbstractString)
    isempty(tag) && throw(ArgumentError("tag must not be empty"))
    occursin(r"[/\\\s]", tag) && throw(ArgumentError("tag must not contain slashes or whitespace"))
    return joinpath(e.model, e.lattice_name, _parameter_dir(e), "beta_$(_fmt(e.beta))_$(tag).h5")
end

"Absolute path of the lattice file belonging to the HDF5 file at `h5path`."
lattice_path(h5path::AbstractString, lattice_name::AbstractString) =
    normpath(joinpath(dirname(abspath(h5path)), "..", lattice_name * ".toml"))

# Zenodo stores files flat; a path is flattened by replacing separators.
const _FLAT_SEP = "__"
flatten_path(rel::AbstractString) = replace(rel, "/" => _FLAT_SEP)
unflatten_path(name::AbstractString) = replace(name, _FLAT_SEP => "/")
