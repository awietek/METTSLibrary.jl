# Layout of the library tree:
#
#   <model>/<lattice_name>/<lattice_name>.toml
#   <model>/<lattice_name>/<parameters>/<sector>/T=<T>_beta=<beta>/<tag>.h5
#
# One file is one METTS run. The tag must identify the run uniquely: across the
# cluster archive the seed alone repeats at different bond dimensions, (seed,
# maxdim) repeats across source trees, and (seed, maxdim, source) repeats across
# X/Z comparison runs -- so the tag is seed<n>_maxm<m>_<basis>_<provenance>.

_fmt(x::Real) = x isa Integer ? string(x) : string(Float64(x))

# "J=0.4_t=3.0_t_prime=-0.3" from the couplings, "n=180" from the sector.
# The `=` matters: keys may contain `_` themselves and values may be negative,
# so gluing them together (`t_prime-0.3`) cannot be read back by eye.
_kv_dir(d) = join(["$k=$(_fmt(v))" for (k, v) in sort(collect(d))], "_")

_parameter_dir(e::Ensemble) = (s = _kv_dir(e.parameters); isempty(s) ? "default" : s)
_sector_dir(e::Ensemble)    = (s = _kv_dir(e.sector);     isempty(s) ? "default" : s)
# Both temperature and inverse temperature, because both are read: T is the
# value runs are specified with and is exactly representable, beta is what the
# file stores. Rounding matters -- 1/beta round-trips to 0.037500000000000006,
# so each label is normalised before it becomes part of a path.
_beta_dir(e::Ensemble) = string("T=", _fmt(round(1 / e.beta; sigdigits = 10)),
                               "_beta=", _fmt(round(e.beta; sigdigits = 6)))

"""
    relpath_for(e::Ensemble; tag) -> String

Location of an ensemble inside the library, relative to its root:
`<model>/<lattice_name>/<parameters>/<sector>/T=<T>_beta=<beta>/<tag>.h5`, e.g.
`tJ/square.L32.W4.cyl/J=0.4_t=3.0_t_prime=-0.3/ndn=56_nup=56/T=0.25_beta=4.0/seed3_maxm2000_X.h5`.
The temperature directory carries both labels so it reads either way: `T` is
exact, `beta` is rounded for legibility, and the file records both. `tag` must
identify the run within its ensemble — seed, bond dimension, collapse basis and
source tree, wherever those repeat.
"""
function relpath_for(e::Ensemble; tag::AbstractString)
    isempty(tag) && throw(ArgumentError("tag must not be empty"))
    occursin(r"[/\\\s]", tag) && throw(ArgumentError("tag must not contain slashes or whitespace"))
    return joinpath(e.model, e.lattice_name, _parameter_dir(e), _sector_dir(e),
                    _beta_dir(e), tag * ".h5")
end

"""
    lattice_path(h5path, lattice_name) -> String

Absolute path of the lattice file belonging to the HDF5 file at `h5path`. The
lattice lives at `<model>/<lattice_name>/<lattice_name>.toml`, so this walks up
from the file until it finds the directory named `lattice_name`. Walking rather
than counting levels keeps it correct if the depth below the lattice changes.
"""
function lattice_path(h5path::AbstractString, lattice_name::AbstractString)
    dir = dirname(abspath(h5path))
    while true
        basename(dir) == lattice_name && return joinpath(dir, lattice_name * ".toml")
        parent = dirname(dir)
        parent == dir && break
        dir = parent
    end
    error("no directory named '$lattice_name' above '$h5path'; the lattice file " *
          "belongs at <model>/<lattice_name>/<lattice_name>.toml")
end

# Zenodo stores files flat; a path is flattened by replacing separators.
const _FLAT_SEP = "__"
flatten_path(rel::AbstractString) = replace(rel, "/" => _FLAT_SEP)
unflatten_path(name::AbstractString) = replace(name, _FLAT_SEP => "/")
