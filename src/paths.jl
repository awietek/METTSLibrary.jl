# Layout of the library tree:
#
#   <model>/<project>/README.md                  who ran it, when, which papers
#   <model>/<project>/<lattice_name>/<lattice_name>.toml
#   <model>/<project>/<lattice_name>/<parameters>/<sector>/T=<T>_beta=<beta>/<tag>.h5
#
# A project is one researcher's body of work and may span several lattices.
# Runs at identical physics in two projects are different data and stay apart;
# finding all data at some couplings is an index query, not a directory listing.
#
# One file is one METTS run, and the tag is its filename: the last thing left
# to tell two runs apart once the path has fixed the physics. So the tag is
# method, not physics -- the algorithm parameters that were varied. See
# TAG_FIELDS and default_tag below.

# name=value path components, and their inverse. "J=0.4_t=3.0_t_prime=-0.3"
# from the couplings, "n=180" from the sector, and the same spelling for the
# tag. The `=` matters: keys may contain `_` themselves and values may be
# negative, so gluing them together (`t_prime-0.3`) could be read back neither
# by eye nor by _kv_parse.
#
# The writer and the reader live here side by side on purpose: the path is the
# authoritative record of a file's model, project, lattice, couplings and
# sector, so these two must stay exact inverses.
_fmt(x::Integer) = string(x)
_fmt(x::Real)    = string(Float64(x))
_fmt(x)          = string(x)

_kv_dir(d) = join(["$k=$(_fmt(v))" for (k, v) in sort(collect(d))], "_")
_kv_dir_or(d, default::AbstractString) = (s = _kv_dir(d); isempty(s) ? default : s)

# Split on the "_" that precedes the next "<name>=", never on every "_", so a
# key like t_prime survives and a value may be negative or non-numeric.
const _KV_RE = r"([A-Za-z][A-Za-z0-9_]*)=(.*?)(?=_[A-Za-z][A-Za-z0-9_]*=|$)"

# "canonical" is the empty SECTOR: no quantum number is fixed, so the ensemble
# is the canonical one at its temperature. "default" is the empty parameter
# set -- a different thing, and essentially never seen, since every coupling
# the lattice declares must have a value. Both read back as empty.
_kv_parse(s::AbstractString) = (s == "canonical" || s == "default") ?
    Pair{String,String}[] :
    [String(m.captures[1]) => String(m.captures[2]) for m in eachmatch(_KV_RE, s)]
# Both temperature and inverse temperature, because both get read: T is the
# value runs are specified with, beta is what the file stores.
#
# Fixed 5+6 field. Zero padding is not decoration: lexicographic and numeric
# order agree only when every number has the same count of integer digits, so
# the padding is what makes `ls` list a temperature sweep in temperature order.
# Six decimals keep the value exact -- and round away the float round-trip,
# since 1/(1/0.0375) is 0.037500000000000006.
_tfield(x::Real) = @sprintf("%012.6f", x)
_beta_dir(e::Ensemble) = string("T=", _tfield(1 / e.beta), "_beta=", _tfield(e.beta))

"""
Algorithm parameters the default tag is built from, in order, plus the
pseudo-field `basis` for the collapse basis (which lives in `collapse_bases`,
not in `algorithm`). A field missing from an ensemble's `algorithm` is skipped,
so a run that recorded less simply gets a shorter name.

This is a default, not a fixed vocabulary. A different time evolution records
different parameters — an expansion order, a Krylov dimension, a number of
sweeps — and naming those is a matter of passing `fields` to `default_tag`, or
of amending this vector once for a whole ingest. What a tag has to achieve is
uniqueness within one ensemble directory; `write_ensemble` refuses a file that
would collide, so a field set that omits what actually varies announces itself
at the first duplicate rather than silently overwriting.
"""
const TAG_FIELDS = ["basis", "maxdim", "tau", "cutoff", "seed"]

"""
    default_tag(e::Ensemble; fields=TAG_FIELDS) -> String

Filename stem for an ensemble, as `name=value` pairs joined by `_` — the same
spelling the parameter and sector directories use, e.g.
`basis=X_maxdim=1000_tau=0.1_cutoff=1.0e-10_seed=3`. Values come from
`e.algorithm`, except `basis`, which is `e.collapse_bases`. Fields absent from
`algorithm` are skipped; pass `fields` to name others.

Nothing is lost by leaving a field out: every value is stored properly inside
the file and is indexed. The tag exists to be read by eye and to be unique.
"""
function default_tag(e::Ensemble; fields=TAG_FIELDS)
    parts = String[]
    for f in fields
        k = String(f)
        v = k == "basis" ? join(e.collapse_bases) : get(e.algorithm, k, nothing)
        v === nothing && continue
        push!(parts, string(k, "=", _fmt(v)))
    end
    isempty(parts) && throw(ArgumentError(
        "cannot build a tag: the ensemble records none of $(join(fields, ", ")). " *
        "Pass fields=[...] naming the parameters that distinguish this run, or tag=\"...\"."))
    return join(parts, "_")
end

"""
    relpath_for(e::Ensemble; tag=default_tag(e)) -> String

Location of an ensemble inside the library, relative to its root:
`<model>/<project>/<lattice_name>/<parameters>/<sector>/T=<T>_beta=<beta>/<tag>.h5`, e.g.
`tJ/superconductors/square.L32.W4.cyl/J=0.4_t=3.0_t_prime=-0.3/ndn=56_nup=56/T=00000.250000_beta=00004.000000/basis=X_maxdim=2000_tau=0.1_seed=3.h5`.
The temperature directory carries both labels so it reads either way, each in a
fixed 5+6 field: zero padded so a directory listing comes out in temperature
order, six decimals so the value is exact. The file stores both as attributes
and both are indexed. `tag` must identify the run within its ensemble; by
default it is `default_tag(e)`.
"""
function relpath_for(e::Ensemble; tag::AbstractString=default_tag(e))
    # All four become single path components, and model, project and
    # lattice_name are read back out of the path, so a slash in any of them
    # would silently change what the file claims to be.
    for (what, s) in ("model" => e.model, "project" => e.project,
                      "lattice_name" => e.lattice_name, "tag" => tag)
        (isempty(s) || occursin(r"[/\\\s]", s)) && throw(ArgumentError(
            "$what must be non-empty and contain no slashes or whitespace, got '$s'"))
    end
    return joinpath(e.model, e.project, e.lattice_name,
                    _kv_dir_or(e.parameters, "default"), _kv_dir_or(e.sector, "canonical"),
                    _beta_dir(e), tag * ".h5")
end

"""
    lattice_dir(h5path) -> String

Directory holding the lattice of the ensemble file at `h5path`. `relpath_for`
puts exactly `<parameters>/<sector>/T=<T>_beta=<beta>` between the lattice
directory and the file, so the lattice directory is four levels up. This is
positional on purpose: the lattice is found by *where it is*, never by a name
stored inside the file, so renaming a lattice is a `mv` and not a rewrite of
every file that uses it.
"""
lattice_dir(h5path::AbstractString) = dirname(dirname(dirname(dirname(abspath(h5path)))))

"""
    lattice_path(h5path) -> String

The lattice file of the ensemble at `h5path`: the single `.toml` in its
`lattice_dir`. A lattice directory holds exactly one lattice, so there is
nothing to disambiguate and no name to match.

The file is not checksummed against the ensemble. A lattice that does not
belong to the data fails `validate` instead — on site count, on couplings the
parameters do not cover, or on failing to parse — which catches the mistakes
that actually happen. Only an edit preserving site count and coupling names
would slip through, and the lattice is append-only like everything else.
"""
function lattice_path(h5path::AbstractString)
    dir = lattice_dir(h5path)
    isdir(dir) || error("no lattice directory above '$h5path' (looked at '$dir')")
    tomls = filter(f -> endswith(f, ".toml"), readdir(dir))
    length(tomls) == 1 ||
        error("expected exactly one .toml in the lattice directory '$dir', found " *
              (isempty(tomls) ? "none" : join(tomls, ", ")))
    return joinpath(dir, tomls[1])
end

"""
    layout_names(h5path) -> (model, project, lattice_name)

The three names `relpath_for` encodes in the path, read back from it:
`<model>/<project>/<lattice_name>/…`. All three are authoritative — none is
stored in the ensemble file.

They are names, not data: they say where a file belongs, and nothing in them
is needed to interpret its contents (that is `site_type`, `local_states`, the
lattice, `parameters`, `sector` and `beta`, all of which *are* stored). Keeping
them out of the file means renaming a model, a project or a lattice is a `mv`
plus `build_index`, and that no stored copy can drift from the directory.
"""
function layout_names(h5path::AbstractString)
    lat  = lattice_dir(h5path)
    proj = dirname(lat)
    return (basename(dirname(proj)), basename(proj), basename(lat))
end
