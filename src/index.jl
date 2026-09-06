const INDEX_FILE = "index.toml"

# Turn metadata + file facts into one TOML-friendly index entry.
function _index_entry(root::AbstractString, rel::AbstractString)
    path = joinpath(root, rel)
    is_lfs_pointer(path) &&
        error("'$rel' is a Git LFS pointer, not content. Run `git lfs pull` before building the index.")
    m = read_metadata(path)
    entry = Dict{String,Any}(
        "path"   => rel,
        "sha256" => sha256_file(path),
        "bytes"  => filesize(path),
    )
    for k in ("model", "site_type", "nsites", "nsamples", "beta", "collapse_bases",
              "observables", "couplings", "lattice_name", "lattice_sha256")
        entry[k] = m[k]
    end
    # lattice file, relative to the library root
    entry["lattice"] = relpath(lattice_path(path, m["lattice_name"]), root)
    startswith(entry["lattice"], "..") &&
        error("lattice file of '$rel' lies outside the library root: $(entry["lattice"])")
    entry["parameters"] = m["parameters"]
    entry["sector"]     = m["sector"]
    # Only scalar algorithm/provenance fields go to the index; TOML cannot hold nothing.
    for grp in ("algorithm", "provenance")
        entry[grp] = Dict{String,Any}(k => v for (k, v) in m[grp] if v isa Union{Real,AbstractString,Bool})
    end
    return entry
end

"""
    build_index(root; write=true) -> Dict

Walk `root` for `*.h5` ensemble files, read their metadata and hashes, and
assemble the index. With `write=true` also writes `<root>/index.toml`.
Hidden directories (e.g. `.git`) are skipped.
"""
function build_index(root::AbstractString; write::Bool=true)
    root = abspath(root)
    rels = String[]
    for (dir, subdirs, files) in walkdir(root)
        filter!(d -> !startswith(d, "."), subdirs)
        for f in files
            endswith(f, ".h5") || continue
            push!(rels, relpath(joinpath(dir, f), root))
        end
    end
    sort!(rels)
    entries = [_index_entry(root, rel) for rel in rels]
    idx = Dict{String,Any}(
        "schema_version" => SCHEMA_VERSION,
        "generated"      => Dates.format(now(UTC), dateformat"yyyy-mm-ddTHH:MM:SS\Z"),
        "ensembles"      => entries,
    )
    if write
        open(joinpath(root, INDEX_FILE), "w") do io
            TOML.print(io, idx; sorted=true)
        end
    end
    return idx
end

"""
    load_index(; source=nothing) -> Dict

Load the index. `source` may be a path to an `index.toml`, a directory
containing one, or `nothing`, in which case the local data root
(`METTSLIBRARY_PATH`) is used if set, and otherwise the index is fetched from
the configured remotes.
"""
function load_index(; source=nothing)
    if source === nothing
        root = data_root()
        if root !== nothing && isfile(joinpath(root, INDEX_FILE))
            return TOML.parsefile(joinpath(root, INDEX_FILE))
        end
        return TOML.parsefile(fetch_file(INDEX_FILE; refresh=true))
    elseif isdir(source)
        return TOML.parsefile(joinpath(source, INDEX_FILE))
    else
        return TOML.parsefile(source)
    end
end

_matches(a::Real, b::Real) = isapprox(a, b; rtol=1e-8, atol=1e-12)
_matches(a, b) = a == b

function _entry_matches(entry, key::String, value)
    if haskey(entry, key)
        return _matches(entry[key], value)
    end
    for grp in ("parameters", "sector", "algorithm", "provenance")
        g = get(entry, grp, nothing)
        g !== nothing && haskey(g, key) && return _matches(g[key], value)
    end
    return false
end

"""
    ensembles(idx=load_index(); kwargs...) -> Vector{Dict}

Filter index entries. Keywords match top-level fields (`model`, `site_type`,
`lattice_name`, `nsites`, `beta`, `nsamples`, ...) or entries of `parameters`,
`sector`, `algorithm`, `provenance`. Numbers are compared approximately.

    ensembles(model="tJ", lattice_name="square.L32.W4.cyl", beta=4.0, t=3.0)
"""
function ensembles(idx::AbstractDict; kwargs...)
    out = Vector{Dict{String,Any}}()
    for entry in idx["ensembles"]
        all(_entry_matches(entry, String(k), v) for (k, v) in kwargs) && push!(out, entry)
    end
    return out
end
ensembles(; kwargs...) = ensembles(load_index(); kwargs...)

"""
    load(entry::AbstractDict) -> Ensemble

Fetch (if needed), verify and read the ensemble described by an index entry,
together with its lattice file.
"""
function load(entry::AbstractDict)
    fetch_file(entry["lattice"]; sha256=entry["lattice_sha256"])
    return read_ensemble(fetch_file(entry["path"]; sha256=entry["sha256"]))
end
