# ---------------------------------------------------------------------------
# HDF5 layout (schema version 2); see docs/src/schema.md
#
# Version 2 dropped /chain: one file is one METTS run, so it was a constant
# column duplicating algorithm/seed. Files written by version 1 still read —
# their /chain dataset is simply ignored.
#
# /                          attributes: schema, schema_version, model, site_type,
#                            local_states, nsites, nsamples, beta, collapse_bases,
#                            lattice_name, lattice_sha256
# ../<lattice_name>.toml     the lattice file, one directory up, shared
# /coordinates               Float64 (dim, nsites), copy of the lattice's Coordinates
# /states                    UInt8   (nsites, nsamples)   0-based into local_states
# /basis                     UInt8   (nsamples,)          0-based into collapse_bases
# /step                      Int32   (nsamples,)
# /parameters                group, one scalar attribute per coupling
# /sector                    group, one integer attribute per quantum number
# /algorithm                 group, scalar attributes
# /provenance                group, scalar attributes
# /observables/<name>        Float64 arrays, last dimension nsamples
# ---------------------------------------------------------------------------

const _ATTR_SCALAR = Union{Integer,AbstractFloat,AbstractString,Bool}

_attr_value(v::_ATTR_SCALAR) = v
_attr_value(v::AbstractVector{<:AbstractString}) = String.(v)
_attr_value(v::AbstractVector{<:Real}) = collect(v)
_attr_value(v) = string(v)

function _write_attr_group!(parent, name::String, d::AbstractDict)
    g = create_group(parent, name)
    for (k, v) in d
        attributes(g)[String(k)] = _attr_value(v)
    end
    return g
end

function _read_attr_group(parent, name::String)
    d = Dict{String,Any}()
    haskey(parent, name) || return d
    g = parent[name]
    for k in keys(attributes(g))
        d[k] = read_attribute(g, k)
    end
    return d
end

sha256_string(s::AbstractString) = bytes2hex(sha256(codeunits(s)))

"""
    write_ensemble(root, e::Ensemble; tag) -> relpath

Write `e` into the library at `root`, at
`<model>/<lattice_name>/<parameters>/beta_<beta>_<tag>.h5`, and return that
relative path. The lattice is written to `<model>/<lattice_name>/<lattice_name>.toml`
if it is not there yet; if it is, it must be identical, since other
ensembles share it. Files are append-only: an existing file at the target
path is an error, use another `tag`. Validates before writing.
"""
function write_ensemble(root::AbstractString, e::Ensemble; tag::AbstractString)
    validate(e)
    rel  = relpath_for(e; tag)
    path = joinpath(root, rel)
    isfile(path) && error("'$rel' already exists. Ensemble files are append-only; use a different tag.")
    mkpath(dirname(path))

    lp = lattice_path(path, e.lattice_name)
    lhash = sha256_string(e.lattice)
    if isfile(lp)
        sha256_file(lp) == lhash ||
            error("lattice file '$lp' exists with different content. It is shared by the ensembles " *
                  "below it and is never overwritten; this is a different lattice and needs a different name.")
    else
        write(lp, e.lattice)
    end

    N, M = size(e.states)
    h5open(path, "w") do f
        a = attributes(f)
        a["schema"]         = SCHEMA_NAME
        a["schema_version"] = SCHEMA_VERSION
        a["model"]          = e.model
        a["site_type"]      = e.site_type
        a["local_states"]   = e.local_states
        a["nsites"]         = N
        a["nsamples"]       = M
        a["beta"]           = e.beta
        a["collapse_bases"] = e.collapse_bases
        a["lattice_name"]   = e.lattice_name
        a["lattice_sha256"] = lhash

        f["coordinates"] = parse_lattice(e.lattice).coordinates
        ds = create_dataset(f, "states", datatype(UInt8), dataspace(e.states);
                            chunk=(N, max(1, min(M, 4096))), deflate=3)
        write(ds, e.states)
        f["basis"] = e.basis
        f["step"]  = e.step

        _write_attr_group!(f, "parameters", e.parameters)
        _write_attr_group!(f, "sector",     e.sector)
        _write_attr_group!(f, "algorithm",  e.algorithm)
        _write_attr_group!(f, "provenance", e.provenance)

        obs = create_group(f, "observables")
        for (k, v) in e.observables
            obs[k] = v
        end
    end
    return rel
end

function _check_schema(f, path)
    a = attributes(f)
    haskey(a, "schema") && read_attribute(f, "schema") == SCHEMA_NAME ||
        error("'$path' is not a $SCHEMA_NAME file (missing or wrong 'schema' attribute).")
    v = read_attribute(f, "schema_version")
    v <= SCHEMA_VERSION ||
        error("'$path' has schema version $v, this package understands up to $SCHEMA_VERSION.")
    return v
end

# Everything except the sample arrays, including the verified lattice text.
function _read_header(f, path)
    _check_schema(f, path)
    d = Dict{String,Any}(k => read_attribute(f, k) for k in
        ("model", "site_type", "local_states", "nsites", "nsamples", "beta", "collapse_bases",
         "lattice_name", "lattice_sha256"))
    d["parameters"]  = Dict{String,Float64}(k => Float64(v) for (k, v) in _read_attr_group(f, "parameters"))
    d["sector"]      = Dict{String,Int}(k => Int(v) for (k, v) in _read_attr_group(f, "sector"))
    d["algorithm"]   = _read_attr_group(f, "algorithm")
    d["provenance"]  = _read_attr_group(f, "provenance")
    d["observables"] = haskey(f, "observables") ? sort!(collect(keys(f["observables"]))) : String[]

    lp = lattice_path(path, d["lattice_name"])
    isfile(lp) || error("lattice file '$lp' belonging to '$path' not found")
    d["lattice"] = read(lp, String)
    sha256_string(d["lattice"]) == d["lattice_sha256"] ||
        error("lattice file '$lp' has been modified: its hash does not match the one recorded in '$path'")
    d["couplings"] = lattice_couplings(d["lattice"])
    return d
end

# Metadata only; used to build the index.
read_metadata(path::AbstractString) = h5open(f -> _read_header(f, path), path, "r")

"""
    read_ensemble(path) -> Ensemble

Read an ensemble file together with its lattice file, verify both, and
validate the result.
"""
function read_ensemble(path::AbstractString)
    h5open(path, "r") do f
        h = _read_header(f, path)
        obs = Dict{String,Array{Float64}}(k => Array{Float64}(read(f["observables"][k])) for k in h["observables"])
        validate(Ensemble(
            model=h["model"], site_type=h["site_type"], local_states=Vector{String}(h["local_states"]),
            lattice=h["lattice"], lattice_name=h["lattice_name"], beta=Float64(h["beta"]),
            parameters=h["parameters"], sector=h["sector"], algorithm=h["algorithm"], provenance=h["provenance"],
            collapse_bases=Vector{String}(h["collapse_bases"]),
            states=Matrix{UInt8}(read(f["states"])), basis=Vector{UInt8}(read(f["basis"])),
            step=Vector{Int32}(read(f["step"])),
            observables=obs))
    end
end

# Hex SHA-256 of a file, streamed. Coincides with the Git LFS OID.
function sha256_file(path::AbstractString)
    ctx = SHA.SHA256_CTX()
    open(path, "r") do io
        buf = Vector{UInt8}(undef, 1 << 20)
        while !eof(io)
            n = readbytes!(io, buf)
            SHA.update!(ctx, view(buf, 1:n))
        end
    end
    return bytes2hex(SHA.digest!(ctx))
end

# True if the file is a Git LFS pointer rather than content (GIT_LFS_SKIP_SMUDGE clones).
function is_lfs_pointer(path::AbstractString)
    filesize(path) < 512 || return false
    return startswith(read(path, String), "version https://git-lfs")
end
