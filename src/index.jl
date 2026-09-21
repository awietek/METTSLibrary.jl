const INDEX_FILE = "index.toml"

# The index is split per project, and each entry stores only what the path
# cannot give.
#
#   <root>/index.toml                     manifest: one line per project
#   <model>/<project>/index.toml          that project's lattices and ensembles
#
# Two reasons. Splitting keeps an ingest from rewriting the whole index: a
# single monolithic file was 1.8 kB per entry, so the full archive of 35,218
# runs would have been a ~63 MB text file rewritten in git on every batch.
# Slimming removes what the layout already carries -- model, project,
# lattice_name, parameters and sector all come back from the path, and
# nsites, couplings and the lattice's own path and hash belong to the lattice,
# not to each of the hundreds of runs that share it.
#
# What an entry keeps is what the path cannot express: the file's hash and
# size, its sample count, its observables, its collapse bases and site type,
# its algorithm block, and `beta`. Beta is stored even though the temperature
# directory names it, because that directory rounds to six decimals for
# legibility and 1/T is non-terminating for most temperatures on the cluster.
#
# `load_index` expands entries back into the full shape, so `ensembles` and
# `load` see exactly what they saw when the index was one flat file.

# name=value pairs out of a directory or tag component. Keys may contain "_"
# and values may be negative or non-numeric, so pairs are split on the "_"
# that precedes the next "<name>=" rather than on every "_".
const _KV_RE = r"([A-Za-z][A-Za-z0-9_]*)=(.*?)(?=_[A-Za-z][A-Za-z0-9_]*=|$)"

function _kv_parse(s::AbstractString)
    s == "default" && return Pair{String,String}[]
    return [String(m.captures[1]) => String(m.captures[2]) for m in eachmatch(_KV_RE, s)]
end

"""
    parse_relpath(rel) -> NamedTuple

Recover from a library-relative ensemble path everything the layout encodes:
`model`, `project`, `lattice_name`, `parameters`, `sector` and `tag`.
"""
function parse_relpath(rel::AbstractString)
    p = splitpath(rel)
    length(p) == 7 || error("'$rel' is not <model>/<project>/<lattice>/<parameters>/" *
                            "<sector>/T=<T>_beta=<beta>/<tag>.h5")
    return (model = p[1], project = p[2], lattice_name = p[3],
            parameters = Dict{String,Float64}(k => parse(Float64, v) for (k, v) in _kv_parse(p[4])),
            sector     = Dict{String,Int}(k => parse(Int, v) for (k, v) in _kv_parse(p[5])),
            tag        = replace(p[7], r"\.h5$" => ""))
end

_scalars(d) = Dict{String,Any}(k => v for (k, v) in d if v isa Union{Real,AbstractString,Bool})

# "1000" -> 1000, "0.2" -> 0.2, "X" -> "X"
function _untag(s::AbstractString)
    v = tryparse(Int, s);     v === nothing || return v
    w = tryparse(Float64, s); w === nothing || return w
    return String(s)
end

"The directory part below a project: <lattice>/<parameters>/<sector>/T=…_beta=…"
function _parse_dir(dir::AbstractString)
    p = splitpath(dir)
    length(p) == 4 || error("'$dir' is not <lattice>/<parameters>/<sector>/T=<T>_beta=<beta>")
    return (lattice_name = p[1],
            parameters = Dict{String,Float64}(k => parse(Float64, v) for (k, v) in _kv_parse(p[2])),
            sector     = Dict{String,Int}(k => parse(Int, v) for (k, v) in _kv_parse(p[3])))
end

# Expand one run of one ensemble back into the flat entry callers expect.
function _expand_run(model, project, ens, run, defaults, lattices, tag_algorithm)
    q   = _parse_dir(ens["dir"])
    lat = get(lattices, q.lattice_name, nothing)
    alg = Dict{String,Any}(get(defaults, "algorithm", Dict{String,Any}()))
    merge!(alg, get(run, "algorithm", Dict{String,Any}()))
    # the fields the filename already carries come back from it
    tagkv = Dict(_kv_parse(run["tag"]))
    for k in tag_algorithm
        haskey(tagkv, k) && (alg[k] = _untag(tagkv[k]))
    end

    out = Dict{String,Any}(
        "path"        => joinpath(model, project, ens["dir"], run["tag"] * ".h5"),
        "model"       => model,
        "project"     => project,
        "lattice_name"=> q.lattice_name,
        "parameters"  => q.parameters,
        "sector"      => q.sector,
        "beta"        => ens["beta"],
        "temperature" => 1 / ens["beta"],
        "nsamples"    => run["nsamples"],
        "bytes"       => run["bytes"],
        "sha256"      => run["sha256"],
        "algorithm"   => alg,
    )
    for k in ("site_type", "collapse_bases", "observables")
        v = get(run, k, get(defaults, k, nothing))
        v === nothing || (out[k] = v)
    end
    if lat !== nothing
        out["nsites"]         = lat["nsites"]
        out["couplings"]      = lat["couplings"]
        out["lattice"]        = lat["path"]
        out["lattice_sha256"] = lat["sha256"]
    end
    return out
end

function _expand_project(pidx)
    model, project = pidx["model"], pidx["project"]
    defaults  = get(pidx, "defaults", Dict{String,Any}())
    lattices  = Dict{String,Any}(l["name"] => l for l in get(pidx, "lattices", []))
    tagalg    = get(pidx, "tag_algorithm", String[])
    out = Dict{String,Any}[]
    for ens in get(pidx, "ensembles", []), run in ens["runs"]
        push!(out, _expand_run(model, project, ens, run, defaults, lattices, tagalg))
    end
    return out
end

_stamp() = Dates.format(now(UTC), dateformat"yyyy-mm-ddTHH:MM:SS\Z")

# Assemble one project's index: constants hoisted into `defaults`, runs grouped
# under the ensemble directory they share, and anything the filename already
# encodes left out.
function _build_project(root, model, project, rels)
    lattices = Dict{String,Any}()
    raw = Tuple{String,String,Dict{String,Any}}[]     # dir, tag, facts
    for rel in rels
        path = joinpath(root, rel)
        is_lfs_pointer(path) &&
            error("'$rel' is a Git LFS pointer, not content. Run `git lfs pull` before building the index.")
        m = read_metadata(path)
        name = m["lattice_name"]
        if !haskey(lattices, name)
            lp = lattice_path(path)
            lattices[name] = Dict{String,Any}(
                "name" => name, "path" => relpath(lp, root), "sha256" => sha256_file(lp),
                "nsites" => m["nsites"], "couplings" => m["couplings"])
        end
        p = splitpath(rel)
        push!(raw, (joinpath(p[3:6]...), replace(p[7], r"\.h5$" => ""),
                    Dict{String,Any}(
                        "sha256" => sha256_file(path), "bytes" => filesize(path),
                        "beta" => m["beta"], "nsamples" => m["nsamples"],
                        "site_type" => m["site_type"], "collapse_bases" => m["collapse_bases"],
                        "observables" => m["observables"],
                        "algorithm" => _scalars(m["algorithm"]))))
    end

    # Which algorithm entries does the filename already carry? Drop exactly
    # those, and only where the tag really spells the same value.
    tagalg = Set{String}()
    for (_, tag, f) in raw
        tagkv = Dict(_kv_parse(tag))
        for (k, v) in f["algorithm"]
            haskey(tagkv, k) && _untag(tagkv[k]) == v && push!(tagalg, k)
        end
    end
    for (_, _, f) in raw
        f["algorithm"] = Dict{String,Any}(k => v for (k, v) in f["algorithm"] if !(k in tagalg))
    end

    # Fields with one value across the whole project belong in the header, not
    # on every run. site_type, collapse_bases and observables are typically
    # constant, and so are most of the algorithm's remaining entries.
    defaults = Dict{String,Any}()
    for k in ("site_type", "collapse_bases", "observables")
        vs = unique(f[k] for (_, _, f) in raw)
        length(vs) == 1 && (defaults[k] = only(vs); foreach(t -> delete!(t[3], k), raw))
    end
    algkeys = union((keys(f["algorithm"]) for (_, _, f) in raw)...)
    dalg = Dict{String,Any}()
    for k in algkeys
        vs = unique(get(f["algorithm"], k, nothing) for (_, _, f) in raw)
        length(vs) == 1 && only(vs) !== nothing || continue
        dalg[k] = only(vs)
        foreach(t -> delete!(t[3]["algorithm"], k), raw)
    end
    isempty(dalg) || (defaults["algorithm"] = dalg)

    # Group the runs under the ensemble directory they share, so a ~110
    # character path prefix is written once instead of once per chain.
    order = String[]; byens = Dict{String,Vector{Any}}()
    betas = Dict{String,Float64}()
    for (dir, tag, f) in raw
        haskey(byens, dir) || (push!(order, dir); byens[dir] = [])
        betas[dir] = f["beta"]
        run = Dict{String,Any}("tag" => tag, "nsamples" => f["nsamples"],
                               "bytes" => f["bytes"], "sha256" => f["sha256"])
        for k in ("site_type", "collapse_bases", "observables")
            haskey(f, k) && (run[k] = f[k])          # only where it differs
        end
        isempty(f["algorithm"]) || (run["algorithm"] = f["algorithm"])
        push!(byens[dir], run)
    end

    return Dict{String,Any}(
        "schema_version" => SCHEMA_VERSION,
        "generated"      => _stamp(),
        "model"          => model,
        "project"        => project,
        "tag_algorithm"  => sort(collect(tagalg)),
        "defaults"       => defaults,
        "lattices"       => [lattices[k] for k in sort(collect(keys(lattices)))],
        "ensembles"      => [Dict{String,Any}("dir" => d, "beta" => betas[d],
                                              "runs" => byens[d]) for d in order],
    )
end

function _toml_string(d)
    io = IOBuffer()
    TOML.print(io, d; sorted=true)
    return String(take!(io))
end

"""
    build_index(root; write=true) -> Dict

Walk `root` for `*.h5` ensemble files and assemble the index. With
`write=true` also writes one `index.toml` per project plus the root manifest.

Returns the expanded, flat view — every entry as if the index were still a
single file — so callers do not have to know the index is split.
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

    groups = Dict{Tuple{String,String},Vector{String}}()
    for rel in rels
        q = parse_relpath(rel)
        push!(get!(groups, (q.model, q.project), String[]), rel)
    end

    projects = Dict{String,Any}[]
    flat     = Dict{String,Any}[]
    for (model, project) in sort(collect(keys(groups)))
        pidx = _build_project(root, model, project, groups[(model, project)])
        rel_idx = joinpath(model, project, INDEX_FILE)
        text    = _toml_string(pidx)
        write && Base.write(joinpath(root, rel_idx), text)
        nfiles = sum(length(e["runs"]) for e in pidx["ensembles"]; init=0)
        push!(projects, Dict{String,Any}(
            "model"     => model,
            "project"   => project,
            "index"     => rel_idx,
            "sha256"    => sha256_string(text),
            "nfiles"    => nfiles,
            "nensembles"=> length(pidx["ensembles"]),
            "nsamples"  => sum(r["nsamples"] for e in pidx["ensembles"] for r in e["runs"]; init=0),
            "lattices"  => [l["name"] for l in pidx["lattices"]],
        ))
        append!(flat, _expand_project(pidx))
    end

    manifest = Dict{String,Any}(
        "schema_version" => SCHEMA_VERSION,
        "generated"      => _stamp(),
        "projects"       => projects,
    )
    write && Base.write(joinpath(root, INDEX_FILE), _toml_string(manifest))
    return Dict{String,Any}("schema_version" => SCHEMA_VERSION,
                            "generated" => manifest["generated"],
                            "projects" => projects,
                            "ensembles" => flat)
end

# Read one project index, locally if we have a root, otherwise from a remote.
function _read_project_index(root, rel, sha)
    root !== nothing && isfile(joinpath(root, rel)) && return TOML.parsefile(joinpath(root, rel))
    return TOML.parsefile(fetch_file(rel; sha256=sha))
end

"""
    load_index(; source=nothing) -> Dict

Load the index and expand it into the flat view: `idx["ensembles"]` is a
vector of entries carrying everything, including the fields recovered from
each path. `source` may be a library root, a path to the root manifest, or
`nothing` — in which case `METTSLIBRARY_PATH` is used if set, and otherwise
the manifest and the project indexes are fetched from the configured remotes.

Only the projects actually present are read, so a remote query costs one small
manifest plus one file per project rather than one large file for everything.
"""
function load_index(; source=nothing)
    local root, manifest
    if source === nothing
        root = data_root()
        manifest = root !== nothing && isfile(joinpath(root, INDEX_FILE)) ?
            TOML.parsefile(joinpath(root, INDEX_FILE)) :
            TOML.parsefile(fetch_file(INDEX_FILE; refresh=true))
    elseif isdir(source)
        root = source
        manifest = TOML.parsefile(joinpath(source, INDEX_FILE))
    else
        root = dirname(abspath(source))
        manifest = TOML.parsefile(source)
    end

    flat = Dict{String,Any}[]
    for p in get(manifest, "projects", Dict{String,Any}[])
        append!(flat, _expand_project(_read_project_index(root, p["index"], p["sha256"])))
    end
    out = Dict{String,Any}(manifest)
    out["ensembles"] = flat
    return out
end

_matches(a::Real, b::Real) = isapprox(a, b; rtol=1e-8, atol=1e-12)
_matches(a, b) = a == b

function _entry_matches(entry, key::String, value)
    if haskey(entry, key)
        return _matches(entry[key], value)
    end
    for grp in ("parameters", "sector", "algorithm")
        g = get(entry, grp, nothing)
        g !== nothing && haskey(g, key) && return _matches(g[key], value)
    end
    return false
end

"""
    ensembles(idx=load_index(); kwargs...) -> Vector{Dict}

Filter index entries. Keywords match top-level fields (`model`, `project`,
`site_type`, `lattice_name`, `nsites`, `beta`, `temperature`, `nsamples`, ...)
or entries of `parameters`, `sector` and `algorithm`. Numbers are compared
approximately.

    ensembles(model="tJ", lattice_name="square.L32.W4.cyl", beta=4.0, t=3.0)

`provenance` is not indexed: it is per-run bookkeeping that nobody queries and
it was two thirds of the index's size. Read it from the file with
`read_ensemble` when you need it.
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
