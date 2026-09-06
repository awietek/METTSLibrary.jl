# ---------------------------------------------------------------------------
# Zenodo: reading a record's file list (for ZenodoRemote) and publishing a
# set of library files as a new record (for sharing frozen subsets).
# Uses the REST API; token from the `ZENODO_TOKEN` environment variable
# (`ZENODO_SANDBOX_TOKEN` for the sandbox) unless passed explicitly.
# ---------------------------------------------------------------------------

zenodo_api(sandbox::Bool) = sandbox ? "https://sandbox.zenodo.org/api" : "https://zenodo.org/api"

function zenodo_token(sandbox::Bool)
    t = get(ENV, sandbox ? "ZENODO_SANDBOX_TOKEN" : "ZENODO_TOKEN", "")
    return isempty(t) ? nothing : strip(t)
end

function _zenodo_request(method::String, url::String; token=nothing, body=nothing, input=nothing,
                         content_type="application/json")
    headers = Pair{String,String}[]
    token === nothing || push!(headers, "Authorization" => "Bearer $token")
    body === nothing || push!(headers, "Content-Type" => content_type)
    out = IOBuffer()
    inp = body !== nothing ? IOBuffer(body) : input
    resp = Downloads.request(url; method, headers, input=inp, output=out, throw=false)
    resp isa Downloads.Response ||
        error("request to $url failed: " * sprint(showerror, resp))
    text = String(take!(out))
    200 <= resp.status < 300 ||
        error("Zenodo returned HTTP $(resp.status) for $method $url\n$text")
    return isempty(text) ? nothing : JSON.parse(text)
end

"""
    zenodo_files(record; token=nothing, sandbox=false) -> Dict{String,String}

Map library-relative paths to download URLs for the files of a Zenodo record.
File names in the record are flattened paths (see `flatten_path`).
"""
function zenodo_files(record::Integer; token=nothing, sandbox::Bool=false)
    token = token === nothing ? zenodo_token(sandbox) : token
    rec = _zenodo_request("GET", "$(zenodo_api(sandbox))/records/$record"; token)
    return Dict{String,String}(unflatten_path(f["key"]) => f["links"]["self"] for f in rec["files"])
end

"""
    publish_zenodo(relpaths, root; title, description, creators, access="restricted",
                   license="cc-by-4.0", keywords=String[], sandbox=false, token=nothing,
                   publish=true) -> Dict

Create a Zenodo record holding the given ensemble files (paths relative to
`root`), the lattice files they reference, and a matching `index.toml`, and
publish it. Returns the record
JSON including `id`, `doi` and `links`.

- `access`: "open", "restricted", or "embargoed".
- `creators`: vector of Dicts like `Dict("name" => "Wietek, Alexander", "affiliation" => "...", "orcid" => "...")`.
- `sandbox=true` targets sandbox.zenodo.org; use it first.
- `publish=false` leaves the record as an editable draft.

A published record cannot be deleted. A record may hold at most 100 files
and 50 GB.
"""
function publish_zenodo(relpaths::AbstractVector{<:AbstractString}, root::AbstractString;
                        title::AbstractString, description::AbstractString,
                        creators::AbstractVector, access::AbstractString="restricted",
                        license::AbstractString="cc-by-4.0", keywords=String[],
                        sandbox::Bool=false, token=nothing, publish::Bool=true)
    token = token === nothing ? zenodo_token(sandbox) : token
    token === nothing && error("no Zenodo token: set ZENODO_TOKEN (or ZENODO_SANDBOX_TOKEN)")
    api = zenodo_api(sandbox)

    # subset index, written to a temporary directory
    full = build_index(root; write=false)
    subset = filter(e -> e["path"] in relpaths, full["ensembles"])
    length(subset) == length(relpaths) ||
        error("some requested paths are not ensemble files under '$root'")
    lattices = unique([e["lattice"] for e in subset])
    uploads = vcat(collect(relpaths), lattices)
    length(uploads) + 1 <= 100 || error("Zenodo allows at most 100 files per record ($(length(uploads)) + index requested)")
    tmpdir = mktempdir()
    idxfile = joinpath(tmpdir, INDEX_FILE)
    open(idxfile, "w") do io
        TOML.print(io, Dict("schema_version" => SCHEMA_VERSION, "generated" => full["generated"],
                            "ensembles" => subset); sorted=true)
    end

    dep = _zenodo_request("POST", "$api/deposit/depositions"; token, body="{}")
    id = dep["id"]
    bucket = dep["links"]["bucket"]

    for rel in uploads
        open(joinpath(root, rel), "r") do io
            _zenodo_request("PUT", "$bucket/$(flatten_path(rel))"; token, input=io,
                            content_type="application/octet-stream")
        end
    end
    open(idxfile, "r") do io
        _zenodo_request("PUT", "$bucket/$INDEX_FILE"; token, input=io, content_type="application/octet-stream")
    end

    meta = Dict("metadata" => Dict(
        "title"       => title,
        "upload_type" => "dataset",
        "description" => description,
        "creators"    => creators,
        "access_right" => access,
        "license"     => license,
        "keywords"    => vcat(["METTS", "product states", "quantum many-body"], keywords),
    ))
    _zenodo_request("PUT", "$api/deposit/depositions/$id"; token, body=JSON.json(meta))

    if publish
        return _zenodo_request("POST", "$api/deposit/depositions/$id/actions/publish"; token)
    else
        return _zenodo_request("GET", "$api/deposit/depositions/$id"; token)
    end
end
