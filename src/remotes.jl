# ---------------------------------------------------------------------------
# Where files come from, in order of preference:
#   1. the local data root      (ENV METTSLIBRARY_PATH, e.g. a clone on the cluster)
#   2. the local download cache (ENV METTSLIBRARY_CACHE, default ~/.cache/mettslibrary)
#   3. remotes: Hugging Face by default, or anything added with add_remote!
# ---------------------------------------------------------------------------

const DEFAULT_HF_REPO = "awietek/mettslibrary"
const DEFAULT_HF_BASE = "https://huggingface.co/datasets/$DEFAULT_HF_REPO/resolve/main"

abstract type Remote end

"Any HTTPS tree where `base/<relpath>` is the file. Optional bearer token."
struct HTTPRemote <: Remote
    base::String
    token::Union{Nothing,String}
end

"A Zenodo record; files are flat, names are flattened paths."
struct ZenodoRemote <: Remote
    record::Int
    files::Dict{String,String}      # relpath => download url
    token::Union{Nothing,String}
    sandbox::Bool
end

const REMOTES = Remote[]

"""
    data_root() -> Union{Nothing,String}

Local copy of the library from the environment variable `METTSLIBRARY_PATH`,
or `nothing` if unset. Files found there are used directly, without network.
"""
function data_root()
    p = get(ENV, "METTSLIBRARY_PATH", "")
    return isempty(p) ? nothing : expanduser(p)
end

"""
    cache_dir() -> String

Where downloaded files are kept: `METTSLIBRARY_CACHE`, default `~/.cache/mettslibrary`.
"""
cache_dir() = expanduser(get(ENV, "METTSLIBRARY_CACHE", joinpath(homedir(), ".cache", "mettslibrary")))

"""
    hf_token() -> Union{Nothing,String}

Hugging Face token from `HF_TOKEN`, or from the file written by
`huggingface-cli login` (`\$HF_HOME/token`, default `~/.cache/huggingface/token`).
"""
function hf_token()
    t = get(ENV, "HF_TOKEN", "")
    isempty(t) || return strip(t)
    home = get(ENV, "HF_HOME", joinpath(homedir(), ".cache", "huggingface"))
    f = joinpath(expanduser(home), "token")
    isfile(f) && return strip(read(f, String))
    return nothing
end

default_remotes() = Remote[HTTPRemote(DEFAULT_HF_BASE, hf_token())]

active_remotes() = isempty(REMOTES) ? default_remotes() : REMOTES

"""
    add_remote!(base::AbstractString; token=nothing)

Add an HTTPS tree as a source, e.g. another Hugging Face repository
(`.../resolve/main`), an institute web server, or a bucket. Remotes are tried
in the order added, before the default Hugging Face repository.
"""
function add_remote!(base::AbstractString; token=nothing)
    r = HTTPRemote(String(rstrip(base, '/')), token)
    push!(REMOTES, r)
    return r
end

"""
    add_remote!(; zenodo::Integer, token=nothing, sandbox=false)

Add a Zenodo record as a source. The record's file list is fetched once.
`token` is a personal access token, needed only for restricted records.
"""
function add_remote!(; zenodo::Integer, token=nothing, sandbox::Bool=false)
    r = ZenodoRemote(Int(zenodo), zenodo_files(zenodo; token, sandbox), token, sandbox)
    push!(REMOTES, r)
    return r
end

"Forget all added remotes; the default Hugging Face remote is used again."
clear_remotes!() = (empty!(REMOTES); nothing)

_headers(token) = token === nothing ? Pair{String,String}[] : ["Authorization" => "Bearer $token"]

url_for(r::HTTPRemote, rel::AbstractString)   = r.base * "/" * rel
url_for(r::ZenodoRemote, rel::AbstractString) = get(r.files, rel, nothing)

# Download to a temporary name and rename, so an interrupted download leaves no half file.
function _download(url::AbstractString, dest::AbstractString, token)
    tmp = dest * ".part"
    try
        Downloads.download(url, tmp; headers=_headers(token))
        mv(tmp, dest; force=true)
    catch err
        rm(tmp; force=true)
        rethrow(err)
    end
end

const _ACCESS_HINT = ". The repository is private or restricted: set HF_TOKEN (or run `huggingface-cli login`) " *
                     "with a token that has read access, or set METTSLIBRARY_PATH to a local copy."

function _explain(err, url)
    err isa Downloads.RequestError || return sprint(showerror, err)
    s = err.response.status
    return "HTTP $s for $url" * (s in (401, 403, 404) ? _ACCESS_HINT : "")
end

"""
    fetch_file(relpath; sha256=nothing, refresh=false) -> local path

Return a local path for the library file `relpath`. Looks in the data root,
then the cache, then downloads from the remotes in order. If `sha256` is
given the file is verified after download and a cached copy with a wrong hash
is replaced. `refresh=true` forces a new download (used for `index.toml`).
"""
function fetch_file(rel::AbstractString; sha256=nothing, refresh::Bool=false)
    root = data_root()
    if root !== nothing
        p = joinpath(root, rel)
        if isfile(p) && !is_lfs_pointer(p)
            sha256 === nothing || sha256_file(p) == sha256 ||
                @warn "local file does not match the index hash" path=p
            return p
        end
    end

    cached = joinpath(cache_dir(), rel)
    if isfile(cached) && !refresh
        if sha256 === nothing || sha256_file(cached) == sha256
            return cached
        end
        @warn "cached file has wrong hash, re-downloading" path=cached
        rm(cached)
    end
    mkpath(dirname(cached))

    failures = String[]
    for r in active_remotes()
        url = url_for(r, rel)
        url === nothing && continue
        try
            _download(url, cached, r.token)
            if sha256 !== nothing && sha256_file(cached) != sha256
                rm(cached)
                push!(failures, "hash mismatch from $url")
                continue
            end
            return cached
        catch err
            push!(failures, _explain(err, url))
        end
    end
    error("could not fetch '$rel':\n  " * join(failures, "\n  "))
end
