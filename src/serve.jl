# ---------------------------------------------------------------------------
# A browsable catalogue of a local library, served over HTTP on localhost.
#
#     using METTSLibrary; serve()
#
# Built on the `Sockets` stdlib rather than HTTP.jl on purpose: this serves a
# handful of routes to one person on 127.0.0.1, and HTTP.jl would pull a
# substantial dependency tree into a package whose surface is deliberately
# small. What is needed is request-line parsing, a query string, and four JSON
# endpoints.
#
# Everything the page shows already exists: the index carries parameters,
# sector, temperature, sample counts and algorithm settings, and the lattice
# files carry coordinates and bonds. Nothing here recomputes physics; the one
# endpoint that touches an .h5 is the per-sample series.
# ---------------------------------------------------------------------------

using Sockets

"""
    serve(root=nothing; port=8801, host=Sockets.localhost, verbose=true)

Serve a browsable catalogue of the library at `root` (default: the configured
`METTSLIBRARY_PATH`) on `http://localhost:<port>`. Blocks until interrupted.

The page lists projects, their lattices and ensembles, draws any lattice from
its coordinates and bonds, and plots the stored per-sample observables of a run.

Binds to localhost only. Paths arriving from the browser are never used to read
the filesystem directly — they are looked up in the index first, so a request
cannot escape the library root.
"""
function serve(root=nothing; port::Integer=8801, host=Sockets.localhost, verbose::Bool=true)
    r = root === nothing ? data_root() : String(root)
    r === nothing && error("no library root: pass one, or set METTSLIBRARY_PATH")
    isdir(r) || error("not a directory: '$r'")
    isfile(joinpath(r, INDEX_FILE)) ||
        error("no $(INDEX_FILE) in '$r' — run build_index first")

    ctx = _ServeContext(r)
    verbose && @info "METTSLibrary catalogue" root=r files=length(ctx.entries) url="http://localhost:$port"

    server = Sockets.listen(host, port)
    try
        while true
            sock = Sockets.accept(server)
            @async try
                _serve_one(sock, ctx)
            catch err
                err isa Base.IOError || @debug "request failed" exception=err
            finally
                close(sock)
            end
        end
    catch err
        err isa InterruptException || rethrow()
        verbose && @info "catalogue stopped"
    finally
        close(server)
    end
    return nothing
end

# Index is read once; a catalogue of a static library does not need to re-read
# it per request, and reading it per request would make every click hash files.
struct _ServeContext
    root::String
    entries::Vector{Dict{String,Any}}
    bypath::Dict{String,Dict{String,Any}}
end

function _ServeContext(root::AbstractString)
    entries = ensembles(load_index(source=root))
    bypath = Dict{String,Dict{String,Any}}(e["path"] => e for e in entries)
    return _ServeContext(String(root), entries, bypath)
end

# --- tiny HTTP ------------------------------------------------------------

function _serve_one(sock, ctx)
    line = readline(sock)
    isempty(line) && return
    parts = split(line)
    length(parts) >= 2 || return
    method, target = parts[1], parts[2]
    while true                      # drain headers; no request bodies are read
        h = readline(sock)
        isempty(h) && break
    end
    if method != "GET"
        return _respond(sock, "405 Method Not Allowed", "text/plain", "GET only")
    end
    path, query = _split_target(target)
    try
        _route(sock, path, query, ctx)
    catch err
        _respond(sock, "500 Internal Server Error", "application/json",
                 JSON.json(Dict("error" => sprint(showerror, err))))
    end
end

function _respond(sock, status, ctype, body)
    b = body isa String ? codeunits(body) : body
    write(sock, "HTTP/1.1 $status\r\n",
                "Content-Type: $ctype\r\n",
                "Content-Length: $(length(b))\r\n",
                "Cache-Control: no-store\r\n",
                "Connection: close\r\n\r\n")
    write(sock, b)
    return nothing
end

_json(sock, x) = _respond(sock, "200 OK", "application/json; charset=utf-8", JSON.json(x))

function _split_target(t)
    i = findfirst('?', t)
    i === nothing && return (_percent(t), Dict{String,String}())
    return (_percent(t[1:i-1]), _query(t[i+1:end]))
end

function _query(s)
    out = Dict{String,String}()
    for kv in split(s, '&'; keepempty=false)
        j = findfirst('=', kv)
        j === nothing ? (out[_percent(kv)] = "") :
                        (out[_percent(kv[1:j-1])] = _percent(kv[j+1:end]))
    end
    return out
end

function _percent(s)
    s = replace(s, '+' => ' ')
    occursin('%', s) || return s
    io = IOBuffer(); i = 1; b = codeunits(s)
    while i <= length(b)
        if b[i] == UInt8('%') && i + 2 <= length(b)
            write(io, parse(UInt8, String(b[i+1:i+2]); base=16)); i += 3
        else
            write(io, b[i]); i += 1
        end
    end
    return String(take!(io))
end

# --- routes ---------------------------------------------------------------

function _route(sock, path, q, ctx)
    if path == "/" || path == "/index.html"
        return _respond(sock, "200 OK", "text/html; charset=utf-8", CATALOGUE_HTML)
    elseif path == "/api/browse"
        return _json(sock, _browse(ctx, get(q, "path", "")))
    elseif path == "/api/lattice"
        return _json(sock, _lattice_for(ctx, get(q, "path", "")))
    elseif path == "/api/series"
        return _json(sock, _series_for(ctx, get(q, "path", ""), get(q, "obs", "")))
    elseif path == "/api/readme"
        return _json(sock, _readme_for(ctx, get(q, "path", "")))
    end
    return _respond(sock, "404 Not Found", "text/plain", "no such route: $path")
end

# The library's own layout is the navigation:
#   <model>/<project>/<lattice>/<parameters>/<sector>/T=<T>_beta=<beta>/<tag>.h5
# so browsing is grouping index entries by path prefix. Nothing else scales:
# a flat list of every project's lattices is unreadable past a few projects,
# and the archive has dozens more to come.
const _LEVELS = ["model", "project", "lattice", "parameters", "sector",
                 "temperature", "run"]

_parts(p) = isempty(p) ? String[] : String.(split(strip(p, '/'), '/'))

function _under(entrypath, pre)
    ep = split(entrypath, '/')
    length(ep) > length(pre) || return false
    for (i, c) in enumerate(pre)
        ep[i] == c || return false
    end
    return true
end

"""
    _browse(ctx, path) -> Dict

One level of the library, as directories do it. `path` is a prefix of an
entry's path; the response lists the distinct next components with their run
and sample counts, plus whatever is worth showing at that level (site counts
and couplings for a lattice, the parsed values for a parameter or sector
directory, T and beta for a temperature directory). At the deepest level the
children are runs and carry the detail the table shows.
"""
function _browse(ctx, path)
    pre = _parts(path)
    d = length(pre)
    d <= 6 || return Dict("error" => "path is deeper than the layout")
    sel = [e for e in ctx.entries if _under(e["path"], pre)]

    groups = Dict{String,Vector{Dict{String,Any}}}()
    order = String[]
    for e in sel
        name = split(e["path"], '/')[d+1]
        haskey(groups, name) || push!(order, name)
        push!(get!(groups, name, Dict{String,Any}[]), e)
    end

    kids = Dict{String,Any}[]
    for name in order
        g = groups[name]
        e1 = g[1]
        k = Dict{String,Any}(
            "name" => name,
            "path" => isempty(path) ? name : "$(rstrip(path,'/'))/$name",
            "nfiles" => length(g),
            "nsamples" => sum(x["nsamples"] for x in g))
        if d == 0                                   # a model: how many projects
            k["nprojects"] = length(unique(x["project"] for x in g))
        elseif d == 1                               # a project: how many lattices
            k["nlattices"] = length(unique(x["lattice_name"] for x in g))
            k["ntemperatures"] = length(unique(x["temperature"] for x in g))
        elseif d == 2                               # a lattice
            k["nsites"] = get(e1, "nsites", nothing)
            k["couplings"] = get(e1, "couplings", String[])
            k["ntemperatures"] = length(unique(x["temperature"] for x in g))
            k["bases"] = sort(unique(vcat((get(x, "collapse_bases", String[]) for x in g)...)))
        elseif d == 3                               # a parameter directory
            k["parameters"] = e1["parameters"]
            k["nsectors"] = length(unique(split(x["path"], '/')[5] for x in g))
        elseif d == 4                               # a sector directory
            k["sector"] = e1["sector"]
            k["ntemperatures"] = length(unique(x["temperature"] for x in g))
        elseif d == 5                               # a temperature directory
            k["temperature"] = e1["temperature"]
            k["beta"] = e1["beta"]
            # the sweep view needs the observable names to offer a choice
            k["observables"] = get(e1, "observables", String[])
        else                                        # a run
            alg = get(e1, "algorithm", Dict{String,Any}())
            k["name"] = replace(name, r"\.h5$" => "")
            k["runpath"] = e1["path"]
            k["temperature"] = e1["temperature"]
            k["parameters"] = e1["parameters"]
            k["sector"] = e1["sector"]
            k["collapse_bases"] = get(e1, "collapse_bases", String[])
            k["observables"] = get(e1, "observables", String[])
            k["bytes"] = get(e1, "bytes", 0)
            for f in ("maxdim", "seed", "tau", "cutoff", "nwarm", "nmetts")
                k[f] = get(alg, f, nothing)
            end
        end
        push!(kids, k)
    end

    # Sort each level the way a reader wants it, not by string accident.
    # `d` is the depth of the PREFIX, so the children are one level deeper:
    # d==5 lists temperature directories, d==6 lists the runs themselves.
    if d == 6
        sort!(kids; by = k -> (something(k["maxdim"], 0), something(k["seed"], 0)))
    elseif d == 5
        sort!(kids; by = k -> k["temperature"])
    else
        sort!(kids; by = k -> k["name"])
    end

    crumbs = [Dict("name" => pre[i], "path" => join(pre[1:i], "/")) for i in 1:d]
    return Dict("path" => path,
                "childLevel" => _LEVELS[min(d + 1, 7)],
                "crumbs" => crumbs,
                "nfiles" => length(sel),
                "nsamples" => sum(e["nsamples"] for e in sel; init = 0),
                "lattice" => d >= 3 ? join(pre[1:3], "/") : nothing,
                "children" => kids,
                # library totals travel with every response so the header does
                # not need a second round trip at each level
                "root" => ctx.root,
                "libfiles" => length(ctx.entries),
                "libsamples" => sum(e["nsamples"] for e in ctx.entries; init = 0))
end

"""
    _readme_for(ctx, path) -> Dict

The project's own `README.md`, as markdown for the page to render. Every
project carries one saying who computed the runs, with what code, and what is
known about their quality — which is exactly what someone browsing needs
before using the data, and it is otherwise invisible from inside the library.

`path` is `model/project`; it is checked against the index rather than joined
onto the root.
"""
function _readme_for(ctx, path)
    pre = _parts(path)
    length(pre) >= 2 || return Dict("error" => "need model/project")
    any(e -> _under(e["path"], pre[1:2]), ctx.entries) ||
        return Dict("error" => "no such project")
    f = joinpath(ctx.root, pre[1], pre[2], "README.md")
    isfile(f) || return Dict("path" => path, "markdown" => nothing)
    return Dict("path" => path, "markdown" => read(f, String))
end

"Coordinates and bonds of one lattice, ready to draw. `path` is model/project/lattice."
function _lattice_for(ctx, path)
    pre = _parts(path)
    length(pre) >= 3 || return Dict("error" => "need model/project/lattice")
    i = findfirst(e -> _under(e["path"], pre[1:3]), ctx.entries)
    i === nothing && return Dict("error" => "no such lattice")
    lattice = pre[3]
    rel = get(ctx.entries[i], "lattice", nothing)
    rel === nothing && return Dict("error" => "index entry carries no lattice path")
    lat = parse_lattice(read(joinpath(ctx.root, rel), String))
    coords = [[lat.coordinates[d, s] for d in 1:size(lat.coordinates, 1)]
              for s in 1:nsites(lat)]
    bonds = [Dict("coupling" => c, "type" => t, "sites" => s)
             for (c, t, s) in lat.interactions]
    return Dict("name" => lattice, "nsites" => nsites(lat),
                "dim" => size(lat.coordinates, 1),
                "coordinates" => coords, "bonds" => bonds,
                "couplings" => lattice_couplings(read(joinpath(ctx.root, rel), String)))
end

"""
Per-sample observables of one run.

`obs` restricts the reply to a single observable. That matters: the run level
overlays every run in a temperature directory, which can be 120 of them, and
returning all twelve observables for each would be tens of megabytes to send
one curve.

The requested path is resolved through the index rather than joined onto the
root, so a crafted path cannot read outside the library.
"""
function _series_for(ctx, path, obs_wanted="")
    e = get(ctx.bypath, path, nothing)
    e === nothing && return Dict("error" => "no such run in the index")
    ens = read_ensemble(joinpath(ctx.root, path))
    obs = Dict{String,Any}()
    for (k, v) in ens.observables
        ndims(v) == 1 || continue
        isempty(obs_wanted) || k == obs_wanted || continue
        obs[k] = collect(v)
    end
    return Dict("path" => path, "nsamples" => nsamples(ens), "nsites" => nsites(ens),
                "observables" => obs,
                # every name this run holds, so a filtered reply can still
                # populate the observable chooser
                "available" => sort([k for (k, v) in ens.observables if ndims(v) == 1]),
                "collapse_bases" => ens.collapse_bases,
                "parameters" => ens.parameters, "sector" => ens.sector,
                "temperature" => 1 / ens.beta, "beta" => ens.beta,
                "algorithm" => ens.algorithm)
end
