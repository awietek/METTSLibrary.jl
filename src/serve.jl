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
    elseif path == "/api/summary"
        return _json(sock, _summary(ctx))
    elseif path == "/api/ensembles"
        return _json(sock, _ensembles_for(ctx, get(q, "project", ""), get(q, "lattice", "")))
    elseif path == "/api/lattice"
        return _json(sock, _lattice_for(ctx, get(q, "project", ""), get(q, "lattice", "")))
    elseif path == "/api/series"
        return _json(sock, _series_for(ctx, get(q, "path", "")))
    end
    return _respond(sock, "404 Not Found", "text/plain", "no such route: $path")
end

"Projects, their lattices and temperature ranges — enough to navigate."
function _summary(ctx)
    projects = Dict{String,Any}()
    for e in ctx.entries
        key = string(e["model"], "/", e["project"])
        p = get!(projects, key, Dict{String,Any}(
            "model" => e["model"], "project" => e["project"],
            "nfiles" => 0, "nsamples" => 0, "lattices" => Dict{String,Any}()))
        p["nfiles"] += 1
        p["nsamples"] += e["nsamples"]
        l = get!(p["lattices"], e["lattice_name"], Dict{String,Any}(
            "name" => e["lattice_name"], "nfiles" => 0, "nsamples" => 0,
            "nsites" => get(e, "nsites", nothing),
            "couplings" => get(e, "couplings", String[]),
            "temperatures" => Set{Float64}(), "bases" => Set{String}()))
        l["nfiles"] += 1
        l["nsamples"] += e["nsamples"]
        push!(l["temperatures"], e["temperature"])
        for b in get(e, "collapse_bases", String[]); push!(l["bases"], b); end
    end
    for (_, p) in projects, (_, l) in p["lattices"]
        l["ntemperatures"] = length(l["temperatures"])
        ts = sort(collect(l["temperatures"]))
        l["tmin"] = isempty(ts) ? nothing : first(ts)
        l["tmax"] = isempty(ts) ? nothing : last(ts)
        l["temperatures"] = ts
        l["bases"] = sort(collect(l["bases"]))
    end
    return Dict("root" => ctx.root,
                "nfiles" => length(ctx.entries),
                "nsamples" => sum(e["nsamples"] for e in ctx.entries; init=0),
                "projects" => [p for (_, p) in sort(collect(projects); by=first)])
end

"Every run of one lattice, trimmed to what the table shows."
function _ensembles_for(ctx, project, lattice)
    out = Dict{String,Any}[]
    for e in ctx.entries
        e["project"] == project || continue
        isempty(lattice) || e["lattice_name"] == lattice || continue
        alg = get(e, "algorithm", Dict{String,Any}())
        push!(out, Dict(
            "path" => e["path"], "temperature" => e["temperature"], "beta" => e["beta"],
            "parameters" => e["parameters"], "sector" => e["sector"],
            "nsamples" => e["nsamples"], "bytes" => get(e, "bytes", 0),
            "observables" => get(e, "observables", String[]),
            "collapse_bases" => get(e, "collapse_bases", String[]),
            "maxdim" => get(alg, "maxdim", nothing), "seed" => get(alg, "seed", nothing),
            "tau" => get(alg, "tau", nothing), "cutoff" => get(alg, "cutoff", nothing)))
    end
    sort!(out; by = r -> (r["temperature"], something(r["maxdim"], 0), something(r["seed"], 0)))
    return Dict("project" => project, "lattice" => lattice, "runs" => out)
end

"Coordinates and bonds of one lattice, ready to draw."
function _lattice_for(ctx, project, lattice)
    i = findfirst(e -> e["project"] == project && e["lattice_name"] == lattice, ctx.entries)
    i === nothing && return Dict("error" => "no such lattice")
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

The requested path is resolved through the index rather than joined onto the
root, so a crafted path cannot read outside the library.
"""
function _series_for(ctx, path)
    e = get(ctx.bypath, path, nothing)
    e === nothing && return Dict("error" => "no such run in the index")
    ens = read_ensemble(joinpath(ctx.root, path))
    obs = Dict{String,Any}()
    for (k, v) in ens.observables
        ndims(v) == 1 && (obs[k] = collect(v))
    end
    return Dict("path" => path, "nsamples" => nsamples(ens), "nsites" => nsites(ens),
                "observables" => obs,
                "collapse_bases" => ens.collapse_bases,
                "parameters" => ens.parameters, "sector" => ens.sector,
                "temperature" => 1 / ens.beta, "beta" => ens.beta,
                "algorithm" => ens.algorithm)
end
