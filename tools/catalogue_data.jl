#  catalogue.tsv -> JSON grouped into parameter families.
#  A family is one physical system (model, lattice, couplings, sector); the
#  temperature sweep lives inside it, so the review unit is the family.
const W = get(ENV, "METTS_SCAN_DIR", joinpath(pwd(), "scan"))

esc(s) = replace(String(s), '\\' => "\\\\", '"' => "\\\"")

pool = String[]; pos = Dict{String,Int}()
function intern(s)
    s = String(s)
    haskey(pos, s) && return pos[s]
    push!(pool, s); pos[s] = length(pool) - 1
    return pos[s]
end

mutable struct Fam
    project::String; model::String; site::String; lattice::String
    params::String; sector::String; nsites::Int
    fmts::Set{String}; bases::Set{String}; flags::Set{String}
    provs::Set{String}; maxdims::Set{String}
    nruns::Int; nsamples::Int; raw::Int; est::Int
    temps::Vector{Tuple{Float64,Float64,Int,Int,String,String,String}}
end

lines = readlines(joinpath(W, "catalogue.tsv"))
hdr = split(lines[1], '\t'); idx = Dict(String(h) => i for (i, h) in enumerate(hdr))

fams = Dict{String,Fam}(); order = String[]
for l in lines[2:end]
    f = split(l, '\t'); length(f) < 19 && continue
    g(n) = String(f[idx[n]])
    k = join([g("project"), g("model"), g("lattice"), g("params"), g("sector")], "|")
    if !haskey(fams, k)
        fams[k] = Fam(g("project"), g("model"), g("site_type"), g("lattice"),
                      g("params"), g("sector"), parse(Int, g("nsites")),
                      Set{String}(), Set{String}(), Set{String}(),
                      Set{String}(), Set{String}(), 0, 0, 0, 0,
                      Tuple{Float64,Float64,Int,Int,String,String,String}[])
        push!(order, k)
    end
    fm = fams[k]
    push!(fm.fmts, g("format"))
    for x in split(g("provenance"), '+');  isempty(x) || push!(fm.provs, String(x));   end
    for x in split(g("maxdims"), ',');     isempty(x) || push!(fm.maxdims, String(x)); end
    for b in split(g("bases"), ','); isempty(b) || push!(fm.bases, String(b)); end
    for x in split(g("flags"), ','); isempty(x) || push!(fm.flags, String(x)); end
    nr = parse(Int, g("nruns")); ns = parse(Int, g("nsamples"))
    fm.nruns += nr; fm.nsamples += ns
    fm.raw += parse(Int, g("raw_bytes")); fm.est += parse(Int, g("est_bytes"))
    push!(fm.temps, (parse(Float64, g("T")), parse(Float64, g("beta")), nr, ns,
                     g("bases"), g("flags"), g("seeds")))
end

rows = String[]
for k in order
    fm = fams[k]
    sort!(fm.temps, by = first)
    Ts = [t[1] for t in fm.temps]
    temps = join(["[" * string(t[1]) * "," * string(round(t[2], digits=4)) * "," *
                  string(t[3]) * "," * string(t[4]) * "," *
                  string(intern(t[5])) * "," * string(intern(t[6])) * ",\"" * esc(t[7]) * "\"]"
                  for t in fm.temps], ",")
    push!(rows, string("[",
        intern(fm.project), ",", intern(fm.model), ",", intern(fm.site), ",",
        intern(fm.lattice), ",", intern(fm.params), ",", intern(fm.sector), ",",
        fm.nsites, ",", length(fm.temps), ",", minimum(Ts), ",", maximum(Ts), ",",
        fm.nruns, ",", fm.nsamples, ",", fm.raw, ",", fm.est, ",",
        intern(join(sort(collect(fm.bases)), ",")), ",",
        intern(join(sort(collect(fm.flags)), ",")), ",",
        intern(join(sort(collect(fm.fmts)), ",")), ",",
        intern(join(sort(collect(fm.provs)), " + ")), ",",
        intern(join(sort(collect(fm.maxdims), by = x -> parse(Int, x)), ", ")), ",",
        "[", temps, "]]"))
end

open(joinpath(W, "data.json"), "w") do io
    print(io, "{\"pool\":[", join(["\"" * esc(p) * "\"" for p in pool], ","),
          "],\"fams\":[", join(rows, ","), "]}")
end
println("families: ", length(rows), "  ensembles: ", length(lines) - 1,
        "  pool: ", length(pool), "  bytes: ", filesize(joinpath(W, "data.json")))
