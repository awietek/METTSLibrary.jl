# ---------------------------------------------------------------------------
# Lattice files. The canonical form is the TOML model file used by XDiag and
# the user's lattice generators:
#
#   Coordinates  = [[x, y], ...]                 one entry per site, MPS order
#   Interactions = [[coupling, type, sites...], ...]   0-based site indices
#
# Further keys (Symmetries, irreps, comments) are kept verbatim but ignored.
# Two older plain-text `.lat` variants are converted to this TOML form:
#   [Coordinates] / [Interactions]        (C++ metts code)
#   [Dimension]=d / [Sites]=N / [Interactions]=n   (XDiag legacy)
# both with interaction lines `TYPE COUPLING sites...`.
# ---------------------------------------------------------------------------

"""
    Lattice

Parsed view of a lattice file: `coordinates` is `(dim, nsites)`,
`interactions` is a vector of `(coupling, type, sites)` with **0-based** site
indices exactly as in the file.
"""
struct Lattice
    coordinates::Matrix{Float64}
    interactions::Vector{Tuple{String,String,Vector{Int}}}
end

nsites(l::Lattice) = size(l.coordinates, 2)

"""
    parse_lattice(toml::AbstractString) -> Lattice

Parse the TOML text of a lattice file. Requires the `Coordinates` key;
`Interactions` may be absent (empty).
"""
function parse_lattice(toml::AbstractString)
    d = TOML.parse(String(toml))
    haskey(d, "Coordinates") || error("lattice file has no 'Coordinates' key")
    cs = d["Coordinates"]
    isempty(cs) && error("lattice file has empty 'Coordinates'")
    dim = length(cs[1])
    all(length(c) == dim for c in cs) || error("inconsistent coordinate dimension in lattice file")
    coords = Matrix{Float64}(undef, dim, length(cs))
    for (i, c) in enumerate(cs)
        coords[:, i] .= Float64.(c)
    end
    inter = Tuple{String,String,Vector{Int}}[]
    for (k, row) in enumerate(get(d, "Interactions", []))
        length(row) >= 2 || error("interaction $k must have at least a coupling and a type")
        (row[1] isa AbstractString && row[2] isa AbstractString) ||
            error("interaction $k: coupling and type must be strings, got $(row[1:2])")
        sites = Int[]
        for s in row[3:end]
            s isa Integer || error("interaction $k: site indices must be integers, got $s")
            push!(sites, s)
        end
        push!(inter, (String(row[1]), String(row[2]), sites))
    end
    return Lattice(coords, inter)
end

"Unique coupling names appearing in the interactions, in order of appearance."
lattice_couplings(l::Lattice) = unique(first.(l.interactions))
lattice_couplings(toml::AbstractString) = lattice_couplings(parse_lattice(toml))

"""
    read_lattice(path) -> String

Read a lattice file and return its TOML text. A `.toml` file is returned
verbatim; a `.lat` file (either legacy variant) is converted with `lat_to_toml`.
"""
function read_lattice(path::AbstractString)
    isfile(path) || error("lattice file not found: '$path'")
    ext = lowercase(splitext(path)[2])
    txt = read(path, String)
    ext == ".toml" && return (parse_lattice(txt); txt)     # fail early on a bad file
    ext == ".lat"  && return lat_to_toml(txt)
    error("lattice file must end in .toml or .lat: '$path'")
end

"""
    lat_to_toml(text) -> String

Convert a plain-text `.lat` lattice file to the TOML form (used by `read_lattice`). Leading comment
lines are preserved. Interaction lines `TYPE COUPLING sites...` become
`[COUPLING, TYPE, sites...]`. Symmetry and irrep sections are dropped.
"""
function lat_to_toml(text::AbstractString)
    comments = String[]
    coords = Vector{Vector{Float64}}()
    inter = Vector{Tuple{String,String,Vector{Int}}}()
    section = :none
    for raw in split(text, '\n')
        line = strip(raw)
        isempty(line) && continue
        if startswith(line, "#")
            section == :none && push!(comments, String(line))
            continue
        end
        if startswith(line, "[")
            tag = lowercase(strip(split(line, "=")[1], ['[', ']', ' ']))
            section = tag in ("coordinates", "sites") ? :coordinates :
                      tag == "interactions"            ? :interactions :
                      tag == "dimension"               ? section : :other
            continue
        end
        if section == :coordinates
            push!(coords, parse.(Float64, split(line)))
        elseif section == :interactions
            f = split(line)
            length(f) >= 2 || error("bad interaction line: '$line'")
            push!(inter, (String(f[2]), String(f[1]), parse.(Int, f[3:end])))
        end
    end
    isempty(coords) && error("no coordinates found in .lat text")
    return _lattice_toml(comments, coords, inter)
end

function _lattice_toml(comments, coords, inter)
    io = IOBuffer()
    for c in comments
        println(io, c)
    end
    isempty(comments) || println(io)
    println(io, "Coordinates = [")
    for (i, c) in enumerate(coords)
        print(io, "  [", join(repr.(Float64.(c)), ", "), "]")
        println(io, i < length(coords) ? "," : "")
    end
    println(io, "]")
    println(io)
    println(io, "Interactions = [")
    for (i, (cpl, typ, sites)) in enumerate(inter)
        print(io, "  ['", cpl, "', '", typ, "'")
        for s in sites
            print(io, ", ", s)
        end
        print(io, "]")
        println(io, i < length(inter) ? "," : "")
    end
    println(io, "]")
    return String(take!(io))
end

"""
    square_lattice_toml(Lx, Ly; xperiodic=false, yperiodic=true,
                        bonds=[("t", "HOP", :nn), ("J", "HB", :nn)]) -> String

Generate a lattice file for an `Lx` x `Ly` square lattice in ITensorMPS
`square_lattice` ordering. `bonds` lists `(coupling, type, range)` with range
`:nn` (nearest neighbours) or `:nnn` (both diagonals). Site indices are
0-based, as in the TOML format. Use this for simulations that built their
Hamiltonian in code rather than from a lattice file.
"""
function square_lattice_toml(Lx::Integer, Ly::Integer; xperiodic::Bool=false, yperiodic::Bool=true,
                             bonds=[("t", "HOP", :nn), ("J", "HB", :nn)])
    N = Lx * Ly
    coords = [[Float64((n - 1) ÷ Ly), Float64((n - 1) % Ly)] for n in 1:N]
    idx(x, y) = x * Ly + y                   # 0-based
    function wrap(x, y)
        if xperiodic; x = mod(x, Lx); elseif !(0 <= x < Lx); return nothing; end
        if yperiodic; y = mod(y, Ly); elseif !(0 <= y < Ly); return nothing; end
        return idx(x, y)
    end
    nn  = Tuple{Int,Int}[]
    nnn = Tuple{Int,Int}[]
    for x in 0:Lx-1, y in 0:Ly-1
        i = idx(x, y)
        for (dx, dy, list) in ((1, 0, nn), (0, 1, nn), (1, 1, nnn), (1, -1, nnn))
            j = wrap(x + dx, y + dy)
            j === nothing && continue
            j == i && continue
            push!(list, (i, j))
        end
    end
    # periodic wrap on Ly = 2 (or Lx = 2) produces each bond twice; keep one
    unique!(b -> minmax(b...), nn); unique!(b -> minmax(b...), nnn)
    inter = Tuple{String,String,Vector{Int}}[]
    for (cpl, typ, range) in bonds
        list = range == :nn ? nn : range == :nnn ? nnn : error("bond range must be :nn or :nnn")
        for (i, j) in list
            push!(inter, (String(cpl), String(typ), [i, j]))
        end
    end
    comments = ["# Generated by METTSLibrary.square_lattice_toml",
                "# Square lattice Lx=$Lx Ly=$Ly, xperiodic=$xperiodic, yperiodic=$yperiodic",
                "# Site ordering: ITensorMPS.square_lattice (n = x*Ly + y, 0-based)"]
    return _lattice_toml(comments, coords, inter)
end

