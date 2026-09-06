const SCHEMA_NAME = "mettslibrary"
const SCHEMA_VERSION = 1

"""
Ordered local basis labels per ITensors site type. The position in this
vector (0-based in files, 1-based in Julia) is the integer stored for a site.
The order matches ITensors' `val(s, name)` so that `stored + 1` is directly
usable as an ITensors state index.
"""
const LOCAL_STATES = Dict{String,Vector{String}}(
    "S=1/2"    => ["Up", "Dn"],
    "tJ"       => ["Emp", "Up", "Dn"],
    "Electron" => ["Emp", "Up", "Dn", "UpDn"],
)

"""
    Ensemble

One set of collapsed product states for one Hamiltonian, lattice, sector and
temperature, with per-sample observables. Maps one-to-one onto a single HDF5
file. Construct with keyword arguments. Required: `model`, `site_type`,
`lattice`, `lattice_name`, `beta`, `states`; everything else has a default.

Fields
- `model`          : physical model name, e.g. "tJ", "Hubbard", "Heisenberg"
- `site_type`      : ITensors site type, e.g. "tJ", "Electron", "S=1/2"
- `local_states`   : ordered local basis labels (see `LOCAL_STATES`)
- `lattice`        : TOML text of the lattice file (Coordinates, Interactions);
                     see `read_lattice`, `square_lattice_toml`. The lattice defines
                     the geometry, boundary conditions and site ordering.
- `lattice_name`   : short name of the lattice, e.g. "shastry.16.HB.J.Jd.fsl" or
                     "square.L32.W4.cyl". Names the lattice directory in the
                     library, which holds the lattice file `<lattice_name>.toml`.
- `beta`           : inverse temperature of the ensemble
- `parameters`     : coupling values; must cover every coupling named in `lattice`
- `sector`         : conserved quantum numbers, e.g. "nup" => 56, "ndn" => 56
- `algorithm`      : how the samples were made: code, code_version, tau, maxdim, cutoff, seed, ...
- `provenance`     : where the file came from: created, creator, source_format, source_path, ...
- `collapse_bases` : labels of collapse bases used, e.g. ["Z"] or ["Z", "X"]
- `states`         : (nsites, nsamples) UInt8, 0-based index into `local_states`
- `basis`          : (nsamples,) UInt8, 0-based index into `collapse_bases`
- `chain`, `step`  : Markov chain id and step within chain, per sample
- `observables`    : name => array whose last dimension is nsamples
"""
Base.@kwdef struct Ensemble
    model::String
    site_type::String
    local_states::Vector{String} = LOCAL_STATES[site_type]
    lattice::String
    lattice_name::String
    beta::Float64
    parameters::Dict{String,Float64} = Dict{String,Float64}()
    sector::Dict{String,Int} = Dict{String,Int}()
    algorithm::Dict{String,Any} = Dict{String,Any}()
    provenance::Dict{String,Any} = Dict{String,Any}()
    collapse_bases::Vector{String} = ["Z"]
    states::Matrix{UInt8}
    basis::Vector{UInt8} = zeros(UInt8, size(states, 2))
    chain::Vector{Int32} = ones(Int32, size(states, 2))
    step::Vector{Int32} = Int32.(1:size(states, 2))
    observables::Dict{String,Array{Float64}} = Dict{String,Array{Float64}}()
end

"Number of MPS sites of an ensemble or lattice."
nsites(e::Ensemble) = size(e.states, 1)
"Number of stored samples."
nsamples(e::Ensemble) = size(e.states, 2)

"Parsed lattice of an ensemble: coordinates and interactions."
lattice(e::Ensemble) = parse_lattice(e.lattice)

function Base.show(io::IO, e::Ensemble)
    print(io, "Ensemble(", e.model, ", ", e.lattice_name, " (", nsites(e), " sites)")
    print(io, ", beta=", e.beta, ", ", nsamples(e), " samples")
    isempty(e.parameters) || print(io, ", ", join(["$k=$v" for (k, v) in sort(collect(e.parameters))], " "))
    isempty(e.sector)     || print(io, ", ", join(["$k=$v" for (k, v) in sort(collect(e.sector))], " "))
    print(io, ")")
end
