# ---------------------------------------------------------------------------
# Converter for the C++ `metts` code (ITensor v3 + lime, 2020-2022).
#
# Output file `<outfile>.dump.h5` holds one dataset per observable with the
# METTS step as first (C-order) dimension, plus
#   /ProductState  float64 (nsteps, nsites)   0 Emp, 1 Up, 2 Dn, 3 UpDn
# No metadata is stored in the file; it lives in the run script, the lattice
# file (`.lat`, [Coordinates]/[Interactions]) and the coupling file.
# The collapse basis was chosen with `--updates x|z`. In the X basis, "Up"/"Dn"
# denote spin-x eigenstates on occupied sites; holes are unaffected.
# ---------------------------------------------------------------------------

const LEGACY_CPP_LABELS = ["Emp", "Up", "Dn", "UpDn"]   # value 0, 1, 2, 3

# Known per-step scalar datasets and their library names.
const LEGACY_CPP_OBSERVABLES = Dict(
    "H"          => "energy",
    "H2"         => "energy2",
    "Entropy"    => "entropy",
    "MaxBondDim" => "maxdim",
    "Number"     => "n",
    "Number2"    => "n2",
    "TotalSz"    => "sz",
    "TotalSz2"   => "sz2",
    "DoubleOcc"  => "double_occupancy",
    "DoubleOcc2" => "double_occupancy2",
)

"""
    from_legacy_cpp(dump_h5; lattice, model, site_type, beta=nothing, temperature=nothing,
                    basis="Z", parameters=Dict(), sector=Dict(), algorithm=Dict(),
                    lattice_name=nothing) -> Ensemble

Convert a `*.dump.h5` file of the legacy C++ METTS code. The file holds no
metadata, so everything must be supplied: the lattice file used for the run
(path to the `.lat`, or TOML text), the model and site type, either `beta` or
`temperature`, the collapse basis, and a value for every coupling named in
the lattice. `lattice_name` defaults to the lattice file's basename.
Every per-step scalar dataset is kept as an observable; known ones are
renamed (see `LEGACY_CPP_OBSERVABLES`), others keep their lowercased name.
"""
function from_legacy_cpp(dump_h5::AbstractString;
                         lattice::AbstractString, model::AbstractString, site_type::AbstractString,
                         beta=nothing, temperature=nothing, basis::AbstractString="Z",
                         parameters=Dict{String,Float64}(), sector=Dict{String,Int}(),
                         algorithm=Dict{String,Any}(), lattice_name=nothing)
    local_states = LOCAL_STATES[site_type]
    # legacy value -> 0-based library index, or -1 if the label does not exist for this site type
    lut = [(i = findfirst(==(l), local_states); i === nothing ? -1 : i - 1) for l in LEGACY_CPP_LABELS]

    states, obs = h5open(dump_h5, "r") do f
        haskey(f, "ProductState") || error("'$dump_h5' has no /ProductState dataset")
        ps = read(f["ProductState"])          # Julia sees (nsites, nsteps)
        N, M = size(ps)
        st = Matrix{UInt8}(undef, N, M)
        for j in 1:M, i in 1:N
            v = Int(round(ps[i, j]))
            0 <= v <= 3 || error("unexpected ProductState value $v")
            lut[v + 1] >= 0 || error("label '$(LEGACY_CPP_LABELS[v + 1])' not valid for site type '$site_type'")
            st[i, j] = lut[v + 1]
        end
        obs = Dict{String,Array{Float64}}()
        for k in keys(f)
            d = f[k]
            # per-step scalars are stored as (nsteps, 1) in C order -> (1, nsteps) in Julia
            k != "ProductState" && d isa HDF5.Dataset && size(d) == (1, M) || continue
            obs[get(LEGACY_CPP_OBSERVABLES, k, lowercase(k))] = vec(Float64.(read(d)))
        end
        st, obs
    end

    return _converted(states, 1:size(states, 2);
        source=dump_h5, source_format="legacy_cpp_dump_h5", code="metts (C++)",
        lattice, lattice_name, model, site_type, beta, temperature, basis,
        parameters, sector, algorithm, observables=obs)
end
