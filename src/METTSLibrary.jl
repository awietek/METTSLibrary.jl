"""
    METTSLibrary

A library of classical product states (collapsed samples) produced by METTS
simulations, together with simple per-sample observables.

The package knows how to write, read, validate, index, fetch, and convert
ensemble files. It knows nothing about which ensembles exist: that is the
content of the data repository (`mettslibrary`), whose `index.toml` this
package reads.

See `docs/src/schema.md` for the file format.
"""
module METTSLibrary

using HDF5
using TOML
using SHA
using JSON
using Dates
using Random
using Downloads

include("schema.jl")
include("lattice.jl")
include("paths.jl")
include("io.jl")
include("validate.jl")
include("index.jl")
include("remotes.jl")
include("states.jl")
include("convert/common.jl")
include("convert/legacy_cpp.jl")
include("convert/metts_jl.jl")
include("zenodo.jl")

export Ensemble, nsites, nsamples, lattice, validate
export Lattice, parse_lattice, read_lattice, lattice_couplings, square_lattice_toml
export write_ensemble, read_ensemble
export build_index, load_index, ensembles, load
export add_remote!, clear_remotes!
export initial_states, state_labels
export from_legacy_cpp, from_samples_txt, from_ttj_run
export publish_zenodo

end # module
