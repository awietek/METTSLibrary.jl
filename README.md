# METTSLibrary.jl

Tools for a library of classical product states sampled by METTS
(minimally entangled typical thermal states) simulations. Building a METTS
chain is expensive; the collapsed product states it yields are tiny. Storing
them lets new simulations start from thermalized states, and lets new
observables be computed from a stored ensemble with one imaginary-time
evolution per sample instead of a new chain.

The package handles the file format, validation, indexing, fetching and
conversion. The data lives in a separate repository, `mettslibrary`, whose
`index.toml` the package reads. See `docs/schema.md` for the file format.

## Using samples

```julia
using METTSLibrary

# find ensembles by any metadata field, parameter or sector
found = ensembles(model = "tJ", lattice_name = "square.L32.W4.cyl", beta = 4.0, t = 3.0)

# fetch, verify and read one
e = load(found[1])

# product states as 1-based ITensors indices, ready for METTS.jl
σs = initial_states(e, 100; basis = "Z", thin = 5)

# or one sample as state names for MPS(sites, ...)
names = state_labels(e, 1)
```

Where files come from, in order:

1. `METTSLIBRARY_PATH`: a local copy of the library. On the institute
   cluster put `export METTSLIBRARY_PATH=/data/condmat/awietek/Data/mettslibrary`
   into your `~/.bashrc`. No network is used if the file is there.
2. `METTSLIBRARY_CACHE` (default `~/.cache/mettslibrary`): previously
   downloaded files.
3. Remotes: the Hugging Face dataset repository by default. While it is
   private you need a token with read access in `HF_TOKEN` or from
   `huggingface-cli login`. Add other sources with `add_remote!`, including a
   Zenodo record: `add_remote!(zenodo = 1234567)`.

## Adding samples

```julia
e = Ensemble(model = "tJ", site_type = "tJ",
             lattice = read_lattice("square.L32.W4.cyl.toml"),  # required: Coordinates + Interactions
             lattice_name = "square.L32.W4.cyl",                 # names the lattice directory
             beta = 4.0, parameters = Dict("t" => 3.0, "J" => 0.4, "t_prime" => -0.3),
             sector = Dict("nup" => 56, "ndn" => 56), collapse_bases = ["Z"],
             states = states,            # (nsites, nsamples) UInt8, 0-based
             observables = Dict("energy" => energies))

root = "/path/to/mettslibrary"
write_ensemble(root, e; tag = "chain01")   # -> tJ/square.L32.W4.cyl/.../beta_4.0_chain01.h5
build_index(root)
```

Then commit and push the new file and `index.toml`. Files are append-only:
never modify a written ensemble, add a new file with a new tag.

A lattice file is mandatory. It is the TOML model file with `Coordinates`
and `Interactions` (the format XDiag reads); older `.lat` files are converted
by `read_lattice`. Every coupling named in the lattice must have a value in
`parameters`. For simulations that built the Hamiltonian in code,
`square_lattice_toml` generates a matching file.

The lattice is written as a separate file `<lattice_name>.toml` in the
lattice directory, one level above the HDF5 files, so every parameter set on
that lattice shares it. Geometry, size and boundary conditions are read from
it and not stored separately. The HDF5 records the relative path and the file's
hash, and refuses to load if the lattice file is missing or changed.

## Converting legacy output

```julia
# C++ metts code: metadata comes from the run script and the .lat file
e = from_legacy_cpp("outfile.dump.h5"; lattice = "chain.nx.16.ny.1.hubbard.lat",
                    model = "Hubbard", site_type = "Electron",
                    temperature = 0.2, basis = "X",
                    parameters = Dict("T" => 1.0, "Tp" => 0.0, "U" => 10.0), sector = Dict("n" => 14))

# METTS.jl ttJ_metts.jl driver: everything is read from the output path
e = from_ttj_run(".../L.32.W.4/J.0.4000/t.3/t_prime.-0.3000/filling.0.87500/T.0.25000/D.1000/tau.0.10.cutoff.1.0e-08.seed.1")
```

## Sharing a subset through Zenodo

```julia
publish_zenodo(["tJ/cylinder_W4_L32/.../beta_4.0_chain01.h5"], root;
               title = "...", description = "...",
               creators = [Dict("name" => "Wietek, Alexander")],
               access = "restricted", sandbox = true)
```

Test with `sandbox = true` first; a published Zenodo record cannot be deleted.
