# METTSLibrary file format, schema version 1

One HDF5 file holds one *ensemble*: the classical product states obtained by
collapsing METTS for one Hamiltonian on one lattice, one quantum-number
sector and one temperature, together with cheap per-sample observables. Files are
self-describing: everything needed to interpret the states is stored as
attributes, so a file is complete on its own without this package.

Files are **append-only**. A file, once written, is never modified. New
samples for the same parameters go into a new file distinguished by its tag.

## Root attributes

| attribute        | type            | meaning |
|------------------|-----------------|---------|
| `schema`         | string          | always `"mettslibrary"` |
| `schema_version` | int             | `1` |
| `model`          | string          | physical model, e.g. `"tJ"`, `"Hubbard"`, `"Heisenberg"` |
| `site_type`      | string          | ITensors site type, e.g. `"tJ"`, `"Electron"`, `"S=1/2"` |
| `local_states`   | string array    | ordered local basis labels, e.g. `["Emp","Up","Dn"]` |
| `nsites`         | int             | number of MPS sites `N` |
| `nsamples`       | int             | number of samples `M` |
| `beta`           | float           | inverse temperature |
| `collapse_bases` | string array    | labels of the collapse bases used, e.g. `["Z"]` or `["Z","X"]` |
| `lattice_name`   | string          | name of the lattice, e.g. `"shastry.16.HB.J.Jd.fsl"`; the lattice file is `../<lattice_name>.toml` relative to this file |
| `lattice_sha256` | string          | SHA-256 of that file, so a modified lattice is detected |

## Datasets

| dataset        | type    | shape (Julia order) | meaning |
|----------------|---------|---------------------|---------|
| `coordinates`  | float64 | `(dim, N)`          | copy of the lattice file's `Coordinates`, for readers that do not want to parse TOML |
| `states`       | uint8   | `(N, M)`            | local state of site `i` in sample `j`, **0-based index into `local_states`** |
| `basis`        | uint8   | `(M,)`              | 0-based index into `collapse_bases` for each sample |
| `chain`        | int32   | `(M,)`              | Markov chain id of each sample |
| `step`         | int32   | `(M,)`              | step within its chain |
| `observables/<name>` | float64 | `(..., M)`    | per-sample observables; last dimension is the sample |

In C or Python (h5py) the shapes appear transposed: `states` is `(M, N)`.

## Groups holding scalar attributes

| group        | content |
|--------------|---------|
| `parameters` | Hamiltonian couplings as float attributes, e.g. `t = 3.0`, `J = 0.4`, `t_prime = -0.3` |
| `sector`     | conserved quantum numbers as int attributes: `nup`, `ndn`, `n` (total), `sz2` (twice Sz) |
| `algorithm`  | how the samples were produced: `code`, `code_version`, `tau`, `maxdim`, `cutoff`, `seed`, ... |
| `provenance` | `created` (ISO 8601 UTC), `creator`, and for converted files `source_format`, `source_path`, `source_sha256`, `converter` |

## The lattice file

Every ensemble references the lattice file that defines the Hamiltonian's
geometry and bond structure. The lattice is a separate TOML file stored in the
library tree at the lattice level, `<lattice_name>.toml` one directory above
the HDF5 file, so that every parameter set, temperature and chain of that
lattice shares it and is guaranteed to adhere to it. Geometry, system size
and boundary conditions are not stored separately: the lattice file is the
single source for them. It
is never overwritten: an ensemble whose lattice differs from the existing
file must reference a different file name. The format is
the TOML model file:

```toml
Coordinates = [
  [0.0, 0.0],
  [1.0, 0.0],
]
Interactions = [
  ['J', 'HB', 0, 1],       # [coupling, type, sites...], 0-based site indices
]
```

Further keys such as `Symmetries` and irreducible representations are kept
verbatim and ignored. Site `i` of the lattice file is MPS site `i+1`, so the
file fixes the site ordering. Every coupling name appearing in
`Interactions` must have a value in `parameters`; the package refuses to
write or read a file where this is not the case. Plain-text `.lat` files of
the two older formats are converted to this form on ingest.

## Conventions

- **Integer encoding.** `states` holds indices into `local_states`, 0-based in
  the file. The order of `local_states` is the ITensors order, so in Julia
  `stored + 1` is the ITensors state index and `local_states[stored + 1]` is
  the state name accepted by `MPS(sites, names)`.
- **Collapse basis.** A sample's meaning depends on its basis. In the `"Z"`
  basis the labels are eigenstates of Sz (and occupation). In the `"X"` basis,
  `"Up"`/`"Dn"` denote spin-x eigenstates on occupied sites; `"Emp"` and
  `"UpDn"` are unaffected. Particle numbers are meaningful in both bases,
  magnetization only in Z.
- **Sector check.** If `sector` is present, every Z-basis sample must have the
  declared particle numbers; the package validates this on read and write.
- **Site ordering.** `coordinates` records where MPS site `i` sits, so the
  snake or column ordering used by the simulation is recoverable.
- **Standard observable names.** `energy`, `energy2`, `entropy` (bipartite
  von Neumann entropy; a vector per sample for all cuts, or a scalar for the
  central cut), `maxdim`, `n`, `sz`, `double_occupancy`. Others are free.

## Location in the library

```
<model>/<lattice_name>/<lattice_name>.toml
<model>/<lattice_name>/<parameters and sector>/beta_<beta>_<tag>.h5

tJ/square.L32.W4.cyl/square.L32.W4.cyl.toml
tJ/square.L32.W4.cyl/J0.4_t3.0_t_prime-0.3_ndn56_nup56/beta_4.0_chain01.h5
tJ/square.L32.W4.cyl/J0.4_t3.0_t_prime-0.3_ndn56_nup56/beta_8.0_chain01.h5
tJ/square.L32.W4.cyl/J0.4_t3.0_t_prime-0.2_ndn56_nup56/beta_4.0_chain01.h5
Heisenberg/shastry.16.HB.J.Jd.fsl/shastry.16.HB.J.Jd.fsl.toml
Heisenberg/shastry.16.HB.J.Jd.fsl/J21.0_J2p1.0_J31.0_J3p1.0_Jd1.0/beta_1.0_seed1.h5
```

`index.toml` at the library root lists every ensemble with its path, SHA-256,
size, metadata, and the root-relative path and hash of its lattice file. It
is regenerated from the files by `build_index`.
