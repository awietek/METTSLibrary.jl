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

The lattice's **name and checksum are deliberately not stored**. The lattice is
the single `.toml` in the lattice directory, four levels above the file, and
its name is that directory's name. Finding it by position rather than by a
stored name is what makes renaming a lattice a `mv` plus `build_index`, instead
of a rewrite of every file that references it. Whether the lattice belongs to
the data is settled by [`validate`](@ref) — site count, couplings covered by
`parameters`, and that it parses at all — rather than by a checksum.

## Datasets

| dataset        | type    | shape (Julia order) | meaning |
|----------------|---------|---------------------|---------|
| `coordinates`  | float64 | `(dim, N)`          | copy of the lattice file's `Coordinates`, for readers that do not want to parse TOML |
| `states`       | uint8   | `(N, M)`            | local state of site `i` in sample `j`, **0-based index into `local_states`** |
| `basis`        | uint8   | `(M,)`              | 0-based index into `collapse_bases` for each sample |
| `step`         | int32   | `(M,)`              | METTS step this sample was taken at |
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
<model>/<project>/README.md
<model>/<project>/<lattice_name>/<lattice_name>.toml
<model>/<project>/<lattice_name>/<parameters>/<sector>/T=<T>_beta=<beta>/<tag>.h5
```

Abbreviating the temperature directory to `T=0.25` for readability, a
project looks like this:

```
tJ/superconductors/README.md
tJ/superconductors/square.L32.W4.cyl/square.L32.W4.cyl.toml
tJ/superconductors/square.L32.W4.cyl/J=0.4_t=3.0_t_prime=-0.3/ndn=56_nup=56/T=0.25/basis=X_maxdim=2000_tau=0.1_seed=1.h5
tJ/superconductors/square.L32.W4.cyl/J=0.4_t=3.0_t_prime=-0.3/ndn=56_nup=56/T=0.25/basis=X_maxdim=2000_tau=0.1_seed=2.h5
tJ/superconductors/square.L32.W4.cyl/J=0.4_t=3.0_t_prime=-0.3/ndn=56_nup=56/T=0.125/basis=X_maxdim=2000_tau=0.1_seed=1.h5
tJ/superconductors/square.L32.W4.cyl/J=0.4_t=3.0_t_prime=-0.2/ndn=56_nup=56/T=0.25/basis=X_maxdim=2000_tau=0.1_seed=1.h5
Heisenberg/frustration/shastry.16.HB.J.Jd.fsl/shastry.16.HB.J.Jd.fsl.toml
Heisenberg/frustration/shastry.16.HB.J.Jd.fsl/J2=1.0_Jd=1.0/default/T=1.0/basis=Z_maxdim=512_seed=1.h5
```

Everything above the filename is physics: model, project, lattice, couplings,
sector, temperature. Two runs that share all of it belong to the same
ensemble and land in the same directory, so the tag — the filename — is
method: the algorithm parameters that were varied. See
[`default_tag`](@ref) for how it is built. A sector with no conserved
quantum numbers gives `default`. All datasets are stored chunked and
deflated.

Couplings and quantum numbers are written `name=value`. The `=` is not
decoration: keys may contain `_` and values may be negative, so `t_prime=-0.3`
is readable where `t_prime-0.3` is not. These names are for humans — the
authoritative values are the `parameters` and `sector` attributes in the file.

`index.toml` at the library root lists every ensemble with its path, SHA-256,
size, metadata, and the root-relative path and hash of its lattice file. It
is regenerated from the files by `build_index`.
