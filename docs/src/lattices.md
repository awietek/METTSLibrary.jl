# Lattice files

Every ensemble references a lattice file. It defines where the MPS sites sit
and which bonds the Hamiltonian contains, and therefore the site ordering
that the stored states refer to. Without it a product state is a string of
integers with no meaning.

## Format

The TOML model file, the same format XDiag reads:

```toml
# comments are allowed and preserved
Coordinates = [
  [0.0, 0.0],
  [0.0, 1.0],
  [1.0, 0.0],
  [1.0, 1.0]
]

Interactions = [
  ['J',  'HB', 0, 1],       # [coupling, type, sites...]
  ['J',  'HB', 2, 3],
  ['Jd', 'HB', 0, 3]
]
```

- Site `i` in the file (0-based) is MPS site `i+1`.
- Each interaction is a coupling name, an interaction type, and the sites it
  acts on. Any number of sites is allowed, including none for global terms.
- Further keys, such as `Symmetries` and irreducible representations, are
  kept verbatim and ignored.

Every coupling name that appears in `Interactions` must have a value in the
ensemble's `parameters`. Couplings that do not appear in the lattice, such
as an on-site `U` that the code handles internally, may still be listed in
`parameters`.

## Reading and converting

```julia
txt = read_lattice("shastry.16.HB.J.Jd.fsl.toml")   # verbatim
txt = read_lattice("shastry.16.HB.J.Jd.fsl.lat")    # converted to TOML
l   = parse_lattice(txt)
nsites(l), l.coordinates, l.interactions, lattice_couplings(l)
```

[`read_lattice`](@ref) accepts the two older plain-text formats and converts
them: the C++ metts `[Coordinates]` / `[Interactions]` layout and the XDiag
`[Dimension]` / `[Sites]` / `[Interactions]` layout. Both list interactions
as `TYPE COUPLING sites`, which becomes `[coupling, type, sites]`.

## Generating one

Simulations that built their Hamiltonian in code have no lattice file.
[`square_lattice_toml`](@ref) generates one for a square lattice in the site
ordering of `ITensorMPS.square_lattice`:

```julia
lat = square_lattice_toml(32, 4; yperiodic = true,
                          bonds = [("t", "HOP", :nn), ("J", "HB", :nn), ("t_prime", "HOP", :nnn)])
```

`:nn` are nearest-neighbour bonds, `:nnn` both diagonals. Check that the
generated bond list matches what your driver actually did.

## Names

Every ensemble carries a `lattice_name`, a short identifier such as
`shastry.16.HB.J.Jd.fsl` or `square.L32.W4.cyl`. It names the directory in
the library that holds all ensembles on that lattice, and the lattice file
inside it. Choose names that identify the lattice completely, including
size and boundary conditions, since these are no longer stored anywhere
else. Using the basename of your lattice file is the natural choice, and the
converters do that by default.

## Where lattice files live

The lattice sits in its own directory, one level above the parameter
directories, and every ensemble below it refers to it:

```
tJ/superconductors/square.L32.W4.cyl/square.L32.W4.cyl.toml
tJ/superconductors/square.L32.W4.cyl/J=0.4_t=3.0/ndn=56_nup=56/T=…_beta=…/basis=X_seed=1.h5
```

An ensemble's HDF5 file records neither the lattice's name nor a checksum of
it. The lattice is the single `.toml` in the lattice directory — four levels
above the HDF5 file, since `<parameters>/<sector>/T=<T>_beta=<beta>` sit
between them — and the lattice's name is that directory's name.

Locating it by position is what keeps a rename cheap: move the directory and
the `.toml`, rerun [`build_index`](@ref), done, with no file rewritten. A name
stored inside every file would have to be rewritten everywhere instead, and
could silently disagree with the directory it sits in.

A lattice that does not belong to the data is caught by [`validate`](@ref),
which checks the site count against the states, that every coupling the lattice
names has a value in `parameters`, and that the file parses at all. That covers
substituting the wrong lattice or truncating one. It does not catch an edit
preserving the site count and coupling names while changing the bonds — the
lattice, like everything else here, is append-only.

The lattice file is written when the first ensemble on that lattice is
written and never overwritten. If you write an ensemble whose lattice differs
from the existing file of the same name, `write_ensemble` refuses: it is a
different lattice and needs a different name.
