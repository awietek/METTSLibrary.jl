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
tJ/square.L32.W4.cyl/square.L32.W4.cyl.toml
tJ/square.L32.W4.cyl/J0.4_t3.0_t_prime-0.3_ndn56_nup56/beta_4.0_seed1.h5
tJ/square.L32.W4.cyl/J0.4_t3.0_t_prime-0.2_ndn56_nup56/beta_4.0_seed1.h5
```

An ensemble's HDF5 file records the lattice name and the file's SHA-256; the
file is always `../<lattice_name>.toml` relative to the HDF5 file. Reading
checks the hash, so a lattice file that was edited after the fact is detected.

The lattice file is written when the first ensemble on that lattice is
written and never overwritten. If you write an ensemble whose lattice differs
from the existing file of the same name, `write_ensemble` refuses: it is a
different lattice and needs a different name.
