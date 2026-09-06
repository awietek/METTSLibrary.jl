# METTSLibrary.jl

A library of classical product states sampled by METTS (minimally entangled
typical thermal states) simulations, and the Julia tools to use it.

## Why

A METTS simulation alternates imaginary-time evolution with a collapse onto
a classical product state. The chain of these collapsed states is expensive
to produce, but each state is tiny: one byte per site. Keeping them has two
payoffs.

- **Warm starts.** A new simulation for the same Hamiltonian can start from
  stored, thermalized states and skip the burn-in of the Markov chain.
- **New observables without a new chain.** Every stored state ``\sigma``
  defines a METTS ``|\psi_\sigma\rangle \propto e^{-\beta H/2}|\sigma\rangle``.
  Any observable can be evaluated on a stored ensemble with one imaginary-time
  evolution per sample, embarrassingly parallel, with no Markov chain at all.

## The pieces

| piece | what it is |
|---|---|
| **ensemble** | the collapsed states of one Hamiltonian on one lattice, in one quantum-number sector, at one temperature, plus cheap per-sample observables. One HDF5 file. |
| **lattice file** | a TOML file with the site coordinates and the interaction bonds. It is the single source for geometry, system size, boundary conditions and site ordering. Shared by all ensembles on that lattice, and referred to by a short `lattice_name`. |
| **library** | a directory tree of ensembles and lattice files, plus an `index.toml` listing every ensemble with its metadata and hash. Kept in a Git repository with LFS. |
| **METTSLibrary.jl** | this package. Reads, writes, validates, indexes, fetches and converts ensembles. Knows nothing about which ensembles exist. |

Code and data are separate. The package is public. The data repository
`mettslibrary` is private for now and lives in three places: a Hugging Face
dataset repository, a clone on the institute cluster that your group reads
from, and a clone on your laptop.

## In one screen

```julia
using METTSLibrary

# which ensembles exist for the t-J cylinder at beta = 4?
found = ensembles(model = "tJ", lattice_name = "square.L32.W4.cyl", beta = 4.0)

# fetch, verify and read one
e = load(found[1])

# thermalization check
using Statistics
mean(e.observables["energy"]), std(e.observables["energy"])

# 100 product states as ITensors state indices, thinned along the chain
σs = initial_states(e, 100; basis = "Z", thin = 5)
```

Continue with [Getting started](@ref).
