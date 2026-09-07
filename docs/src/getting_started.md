# Getting started

## Install the package

```julia
using Pkg
Pkg.add(url = "https://github.com/awietek/METTSLibrary.jl")
```

Once the package is registered this becomes `Pkg.add("METTSLibrary")`.

## Tell it where the data is

The package looks for files in this order and stops at the first hit:

1. **A local copy of the library**, given by the environment variable
   `METTSLIBRARY_PATH`. On the institute cluster this is the shared clone.
   No network is used.
2. **The download cache**, `METTSLIBRARY_CACHE`, default `~/.cache/mettslibrary`.
3. **Remotes.** By default the Hugging Face dataset repository. Others can be
   added, see [Sharing](@ref).

On the institute cluster the group's copy lives at
`/data/condmat/awietek/Data/mettslibrary`. Add one line to your `~/.bashrc`
(or `~/.zshrc`, whichever shell you use) and open a new shell:

```bash
export METTSLIBRARY_PATH=/data/condmat/awietek/Data/mettslibrary
```

That is all the configuration there is on the cluster: no account, no
token, no network. Check that Julia sees it with

```julia
ENV["METTSLIBRARY_PATH"]
```

Jobs submitted through the batch system inherit the variable from your login
environment if `.bashrc` is sourced; if not, put the same `export` line into
the job script.

On your own laptop point it at your clone in the same way, e.g.
`export METTSLIBRARY_PATH=$HOME/Research/Data/mettslibrary`.

## Access to the private repository

Away from a local copy the package downloads from Hugging Face. While the
repository is private this needs a token with read access, obtained in two
steps:

1. Create a free account at huggingface.co and ask the library maintainer to
   add you as a collaborator on the dataset.
2. In your account settings create a read token, then either run
   `huggingface-cli login` once, or put it in the environment variable
   `HF_TOKEN`. The package reads both.

If the token is missing the first download fails with a message that says
exactly this. Once the repository is public this section becomes irrelevant.

## First query

```julia
using METTSLibrary

idx = load_index()                 # index.toml from the local copy or the remote
length(idx["ensembles"])

ensembles(idx; model = "tJ")       # filter by any field
ensembles(idx; lattice_name = "square.L32.W4.cyl", beta = 4.0, t = 3.0, t_prime = -0.3)
ensembles(idx; nsites = 128)
```

Each result is a dictionary with the ensemble's path in the library, its
hash and size, and all metadata: model, site type, lattice name and number
of sites, temperature, couplings, quantum numbers, the collapse bases used,
the number of samples, and the names of the stored observables.

```julia
e = load(ensembles(idx; model = "tJ", beta = 4.0)[1])
```

`load` fetches the file and its lattice if they are not local, checks the
hashes, caches them, and returns an [`Ensemble`](@ref). Later calls for the
same file are served from the cache.
