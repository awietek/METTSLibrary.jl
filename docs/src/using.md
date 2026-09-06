# Using samples

## Finding ensembles

[`ensembles`](@ref) filters the index by keyword. Keys match the top-level
fields of an index entry (`model`, `site_type`, `lattice_name`, `nsites`,
`beta`, `nsamples`), or an entry of `parameters`, `sector`, `algorithm` or
`provenance`. Numbers are compared approximately, strings exactly.

```julia
ensembles(model = "tJ", lattice_name = "square.L32.W4.cyl")   # one lattice
ensembles(model = "tJ", t = 3.0, J = 0.4, nup = 56)           # by coupling and sector
ensembles(code = "METTS.jl", beta = 8.0)                      # by producing code
```

Geometry, size and boundary conditions are properties of the lattice, so you
select them through `lattice_name` or `nsites`, and read the details from
the lattice file itself.

The result is a vector of dictionaries. Two entries worth knowing:

- `entry["path"]` is the location in the library, e.g.
  `tJ/square.L32.W4.cyl/J0.4_t3.0_t_prime-0.3_ndn56_nup56/beta_4.0_chain01.h5`
- `entry["collapse_bases"]` tells you which bases the samples were collapsed
  in, which matters below.

## Loading

```julia
e = load(entry)            # from an index entry
```

Reading validates the file: the lattice file must be present and unchanged,
every state index must be a valid label, and if a sector is declared, every
Z-basis sample must have the declared particle numbers.

An [`Ensemble`](@ref) holds:

| field | content |
|---|---|
| `states` | `(nsites, nsamples)` matrix of `UInt8`, 0-based indices into `local_states` |
| `local_states` | the label table, e.g. `["Emp", "Up", "Dn"]`, in ITensors order |
| `basis`, `collapse_bases` | per-sample basis code and its label table, e.g. `["Z", "X"]` |
| `chain`, `step` | which Markov chain a sample came from, and where in it |
| `observables` | `Dict` of arrays whose last dimension is the sample, e.g. `"energy"`, `"entropy"` |
| `parameters`, `sector`, `beta` | the Hamiltonian couplings, quantum numbers, inverse temperature |
| `lattice`, `lattice_name` | the lattice file's TOML text and its name; `lattice(e)` parses it into coordinates and interactions |

## Checking thermalization

The per-sample observables exist so you can judge an ensemble before using it.

```julia
using Statistics
E = e.observables["energy"]
mean(E), std(E) / sqrt(length(E))

# per chain
for c in unique(e.chain)
    sel = e.chain .== c
    println("chain $c: ", mean(E[sel]), " over ", count(sel), " samples")
end
```

The `step` field lets you drop the beginning of each chain if the run did
not discard a burn-in itself.

## Product states for METTS.jl

[`initial_states`](@ref) returns states as vectors of 1-based ITensors state
indices, which is what `METTS.random_product_state` returns as well, so it
is a drop-in replacement for the start of a chain.

```julia
σs = initial_states(e, 100; basis = "Z", thin = 5, rng = MersenneTwister(1))
```

- `n` states are drawn without replacement from the eligible samples.
  Without `n` you get all eligible samples in file order.
- `basis` restricts to samples collapsed in that basis.
- `thin` keeps only every `thin`-th step along each chain, which reduces the
  autocorrelation between the states you draw.

Turning a state into an MPS with ITensors:

```julia
using ITensors, ITensorMPS, METTS
sites = siteinds("tJ", nsites(e); conserve_qns = true)
σ = σs[1]
names = [local_state_string(sites[i], σ[i]) for i in 1:nsites(e)]   # METTS.jl helper
psi = MPS(sites, names)
```

Or directly with the label strings the library stores:

```julia
psi = MPS(sites, state_labels(e, j))
```

## The collapse basis matters

A sample's labels only mean what its basis says. In the `"Z"` basis, `"Up"`
and `"Dn"` are eigenstates of ``S^z``. In the `"X"` basis they are
eigenstates of ``S^x`` on occupied sites, and holes are unaffected. The
ITensors names `"Up"` and `"Dn"` build Z-basis product states, so an X-basis
sample must be turned into ``(|\!\uparrow\rangle \pm |\!\downarrow\rangle)/\sqrt2``
on each occupied site by hand, and such a state is not an ``S^z`` eigenstate
and cannot be built with `conserve_qns = true`.

Use `basis = "Z"` in `initial_states` unless you know what you are doing,
and check `entry["collapse_bases"]` before choosing an ensemble. Particle
numbers are meaningful in both bases, magnetization only in Z.

## Recomputing observables from stored states

Given a stored ensemble at inverse temperature ``\beta``, the METTS
``|\psi_\sigma\rangle \propto e^{-\beta H/2}|\sigma\rangle`` can be regenerated
for each state with a single imaginary-time evolution, and any observable
averaged over the ensemble. The samples are already distributed according to
the diagonal Gibbs weights ``\langle\sigma|e^{-\beta H}|\sigma\rangle / Z``,
so no Markov chain is needed. This is the main payoff of the library.

```julia
using METTS
σs = initial_states(e; basis = "Z")
results = map(σs) do σ
    psi = MPS(sites, [local_state_string(sites[i], σ[i]) for i in eachindex(σ)])
    psi = timeevo_tdvp_extend(H, psi, -e.beta / 2; tau = 0.1, normalize = true, maxm = 1000)
    inner(psi', O, psi)
end
```

Two caveats. Samples within one chain are correlated, so thin or block
before estimating errors. And a state sampled at ``\beta`` is a warm start at
another temperature, not a sample of it; averaging over it there requires
running the chain.
