# Adding samples

Adding samples means building an [`Ensemble`](@ref) in memory, writing it
into the library tree, rebuilding the index, and committing.

## Rules

- **Files are append-only.** An ensemble file, once written, is never
  modified. New samples for the same parameters go into a new file with a
  new tag. This keeps Git LFS storage from duplicating a growing file and
  keeps hashes and citations stable.
- **A lattice file is required**, and every coupling it names must have a
  value in `parameters`. See [Lattice files](@ref).
- **Metadata is complete.** Everything needed to interpret the states is in
  the file: model, site type, lattice, couplings, sector, temperature,
  collapse basis per sample, step, how the samples were produced.
- **One file is one METTS run.** A run's samples go into their own file,
  tagged with its seed, rather than being merged with other runs at the same
  parameters: if one run turns out to be flawed, the others stay usable.

## Collecting states during a METTS run

Your driver already has the collapsed state after each step. Collect it
together with the cheap observables you compute anyway.

```julia
using METTS, METTSLibrary

N = length(sites)
states  = UInt8[]                 # appended column by column
energy  = Float64[]
entropy = Float64[]
steps   = Int32[]

for step in 1:nmetts
    psi = evolve(psi, beta / 2)                    # your imaginary-time evolution
    push!(energy,  real(inner(psi', H, psi)))
    push!(entropy, entropy_von_neumann(psi, N ÷ 2))
    σ = collapse_with_qn(psi, "Z")                 # Vector{Int}, 1-based ITensors indices
    step > nwarm && (append!(states, UInt8.(σ .- 1)); push!(steps, step))
    psi = MPS(sites, [local_state_string(sites[i], σ[i]) for i in 1:N])
end
states = reshape(states, N, :)
```

The library stores 0-based indices into the label table, hence the `.- 1`.
`METTS.local_state_strings(sites[1])` gives the label table in the same
order as `METTSLibrary.LOCAL_STATES[site_type]`.

## Building the ensemble

```julia
e = Ensemble(
    model      = "tJ",
    site_type  = "tJ",
    lattice    = read_lattice("square.L32.W4.cyl.toml"),     # TOML text, see Lattice files
    lattice_name = "square.L32.W4.cyl",                      # names the directory and lattice file
    beta       = beta,
    parameters = Dict("t" => 3.0, "J" => 0.4, "t_prime" => -0.3),
    sector     = Dict("nup" => 56, "ndn" => 56),
    algorithm  = Dict{String,Any}("code" => "METTS.jl", "code_version" => "1.0.0",
                                  "tau" => 0.1, "maxdim" => 1000, "cutoff" => 1e-8, "seed" => seed),
    provenance = Dict{String,Any}("creator" => "A. Wietek", "project" => "pseudogap"),
    collapse_bases = ["Z"],
    states     = states,
    step       = steps,
    observables = Dict{String,Array{Float64}}("energy" => energy[nwarm+1:end],
                                              "entropy" => entropy[nwarm+1:end]),
)
```

If several bases were used, `collapse_bases = ["Z", "X"]` and `basis` holds a
0 or 1 per sample. `validate(e)` runs the consistency checks and is also
called by `write_ensemble`.

Observables can be any array whose last dimension is the sample: a scalar
energy per sample, or the entropy at every cut as an `(N-1, nsamples)`
matrix. Standard names are `energy`, `energy2`, `entropy`, `maxdim`, `n`,
`sz`, `double_occupancy`.

## Writing into the library

```julia
root = ENV["METTSLIBRARY_PATH"]                     # the clone
rel  = write_ensemble(root, e; tag = "seed$(seed)")
# "tJ/square.L32.W4.cyl/J=0.4_t=3.0_t_prime=-0.3/ndn=56_nup=56/beta=4.0/seed1.h5"
build_index(root)
```

[`write_ensemble`](@ref) derives the location from the metadata and returns
it. The `tag` distinguishes files with otherwise identical metadata, such as
different chains or batches, and is mandatory. An existing file at that
location is an error, never overwritten. The lattice file is written to
`tJ/square.L32.W4.cyl/square.L32.W4.cyl.toml` if it is not there yet, and a
different lattice under the same name is refused.

[`build_index`](@ref) walks the whole tree, reads each file's metadata and
hash, and rewrites `index.toml`. It takes a few seconds per gigabyte.

## Committing

```bash
cd $METTSLIBRARY_PATH
git add tJ/ index.toml
git commit -m "tJ W4 L32 t3 J0.4 tp-0.3 n0.875: beta 4 and 8, seeds 1-4"
git push
```

Git LFS handles the HDF5 files transparently. Commit a batch at a time rather
than one commit per file. Other clones receive the new ensembles with
`git pull`.

## Ingesting on the cluster

The intended flow for a group is that jobs write finished ensembles into an
inbox directory, and one person or a scheduled job ingests them:

```julia
for f in readdir(inbox; join = true)
    e = read_ensemble(f)                                  # validates
    write_ensemble(root, e; tag = tag_for(f))
    rm(f)
end
build_index(root)
```

Note that `read_ensemble` needs the lattice file next to the inbox file's
directory as well, so jobs should write into the inbox with
`write_ensemble(inbox, e; tag = ...)`, which places both.

followed by the commit above. This keeps a single identity committing and
keeps the tree clean.
