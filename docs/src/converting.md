# Converting legacy output

Samples from earlier codes can be brought into the library. Each converter
returns a validated [`Ensemble`](@ref) which you then write and index as in
[Adding samples](@ref). All converters record where the data came from in
`provenance`: source format, source path, hash of the source file.

## C++ metts code (2020 to 2022)

The `outfile.dump.h5` files hold a `/ProductState` dataset with one row per
step, encoded 0 Emp, 1 Up, 2 Dn, 3 UpDn, plus per-step scalars such as `H`
and `Entropy`. Nothing else is in the file: the parameters live in the run
script, the geometry in the `.lat` file. So everything must be passed in.

```julia
e = from_legacy_cpp("outfile_restart.dump.h5";
        lattice     = "chain.nx.16.ny.1.hubbard.lat",   # the run's lattice file
        model       = "Hubbard", site_type = "Electron",
        temperature = 0.2,                              # or beta
        basis       = "X",                              # from --updates x|z
        parameters  = Dict("T" => 1.0, "Tp" => 0.0, "U" => 10.0),
        sector      = Dict("n" => 14),
        algorithm   = Dict{String,Any}("maxdim" => 200, "cutoff" => 1e-7, "seed" => 1))
```

The lattice file names couplings `T` and `Tp`, so both need values even if
`Tp` was zero. The `lattice_name` defaults to the lattice file's basename,
here `chain.nx.16.ny.1.hubbard`; pass `lattice_name = ...` to choose another. Per-step scalars are kept as observables under library names
(`H` becomes `energy`, `Entropy` becomes `entropy`, and so on).

## METTS.jl `ttJ_metts.jl` driver

The example driver writes `samples.txt` with one line per measurement,
`meas: [1, 2, 3, ...]`, in 1-based ITensors indices for the tJ site type,
and encodes the parameters in the output path. Samples are collapsed in the
X basis.

```julia
e = from_ttj_run(".../L.32.W.4/J.0.4000/t.3/t_prime.-0.3000/filling.0.87500/T.0.25000/D.1000/tau.0.10.cutoff.1.0e-08.seed.1")
```

Everything is read from the path. Since the driver built the Hamiltonian in
code, a lattice is generated with `square_lattice_toml` using couplings `t`,
`J` and `t_prime`, named `square.L<L>.W<W>.cyl`. Pass `lattice = ...` and
`lattice_name = ...` to use a real one instead.

For other METTS.jl drivers with the same `samples.txt` format use
[`from_samples_txt`](@ref) and supply the metadata yourself.

## What cannot be recovered

The `metts_heisenberg` code stored energies and correlators but no collapsed
states, and the `metts_checker` code kept only the last state as a restart
point. Those runs cannot be converted.

## Encoding pitfalls

The two legacy families disagree on integers: the C++ files are 0-based with
Emp = 0, METTS.jl uses ITensors' 1-based indices with Emp = 1. The converters
map through label names, never raw integers, so a hole never turns into a
spin. If you write a converter for another format, do the same.
