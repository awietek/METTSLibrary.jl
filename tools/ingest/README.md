# Ingest drivers

One driver per cluster project, plus the SLURM array script that runs it.
These are the record of how the library's contents were produced: every file's
`provenance.source_path` points at a dump these scripts converted, and
`provenance.time_evolution_source` names the run script the tau values came
from. Keep them here rather than in a home directory, or the data stops being
reproducible the moment that directory is tidied.

| driver | project on the cluster | in the library |
|---|---|---|
| `ingest_optical.jl` | `hubbard.optical.lattice` (t-J despite the name) | `tJ/tj.mixed.dimension` |
| `ingest_sc_tj.jl` | `superconductors/tj` + `tj_bkup` | `tJ/superconductors` |
| `ingest_sc_hubbard.jl` | `superconductors/hubbard` + `hubbard_bkup` | `Hubbard/superconductors` |
| `ingest_kagome.jl` | `kagome.superconductors/hubbard` (METTS.jl, `samples.txt`) | **not ingested — see below** |
| `ingest_triangular.jl` | `hubbard.triangular.metts.v2` | `Hubbard/hubbard.triangular.metts` |

Each one enumerates the dumps, parses the metadata out of the paths, resolves
the lattice, deduplicates, checks tag uniqueness, and converts. They share a
shape and differ only where the data does.

`ingest_kagome.jl` is the odd one twice over. Its input is METTS.jl text output
rather than a legacy C++ dump, so it converts a lattice file as it goes and is
not an array job. And **its project was deliberately skipped**: the runs are
too short to be thermal (median 5 samples per run at the coldest temperature,
51 at the warmest) and the driver restarts the chain through sloppy DMRG on
every resubmission instead of continuing it. The reasoning and the numbers are
in the data repo's `CLAUDE.md` under "Projects deliberately skipped". The
driver is kept, working and dry-run clean, in case that call is revisited —
so do not treat its presence here as evidence the data was ingested.

## Things that are not obvious

**Run the dry run first** (`DRYRUN=1 MAXFILES=4`). It prints the target path
and the occupancy computed from the states against the declared sector, which
is the check that catches a misparsed sector or the wrong lattice. It also
runs the tag-uniqueness gate over every ensemble before anything is written —
that gate caught 120 colliding Hubbard ensembles that would otherwise have
aborted an array job partway through.

**Deduplicate on parsed values, never on path strings.** The live and backup
trees write the same run with cosmetically different spellings (`U.10` against
`U.10.00`), which are equal as Float64 and therefore produce the same library
path. A string key leaves them looking distinct and they then collide.

**`--constraint=icelake` is required**, not a preference: the Julia depot's
precompiled images were built on icelake, and the `medium` partition is mostly
zen5 and cascadelake, where loading HDF5 fails with *"Unable to find compatible
target in cached code image"* before any work starts.

**`NSHARDS` must be the full array width**, not `SLURM_ARRAY_TASK_COUNT`.
Requeueing a single shard with `--array=18` sets that variable to 1, which
silently disables sharding and makes one task process every run.

**Lattice files must exist before the array fans out.** `write_ensemble`
writes one only if absent, so parallel tasks starting from nothing race on the
same file. Create each lattice serially first (`LATMATCH=...` with
`MAXFILES=1`) unless the library already has it.

**`RESUME=1`** makes a task skip runs already written, so a requeue after a
timeout or a failed shard costs only what is missing.

**Not every ingest wants an array.** `ingest_kagome.jl` runs as a single serial
task on purpose. Its inputs are ~70-line `samples.txt` files, not seek-bound
multi-GB dumps, so the whole project converts in one job — and serial makes the
lattice-file race impossible, which in turn removes the `MAXFILES=1` pre-create
that leaked into the Hubbard array through `sbatch`'s default `--export=ALL`.
Reach for the array when the per-run cost justifies it, not by default.

**Parse the filename when the directory layout is inconsistent.**
`hubbard.triangular.metts.v2` writes its four METTS subtrees in **six**
different directory shapes — `t` sometimes in the lattice directory and
sometimes in the parameter directory, `T` and `U` swapping levels, and two more
shapes nested under `outfiles.metts.fixedt.save/outfiles.metts.fixedu/`. The
filename is uniform across all of them and carries the complete parameter set,
so `ingest_triangular.jl` parses that and never looks at the path. Encoding six
directory shapes would have been six chances to be wrong.

**Some lattices no longer exist and must be regenerated.** The
`hubbard.triangular.metts.v2` lattice files lived in the Flatiron *home* tree
(`/mnt/home/...`), and only `ceph` was copied to this cluster, so
`/data/condmat/awietek/flatiron` holds `ceph` and nothing else. They are rebuilt
by `triangular/make_lattices.py` on top of the original `latt.py` and
`create_triangular.py`, kept alongside it. Do not trust a regenerated lattice
without a check: that script's geometry reproduces all 37 surviving triangular
lattices byte-for-byte, and its bond assignment was recovered from the stored
energies (`triangular/solve_tp3.py`, `triangular/geom_check.py`) rather than
inferred from the lattice name.

**Rename a lattice when the cluster name does not say what it is.** The
triangular cylinders are `op` and `rect.op` on the cluster; they are the
standard **YC** and **XC** cylinders, so the library calls them that. `rect`
carries no information a reader can use, and the distinction is physical — the
two are not isomorphic (different degree profiles and adjacency spectra), not
merely two drawings of one lattice. `LIBRARY_LATTICE` in
`ingest_triangular.jl` holds the map. Renaming costs a `mv` of the directory
and its `.toml` plus a `build_index`, because nothing inside an `.h5` names its
lattice — which is the whole point of "the path carries names, the file carries
data".

**Observables are only attached when they provably line up.** The kagome driver
writes the energy per step, appends it to an extensible dataset, and saves the
sample *afterwards*, so a run killed in between leaves one extra energy — and a
run that was killed there and then resumed leaves that extra buried in the
middle, shifting everything after it. The two look identical in the data (on
resume the driver re-runs DMRG, so the repeat is not the same number); what
separates them is which file was written last. An unexplained length mismatch
drops the observable rather than storing a misaligned one.

## Environment

`OUT` (library root), `PROJECT`, `DRYRUN`, `RESUME`, `MAXFILES` (0 = all),
`SHARD`/`NSHARDS`, and `LATMATCH` where the driver supports it.
