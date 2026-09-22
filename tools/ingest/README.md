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

Each one enumerates the dumps, parses the metadata out of the paths, resolves
the lattice, deduplicates, checks tag uniqueness, and converts. They share a
shape and differ only where the data does.

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

## Environment

`OUT` (library root), `PROJECT`, `DRYRUN`, `RESUME`, `MAXFILES` (0 = all),
`SHARD`/`NSHARDS`, and `LATMATCH` where the driver supports it.
