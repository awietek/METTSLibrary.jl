# tools — surveying a cluster archive of METTS runs

These scripts build an ensemble-level catalogue of METTS runs scattered across
a filesystem, so you can see what exists before converting any of it. They are
not part of the package: nothing here is loaded by `using METTSLibrary`, and
they are excluded from the tests.

They were written against `/data/condmat/awietek` (347,776 HDF5 files, 9.3 TB of
product states) and carry the conventions of that archive in their path parsers.
Point them at another archive and the parsers in `build_catalogue.jl` will need
extending — that is the part that knows how run directories are named.

Everything writes to one scan directory, set with `METTS_SCAN_DIR`
(default `./scan`). Each stage reads the previous stage's output, so run them
in order; each is restartable and none touches the source archive.

```bash
export METTS_SCAN_DIR=$PWD/scan
mkdir -p "$METTS_SCAN_DIR"

bash  tools/scan_archive.sh   /data/condmat/awietek   # -> all_files.list
bash  tools/classify_hdf5.sh                         # -> classified.tsv
bash  tools/find_state_files.sh                      # -> mates2.tsv, names.paths
bash  tools/infer_basis.sh                           # -> sz.tsv
julia tools/build_catalogue.jl                       # -> catalogue.tsv, collisions.tsv
julia tools/catalogue_data.jl                        # -> data.json
bash  tools/build_page.sh                            # -> catalogue.html
```

`classify_hdf5.sh` and `infer_basis.sh` split their work across
`METTS_SCAN_JOBS` processes (default 12) and call the matching `.jl` on one
chunk each; run the `.jl` directly only for a single chunk. The work is HDF5
metadata over NFS, so it is latency-bound — more workers than cores helps.

## What each stage does

**`scan_archive.sh`** — one traversal of the archive collecting every `.h5`,
`.hdf5`, `.txt` and anything named `*chkpt*`, `*sample*` or `*state*`. One
traversal, because an NFS walk of this size costs 20–40 minutes.

**`classify_hdf5.sh`** — opens every HDF5 file and records whether it holds
`/ProductState`, its site and step counts, and how many per-step scalars it
carries. Runs in 12 parallel processes; ~25 minutes for 350k files. Classify
per *file*: a group of files that share a name pattern does **not** share a
structure, and sampling one representative per group gives wrong answers.

**`find_state_files.sh`** — searches non-HDF5 files near stateless HDF5 for
collapsed states. Five encodings are tested, because runs used all of them:
`12: [1, 2, 3]` (METTS.jl), bare `[1, 2, 3]`, whitespace-separated `0..3` or
`1..4`, and **label names** `Up Dn Emp UpDn` (checkpoint files). An
integers-only search misses the last of these entirely.

**`infer_basis.sh`** — recovers the collapse basis where the filename does not
record it, by reading up to 400 samples per run and counting distinct total
S<sub>z</sub>: X-basis collapse does not conserve it, Z-basis does. It also
scans files whose basis *is* known, as a control — check that those come back
correct before trusting any inference.

**`build_catalogue.jl`** — parses run paths into (model, lattice, couplings,
sector, temperature, seed), groups them into ensembles, propagates the basis
within parameter families, drops checkpoints whose chain is already stored, and
assigns each run a library tag. Writes `catalogue.tsv`, one row per ensemble,
plus `collisions.tsv` for any tag it could not make unique.

**`catalogue_data.jl`** and **`build_page.sh`** — collapse ensembles into
parameter families and render the browsable page.

## Things this archive taught us, which a new one may repeat

- **Directory names lie.** A project called `superconductors` held both t-J and
  Hubbard runs; `hubbard.optical.lattice` was entirely t-J; and 4,271 Kanamori
  square-lattice runs sat under `kagome.superconductors`. Derive the model from
  the couplings and confirm it against the datasets in the file.
- **Backups are not copies.** `hubbard_bkup` was a pure subset of `hubbard`, but
  `tj_bkup` held 580 runs that existed nowhere else. 22.5% of all samples were
  duplicates.
- **Match runs on parsed values, never on filename text.** A dump and its own
  checkpoint format the same numbers differently
  (`t.1.00...T.0.01250` vs `t.1...T.0.0125`), so string substitution silently
  fails to pair them.
- **Check the failure signature, not the total.** Two parser bugs here produced
  entirely plausible aggregate numbers. `build_catalogue.jl` prints counts of
  impossible states (zero sectors, empty parameters, lattice equal to project
  name) for exactly this reason; they should all be 0.
