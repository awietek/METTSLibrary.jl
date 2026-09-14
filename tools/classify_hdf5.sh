#!/bin/bash
# Classify every HDF5 file in the archive, in parallel.
#
# classify_hdf5.jl handles one chunk; this splits the list and runs N of them.
# Serially this takes hours on a few hundred thousand files; at 12-way it is
# about 25 minutes, and the work is metadata-only so it is latency-bound, not
# CPU-bound -- more workers than cores is fine.
#
#   -> classified.tsv   path, status, bytes, nsites, nsteps, nscalars, ndatasets, names
set -eu
W=${METTS_SCAN_DIR:-$PWD/scan}
N=${METTS_SCAN_JOBS:-12}
HERE=$(cd "$(dirname "$0")" && pwd)

[ -s "$W/all_files.list" ] || { echo "all_files.list missing - run scan_archive.sh"; exit 1; }

cut -f2 "$W/all_files.list" | grep -E '\.(h5|hdf5)$' > "$W/h5_only.list"
echo "HDF5 files to classify: $(wc -l < "$W/h5_only.list")  in $N processes"

rm -f "$W"/cls_*.tsv "$W"/chunk_*
split -n "l/$N" -d "$W/h5_only.list" "$W/chunk_"
date +%H:%M:%S
for c in "$W"/chunk_*; do
    julia "$HERE/classify_hdf5.jl" "$c" "$W/cls_$(basename "$c").tsv" &
done
wait
date +%H:%M:%S

cat "$W"/cls_*.tsv > "$W/classified.tsv"
rm -f "$W"/cls_*.tsv "$W"/chunk_*
echo "classified: $(wc -l < "$W/classified.tsv")"
cut -f2 "$W/classified.tsv" | sort | uniq -c | sort -rn | sed 's/^/   /'
