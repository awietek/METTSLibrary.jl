#!/bin/bash
# Search non-HDF5 files for collapsed product states.
#
# Runs that wrote no /ProductState sometimes stored their states beside the
# HDF5 as text. Five encodings appear in practice and all are tested here --
# in particular the checkpoint files use *label names*, which an integers-only
# search misses completely.
#
#   mates2.tsv    kind <TAB> path, for every file holding states
#   names.paths   just the label-name files (checkpoints)
#
# Needs: classified.tsv (classify_hdf5.jl), all_files.list (scan_archive.sh)
set -u
W=${METTS_SCAN_DIR:-$PWD/scan}
AWK=$(command -v gawk || command -v awk)

echo "scan dir: $W"
[ -s "$W/classified.tsv" ] || { echo "classified.tsv missing - run classify_hdf5.jl"; exit 1; }

# Directories holding a stateless or unreadable HDF5, plus their parents: a
# state file may sit one level up from the run it belongs to.
awk -F'\t' '$2=="NOSTATES" || $2=="ERROR" {p=$1; sub(/\/[^\/]*$/,"",p); print p}' \
  "$W/classified.tsv" | sort -u > "$W/nodirs.list"
sed 's|/[^/]*$||' "$W/nodirs.list" | sort -u > "$W/noparents.list"
cat "$W/nodirs.list" "$W/noparents.list" | sort -u > "$W/alldirs.list"
echo "directories to search: $(wc -l < "$W/alldirs.list")"

# Every non-HDF5 regular file in them, whatever it is called. NB the
# directories must precede find's expression, hence sh -c "$@".
: > "$W/cand.list"
xargs -a "$W/alldirs.list" -d '\n' -n 300 sh -c \
  'find "$@" -maxdepth 1 -type f ! -name "*.h5" ! -name "*.hdf5" -size -64M -printf "%s\t%p\n" 2>/dev/null' _ \
  >> "$W/cand.list"
cut -f2 "$W/cand.list" | grep -v '^$' > "$W/cand.paths"
echo "candidate files: $(wc -l < "$W/cand.paths")"

# First 40 lines of each; stop at the first hit.
xargs -a "$W/cand.paths" -d '\n' -n 400 "$AWK" '
  FNR>40 { nextfile }
  /^[0-9]+ *: *\[/                                 { print "LABELLED\t" FILENAME; nextfile }
  /\[[0-9]+ *(, *[0-9]+)+\]/                       { print "BRACKET\t"  FILENAME; nextfile }
  /^(Up|Dn|Emp|UpDn)([ \t]+(Up|Dn|Emp|UpDn)){7,}/  { print "NAMES\t"    FILENAME; nextfile }
  /^[0-3]([ ,\t]+[0-3]){9,}[ \t]*$/                { print "INTS0\t"    FILENAME; nextfile }
  /^[1-4]([ ,\t]+[1-4]){9,}[ \t]*$/                { print "INTS1\t"    FILENAME; nextfile }
' 2>/dev/null > "$W/mates2.tsv"

awk -F'\t' '$1=="NAMES" {print $2}' "$W/mates2.tsv" | sort > "$W/names.paths"

echo "state-bearing files: $(wc -l < "$W/mates2.tsv")"
cut -f1 "$W/mates2.tsv" | sort | uniq -c | sort -rn | sed 's/^/   /'
echo "label-name files (checkpoints): $(wc -l < "$W/names.paths")"
