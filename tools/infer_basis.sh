#!/bin/bash
# Recover the collapse basis from the states, in parallel, with a control set.
#
# X-basis collapse does not conserve total Sz; Z-basis does. Files whose name
# records the basis (updates.x / updates.z) are scanned too, as calibration:
# if those do not come back correct, the inference cannot be trusted for the
# files that lack the token, and you should stop rather than guess.
#
#   -> sz.tsv   path, label, status, nsites, nsteps, nread, ndistinct_Sz, min, max, modecount
set -eu
W=${METTS_SCAN_DIR:-$PWD/scan}
N=${METTS_SCAN_JOBS:-12}
HERE=$(cd "$(dirname "$0")" && pwd)

[ -s "$W/classified.tsv" ] || { echo "classified.tsv missing - run classify_hdf5.sh"; exit 1; }

awk -F'\t' '$2=="STATES" {print $1}' "$W/classified.tsv" > "$W/st.paths"
grep 'updates\.x' "$W/st.paths" | awk 'NR%97==1' | head -80  > "$W/cal.paths" || true
grep 'updates\.z' "$W/st.paths" | awk 'NR%7==1'  | head -80 >> "$W/cal.paths" || true
grep -v 'updates\.' "$W/st.paths" > "$W/unk.paths" || true
cat "$W/cal.paths" "$W/unk.paths" > "$W/sz.paths"
echo "to scan: $(wc -l < "$W/sz.paths")  (calibration: $(wc -l < "$W/cal.paths"))"

rm -f "$W"/sz_szc_*.tsv "$W"/szc_*
split -n "l/$N" -d "$W/sz.paths" "$W/szc_"
for c in "$W"/szc_*; do
    julia "$HERE/infer_basis.jl" "$c" "$W/sz_$(basename "$c").tsv" &
done
wait
cat "$W"/sz_szc_*.tsv > "$W/sz.tsv"
rm -f "$W"/sz_szc_*.tsv "$W"/szc_*

echo
echo "CALIBRATION -- these must agree or the inference is not usable:"
awk -F'\t' '$2!="unknown" && $3=="OK" {
  print "  updates." $2 "  " ($7==1 ? "Sz fixed" : "Sz varies")
}' "$W/sz.tsv" | sort | uniq -c
echo
echo "inferred for unlabelled files:"
awk -F'\t' '$2=="unknown" && $3=="OK" {
  print "  " ($7>1 ? "Sz varies -> X" : ($6>=20 ? "Sz fixed -> Z" : "too few samples -> undetermined"))
}' "$W/sz.tsv" | sort | uniq -c
