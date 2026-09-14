#!/bin/bash
# Inventory every candidate METTS artefact under /data/condmat/awietek:
#   - HDF5 files, classified by whether they hold /ProductState
#   - text files, classified by whether they hold collapsed product states
#   - HDF5 without ProductState that has a state-bearing text file beside it
# One traversal; groups by (project, filename pattern with digits normalised)
# and probes a representative of each group.
export PATH=/usr/local/hdf5/1.10.7/gnu-serial/bin:$PATH
W=${METTS_SCAN_DIR:-$PWD/scan}
ALL=$W/all_files.list
GRP_H5=$W/grp_h5.list
GRP_TX=$W/grp_tx.list
PROBE=$W/probe2.tsv

echo "############ 1. one traversal for HDF5 + text candidates"
date +%H:%M:%S
timeout 3000 find /data/condmat/awietek \
     -path /data/condmat/awietek/Data/mettslibrary -prune -o \
     -type f \( -name '*.h5' -o -name '*.hdf5' -o -name '*.txt' \
                -o -iname '*chkpt*' -o -iname '*sample*' -o -iname '*state*' \) \
     -printf '%s\t%p\n' 2>/dev/null > "$ALL"
RC=$?
date +%H:%M:%S
[ $RC -eq 124 ] && echo "WARNING: find hit the 3000s timeout - list may be incomplete"
echo "candidate files: $(wc -l < "$ALL")"
awk -F'\t' '$2 ~ /\.(h5|hdf5)$/ {n++; s+=$1} END {printf "  hdf5: %d files, %.2f GB\n", n, s/1024/1024/1024}' "$ALL"
awk -F'\t' '$2 !~ /\.(h5|hdf5)$/ {n++; s+=$1} END {printf "  text: %d files, %.2f GB\n", n, s/1024/1024/1024}' "$ALL"

# ---- group -----------------------------------------------------------------
group () {   # $1 = awk condition
  awk -F'\t' -v cond="$1" '
  {
    isb = ($2 ~ /\.(h5|hdf5)$/)
    if ((cond == "h5" && !isb) || (cond == "tx" && isb)) next
    path=$2
    if (match(path, /\/Projects\//)) { rest=substr(path,RSTART+RLENGTH); split(rest,a,"/"); proj=a[1] }
    else { rest=path; sub(/^\/data\/condmat\/awietek\/?/,"",rest); split(rest,a,"/"); proj=a[1] (a[2]!=""?"/" a[2]:"") }
    n=split(path,parts,"/"); base=parts[n]; gsub(/[0-9]+/,"#",base)
    k=proj"\t"base
    cnt[k]++; byt[k]+=$1; if (!(k in rep)) rep[k]=path
  }
  END { for (k in cnt) print cnt[k]"\t"byt[k]"\t"rep[k]"\t"k }' "$ALL" | sort -rn
}
group h5 > "$GRP_H5"
group tx > "$GRP_TX"
echo "distinct hdf5 groups: $(wc -l < "$GRP_H5")   text groups: $(wc -l < "$GRP_TX")"

: > "$PROBE"

echo
echo "############ 2. probe HDF5 groups for /ProductState"
i=0
while IFS=$'\t' read -r cnt bytes rep proj pat; do
  i=$((i+1)); [ $i -gt 250 ] && break
  info=$(timeout 30 h5ls "$rep" 2>/dev/null)
  if echo "$info" | grep -q '^ProductState '; then
    dims=$(echo "$info" | sed -n 's/^ProductState *Dataset {\([0-9]*\)[^,]*, \([0-9]*\)}.*/\1 \2/p')
    st=${dims% *}; si=${dims#* }
    sc=$(echo "$info" | grep -cE 'Dataset \{[0-9]+(/Inf)?, 1\}')
    printf 'H5_STATES\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$cnt" "$bytes" "$proj" "$pat" "${si:-0}" "${st:-0}" "$sc" "$rep" >> "$PROBE"
  elif [ -z "$info" ]; then
    printf 'H5_UNREADABLE\t%s\t%s\t%s\t%s\t0\t0\t0\t%s\n' "$cnt" "$bytes" "$proj" "$pat" "$rep" >> "$PROBE"
  else
    # no ProductState: does a state-bearing text file sit beside it?
    d=$(dirname "$rep"); mate=""
    for t in $(ls "$d" 2>/dev/null | grep -iE '\.txt$|chkpt|sample|state' | head -6); do
      h=$(head -c 3000 "$d/$t" 2>/dev/null)
      if echo "$h" | grep -qE '^[0-9]+ *: *\[|\[[0-9]+ *(, *[0-9]+)+\]|^[0-3]( +[0-3]){7,}|^[1-4]( +[1-4]){7,}'; then
        mate="$t"; break
      fi
    done
    printf 'H5_NOSTATES\t%s\t%s\t%s\t%s\t0\t0\t0\t%s\n' "$cnt" "$bytes" "$proj" "$pat" "${mate:-none}" >> "$PROBE"
  fi
done < "$GRP_H5"
echo "probed $i hdf5 groups"

echo
echo "############ 3. probe text groups for collapsed states"
i=0
while IFS=$'\t' read -r cnt bytes rep proj pat; do
  i=$((i+1)); [ $i -gt 250 ] && break
  h=$(head -c 4000 "$rep" 2>/dev/null)
  kind=TX_OTHER
  if   echo "$h" | grep -qE '^[0-9]+ *: *\['            ; then kind=TX_SAMPLES_LABELLED   # "12: [1, 2, 3]"
  elif echo "$h" | grep -qE '\[[0-9]+ *(, *[0-9]+)+\]'  ; then kind=TX_SAMPLES_BRACKET    # "[1, 2, 3]"
  elif echo "$h" | grep -qE '^[0-3]( +[0-3]){7,}'       ; then kind=TX_STATES_0BASED      # "0 1 2 0 ..."
  elif echo "$h" | grep -qE '^[1-4]( +[1-4]){7,}'       ; then kind=TX_STATES_1BASED      # "1 2 3 1 ..."
  fi
  printf '%s\t%s\t%s\t%s\t%s\t0\t0\t0\t%s\n' "$kind" "$cnt" "$bytes" "$proj" "$pat" "$rep" >> "$PROBE"
done < "$GRP_TX"
echo "probed $i text groups"

echo
echo "############ 4. summary by class"
awk -F'\t' '{n[$1]+=$2; b[$1]+=$3} END {printf "%-18s %9s %12s\n","class","files","GB";
  for (k in n) printf "%-18s %9d %9.2f GB\n", k, n[k], b[k]/1024/1024/1024}' "$PROBE" | sort -k2 -rn

echo
echo "############ 5. HDF5 WITH product states, by project"
printf '%-32s %8s %10s %6s %7s %5s %11s\n' project files raw_GB sites steps scal est_lib
awk -F'\t' '$1=="H5_STATES" {
  k=$4; n[k]+=$2; b[k]+=$3; si[k]+=$6*$2; st[k]+=$7*$2; sc[k]+=$8*$2
  est[k]+=$2*(0.242*$6*$7 + 4.4*$8*$7 + 4*$7 + 5000)
} END { for (k in n) printf "%-32s %8d %10.2f %6d %7d %5d %8.3f GB\n",
        k, n[k], b[k]/1024/1024/1024, si[k]/n[k], st[k]/n[k], sc[k]/n[k], est[k]/1024/1024/1024 }' "$PROBE" | sort -k2 -rn
awk -F'\t' '$1=="H5_STATES" {n+=$2; b+=$3; e+=$2*(0.242*$6*$7+4.4*$8*$7+4*$7+5000)}
  END {printf "\nTOTAL: %d files, %.2f GB raw -> %.2f GB in the library\n", n, b/1024/1024/1024, e/1024/1024/1024}' "$PROBE"

echo
echo "############ 6. TEXT files that look like collapsed states"
printf '%-30s %8s %10s  %-28s %s\n' project files GB pattern class
awk -F'\t' '$1 ~ /^TX_SAMPLES|^TX_STATES/ {printf "%-30s %8d %10.3f  %-28s %s\n", $4, $2, $3/1024/1024/1024, $5, $1}' "$PROBE" | sort -k2 -rn | head -30
echo "(none listed above = no text-based sample files found)"

echo
echo "############ 7. HDF5 WITHOUT product states - is there a state file beside it?"
printf '%-30s %8s %10s  %-26s %s\n' project files GB pattern mate
awk -F'\t' '$1=="H5_NOSTATES" {printf "%-30s %8d %10.2f  %-26s %s\n", $4, $2, $3/1024/1024/1024, $5, $9}' "$PROBE" | sort -k2 -rn | head -25
echo
echo "--- those WITH a state-bearing neighbour:"
awk -F'\t' '$1=="H5_NOSTATES" && $9!="none" {printf "  %s  (%d files)  mate=%s\n", $4"/"$5, $2, $9}' "$PROBE" | head -20
echo "(none listed = every stateless HDF5 has no sample file next to it)"
