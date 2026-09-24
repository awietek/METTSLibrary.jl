#!/bin/bash
#SBATCH --job-name=metts-sinha-ttj
#SBATCH --partition=medium
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --constraint=icelake
#SBATCH --output=/home/awietek/metts-ingest/logs/%x.%j.log

# Ingest Aritra Sinha's t-t'-J cylinders into tJ/sinha.ttprime.cylinder:
# the ~1,800 runs a params line proves ran at tau = 0.2, on three lattices
# (square.L16.W4, L24.W4, L16.W6, all .cyl.ttpJ).
#
# SERIAL ON PURPOSE.
#
# Not for memory reasons -- that was measured and is a non-issue. These
# outfile.h5 run to 832 GB (112 TB over the selected runs) because of the
# four-fermion correlators, and `h5ls` on one needs 16-23 GB of RSS, which is
# how sixteen parallel readers put a login node into swap. But reading
# /energy and /svn BY NAME is O(1) in the file size: 0.4-4.1 s and a flat
# 0.42 GB maxrss from 40 MB to 832 GB. 32G here is generous.
#
# Serial is to avoid the two traps the earlier array ingests hit: no NSHARDS /
# SLURM_ARRAY_TASK_COUNT mismatch, and no two tasks racing to create the same
# lattice .toml on first write. At ~1-4 s per run, 1,800 runs finish serially
# in well under the walltime.
#
# --constraint=icelake is required: the Julia depot's precompiled images were
# built on icelake, and `medium` is 89 icelake against 10 cascadelake and 168
# zen5, where loading HDF5 dies with "Unable to find compatible target in
# cached code image" before any work starts.

set -u
umask 027          # the clone is the group's read copy; the default here is 077

for f in /etc/profile.d/modules.sh /usr/share/Modules/init/bash; do
    [ -f "$f" ] && . "$f" && break
done
module load julia/1.12.0 2>/dev/null

export JULIA_DEPOT_PATH=/home/awietek/julia-packages
export OUT=/data/condmat/awietek/Data/mettslibrary
export PROJECT=sinha.ttprime.cylinder
export RESUME=1

# Set EVERY variable the driver reads, including those whose default is what we
# want: sbatch exports the submitting shell's environment, so anything left
# unset here is silently inherited from it.
export DRYRUN=0
export MAXFILES=0          # 0 = all selected runs
export NSHARDS=1
export SHARD=0
export WANT_TAU=0.2        # only runs a params line proves ran at this tau
export ENERGY_SHIFT=1      # sample[j] <-> energy[j+1]; see the driver's header
export LARGE_H5_GB=Inf     # no size guard: a named-dataset read is O(1) in size
export REQUIRE_ENERGY=1    # a run with no energy series is not written at all

PKG=/home/awietek/julia-packages/dev/METTSLibrary
DRIVER=$PKG/tools/ingest/ingest_sinha_ttj.jl

echo "sinha t-t'-J ingest on $(hostname) at $(date -Is)"
julia --project="$PKG" "$DRIVER"
rc=$?
echo "exited $rc at $(date -Is)"
exit $rc
