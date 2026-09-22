#!/bin/bash
#SBATCH --job-name=metts-ingest-optical
#SBATCH --partition=short
#SBATCH --time=01:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --array=0-19%20
#SBATCH --output=/home/awietek/metts-ingest/logs/%x.%A_%a.log

# Ingest the optical-lattice t-J runs into the library, one shard per array task.
#
# Conversion is seek-bound rather than CPU-bound, so one core per task is right;
# the memory only has to hold one run's ProductState at a time.
#
# Everything this job reads must live on a SHARED filesystem: the session
# scratchpad under /tmp is local to the submit host and invisible to the compute
# nodes, so the driver is staged to /home/awietek/metts-ingest first.
#
# Both lattice .toml files MUST already exist: write_ensemble writes one only if
# absent, and parallel tasks starting from nothing would race on the same file.
# The driver is run once per lattice with LATMATCH beforehand for exactly that.
#
# RESUME=1 makes tasks skip runs already written, so the array can be requeued
# after a timeout or partial failure without tripping the append-only check.

set -u

# The clone is the group's read copy; the default umask here is 077, which would
# create directories the group cannot enter.
umask 027

for f in /etc/profile.d/modules.sh /usr/share/Modules/init/bash; do
    [ -f "$f" ] && . "$f" && break
done
module load julia/1.12.0 2>/dev/null

export JULIA_DEPOT_PATH=/home/awietek/julia-packages
export OUT=/data/condmat/awietek/Data/mettslibrary
export PROJECT=tj.mixed.dimension
export RESUME=1
export DRYRUN=0
export NSHARDS=${SLURM_ARRAY_TASK_COUNT:-20}
export SHARD=${SLURM_ARRAY_TASK_ID:-0}

DRIVER=/home/awietek/metts-ingest/ingest_optical.jl
PKG=/home/awietek/julia-packages/dev/METTSLibrary

echo "task $SHARD/$NSHARDS on $(hostname) at $(date -Is)"
julia --project="$PKG" "$DRIVER"
rc=$?
echo "task $SHARD exited $rc at $(date -Is)"
exit $rc
