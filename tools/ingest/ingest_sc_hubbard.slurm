#!/bin/bash
#SBATCH --job-name=metts-sc-hub
#SBATCH --partition=medium
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --array=0-19%20
#SBATCH --constraint=icelake
#SBATCH --output=/home/awietek/metts-ingest/logs/%x.%A_%a.log

# Ingest the Hubbard half of `superconductors` into Hubbard/superconductors:
# ~5000 distinct runs across hubbard/metts/{outfiles,outfiles.metts} and
# hubbard_bkup, deduplicated by step count.
#
# --constraint=icelake is required, not a preference: the Julia depot's
# precompiled images were built on icelake, and the `medium` partition is 89
# icelake against 10 cascadelake and 168 zen5, where loading HDF5 dies with
# "Unable to find compatible target in cached code image" before any work.
#
# NSHARDS is the FULL array width, never SLURM_ARRAY_TASK_COUNT: requeueing a
# single shard with --array=N sets that to 1, which would silently disable
# sharding and have one task process every run.
#
# The lattice .toml must exist before the array fans out -- write_ensemble
# creates one only if absent and parallel tasks would race. There is only one
# Hubbard lattice here, so a single serial MAXFILES=1 run creates it.

set -u
umask 027          # the clone is the group's read copy; the default here is 077

for f in /etc/profile.d/modules.sh /usr/share/Modules/init/bash; do
    [ -f "$f" ] && . "$f" && break
done
module load julia/1.12.0 2>/dev/null

export JULIA_DEPOT_PATH=/home/awietek/julia-packages
export OUT=/data/condmat/awietek/Data/mettslibrary
export PROJECT=superconductors
export RESUME=1
export DRYRUN=0
export NSHARDS=${NSHARDS:-20}
export SHARD=${SLURM_ARRAY_TASK_ID:-0}

# Set EVERY variable the driver reads, even the ones whose default is what we
# want. sbatch exports the submitting shell's environment (--export=ALL is the
# default), so anything left unset here is inherited: submitting from a shell
# that had just run the serial lattice pre-create with MAXFILES=1 gave all 20
# tasks MAXFILES=1, and each converted one already-written run and stopped.
# The batch script is the authority; the submitting shell is not.
export MAXFILES=0        # 0 = all runs in this shard
export LATMATCH=""       # no lattice filter

DRIVER=/home/awietek/metts-ingest/ingest_sc_hubbard.jl
PKG=/home/awietek/julia-packages/dev/METTSLibrary

echo "task $SHARD/$NSHARDS on $(hostname) at $(date -Is)"
julia --project="$PKG" "$DRIVER"
rc=$?
echo "task $SHARD exited $rc at $(date -Is)"
exit $rc
