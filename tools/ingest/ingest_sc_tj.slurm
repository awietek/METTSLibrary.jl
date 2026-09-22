#!/bin/bash
#SBATCH --job-name=metts-sc-tj
#SBATCH --partition=medium
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --array=0-19%20
#SBATCH --constraint=icelake
#SBATCH --output=/home/awietek/metts-ingest/logs/%x.%A_%a.log

# --constraint=icelake is required, not a preference. The Julia depot at
# JULIA_DEPOT_PATH holds precompiled images built on an Intel icelake node, and
# Julia refuses a cache whose CPU target the running node cannot provide:
#
#   ERROR: Unable to find compatible target in cached code image
#          .../packages/HDF5/.../api/api.jl
#
# The `medium` partition is 89 icelake, 10 cascadelake and 168 zen5 (AMD), so
# without the constraint most nodes would fail -- task 18 of job 1324673 drew
# cascadelake and died after 12 minutes. 89 nodes is ample for a 20-task array.
# The alternative fix is to rebuild the depot with a multi-architecture
# JULIA_CPU_TARGET so one cache serves every node; until then, constrain.

# Ingest the t-J half of `superconductors`: ~3800 distinct runs across the
# tj/ and tj_bkup/ trees, deduplicated by step count.
#
# `medium` with 8 h, not `short` with 90 min: these dumps average 741 MB and
# reach 1.7 GB, roughly 1.3x the optical-lattice set, and the previous array
# came close to its limit. A task that overruns is not a disaster -- RESUME=1
# makes a requeue skip what is already written -- but the retry costs an hour.
#
# The lattice .toml files must all exist before this runs; write_ensemble
# creates one only if absent and parallel tasks would race. Two are already in
# the library; square.op.nx.4.ny.4.ttpJ is pre-created serially with LATMATCH.

set -u
umask 027          # the clone is the group's read copy; default here is 077

for f in /etc/profile.d/modules.sh /usr/share/Modules/init/bash; do
    [ -f "$f" ] && . "$f" && break
done
module load julia/1.12.0 2>/dev/null

export JULIA_DEPOT_PATH=/home/awietek/julia-packages
export OUT=/data/condmat/awietek/Data/mettslibrary
export PROJECT=superconductors
export RESUME=1
export DRYRUN=0
# NSHARDS must be the FULL array width, not the number of tasks in this
# submission. SLURM_ARRAY_TASK_COUNT reports the latter: requeueing one shard
# with `--array=18` sets it to 1, the driver's `nsh > 1` guard then skips
# sharding entirely, and that single task quietly processes every run instead
# of its twentieth -- racing the other tasks for the same output paths. Fixed
# at 20 to match --array=0-19; override by exporting NSHARDS explicitly.
export NSHARDS=${NSHARDS:-20}
export SHARD=${SLURM_ARRAY_TASK_ID:-0}

DRIVER=/home/awietek/metts-ingest/ingest_sc_tj.jl
PKG=/home/awietek/julia-packages/dev/METTSLibrary

echo "task $SHARD/$NSHARDS on $(hostname) at $(date -Is)"
julia --project="$PKG" "$DRIVER"
rc=$?
echo "task $SHARD exited $rc at $(date -Is)"
exit $rc
