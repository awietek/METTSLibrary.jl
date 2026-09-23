#!/bin/bash
#SBATCH --job-name=metts-tri
#SBATCH --partition=medium
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --array=0-19%20
#SBATCH --constraint=icelake
#SBATCH --output=/home/awietek/metts-ingest/logs/%x.%A_%a.log

# Ingest the `thermo` subtree of hubbard.triangular.metts.v2 into
# Hubbard/hubbard.triangular.metts: 4,613 runs of 5,095, median 245 samples.
# The other three METTS subtrees are a different campaign that never dumped
# collapsed states -- see the header of ingest_triangular.jl.
#
# --constraint=icelake is required, not a preference: the depot's precompiled
# images were built on icelake and `medium` is mostly zen5/cascadelake, where
# loading HDF5 dies with "Unable to find compatible target in cached code
# image" before any work starts.
#
# NSHARDS is the FULL array width, never SLURM_ARRAY_TASK_COUNT: requeueing a
# single shard with --array=N sets that to 1, silently disabling sharding so
# one task processes every run.
#
# The lattice files MUST exist in the library before this array fans out.
# Committing them under tools/ingest/triangular/lattices/ is not enough: that
# is only what the driver READS. write_ensemble still creates the library's own
# copy at <model>/<project>/<lattice>/<lattice>.toml the first time it is
# needed, and 20 parallel tasks would race on it exactly as in the earlier
# ingests. Copy all 13 into place first (see the project README).
#
# Each task repeats the header scan over all 5,095 dumps before sharding,
# because deduplication and the collapse-basis vote need the whole picture.
# That is ~20x the NFS reads but the same wall time, and at 5k files it is
# cheaper than maintaining a separate scan/combine/convert pipeline.

set -u
umask 027          # the clone is the group's read copy; the default here is 077

for f in /etc/profile.d/modules.sh /usr/share/Modules/init/bash; do
    [ -f "$f" ] && . "$f" && break
done
module load julia/1.12.0 2>/dev/null

export JULIA_DEPOT_PATH=/home/awietek/julia-packages
export OUT=/data/condmat/awietek/Data/mettslibrary
export PROJECT=hubbard.triangular.metts
export RESUME=1
export NSHARDS=${NSHARDS:-20}
export SHARD=${SLURM_ARRAY_TASK_ID:-0}

# Set EVERY variable the driver reads, including the ones whose default is what
# we want. sbatch exports the submitting shell's environment, so anything left
# unset here is inherited -- submitting from a shell that had just run a
# MAXFILES=1 smoke test once gave all 20 tasks MAXFILES=1.
export DRYRUN=0
export MAXFILES=0        # 0 = all runs in this shard
export SCANMAX=0         # 0 = scan everything; nonzero is a smoke test only

PKG=/home/awietek/julia-packages/dev/METTSLibrary
DRIVER=$PKG/tools/ingest/ingest_triangular.jl

echo "task $SHARD/$NSHARDS on $(hostname) at $(date -Is)"
julia --project="$PKG" "$DRIVER"
rc=$?
echo "task $SHARD exited $rc at $(date -Is)"
exit $rc
