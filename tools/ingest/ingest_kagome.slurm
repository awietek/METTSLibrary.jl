#!/bin/bash
#SBATCH --job-name=metts-kagome
#SBATCH --partition=medium
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --constraint=icelake
#SBATCH --output=/home/awietek/metts-ingest/logs/%x.%j.log

# Ingest kagome.superconductors (Hubbard on kagome cylinders) into
# Hubbard/kagome.superconductors: 3,161 METTS.jl runs, three lattices.
#
# Deliberately NOT an array job. Each run here is a ~70-line samples.txt plus a
# small companion .h5 -- fractions of a second, against the ~156 s of a
# seek-bound legacy C++ dump. Serial finishes well inside the walltime and buys
# back everything sharding cost on the previous two ingests: no NSHARDS/
# SLURM_ARRAY_TASK_COUNT trap, no parallel tasks racing to create the same
# lattice .toml, and so no serial MAXFILES=1 pre-create to leak MAXFILES into
# the array through sbatch's default --export=ALL.
#
# --constraint=icelake is still required: the Julia depot's precompiled images
# were built on icelake, and `medium` is 89 icelake against 10 cascadelake and
# 168 zen5, where loading HDF5 dies with "Unable to find compatible target in
# cached code image" before any work happens.

set -u
umask 027          # the clone is the group's read copy; the default here is 077

for f in /etc/profile.d/modules.sh /usr/share/Modules/init/bash; do
    [ -f "$f" ] && . "$f" && break
done
module load julia/1.12.0 2>/dev/null

export JULIA_DEPOT_PATH=/home/awietek/julia-packages
export OUT=/data/condmat/awietek/Data/mettslibrary
export PROJECT=kagome.superconductors
export RESUME=1

# Set EVERY variable the driver reads, including the ones whose default is what
# we want: sbatch exports the submitting shell's environment, so anything left
# unset here is silently inherited from it.
export DRYRUN=0
export MAXFILES=0        # 0 = all runs
export NSHARDS=1
export SHARD=0

# Run the tracked copy directly rather than a duplicate under ~/metts-ingest,
# so what ran is exactly what is committed.
PKG=/home/awietek/julia-packages/dev/METTSLibrary
DRIVER=$PKG/tools/ingest/ingest_kagome.jl

echo "kagome ingest on $(hostname) at $(date -Is)"
julia --project="$PKG" "$DRIVER"
rc=$?
echo "exited $rc at $(date -Is)"
exit $rc
