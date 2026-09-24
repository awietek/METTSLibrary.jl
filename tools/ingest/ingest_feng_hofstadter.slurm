#!/bin/bash
#SBATCH --job-name=metts-feng-hof
#SBATCH --partition=medium
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --constraint=icelake
#SBATCH --output=/home/awietek/metts-ingest/logs/%x.%j.log

# Ingest Chunhan Feng's Hofstadter-Hubbard triangular cylinders (MBD2000) into
# Hubbard/feng.hofstadter.triangular: 16,000 runs, four lattices, 10 values of
# U, 40 temperatures, 10 seeds, all at half filling in a fixed pi/2 flux.
#
# SERIAL ON PURPOSE, for two independent reasons.
#
# 1. The per-sample energies live only in checkpoint.bin, a ~9 MB Julia
#    Serialization blob per run, so a full pass reads ~144 GB. That is fine
#    one file at a time and unkind to the filesystem sixteen at a time.
# 2. The lattice .toml is created lazily on the first write for each system
#    size. Parallel tasks would race to create the same file.
#
# RESUME=1, so if 24 h is not enough, resubmitting continues where it stopped
# rather than redoing the work.
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
export PROJECT=feng.hofstadter.triangular
export RESUME=1

# Set EVERY variable the driver reads: sbatch exports the submitting shell's
# environment, so anything left unset here is silently inherited from it.
export DRYRUN=0
export MAXFILES=0          # 0 = all runs
export NSHARDS=1
export SHARD=0
export ENERGY_SHIFT=1      # sample[j] <-> energy[j+1]; see the driver's header
export TMAX=Inf            # all 40 temperatures

PKG=/home/awietek/julia-packages/dev/METTSLibrary
DRIVER=$PKG/tools/ingest/ingest_feng_hofstadter.jl

echo "feng Hofstadter ingest on $(hostname) at $(date -Is)"
julia --project="$PKG" "$DRIVER"
rc=$?
echo "exited $rc at $(date -Is)"
exit $rc
