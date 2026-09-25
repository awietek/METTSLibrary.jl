#!/bin/bash
#SBATCH --job-name=metts-soares-j1j2
#SBATCH --partition=medium
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --constraint=icelake
#SBATCH --output=/home/awietek/metts-ingest/logs/%x.%j.log
#
# Ingest Rafael Soares' J1-J2 Heisenberg runs on triangular cylinders into
# Heisenberg/soares.j1j2.triangular: 1,112 runs on two lattices (12x6 and
# 18x9), J2 in {0, 0.125, 0.5}, 14 inverse temperatures, four bond dimensions.
# The library's first non-fermionic project.
#
# Deliberately NOT an array job, and not sharded. The whole source is 223 MB of
# small HDF5 files -- a few hundred kilobytes each -- so a serial pass finishes
# in minutes. Sharding would buy nothing and would reintroduce the hazard of
# two tasks racing to create the same lattice .toml on first write.
#
# The driver builds and verifies both lattices up front, before writing
# anything, so an unreadable or wrong-shaped lattice file stops the run
# immediately rather than half way through.
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
export PROJECT=soares.j1j2.triangular
export RESUME=1

# Set EVERY variable the driver reads: sbatch exports the submitting shell's
# environment, so anything left unset here is silently inherited from it.
export DRYRUN=0
export MAXFILES=0          # 0 = all runs
export NSHARDS=1
export SHARD=0
export ENERGY_SHIFT=1      # sample[j] <-> energy[j+1]; see the driver's header

PKG=/home/awietek/julia-packages/dev/METTSLibrary
DRIVER=$PKG/tools/ingest/ingest_soares_j1j2.jl

echo "soares J1-J2 ingest on $(hostname) at $(date -Is)"
julia --project="$PKG" "$DRIVER"
rc=$?
echo "exited $rc at $(date -Is)"
exit $rc
