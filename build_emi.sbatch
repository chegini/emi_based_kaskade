#!/usr/bin/env bash
#SBATCH --job-name=build-emi
#SBATCH --partition=small
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=00:30:00
#SBATCH --output=logs/build-emi-%j.out
#SBATCH --error=logs/build-emi-%j.err

set -euo pipefail

BASE=/data/numerik/people/fchegini/project/fatemeh/MicroCard/current_papers/BDDC_petsc_SDC_2023/gatevariables/20Nov2023/BDDC
IMAGE=localhost/kaskade7-mpi-build-environment:latest
IMAGE_ARCHIVE=/data/numerik/people/fchegini/kaskade7-mpi-build-environment.tar
BUILD_TMP="$BASE/emiModelBddcSdc_new_scratch/build-tmp-${SLURM_JOB_ID:-manual}"

cd "$BASE/emiModelBddcSdc_new_scratch"
mkdir -p "$BUILD_TMP"

echo "Job ID: ${SLURM_JOB_ID:-unknown}"
echo "Node: $(hostname)"
echo "Start: $(date)"

if ! podman image exists "$IMAGE"; then
  podman load -i "$IMAGE_ARCHIVE"
fi

echo "MPI compiler: /usr/lib64/mpi/gcc/openmpi4/bin/mpicxx"

podman run --rm --pull=never \
  -e PATH=/usr/lib64/mpi/gcc/openmpi4/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin \
  -e LD_LIBRARY_PATH=/usr/lib64/mpi/gcc/openmpi4/lib64:/usr/lib64/mpi/gcc/openmpi4/lib \
  -e TMPDIR=/BDDC/emiModelBddcSdc_new_scratch/$(basename "$BUILD_TMP") \
  -v /data/numerik/people/fchegini/kaskade7_new:/kaskade7 \
  -v "$BASE:/BDDC" \
  -w /BDDC/emiModelBddcSdc_new_scratch \
  "$IMAGE" \
  bash -lc 'make clean KASKADE7=/kaskade7 MPI=1 MPICXX=/usr/lib64/mpi/gcc/openmpi4/bin/mpicxx && make -j2 KASKADE7=/kaskade7 MPI=1 MPICXX=/usr/lib64/mpi/gcc/openmpi4/bin/mpicxx'

echo "Finished: $(date)"
