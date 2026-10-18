#!/bin/bash
#SBATCH --job-name=emi-mpi8
#SBATCH --partition=small
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=1
#SBATCH --time=02:00:00
#SBATCH --output=logs/emi-mpi8-%j.out
#SBATCH --error=logs/emi-mpi8-%j.err

set -euo pipefail

BASE=/data/numerik/people/fchegini/project/fatemeh/MicroCard/current_papers/BDDC_petsc_SDC_2023/gatevariables/20Nov2023/BDDC
IMAGE=localhost/kaskade7-mpi-build-environment:latest
ARCHIVE=/data/numerik/people/fchegini/kaskade7-mpi-build-environment.tar

cd "$BASE/emiModelBddcSdc_new_scratch"
mkdir -p logs

if ! podman image exists "$IMAGE"; then
  podman load -i "$ARCHIVE"
fi

podman run --rm --pull=never \
  -e PATH=/usr/lib64/mpi/gcc/openmpi4/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin \
  -e LD_LIBRARY_PATH=/usr/lib64/mpi/gcc/openmpi4/lib64:/usr/lib64/mpi/gcc/openmpi4/lib \
  -e OMPI_ALLOW_RUN_AS_ROOT=1 \
  -e OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1 \
  -v /data/numerik/people/fchegini/kaskade7_new:/kaskade7 \
  -v "$BASE:/BDDC" \
  -w /BDDC/emiModelBddcSdc_new_scratch \
  "$IMAGE" \
  bash -lc 'MPI_RANKS=8 MPI_ALLOW_ROOT=1 THREADS=1 BITS=22 LIFTING_MODE=1 EXCLUDE_CORNERS=1 VTK=0 ./run_case_mpi.sh'
