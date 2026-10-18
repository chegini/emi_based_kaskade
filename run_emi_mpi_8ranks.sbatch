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
HOST_MPI=/home/htc/fchegini/opencarp_dependencies/petsc_install

cd "$BASE/emiModelBddcSdc_new_scratch"
mkdir -p logs

if ! podman image exists "$IMAGE"; then
  podman load -i "$ARCHIVE"
fi

podman run --rm --pull=never \
  -e PATH=$HOST_MPI/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin \
  -e LD_LIBRARY_PATH=/kaskade7/libs:$HOST_MPI/lib:$HOST_MPI/lib64:/usr/lib64:/usr/lib:/usr/local/lib64:/usr/local/lib \
  -e MPI_RANKS=8 \
  -e THREADS=1 \
  -e BITS=22 \
  -e LIFTING_MODE=1 \
  -e EXCLUDE_CORNERS=1 \
  -e VTK=0 \
  -v "$HOST_MPI:$HOST_MPI:ro" \
  -v /data/numerik/people/fchegini/kaskade7_new:/kaskade7 \
  -v "$BASE:/BDDC" \
  -w /BDDC/emiModelBddcSdc_new_scratch \
  "$IMAGE" \
  bash -lc './run_case_mpi.sh'
