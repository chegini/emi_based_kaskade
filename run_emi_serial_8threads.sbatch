#!/bin/bash
#SBATCH --job-name=emi-serial8
#SBATCH --partition=small
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=02:00:00
#SBATCH --output=logs/emi-serial8-%j.out
#SBATCH --error=logs/emi-serial8-%j.err

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
  -e PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin \
  -e OMP_NUM_THREADS=8 \
  -v /data/numerik/people/fchegini/kaskade7_new:/kaskade7 \
  -v "$BASE:/BDDC" \
  -w /BDDC/emiModelBddcSdc_new_scratch \
  "$IMAGE" \
  bash -lc 'THREADS=8 BITS=22 LIFTING_MODE=1 EXCLUDE_CORNERS=1 VTK=0 ./run_case.sh'
