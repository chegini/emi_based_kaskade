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
HOST_MPI=/home/htc/fchegini/opencarp_dependencies/petsc_install
BUILD_TMP="$BASE/emiModelBddcSdc_new_scratch/build-tmp-${SLURM_JOB_ID:-manual}"
MPI_CXX="$HOST_MPI/bin/mpicxx"
# Use the real MPICH API and suppress MUMPS' sequential mpi.h shim.
MPI_INCLUDE_FLAGS="-I$HOST_MPI/include -include $HOST_MPI/include/mpi.h -DMUMPS_MPI_H"
MPI_LIBRARY_PATH="$HOST_MPI/lib:$HOST_MPI/lib64"

cd "$BASE/emiModelBddcSdc_new_scratch"
mkdir -p "$BUILD_TMP"

echo "Job ID: ${SLURM_JOB_ID:-unknown}"
echo "Node: $(hostname)"
echo "Start: $(date)"

if ! podman image exists "$IMAGE"; then
  podman load -i "$IMAGE_ARCHIVE"
fi

echo "MPI compiler: $HOST_MPI/bin/mpicxx"

podman run --rm --pull=never \
  -e PATH=$HOST_MPI/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin \
  -e LD_LIBRARY_PATH=$MPI_LIBRARY_PATH:/usr/local/lib64:/usr/local/lib \
  -e CPATH=$HOST_MPI/include \
  -e CPLUS_INCLUDE_PATH=$HOST_MPI/include \
  -e LIBRARY_PATH=$MPI_LIBRARY_PATH \
  -e HOST_MPI=$HOST_MPI \
  -e MPI_CXX=$MPI_CXX \
  -e MPI_INCLUDE_FLAGS="$MPI_INCLUDE_FLAGS" \
  -e TMPDIR=/BDDC/emiModelBddcSdc_new_scratch/$(basename "$BUILD_TMP") \
  -v "$HOST_MPI:$HOST_MPI:ro" \
  -v /data/numerik/people/fchegini/kaskade7_new:/kaskade7 \
  -v "$BASE:/BDDC" \
  -w /BDDC/emiModelBddcSdc_new_scratch \
  "$IMAGE" \
  bash -lc '
    set -euo pipefail
    export PATH="$HOST_MPI/bin:$PATH"
    export LD_LIBRARY_PATH="$HOST_MPI/lib:$HOST_MPI/lib64:/usr/local/lib64:/usr/local/lib:${LD_LIBRARY_PATH:-}"
    export LIBRARY_PATH="$HOST_MPI/lib:$HOST_MPI/lib64:${LIBRARY_PATH:-}"
    test -x "$MPI_CXX"
    test -f "$HOST_MPI/include/mpi.h"
    ldd /kaskade7/tools/remove_install_deps | grep -E "libmpi|not found" || true
    "$MPI_CXX" --version
    make clean KASKADE7=/kaskade7 MPI=1 \
      MPICXX="$MPI_CXX" \
      MPI_INCLUDE_FLAGS="$MPI_INCLUDE_FLAGS"
    make -j2 KASKADE7=/kaskade7 MPI=1 \
      MPICXX="$MPI_CXX" \
      MPI_INCLUDE_FLAGS="$MPI_INCLUDE_FLAGS"
    test -x emiModel
    ls -lh emiModel
  '

echo "Finished: $(date)"
