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
HOST_MPI=/home/htc/fchegini/opencarp_dependencies/petsc_install
THREADS=${THREADS:-8}
BITS=${BITS:-22}
EXCLUDE_CORNERS=${EXCLUDE_CORNERS:-1}
MESH=${MESH:-../input_emi_mesh_blocks/SX_h16_original_mesh.vtu}
EXTRA=${EXTRA:-../input_emi_mesh_blocks/SX_h16_original_list_extracellular.txt}
INTRA=${INTRA:-../input_emi_mesh_blocks/SX_h16_original_list_intracellular.txt}
EXCITED=${EXCITED:-../input_emi_mesh_blocks/SX_h16_original_early_excited.txt}

cd "$BASE/emiModelBddcSdc_new_scratch"
mkdir -p logs

if ! podman image exists "$IMAGE"; then
  podman load -i "$ARCHIVE"
fi

srun --ntasks=1 --cpu-bind=cores --hint=nomultithread podman run --rm --pull=never \
  -e PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin \
  -e LD_LIBRARY_PATH=/kaskade7/libs:$HOST_MPI/lib:$HOST_MPI/lib64:/usr/lib64:/usr/lib:/usr/local/lib64:/usr/local/lib \
  -e OMP_NUM_THREADS="$THREADS" \
  -e OMP_PROC_BIND=close \
  -e OMP_PLACES=cores \
  -e MESH="$MESH" \
  -e EXTRA="$EXTRA" \
  -e INTRA="$INTRA" \
  -e EXCITED="$EXCITED" \
  -v "$HOST_MPI:$HOST_MPI:ro" \
  -v /data/numerik/people/fchegini/kaskade7_new:/kaskade7 \
  -v "$BASE:/BDDC" \
  -w /BDDC/emiModelBddcSdc_new_scratch \
  "$IMAGE" \
  bash -lc "THREADS=$THREADS BITS=$BITS LIFTING_MODE=1 EXCLUDE_CORNERS=$EXCLUDE_CORNERS VTK=0 ./run_case.sh"

if command -v sacct >/dev/null 2>&1; then
  sacct -j "$SLURM_JOB_ID" --format=JobID,State,Elapsed,TotalCPU,AllocCPUS,MaxRSS
fi
