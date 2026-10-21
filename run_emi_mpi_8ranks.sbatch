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
MPI_RANKS=${MPI_RANKS:-8}
THREADS=${THREADS:-1}
BITS=${BITS:-22}
LIFTING_MODE=${LIFTING_MODE:-1}
EXCLUDE_CORNERS=${EXCLUDE_CORNERS:-1}
MESH=${MESH:-../input_emi_mesh_blocks/SX_h16_original_mesh.vtu}
EXTRA=${EXTRA:-../input_emi_mesh_blocks/SX_h16_original_list_extracellular.txt}
INTRA=${INTRA:-../input_emi_mesh_blocks/SX_h16_original_list_intracellular.txt}
EXCITED=${EXCITED:-../input_emi_mesh_blocks/SX_h16_original_early_excited.txt}
WRAPPER_CPUS=${WRAPPER_CPUS:-$MPI_RANKS}

cd "$BASE/emiModelBddcSdc_new_scratch"
mkdir -p logs

echo "Job ID: $SLURM_JOB_ID"
echo "Host: $(hostname)"
echo "Ranks: $MPI_RANKS"
echo "Threads per rank: $THREADS"
echo "Bits: $BITS"
echo "CPU binding: cores"
echo "Wrapper CPUs: $WRAPPER_CPUS"

if ! podman image exists "$IMAGE"; then
  podman load -i "$ARCHIVE"
fi

srun --ntasks=1 --cpus-per-task="$WRAPPER_CPUS" --cpu-bind=cores --hint=nomultithread podman run --rm --pull=never --pids-limit=4096 --ulimit nproc=4096:4096 \
  -e PATH=$HOST_MPI/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin \
  -e LD_LIBRARY_PATH=/kaskade7/libs:$HOST_MPI/lib:$HOST_MPI/lib64:/usr/lib64:/usr/lib:/usr/local/lib64:/usr/local/lib \
  -e MPI_RANKS="$MPI_RANKS" \
  -e THREADS="$THREADS" \
  -e BITS="$BITS" \
  -e LIFTING_MODE="$LIFTING_MODE" \
  -e EXCLUDE_CORNERS="$EXCLUDE_CORNERS" \
  -e ZIGZAG="$ZIGZAG" \
  -e HUFFMAN="$HUFFMAN" \
  -e BITLENGTH="$BITLENGTH" \
  -e HISTOGRAM="$HISTOGRAM" \
  -e MESH="$MESH" \
  -e EXTRA="$EXTRA" \
  -e INTRA="$INTRA" \
  -e EXCITED="$EXCITED" \
  -e MPI_BIND_ARGS="-bind-to core" \
  -e OMP_PROC_BIND=close \
  -e OMP_PLACES=cores \
  -e VTK="${VTK:-0}" \
  -e FINAL_TIME="${FINAL_TIME:-0.1}" \
  -e MIN_INTERFACE_SIZE="${MIN_INTERFACE_SIZE:-0}" \
  -v "$HOST_MPI:$HOST_MPI:ro" \
  -v /data/numerik/people/fchegini/kaskade7_new:/kaskade7 \
  -v "$BASE:/BDDC" \
  -w /BDDC/emiModelBddcSdc_new_scratch \
  "$IMAGE" \
  bash -lc './run_case_mpi.sh'

if command -v sacct >/dev/null 2>&1; then
  sacct -j "$SLURM_JOB_ID" --format=JobID,State,Elapsed,TotalCPU,AllocCPUS,MaxRSS
fi
