#!/bin/bash
#SBATCH --job-name=emi-2n8r
#SBATCH --partition=big
#SBATCH --nodes=2
#SBATCH --ntasks=8
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=1
#SBATCH --time=04:00:00
#SBATCH --output=logs/emi-2nodes-8ranks-%j.out
#SBATCH --error=logs/emi-2nodes-8ranks-%j.err

# Multi-node launch: Slurm starts one container per MPI task. Do not call
# run_case_mpi.sh here because it would start a nested mpirun on one node.
# The application is built against the cluster MPICH installation. Use PMI2,
# which is the Slurm integration supported by this MPICH Hydra runtime.
set -euo pipefail

BASE=/data/numerik/people/fchegini/project/fatemeh/MicroCard/current_papers/BDDC_petsc_SDC_2023/gatevariables/20Nov2023/BDDC
KASKADE7=/data/numerik/people/fchegini/kaskade7_new
IMAGE=localhost/kaskade7-mpi-build-environment:latest
ARCHIVE=/data/numerik/people/fchegini/kaskade7-mpi-build-environment.tar
MPI_IMPL=${MPI_IMPL:-pmi2}
HOST_MPI=/home/htc/fchegini/opencarp_dependencies/petsc_install
MPI_RANKS=${MPI_RANKS:-${SLURM_NTASKS}}
THREADS=${THREADS:-1}
CPUS_PER_RANK=${CPUS_PER_RANK:-1}
BITS=${BITS:-22}
EXCLUDE_CORNERS=${EXCLUDE_CORNERS:-1}
ZIGZAG=${ZIGZAG:-1}
HUFFMAN=${HUFFMAN:-1}
BITLENGTH=${BITLENGTH:-1}
HISTOGRAM=${HISTOGRAM:-0}
FINAL_TIME=${FINAL_TIME:-0.1}
REFINE=${REFINE:-0}
GRAPH_LIFTING=${GRAPH_LIFTING:-1}
GRAPH_LIFTING_AFTER_QUANTIZATION=${GRAPH_LIFTING_AFTER_QUANTIZATION:-1}
BDD_ITERATIONS=${BDD_ITERATIONS:-300}
MIN_INTERFACE_SIZE=${MIN_INTERFACE_SIZE:-0}
THRESHOLD_INTERFACE_TYPES=${THRESHOLD_INTERFACE_TYPES:-2}
INTERFACE_TYPES=${INTERFACE_TYPES:-7}
MESH=${MESH:-../input_emi_mesh_blocks/SX_h16_original_mesh.vtu}
EXTRA=${EXTRA:-../input_emi_mesh_blocks/SX_h16_original_list_extracellular.txt}
INTRA=${INTRA:-../input_emi_mesh_blocks/SX_h16_original_list_intracellular.txt}
EXCITED=${EXCITED:-../input_emi_mesh_blocks/SX_h16_original_early_excited.txt}
RUN_NAME=${RUN_NAME:-multi-node-ranks-${MPI_RANKS}-bits-${BITS}-job-${SLURM_JOB_ID}}
RUN="output/multi-node/${RUN_NAME}"
export BASE KASKADE7 IMAGE ARCHIVE HOST_MPI RUN
export MPI_IMPL BITS EXCLUDE_CORNERS THREADS CPUS_PER_RANK MESH EXTRA INTRA EXCITED
export ZIGZAG HUFFMAN BITLENGTH HISTOGRAM FINAL_TIME
export REFINE GRAPH_LIFTING GRAPH_LIFTING_AFTER_QUANTIZATION BDD_ITERATIONS
export MIN_INTERFACE_SIZE THRESHOLD_INTERFACE_TYPES INTERFACE_TYPES

cd "$BASE/emiModelBddcSdc_new_scratch"
mkdir -p logs "$RUN"

echo "Job ID: $SLURM_JOB_ID"
echo "Nodes: $(scontrol show hostnames "$SLURM_JOB_NODELIST" | tr "\n" " ")"
echo "Tasks: $SLURM_NTASKS"
echo "MPI mode: $MPI_IMPL"
echo "Ranks per node: ${SLURM_NTASKS_PER_NODE:-unknown}"
echo "Bits: $BITS"
echo "Minimum interface size: ${MIN_INTERFACE_SIZE:-0}"
echo "Threshold interface types: ${THRESHOLD_INTERFACE_TYPES:-2}"
echo "Mesh: $MESH"
echo "CPU binding: cores"
echo "CPUs per rank: $CPUS_PER_RANK"

# The image is already installed in the normal Podman store on the ZIB
# system. Do not import the archive inside an sbatch job: the batch
# environment terminates podman load while it is unpacking layers.
if ! podman image exists "$IMAGE"; then
  echo "ERROR: image $IMAGE is not available in the normal Podman store on $(hostname)" >&2
  exit 2
fi

srun --mpi="$MPI_IMPL" --ntasks="$SLURM_NTASKS" \
  --ntasks-per-node="$SLURM_NTASKS_PER_NODE" \
  --cpus-per-task="$CPUS_PER_RANK" \
  --cpu-bind=cores --hint=nomultithread --distribution=block:block \
  --kill-on-bad-exit=0 \
  --output="logs/emi-multi-%j-rank%t.out" \
  --error="logs/emi-multi-%j-rank%t.err" \
  bash -lc '
    set -euxo pipefail
    echo "mpi-container rank=${SLURM_PROCID} node=$(hostname) image=$IMAGE"
    : "${PMI_FD:?PMI_FD was not provided by Slurm}"
    eval "exec 3<&${PMI_FD}"
    exec podman run --rm --cgroups=disabled --pull=never \
      --network=host \
      --ipc=host \
      --preserve-fds=1 \
      --env-host \
      -e PMI_FD=3 \
      -v "$HOST_MPI:$HOST_MPI:ro" \
      -e PATH=$HOST_MPI/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin \
      -e LD_LIBRARY_PATH=/kaskade7/libs:$HOST_MPI/lib:$HOST_MPI/lib64:/usr/local/lib64:/usr/local/lib \
      -e OMP_NUM_THREADS="$THREADS" \
      -e MIN_INTERFACE_SIZE="${MIN_INTERFACE_SIZE:-0}" \
      -e THRESHOLD_INTERFACE_TYPES="${THRESHOLD_INTERFACE_TYPES:-2}" \
      -e OMP_PROC_BIND=close \
      -e OMP_PLACES=cores \
      -v "$KASKADE7:/kaskade7" \
      -v "$BASE:/BDDC" \
      -w /BDDC/emiModelBddcSdc_new_scratch \
      "$IMAGE" \
      ./emiModel \
        --mpi 1 \
        --input "$MESH" \
        --extra_set "$EXTRA" \
        --intra_set "$INTRA" \
        --excited "$EXCITED" \
        --refine "${REFINE:-0}" \
        --bddc 1 \
        --bddcCompression 1 \
        --bddcCompressionBits "$BITS" \
        --bddcCompressionMinInterfaceSize "$MIN_INTERFACE_SIZE" \
        --bddcCompressionThresholdInterfaceTypes "$THRESHOLD_INTERFACE_TYPES" \
        --bddcGraphLifting "$GRAPH_LIFTING" \
        --bddcGraphLiftingAfterQuantization "$GRAPH_LIFTING_AFTER_QUANTIZATION" \
        --bddcZigzag "$ZIGZAG" \
        --bddcHuffman "$HUFFMAN" \
        --bddcBitlength "$BITLENGTH" \
        --bddcExcludeCornersFromCompression "$EXCLUDE_CORNERS" \
        --bddcCompressionReport 1 \
        --bddcHistogram "$HISTOGRAM" \
        --bddcTolerance 1e-7 \
        --bddcIterations "$BDD_ITERATIONS" \
        --bddcVerbose 1 \
        --nThreads "$THREADS" \
        --dt 0.01 \
        --finalTime "$FINAL_TIME" \
        --vtk 0 \
        --dir "/BDDC/emiModelBddcSdc_new_scratch/$RUN"
  '

if command -v sacct >/dev/null 2>&1; then
  sacct -j "$SLURM_JOB_ID" --format=JobID,State,Elapsed,TotalCPU,AllocCPUS,MaxRSS
fi
