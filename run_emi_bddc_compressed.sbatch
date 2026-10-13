#!/bin/bash
#SBATCH --job-name=emi-bddc-comp
#SBATCH --output=/data/numerik/people/fchegini/project/fatemeh/MicroCard/current_papers/BDDC_petsc_SDC_2023/gatevariables/20Nov2023/BDDC/emiModelBddcSdc_new_scratch/logs/session13/slurm-%j.out
#SBATCH --error=/data/numerik/people/fchegini/project/fatemeh/MicroCard/current_papers/BDDC_petsc_SDC_2023/gatevariables/20Nov2023/BDDC/emiModelBddcSdc_new_scratch/logs/session13/slurm-%j.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G

set -euo pipefail

BASE=/data/numerik/people/fchegini/project/fatemeh/MicroCard/current_papers/BDDC_petsc_SDC_2023/gatevariables/20Nov2023/BDDC
WORK="$BASE/emiModelBddcSdc_new_scratch"

IMAGE=localhost/kaskade7-mpi-build-environment:latest
IMAGE_ARCHIVE=/data/numerik/people/fchegini/kaskade7-mpi-build-environment.tar

RUN_HOST="output/session13/sbatch-bddc-compressed-${SLURM_JOB_ID}"
LOG_HOST="logs/session13/sbatch-bddc-compressed-${SLURM_JOB_ID}.log"

cd "$WORK"

mkdir -p "$RUN_HOST" "$(dirname "$LOG_HOST")"

echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
echo "Start: $(date)"

if ! podman image exists "$IMAGE"; then
  echo "Container image missing on compute node. Loading archive..."
  podman load -i "$IMAGE_ARCHIVE"
fi

podman image exists "$IMAGE"

podman run --rm \
  --pull=never \
  -e PATH=/usr/lib64/mpi/gcc/openmpi4/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin \
  -e LD_LIBRARY_PATH=/usr/lib64/mpi/gcc/openmpi4/lib64:/usr/lib64/mpi/gcc/openmpi4/lib \
  -v /data/numerik/people/fchegini/kaskade7_new:/kaskade7 \
  -v "$BASE:/BDDC" \
  -w /BDDC/emiModelBddcSdc_new_scratch \
  "$IMAGE" \
  ./emiModel \
    --mpi 0 \
    --input ../input_emi_mesh/pepe_combi_domi.vtu \
    --extra_set ../input_emi_mesh/pepe_combi_domi_extracellular.txt \
    --intra_set ../input_emi_mesh/pepe_combi_domi_intracellular.txt \
    --excited ../input_emi_mesh/pepe_combi_domi_excited.txt \
    --bddc 1 \
    --bddcCompression 1 \
    --bddcCompressionBits 16 \
    --bddcGraphLifting 1 \
    --bddcZigzag 1 \
    --bddcHuffman 1 \
    --bddcBitlength 1 \
    --bddcCompressionReport 1 \
    --bddcHistogram 1 \
    --bddcTolerance 1e-7 \
    --bddcIterations 20 \
    --bddcVerbose 1 \
    --finalTime 0.15 \
    --dir "/BDDC/emiModelBddcSdc_new_scratch/$RUN_HOST" \
  2>&1 | tee "$LOG_HOST"

echo "Finished: $(date)"
echo "Output: $WORK/$RUN_HOST"
echo "Log: $WORK/$LOG_HOST"
