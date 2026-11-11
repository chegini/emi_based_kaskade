#!/usr/bin/env bash
#SBATCH --partition=big
#SBATCH --constraint=Gold6342
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --job-name=pepe-physical
#SBATCH --output=output/emi-physical-histogram-%j.log
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?Submit from emiModelBddcSdc_new_scratch}"
export LANG=C LC_ALL=C
emi_sif=${emi_sif:-/data/numerik/people/fchegini/containers/kaskade-robin-20261105/kaskade7-robin.sif}
emi_kaskade=${emi_kaskade:-/data/numerik/people/fchegini/kaskade7_new}
# One container, two local MPI ranks. FCG replaces the original histogram study's FGMRES.
# Rebuilds this checkout's emiModel: do not run another build in this checkout concurrently.
srun --mpi=none --nodes=1 --ntasks=1 --cpus-per-task=2 \
  apptainer exec --cleanenv \
  --bind "$emi_kaskade:/kaskade7" \
  --bind "$(dirname "$PWD"):/BDDC" \
  --pwd /BDDC/emiModelBddcSdc_new_scratch "$emi_sif" \
  /bin/bash tests/run_emi_physical_histogram_container.sh "$SLURM_JOB_ID" "${1:-}" </dev/null
