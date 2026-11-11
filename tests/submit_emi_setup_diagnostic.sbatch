#!/usr/bin/env bash
# Submit from the EMI repository root; rebuild and run one raw FCG solve.
#SBATCH --partition=big
#SBATCH --nodes=8
#SBATCH --ntasks=16
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --job-name=emi-setup-check
#SBATCH --output=output/emi-setup-check-%j.log
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?Submit from the EMI repository root}"
export LANG=C LC_ALL=C LC_CTYPE=C PYTHONDONTWRITEBYTECODE=1
emi_sif=${emi_sif:-/data/numerik/people/fchegini/containers/kaskade-robin-20261105/kaskade7-robin.sif}
emi_kaskade=${emi_kaskade:-/data/numerik/people/fchegini/kaskade7_new}
printf 'Setup diagnostic job %s on %s\n' "$SLURM_JOB_ID" "$SLURM_JOB_NODELIST"
srun --mpi=none --nodes=1 --ntasks=1 --ntasks-per-node=1 \
  apptainer exec --cleanenv \
  --bind "$emi_kaskade:/kaskade7" \
  --bind "$(dirname "$PWD"):/BDDC" \
  --pwd "/BDDC/$(basename "$PWD")" "$emi_sif" \
  /bin/bash tests/build_emi_mpi_routing_cluster.sh
exec python3 -u tests/run_emi_solver_study.py \
  --ranks "$SLURM_NTASKS" --expected-nodes "$SLURM_JOB_NUM_NODES" \
  --slurm-container "$emi_sif" --kaskade "$emi_kaskade" \
  --case kermit --solvers fcg --case-minutes 60 \
  --setup-diagnostic-only "$@"
