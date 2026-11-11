#!/usr/bin/env bash
# Submit from the EMI repository root. Each case is a separate bounded srun step.
#SBATCH --partition=big
#SBATCH --nodes=8
#SBATCH --ntasks=16
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=08:00:00
#SBATCH --job-name=pepe-solver-study
#SBATCH --output=output/emi-solver-study-%j.log
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?Submit from the EMI repository root}"
test -f tests/run_emi_solver_study.py
export LANG=C LC_ALL=C LC_CTYPE=C
export PYTHONDONTWRITEBYTECODE=1
emi_sif=${emi_sif:-/data/numerik/people/fchegini/containers/kaskade-robin-20261105/kaskade7-robin.sif}
emi_kaskade=${emi_kaskade:-/data/numerik/people/fchegini/kaskade7_new}
printf 'Solver study job %s: %s nodes, %s ranks on %s\n' \
  "$SLURM_JOB_ID" "$SLURM_JOB_NUM_NODES" "$SLURM_NTASKS" "$SLURM_JOB_NODELIST"
exec python3 -u tests/run_emi_solver_study.py \
  --ranks "$SLURM_NTASKS" --expected-nodes "$SLURM_JOB_NUM_NODES" \
  --slurm-container "$emi_sif" --kaskade "$emi_kaskade" "$@"
