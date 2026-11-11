#!/usr/bin/env bash
#SBATCH --partition=big
#SBATCH --constraint=Gold6342
#SBATCH --nodes=2
#SBATCH --ntasks=4
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --time=03:00:00
#SBATCH --job-name=emi-master-check
#SBATCH --output=output/emi-master-worker-%j.log
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?Submit from the EMI repository root}"
export LANG=C LC_ALL=C LC_CTYPE=C PYTHONDONTWRITEBYTECODE=1
emi_sif=${emi_sif:-/data/numerik/people/fchegini/containers/kaskade-robin-20261105/kaskade7-robin.sif}
emi_kaskade=${emi_kaskade:-/data/numerik/people/fchegini/kaskade7_new}
root="$PWD/output/emi-master-worker-${SLURM_JOB_ID}"
mkdir -p "$root"
for mode in 0 1; do
  python3 -u tests/run_emi_solver_study.py \
    --case weak_32_h16 --solvers fcg --bits 16 --split off \
    --transform-order before --master-worker "$mode" \
    --ranks "$SLURM_NTASKS" --expected-nodes "$SLURM_JOB_NUM_NODES" \
    --iterations 1000 --case-minutes 20 \
    --slurm-container "$emi_sif" --kaskade "$emi_kaskade" \
    --output "$root/mode-$mode"
done
printf 'Validation complete. Compare %s/mode-{0,1}/results.csv\n' "$root"
