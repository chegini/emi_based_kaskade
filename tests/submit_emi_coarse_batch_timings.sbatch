#!/usr/bin/env bash
#SBATCH --partition=big
#SBATCH --constraint=Gold6342
#SBATCH --nodes=8
#SBATCH --ntasks=16
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --mem=64G
#SBATCH --time=03:00:00
#SBATCH --job-name=lego-coarse-timings
#SBATCH --output=output/emi-weak-timings-%j.log
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?Submit from the EMI repository root}"
export LANG=C LC_ALL=C LC_CTYPE=C PYTHONDONTWRITEBYTECODE=1
emi_check_job=${1:?Supply the completed coarse-batch validation job ID}
[[ "$emi_check_job" =~ ^[0-9]+$ ]] || exit 2
emi_root="$PWD/output/emi-coarse-batch-$emi_check_job"
test -f "$emi_root/comparison.json"
test -f "$emi_root/VALIDATION_COMPLETE"
exec /bin/bash tests/submit_emi_weak_timings.sbatch \
  "$emi_root/validation-weak_32_h16" \
  "$emi_root/validation-weak_64_h16" \
  "$emi_root/validation-weak_128_h16"
