#!/usr/bin/env bash
# One eight-node allocation; run the three layouts sequentially on this node pool.
#SBATCH --partition=big
#SBATCH --nodes=8
#SBATCH --ntasks=16
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --mem=64G
#SBATCH --time=03:00:00
#SBATCH --job-name=lego-weak-timings
#SBATCH --output=output/emi-weak-timings-%j.log
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?Submit from the EMI repository root}"
export LANG=C LC_ALL=C LC_CTYPE=C PYTHONDONTWRITEBYTECODE=1
emi_sif=${emi_sif:-/data/numerik/people/fchegini/containers/kaskade-robin-20261105/kaskade7-robin.sif}
emi_kaskade=${emi_kaskade:-/data/numerik/people/fchegini/kaskade7_new}
emi_cpu_model=${emi_cpu_model:-Intel(R) Xeon(R) Gold 6342 CPU @ 2.80GHz}
if (( $# != 3 )); then
  echo 'Usage: sbatch [node selection] tests/submit_emi_weak_timings.sbatch REF_2_NODES REF_4_NODES REF_8_NODES' >&2
  exit 2
fi
printf 'Timing job %s on %s; CPU model: %s\n' "$SLURM_JOB_ID" "$SLURM_JOB_NODELIST" "$emi_cpu_model"
for emi_reference in "$@"; do
  python3 -u tests/run_emi_weak_timings.py \
    --reference "$emi_reference" --slurm-container "$emi_sif" --kaskade "$emi_kaskade" \
    --cpu-model "$emi_cpu_model" --warmups 1 --repetitions 3
done
