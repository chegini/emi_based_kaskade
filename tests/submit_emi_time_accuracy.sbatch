#!/bin/bash
#SBATCH --partition=big
#SBATCH --nodes=2
#SBATCH --ntasks=4
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=04:00:00
#SBATCH --job-name=pepe-time-accuracy
#SBATCH --output=output/emi-time-accuracy-%j.log

# Submit from the EMI repository root. Override nodes/tasks together with sbatch.
# Optional argument --with-timings adds the 36-solve timing campaign.
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?Submit from the EMI repository root}"
export LANG=C LC_ALL=C LC_CTYPE=C
emi_sif=${emi_sif:-/data/numerik/people/fchegini/containers/kaskade-robin-20261105/kaskade7-robin.sif}
emi_kaskade=${emi_kaskade:-/data/numerik/people/fchegini/kaskade7_new}
emi_nodes=${SLURM_JOB_NUM_NODES:?}
emi_ranks=${SLURM_NTASKS:?}
case "$emi_nodes" in 2|4|8) ;; *) echo 'Require 2, 4 or 8 nodes' >&2; exit 2 ;; esac
if (( emi_ranks != 2 * emi_nodes )); then
  echo 'Require exactly two ranks per node; override --nodes and --ntasks together' >&2
  exit 2
fi
if (( $# > 1 )) || [[ ${1:-} != '' && ${1:-} != '--with-timings' ]]; then
  echo 'Usage: sbatch [resource options] tests/submit_emi_time_accuracy.sbatch [--with-timings]' >&2
  exit 2
fi
emi_campaign="$PWD/output/emi-time-campaign-n${emi_nodes}-${SLURM_JOB_ID}"
mkdir "$emi_campaign"
printf 'Job %s: %s nodes, %s ranks; hosts %s\n' "$SLURM_JOB_ID" "$emi_nodes" "$emi_ranks" "$SLURM_JOB_NODELIST"

# Fresh validation in the same allocation supplies both the eight-node gate
# and the exact first-step reference for each longer trajectory.
python3 -u tests/run_emi_mpi_routing_check.py \
  --case pepe --ranks "$emi_ranks" --expected-nodes "$emi_nodes" \
  --slurm-container "$emi_sif" --kaskade "$emi_kaskade" \
  --slurm-btl self,vader,tcp --slurm-vader-single-copy none \
  --check-local-bypass --bits 16 --transform-order before \
  --output "$emi_campaign/one-step-reference"

python3 -u tests/run_emi_time_accuracy_check.py \
  --reference "$emi_campaign/one-step-reference" \
  --slurm-container "$emi_sif" --kaskade "$emi_kaskade" \
  --steps 10 --output "$emi_campaign/accuracy-10steps"

if [[ ${1:-} == '--with-timings' ]]; then
  python3 -u tests/run_emi_routing_timings.py \
    --reference "$emi_campaign/one-step-reference" \
    --slurm-container "$emi_sif" --kaskade "$emi_kaskade" \
    --warmups 1 --repetitions 3 --output "$emi_campaign/timings"
fi
printf 'Campaign complete: %s\n' "$emi_campaign"
