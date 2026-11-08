#!/bin/bash
#SBATCH --partition=big
#SBATCH --nodes=2
#SBATCH --ntasks=4
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=03:00:00
#SBATCH --job-name=pepe-routing-timings
#SBATCH --output=output/emi-routing-timings-%j.log

# Submit from the EMI repository root. Slurm copies this script to its spool,
# so resolve the project through SLURM_SUBMIT_DIR, not BASH_SOURCE.
# For a four-node/eight-rank reference, submit with:
# sbatch --nodes=4 --ntasks=8 --ntasks-per-node=2 tests/submit_emi_routing_timings.sbatch REFERENCE
# The Python runner inherits the layout from REFERENCE and checks the allocation.
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?Submit with sbatch from the EMI repository root}"
test -f tests/run_emi_routing_timings.py
export LANG=C LC_ALL=C LC_CTYPE=C
emi_sif=${emi_sif:-/data/numerik/people/fchegini/containers/kaskade-robin-20261105/kaskade7-robin.sif}
emi_kaskade=${emi_kaskade:-/data/numerik/people/fchegini/kaskade7_new}
emi_reference=${1:-output/emi-mpi-routing-pepe-20261107T120921Z-2561102}
printf 'Timing job %s on %s\n' "$SLURM_JOB_ID" "$SLURM_JOB_NODELIST"
exec python3 -u tests/run_emi_routing_timings.py \
  --reference "$emi_reference" \
  --slurm-container "$emi_sif" \
  --kaskade "$emi_kaskade" \
  --warmups 1 --repetitions 3
