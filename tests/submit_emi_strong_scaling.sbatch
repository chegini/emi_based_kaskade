#!/usr/bin/env bash
#SBATCH --partition=big
#SBATCH --constraint=Gold6342
#SBATCH --nodes=16
#SBATCH --ntasks=32
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --job-name=lego-strong
#SBATCH --output=output/emi-strong-scaling-%j.log

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?Submit from the EMI repository root}"
export LANG=C LC_ALL=C LC_CTYPE=C

emi_sif=/data/numerik/people/fchegini/containers/kaskade-robin-20261105/kaskade7-robin.sif
emi_kaskade=/data/numerik/people/fchegini/kaskade7_new
emi_executable="$PWD/output/emi-coarse-batch-6009589/build/emiModel"

test -x "$emi_executable"

for layout in "2 4" "4 8" "8 16" "16 32"; do
  read -r nodes ranks <<< "$layout"
  echo "=== strong scaling: weak_128_h16, ${nodes} nodes, ${ranks} ranks ==="

  python3 -u tests/run_emi_solver_study.py \
    --executable "$emi_executable" \
    --ranks "$ranks" \
    --expected-nodes "$nodes" \
    --slurm-container "$emi_sif" \
    --kaskade "$emi_kaskade" \
    --case weak_128_h16 \
    --solvers fcg \
    --split on \
    --bits 12 16 \
    --transform-order before \
    --iterations 1000 \
    --case-minutes 60
done
