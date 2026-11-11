#!/bin/bash
#SBATCH --partition=big
#SBATCH --constraint=Gold6342
#SBATCH --nodes=8
#SBATCH --ntasks=16
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --job-name=pepe-fcg-100steps
#SBATCH --output=output/emi-pepe-100steps-%j.log
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?}"
export LANG=C LC_ALL=C PYTHONDONTWRITEBYTECODE=1
emi_sif=${emi_sif:-/data/numerik/people/fchegini/containers/kaskade-robin-20261105/kaskade7-robin.sif}
emi_kaskade=${emi_kaskade:-/data/numerik/people/fchegini/kaskade7_new}
emi_root="$PWD/output/emi-pepe-100steps-$SLURM_JOB_ID"
mkdir "$emi_root"
# Separate build keeps the existing emiModel executable intact.
srun --mpi=none --nodes=1 --ntasks=1 --ntasks-per-node=1 \
  apptainer exec --cleanenv --bind "$emi_kaskade:/kaskade7" \
  --bind "$PWD:$PWD" --pwd "$PWD" "$emi_sif" \
  /bin/bash tests/build_emi_coarse_batch_cluster.sh "$emi_root/build"
python3 -u tests/run_emi_pepe_trajectory.py \
  --case pepe --solvers fcg --split on --bits 16 --transform-order before \
  --ranks 16 --expected-nodes 8 --iterations 1000 --case-minutes 150 \
  --executable "$emi_root/build/emiModel" \
  --slurm-container "$emi_sif" --kaskade "$emi_kaskade" \
  --output "$emi_root/trajectory"
