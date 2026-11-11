#!/usr/bin/env bash
#SBATCH --partition=big
#SBATCH --constraint=Gold6342
#SBATCH --nodes=2
#SBATCH --ntasks=4
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --job-name=emi-pack-solve
#SBATCH --output=output/emi-packing-solve-%j.log
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?Submit from the EMI repository root}"
export LANG=C LC_ALL=C PYTHONDONTWRITEBYTECODE=1
emi_sif=${emi_sif:-/data/numerik/people/fchegini/containers/kaskade-robin-20261105/kaskade7-robin.sif}
emi_kaskade=${emi_kaskade:-/data/numerik/people/fchegini/kaskade7_new}
emi_root="$PWD/output/emi-packing-solve-$SLURM_JOB_ID"
mkdir "$emi_root"
# Freeze hashes before compilation; the builder copies application sources.
sha256sum EMI_model.cpp ./*.hh tests/run_emi_lossless_fullsolve.py \
  "$emi_kaskade"/mg/bddc.hh "$emi_kaskade"/mg/bddc.hpp \
  "$emi_kaskade"/mg/bddcSpaceTransfer.hh "$emi_kaskade"/mg/bddcSpaceTransfer.hpp \
  "$emi_kaskade"/mg/bddcPacketPacking.hh > "$emi_root/source-sha256.txt"
srun --mpi=none --nodes=1 --ntasks=1 --ntasks-per-node=1 --cpus-per-task=1 \
  apptainer exec --cleanenv --bind "$emi_kaskade:/kaskade7" \
  --bind "$PWD:$PWD" --pwd "$PWD" "$emi_sif" \
  /bin/bash tests/build_emi_coarse_batch_cluster.sh "$emi_root/build"
sha256sum -c "$emi_root/source-sha256.txt"
python3 -u tests/run_emi_lossless_fullsolve.py \
  --executable "$emi_root/build/emiModel" --output "$emi_root/validation" \
  --kaskade "$emi_kaskade" --slurm-container "$emi_sif" \
  --ranks "$SLURM_NTASKS" --expected-nodes "$SLURM_JOB_NUM_NODES" "$@"
sha256sum -c "$emi_root/source-sha256.txt"
echo "Results: $emi_root/validation"
