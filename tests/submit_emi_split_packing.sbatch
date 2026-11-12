#!/usr/bin/env bash
#SBATCH --partition=big
#SBATCH --constraint=Gold6342
#SBATCH --nodes=8
#SBATCH --ntasks=16
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --hint=nomultithread
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --job-name=emi-split-pack
#SBATCH --output=output/emi-split-packing-%j.log
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?Submit from the EMI repository root}"
export LANG=C LC_ALL=C PYTHONDONTWRITEBYTECODE=1
emi_sif=${emi_sif:-/data/numerik/people/fchegini/containers/kaskade-robin-20261105/kaskade7-robin.sif}
emi_kaskade=${emi_kaskade:-/data/numerik/people/fchegini/kaskade7_new}
emi_root="$PWD/output/emi-split-packing-$SLURM_JOB_ID"
mkdir "$emi_root"
scontrol show job "$SLURM_JOB_ID" > "$emi_root/allocation.txt"
# Use the checked-out sources directly. Build a new executable for this job.
# Do not edit these sources while the campaign runs.
sha256sum EMI_model.cpp ./*.hh Makefile tests/*.py tests/*.cpp tests/*.sh \
  tests/submit_emi_split_packing.sbatch tests/table5_fcg_previous.csv \
  "$emi_kaskade"/mg/bddc*.hh "$emi_kaskade"/mg/bddc*.hpp \
  > "$emi_root/source-sha256.txt"
srun --mpi=none --nodes=1 --ntasks=1 --ntasks-per-node=1 --cpus-per-task=1 \
  apptainer exec --cleanenv --bind "$emi_kaskade:/kaskade7" \
  --bind "$PWD:$PWD" --pwd "$PWD" "$emi_sif" \
  /bin/bash tests/build_emi_coarse_batch_cluster.sh "$emi_root/build"
sha256sum -c "$emi_root/source-sha256.txt"

# No precision sweep until matched 16-bit OFF/ON full solves preserve exact bytes.
for split in on off; do
  python3 -u tests/run_emi_lossless_fullsolve.py \
    --repo "$PWD" --executable "$emi_root/build/emiModel" \
    --kaskade "$emi_kaskade" --slurm-container "$emi_sif" \
    --case pepe --split "$split" --bits 16 --case-minutes 40 \
    --ranks "$SLURM_NTASKS" --expected-nodes "$SLURM_JOB_NUM_NODES" \
    --output "$emi_root/parity-split-$split"
done
failed=0
for split in on off; do
  mkdir "$emi_root/study-split-$split"
  for order in before after; do
    python3 -u tests/run_emi_solver_study.py \
      --repo "$PWD" --executable "$emi_root/build/emiModel" \
      --kaskade "$emi_kaskade" --slurm-container "$emi_sif" \
      --case pepe --solvers fcg --split "$split" --lossless-packing 1 \
      --bits 8 12 16 24 32 --transform-order "$order" --case-minutes 40 \
      --ranks "$SLURM_NTASKS" --expected-nodes "$SLURM_JOB_NUM_NODES" \
      --output "$emi_root/study-split-$split/$order" || failed=1
  done
done
cp tests/table5_fcg_previous.csv "$emi_root/study-split-on/previous.csv"
python3 tests/compare_emi_table5_fcg.py "$emi_root/study-split-on"
python3 tests/summarize_split_packing.py "$emi_root"
sha256sum -c "$emi_root/source-sha256.txt"
echo "Completed all requested cases; failure flag=$failed (inspect retained failed configurations)."
echo "Results: $emi_root"
exit "$failed"
