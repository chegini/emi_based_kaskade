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
#SBATCH --job-name=emi-table5-fcg
#SBATCH --output=output/emi-table5-fcg-%j.log
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?Submit from the emi repository}"
export LANG=C LC_ALL=C PYTHONDONTWRITEBYTECODE=1
bench_sif=${bench_sif:-/data/numerik/people/fchegini/containers/kaskade-robin-20261105/kaskade7-robin.sif}
bench_kaskade=${bench_kaskade:-/data/numerik/people/fchegini/kaskade7_new}
bench_root="$PWD/output/emi-table5-fcg-$SLURM_JOB_ID"
mkdir "$bench_root"
scontrol show job "$SLURM_JOB_ID" > "$bench_root/allocation.txt"
sha256sum ./*.cpp ./*.hh Makefile tests/submit_emi_table5_fcg.sbatch tests/table5_fcg_previous.csv tests/compare_emi_table5_fcg.py \
 "$bench_kaskade"/mg/bddc.hh "$bench_kaskade"/mg/bddc.hpp \
 "$bench_kaskade"/mg/bddcSpaceTransfer.hh "$bench_kaskade"/mg/bddcSpaceTransfer.hpp \
 "$bench_kaskade"/mg/bddcPacketPacking.hh "$bench_kaskade"/tools/bddc_packing_benchmark.py \
 > "$bench_root/source-sha256.txt"
srun --mpi=none --nodes=1 --ntasks=1 --ntasks-per-node=1 --cpus-per-task=1 \
 apptainer exec --cleanenv --bind "$bench_kaskade:/kaskade7" \
 --bind "$(dirname "$PWD"):$(dirname "$PWD")" --pwd "$PWD" "$bench_sif" \
 /bin/bash tests/build_emi_coarse_batch_cluster.sh "$bench_root/build"
sha256sum -c "$bench_root/source-sha256.txt"
cp tests/table5_fcg_previous.csv "$bench_root/previous.csv"
failed=0
for order in before after; do
 python3 -u tests/run_emi_solver_study.py \
  --repo "$PWD" --executable "$bench_root/build/emiModel" \
  --case pepe --solvers fcg --bits 8 12 16 24 32 --split on \
  --master-worker 0 --lossless-packing 0 --transform-order "$order" \
  --iterations 1000 --case-minutes 40 \
  --ranks "$SLURM_NTASKS" --expected-nodes "$SLURM_JOB_NUM_NODES" \
  --slurm-container "$bench_sif" --kaskade "$bench_kaskade" \
  --output "$bench_root/$order" || failed=1
done
python3 tests/compare_emi_table5_fcg.py "$bench_root"
sha256sum -c "$bench_root/source-sha256.txt"
echo "Results: $bench_root/comparison.csv"

exit "$failed"
