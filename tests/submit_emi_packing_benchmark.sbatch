#!/usr/bin/env bash
#SBATCH --partition=big
#SBATCH --constraint=Gold6342
#SBATCH --nodes=2
#SBATCH --ntasks=4
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --hint=nomultithread
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --job-name=emi-pack-bench
#SBATCH --output=output/emi-packing-benchmark-%j.log
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?Submit from the emi repository}"
export LANG=C LC_ALL=C PYTHONDONTWRITEBYTECODE=1
bench_sif=${bench_sif:-/data/numerik/people/fchegini/containers/kaskade-robin-20261105/kaskade7-robin.sif}
bench_kaskade=${bench_kaskade:-/data/numerik/people/fchegini/kaskade7_new}
bench_root="$PWD/output/emi-packing-benchmark-$SLURM_JOB_ID"
mkdir "$bench_root"
scontrol show job "$SLURM_JOB_ID" > "$bench_root/allocation.txt"
sha256sum ./*.cpp ./*.hh Makefile tests/*packing_benchmark* \
 "$bench_kaskade"/mg/bddc.hh "$bench_kaskade"/mg/bddc.hpp \
 "$bench_kaskade"/mg/bddcSpaceTransfer.hh "$bench_kaskade"/mg/bddcSpaceTransfer.hpp \
 "$bench_kaskade"/mg/bddcPacketPacking.hh "$bench_kaskade"/tools/bddc_packing_benchmark.py \
 > "$bench_root/source-sha256.txt"
srun --mpi=none --nodes=1 --ntasks=1 --ntasks-per-node=1 --cpus-per-task=1 \
 apptainer exec --cleanenv --bind "$bench_kaskade:/kaskade7" \
 --bind "$(dirname "$PWD"):$(dirname "$PWD")" --pwd "$PWD" "$bench_sif" \
 /bin/bash tests/build_emi_coarse_batch_cluster.sh "$bench_root/build"
sha256sum -c "$bench_root/source-sha256.txt"
python3 -u "$bench_kaskade/tools/bddc_packing_benchmark.py" \
 --model emi --repo "$PWD" --executable "$bench_root/build/emiModel" \
 --output "$bench_root/results" --kaskade "$bench_kaskade" --container "$bench_sif" \
 --nodes "$SLURM_JOB_NUM_NODES" --ranks "$SLURM_NTASKS" "$@"
sha256sum -c "$bench_root/source-sha256.txt"
echo "Results: $bench_root/results"
