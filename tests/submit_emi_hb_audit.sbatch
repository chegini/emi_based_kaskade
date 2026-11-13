#!/usr/bin/env bash
#SBATCH --partition=big
#SBATCH --constraint=Gold6342
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --job-name=emi-hb-audit
#SBATCH --output=output/emi-hb-audit-%j.log
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?}"
export LANG=C LC_ALL=C PYTHONDONTWRITEBYTECODE=1
root="$PWD/output/emi-hb-audit-$SLURM_JOB_ID"
sif=/data/numerik/people/fchegini/containers/kaskade-robin-20261105/kaskade7-robin.sif
kaskade=/data/numerik/people/fchegini/kaskade7_new
python3 tests/select_emi_hb_inputs.py \
 --source "${HB_CAPTURE_INPUTS:-$PWD/output/emi-hb-audit-6143753/inputs}" \
 --output "$root" --calls "${HB_AUDIT_CALLS:-3}"
sha256sum tests/compare_emi_hb_frozen.cpp tests/summarize_emi_hb_audit.py \
 tests/select_emi_hb_inputs.py "$kaskade"/mg/bddc*.hh "$kaskade"/mg/bddc*.hpp > "$root/source-sha256.txt"
# Compile and replay inside the same Dune/MPI container as production.
srun --mpi=none --nodes=1 --ntasks=1 --ntasks-per-node=1 \
 apptainer exec --cleanenv --bind "$kaskade:/kaskade7" --bind "$PWD:$PWD" \
 --pwd "$PWD" "$sif" bash -c '
 set -euo pipefail
 export PATH=/usr/lib64/mpi/gcc/openmpi4/bin:$PATH
 export LD_LIBRARY_PATH=/usr/lib64/mpi/gcc/openmpi4/lib64:/usr/lib64/mpi/gcc/openmpi4/lib
 mpicxx -std=c++20 -O2 -I /kaskade7 tests/compare_emi_hb_frozen.cpp -o "$1/replay"
 for bits in 8 16; do
  for depth in 1 2 3 0; do
   echo "Starting bits=$bits depth=$depth"
   "$1/replay" "$1/inputs" "$1/replay-b$bits-d$depth" "$bits" "$depth"
  done
 done
 ' audit "$root"
sha256sum -c "$root/source-sha256.txt"
python3 tests/summarize_emi_hb_audit.py "$root"
echo "Audit complete: $root"
