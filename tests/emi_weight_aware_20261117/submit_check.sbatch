#!/usr/bin/env bash
# Fresh-build validation only. This file never submits another SLURM job.
#SBATCH --partition=big
#SBATCH --constraint=Gold6342
#SBATCH --nodes=2
#SBATCH --ntasks=4
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --hint=nomultithread
#SBATCH --mem=64G
#SBATCH --time=00:30:00
#SBATCH --job-name=pepe-split
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?Submit from emiModelBddcSdc_new_scratch}"
[[ ${SLURM_JOB_NUM_NODES:?} == 2 && ${SLURM_NTASKS:?} == 4 ]] || exit 2
bits=${1:-16}
case "$bits" in 12|16) ;; *) echo 'Use precision 12 or 16' >&2; exit 2 ;; esac
export LANG=C LC_ALL=C PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
kaskade=${emi_kaskade:-/data/numerik/people/fchegini/kaskade7_new}
sif=${emi_sif:-/data/numerik/people/fchegini/containers/kaskade-robin-20261105/kaskade7-robin.sif}
out=${EMI_WEIGHT_AWARE_OUTPUT_ROOT:-/scratch/htc/$USER/emi-weight-aware}
root="$out/check-$SLURM_JOB_ID"
package="$PWD/tests/emi_weight_aware_20261117"
mkdir -p "$out"
mkdir "$root"
[[ -f "$kaskade/mg/bddcConstantSplit.hh" ]] || {
  echo 'The corrected weight-aware splitting header is missing' >&2; exit 2;
}
scontrol show job "$SLURM_JOB_ID" > "$root/allocation.txt"
git rev-parse HEAD > "$root/emi-git-head.txt"
git -C "$kaskade" rev-parse HEAD > "$root/kaskade-git-head.txt"
# Include the new policy header, all BDDC headers, sources, scripts and linked libraries.
sha256sum EMI_model.cpp ./*.hh Makefile tests/*.py tests/*.cpp \
  tests/build_emi_coarse_batch_cluster.sh "$package"/*.py "$package"/*.sbatch \
  "$kaskade"/mg/bddc*.hh "$kaskade"/mg/bddc*.hpp \
  "$kaskade"/libs/*.o "$kaskade"/libs/libkaskade.a > "$root/source-sha256.txt"
trap 'sha256sum --quiet -c "$root/source-sha256.txt" || exit 1' EXIT
srun --mpi=none --nodes=1 --ntasks=1 --ntasks-per-node=1 --cpu-bind=cores \
  apptainer exec --cleanenv --bind "$kaskade:/kaskade7" --bind "$PWD:$PWD" \
  --bind "$out:$out" --pwd "$PWD" "$sif" \
  /bin/bash tests/build_emi_coarse_batch_cluster.sh "$root/build"
# Run native checks with real Dune/MPI headers in the same production container.
srun --mpi=none --nodes=1 --ntasks=1 --ntasks-per-node=1 --cpu-bind=cores \
  apptainer exec --cleanenv --bind "$kaskade:/kaskade7" --bind "$PWD:$PWD" \
  --bind "$out:$out" --pwd "$PWD" "$sif" /bin/bash -c '
    set -euo pipefail
    export PATH=/usr/lib64/mpi/gcc/openmpi4/bin:$PATH
    export LD_LIBRARY_PATH=/usr/lib64/mpi/gcc/openmpi4/lib64:/usr/lib64/mpi/gcc/openmpi4/lib
    CXX=mpicxx KASKADE7=/kaskade7 python3 -m unittest discover \
      -s tests/emi_weight_aware_20261117 -p test_check.py
    mpicxx -std=c++20 -O2 -I /kaskade7 tests/test_bddc_constant_split.cpp -o "$1/constant-split"
    "$1/constant-split"
    mpicxx -std=c++20 -O2 -I /kaskade7 tests/test_bddc_split_recenter.cpp -o "$1/split-recenter"
    KASKADE_LIFTING_OPTIMIZATION=4 python3 tests/check_emi_split_recenter.py \
      --executable "$1/split-recenter" --output "$1/recenter-parity" --bits "$2"
  ' native-checks "$root/build" "$bits"
sha256sum --quiet -c "$root/source-sha256.txt"
set +e
python3 -u "$package/run_check.py" --executable "$root/build/emiModel" \
  --kaskade "$kaskade" --slurm-container "$sif" --output "$root/results" --bits "$bits" --execute
status=$?
set -e
echo "Results: $root/results; Pepe weight-aware validation exit=$status"
exit "$status"
