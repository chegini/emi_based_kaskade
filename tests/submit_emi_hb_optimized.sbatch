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
#SBATCH --time=24:00:00
#SBATCH --job-name=emi-hb-opt
#SBATCH --output=output/emi-hb-opt-%j.log
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?Submit from the EMI repository root}"
export LANG=C LC_ALL=C PYTHONDONTWRITEBYTECODE=1
emi_sif=${emi_sif:-/data/numerik/people/fchegini/containers/kaskade-robin-20261105/kaskade7-robin.sif}
emi_kaskade=${emi_kaskade:-/data/numerik/people/fchegini/kaskade7_new}
emi_mode=${1:-pilot}
case "$emi_mode" in
  pilot) emi_case=${2:-pepe} ;;
  strong) emi_case=${2:-weak_32_h16} ;;
  weak) emi_case="weak_$((2*SLURM_NTASKS))_h16" ;;
  *) echo 'Use pilot [mesh], strong [fixed-mesh], or weak' >&2; exit 2 ;;
esac
emi_output=${EMI_HB_OUTPUT_ROOT:-/scratch/htc/$USER/emi-hb}
mkdir -p "$emi_output"
emi_root="$emi_output/emi-hb-opt-$SLURM_JOB_ID"
mkdir "$emi_root"
scontrol show job "$SLURM_JOB_ID" > "$emi_root/allocation.txt"
# Retain exact lifting sources alongside the executable for reproducibility.
tar -cf "$emi_root/lifting-sources.tar" -C "$emi_kaskade" \
  mg/bddcLiftingPlan.hh mg/bddcLiftingPerformance.hh \
  mg/bddcSpaceTransfer.hh mg/bddcSpaceTransfer.hpp \
  mg/bddcTransferTiming.hh mg/bddcRecenter.hh
sha256sum EMI_model.cpp ./*.hh Makefile tests/*.py tests/*.cpp tests/*.sh \
  tests/submit_emi_hb_optimized.sbatch \
  "$emi_kaskade"/mg/bddc*.hh "$emi_kaskade"/mg/bddc*.hpp \
  > "$emi_root/source-sha256.txt"
srun --mpi=none --nodes=1 --ntasks=1 --ntasks-per-node=1 \
  apptainer exec --cleanenv --bind "$emi_kaskade:/kaskade7" \
  --bind "$emi_output:$emi_output" --bind "$PWD:$PWD" --pwd "$PWD" "$emi_sif" \
  /bin/bash tests/build_emi_coarse_batch_cluster.sh "$emi_root/build"
# Native recorder and production-plan tests are independent of the PDE solve.
srun --mpi=none --nodes=1 --ntasks=1 --ntasks-per-node=1 \
  apptainer exec --cleanenv --bind "$emi_kaskade:/kaskade7" \
  --bind "$emi_output:$emi_output" --bind "$PWD:$PWD" --pwd "$PWD" "$emi_sif" /bin/bash -c '
    set -euo pipefail
    export PYTHONDONTWRITEBYTECODE=1
    export PATH=/usr/lib64/mpi/gcc/openmpi4/bin:$PATH
    export LD_LIBRARY_PATH=/usr/lib64/mpi/gcc/openmpi4/lib64:/usr/lib64/mpi/gcc/openmpi4/lib
    CXX=mpicxx KASKADE7=/kaskade7 python3 tests/test_emi_hb_time_study.py
    CXX=mpicxx KASKADE_TRANSFER_TIMING=1 python3 tests/test_multilevel_hb_transfer.py --kaskade-root /kaskade7
    python3 tests/test_emi_hb_time_analysis.py
    python3 tests/test_lifting_performance.py --kaskade /kaskade7 --output "$1/plan-regression"
    python3 tests/test_emi_hb_time_mpi.py --kaskade /kaskade7 --executable "$1/hb-time-mpi-test" --build-only
  ' hb-test "$emi_root/build"
srun --mpi=pmix_v5 --nodes="$SLURM_JOB_NUM_NODES" --ntasks="$SLURM_NTASKS" --ntasks-per-node=2 \
  --cpu-bind=cores --kill-on-bad-exit=1 apptainer exec \
  --bind "$emi_kaskade:/kaskade7" --bind "$emi_output:$emi_output" --bind "$PWD:$PWD" --pwd "$PWD" "$emi_sif" /bin/bash -c '
    set -euo pipefail
    export PYTHONDONTWRITEBYTECODE=1
    export PATH=/usr/lib64/mpi/gcc/openmpi4/bin:$PATH
    export LD_LIBRARY_PATH=/usr/lib64/mpi/gcc/openmpi4/lib64:/usr/lib64/mpi/gcc/openmpi4/lib
    export OMPI_MCA_pml=ob1 OMPI_MCA_btl=self,vader,tcp
    export OMPI_MCA_btl_vader_single_copy_mechanism=none
    export KASKADE_TRANSFER_TIMING=2 KASKADE_TIMESTEP_BYTES=1
    exec "$1/hb-time-mpi-test"
  ' hb-mpi "$emi_root/build" > "$emi_root/mpi-preflight.log" 2>&1
srun --mpi=none --nodes="$SLURM_JOB_NUM_NODES" --ntasks="$SLURM_JOB_NUM_NODES" --ntasks-per-node=1 \
  /bin/bash -c 'hostname; lscpu' > "$emi_root/cpu-inventory.log"
sha256sum --quiet -c "$emi_root/source-sha256.txt"
read -r -a emi_bits <<< "${EMI_HB_BITS:-16}"
read -r -a emi_tiers <<< "${EMI_HB_TIERS:-0 1}"
read -r -a emi_optimizations <<< "${EMI_HB_OPTIMIZATIONS:-0 1 4 5}"
emi_diagnostics=()
if [[ ${EMI_HB_PROFILE_SETUP:-0} == 1 ]]; then emi_diagnostics+=(--hierarchy-profile); fi
set +e
python3 -u tests/run_emi_hb_time_study.py \
  --executable "$emi_root/build/emiModel" --kaskade "$emi_kaskade" \
  --slurm-container "$emi_sif" --output "$emi_root/results" \
  --case "$emi_case" --scaling "$emi_mode" \
  --nodes "$SLURM_JOB_NUM_NODES" --ranks "$SLURM_NTASKS" \
  --steps "${EMI_HB_STEPS:-3}" --bits "${emi_bits[@]}" --depths 1 4 \
  --lifting-optimizations "${emi_optimizations[@]}" "${emi_diagnostics[@]}" \
  --repeats "${EMI_HB_REPEATS:-3}" --detail-repeats "${EMI_HB_DETAIL_REPEATS:-3}" \
  --tiers "${emi_tiers[@]}" --warmups 1 --case-minutes "${EMI_HB_CASE_MINUTES:-60}"
emi_status=$?
set -e
sha256sum --quiet -c "$emi_root/source-sha256.txt"
echo "Optimized fixed-topology EMI study complete: $emi_root/results; exit=$emi_status"
exit "$emi_status"
