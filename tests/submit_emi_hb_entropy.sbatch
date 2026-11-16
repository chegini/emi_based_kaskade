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
#SBATCH --time=24:00:00
#SBATCH --job-name=emi-hb-entropy
#SBATCH --output=/scratch/htc/%u/emi-hb/emi-hb-entropy-%j.log
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?Submit from EMI repository root}"
export LANG=C LC_ALL=C PYTHONDONTWRITEBYTECODE=1
emi_sif=${emi_sif:-/data/numerik/people/fchegini/containers/kaskade-robin-20261105/kaskade7-robin.sif}
emi_kaskade=${emi_kaskade:-/data/numerik/people/fchegini/kaskade7_new}
export emi_sif emi_kaskade
emi_output=${EMI_HB_OUTPUT_ROOT:-/scratch/htc/$USER/emi-hb}
export EMI_HB_OUTPUT_ROOT="$emi_output"
mkdir -p "$emi_output"
case "${EMI_HB_BITS:-16}" in 8|12|16|24|32) ;; *) echo 'Set EMI_HB_BITS to ONE of 8,12,16,24,32' >&2;exit 2 ;; esac
case "${EMI_HB_PREFLIGHT_ONLY:-0}" in 0|1) ;; *) echo 'EMI_HB_PREFLIGHT_ONLY must be 0 or 1' >&2;exit 2 ;; esac
if [[ "${EMI_HB_PREFLIGHT_ONLY:-0}" == 1 ]];then
  if [[ "$SLURM_JOB_NUM_NODES" != 1 || "$SLURM_NTASKS" != 1 ]];then
    echo 'Preflight-only mode requires 1 node and 1 task.' >&2;exit 2
  fi
elif [[ "$SLURM_JOB_NUM_NODES" != 8 || "$SLURM_NTASKS" != 16 ]];then
  echo 'This entropy campaign requires 8 nodes, 16 ranks (2/node).' >&2;exit 2
fi
emi_root="$emi_output/emi-hb-entropy-$SLURM_JOB_ID"
mkdir "$emi_root"
# Never inherit entropy instrumentation into the ordinary benchmark executable.
unset EMI_MESSAGE_ENTROPY_PREFIX APPTAINERENV_EMI_MESSAGE_ENTROPY_PREFIX
sha256sum tests/*emi*entropy*.py tests/emi_message_entropy.hh tests/submit_emi_hb_entropy.sbatch \
  "$emi_kaskade"/mg/bddc*.hh "$emi_kaskade"/mg/bddc*.hpp > "$emi_root/source-sha256.txt"
python3 tests/test_emi_entropy_workflow.py
python3 tests/prepare_emi_entropy_overlay.py --source "$emi_kaskade" --output "$emi_root/kaskade-observer"
# Test the observer before spending hours on the repeated PDE study.
srun --mpi=none --nodes=1 --ntasks=1 --ntasks-per-node=1 \
  apptainer exec --cleanenv --bind "$emi_kaskade:$emi_kaskade" \
  --bind "$emi_root/kaskade-observer:/kaskade7" \
  --bind "$emi_output:$emi_output" --bind "$PWD:$PWD" --pwd "$PWD" "$emi_sif" /bin/bash -c '
    set -euo pipefail
    export PYTHONDONTWRITEBYTECODE=1
    export PATH=/usr/lib64/mpi/gcc/openmpi4/bin:$PATH
    export LD_LIBRARY_PATH=/usr/lib64/mpi/gcc/openmpi4/lib64:/usr/lib64/mpi/gcc/openmpi4/lib
    CXX=mpicxx python3 tests/test_emi_message_entropy.py
    CXX=mpicxx python3 tests/check_emi_entropy_observer.py --original "$2" \
      --overlay /kaskade7 --output "$1/observer-parity"
    /bin/bash tests/build_emi_coarse_batch_cluster.sh "$1/build"
  ' entropy-preflight "$emi_root" "$emi_kaskade"
if [[ "${EMI_HB_PREFLIGHT_ONLY:-0}" == 1 ]];then
  sha256sum --quiet -c "$emi_root/source-sha256.txt"
  echo "PREFLIGHT PASS: native entropy tests, original/observer packet parity, and full EMI build; $emi_root"
  echo 'No performance or entropy simulation trajectories launched.'
  exit 0
fi
# Original driver, original sources, original executable: entropy disabled for all repeats.
set +e
/bin/bash tests/submit_emi_hb_depth_split.sbatch "${1:-pilot}" "${2:-pepe}"
emi_performance_status=$?
set -e
emi_performance="$emi_output/emi-hb-depth-split-$SLURM_JOB_ID/results"
if [[ ! -f "$emi_performance/combined-comparison.csv" ]];then
  echo 'Performance study did not complete; do not diagnose an incomplete baseline.' >&2
  exit 1
fi
# Low-bit configurations can fail convergence; retain those failures and only
# diagnose configurations with complete, exact-validated trajectories.
printf '%s\n' "$emi_performance" > "$emi_root/performance-results.txt"
set +e
python3 -u tests/run_emi_message_entropy.py --performance-results "$emi_performance" \
  --executable "$emi_root/build/emiModel" --overlay "$emi_root/kaskade-observer" \
  --output "$emi_root/results"
emi_entropy_status=$?
set -e
sha256sum --quiet -c "$emi_root/source-sha256.txt"
echo "Runtime results (entropy OFF): $emi_performance"
echo "Message entropy results: $emi_root/results"
echo "Performance exit=$emi_performance_status; entropy exit=$emi_entropy_status"
if (( emi_performance_status != 0 || emi_entropy_status != 0 ));then exit 1;fi
