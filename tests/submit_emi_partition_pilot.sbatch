#!/usr/bin/env bash
# Select preflight, smoke, smoke-raw, or pilot; this script never calls sbatch itself.
#SBATCH --partition=big
#SBATCH --constraint=Gold6342
#SBATCH --nodes=2
#SBATCH --ntasks=4
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --hint=nomultithread
#SBATCH --mem=64G
#SBATCH --time=01:00:00
#SBATCH --job-name=emi-partition
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?Submit from the EMI repository root}"
export LANG=C LC_ALL=C PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
emi_mode=${1:-preflight}
case "${EMI_MPI_BYTE_ENTROPY:-0}" in
  0) emi_entropy=() ;;
  1) emi_entropy=(--byte-entropy) ;;
  *) echo 'EMI_MPI_BYTE_ENTROPY must be 0 or 1' >&2; exit 2 ;;
esac
case "$emi_mode" in
  preflight) ;;
  smoke|smoke-raw)
    [[ $SLURM_JOB_NUM_NODES == 2 && $SLURM_NTASKS == 4 ]] || {
      echo 'Smoke check requires exactly 2 nodes and 4 ranks' >&2; exit 2;
    } ;;
  pilot)
    [[ $SLURM_JOB_NUM_NODES == 8 && $SLURM_NTASKS == 16 ]] || {
      echo 'Pilot requires exactly 8 nodes and 16 ranks' >&2; exit 2;
    } ;;
  *) echo 'Use preflight, smoke, smoke-raw, or pilot; no full sweep is implemented' >&2; exit 2 ;;
esac
emi_sif=${emi_sif:-/data/numerik/people/fchegini/containers/kaskade-robin-20261105/kaskade7-robin.sif}
emi_kaskade=${emi_kaskade:-/data/numerik/people/fchegini/kaskade7_new}
emi_output=${EMI_HB_OUTPUT_ROOT:-/scratch/htc/$USER/emi-hb}
mkdir -p "$emi_output"
emi_root="$emi_output/emi-partition-$emi_mode-$SLURM_JOB_ID"
mkdir "$emi_root"
scontrol show job "$SLURM_JOB_ID" > "$emi_root/allocation.txt"
sha256sum EMI_model.cpp ./*.hh Makefile tests/run_emi_partition_pilot.py tests/run_emi_partition_smoke.py \
  tests/test_emi_serialized_byte_entropy.py tests/test_bddc_serialized_byte_entropy.cpp \
  tests/test_emi_partition* tests/submit_emi_partition_pilot.sbatch \
  "$emi_kaskade"/mg/bddc*.hh "$emi_kaskade"/mg/bddc*.hpp > "$emi_root/source-sha256.txt"
# This is a fresh separate build; production executables and original results stay intact.
srun --mpi=none --nodes=1 --ntasks=1 --ntasks-per-node=1 \
  apptainer exec --cleanenv --bind "$emi_kaskade:/kaskade7" \
  --bind "$emi_output:$emi_output" --bind "$PWD:$PWD" --pwd "$PWD" "$emi_sif" \
  /bin/bash tests/build_emi_coarse_batch_cluster.sh "$emi_root/build"
srun --mpi=none --nodes=1 --ntasks=1 --ntasks-per-node=1 \
  apptainer exec --cleanenv --bind "$emi_kaskade:/kaskade7" \
  --bind "$emi_output:$emi_output" --bind "$PWD:$PWD" --pwd "$PWD" "$emi_sif" /bin/bash -c '
    set -euo pipefail
    export PYTHONDONTWRITEBYTECODE=1 PATH=/usr/lib64/mpi/gcc/openmpi4/bin:$PATH
    export LD_LIBRARY_PATH=/usr/lib64/mpi/gcc/openmpi4/lib64:/usr/lib64/mpi/gcc/openmpi4/lib
    mpicxx -std=c++20 -O2 -Wall -Wextra -Werror -I . -I /kaskade7 tests/test_emi_partitioning.cpp -o "$1/partition-policy"
    "$1/partition-policy"
    mpicxx -std=c++20 -O2 -Wall -Wextra -Werror -DKASKADE_HAVE_MPI -I . -I /kaskade7 tests/test_emi_partition_runtime.cpp -o "$1/partition-runtime"
    mpicxx -std=c++20 -O2 -Wall -Wextra -Werror -DKASKADE_HAVE_MPI -I . -I /kaskade7 tests/test_emi_interface_exchange.cpp -o "$1/interface-owner-routing"
    CXX=mpicxx python3 tests/test_emi_partition_mpi.py --kaskade /kaskade7 --output "$1/owner-routing" --build-only
    CXX=mpicxx KASKADE7=/kaskade7 python3 tests/test_emi_partition_report.py
    python3 tests/test_emi_partition_pilot.py
    python3 tests/test_emi_partition_smoke.py
    CXX=mpicxx KASKADE7=/kaskade7 python3 tests/test_emi_partition_options.py
    CXX=mpicxx KASKADE7=/kaskade7 python3 tests/test_emi_serialized_byte_entropy.py
    CXX=mpicxx KASKADE7=/kaskade7 python3 -m unittest discover -s tests -p test_emi_mpi_routing_check.py
    CXX=mpicxx KASKADE7=/kaskade7 python3 -m unittest discover -s tests -p test_emi_histogram_start_check.py
    mpicxx -std=c++20 -O2 -I /kaskade7 tests/test_bddc_split_recenter.cpp -o "$1/partition-codec"
    KASKADE_LIFTING_OPTIMIZATION=4 python3 tests/check_emi_split_recenter.py --executable "$1/partition-codec" --output "$1/partition-codec-check" --bits 16
  ' partition-tests "$emi_root/build"
srun --mpi=pmix_v5 --nodes="$SLURM_JOB_NUM_NODES" --ntasks="$SLURM_NTASKS" --ntasks-per-node=2 \
  --distribution=block:block --cpu-bind=cores --kill-on-bad-exit=1 \
  apptainer exec --bind "$emi_kaskade:/kaskade7" --bind "$emi_output:$emi_output" \
  --bind "$PWD:$PWD" --pwd "$PWD" "$emi_sif" /bin/bash -c '
    set -euo pipefail
    export PATH=/usr/lib64/mpi/gcc/openmpi4/bin:$PATH
    export LD_LIBRARY_PATH=/usr/lib64/mpi/gcc/openmpi4/lib64:/usr/lib64/mpi/gcc/openmpi4/lib
    export OMPI_MCA_pml=ob1 OMPI_MCA_btl=self,vader,tcp OMPI_MCA_btl_vader_single_copy_mechanism=none
    "$1/partition-runtime" "$2/ownership-preflight" "$3"
    "$1/interface-owner-routing"
    "$1/owner-routing/ownership-mpi"
    exec "$1/owner-routing/ownership-mpi" "$2/serialized-byte-entropy-preflight"
  ' partition-mpi "$emi_root/build" "$emi_root" "$SLURM_JOB_NUM_NODES" > "$emi_root/mpi-preflight.log" 2>&1
sha256sum --quiet -c "$emi_root/source-sha256.txt"
if [[ $emi_mode == preflight ]]; then
  echo "PREFLIGHT complete: $emi_root; no EMI trajectories executed"
  exit 0
fi
srun --mpi=none --nodes="$SLURM_JOB_NUM_NODES" --ntasks="$SLURM_JOB_NUM_NODES" --ntasks-per-node=1 \
  /bin/bash -c 'hostname; lscpu' > "$emi_root/cpu-inventory.log"
emi_pools=()
[[ -z ${EMI_PARTITION_INTRA_NODES:-} ]] || emi_pools+=(--intra-nodes "$EMI_PARTITION_INTRA_NODES")
[[ -z ${EMI_PARTITION_EXTRA_NODES:-} ]] || emi_pools+=(--extra-nodes "$EMI_PARTITION_EXTRA_NODES")
set +e
if [[ $emi_mode == smoke || $emi_mode == smoke-raw ]]; then
emi_raw=()
[[ $emi_mode != smoke-raw ]] || emi_raw+=(--uncompressed)
python3 -u tests/run_emi_partition_smoke.py \
  --executable "$emi_root/build/emiModel" --kaskade "$emi_kaskade" \
  --slurm-container "$emi_sif" --output "$emi_root/results" "${emi_raw[@]}" "${emi_entropy[@]}"
else
python3 -u tests/run_emi_partition_pilot.py \
  --executable "$emi_root/build/emiModel" --kaskade "$emi_kaskade" \
  --slurm-container "$emi_sif" --output "$emi_root/results" \
  --nodes 8 --ranks 16 --steps 3 --repeats 3 --detail-repeats 3 --warmups 1 \
  "${emi_pools[@]}" "${emi_entropy[@]}"
fi
emi_status=$?
set -e
sha256sum --quiet -c "$emi_root/source-sha256.txt"
echo "Controlled EMI partition $emi_mode complete: $emi_root/results; exit=$emi_status"
exit "$emi_status"
