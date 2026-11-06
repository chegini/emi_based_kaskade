#!/usr/bin/env bash
# Single-node validation inside the Open MPI Kaskade container.
set -euo pipefail
cd "$(dirname "$0")/.."
python3 tests/run_emi_histogram_start_check.py "$@" --plan-only
kaskade_dir=${KASKADE7:-/kaskade7}
mpi_dir=/usr/lib64/mpi/gcc/openmpi4
if [[ ! -x "$mpi_dir/bin/mpicxx" ]]; then
  echo 'Run inside the Kaskade Open MPI container.' >&2
  exit 2
fi
export PATH="$mpi_dir/bin:/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin"
export LD_LIBRARY_PATH="$mpi_dir/lib64:$mpi_dir/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 LANG=C LC_ALL=C
export OMPI_MCA_pml=ob1 OMPI_MCA_btl=self,vader,tcp
export PYTHONDONTWRITEBYTECODE=1
unset KASKADE_FROZEN_TRANSFER_DIR KASKADE_RESTRICTION_CAPTURE KASKADE_TRANSFER_AUDIT_DIR
check_dir=$(mktemp -d /tmp/emi-histogram-start.XXXXXX)
mpicxx -std=c++20 -O2 tests/test_emi_bddc_step_schedule.cpp -o "$check_dir/step-schedule"
"$check_dir/step-schedule"
mpicxx -std=c++20 -O2 -I "$kaskade_dir" \
  "$kaskade_dir/tests/bddcTransferOrder/bddcTransferOrder.cpp" -o "$check_dir/transfer-order"
"$check_dir/transfer-order"
python3 tests/test_emi_histogram_start_check.py
python3 tests/test_emi_compression_report.py
if [[ -f EMI_model.d ]]; then
  cp -- EMI_model.d "$check_dir/EMI_model.d"
else
  touch EMI_model.d
fi
if ! make -W EMI_model.cpp -o "$kaskade_dir/libs/libkaskade.a" \
    KASKADE7="$kaskade_dir" MPI=1 MPICXX=mpicxx \
    MPI_INCLUDE_FLAGS="-include $mpi_dir/include/mpi.h -DMUMPS_MPI_H"; then
  if [[ -f gccerr.txt ]]; then tail -80 gccerr.txt; fi
  exit 1
fi
python3 tests/run_emi_histogram_start_check.py "$@"
