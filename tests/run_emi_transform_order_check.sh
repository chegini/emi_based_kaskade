#!/usr/bin/env bash
# Run inside the Open MPI container with /kaskade7 and /BDDC mounted.
set -euo pipefail
cd "$(dirname "$0")/.."
kaskade_dir=${KASKADE7:-/kaskade7}
mpi_dir=/usr/lib64/mpi/gcc/openmpi4
if [[ ! -x "$mpi_dir/bin/mpicxx" ]]; then
  echo 'Run inside the Kaskade Open MPI container.' >&2
  exit 2
fi
export PATH="$mpi_dir/bin:/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin"
export LD_LIBRARY_PATH="$mpi_dir/lib64:$mpi_dir/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 LANG=C LC_ALL=C
# Single-node validation only.
export OMPI_MCA_pml=ob1 OMPI_MCA_btl=self,vader,tcp
unset KASKADE_FROZEN_TRANSFER_DIR KASKADE_RESTRICTION_CAPTURE KASKADE_TRANSFER_AUDIT_DIR
check_dir=$(mktemp -d /tmp/emi-transform-order.XXXXXX)
mpicxx -std=c++20 -O2 -I "$kaskade_dir" \
  "$kaskade_dir/tests/bddcTransferOrder/bddcTransferOrder.cpp" -o "$check_dir/transfer-order"
mpicxx -std=c++20 -O2 -I "$kaskade_dir" tests/test_bddc_constant_split.cpp \
  -o "$check_dir/constant-split"
"$check_dir/constant-split"
mpicxx -std=c++20 -O2 tests/test_emi_cg_split_options.cpp -o "$check_dir/split-options"
"$check_dir/split-options"
"$check_dir/transfer-order"
PYTHONDONTWRITEBYTECODE=1 python3 tests/test_emi_quantization_diagnostic.py
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
# Same original-residual CG scalar convention as the previous Pepe campaigns.
# Defaults retain CG/off; --solvers and --split select the full experiment matrix.
python3 tests/run_emi_quantization_diagnostic.py \
  --case pepe --ranks 2 --refine 0 --bits 12 16 24 \
  --iterations 1000 --cg-original-residual 1 --transform-orders before after "$@"
