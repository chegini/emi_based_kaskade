#!/usr/bin/env bash
set -euo pipefail
cd /BDDC/emiModelBddcSdc_new_scratch
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python3 tests/test_emi_transfer_roundtrip_audit.py
audit_build=$(mktemp -d /tmp/emi-transfer-audit-build.XXXXXX)
mpicxx -std=c++20 -O2 -I/kaskade7 tests/test_bddc_transfer_audit.cpp -o "$audit_build/test_capture"
"$audit_build/test_capture" "$audit_build/capture"
mpicxx -std=c++20 -O2 -Wall -Wextra -Werror -I/kaskade7 \
  tests/test_bddc_lifting_plan.cpp -o "$audit_build/test_lifting"
"$audit_build/test_lifting"
mpicxx -std=c++20 -O2 -I/kaskade7 tests/test_bddc_transfer_consistency.cpp \
  -o "$audit_build/test_transfer"
if ! "$audit_build/test_transfer" > "$audit_build/transfer.log" 2>&1; then
  tail -n 80 "$audit_build/transfer.log"
  exit 1
fi
tail -n 1 "$audit_build/transfer.log"
if [[ -f EMI_model.d ]]; then cp -- EMI_model.d "$audit_build/EMI_model.d"; fi
: > EMI_model.d
make -W EMI_model.cpp -o /kaskade7/libs/libkaskade.a \
  MPI=1 MPICXX=mpicxx \
  MPI_INCLUDE_FLAGS="-include /usr/lib64/mpi/gcc/openmpi4/include/mpi.h -DMUMPS_MPI_H"
python3 tests/run_emi_transfer_roundtrip_audit.py "$@"
