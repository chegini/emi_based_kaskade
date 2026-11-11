#!/usr/bin/env bash
# Build only, inside one compute-node container. No nested MPI launcher.
set -euo pipefail
cd "$(dirname "$0")/.."
export PATH=/usr/lib64/mpi/gcc/openmpi4/bin:/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin
export LD_LIBRARY_PATH=/usr/lib64/mpi/gcc/openmpi4/lib64:/usr/lib64/mpi/gcc/openmpi4/lib
export LANG=C LC_ALL=C OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
python3 tests/test_emi_mpi_routing_check.py
check_dir=$(mktemp -d /tmp/emi-routing-cluster-build.XXXXXX)
mpicxx -std=c++20 -O2 -Wall -Wextra -Werror \
  tests/test_emi_local_support_diagnostic.cpp -o "$check_dir/local-support"
"$check_dir/local-support"
CXX=mpicxx python3 tests/test_emi_local_support_extraction.py
mpicxx -std=c++20 -O2 -I /kaskade7 \
  /kaskade7/tests/bddcTransferOrder/bddcTransferOrder.cpp -o "$check_dir/transfer-order"
(cd "$check_dir" && ./transfer-order)
if [[ -f EMI_model.d ]]; then
  cp -- EMI_model.d "$check_dir/EMI_model.d"
else
  touch EMI_model.d
fi
if ! make -W EMI_model.cpp -o /kaskade7/libs/libkaskade.a \
    KASKADE7=/kaskade7 MPI=1 MPICXX=mpicxx \
    MPI_INCLUDE_FLAGS="-include /usr/lib64/mpi/gcc/openmpi4/include/mpi.h -DMUMPS_MPI_H"; then
  if [[ -f gccerr.txt ]]; then tail -80 gccerr.txt; fi
  exit 1
fi
echo 'Cluster EMI build passed'
