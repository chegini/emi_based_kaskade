#!/usr/bin/env bash
# First-stage, single-node validation inside the Open MPI Kaskade container.
set -euo pipefail
cd "$(dirname "$0")/.."
python3 tests/run_emi_mpi_routing_check.py "$@" --plan-only
kaskade_dir=${KASKADE7:-/kaskade7}
mpi_dir=/usr/lib64/mpi/gcc/openmpi4
test -x "$mpi_dir/bin/mpicxx"
export PATH="$mpi_dir/bin:/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin"
export LD_LIBRARY_PATH="$mpi_dir/lib64:$mpi_dir/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 LANG=C LC_ALL=C PYTHONDONTWRITEBYTECODE=1
export OMPI_MCA_pml=ob1 OMPI_MCA_btl=self,vader,tcp
export OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1
unset KASKADE_FROZEN_TRANSFER_DIR KASKADE_RESTRICTION_CAPTURE KASKADE_TRANSFER_AUDIT_DIR
check_dir=$(mktemp -d /tmp/emi-mpi-routing.XXXXXX)
mpicxx -std=c++20 -O2 -Wall -Wextra -Werror -DKASKADE_HAVE_MPI -I "$kaskade_dir" \
  "$kaskade_dir/tests/bddcMpiRouting/bddcMpiRouting.cpp" -o "$check_dir/routing"
mpirun -np 2 "$check_dir/routing" 1
# Stay within the requested two CPUs; larger ranks are tested separately on the cluster.
MPIEXEC=mpirun MPICXX=mpicxx python3 tests/test_bddc_mpi_exchange.py \
  --kaskade "$kaskade_dir" --ranks 1 2
python3 tests/test_emi_mpi_routing_check.py
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
python3 tests/run_emi_mpi_routing_check.py "$@"
