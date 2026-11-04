#!/usr/bin/env bash
set -euo pipefail
cd /BDDC/emiModelBddcSdc_new_scratch
test "$(git -C /kaskade7 branch --show-current)" = master
test "$(git branch --show-current)" = experiment/cg-fcg-constant-split
python3 tests/run_emi_cg_split_check.py "$@" --validate-only
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
check_dir=$(mktemp -d /tmp/emi-cg-split-build.XXXXXX)
mpicxx -std=c++20 -O2 tests/test_emi_cg_split_options.cpp -o "$check_dir/test_options"
"$check_dir/test_options"
python3 tests/test_emi_cg_split_commands.py
# Preserve and regenerate branch-dependent generated dependencies.
if [[ -f EMI_model.d ]]; then
  cp -- EMI_model.d "$check_dir/EMI_model.d"
  printf 'Previous dependency file saved in %s/EMI_model.d\n' "$check_dir"
fi
# gccmakedep requires an existing file; empty it only after preserving it.
: > EMI_model.d
make -W EMI_model.cpp -o /kaskade7/libs/libkaskade.a \
  MPI=1 MPICXX=mpicxx \
  MPI_INCLUDE_FLAGS="-include /usr/lib64/mpi/gcc/openmpi4/include/mpi.h -DMUMPS_MPI_H"
python3 tests/run_emi_cg_split_check.py "$@"
