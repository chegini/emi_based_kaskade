#!/usr/bin/env bash
# Run inside one compute-node container. Build into a fresh separate directory.
set -euo pipefail
emi_source=$PWD
emi_build=${1:?Supply a fresh absolute build directory}
export PATH=/usr/lib64/mpi/gcc/openmpi4/bin:/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin
export LD_LIBRARY_PATH=/usr/lib64/mpi/gcc/openmpi4/lib64:/usr/lib64/mpi/gcc/openmpi4/lib
export LANG=C LC_ALL=C OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
mkdir "$emi_build"
cp "$emi_source"/EMI_model.cpp "$emi_source"/*.hh "$emi_source"/Makefile "$emi_build/"
mkdir "$emi_build/source-mg"
cp /kaskade7/mg/bddc.hh /kaskade7/mg/bddc.hpp \
  /kaskade7/mg/bddcCoarseReduction.hh "$emi_build/source-mg/"
mpicxx -std=c++20 -O2 -Wall -Wextra -Werror -I /kaskade7 \
  tests/test_bddc_coarse_reduction.cpp -o "$emi_build/coarse-reduction-test"
CXX=mpicxx python3 tests/test_bddc_flexible_cg.py --kaskade-root /kaskade7
CXX=mpicxx python3 tests/test_emi_local_support_extraction.py
cd "$emi_build"
touch EMI_model.d
if ! make -W EMI_model.cpp -o /kaskade7/libs/libkaskade.a \
    KASKADE7=/kaskade7 MPI=1 MPICXX=mpicxx \
    MPI_INCLUDE_FLAGS="-include /usr/lib64/mpi/gcc/openmpi4/include/mpi.h -DMUMPS_MPI_H"; then
  if [[ -f gccerr.txt ]]; then tail -80 gccerr.txt; fi
  exit 1
fi
python3 "$emi_source/tests/verify_emi_coarse_build.py" \
  --build "$emi_build" --repo "$emi_source" --kaskade /kaskade7 --write-manifest
echo "Separate coarse-batch executable built: $emi_build/emiModel"
