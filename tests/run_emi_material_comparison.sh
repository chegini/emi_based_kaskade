#!/usr/bin/env bash
# Run manually inside the MPI-enabled Kaskade container.
set -euo pipefail
repo_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
cd "$repo_dir"

mode=${1:-audit}
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
emi_mpiexec=${MPIEXEC:-mpirun}
emi_mpicxx=${MPICXX:-mpicxx}
emi_kaskade=${KASKADE_ROOT:-/kaskade7}

physics=(
  --mpi 1 --sdc 0
  --input ../input_emi_mesh/pepe_combi_domi.vtu
  --extra_set ../input_emi_mesh/pepe_combi_domi_extracellular.txt
  --intra_set ../input_emi_mesh/pepe_combi_domi_intracellular.txt
  --excited ../input_emi_mesh/pepe_combi_domi_excited.txt
  --C_m 1 --R 0.1 --dt 0.01
  --nThreads 1 --profile 1 --diagnostics 0
)
bddc=(
  --bddc 1 --direct 0 --bddcCompression 0 --masterWorker 0
  --bddcSplitInterfaces 0 --bddcInterfaceTypes 7
  --bddcUseCg 1 --bddcRelativeTolerance 1e-7 --bddcVerbose 0
  --bddcInterfaceAssembly 2 --bddcInterfaceAssemblyCheck 0
)

run_case()
{
  local name=$1 ranks=$2 vtk=$3 steps=$4
  shift 4
  local final_time=0.1
  if [[ "$steps" == 1 ]]; then final_time=0.01; fi
  local out="output/$name"
  mkdir -p "$out"
  "$emi_mpiexec" --allow-run-as-root -np "$ranks" \
    -x OPENBLAS_NUM_THREADS -x OMP_NUM_THREADS \
    ./emiModel "${physics[@]}" --vtk "$vtk" \
    --finalTime "$final_time" --maximumNumberOfTimeSteps "$steps" \
    "$@" --dir "$out" 2>&1 | tee "$out/run.log"
}

unit_tests()
{
  local unit_dir
  unit_dir=$(mktemp -d /tmp/emi-material-tests.XXXXXX)
  # Keep executables available for inspection/re-running; no repo artifacts.
  c++ -std=c++20 -O2 -Wall -Wextra -pedantic -I"$emi_kaskade" \
    tests/test_emi_material_assembly.cpp -o "$unit_dir/material_serial"
  "$unit_dir/material_serial"
  "$emi_mpicxx" -std=c++20 -O2 -Wall -Wextra -pedantic \
    -DKASKADE_HAVE_MPI -I"$emi_kaskade" \
    tests/test_emi_material_assembly.cpp -o "$unit_dir/material_mpi"
  local ranks
  for ranks in 1 2 8; do
    "$emi_mpiexec" --allow-run-as-root -np "$ranks" "$unit_dir/material_mpi"
  done
  printf 'Standalone test executables: %s\n' "$unit_dir"
}

case "$mode" in
  unit)
    unit_tests
    ;;
  audit)
    unit_tests
    for ranks in 1 8; do
      run_case "pepe-material-owner-audit-np${ranks}-one-step" "$ranks" 0 1 \
        "${bddc[@]}" --bddcMaterialAssembly 1 --bddcMaterialAssemblyCheck 1
    done
    ;;
  accuracy)
    direct_reference=output/pepe-direct-interface-check-10steps/emiLast.vtu
    for reference in "$direct_reference" \
      output/pepe-bddc-interface-check-np1-10steps/emiBDDCLast.vtu \
      output/pepe-bddc-interface-check-np8-10steps/emiBDDCLast.vtu; do
      if [[ ! -f "$reference" ]]; then
        printf 'Missing saved ten-step reference: %s\n' "$reference" >&2
        exit 1
      fi
    done
    for ranks in 1 8; do
      name="pepe-material-owner-accuracy-np${ranks}-10steps"
      run_case "$name" "$ranks" 1 10 \
        "${bddc[@]}" --bddcMaterialAssembly 1 --bddcMaterialAssemblyCheck 0
      python3 compare_emi_vtu.py "output/$name/emiBDDCLast.vtu" "$direct_reference" \
        | tee "output/$name/compare-direct.log"
      python3 compare_emi_vtu.py "output/$name/emiBDDCLast.vtu" \
        "output/pepe-bddc-interface-check-np${ranks}-10steps/emiBDDCLast.vtu" \
        | tee "output/$name/compare-replicated.log"
    done
    ;;
  benchmark)
    # Identical ten-step physics/output settings. No VTU comparison in this
    # mode: accuracy is checked separately against the saved VTK references.
    run_case pepe-material-benchmark-direct-np1-10steps 1 0 10 --bddc 0 --direct 1
    for material_mode in 0 1; do
      for ranks in 1 8; do
        run_case "pepe-material-benchmark-mode${material_mode}-np${ranks}-10steps" "$ranks" 0 10 \
          "${bddc[@]}" --bddcMaterialAssembly "$material_mode" --bddcMaterialAssemblyCheck 0
      done
    done
    ;;
  *)
    printf 'Usage: bash tests/run_emi_material_comparison.sh {unit|audit|accuracy|benchmark}\n' >&2
    exit 2
    ;;
esac
