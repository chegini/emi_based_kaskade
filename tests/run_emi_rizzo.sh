#!/usr/bin/env bash
# Run manually inside the MPI-enabled Kaskade container after rebuilding emiModel.
set -euo pipefail
repo_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
cd "$repo_dir"
quiet=${QUIET:-0}
if [[ "$quiet" != 0 && "$quiet" != 1 ]]; then
  printf 'QUIET must be 0 or 1.\n' >&2
  exit 2
fi

convergence_report=${CONVERGENCE_REPORT:-1}
if [[ "$convergence_report" != 0 && "$convergence_report" != 1 ]]; then
  printf 'CONVERGENCE_REPORT must be 0 or 1.\n' >&2
  exit 2
fi

steps=${1:-1}
if [[ ! "$steps" =~ ^[1-9][0-9]*$ ]]; then
  printf 'Usage: bash tests/run_emi_rizzo.sh [positive number of steps]\n' >&2
  exit 2
fi
read -r -a rank_counts <<< "${RANKS:-1 8}"
if (( ${#rank_counts[@]} == 0 )); then
  printf 'RANKS must contain at least one positive rank count.\n' >&2
  exit 2
fi
for ranks in "${rank_counts[@]}"; do
  if [[ ! "$ranks" =~ ^[1-9][0-9]*$ ]]; then
    printf 'Invalid MPI rank count: %s\n' "$ranks" >&2
    exit 2
  fi
done
vtk=${VTK:-0}
if [[ "$vtk" != 0 && "$vtk" != 1 ]]; then
  printf 'VTK must be 0 or 1.\n' >&2
  exit 2
fi
mesh=../input_emi_mesh/rizzo_combi_domi
for input in "$mesh.vtu" "${mesh}_extracellular.txt" \
             "${mesh}_intracellular.txt" "${mesh}_excited.txt"; do
  if [[ ! -r "$input" ]]; then
    printf 'Missing Rizzo input: %s\n' "$input" >&2
    exit 1
  fi
done
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
final_time=$(python3 -c 'import sys; print(format(int(sys.argv[1])*0.01, ".17g"))' "$steps")
for ranks in "${rank_counts[@]}"; do
  out="output/rizzo-bddc-optimized-cond-np${ranks}-${steps}steps"
  mkdir -p "$out"
  "${MPIEXEC:-mpirun}" --allow-run-as-root -np "$ranks" \
    -x OPENBLAS_NUM_THREADS -x OMP_NUM_THREADS \
    ./emiModel \
    --mpi 1 --sdc 0 \
    --input "$mesh.vtu" \
    --extra_set "${mesh}_extracellular.txt" \
    --intra_set "${mesh}_intracellular.txt" \
    --excited "${mesh}_excited.txt" \
    --sigma_i 3 --sigma_e 20 --C_m 1 --R 0.1 \
    --dt 0.01 --finalTime "$final_time" --maximumNumberOfTimeSteps "$steps" \
    --nThreads 1 --profile 1 --diagnostics 0 --vtk "$vtk" --quiet "$quiet" \
    --bddc 1 --direct 0 --bddcCompression 0 --masterWorker 0 \
    --bddcSplitInterfaces 0 --bddcInterfaceTypes 7 \
    --bddcUseCg 1 --bddcRelativeTolerance 1e-7 --bddcVerbose 1 \
    --bddcConditionEstimate 1 --bddcConvergenceReport "$convergence_report" \
    --bddcInterfaceAssembly 2 --bddcInterfaceAssemblyCheck 0 \
    --bddcMaterialAssembly 1 --bddcMaterialAssemblyCheck 0 \
    --dir "$out" 2>&1 | tee "$out/run.log"
done
