#!/usr/bin/env bash
# Global reference solve, manually run inside the MPI-enabled container.
set -euo pipefail
repo_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
cd "$repo_dir"
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
mesh=../input_emi_mesh/rizzo_combi_domi
out=output/rizzo-direct-mumps-spd-10steps
for input in "$mesh.vtu" "${mesh}_extracellular.txt" \
             "${mesh}_intracellular.txt" "${mesh}_excited.txt"; do
  if [[ ! -r "$input" ]]; then
    printf 'Missing Rizzo input: %s\n' "$input" >&2
    exit 1
  fi
done
mkdir -p "$out"
timing=()
if [[ -x /usr/bin/time ]]; then timing=(/usr/bin/time -v); fi
"${MPIEXEC:-mpirun}" --allow-run-as-root -np 1 \
  -x OPENBLAS_NUM_THREADS -x OMP_NUM_THREADS \
  "${timing[@]}" ./emiModel \
  --mpi 1 --sdc 0 --bddc 0 --direct 1 --directBackend mumps-spd \
  --input "$mesh.vtu" \
  --extra_set "${mesh}_extracellular.txt" \
  --intra_set "${mesh}_intracellular.txt" \
  --excited "${mesh}_excited.txt" \
  --sigma_i 3 --sigma_e 20 --C_m 1 --R 0.1 \
  --dt 0.01 --finalTime 0.1 --maximumNumberOfTimeSteps 10 \
  --nThreads 1 --vtk 1 --diagnostics 0 --profile 1 \
  --dir "$out" 2>&1 | tee "$out/run.log"

bddc=output/rizzo-bddc-optimized-cond-np8-10steps/emiBDDCLast.vtu
if [[ -f "$bddc" ]]; then
  python3 compare_emi_vtu.py "$bddc" "$out/emiLast.vtu" \
    | tee "$out/compare-bddc.log"
fi
