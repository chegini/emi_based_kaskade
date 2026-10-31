#!/usr/bin/env bash
# Run from the EMI repository inside the container. No SSH is used.
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
bits=${BITS:-24}
ranks=${RANKS:-1}
tolerance=${TOLERANCE:-1e-7}
use_cg=${USE_CG:-1}
case "$use_cg" in
  0|1) ;;
  *) printf 'USE_CG must be 0 or 1\n' >&2; exit 2 ;;
esac
# The detailed convergence CSV currently requires CG. Both solver variants
# still use and validate the same true global residual stopping criterion.
convergence_report=$use_cg
python3 - "$tolerance" <<'CHECK_TOLERANCE'
import math
import sys
try:
    valid = math.isfinite(float(sys.argv[1])) and float(sys.argv[1]) > 0
except ValueError:
    valid = False
if not valid:
    raise SystemExit("TOLERANCE must be finite and positive")
CHECK_TOLERANCE
base=./input/10Cells3d_10extra
root=output/10cells-compression-check-b${bits}-np${ranks}-cg${use_cg}-tol${tolerance}-$(date -u +%Y%m%dT%H%M%SZ)-$$
for file in "${base}_mesh.vtu" "${base}_list_extracellular.txt" \
            "${base}_list_intracellular.txt" "${base}_early_excited.txt"; do
  if [[ ! -r "$file" ]]; then printf 'Missing input: %s\n' "$file" >&2; exit 2; fi
done
mkdir -p "$root"
printf 'Output: %s\nRelative tolerance: %s\nUse CG: %s\n' "$root" "$tolerance" "$use_cg"
failed=0
for stage in uncompressed quantization full; do
  compression=1; graph=0; coding=0
  if [[ "$stage" == uncompressed ]]; then compression=0; fi
  if [[ "$stage" == full ]]; then graph=1; coding=1; fi
  out="$root/$stage"
  mkdir -p "$out"
  command=("${MPIEXEC:-mpirun}" --allow-run-as-root -np "$ranks"
    -x OPENBLAS_NUM_THREADS -x OMP_NUM_THREADS ./emiModel
    --mpi 1 --sdc 0 --bddc 1 --direct 0 --refine 0
    --input "${base}_mesh.vtu" --extra_set "${base}_list_extracellular.txt"
    --intra_set "${base}_list_intracellular.txt" --excited "${base}_early_excited.txt"
    --sigma_i 3 --sigma_e 20 --C_m 1 --R 0.1
    --dt 0.01 --finalTime 0.01 --maximumNumberOfTimeSteps 1
    --nThreads 1 --quiet 0 --vtk 1 --diagnostics 0 --profile 1
    --masterWorker 0 --bddcSplitInterfaces 0 --bddcInterfaceTypes 7
    --bddcUseCg "$use_cg" --bddcRelativeTolerance "$tolerance" --bddcIterations 300
    --bddcVerbose 0 --bddcConvergenceReport "$convergence_report" --bddcConditionEstimate 0
    --bddcMaterialAssembly 1 --bddcInterfaceAssembly 2
    --bddcCompression "$compression" --bddcCompressionBits "$bits"
    --bddcGraphLifting "$graph" --bddcGraphLiftingAfterQuantization "$graph"
    --bddcZigzag "$coding" --bddcBitlength "$coding" --bddcHuffman "$coding"
    --bddcExcludeCornersFromCompression 1 --bddcCompressionMinInterfaceSize 0
    --bddcCompressionThresholdInterfaceTypes 4 --bddcCompressionReport 1
    --bddcHistogram 0 --dir "$out")
  printf '%q ' "${command[@]}" > "$out/command.sh"
  printf '\n' >> "$out/command.sh"
  printf '\nRunning %s: bits=%s ranks=%s\nLog: %s/run.log\n' "$stage" "$bits" "$ranks" "$out"
  if "${command[@]}" > "$out/run.log" 2>&1; then rc=0; else rc=$?; fi
  if python3 tests/collect_emi_scaling_smoke.py "$out/run.log" \
       --steps 1 --ranks "$ranks" --tolerance "$tolerance" --exit-code "$rc"; then
    if [[ "$stage" != uncompressed ]]; then
      python3 compare_emi_vtu.py "$out/emiBDDCLastCompression.vtu" \
        "$root/uncompressed/emiBDDCLast.vtu" > "$out/comparison.log"
      grep -E '^(same coordinates|relative L2 error|max absolute error):' "$out/comparison.log"
    fi
  else
    failed=$((failed+1))
    grep -E 'BDDC global residual|failed after|target not reached' "$out/run.log" | head -6 || true
    if [[ "$stage" == uncompressed ]]; then
      printf 'Uncompressed baseline failed; stopping.\n' >&2
      exit 1
    fi
  fi
done
printf '\nFinished with %s failed stage(s). Logs and any convergence CSVs: %s\n' "$failed" "$root"
if (( failed>0 )); then exit 1; fi
