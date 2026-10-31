#!/usr/bin/env bash
# Container validation before timing/scaling runs. No Slurm jobs are submitted.
set -euo pipefail
repo_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
cd "$repo_dir"
ranks=${RANKS:-8}
steps=${STEPS:-1}
bits=${BITS:-12}
max_iterations=${MAX_ITERATIONS:-300}
read -r -a cases <<< "${CASES:-weak_8_h16 strong_h16}"
for value in "$ranks" "$steps" "$bits" "$max_iterations"; do
  if [[ ! "$value" =~ ^[1-9][0-9]*$ ]]; then
    printf 'RANKS, STEPS, BITS and MAX_ITERATIONS must be positive integers.\n' >&2
    exit 2
  fi
done
if (( ${#cases[@]} == 0 )); then
  printf 'CASES must not be empty.\n' >&2
  exit 2
fi
mesh_dir=${MESH_DIR:-../input_emi_mesh_blocks/weak_strong_scaling}
for case_name in "${cases[@]}"; do
  if [[ ! "$case_name" =~ ^weak_[1-9][0-9]*_h16$ && "$case_name" != strong_h16 ]]; then
    printf 'Unsupported smoke-test case: %s\n' "$case_name" >&2
    exit 2
  fi
  base="$mesh_dir/${case_name}_kaskade"
  for input in "${base}_mesh.vtu" "${base}_list_extracellular.txt" \
               "${base}_list_intracellular.txt" "${base}_early_excited.txt"; do
    if [[ ! -r "$input" ]]; then
      printf 'Missing input: %s\n' "$input" >&2
      exit 2
    fi
  done
done
test -x ./emiModel
out=${OUT_ROOT:-output/scaling-full-smoke-$(date -u +%Y%m%dT%H%M%SZ)-$$}
if [[ -e "$out" ]]; then
  printf 'Refusing to reuse an existing output directory: %s\n' "$out" >&2
  exit 2
fi
mkdir -p "$out"
printf 'Container validation output: %s\n' "$out"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
final_time=$(python3 -c 'import sys; print(format(int(sys.argv[1])*0.01,".17g"))' "$steps")
{
  printf 'purpose=container validation; VTK enabled; not scaling timing data\n'
  printf 'ranks=%s\nsteps=%s\nbits=%s\nmaximum_iterations=%s\n' \
    "$ranks" "$steps" "$bits" "$max_iterations"
  printf 'compression=quantization,integer graph lifting,zigzag,bitlength,Huffman\n'
  printf 'corner_dofs=exact\ntrue_relative_tolerance=1e-7\n'
  python3 - <<'PY'
import hashlib
import os
from pathlib import Path
h=hashlib.sha256()
print("hostname="+os.uname().nodename)
with Path("emiModel").open("rb") as stream:
    for chunk in iter(lambda: stream.read(1024*1024),b""):
        h.update(chunk)
print("executable_sha256="+h.hexdigest())
try:
    print("visible_cpus="+str(len(os.sched_getaffinity(0))))
except (AttributeError,OSError):
    print("visible_cpus=unavailable")
PY
  git rev-parse HEAD || true
  git -C /kaskade7 rev-parse HEAD || true
  "${MPIEXEC:-mpirun}" --version || true
} > "$out/environment.txt" 2>&1
for case_name in "${cases[@]}"; do
  base="$mesh_dir/${case_name}_kaskade"
  for compression in 0 1; do
    mode=uncompressed
    if [[ "$compression" == 1 ]]; then mode="full-bits${bits}"; fi
    run="$out/$case_name/$mode"
    mkdir -p "$run"
    command=("${MPIEXEC:-mpirun}" --allow-run-as-root -np "$ranks"
      -x OPENBLAS_NUM_THREADS -x OMP_NUM_THREADS ./emiModel
      --mpi 1 --sdc 0 --bddc 1 --direct 0 --refine 0
      --input "${base}_mesh.vtu"
      --extra_set "${base}_list_extracellular.txt"
      --intra_set "${base}_list_intracellular.txt"
      --excited "${base}_early_excited.txt"
      --sigma_i 3 --sigma_e 20 --C_m 1 --R 0.1
      --dt 0.01 --finalTime "$final_time" --maximumNumberOfTimeSteps "$steps"
      --nThreads 1 --quiet 0 --vtk 1 --diagnostics 0 --profile 1
      --masterWorker 0 --bddcSplitInterfaces 0 --bddcInterfaceTypes 7
      --bddcUseCg 1 --bddcRelativeTolerance 1e-7 --bddcIterations "$max_iterations"
      --bddcVerbose 0 --bddcConvergenceReport 0 --bddcConditionEstimate 0
      --bddcMaterialAssembly 1 --bddcInterfaceAssembly 2
      --bddcCompression "$compression" --bddcCompressionBits "$bits"
      --bddcGraphLifting 1 --bddcGraphLiftingAfterQuantization 1
      --bddcZigzag 1 --bddcBitlength 1 --bddcHuffman 1
      --bddcExcludeCornersFromCompression 1 --bddcCompressionMinInterfaceSize 0
      --bddcCompressionThresholdInterfaceTypes 4 --bddcCompressionReport 1
      --bddcHistogram 0 --dir "$run")
    printf '%q ' "${command[@]}" > "$run/command.sh"
    printf '\n' >> "$run/command.sh"
    printf 'Running %s, %s, %s rank(s), %s step(s)...\n' "$case_name" "$mode" "$ranks" "$steps"
    printf 'Log: %s\n' "$run/run.log"
    if [[ -r /sys/fs/cgroup/memory.events ]]; then
      cp /sys/fs/cgroup/memory.events "$run/cgroup-memory-events-before.txt"
    fi
    set +e
    "${command[@]}" > "$run/run.log" 2>&1
    rc=$?
    set -e
    if [[ -r /sys/fs/cgroup/memory.events ]]; then
      cp /sys/fs/cgroup/memory.events "$run/cgroup-memory-events-after.txt"
    fi
    if ! python3 tests/collect_emi_scaling_smoke.py "$run/run.log" \
          --steps "$steps" --ranks "$ranks" --exit-code "$rc"; then
      printf 'Stopped. Inspect %s and metrics.json; no tolerance was relaxed.\n' "$run/run.log" >&2
      tail -n 20 "$run/run.log" >&2
      exit 1
    fi
  done
  comparison="$out/$case_name/comparison.log"
  python3 compare_emi_vtu.py \
    "$out/$case_name/full-bits${bits}/emiBDDCLastCompression.vtu" \
    "$out/$case_name/uncompressed/emiBDDCLast.vtu" > "$comparison"
  python3 - "$comparison" <<'PY'
import math
from pathlib import Path
import re
import sys
text=Path(sys.argv[1]).read_text()
for label in ("relative L2 error", "max absolute error"):
    match=re.search(r"^"+re.escape(label)+r":\s+(\S+)\s*$",text,re.MULTILINE)
    if not match or not math.isfinite(float(match.group(1))):
        raise SystemExit("Invalid solution comparison: "+label)
    print(label+": "+match.group(1))
PY
done
printf 'Container checks completed. Review solution differences before the timing campaign.\n'
printf 'Full logs, exact commands and metrics: %s\n' "$out"
