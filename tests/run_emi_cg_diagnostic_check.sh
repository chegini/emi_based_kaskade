#!/usr/bin/env bash
# Validation only: unchanged 10-cell physical problem, one step, tolerance 1e-7.
set -euo pipefail
repo=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
cd "$repo"
audit_out="${OUTPUT:-output/emi-cg-audit-$(date -u +%Y%m%dT%H%M%SZ)-$$}"
mkdir -p "$audit_out"
common=(--ranks "${RANKS:-8}" --refine "${REFINE:-2}" --tolerance 1e-7)
runner=tests/run_emi_quantization_diagnostic.py

# First verify that observing a convergent run does not change its trajectory.
python3 "$runner" "${common[@]}" --bits 16 --restriction-passes 1 \
  --output "$audit_out/b16-off"
python3 "$runner" "${common[@]}" --bits 16 --restriction-passes 1 \
  --cg-diagnostics --output "$audit_out/b16-on"
python3 - "$audit_out" <<'PY'
from pathlib import Path
import sys
root=Path(sys.argv[1])
for case in ('uncompressed-cg','quantization-16'):
    name='bddc-convergence-step1.csv'
    a=(root/'b16-off'/case/name).read_bytes()
    b=(root/'b16-on'/case/name).read_bytes()
    if a!=b:
        raise SystemExit(f'STOP: {case} convergence history changed with diagnostics; inspect both runs.')
print('Observer parity passed: uncompressed and 16-bit convergence CSVs are identical.')
(root/'observer-parity.txt').write_text('PASS: uncompressed and 16-bit convergence CSVs are identical.\n')
PY

# A numerical failure here is evidence to retain, not a reason to loosen tolerance.
if python3 "$runner" "${common[@]}" --bits 12 --restriction-passes 1 \
    --cg-diagnostics --output "$audit_out/b12-passes1"; then
  printf 'Single-pass 12-bit run passed.\n'
else
  printf 'Single-pass 12-bit diagnostic returned failure; checking reference runs before continuing.\n'
fi
python3 - "$audit_out/b12-passes1/summary.json" <<'PY'
import json,sys
rows=json.load(open(sys.argv[1]))
if len(rows)!=3 or any(r['status']!='PASS' for r in rows[:2]):
    raise SystemExit('STOP: direct or uncompressed reference failed.')
PY
python3 "$runner" "${common[@]}" --bits 12 --restriction-passes 2 \
  --cg-diagnostics --output "$audit_out/b12-passes2"
printf '\nCG audit results: %s\n' "$audit_out"
printf 'Inspect b16-on, b12-passes1 and b12-passes2 with tests/summarize_emi_cg_diagnostics.py.\n'
