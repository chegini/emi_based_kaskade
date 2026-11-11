#!/usr/bin/env bash
# A real script file keeps MPI stdin forwarding away from the shell program.
set -euo pipefail
export PATH=/usr/lib64/mpi/gcc/openmpi4/bin:/usr/local/bin:/usr/bin:/bin
export LD_LIBRARY_PATH=/usr/lib64/mpi/gcc/openmpi4/lib64:/usr/lib64/mpi/gcc/openmpi4/lib
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
export OMPI_MCA_pml=ob1 OMPI_MCA_btl=self,vader,tcp
export PYTHONDONTWRITEBYTECODE=1
check=$(mktemp -d /tmp/emi-physical.XXXXXX)
mpicxx -std=c++20 -O2 -I/kaskade7 /kaskade7/tests/bddcTransferOrder/physicalHistogram.cpp -o "$check/physical"
"$check/physical"
resume=${2:-}
result="$PWD/output/emi-physical-histogram-${resume:-$1}"
printf 'Campaign results: %s\n' "$result"
if [[ -n "$resume" ]]; then
  [[ "$resume" =~ ^[0-9]+$ ]] || { echo 'Resume argument must be a job number' >&2; exit 2; }
  python3 - "$result" <<'VERIFY'
import hashlib,json,sys
from pathlib import Path
root=Path(sys.argv[1]); previous=root/'onset-1'
s=json.loads((previous/'summary.json').read_text())
e=json.loads((previous/'environment.json').read_text())
assert not s['errors'], 'Previous campaign has errors'
assert s['checks']['observer']=='EXACT_SNAPSHOTS_TRACES_ITERATIONS_AND_BYTES'
assert s['checks']['physical_capture']=='FIRST_CALL_BOTH_PHASES_FINITE_AND_ORDERED'
assert (e['solver'],e['split'],e['bits'],e['case'],e['steps'],e['ranks'])==('fcg','on',16,'pepe',1,2)
assert hashlib.sha256(Path('emiModel').read_bytes()).hexdigest()==e['executable_sha256'], 'Executable changed; submit a fresh campaign'
for name,digest in e['source_sha256'].items():
    assert hashlib.sha256(Path(name).read_bytes()).hexdigest()==digest, 'Source changed: '+name
assert not (root/'onset-6').exists(), 'onset-6 already exists; inspect it before resuming'
print('Verified completed onset-1 and unchanged executable/source hashes',flush=True)
VERIFY
else
  mkdir "$result"
fi
# Existing validation checks exact observer-off/on snapshots, solver traces and bytes.
if [[ -z "$resume" ]]; then
bash tests/run_emi_histogram_start_check.sh \
  --case pepe --ranks 2 --steps 1 --uncompressed-steps 0 --histogram-steps 1 \
  --bits 16 --solver fcg --split on --transform-order before --physical-values \
  --output "$result/onset-1"
fi
python3 tests/run_emi_histogram_start_check.py \
  --case pepe --ranks 2 --steps 6 --uncompressed-steps 5 --histogram-steps 6 \
  --bits 16 --solver fcg --split on --transform-order before --physical-values \
  --output "$result/onset-6"
python3 tests/summarize_emi_physical_histograms.py "$result"
printf 'Physical EMI capture complete: %s\n' "$result"
