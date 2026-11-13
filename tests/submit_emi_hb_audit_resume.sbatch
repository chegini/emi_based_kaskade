#!/usr/bin/env bash
#SBATCH --partition=big
#SBATCH --constraint=Gold6342
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --job-name=emi-hb-resume
#SBATCH --output=output/emi-hb-resume-%j.log
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?Submit from EMI root}"
export LANG=C LC_ALL=C PYTHONDONTWRITEBYTECODE=1
root="$PWD/output/emi-hb-audit-6144259"
resume="$PWD/output/emi-hb-resume-$SLURM_JOB_ID"
sif=/data/numerik/people/fchegini/containers/kaskade-robin-20261105/kaskade7-robin.sif
kaskade=/data/numerik/people/fchegini/kaskade7_new
# Validate completion from the original log, not merely from existing CSV files.
python3 - <<'PY'
from pathlib import Path
import re
log=Path('output/emi-hb-audit-6144259.log').read_text()
completed=set();current=None
for line in log.splitlines():
    m=re.fullmatch(r'Starting bits=(\d+) depth=(\d+)',line)
    if m:current=tuple(map(int,m.groups()))
    if line.startswith('Frozen replay passed:'):
        assert current is not None
        completed.add(current)
expected={(b,d) for b in [8,16] for d in [1,2,3,0]}
assert expected-completed=={(16,0)},('Unexpected unfinished cases',expected-completed)
for b,d in completed:
    folder=Path(f'output/emi-hb-audit-6144259/replay-b{b}-d{d}')
    assert all((folder/f).is_file() for f in ['packets.csv','level-stats.csv','codebooks.csv'])
print('Verified seven completed replays; resuming only 16-bit recursive HB.',flush=True)
PY
sha256sum -c "$root/source-sha256.txt"
mkdir "$resume"
sha256sum "$root/replay" > "$resume/replay-sha256.txt"
srun --mpi=none --nodes=1 --ntasks=1 --ntasks-per-node=1 \
 apptainer exec --cleanenv --bind "$kaskade:/kaskade7" --bind "$PWD:$PWD" \
 --pwd "$PWD" "$sif" bash -c '
 set -euo pipefail
 export PATH=/usr/lib64/mpi/gcc/openmpi4/bin:$PATH
 export LD_LIBRARY_PATH=/usr/lib64/mpi/gcc/openmpi4/lib64:/usr/lib64/mpi/gcc/openmpi4/lib
 "$1/replay" "$1/inputs" "$2/replay-b16-d0" 16 0
 ' audit "$root" "$resume"
sha256sum -c "$resume/replay-sha256.txt"
sha256sum -c "$root/source-sha256.txt"
# Keep the timed-out files; install the replacement only after all checks pass.
mkdir -p "$root/incomplete"
if [[ -d "$root/replay-b16-d0" ]]; then
 mv "$root/replay-b16-d0" "$root/incomplete/replay-b16-d0-$SLURM_JOB_ID"
fi
mv "$resume/replay-b16-d0" "$root/replay-b16-d0"
python3 tests/summarize_emi_hb_audit.py "$root"
echo "All eight replays complete. Summary tables: $root"
