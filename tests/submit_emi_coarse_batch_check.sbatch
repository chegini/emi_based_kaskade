#!/usr/bin/env bash
# One CPU family and one node pool for all three validation layouts.
#SBATCH --partition=big
#SBATCH --constraint=Gold6342
#SBATCH --nodes=8
#SBATCH --ntasks=16
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --mem=64G
#SBATCH --time=03:00:00
#SBATCH --job-name=lego-coarse-check
#SBATCH --output=output/emi-coarse-batch-check-%j.log
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?Submit from the EMI repository root}"
export LANG=C LC_ALL=C LC_CTYPE=C PYTHONDONTWRITEBYTECODE=1
emi_sif=${emi_sif:-/data/numerik/people/fchegini/containers/kaskade-robin-20261105/kaskade7-robin.sif}
emi_kaskade=${emi_kaskade:-/data/numerik/people/fchegini/kaskade7_new}
emi_reuse_build=
if (( $# )); then
  if (( $# != 2 )) || [[ "$1" != --reuse-build ]]; then
    echo 'Usage: sbatch tests/submit_emi_coarse_batch_check.sbatch [--reuse-build ABSOLUTE_BUILD_DIR]' >&2
    exit 2
  fi
  emi_reuse_build=$2
  [[ "$emi_reuse_build" == /* ]] || { echo 'Reuse build path must be absolute' >&2; exit 2; }
fi
emi_root="$PWD/output/emi-coarse-batch-$SLURM_JOB_ID"
mkdir "$emi_root"
python3 - <<'PY'
import hashlib, json
from pathlib import Path
reference = Path('output/emi-solver-study-weak_32_h16-before-20261108T203929Z-2844236')
expected = json.loads((reference/'environment.json').read_text())['executable_sha256']
actual = hashlib.sha256(Path('emiModel').read_bytes()).hexdigest()
if actual != expected:
    raise SystemExit('Root emiModel differs from the recorded scalar baseline; inspect before building')
PY
cp -- emiModel "$emi_root/emiModel-baseline"
sha256sum emiModel > "$emi_root/baseline-sha256.txt"
printf 'Coarse-batch validation job %s on %s\nOutput: %s\n' \
  "$SLURM_JOB_ID" "$SLURM_JOB_NODELIST" "$emi_root"
emi_build="$emi_root/build"
if [[ -n "$emi_reuse_build" ]]; then
  emi_build="$emi_reuse_build"
  printf '%s\n' "$emi_build" > "$emi_root/reused-build.txt"
else
  srun --mpi=none --nodes=1 --ntasks=1 --ntasks-per-node=1 \
    apptainer exec --cleanenv --bind "$emi_kaskade:/kaskade7" \
    --bind "$PWD:$PWD" --pwd "$PWD" "$emi_sif" \
    /bin/bash tests/build_emi_coarse_batch_cluster.sh "$emi_build"
fi
# Reuse reads the original build directly: no copy, chmod, ACL or xattr changes.
python3 tests/verify_emi_coarse_build.py \
  --build "$emi_build" --repo "$PWD" --kaskade "$emi_kaskade" \
  --manifest "$emi_root/verified-build.json" --write-manifest
srun --mpi=pmix_v5 --nodes=2 --ntasks=4 --ntasks-per-node=2 --cpus-per-task=1 \
  --cpu-bind=cores --distribution=block:block --kill-on-bad-exit=1 --time=00:05:00 \
  apptainer exec --bind "$PWD:$PWD" --pwd "$PWD" "$emi_sif" \
  /bin/bash --noprofile --norc -c '
export PATH=/usr/lib64/mpi/gcc/openmpi4/bin:/usr/local/bin:/usr/bin:/bin
export LD_LIBRARY_PATH=/usr/lib64/mpi/gcc/openmpi4/lib64:/usr/lib64/mpi/gcc/openmpi4/lib
export OMPI_MCA_pml=ob1 OMPI_MCA_btl=self,vader,tcp OMPI_MCA_btl_vader_single_copy_mechanism=none
exec "$1"
' coarse-test "$emi_build/coarse-reduction-test"
emi_verify_sources()
{
  python3 tests/verify_emi_coarse_build.py \
    --build "$emi_build" --repo "$PWD" --kaskade "$emi_kaskade" \
    --manifest "$emi_root/verified-build.json"
}
for emi_size in 32 64 128; do
  emi_verify_sources
  emi_nodes=$((emi_size / 16))
  python3 -u tests/run_emi_solver_study.py \
    --case "weak_${emi_size}_h16" --solvers fcg --split on --bits 12 16 \
    --transform-order before --iterations 1000 --case-minutes 20 \
    --ranks "$((2 * emi_nodes))" --expected-nodes "$emi_nodes" \
    --executable "$emi_build/emiModel" \
    --slurm-container "$emi_sif" --kaskade "$emi_kaskade" \
    --output "$emi_root/validation-weak_${emi_size}_h16"
done
emi_verify_sources
python3 tests/compare_emi_coarse_batch.py --campaign "$emi_root" \
  | tee "$emi_root/comparison.txt"
cmp -- emiModel "$emi_root/emiModel-baseline"
touch "$emi_root/VALIDATION_COMPLETE"
echo 'Validation complete. Review comparison.txt before submitting timing repeats.'
