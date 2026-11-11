Coarse MPI batching experiment
=============================

Scope: scalar coarse-matrix and coarse-RHS sums are packed in deterministic
storage order and reduced with MPI_SUM in contiguous double buffers. Counts
larger than INT_MAX are chunked. Structural zeros are retained. No matrix,
constraint, solver, tolerance, splitting or codec settings are changed.
The MPI summation order can change; floating-point/trajectory differences must
be checked. This is not a claim of bitwise equivalence or measured speedup.

The default executable emiModel is not rebuilt by this workflow. The original
timing executable and sources were preserved before editing in:
output/emi-coarse-batch-baseline-20261108T213215Z/
Each new validation job also copies the baseline executable into its own output.

From the EMI root on the cluster:
mkdir -p output
emi_coarse_check_job=$(sbatch --parsable tests/submit_emi_coarse_batch_check.sbatch)
printf '%s\n' "$emi_coarse_check_job" | tee output/emi-coarse-batch-check-job.txt
tail -n 60 -F "output/emi-coarse-batch-check-${emi_coarse_check_job}.log"

The job requests one exclusive eight-node Gold6342 allocation. It compiles the
production EMI translation unit into output/emi-coarse-batch-JOB/build/emiModel,
runs a four-rank/two-node MPI helper test and validates weak_32/64/128_h16 on
2/4/8 nodes (two ranks per node). Each layout runs direct, uncompressed FCG,
FCG+split before 12 bits, and FCG+split before 16 bits. Snapshots are enabled.
No performance repeats run in this first job.

After successful completion:
emi_coarse_check_job=$(cat output/emi-coarse-batch-check-job.txt)
cat "output/emi-coarse-batch-${emi_coarse_check_job}/comparison.txt"

The comparison checks inputs/layout, convergence and complete finite snapshots,
and reports raw/centered direct-reference errors, changed solution norms,
iterations and residuals against the archived scalar validation. It imposes no
scientific solution-error acceptance threshold. Review changes before timing.
The full comparison.json also records direct-solution changes and payloads.

Then submit the matched timing campaign:
emi_coarse_timing_job=$(sbatch --parsable tests/submit_emi_coarse_batch_timings.sbatch "$emi_coarse_check_job")
printf '%s\n' "$emi_coarse_timing_job" | tee output/emi-coarse-batch-timing-job.txt
tail -n 60 -F "output/emi-weak-timings-${emi_coarse_timing_job}.log"

Timing uses the candidate executable from its new validation manifests, checks
its SHA256, and retains the existing one warm-up/three measured protocol,
Gold6342 CPU check, two ranks/node, snapshots/histograms/VTK off, profiles and
true residual checks on. It rejects timing/validation trajectory mismatches.
Compare new medians/spreads with scalar timing job 6009538; node pool differs,
so small time differences alone are not proof of improvement.

Local checks (MPI launcher must be permitted to open local sockets):
mpicxx -std=c++20 -O2 -Wall -Wextra -Werror -I "$KASKADE7" tests/test_bddc_coarse_reduction.cpp -o /tmp/coarse-reduction-test
mpirun -n 2 /tmp/coarse-reduction-test
mpirun -n 4 /tmp/coarse-reduction-test
python3 -m unittest discover -s tests -p 'test_emi_solver_study.py'
python3 -m unittest discover -s tests -p 'test_emi_weak_timings.py'

The MPI helper test exercises production packing/reduction against independent
scalar PMPI sums, checks collective counts and forced chunk boundaries, and
covers empty systems, structural zeros, inactive owners and repeated RHS data.
The full application build and multi-node validation require the cluster image;
local helper and runner checks do not replace them.

New campaigns use 12/16 bits. Existing 12/24-bit references remain supported.
The archived scalar baseline has no 16-bit case: its old/new fields are marked
unavailable. The new 16-bit case is checked against its own direct reference
and raw FCG control; no scalar-vs-batched 16-bit speedup can be claimed.

Reuse after a build-only failure:
sbatch tests/submit_emi_coarse_batch_check.sbatch --reuse-build "$PWD/output/emi-coarse-batch-6009589/build"
The reused build is read in place. No recursive archive copy or permission/ACL
preservation is attempted. Verification hashes and new validation outputs go
into the new job directory. Keep the original build for subsequent timings.
