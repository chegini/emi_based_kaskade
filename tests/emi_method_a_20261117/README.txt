Pepe Method A depth/splitting comparison, prepared 2026-11-17

Default submission: 12 bits, 2 physical nodes, 4 MPI ranks, 3 timesteps.
Nine trajectories: matched raw FCG, then Method A at depths 1,2,3,4 with
prolongation splitting OFF and ON. Both R/P use the specified bit depth and
lifting cap. Exact corners; full coding/packet packing; corrected weight-aware
splitting; P2 baseline ownership; true tolerance 1e-7. Method A quantizes first
and uses reversible integer lifting; it has no separate recentering switch.
Validated lifting optimization mask A=5 (raw=0). No solver code modifications.

Submit from emiModelBddcSdc_new_scratch, after these files are visible on HTC:

out=/scratch/htc/$USER/emi-method-a
mkdir -p "$out"
sbatch --nodes=2 --ntasks=4 --ntasks-per-node=2 \
  --time=01:00:00 --job-name=pepe-A12 \
  --output="$out/method-A-%j.log" --error="$out/method-A-%j.log" \
  tests/emi_method_a_20261117/submit_check.sbatch 12

The wrapper creates a fresh isolated build and runs native checks before the
trajectories. It never submits additional jobs. Existing B data are not overwritten.
The one-hour limit is an allocation limit, not a prediction of queue/run duration.
Optional 16-bit matched A comparison uses the same script with argument 16;
it is not required for the initial 12-bit investigation.

Results: /scratch/htc/$USER/emi-method-a/check-JOB_ID/results/
- report.txt: configuration status
- summary.csv: iterations, raw-relative errors, bytes, savings, bytes/iteration,
  diagnostic application/FCG timings, entropy, and integer-depth trajectory parity
- per-step.csv: actual iterations/residual/error and separate R/P bytes/applications
- serialized-message-byte-entropy.csv: complete serialized send-buffer entropy,
  pooled across ranks and messages, separately R/P and communication scope
- per-case hb-time CSVs, stage timers, partition transport CSVs, state snapshots,
  convergence CSVs and run.log: retained for later A/B analysis
- manifest/source hashes and allocation metadata: provenance

All timing is diagnostic and observer-instrumented, with no repeated timing runs.
Do not infer speedup. Compare inter_node with inter_node; inter_rank includes
same_node and inter_node. Do not sum overlapping scopes or rank-max timers.
Failed trajectories are FAILED; absent raw reference causes compressed runs to
be SKIPPED. Existing raw snapshots are not an independent exact-solution error.
Depth>1 state and residual histories are compared exactly against depth 1 for
the same splitting mode; a discrepancy is a validation failure, not suppressed.

Offline preparation:
PYTHONDONTWRITEBYTECODE=1 python3 tests/emi_method_a_20261117/run_check.py --bits 12
PYTHONDONTWRITEBYTECODE=1 python3 -m unittest discover \
  -s tests/emi_method_a_20261117 -p test_check.py
No solver execution occurs without --execute and a valid SLURM allocation.
