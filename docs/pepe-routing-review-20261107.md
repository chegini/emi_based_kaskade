# Pepe: two-node routing and local-bypass review

Source: `output/emi-mpi-routing-pepe-20261107T120921Z-2561102`. All 16 validation cases passed.

Configuration: 4 ranks on 2 nodes (htc-cmp011/012), 20 subdomains, 615043 cells, 120839 DOFs; 16-bit transform-first compression with constant splitting; one physical step; histograms off; Open MPI self,vader,tcp with Vader single-copy none.

These are single validation runs with solution snapshots enabled. They do not establish a speedup. Solver time is the maximum solver-wrapper time across ranks, including setup. Application time is the reported application wall time.

## Iterations and timings

| Case | Iterations | True relative residual | Solver s | Application s | Solver change vs matched raw destination |
|---|---:|---:|---:|---:|---:|
| cg-uncompressed-legacy | 41 | 7.104e-08 | 34.073 | 100.000 | -1.36% |
| cg-uncompressed-destination | 41 | 7.104e-08 | 34.542 | 99.500 | +0.00% |
| cg-compressed-split-legacy | 44 | 9.376e-08 | 37.060 | 103.460 | +7.29% |
| cg-compressed-split-destination | 44 | 9.376e-08 | 37.270 | 103.920 | +7.90% |
| fcg-uncompressed-legacy | 41 | 7.104e-08 | 34.776 | 101.250 | -0.94% |
| fcg-uncompressed-destination | 41 | 7.104e-08 | 35.104 | 102.530 | +0.00% |
| fcg-compressed-split-legacy | 45 | 7.766e-08 | 36.927 | 103.200 | +5.19% |
| fcg-compressed-split-destination | 45 | 7.766e-08 | 37.183 | 102.560 | +5.92% |
| fgmres-uncompressed-legacy | 43 | 8.345e-08 | 39.831 | 107.000 | +0.60% |
| fgmres-uncompressed-destination | 43 | 8.345e-08 | 39.595 | 107.420 | +0.00% |
| fgmres-compressed-split-legacy | 43 | 8.356e-08 | 40.694 | 107.860 | +2.77% |
| fgmres-compressed-split-destination | 43 | 8.356e-08 | 44.078 | 112.050 | +11.32% |
| cg-compressed-split-local-bypass | 42 | 8.629e-08 | 39.376 | 116.110 | +13.99% |
| fcg-compressed-split-local-bypass | 42 | 8.580e-08 | 37.719 | 114.170 | +7.45% |
| fgmres-compressed-split-local-bypass | 43 | 8.345e-08 | 42.168 | 114.150 | +6.50% |

## Logical payloads

MB means 1,000,000 bytes. Same-node excludes same-rank in the table; both are local. Restriction and prolongation are summed. Savings use the total bytes of the separately converged matched raw destination solve, including any iteration-count changes. These are not hardware network counters.

| Case | Same-rank MB | Same-node MB | Inter-node MB | Total MB | Inter-node saved | Total saved |
|---|---:|---:|---:|---:|---:|---:|
| cg-uncompressed-destination | 10.124704 | 29.799456 | 12.949440 | 52.873600 | 0.00% | 0.00% |
| cg-compressed-split-destination | 2.886008 | 8.345864 | 3.841328 | 15.073200 | 70.34% | 71.49% |
| fcg-uncompressed-destination | 10.124704 | 29.799456 | 12.949440 | 52.873600 | 0.00% | 0.00% |
| fcg-compressed-split-destination | 2.950000 | 8.534976 | 3.927040 | 15.412016 | 69.67% | 70.85% |
| fgmres-uncompressed-destination | 10.618592 | 31.253088 | 13.581120 | 55.452800 | 0.00% | 0.00% |
| fgmres-compressed-split-destination | 2.818752 | 8.153520 | 3.753432 | 14.725704 | 72.36% | 73.44% |
| cg-compressed-split-local-bypass | 10.371648 | 30.526272 | 7.529480 | 48.427400 | 41.85% | 8.41% |
| fcg-compressed-split-local-bypass | 10.371648 | 30.526272 | 7.529872 | 48.427792 | 41.85% | 8.41% |
| fgmres-compressed-split-local-bypass | 10.618592 | 31.253088 | 7.714280 | 49.585960 | 43.20% | 10.58% |

## Solution errors

Raw and centered relative Euclidean errors are both retained. Centering removes the arithmetic mean error; it is a diagnostic, not an automatic scientific justification for ignoring offsets. No solution-error acceptance threshold has been agreed.

| Case | Raw relative L2/direct | Centered relative L2/direct | Max absolute/direct | Mean error/direct | Raw relative L2/raw BDDC |
|---|---:|---:|---:|---:|---:|
| cg-uncompressed-legacy | 7.0707e-10 | 7.0104e-10 | 1.9036e-09 | 6.1209e-12 | -- |
| cg-uncompressed-destination | 7.0707e-10 | 7.0104e-10 | 1.9036e-09 | 6.1209e-12 | -- |
| cg-compressed-split-legacy | 5.1709e-06 | 1.1215e-09 | 3.4475e-07 | -3.4331e-07 | 5.1710e-06 |
| cg-compressed-split-destination | 5.1709e-06 | 1.1215e-09 | 3.4475e-07 | -3.4331e-07 | 5.1710e-06 |
| fcg-uncompressed-legacy | 7.0742e-10 | 7.0104e-10 | 1.9037e-09 | 6.2908e-12 | -- |
| fcg-uncompressed-destination | 7.0742e-10 | 7.0104e-10 | 1.9037e-09 | 6.2908e-12 | -- |
| fcg-compressed-split-legacy | 4.3227e-06 | 9.5715e-10 | 2.8822e-07 | -2.8700e-07 | 4.3228e-06 |
| fcg-compressed-split-destination | 4.3227e-06 | 9.5715e-10 | 2.8822e-07 | -2.8700e-07 | 4.3228e-06 |
| fgmres-uncompressed-legacy | 1.0936e-09 | 1.0889e-09 | 2.3298e-09 | 6.7524e-12 | -- |
| fgmres-uncompressed-destination | 1.0936e-09 | 1.0889e-09 | 2.3298e-09 | 6.7524e-12 | -- |
| fgmres-compressed-split-legacy | 4.3673e-05 | 1.1163e-09 | 2.9014e-06 | 2.8995e-06 | 4.3673e-05 |
| fgmres-compressed-split-destination | 4.3673e-05 | 1.1163e-09 | 2.9014e-06 | 2.8995e-06 | 4.3673e-05 |
| cg-compressed-split-local-bypass | 1.6461e-06 | 9.6226e-10 | 1.1190e-07 | 1.0929e-07 | 1.6461e-06 |
| fcg-compressed-split-local-bypass | 8.9795e-08 | 9.0013e-10 | 8.3915e-09 | 5.9613e-09 | 8.9698e-08 |
| fgmres-compressed-split-local-bypass | 6.6406e-07 | 1.0905e-09 | 4.5953e-08 | 4.4088e-08 | 6.6395e-07 |

## Interpretation

Destination and legacy routing produced identical snapshots, solver traces, iterations and total payloads for the paired configurations. Their one-run timing differences cannot establish a routing speedup.

Local bypass saves 41.85% of inter-node bytes for CG/FCG and 43.20% for FGMRES. Total logical-byte savings are only about 8.41% (CG/FCG) and 10.58% (FGMRES), because local data remains raw and consistency requires some remote DOFs to remain exact.

FCG with bypass has the smallest raw error among bypass cases (8.98e-8 relative to direct). It is 7.45% slower in solver-wrapper time than its raw destination baseline in this single run. This is a candidate for repeated measurements, not a demonstrated performance improvement.

## Repeated timing campaign

Use tests/run_emi_routing_timings.py with this run as --reference. It compares nine destination-routed cases: raw, all-interface compression+split, and local bypass, for each solver. Each case gets one excluded warm-up and three measured repetitions in reproducibly shuffled order (36 solves).

The timing runner disables potentialDiagnostics snapshot output only. It retains physical/solver parameters, residual checks, profiles, byte counters and numerical traces; iteration and payload parity with the reference are checked. Each repetition starts a new process and rebuilds solver setup. Codebook work remains included; no one-time training assumption is made.

It records median/min/max solver, application and launcher times, iterations/residuals and three route-byte categories in summary.json and writes timings.md on completion. Reference solution errors are not presented as newly measured errors in the timing runs.

From the host allocation shell, in the EMI repository, run:

```bash
python3 tests/run_emi_routing_timings.py \
  --reference output/emi-mpi-routing-pepe-20261107T120921Z-2561102 \
  --slurm-container "$emi_sif" \
  --kaskade "$emi_kaskade" \
  --warmups 1 --repetitions 3
```

Use exactly two allocated nodes and four MPI tasks, two per node, with one CPU
per task and 24 GB per node. Allow roughly 1–2 hours based on the validation
durations; request three hours for margin if creating a fresh allocation.
The runner uses the reference precision and transform order and explicitly
sets `self,vader,tcp` and Vader single-copy `none`. Its executable hash must
match the reference; if the binary changed, perform a fresh full validation
and use that new reference. Do not run another EMI campaign concurrently on
the same allocated CPUs. Preserve the entire output directory, including
warm-up logs and environment.json. Add `--plan-only` to inspect the schedule
without launching jobs.

For unattended execution, use a fresh SSH login shell (outside any interactive
allocation), enter the EMI repository root, and submit:

```bash
mkdir -p output
sbatch tests/submit_emi_routing_timings.sbatch
```

The batch script requests its own two-node, three-hour allocation and records
its console output in `output/emi-routing-timings-JOB_ID.log`. Do not wrap it in
`salloc`. It defaults to the reference above; pass another validated reference
directory as its first argument if needed. Export `emi_sif` and `emi_kaskade`
before submission to override the supplied cluster paths.

The first timing attempt, `emi-routing-timings-pepe-20261107T125605Z-2784271`,
was interrupted after five successful warm-ups, during the sixth. No measured
repetitions completed. Its manifest records job 5996341; the terminal reports
expiry of the older job 5996261 followed by Hangup, consistent with a nested
interactive allocation being terminated through its parent shell. Preserve
that incomplete campaign and start a fresh batch campaign; do not combine
its warm-up samples with measured results.
