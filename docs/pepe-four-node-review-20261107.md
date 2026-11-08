# Pepe four-node validation: 7 November 2026

Reference: `output/emi-mpi-routing-pepe-20261107T182241Z-76236`. All 16 cases passed. Eight ranks on four nodes, two ranks per node; allocated hosts htc-cmp301,302,307,308. Sixteen-bit transform-first compression, constant splitting in compressed cases; histograms off. Transport self,vader,tcp; Vader single-copy none.

The following destination-routed cases are the reference for repeated timings. MB is 1,000,000 bytes. Local bytes combine same-rank and same-node routes. Inter-node savings compare total payload with the matched separately converged raw destination solve.

| Case | Iterations | Solver s | Application s | Local MB | Inter-node MB | Inter-node saved | Raw relative error/direct |
|---|---:|---:|---:|---:|---:|---:|---:|
| cg-uncompressed-destination | 41 | 19.841 | 68.040 | 30.948768 | 21.924832 | 0.00% | 7.1091e-10 |
| cg-compressed-split-destination | 45 | 20.377 | 67.280 | 8.821944 | 6.595072 | 69.92% | 4.1407e-06 |
| fcg-uncompressed-destination | 41 | 18.232 | 66.360 | 30.948768 | 21.924832 | 0.00% | 7.1253e-10 |
| fcg-compressed-split-destination | 45 | 18.861 | 66.920 | 8.825744 | 6.598152 | 69.91% | 4.1330e-06 |
| fgmres-uncompressed-destination | 43 | 18.739 | 69.780 | 32.458464 | 22.994336 | 0.00% | 1.0935e-09 |
| fgmres-compressed-split-destination | 43 | 20.229 | 70.930 | 8.426040 | 6.299640 | 72.60% | 4.8053e-05 |
| cg-compressed-split-local-bypass | 43 | 19.435 | 71.540 | 32.458464 | 11.643944 | 46.89% | 1.7076e-06 |
| fcg-compressed-split-local-bypass | 43 | 19.278 | 69.710 | 32.458464 | 11.647704 | 46.87% | 2.0780e-06 |
| fgmres-compressed-split-local-bypass | 43 | 21.234 | 70.490 | 32.458464 | 11.637384 | 49.39% | 2.7134e-06 |

These single validation timings include snapshot diagnostics. The two-node measurements used different hosts (htc-cmp011/012). Do not attribute the timing difference solely to node/rank count without controlling hardware and collecting repeated timings.

Local bypass FCG error/direct increased from approximately 8.98e-8 in the two-node reference to 2.08e-6 here. The local-versus-remote ownership classification and numerical reduction order can change with placement. Passing the residual check does not imply the solution errors are identical or scientifically acceptable; no solution-error threshold has been agreed.

## Repeat timings

The timing runner now inherits rank/node counts from the reference. It retains the executable hash, iteration and payload consistency checks and disables only potential snapshot diagnostics. No solver rebuild is needed.

From a fresh login shell, in the repository root, submit:

```bash
sbatch --nodes=4 --ntasks=8 --ntasks-per-node=2 \
  tests/submit_emi_routing_timings.sbatch \
  output/emi-mpi-routing-pepe-20261107T182241Z-76236
```

The batch script requests three hours. Its reference argument is required here: the default reference is the older two-node run. Release the completed interactive allocation before requesting these resources; do not nest another salloc. Preserve the new node list, since the scheduler may choose other hosts.

Alternatively, in the existing four-node allocation if sufficient time remains:

```bash
python3 tests/run_emi_routing_timings.py \
  --reference output/emi-mpi-routing-pepe-20261107T182241Z-76236 \
  --slurm-container "$emi_sif" --kaskade "$emi_kaskade" \
  --warmups 1 --repetitions 3
```

Both paths perform 36 solves. Run only one path. Inspect summary.json and timings.md after completion; measured repetitions exclude warm-ups.
