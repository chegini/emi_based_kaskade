Pepe check of weight-aware prolongation splitting
================================================

Prepared 2026-11-17. No cluster job has been submitted by this workflow.
No production solver, connectivity, assembly or compression code is changed here.

Why EMI needs its own check
--------------------------
Kaskade commit d8c532f changes the split-ON constant selected by the shared
SpaceTransferDataCompression::preProlongate method used by both Poisson and EMI.
It does NOT change or normalize the averaging weights themselves.
Earlier Pepe results predate that correction and do not establish current behavior.

For a producer s and shared DOF i, let w_si be its existing diagonal averaging
weight. Splitting represents w_si*u_si as w_si*(u_si-c_s) + w_si*c_s.
After inverse lifting, each receiver restores the producer's w_si*c_s using
the original producer weights exchanged during setup. The accumulated weighted
sum is divided by the original sum of weights. This identity allows different
constants and weights for every producer, including DOFs with several I/E
neighbors. It does not require one myocyte per subdomain or a bipartite graph.
Exact corners remain excluded. Interface graphs and producer plans are unchanged.

The new selector tests zero, the arithmetic mean and the rho-squared weighted
mean, after rounding each candidate to its transmitted float32 constant.
It evaluates the existing graph transform of the weighted remainder and picks
the narrowest coefficient range (smaller magnitude breaks range ties).
The zero candidate prevents a range larger than the unsplit remainder on the
SAME input. This is not a guarantee of smaller decoded error or convergence.
Splitting ON can change quantized values, solution errors and iteration counts.
Reversible recentering is separately checked to preserve quantized reconstruction.

Controlled validation
---------------------
Pepe; P2 baseline; 2 physical nodes; 4 MPI ranks; 3 timesteps; unchanged dt=0.01,
FCG and true relative tolerance 1e-7. Compressed runs use 12- or 16-bit transform-first
Method B, one lifting level, Huffman/bit-length coding and compact packing,
with exact corners. Four trajectories, with matched current raw reference:

  raw-fcg                  No compression
  B-recenter-split-off      Full B+R; splitting OFF
  B-split-on                Full B; corrected splitting ON; recentering OFF
  B-recenter-split-on       Full B+R; corrected splitting ON; recentering ON

Each solve must pass actual true residual checks, completion records on every
rank, topology/ownership consistency, initial-state and operator fingerprints,
and serialized-byte accounting. State errors are measured for every timestep
against matched converged raw FCG. Recoding OFF/ON must give identical state
snapshots and residual histories for the split-ON pair.

Native tests additionally exercise weighted multi-producer assembly, different
local numbering, exact corners, invalid headers and split/recenter packet parity.
They are synthetic checks, not independent Pepe PDE solves.

Usage on cluster, from emiModelBddcSdc_new_scratch
------------------------------------------------
  out=/scratch/htc/$USER/emi-weight-aware
  mkdir -p "$out"
  sbatch --output="$out/check-%j.log" --error="$out/check-%j.log" \
    tests/emi_weight_aware_20261117/submit_check.sbatch

The default remains 16 bits. To repeat the SAME four cases at 12 bits:
  sbatch --job-name=pepe-split-b12 \
    --output="$out/check-%j.log" --error="$out/check-%j.log" \
    tests/emi_weight_aware_20261117/submit_check.sbatch 12

The selected precision applies to restriction and prolongation and is recorded
in the manifest, report and summary. Native recentering checks use that precision.
Each job has a new result directory, preserving job 6249125 and all old measurements.
No production numerical implementation changes are needed to select 12 bits.

The script ALWAYS builds a fresh isolated executable against the corrected
shared headers. It does not reuse the old EMI binaries, stale Poisson builds,
or overwrite results. Source, library and executable hashes are recorded/checked.
The script requests a 30-minute limit; expected work is one build and four small
three-timestep runs, usually several minutes after allocation (queue wait unknown).

After completion, substitute the actual submitted ID:
  job=YOUR_JOB_ID
  sacct -X -j "$job" --format=JobID,State,ExitCode,Elapsed
  r="$out/check-$job/results"
  cat "$r/report.txt"
  cat "$r/summary.csv"
  cat "$r/per-step.csv"

Per-run files preserve residual histories and state snapshots, R/P phase byte
counts, topology/ownership data and serialized-message byte entropy for all
inter-rank/inter-node/same-node scopes. Entropy is pooled bits/byte, not Huffman
training entropy. Observer timings are diagnostic, not repeated speed benchmarks.
No solution-error threshold is invented; measured errors require scientific review.
PASS refers to convergence and structural checks, not a universal accuracy claim.

Offline plan (no simulation, no allocation):
  python3 tests/emi_weight_aware_20261117/run_check.py
  python3 tests/emi_weight_aware_20261117/run_check.py --bits 12

Offline command/invariant tests:
  python3 -m unittest discover -s tests/emi_weight_aware_20261117 -p test_check.py
