# Pepe: short time-integration checks on 2, 4 and 8 nodes

The completed two- and four-node timing results are preserved. A source/result
checkpoint was written before this stage to:

```text
output/checkpoint-before-time-accuracy-20261107T215010Z/
```

It contains both repository revisions/branches, tracked patches, relevant
source archives and checksums, and copies of four validation/timing manifests
and summaries. Full original result directories remain untouched. The image,
mesh data and executable are not duplicated in this source checkpoint. This
is a filesystem checkpoint, not a new Git commit or push.

## Experiment

Each allocation uses two MPI ranks per node:

| Nodes | Ranks | Subdomains |
|---:|---:|---:|
| 2 | 4 | 20 |
| 4 | 8 | 20 |
| 8 | 16 | 20 |

The fixed Pepe mesh is unchanged. With 16 ranks and 20 subdomains, work is not
necessarily balanced evenly, and more nodes need not improve runtime.

The batch script runs these stages sequentially and stops on failure:

1. Fresh full one-step routing/local-bypass validation on the allocated layout.
   This includes direct, CG, FCG and FGMRES cases. It is required for eight
   nodes, which has not yet been validated, and gives an exact first-step
   reference on the actual hosts for every layout.
2. Four ten-step trajectories: direct reference (one rank), raw FGMRES BDDC,
   FGMRES with all-interface compression+split, and FGMRES with local bypass+split.
3. Optional repeated one-step timings, enabled with `--with-timings`.
   Use this for eight nodes to match the completed two-/four-node campaigns.

Every trajectory starts from the same default excited-domain initial state.
Compression begins at physical step 1. The configuration is semi-implicit
Euler, dt=0.01, 16-bit transform-first coding, one restriction pass, FGMRES
restart 30, maximum 1000 iterations per step and true relative tolerance 1e-7.
Histograms and VTK are off; full-precision potential snapshots are on.
Transport is self,vader,tcp with Vader single-copy none.

CG/FCG constant-split experiments remain one-step-only. The existing FGMRES
limit of ten steps remains in force. No solver algorithm, build or numerical
guard is changed by these new scripts.

The first ten-step snapshot is required to match the corresponding one-step
reference exactly, including the exported RHS/correction/residual. If this
fails, inspect the difference before proceeding; do not relax the check without
understanding it. Subsequent RHS values evolve from each run's own state and
are allowed to differ. Steps 1, 5 and 10 are checkpoints of one trajectory,
not three independently restarted simulations.

## Submission

Use a fresh login shell, outside an old interactive allocation, and enter the
EMI repository root. The scripts use the existing image and Kaskade paths;
export `emi_sif` and `emi_kaskade` first if using different locations. Do not
rebuild or change the executable while these jobs are pending or running.

The following submits the layouts in sequence. A failed job prevents the next
job from starting. Pending dependent jobs may then need cancellation/re-submission
after the failure is resolved.

```bash
mkdir -p output
emi_job2=$(sbatch --parsable \
  --nodes=2 --ntasks=4 --ntasks-per-node=2 \
  tests/submit_emi_time_accuracy.sbatch)
echo "Two-node job: $emi_job2"

emi_job4=$(sbatch --parsable --dependency="afterok:$emi_job2" \
  --nodes=4 --ntasks=8 --ntasks-per-node=2 \
  tests/submit_emi_time_accuracy.sbatch)
echo "Four-node job: $emi_job4"

emi_job8=$(sbatch --parsable --dependency="afterok:$emi_job4" \
  --nodes=8 --ntasks=16 --ntasks-per-node=2 \
  tests/submit_emi_time_accuracy.sbatch --with-timings)
echo "Eight-node job: $emi_job8"
```

Only continue a submission if the previous `sbatch` succeeded and returned a
job number. Each job requests four hours, 24 GB per node and one CPU per rank.
Each ten-step solve has a 40-minute Slurm step limit. The batch stage does not
launch mpirun inside srun; it uses the existing host Slurm/Apptainer launcher.

```bash
squeue -u "$USER"
tail -f "output/emi-time-accuracy-${emi_job2}.log"
```

The three jobs have separate logs and output directories:

```text
output/emi-time-campaign-nN-JOB_ID/
  one-step-reference/
  accuracy-10steps/
    environment.json
    summary.json
    accuracy.md
    direct/
    fgmres-uncompressed-destination/
    fgmres-compressed-split-destination/
    fgmres-compressed-split-local-bypass/
  timings/                         only with --with-timings
```

For an already validated layout in an active allocation, the standalone
accuracy runner also accepts an existing matching reference:

```bash
python3 tests/run_emi_time_accuracy_check.py \
  --reference output/emi-mpi-routing-pepe-20261107T182241Z-76236 \
  --slurm-container "$emi_sif" --kaskade "$emi_kaskade" --steps 10
```

That example requires four nodes/eight ranks. Use `--plan-only` to inspect
commands without launching applications. Do not run it alongside the batch
workflow on the same allocated CPUs.

## Interpretation

`summary.json` records errors at every physical step, including raw and centered
relative Euclidean potential errors, mean offsets and maximum absolute errors
against both the direct trajectory and the uncompressed BDDC trajectory. It
also records iterations and true residuals per step and cumulative local/inter-node
payloads. `accuracy.md` highlights steps 1, 5 and 10. These timings include
validation snapshots and must not be mixed with the snapshot-free timing campaigns.

PASS means the configured residual, completion, snapshot and accounting checks
passed. No scientific solution-error threshold has been agreed. Small centered
error is not permission to discard a raw offset. These are global algebraic
potential errors, not membrane-voltage errors or H^(1/2)/H^(-1/2) norms; those
measurements remain separate work.

Record the actual nodes for every job. Different physical machines and different
local/remote compression eligibility confound a pure strong-scaling claim.
Compare each compressed trajectory with its matched raw/direct reference first.
The new scripts have local command/parser tests; the cluster executions are
still required to validate this multi-step, multi-node workflow.
