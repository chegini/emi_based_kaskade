# EMI model with Kaskade7

## Controlled one-step grouped prolongation experiment

`tests/run_emi_grouped_prolongation_check.py` runs **10cells, one time step**
with fresh direct, CG16/16 and FCG16/16 controls, followed by FCG8/8 with
`--bddcProlongateGroupSize 32`. It preserves the saved physical problem,
eight BDDC ranks, refinement level 2, one restriction pass and true tolerance
`1e-7`. Restriction uses the existing exact-zero quantizer. Prolongation uses
the original affine quantizer on consecutive groups of 32 non-corner values in
each producer's union interface order. A group has the same grid for its local
copy and every destination; exact corners retain their existing representation.

Scale and minimum remain **float32**. Packets carry a uint32 run count and a
uint32 scalar count plus float32 scale/minimum for each contiguous group run.
All descriptors, padding and exact-corner trailers are counted. Integer lifting
and Huffman coding remain enabled. Grouping is optional (default zero), and
does not change the CG/FCG recurrence, number of restriction passes or MPI
exchange schedule. Metadata and packing work increase; smaller groups can
therefore cost more despite using fewer bits per quantized value.

After syncing both repositories and rebuilding `emiModel`, run in the container:

```bash
python3 tests/run_emi_grouped_prolongation_check.py \
  --reference-campaign output/emi-sd-cg8-20261101T121744Z-9497
```

For the controlled follow-up with **16-value groups**, add `--group-size 16`.
The default remains 32. This changes only the candidate's prolongation group
size; precision remains 8/8 and restriction remains one pass. No rebuilding is
needed when changing this runtime setting. The fresh controls and actual-byte
accounting are identical, and the chosen size is recorded in the output.

The source campaign is read only. Each invocation creates a new output folder
and saves exact commands, executable digest, logs, metrics and full-precision
raw/centered errors against the fresh direct solve. The fresh CG16 first-step
residual and iteration count must match the saved control. A failed control
stops the experiment. Candidate convergence and total-byte savings against
both CG16 and FCG16 are reported separately. No hidden solution-error threshold
is applied. Byte totals include intra-rank deliveries but exclude MPI framing
and codebook training; validation timings are not benchmark timings.

Positive `r.z` in the earlier frozen-buffer screen is **not** a convergence
guarantee. This runner does not launch a ten-step campaign, add refinement
passes, raise precision, relax tolerance or apply post-solve constant correction.

This project is the clean EMI implementation based on `kaskade7_new`. It
supports EMI time integration, SDC, algebraic adaptivity (AA), BDDC, and
optional compressed BDDC interface transfer.

The old project and old compression tree are separate:

```text
emiModelBddcSdc  -> kaskade7_new
emiModelBddcSdc_new           -> kaskade7_compression
```

The old `kaskade7_compression` tree must not be modified by this project.

## Podman environment

For a first-time setup, including creation of an MPI-enabled image, see
[`docs/mpi-container-setup.md`](docs/mpi-container-setup.md).

Set these two host paths first:

```bash
PATH_TO_KASKADE=/path/to/kaskade7_new
PATH_TO_PARENT_PROJECT=/path/to/BDDC
PATH_TO_BUILD_ENVIRONMENT=/path/to/kaskade7-build-environment.tar.gz
```

Here, `PATH_TO_PARENT_PROJECT` is the parent directory that contains
`emiModelBddcSdc` and `input_emi_mesh`.

Run this command from any directory:

```bash
podman run -it --rm \
  -v "${PATH_TO_KASKADE}:/kaskade7" \
  -v "${PATH_TO_PARENT_PROJECT}:/BDDC" \
  -w /BDDC/emiModelBddcSdc \
  localhost/kaskade7-build-environment:latest
```

If the image has not been loaded:

```bash
gunzip -c "${PATH_TO_BUILD_ENVIRONMENT}" | \
  podman load
```

The complete `BDDC` directory is mounted as `/BDDC`. Therefore the mesh files
in `../input_emi_mesh` are visible from the working directory.

## Build

Inside the container:

```bash
make
```

To clean the executable, dependency files, and generated output:

```bash
make clean
```

### Cluster workflow

The Slurm build script compiles the MPI-enabled executable inside the
container:

```bash
sbatch build_emi.sbatch
```

The build is needed after changing the EMI source, changing Kaskade headers,
or when `emiModel` is missing. Reconnecting to the cluster or running on a
different login or compute node does not by itself require a rebuild because
the project is stored on the shared filesystem.

After a successful build, submit the two benchmark runs independently:

```bash
sbatch run_emi_serial_8threads.sbatch
sbatch run_emi_mpi_8ranks.sbatch
```

The first job runs one MPI process with eight assembly/BDDC threads. The
second runs eight MPI ranks with one thread per rank. Logs are written to
`logs/emi-serial8-<jobid>.out` and `logs/emi-mpi8-<jobid>.out`.

The scenario helper is available through `steps_to_run.mk`:

```bash
make -f steps_to_run.mk help
make -f steps_to_run.mk build
```

MPI is optional and is only enabled for the smoke test at this stage. The
default build does not include `mpi.h` or link against MPI:

```bash
make -f steps_to_run.mk build
```

To build the MPI-enabled executable and verify MPI startup plus an
`MPI_Allreduce`, run inside a container with `mpicxx` and `mpirun` available:

```bash
make -f steps_to_run.mk build-mpi MPICXX=mpicxx
make -f steps_to_run.mk run-mpi-smoke MPI_RANKS=2 MPICXX=mpicxx
```

The smoke test runs the existing EMI model independently on every rank. It
does not distribute BDDC subdomains or exchange BDDC data yet. A successful
run prints one line per rank, for example:

```text
MPI smoke test: rank 0/2, allreduce sum=3
MPI smoke test: rank 1/2, allreduce sum=3
```

After an MPI build, return to the normal non-MPI executable with:

```bash
make clean
make -f steps_to_run.mk build
```

## Input files

The commands use the 10-cell test mesh:

```text
../input_emi_mesh/10Cells3d_10extra_mesh.vtu
../input_emi_mesh/10Cells3d_10extra_list_extracellular.txt
../input_emi_mesh/10Cells3d_10extra_list_intracellular.txt
../input_emi_mesh/10Cells3d_10extra_early_excited.txt
```

The input directory is outside the scratch project so it can be shared by all
solver configurations.

The VTU file must contain cell data named `domain`, whose integer-valued tags
identify the material region. Each text tag file starts with the number of
tags, followed by that many integer tags. Keep extracellular and intracellular
tag lists disjoint; excited tags should identify the intracellular regions to
stimulate. The initial potential is `0.5` on excited tags and `0` elsewhere.

The mapper uses one shared potential component for extracellular tags and a
separate disconnected component for each intracellular tag. Thus material
labels affect both the physical interface terms and the global DOF numbering.

## Solver configurations

Each scenario has a target in `steps_to_run.mk`:

| Scenario | Target | Main options |
|---|---|---|
| EMI only | `run-emi` | semi-implicit Euler |
| EMI + SDC | `run-sdc` | `--sdc 1` |
| EMI + SDC + AA | `run-sdc-aa` | `--sdc 1 --algebraicAdaptivity 1` |
| EMI + BDDC | `run-bddc` | `--bddc 1 --bddcCompression 0` |
| EMI + BDDC + compression | `run-bddc-compression` | `--bddc 1 --bddcCompression 1` |
| EMI + SDC + BDDC | `run-sdc-bddc` | `--sdc 1 --bddc 1 --bddcCompression 0` |

MPI-BDDC uses owner-only subdomain storage. With `--mpi 1 --bddc 1`, each rank
constructs and solves only its deterministic round-robin subdomains, exchanges
remote interface payloads with one collective for restriction and one for
prolongation, and contributes coarse and global corrections through MPI
reductions. The old replicated MPI path is removed. Owner-only storage is
automatic whenever MPI-BDDC is enabled.
| EMI + SDC + BDDC + AA | `run-sdc-bddc-aa` | adds AA to SDC+BDDC |
| EMI + SDC + BDDC + compression | `run-sdc-bddc-compression` | compressed BDDC transfer |
| EMI + SDC + BDDC + AA + compression | `run-sdc-bddc-aa-compression` | compressed transfer plus AA |

For example:

```bash
make -f steps_to_run.mk run-sdc-bddc-aa-compression
```

## Command-line parameters

The executable accepts all runtime options in the form
`./emiModel --option value`. Use `./emiModel --help` for the generated option
list. The defaults below are defined in `EMI_model.cpp`.

### Input, output, and discretization

| Option | Default | Meaning |
| --- | ---: | --- |
| `--input` | built-in 2-cell VTU | VTU mesh with domain cell data |
| `--extra_set` | built-in 2-cell file | Extracellular material tags |
| `--intra_set` | built-in 2-cell file | Intracellular material tags |
| `--excited` | built-in 2-cell file | Initially excited material tags |
| `--dir` | `output` | Output directory |
| `--refine` | `0` | Number of uniform mesh refinements |
| `--order` | `1` | Finite-element polynomial order |
| `--nThreads` | `2` | Assembly/setup task limit and Kaskade global worker limit per rank; `1` means sequential Kaskade tasks, `0` means automatic |
| `--vtk` | `1` | Write VTK output: `0` no, `1` yes |
| `--diagnostics` | `0` | Print per-step vector diagnostics |
| `--profile` | `0` | Print solver timings plus an end-to-end phase table with MPI min/mean/max |

### Time integration

| Option | Default | Meaning |
| --- | ---: | --- |
| `--finalTime` | `0.01` | Final simulation time |
| `--dt` | `0.01` | Time-step size |
| `--maximumNumberOfTimeSteps` | `0` | Maximum number of steps; `0` uses `finalTime/dt` |

### Optional compact output

Add `--quiet 1` to a non-SDC BDDC or direct run for compact console output.
The default, `--quiet 0`, retains the existing output. Quiet mode hides routine
assembly messages, timing/compression tables and MUMPS informational output.
The true-residual iteration table is enabled by default for supported BDDC CG
runs (see the convergence-report section below). When that report is disabled,
`--bddcVerbose 1` enables the internal iteration table even in quiet mode;
`--bddcVerbose 0` shows only per-step summaries. The Pepe and Rizzo runners set
`--bddcVerbose 1` so the internal table remains available as an alternative.
During BDDC setup, rank zero lists every local submatrix
with its material tag and dimensions once material assembly/exchange and local
matrix extraction have completed, including submatrices assembled by other MPI
ranks. The coarse-interface counts (corners, edges and faces) are printed once
on rank zero during setup, including in quiet mode. The true-residual stopping
rule is printed once on rank zero before the
first solve, including in quiet mode. When the internal table is selected,
a legend follows the stopping rule: `||du||_A` is the energy norm of the unscaled search
direction, and `p[0]` is the internal residual/preconditioned-residual inner
product. The table values precede the update; the stopping rule uses the
accumulated correction and the separately measured global residual. The RHS
norm is fixed during each timestep solve. Quiet mode retains mesh information,
per-step progress, failures,
explicitly requested diagnostics/audits, and final wall time. Solver settings,
stopping checks, MPI reductions and VTK output are unchanged.

With `--bddcConditionEstimate 1`, every BDDC step includes `condition_estimate`:

```text
BDDC step 1/10: iterations=44 rel=6.66935e-08 status=converged condition_estimate=39.6139
```

This is the RHS-dependent CG/Lanczos estimate for the preconditioned operator
after interior elimination, not the exact condition number of the full matrix.
An unavailable estimate is reported explicitly. Without true-residual stopping
or diagnostics, quiet mode labels the existing internal norm as `internal_norm`.
Direct solvers have no iteration count; MUMPS SPD prints the checked relative
residual for every step. Quiet mode does not add a condition estimator to direct
solves or enable BDDC condition estimation automatically.

After recompiling in the container, the existing Rizzo runners accept `QUIET=1`:

```bash
# Optimized BDDC, ten steps, eight ranks, with per-step condition estimates.
QUIET=1 RANKS=8 VTK=1 bash tests/run_emi_rizzo.sh 10

# MUMPS SPD direct reference, ten steps, one rank.
QUIET=1 bash tests/run_emi_rizzo_direct.sh
```

For matching one-step BDDC and direct solves on the Pepe mesh:

```bash
# Optimized BDDC on eight ranks, with VTK output for the comparison.
QUIET=1 RANKS=8 VTK=1 bash tests/run_emi_pepe.sh 1

# MUMPS SPD direct reference on one rank, followed by the error comparison.
QUIET=1 bash tests/run_emi_pepe_direct.sh 1
```

Both Pepe runners accept a positive step count as their argument; use `10` in
both commands for ten steps. Without an argument, the BDDC runner defaults to
one step and the direct runner retains its previous ten-step default. The
direct runner always uses one MPI rank and computes `finalTime=0.01*steps`.
It saves the MUMPS SPD result under `output/pepe-direct-mumps-spd-<STEPS>steps`
and compares it with the matching optimized BDDC solution
`output/pepe-bddc-optimized-cond-np8-<STEPS>steps/emiBDDCLast.vtu` when present.
Set `BDDC_REFERENCE=/path/to/emiBDDCLast.vtu` to compare against a different
matching Pepe run. A missing reference is reported explicitly; the comparison
report is saved as `compare-bddc.log` in the direct output directory.

Omit `QUIET=1` for the original output. The direct runner still prints its final
VTU comparison when a BDDC reference file exists; its verbose memory/timing
report is disabled in quiet mode. Both runners keep their existing log paths.

### Optional true-residual convergence report

The default report shows the accepted CG update count, step length,
`||b-A*du||_2`, `||b||_2`, and their ratio for non-SDC BDDC CG runs with
`--bddcRelativeTolerance` greater than zero. The application setting defaults
to `--bddcConvergenceReport -1` (automatic): it enables the report for those
runs and keeps it off for other solver modes. This changes reporting only;
it does not enable CG or change any stopping tolerance.

`--bddcConvergenceReport 1` explicitly enables the report and requires the
supported settings above. `--bddcConvergenceReport 0` disables it and restores
the original internal table when `--bddcVerbose 1` is selected. The final
residual and requested condition estimate remain visible in both modes.

Row zero measures the initial accepted correction after interior elimination
and has no step length. Subsequent rows measure the accumulated correction
AFTER applying the listed alpha. Rejected corrections produce no additional
row. The RHS norm is fixed throughout a timestep's linear solve; for zero RHS,
the `relative_residual` field uses the documented absolute-residual fallback.
The report observes existing stopping checks: it adds no residual evaluations
and changes no matrix, correction or convergence decision. It adds optional
console/file I/O; file-open/write failure is propagated to all MPI ranks.

Rank zero writes `bddc-convergence-step<N>.csv` in the run's output directory,
with 17-digit values. A failing solve retains its available history. For Pepe:

```bash
QUIET=1 RANKS=8 VTK=1 bash tests/run_emi_pepe.sh 1
```

Both Pepe and Rizzo runners enable the report by default. Set
`CONVERGENCE_REPORT=0` to select their original internal iteration table.
For plotting, run the standalone helper on a machine with Matplotlib installed:

```bash
python3 plot_emi_bddc_convergence.py \
  output/pepe-bddc-optimized-cond-np8-1steps/bddc-convergence-step1.csv
```

The saved PNG has two panels: true relative residual on a logarithmic axis
with the stopping target, and CG step length. Use `--absolute` to plot the
absolute residual, fixed RHS norm and corresponding absolute stopping target.
`--format pdf` or `--format svg` provides vector output. Exact zero residuals
are identified in a note rather than assigned artificial positive values on
the logarithmic axis. Multiple CSV paths can be supplied to plot several steps.
CSV generation has no Python dependency; the figure helper needs Matplotlib
(`python3 -m pip install matplotlib` in your chosen plotting environment).
The previously exported container image must be rebuilt to include these files,
or compile the updated host checkout mounted into the development container.

### PCG and direct linear solvers

The default EMI linear solver is PCG. Set `--direct 1` to use the direct
solver instead.

| Option | Default | Meaning |
| --- | ---: | --- |
| `--direct` | `0` | `0` PCG, `1` direct solver |
| `--directBackend` | `umfpack` | Global reference factorization: `umfpack`, or opt-in `mumps-spd` for a single-rank, non-SDC direct solve |
| `--cgTol` | `1e-8` | PCG stopping tolerance |
| `--maxCGIter` | `10000` | Maximum PCG iterations |
| `--CG_shift` | `1` | Apply constant-shift correction after PCG |

Example:

```bash
./emiModel --direct 0 --cgTol 1e-8 --maxCGIter 10000
```

If UMFPACK's global LU factorization runs out of memory, the optional
`--direct 1 --directBackend mumps-spd --bddc 0 --sdc 0` path uses the already
linked sequential MUMPS solver in positive-definite mode. It checks the full
assembled scalar matrix for finite entries, positive diagonal entries and
`max|A-A^T| <= 1e-12*max|A|` before submitting only the nonzero lower triangle.
This input check is not an SPD proof; a failed MUMPS SPD factorization stops
the run. Only exact zeros and the redundant upper triangle are omitted.
No diagonal shift, coefficient change or dropping tolerance is applied.
Changing the backend can change floating-point roundoff, and reduced memory
use is a possibility rather than a guarantee. UMFPACK remains the default;
BDDC's local/coarse solvers are unaffected.

MUMPS factors the fixed global matrix once, reuses it for all steps, and prints
the absolute and relative residual of each correction against the original
full matrix before advancing the state. Nonfinite residuals stop the run.
There is no iterative stopping tolerance in this direct path; check that its
reported residuals are well below the BDDC target before treating the output
as a reference. The path requires `-np 1`; increasing MPI ranks does not
distribute this sequential factorization.

After rebuilding the application in the MPI-enabled container:

```bash
bash tests/run_emi_rizzo_direct.sh
```

This runs the same ten-step Rizzo computational problem with VTK enabled,
records output in `output/rizzo-direct-mumps-spd-10steps/run.log`, and compares
the final solution with the existing eight-rank BDDC VTU when present. If
`/usr/bin/time` is available, it measures the application inside the MPI
launch, including its maximum resident set size. The direct-matrix input
checks can be tested independently of Dune/MUMPS:

```bash
c++ -std=c++20 -O2 -Wall -Wextra -pedantic \
  tests/test_direct_matrix.cpp -o /tmp/test_direct_matrix
/tmp/test_direct_matrix
```

### Spectral deferred correction

Enable SDC with `--sdc 1`.

| Option | Default | Meaning |
| --- | ---: | --- |
| `--sdc` | `0` | Enable SDC: `0` no, `1` yes |
| `--sdcCollocationPoints` | `3` | Target collocation-point count |
| `--sdcStartCollocationPoints` | `3` | Initial collocation-point count |
| `--minimumSdcSweeps` | `3` | Minimum SDC sweeps |
| `--maximumSdcSweeps` | `5` | Maximum SDC sweeps |
| `--sdcSweepType` | `1` | `0` Euler sweep, `1` LU sweep |
| `--sdcTolerance` | `1e-6` | Correction-norm stopping tolerance |
| `--sdcAbsoluteTolerance` | `1e-12` | Absolute SDC error-estimate tolerance |
| `--sdcInitialContraction` | `0.2` | Initial contraction estimate |

### Algebraic adaptivity

In this implementation, `AA` means algebraic adaptivity. It selects active
degrees of freedom or complete BDDC owner subdomains for later SDC sweeps.

| Option | Default | Meaning |
| --- | ---: | --- |
| `--algebraicAdaptivity` | `0` | Enable AA: `0` no, `1` yes |
| `--algebraicAdaptivityTolerance` | `0.0` | AA selection tolerance; `0` disables selection |

Example:

```bash
./emiModel --sdc 1 --algebraicAdaptivity 1 \
  --algebraicAdaptivityTolerance 1e-4
```

### BDDC

| Option | Default | Meaning |
| --- | ---: | --- |
| `--bddc` | `0` | Enable BDDC |
| `--bddcIterations` | `3000` | Maximum BDDC iterations per time step |
| `--bddcTolerance` | `1e-8` | Legacy internal search-direction energy-norm tolerance |
| `--bddcRelativeTolerance` | `0` | True global relative residual target for non-SDC BDDC; `0` keeps the legacy stopping rule |
| `--bddcInterfaceTypes` | `7` | Interface flags: `1` corner, `2` edge, `4` face, `7` all |
| `--bddcVerbose` | `0` | Show the internal iteration table when the convergence report is off, including in quiet mode; enabled by the Pepe/Rizzo runners |
| `--bddcUseCg` | `1` | Use CG in the BDDC coarse-solve path |
| `--bddcConvergenceReport` | `-1` | Automatic true-residual table/CSV for non-SDC CG with relative-residual stopping; `0` disables, `1` explicitly enables |
| `--bddcConditionEstimate` | `0` | Print an RHS-dependent CG/Lanczos Ritz estimate of the BDDC-preconditioned spectral condition number; requires uncompressed, non-SDC CG |
| `--bddcAudit` | `0` | Audit the constant mode and coarse interface groups; requires `--bddc 1 --sdc 0` |
| `--bddcMaterialAssembly` | `0` | Non-SDC material matrix assembly: `0` replicated legacy path; `1` assemble on each subdomain's MPI owner, then exchange sparse nonzeros |
| `--bddcMaterialAssemblyCheck` | `0` | Audit owner-assembled material matrices against the replicated assembler before interface redistribution; requires material mode `1` |
| `--bddcInterfaceAssembly` | `0` | Non-SDC interface mass assembly: `0` legacy pair-by-pair, `1` shared traversal on every rank, `2` shared traversal on source MPI owners |
| `--bddcInterfaceAssemblyCheck` | `0` | Compare optimized directed mass blocks against legacy pair assembly before redistribution; requires assembly mode `1` or `2` |
| `--bddcSplitInterfaces` | `0` | Separate coarse averages by connected electrical component; requires `--bddc 1 --sdc 0` |
| `--bddcOperatorAudit` | `0` | Diagnose uncompressed BDDC operators, then exit without timestep solves; requires `--bddc 1 --sdc 0 --bddcCompression 0` |

With `--bddcConditionEstimate 1`, each non-SDC timestep also prints
`lambda_min`, `lambda_max`, `condition_estimate` and `krylov_dim`. The existing accepted
CG alpha/beta coefficients reconstruct the small Lanczos tridiagonal, whose
extreme eigenvalues are computed by scaled Sturm bisection. This adds no
operator applications or MPI reductions and leaves the solver/stopping rule
unchanged. The Kaskade headers must include `BDDCSolver::lastCgStep()`.
The estimate assumes a fixed SPD preconditioner and describes the spectrum
of the symmetrically preconditioned operator on the sampled Krylov space
after interior elimination. It is not the exact condition number of the full
matrix `A`, nor the Euclidean singular-value condition number of `B*A`.
It depends on the RHS, can miss eigenmodes, and is affected by floating-point
loss of CG orthogonality. A fresh history is used per RHS; rejected updates
are excluded. Fewer than two accepted iterations or invalid spectral data
produce `condition_estimate=unavailable` with an explanatory status.
The option requires `--bddc 1 --sdc 0 --bddcUseCg 1 --bddcCompression 0
--bddcOperatorAudit 0`. See the
[Netlib Lanczos connection](https://www.netlib.org/linalg/html_templates/node118.html).

The standalone estimator check needs only a C++20 compiler:

```bash
c++ -std=c++20 -O2 -Wall -Wextra -pedantic \
  tests/test_bddc_condition.cpp -o /tmp/test_bddc_condition
/tmp/test_bddc_condition
```

The existing Rizzo inputs are `../input_emi_mesh/rizzo_combi_domi.vtu` and
the accompanying `_extracellular.txt`, `_intracellular.txt`, `_excited.txt`
files: 20 extracellular tags, 20 intracellular tags and excitation tag 147.
After recompiling the application in the MPI-enabled container:

```bash
# One-step smoke runs at 1 and 8 ranks; optimized assembly and Ritz estimates.
bash tests/run_emi_rizzo.sh 1

# Ten steps with VTK output, on eight ranks only.
RANKS=8 VTK=1 bash tests/run_emi_rizzo.sh 10
```

The script retains the Pepe computational benchmark coefficients
`sigma_i=3`, `sigma_e=20`, `C_m=1`, `R=0.1`, `dt=0.01`, the weak Robin
boundary and the existing initial-excitation mechanism. It does not claim
physical equivalence to an openCARP setup or introduce a units correction.
Default VTK output is off for timing. Logs are written to
`output/rizzo-bddc-optimized-cond-np<RANKS>-<STEPS>steps/run.log`.

With `--bddcRelativeTolerance 1e-7`, the non-SDC BDDC solve stops when the
assembled, unpreconditioned correction residual satisfies
`||b-A*du||_2 / ||b||_2 <= 1e-7`. The denominator is the original global RHS norm,
corresponding to a zero initial global correction, even though the solver first
eliminates interior residuals locally. It is not normalized by the remaining
residual after that elimination. For a zero RHS, the criterion uses the absolute
residual instead. The initial accepted correction and every accepted update are
checked; an already solved system needs no CG calls.

This option overrides `bddcTolerance`. It does not change the matrix, physical
BCs, scaling or coarse constraints. Each check assembles the global correction
and performs a matrix-vector product on rank zero. The two norms are broadcast
so all MPI ranks make the same stopping decision; replicated norms are not
summed. Profiling reports checking and correction-assembly costs separately.
These extra checking collectives are outside the BDDC transfer compression
report; total wall time includes them.
The global residual and convergence status are printed even with diagnostics
disabled. If the iteration limit, stagnation or breakdown prevents meeting the
target, the run fails before advancing that timestep or writing a final
solution. The option requires `--bddc 1 --sdc 0 --bddcOperatorAudit 0`.

Using a true residual can change iteration counts and the final solution within
iterative-solver accuracy. The legacy rule remains the default. Comparing to
other solvers still requires matching the mesh, physical BCs, coefficients,
discretization, initial-guess policy and norm reference. In particular, this
application solves for a timestep correction; a solver working on the full
potential can have a different RHS and reference norm even for equivalent
physical equations.

The stopping-policy regression check has no Dune dependencies:

```bash
c++ -std=c++20 -Wall -Wextra -pedantic tests/test_bddc_stopping.cpp -o /tmp/test_bddc_stopping
/tmp/test_bddc_stopping
```

With `--bddcAudit 1`, rank zero reports the current LHS applied to the nodal
constant field, its normalized residual and Rayleigh quotient, and the initial
RHS constant component. The Robin boundary term remains active. A small
`||A*1||_inf / ||A||_inf` identifies a candidate near-null mode; it does not by
itself establish poor BDDC conditioning. The RHS compatibility ratio is
informational for the current Robin problem, not a pure-Neumann solvability test.

The audit writes `bddc_audit/interface_groups.csv` and
`bddc_audit/interface_components.csv` below the output directory. Groups match
the subdomain-ID grouping in `InterfaceAverages`. Cell-connected components join
only DOFs sharing an original volume element within the group, so duplicated
membrane traces remain separate. Electrical label `0` is the shared extracellular
space; other labels are intracellular material tags, matching the mapper.
Multiple components or labels identify candidates for additional constraints,
not proven errors. The hypothetical constraint count does not change the solver.
The `C*1` check verifies agreement of the actual constraint row sums across
subdomains; it does not test near-null-mode reproduction by the coarse basis.

Material assembly mode `1` assigns each material's existing filtered assembler
call to its deterministic BDDC owner, including `--masterWorker`. For 20
subdomains and eight ordinary ranks, each rank assembles two or three matrices.
The quadrature, global numbering, symmetric extraction, physical coefficients,
Robin BC, initial excitation and solver tolerances are unchanged. No new
stimulation or Dirichlet configuration is introduced by this option.

Each owner broadcasts its nonzero entries as `(row, column, double-bits)` words,
without a floating-point reduction. Receivers copy the full global matrix
pattern and fill the received entries; retaining replicated matrices preserves
the current interface redistribution, local extraction, weights, dumps and
audits. One material is exchanged at a time to bound temporary storage, and
large messages are chunked within MPI's integer count limit. Empty messages
produce zero matrices; nonfinite entries, invalid indices and writes outside
the pattern fail before the matrices are used. Collective error checks prevent
peers proceeding to the next broadcast after an owner/receiver failure.

With one rank or MPI disabled, material mode `1` assembles every subdomain
locally and needs no exchange. Mode `0` remains the default and adds no material
exchange collectives. With `--profile 1`, `material subdomain assembly` and
`material matrix exchange` are separate rows. Waiting for other owners is
included in assembly, and allocation/initialization of replicas is included in
exchange; compare total time as well as these phases. Saved assembly work does
not guarantee a speedup once transfer and copying are included.

`--bddcMaterialAssemblyCheck 1` additionally reassembles every material on every
rank and compares each entire sparse operator before interface redistribution,
including entries absent from either pattern. The tolerance is
`maxAbsDifference <= 1e-14 + 1e-11*maxReference`. Audit cost is charged to
`diagnostics and interface splitting`; keep this option off for timings.

After rebuilding **only the application**, run these commands manually in the
MPI-enabled container, in this order:

```bash
# Serial codec tests + real MPI tests at 1, 2, 8 ranks; then one-step
# material matrix audits and true-residual solves at 1 and 8 ranks.
bash tests/run_emi_material_comparison.sh audit

# Ten-step VTK runs at 1/8 ranks, compared against the saved direct and
# replicated-material BDDC results from the previous validated experiment.
bash tests/run_emi_material_comparison.sh accuracy

# Ten-step direct run and both material modes at 1/8 ranks. Every run
# disables VTK, diagnostics and matrix audits; all use one thread per rank.
bash tests/run_emi_material_comparison.sh benchmark
```

The script uses the Pepe mesh, `C_m=1`, `R=0.1`, `dt=0.01`, shared interface
assembly mode `2`, unchanged coarse constraints and true relative tolerance
`1e-7`. It does not depend on a pre-existing shell argument array. Accuracy
mode requires these previously validated files:

- `output/pepe-direct-interface-check-10steps/emiLast.vtu`
- `output/pepe-bddc-interface-check-np1-10steps/emiBDDCLast.vtu`
- `output/pepe-bddc-interface-check-np8-10steps/emiBDDCLast.vtu`

Each case records its full output under its own `output/pepe-material-...`
directory as `run.log`, and accuracy mode records both solution comparisons.
`unit` runs only the standalone tests. The test matrix fixture matches the
scalar-block/proxy and prefix-only iterator API of `NumaBCRSMatrix`; the
production-size audit in the full executable remains necessary.

The Pepe container validation passed the standalone tests at 1, 2 and 8 MPI
ranks. Full material-matrix audits at 1 and 8 ranks reported zero maximum
absolute and relative L2 differences from replicated assembly. After ten
steps, the exported final solutions matched the previous same-rank BDDC
outputs exactly. Relative L2 errors against the direct solution were
`1.53019e-8` (1 rank) and `1.52973e-8` (8 ranks), with maximum absolute error
approximately `1e-7`. Both BDDC modes required 41 iterations per step and met
the unchanged true relative residual tolerance `1e-7` at every step. These
checks establish numerical agreement for this case, rather than increased
accuracy or a discretization-error estimate.

The matched ten-step benchmark disabled VTK, diagnostics and matrix audits:

| Solver / material mode | MPI ranks | Wall time (s) |
| --- | ---: | ---: |
| Direct | 1 | 147.70 |
| BDDC, replicated material assembly (`0`) | 1 | 233.16 |
| BDDC, owner material assembly (`1`) | 1 | 235.69 |
| BDDC, replicated material assembly (`0`) | 8 | 277.40 |
| BDDC, owner material assembly (`1`) | 8 | 268.05 |

At eight ranks, mean material assembly plus exchange decreased from
`15.724453` to `2.902863 + 0.913891 = 3.816754` seconds, about 76% less time
for these phases. Overall wall time decreased by 3.4%; the solver wrapper
still took about 227 seconds. At one rank, owner mode was 1.1% slower in this
measurement. These are single measurements, so repeated runs are needed to
quantify timing variability. Keep owner assembly as an optional MPI setup
optimization: it preserves the tested solutions, but does not reduce iteration
counts, remove replicated matrix storage, or make BDDC faster than the direct
solver on this benchmark. Material mode `0` remains the default.

Interface assembly modes `1` and `2` replace one full assembler call per
material-interface direction with a single grid traversal. Each directed face
uses the existing face integration order, Kaskade basis evaluator (including
mapper transformations), and linearized mass cache. Same-material and
extracellular/extracellular faces contribute no interface mass. Intracellular
interfaces retain their existing capacitance terms. Both face orientations are
integrated independently using their own inside/outside coordinate maps.

The legacy EMI assembler declares the scalar matrix block symmetric. It
retains only global row >= column during scatter and mirrors those entries
when extracting each directed pair matrix. The shared traversal follows this
convention too; it does not replace a directed pair by an unsymmetrized jump
block. Each assembled directed block is redistributed in the existing pair
order, subtracting half from its source material and adding half to its
neighbor. The global matrix, weak outer Robin term, material assembly, coarse
constraints, weights formula, and solver stopping criterion are unchanged.
Face-wise accumulation may produce floating-point roundoff differences.

Mode `2` uses the same deterministic source-subdomain ownership as the BDDC
solver, including `--masterWorker`. Each directed block is integrated once by
its source owner. Sparse entries are exchanged using integer words carrying
exact double bits, without a floating-point MPI reduction. The current
`BddcData` still retains replicated local matrices for diagonal weights, dumps
and audits, so the sparse setup blocks are gathered to all ranks, including
both owners. This removes duplicated integration, while leaving the existing
matrix storage and fine/coarse solver communication unchanged. With one rank
or MPI disabled, mode `2` is the shared serial traversal. The traversal itself
is sequential within each rank; `--nThreads` continues to control the existing
assembler paths. Legacy mode `0` remains the default until the new path has
been validated on the target meshes.

Use `--bddcInterfaceAssemblyCheck 1` for numerical validation, not performance
measurement: it additionally runs the old assembler for every directed pair,
compares all nonzero entries (including missing or extra entries), and fails
collectively before the solve if `maxAbsDifference > 1e-14 + 1e-11*maxReference`.
It reports the maximum absolute difference and aggregate relative Frobenius
error. Reference assembly is charged to the diagnostics phase. Then rerun
with `--bddcInterfaceAssemblyCheck 0` for timings, using identical BCs,
convergence tolerance, audit and VTK settings in the baseline and new runs.
Compare the true global residual and exported solutions as well as timings.

The sparse assembly/redistribution test needs only C++20:

```bash
c++ -std=c++20 -Wall -Wextra -pedantic tests/test_emi_interface_blocks.cpp -o /tmp/test_emi_interface_blocks
/tmp/test_emi_interface_blocks
```

The real MPI transport test uses Kaskade's lightweight ownership header,
without Dune or the Kaskade library:

```bash
mpicxx -std=c++20 -O2 -Wall -Wextra -pedantic -DKASKADE_HAVE_MPI -I/kaskade7 \
  tests/test_emi_interface_exchange.cpp -o /tmp/test_emi_interface_exchange
for ranks in 1 2 8; do
  mpirun --allow-run-as-root -np "$ranks" /tmp/test_emi_interface_exchange || break
done
```

With `--bddcSplitInterfaces 1`, each owner group is partitioned using that same
original-cell graph. Each connected component receives its own coarse average.
The mapper must give each shared DOF exactly one electrical label; ambiguous
labels cause an error. Shared extracellular potentials keep label `0`, and
duplicated membrane traces stay separate. The global matrix, RHS, Robin boundary
term, diagonal scaling and convergence tolerance are unchanged. Only the coarse
preconditioner changes. This requires the Kaskade `InterfaceAverages` overload
accepting partition labels.

Use `--bddcInterfaceTypes 7` for the experiment: each component is classified
individually, so splitting an edge can produce corners. The audit CSV groups still
describe the original owner groups, while its enabled-constraint count and the
local-data dump reflect the actual split constraints. For the audited Pepe case,
expect 493 constraints instead of 202. Compare iterations, true global residual,
solution error and total wall time; fewer iterations can come with more expensive
coarse setup. Compare timings with the same audit and VTK settings.

The component graph has a standalone regression check (no Dune dependencies):

```bash
c++ -std=c++20 -Wall -Wextra -pedantic tests/test_bddc_components.cpp -o /tmp/test_bddc_components
/tmp/test_bddc_components
```

`--bddcOperatorAudit 1` constructs the fixed local and coarse factors once and
checks the actual uncompressed transfer and local-solve kernels. It reports:

- The external diagonal weights' signs, partition of unity and agreement with
  the normalized local diagonals used by the internal transfer.
- Global and local matrix symmetry, measured by maximum entrywise skew relative
  to the maximum absolute matrix entry.
- Transfer adjointness on shared DOFs and preservation of a constant trace.
- Local constrained-basis reconstruction of a constant from `C_i*1`, including
  its constraint defect. Weak Robin terms can prevent exact reconstruction.
- Sampled reciprocity and quadratic forms (positivity checks) of the complete
  linear preconditioner, its additivity and repeatability after reset.
- The response `B*A*1` of a fresh preconditioner application to the nearly
  constant mode. No assertion that this must equal `1` is imposed.

Here `B` includes the initial interior Dirichlet solve and one unscaled BDDC
correction, before any CG step length is applied. The two deterministic probes
have unit Euclidean norm and zero mean to reduce domination by the near-null
constant mode. Their responses' means are also printed. Norm-scaled and
pair-scaled reciprocity defects use different normalizations; near-zero cross
products can make the latter misleading. A few sampled checks cannot prove
positive definiteness or establish spectral bounds.

The audit exits before timestep solves and final solution VTK output. Use
`--vtk 0` and a separate output directory. It can also be used with split coarse
interfaces or MPI ownership. The inexpensive transfer check simulates message
delivery on rank zero using replicated metadata; the complete preconditioner
probes use the actual configured MPI exchange. This is a diagnostic run, not a
timing benchmark. Ordinary solves are unchanged when the option is zero.

### BDDC compression

| Option | Default | Meaning |
| --- | ---: | --- |
| `--bddcCompression` | `0` | Enable quantized BDDC transfer |
| `--bddcCompressionBits` | `16` | Quantization precision in bits |
| `--bddcCompressionMinInterfaceSize` | `0` | Minimum scalar interface size for thresholded fallback; `0` disables it |
| `--bddcCompressionThresholdInterfaceTypes` | `2` | Interface type receiving the threshold: `2=edge`, `4=face`, `0=all` |
| `--bddcGraphLifting` | `1` | Enable graph-lifting transform |
| `--bddcHuffman` | `1` | Enable Huffman coding |
| `--bddcBitlength` | `1` | Enable bit-length/tail encoding |
| `--bddcCompressionReport` | `1` | Print aggregate transfer-volume statistics |

The full compressed transfer pipeline is enabled with:

```bash
--bddcCompression 1 \
--bddcCompressionBits 16 \
--bddcGraphLifting 1 \
--bddcHuffman 1 \
--bddcBitlength 1 \
--bddcCompressionReport 1
```

The compression stages are quantization, optional graph lifting, zigzag and
bit-length coding, Huffman coding, and payload packing. The decoder applies
the inverse operations in reverse order. Test the accuracy/communication
trade-off by changing `--bddcCompressionBits`, for example to `8`, `12`, `16`,
or `24`.

### MPI

| Option | Default | Meaning |
| --- | ---: | --- |
| `--mpi` | `0` | Initialize MPI and run the MPI communication check |

Serial execution uses `--mpi 0` and does not require `mpirun`:

```bash
./emiModel --mpi 0
```

MPI execution requires the MPI-enabled build and an MPI launcher:

```bash
make -f steps_to_run.mk build-mpi MPICXX=mpicxx
mpirun --allow-run-as-root -np 2 ./emiModel --mpi 1
```

The current BDDC path distributes local factorizations and iterative subdomain
work by ownership. Mesh data, global operators and material/interface assembly
are still replicated, as is the coarse factorization. Consequently, faster
iterations do not guarantee a shorter full run. Compare total wall time and
individual phases when reporting MPI performance.

The driver initializes Kaskade's global thread pool from `--nThreads` before
grid/FE allocation. The first initialization fixes its configuration for the
run. With `--nThreads 1`, global and node-local Kaskade tasks execute inline;
the matching library change also avoids creating unused workers in this mode.
Pool initialization remains inside the reported wall-time interval. For
values greater than one, the value limits global workers; NUMA node workers
retain the library's separate sizing policy. BLAS/OpenMP threading is controlled
separately, not by this application option.

For an MPI benchmark with one computational thread per rank, set the BLAS
environment before launching the executable:

```bash
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
mpirun --allow-run-as-root -np 8 \
  -x OPENBLAS_NUM_THREADS -x OMP_NUM_THREADS \
  ./emiModel --mpi 1 --nThreads 1
```

OpenBLAS uses `OPENBLAS_NUM_THREADS` in pthread builds and `OMP_NUM_THREADS`
in OpenMP builds; see its
[runtime variable documentation](https://www.openmathlib.org/OpenBLAS/docs/runtime_variables/).
The startup log prints the requested Kaskade setting, actual sequential mode,
and these environment values. Environment values describe the request rather
than measuring the BLAS runtime's actual thread count.

Use the same environment, mesh, solver settings and VTK settings for both
`-np 1` and `-np 8`. Check solution agreement and the true global residual
target before interpreting timing changes.

The threading change includes a compiled library source file. Rebuild the
library before relinking the application; skipping the archive rebuild would
retain the old worker-creation behavior. A library-only build can omit the
test/tutorial makefile modules and their generated dependency files:

```bash
make -C /kaskade7 TESTMODULES= TUTORIALMODULES= kasklib
```

With `--profile 1`, the final `EMI end-to-end phases` table separates mesh/grid
loading, material tags/refinement, FE spaces/initial state/DOF support, initial
VTK output, global assembly, BDDC partition metadata, material-local assembly,
material matrix exchange, interface-neighbor discovery, interface mass-block assembly, mass redistribution,
local matrix extraction, weights, diagnostics/splitting, the solver wrapper and
timestep/final output.

The intervals are exclusive within each rank. Phase switches do not add MPI
barriers; three vector reductions at the end collect minimum, mean and maximum
wall seconds, and rank zero prints one summary. Waiting is included in the phase
where it occurs. Different phases can have different slowest ranks, so do not
add their maximum column: use `accounted total`, reduced from each rank's sum.
The accounted interval starts just inside the existing total timer, before
thread-pool initialization, and ends before reporting collectives/printing.

The solver-wrapper row includes the setup and solve phases already reported
inside BDDC; those inner timings are a breakdown, not additional time to add.
SDC wrappers also include any intermediate VTK writes performed internally.
Direct semi-implicit Euler timestep output is separated explicitly. With
`--profile 0`, the new profiler performs no clock queries, MPI reductions or
output. This instrumentation changes no numerical operators or stopping rules.

The timing-accounting test is independent of Dune:

```bash
c++ -std=c++20 -Wall -Wextra -pedantic tests/test_emi_phase_profile.cpp -o /tmp/test_emi_phase_profile
/tmp/test_emi_phase_profile
```

It can also check real MPI min/mean/max accounting with deterministic synthetic
rank times:

```bash
mpicxx -std=c++20 -DKASKADE_HAVE_MPI tests/test_emi_phase_profile.cpp -o /tmp/test_emi_phase_profile_mpi
mpirun --allow-run-as-root -np 8 /tmp/test_emi_phase_profile_mpi
```

### Physical parameters

| Option | Default | Meaning |
| --- | ---: | --- |
| `--penalty` | `1e6` | Boundary penalty parameter |
| `--sigma_i` | `3.0` | Intracellular conductivity |
| `--sigma_e` | `20.0` | Extracellular conductivity |
| `--C_m` | `1.0` | Membrane capacitance |
| `--R` | `0.1` | Gap-junction conductance |

### Complete example

```bash
./emiModel \
  --mpi 0 \
  --sdc 1 \
  --bddc 1 \
  --bddcCompression 1 \
  --bddcCompressionBits 16 \
  --bddcGraphLifting 1 \
  --bddcHuffman 1 \
  --bddcBitlength 1 \
  --bddcCompressionReport 1 \
  --algebraicAdaptivity 1 \
  --algebraicAdaptivityTolerance 1e-4 \
  --direct 0 \
  --cgTol 1e-8 \
  --maxCGIter 10000 \
  --minimumSdcSweeps 2 \
  --maximumSdcSweeps 4 \
  --bddcIterations 100 \
  --bddcTolerance 1e-8 \
  --nThreads 2 \
  --finalTime 0.01 \
  --dt 0.01 \
  --maximumNumberOfTimeSteps 1 \
  --input ../input_emi_mesh/10Cells3d_10extra_mesh.vtu \
  --extra_set ../input_emi_mesh/10Cells3d_10extra_list_extracellular.txt \
  --intra_set ../input_emi_mesh/10Cells3d_10extra_list_intracellular.txt \
  --excited ../input_emi_mesh/10Cells3d_10extra_early_excited.txt \
  --vtk 0 \
  --profile 1 \
  --dir output/emi_sdc_bddc_aa_compression

The helper defaults to two threads and one time step for quick checks. `--nThreads`
limits parallel assembly, BDDC subdomain construction and global Kaskade workers.
Override the thread
count with `NTHREADS`, for example `make -f steps_to_run.mk run-sdc-bddc
NTHREADS=4`. The one-step limit is a test setting in
`steps_to_run.mk`; increase `TIME_ARGS` there for longer runs. For direct
executable runs, `--nThreads` has the same two-thread default and
`--maximumNumberOfTimeSteps 0` means use all steps implied by `finalTime` and
`dt`.

To run all nine configurations sequentially:

```bash
make -f steps_to_run.mk run-all
```

For an SDC+BDDC+AA timing run, use:

```bash
make -f steps_to_run.mk run-sdc-bddc-aa \
  PROFILE_ARGS="--profile 1" VTK=0 \
  TIME_ARGS="--finalTime 1 --dt 0.01 --maximumNumberOfTimeSteps 100"
```

The report separates local operator extraction, residual assembly, local
collocation matrix/RHS construction, interface setup, subdomain construction
and factorization, transfer configuration, BDDC solver/coarse setup, RHS setup,
BDDC solve calls and time, SDC updates, and AA selection. `PROFILE_ARGS` and
`VTK` are also available on the other SDC+BDDC targets. Profiling is off by
default and currently instruments the SDC+BDDC path. It also reports the number
of subdomain constructions. Although interval operators repeat between time
steps, caching their factorizations is deferred: BDDC subdomains retain mutable
solution and transfer state, so reuse first needs a well-defined reset operation
in the BDDC layer. Subdomain construction is parallelized across independent
subdomains, with results restored to subdomain-ID order before creating the
BDDC solver.

## Compression controls

Compression is disabled with:

```bash
--bddcCompression 0
```

Compression is enabled with:

```bash
--bddcCompression 1 --bddcCompressionBits 16
```

The compressed transfer uses the updated implementation in `kaskade7_new`:

```text
kaskade7_new/mg/bddcSpaceTransfer.hh
kaskade7_new/mg/bddcSpaceTransfer.hpp
kaskade7_new/mg/prefixcoder.hh
kaskade7_new/mg/prefixcoder.hpp
kaskade7_new/mg/alphabet.hh
```

The compressed path configures quantization, bit-length coding, and Huffman
coding for both restriction and prolongation. The transform is selected with:

```bash
--bddcGraphLifting 0   # no transform
--bddcGraphLifting 1   # connectivity-based graph lifting
```

The compressed stages can be isolated for volume studies:

```bash
--bddcCompression 1 --bddcHuffman 0 --bddcBitlength 0  # quantization only
--bddcCompression 1 --bddcHuffman 1 --bddcBitlength 0  # quantization + Huffman
--bddcCompression 1 --bddcHuffman 1 --bddcBitlength 1  # full coded pipeline
```

`--bddcCompressionReport 1` prints aggregate logical and transmitted byte
counts, bits saved, savings percentage, compression ratio, and the estimated
shared codebook size.
For SDC+BDDC+AA the report covers all collocation intervals, sweeps, and time
steps. The transmitted count includes the packed per-message metadata and
payload; the shared process-local Huffman codebook is not repeated per
message. Use `--bddcCompressionReport 0` to suppress the report.

Graph lifting is initialized from the local BDDC matrix sparsity, so it uses
finite-element connectivity rather than the arbitrary ordering of an
unstructured interface vector. If compression is disabled, the standard
`SpaceTransfer` implementation is used.

The compressed transfer pipeline is:

```text
values -> optional graph lifting -> quantization -> zigzag/bit-length coding
       -> Huffman encoding -> compressed payload
```

Decoding applies the exact inverse order:

```text
compressed payload -> Huffman decoding -> bit-length/zigzag reconstruction
                    -> dequantization -> inverse graph lifting -> values
```

The restriction and prolongation codebooks are built from the first transfer
data encountered by a process and then reused by later BDDC iterations and
time steps. Later transfers rebuild only when the transfer layout changes; the
normal per-iteration work is encoding/decoding the payload with the existing
codebook. The shared codebook is process-local, so this policy does not yet
describe MPI-rank sharing.

## Output files

Each helper target uses a separate output directory:

```text
output/emi/
output/sdc/
output/sdc_aa/
output/bddc/
output/sdc_bddc/
output/sdc_bddc_aa/
output/compression/bddc/
output/compression/sdc_bddc/
output/compression/sdc_bddc_aa/
```

Typical final files are:

```text
emiSDCLast.vtu
emiSDCBDDCLast.vtu
emiSDCBDDCLastCompression.vtu
emiBDDCLast.vtu
emiBDDCLastCompression.vtu
```

The compressed SDC+BDDC run writes:

```text
output/compression/sdc_bddc/emiSDCBDDCLastCompression.vtu
```

The compressed SDC+BDDC+AA run writes:

```text
output/compression/sdc_bddc_aa/emiSDCBDDCLastCompression.vtu
```

The current scratch program writes VTK solution output and solver progress to
standard output. Extended compression CSV, codebook, and ASCII exchange files
are not enabled by default yet.

## AA behavior

In ordinary SDC, AA selects active DOFs for later sweeps and expands that set
through the mass/stiffness matrix neighborhood. In SDC+BDDC, the selection is
converted to complete owner subdomains, since BDDC solves subdomain-local
problems. The residual assembly then visits cells incident to the active DOFs
or subdomains, preserving contributions to their residual rows. With AA off,
all DOFs or all subdomains remain active on every sweep.

## Numerical method and inputs

See [docs/emi-numerics.md](docs/emi-numerics.md) for the implemented weak-form
terms, SDC collocation equations and stopping rules, AA selection estimate,
BDDC partitioning, assumptions, output names, and validation recommendations.

The implementation uses the Aliev-Panfilov membrane model. MPI communication
is not implemented by this EMI driver; extending communication in Kaskade's
original `mg/bddc.hpp` is a separate planned task. The compression option
currently configures compressed transfer; extended compression logs and
exchange-file output are not enabled by default.
Histogram output can be enabled for compressed BDDC transfers with
`--bddcHistogram 1`. The CSV files are written below
`<output-directory>/histograms` and contain the raw, zigzag, and bit-length
symbol distributions for restriction and prolongation. For example:

```bash
./emiModel \
  --bddc 1 \
  --bddcCompression 1 \
  --bddcCompressionBits 16 \
  --bddcGraphLifting 1 \
  --bddcHuffman 1 \
  --bddcBitlength 1 \
  --bddcHistogram 1 \
  --dir output/emi-histogram
```

### Container checks before full-compression EMI scaling

Run the corrected uncompressed solver and full compression on the weak-8 and
fixed strong meshes before submitting a new scaling campaign:

```bash
CASES=weak_8_h16 BITS=16 RANKS=8 STEPS=1 bash tests/run_emi_scaling_smoke.sh
CASES=strong_h16 BITS=16 RANKS=8 STEPS=1 bash tests/run_emi_scaling_smoke.sh
```

The corrected CG path has passed ten-step checks at 16 and 24 bits on both
meshes at eight ranks and a true relative tolerance of `1e-7`. The iteration,
comparison and provenance records are in
`docs/emi_compression_validation_20261031.txt` and its companion JSON.
The 12-bit CG configuration remains nonconverged in the tested cases; use an
explicit `BITS=16` or `BITS=24` for these validated checks. The runner's
`BITS=12` default remains available as a diagnostic, not a validated setting.

For the small 10-cell mesh, run the three-stage uncompressed, quantization-only,
and full-compression check inside the container:

```bash
USE_CG=1 TOLERANCE=1e-7 BITS=16 RANKS=8 bash run_emi_10cells_compression_check.sh
```

This small-case runner defaults to `USE_CG=1`, `TOLERANCE=1e-7`, `BITS=24`, and
`RANKS=1`. `USE_CG=0` selects the existing steepest-descent iteration as a
separate diagnostic; it disables the CG-only convergence CSV but retains the
same true-residual acceptance check. `TOLERANCE` is forwarded to both the solver
and log validator. Alternative tolerances are separate experiments and are not
mixed into the `1e-7` campaign.

Each invocation creates a fresh output directory. The runner enables
quantization, integer graph lifting after quantization, zigzag mapping,
bit-length/tail coding and Huffman coding; corner DOFs remain exact.
It uses the same coefficients, diagonal weights, constraints and true global
relative residual target of `1e-7` in the paired runs. It records exact commands,
executable SHA256, source revisions, transfer counters, iterations, residuals,
phase timings and final exported-solution comparisons. Full logs stay in files;
the console prints short summaries. A failed MPI run or an unmet residual target
stops the test and writes a failed `metrics.json`; it is not treated as a timing
success. Lower precision is not made to converge by relaxing the tolerance.

`CASES` defaults to `weak_8_h16 strong_h16`, `RANKS=8`, `STEPS=1`, `BITS=12`, and
`MAX_ITERATIONS=300`. These are VTK-enabled validation runs, not paper-quality
timing measurements. The condition estimator is disabled because the compressed
path does not satisfy its fixed-SPD-preconditioner assumption. All profiler and
compression output is retained in `run.log` (`quiet=0`), while the legacy CG
table and per-update CSV are off in both runs. `STEPS=10` performs the follow-up
ten-step validation, with final time `0.1`; it still includes VTK output.

After validating, use controlled allocations, VTK off, repeated timings, and
paired uncompressed/full runs for the scaling tables. The report's existing weak
mesh/rank pairing is 8/8, 16/16, 32/32, 64/64 and 128/128. That family has two
material subdomains per rank, not one. Strong scaling keeps `strong_h16` fixed.
Changing only the rank count on weak-8 is not a weak-scaling experiment.
The existing historical Slurm scaling launcher uses the legacy internal
`bddcTolerance`; it needs a separate update to the current true-residual protocol
before launching this new campaign. Keep old measurements separate from new
runs made with the corrected material assembly and empty-interior handling.

Reported transfer bytes are the solver's aggregate counters; they are not a
hardware measurement of network traffic. `cgroup-memory-events-before/after.txt`
are cumulative OOM-counter snapshots, not per-run peak-memory measurements.
Review solution differences before deciding the precision for the final campaign.

The runner and log-validation checks can be tested without running EMI:

```bash
python3 -m unittest discover -s tests -p 'test_emi_scaling_smoke.py' -v
```

Regression checks for the numerical and reporting corrections:

```bash
python3 tests/test_bddc_residual_update.py
python3 tests/test_emi_compression_report.py
mpicxx -std=c++20 -O2 -I/kaskade7 tests/test_bddc_transfer_consistency.cpp -o /tmp/test_bddc_transfer_consistency
/tmp/test_bddc_transfer_consistency
```

The Python residual/report checks compile isolated production methods against
small matrix or MPI-rank stand-ins. The C++ transfer check uses the actual
Kaskade/Dune headers and checks shared-value consistency, unquantized corner
residuals, and byte counters. Container simulations remain the integration checks.
All ranks must use a rebuilt executable after the Kaskade corner-transfer fix.
The uncompressed MPI byte-reporting fix contributes globally aggregated domain
traffic only once; it changes reported baseline bytes, not solver arithmetic.

### Directional precision diagnostics for low-bit failures

Optional `--bddcRestrictCompressionBits` and `--bddcProlongateCompressionBits`
override `--bddcCompressionBits` in the corresponding transfer direction. Each
option defaults to `-1` (inherit the common setting), preserving existing runs;
explicit overrides must be in `1..32`, matching EMI's uint32 quantized symbols.
The setters are used consistently in non-SDC and SDC transfer configuration.

The 10-cell runner exposes these as `RESTRICT_BITS` and `PROLONGATE_BITS`.
Unset values inherit `BITS`. Each effective precision is printed and recorded in
the output path and exact command. For example:

```bash
USE_CG=0 BITS=8 RESTRICT_BITS=32 PROLONGATE_BITS=8 RANKS=1 TOLERANCE=1e-7 bash run_emi_10cells_compression_check.sh
USE_CG=0 BITS=8 RESTRICT_BITS=8 PROLONGATE_BITS=32 RANKS=1 TOLERANCE=1e-7 bash run_emi_10cells_compression_check.sh
```

These are mixed-precision diagnostic experiments, not uniform 8-bit validation.
32-bit quantization still differs from unquantized double transmission. This
option does not fix or establish convergence of uniform 8-bit or 12-bit CG.
Rebuild EMI after updating the options and precision helper. The helper itself
can be checked without Dune or MPI:

```bash
c++ -std=c++20 -O2 -I. tests/test_bddc_precision.cpp -o /tmp/test_bddc_precision
/tmp/test_bddc_precision
```

### One direct-reference and quantization-only sweep

For a fixed 10-cell, one-timestep diagnostic, run:

```bash
python3 tests/run_emi_quantization_diagnostic.py
```

This runs sequential MUMPS once, uncompressed BDDC+CG once, then uniform
quantization at 8, 12, 16, 24 and 32 bits. Both transfer directions use the listed
bit depth; graph lifting, zigzag, bit-length and Huffman stages are disabled.
The default is one MPI rank, one thread, 300 solver calls and true relative
tolerance `1e-7`. Existing directional-precision environment variables do not
override this uniform sweep. `--bits`, `--ranks`, `--iterations` and `--tolerance`
provide explicit diagnostic overrides. No solver default or physics is changed.

`--refine N` applies `N` uniform mesh refinement levels to every case, including
the direct reference; the default is `0`. This changes the mesh resolution and
is separate from `--restriction-passes`, which refines quantized residual updates.
The selected mesh refinement level is saved in `environment.json` and each
case's command. Mesh refinement increases assembly and factorization costs.

All converged BDDC results are compared with the direct reference. Quantized
results are additionally compared with the uncompressed BDDC baseline. The raw
exported-u differences are the primary displayed comparison; offset-removed
statistics remain separate diagnostics. PASS verifies residual/completion and
successful comparisons, without silently imposing a new solution-error budget.

A failed direct or uncompressed baseline stops the sweep. Failed quantized cases
are recorded while later bit depths continue, with a nonzero final exit code if
any case failed. One summary.txt/summary.json, exact commands, environment and
binary hash, complete run logs and convergence CSVs are retained below a fresh
output directory. This is a convergence diagnostic with VTK, not a timing study.
The runner's launch/validation control can be checked without real EMI or MPI:

```bash
python3 -m unittest discover -s tests -p 'test_emi_quantization_diagnostic.py' -v
```

### Experimental refinement of low-bit restriction updates

`--bddcRestrictionPasses` defaults to `1`. Values `2..16` enable additional
restriction exchanges before a compressed non-SDC BDDC update. Between these
exchanges, the original local residual and solution stay fixed. The existing
published-residual feedback sends the remaining communication error each time,
so additional low-bit updates refine the restricted residual before the coarse
solve. No matrix, boundary coefficient or stopping tolerance is changed.

This is a separate communication/solver configuration: each exchange quantizes
values at the configured bit depth, but additional messages carry additional
information. The quantization-only diagnostic still uses its existing unencoded
transport format; quantization bits do not specify the physical wire width.
Traffic counters and restriction timing include every exchange. An 8-bit run
with four passes must be reported as 8-bit restriction refinement, not as the
original one-pass experiment. Neither speedup nor net byte savings is implied.
The application currently enables extra passes only for compressed non-SDC runs.

A focused experiment, after rebuilding the updated Kaskade/EMI code, is:

```bash
python3 tests/run_emi_quantization_diagnostic.py --bits 8 --restriction-passes 4
```

The direct and uncompressed CG baselines remain unchanged and are run once each.
The quantized case uses 8 bits in BOTH transfer directions, with four restriction
passes, and no graph or entropy coding. The summary explicitly records `Rpasses`.
The default driver still uses one pass. Passing the transfer-kernel regression
below does not establish EMI convergence; the container experiment is required:

```bash
mpicxx -std=c++20 -O2 -I/kaskade7 tests/test_bddc_restriction_refinement.cpp -o /tmp/test_bddc_restriction_refinement
/tmp/test_bddc_restriction_refinement
```

The diagnostic driver can optionally add the full coding pipeline after each
quantization-only case, reusing the same direct and uncompressed references:

```bash
python3 tests/run_emi_quantization_diagnostic.py --bits 8 --restriction-passes 4 --full-compression --ranks 8
```

`--full-compression` enables integer graph lifting after quantization, zigzag,
bit-length/tail coding and Huffman in the `full-*` rows only. `quantization-*`
rows keep all those stages disabled. The default remains quantization only.
The full case is attempted even if its quantization-only counterpart fails.

The `net bytes saved` column is `100*(1-case_bytes/uncompressed_CG_bytes)` using
total logical transferred bytes over the entire solve. It includes changed
iteration counts and every additional restriction pass. This differs from the
codec's saving against its own equivalent raw stream; a positive per-message
saving can coexist with negative net savings when refinement adds exchanges.
These remain application counters including intra-rank transfers, not measured
network traffic. Accuracy, net byte reduction and wall time are separate results.

The one-timestep 10-cell validation of this four-pass 8-bit variant passed with
CG on eight MPI ranks. Uncompressed CG used 10 iterations; quantization-only and
full coding each used 11, with a true relative residual of `4.56e-8` and raw
exported-u relative difference `1.04e-9` versus sequential MUMPS. Full coding
transferred 268,112 logical bytes versus 442,240 for uncompressed CG: **39.37%
net savings**, including all extra exchanges. The codec's 77.95% saving against
its own raw stream is a different metric. This validates one timestep on this
small mesh; larger meshes, multiple timesteps and speedup remain unverified for
this variant. The original one-pass 8-bit failure is retained separately.
Detailed provenance is in `docs/emi_8bit_refinement_validation_20261031.txt` and
its companion JSON.

## Optional EMI CG communication diagnostics

Use `--bddcCgDiagnostics 1` (default 0), or pass `--cg-diagnostics` to
`tests/run_emi_quantization_diagnostic.py`. This is available for non-SDC BDDC CG
with the original true-residual stopping rule. It does not change quantization,
restriction passes, alpha/beta, the preconditioner, or the convergence tolerance.
Rebuild EMI against the matching Kaskade headers before enabling the option.

Each timestep writes `bddc-cg-diagnostics-stepN.csv`. Restriction rows are
recorded after **every pass**, before the coarse solve. Search rows are recorded
before updating the solution, with an `applied` flag filled after the solve call.
Solver calls are one-based; their true residual is the **pre-update** residual.
Use the existing convergence CSV for accepted post-update residuals.

Definitions:

- Original residual: `r = b - A*x`, using the original assembled EMI matrix,
  current RHS, and weighted global correction. Root evaluates global norms once.
- Restriction error: the Euclidean norm of `sum_s rRes_s - r`. Each owned
  subdomain contributes once; the restricted dual residual is not weighted a
  second time. Relative error divides by `||r||_2`; for zero residual it is absolute.
- Reference CG scalar: `rho_true = r^T*z`, where `z` is the actual prolonged
  preconditioned correction, including its configured quantization. This is
  compared with the solver's existing `rho_internal`, not substituted for it.
  Relative discrepancy is blank if the reference scalar is too small
  (`|rho_true| <= 64*epsilon*||r||*||z||`) or invalid; the absolute discrepancy
  and `rho_true_near_zero` flag remain available.
- Adjacent conjugacy defect: `|q_k^T A q_previous| /
  sqrt((q_k^T A q_k)*(q_previous^T A q_previous))`. The previous direction is
  the last **accepted** direction. The first sample and nonpositive energies
  yield a blank value, not a false zero. This checks adjacent directions only.
- Internal curvature is the value used by CG. Global curvature is independently
  evaluated with the original assembled matrix.
- Shared-copy jumps for `x`, `z`, and `q` are the owner-local Euclidean norm
  of deviations from their weighted global reconstruction, divided by
  `max(1,||global vector||)`. Nonzero jumps must be considered before interpreting
  scalar/conjugacy measurements of the reconstructed vector.
- `solution_change_since_restriction` is the owner-local Euclidean norm of
  solution changes since the preceding restriction observation in the same call.
  The first pass defines the reference. It should remain exactly zero.

All diagnostic observations are read-only. They add two sum reductions per event,
root matrix-vector products, memory, and CSV output. File-status checks also
broadcast on MPI runs. These costs are **not** included in transfer byte counters,
and diagnostic timings are unsuitable for compression-performance claims.
Default-off execution adds no diagnostic reductions or matrix-vector products.

Start with the same 10-cell EMI input, refinement 2, 8 ranks, one timestep,
unchanged physical parameters, CG, and tolerance `1e-7`:

`RANKS=8 REFINE=2 bash tests/run_emi_cg_diagnostic_check.sh` runs the complete
comparison. It first requires identical convergence CSVs with diagnostics off/on
for uncompressed CG and 16-bit quantization. It then retains a failing one-pass
12-bit run and proceeds to two passes, provided both reference solves succeeded.
The script does not compile; rebuild first. Alternatively, run the cases directly:

```bash
python3 tests/run_emi_quantization_diagnostic.py \
  --bits 12 16 --restriction-passes 1 --ranks 8 --refine 2 \
  --tolerance 1e-7 --cg-diagnostics
python3 tests/run_emi_quantization_diagnostic.py \
  --bits 12 --restriction-passes 2 --ranks 8 --refine 2 \
  --tolerance 1e-7 --cg-diagnostics
python3 tests/summarize_emi_cg_diagnostics.py output/<diagnostic-output-directory>
```

The first command is expected to exit nonzero if the single-pass 12-bit case
still fails; its diagnostics are retained and the 16-bit case is still attempted.
Graph lifting and entropy coding stay disabled unless `--full-compression` is
explicitly requested. Before drawing conclusions, compare an instrumented
16-bit run's convergence CSV with an otherwise identical default-off run.

Lightweight checks (not a replacement for container MPI validation):

```bash
python3 tests/test_bddc_cg_observer.py --kaskade-root /kaskade7
c++ -std=c++20 -O2 tests/test_bddc_cg_diagnostics.cpp -o /tmp/test_bddc_cg_diagnostics
diagnostic_test_dir=$(mktemp -d /tmp/emi-cg-diagnostics.XXXXXX)
/tmp/test_bddc_cg_diagnostics "$diagnostic_test_dir"
```

### Experimental single-pass, zero-preserving restriction

`--bddcRestrictZeroPreserving 1` selects an **optional restriction-only** quantizer.
The default (`0`) retains the existing quantization arithmetic. CG scalars,
recurrences, stopping rule, restriction pass count and prolongation quantization
are unchanged. This experiment does not establish that 12-bit CG converges.

The candidate uses an unsigned b-bit alphabet with an implicit integer zero point
`Z=2^(b-1)` and symmetric limit `L=Z-1`. For each subdomain's residual increments,
choose `s <= L/max(abs(increment))`, round the scale down to float32, encode
`q=Z+round(s*increment)` and decode `(q-Z)/s`. Signed values are rounded before
adding Z. The code range is `1..2^b-1`, leaving one code unused, so exact zero and
sign symmetry are preserved. Sender feedback and received values use exactly the
same float32 scale. The existing two-float header contains `(s,0)`; there is no
additional metadata or exchange. All ranks must use the same rebuilt executable,
bit depth and quantizer mode. Supported depths are 2..32 bits. Nodal quantization
and optional integer lifting **after** quantization are supported; transforms
before quantization are rejected for this experiment. Extreme data requiring a
scale outside the supported float32 range are rejected rather than silently
changing precision.

This symmetric choice can have coarser spacing than the original affine min–max
grid on an asymmetric range. It does not preserve weighted sums, change the
quantization group, or add a second error-correction payload. Fully encoded size
may improve or worsen; quantization-only still uses double-sized symbol slots.

Run the controlled experiment in the existing container after synchronizing both
this repository and `kaskade7_new`:

```bash
(
set -e
mpicxx -std=c++20 -O2 -I/kaskade7 tests/test_bddc_zero_quantization.cpp \
  -o /tmp/test_bddc_zero_quantization
/tmp/test_bddc_zero_quantization
mpicxx -std=c++20 -O2 -I/kaskade7 tests/test_bddc_transfer_consistency.cpp \
  -o /tmp/test_bddc_transfer_consistency
/tmp/test_bddc_transfer_consistency
make -W EMI_model.cpp -o /kaskade7/libs/libkaskade.a \
  MPI=1 MPICXX=mpicxx \
  MPI_INCLUDE_FLAGS="-include /usr/lib64/mpi/gcc/openmpi4/include/mpi.h -DMUMPS_MPI_H"
python3 tests/run_emi_zero_preserving_check.py
)
```

The driver fixes the original 10-cell mesh with refinement 2, 8 BDDC ranks,
one timestep and relative tolerance `1e-7`. It runs matched direct/uncompressed
references, legacy 16-bit/one-pass, and candidate 12-bit/one-pass, with both
quantization-only and full coding. The candidate is 12 bits in **both** directions,
but only restriction uses the new quantizer. Output is placed in a new
`output/emi-zero-preserving-*` directory; existing runs are never overwritten.
The direct and uncompressed cases never receive the new option.

`comparison.json` declares the objective met only when the candidate converges,
its raw exported-u relative Euclidean and maximum absolute errors versus direct
are no worse than the matched 16-bit baseline, and its total fully encoded logical
bytes beat both that baseline and the historical **2,430,000-byte** target.
Quantization/full convergence histories must agree. No mean offset is removed
for this criterion. CG diagnostics remain enabled to explain success or failure;
their communication is outside the codec counters, and their times are **not**
performance measurements. The driver does not silently increase bit depth,
restriction passes, iteration limit or tolerance after a failure.

For the candidate alone, the existing runner also accepts:

```bash
python3 tests/run_emi_quantization_diagnostic.py --bits 12 \
  --restriction-passes 1 --zero-preserving-restriction \
  --full-compression --cg-diagnostics --ranks 8 --refine 2 --tolerance 1e-7
```

### Read-only restriction range analysis before grouped quantization

Existing aggregate diagnostics do not retain buffer values. After rebuilding,
capture the existing 16-bit and zero-preserving 12-bit trajectories with:

```bash
python3 tests/run_emi_restriction_range_check.py \
  --reference output/emi-zero-preserving-20261101T081725Z
```

This runs quantization-only diagnostics (with direct and uncompressed references)
at the same 10-cell mesh, refinement 2, eight ranks, one time step, physics, one
restriction pass, 300-call limit and true relative tolerance 1e-7. The known
12-bit failure is retained. No new quantizer is applied, and no full-codec rerun
is needed to observe these increments. Every convergence CSV must match its
saved counterpart byte-for-byte before range results are interpreted.

The opt-in `--restriction-range-capture` runner flag exports
`KASKADE_RESTRICTION_CAPTURE` separately for each compressed case. Kaskade writes
the actual pre-quantization increment (local residual minus published feedback),
exact-corner exclusions and ordered neighbor membership. Capture occurs before
encoding/feedback changes, without MPI, at restriction invocations 1..32, 50,
100, 200 and 300. This is sampled data, not a complete failed-run byte count.
Files are distinct per owned subdomain and invocation; existing files are never
overwritten. Capture I/O invalidates timing comparisons. The low-level capture
counter is per transfer object and does not reset at time steps.

`range-analysis/buffers.csv` and `groups.csv` measure ranges for groups of
32, 64, 128, 256 and 512 compressible entries in producer interface order.
One producer DOF has one group, reused across all neighbors and its local copy.
Naively quantizing independently per neighbor would violate this requirement.
These groups use index order, not inferred spatial connectivity.

The proposed metadata model is 12 bytes per contiguous group segment (uint32
length and two float32 scale/offset fields), 4 bytes per packet for segment
count, the existing 8 bytes of stream lengths, and 8-byte packet alignment.
Receiver DOF order is already known; no per-value group IDs are assumed.
Fragmented groups pay repeated segment descriptors. Exact-corner doubles remain
unchanged. The analyzer explicitly counts these segments from captured topology.

For full groups without fragmentation, descriptor overhead is 0.375, 0.1875,
0.09375, 0.046875 and 0.0234375 bytes/value respectively. Reducing a fixed-width
payload from 16 to 12 bits offers only 0.5 bytes/value before metadata. Small
or fragmented groups can therefore erase that saving. The reported packed-bit
model is **not** a prediction of integer lifting/Huffman sizes: changing scales
may change symbol smoothness and entropy. Quantization-only double slots do not
shrink when fewer quantization bits are selected.

The numerical indicator compares the ideal symmetric 12-bit group step
`maxabs(group)/2047` with the original affine 16-bit whole-buffer step
`span(buffer)/65535` on the same captured values. It is a screening bound,
not a measured quantization error or a CG convergence guarantee. Constant
whole buffers have undefined relative step ratios. Offline replay cannot predict
how changed quantization will alter future feedback and search directions.

After reviewing these measurements, select one group size for a controlled
experiment: grouped zero-preserving **12-bit restriction**, original **16-bit
prolongation**, one pass, same CG/physics/tolerance. Keep exact corners and full
coding settings fixed. Do not implement this candidate until the range and
metadata evidence supports it. Judge total logical bytes against the matched
16/16 baseline (historically 2,430,000), plus convergence and raw solution
differences against direct. Per-packet savings alone are insufficient.

### Controlled grouped restriction experiment (128 values, R12/P16)

The optional `--bddcRestrictGroupSize 128` implementation now enables this
experiment with `--bddcRestrictZeroPreserving 1`. Zero remains the default group
size, retaining the previous path. Groups count compressible DOFs in sorted
producer interface order; exact corners are skipped and remain unquantized.
Per-group maximum magnitudes determine symmetric float32 scales. Own feedback
and every neighbor use the same code and scale for a producer value.

Neighbor packets retain their existing DOF order. A uint32 run count followed
by (uint32 scalar count, float32 scale, float32 zero offset) descriptors replaces
the previous scale/min header. The descriptor header is independently padded
to an 8-byte slot boundary; this can cost up to four bytes more per packet than
the earlier combined-padding model. Encoded stream lengths and payloads follow.
Integer lifting/Huffman operate on the grouped codes and remain reversible;
changing scales can change their entropy. Exact-corner trailers are unchanged.
Sent/received byte counters include all descriptors, padding and payloads.
There are no additional MPI exchanges. Packet repacking and receiver scale
expansion add allocation, copying and computation overhead.

After rebuilding, run:

```bash
python3 tests/run_emi_grouped_restriction_check.py
```

This preserves the fixed 10-cell/refinement-2 problem, eight BDDC ranks, one step,
one restriction pass, 300-iteration limit, physics and tolerance 1e-7. It runs
direct/uncompressed references, quantization-only and full-codec variants for:

- Baseline: original restriction/prolongation quantizers, 16/16 bits, no groups.
- Candidate: grouped exact-zero restriction at 12 bits, original prolongation
  at 16 bits, group size 128.

The baseline must reproduce the preserved convergence CSVs in
`output/emi-zero-preserving-20261101T081725Z/baseline16`; use `--reference`
to specify the corresponding preserved baseline directory elsewhere.
Outputs go to a fresh `output/emi-grouped-restriction-*` directory. No existing
results are changed. `comparison.json` requires true-residual convergence,
raw exported-u relative Euclidean and maximum absolute errors no worse than the
matched full-16/direct comparison, and fewer total logical bytes than both the
matched baseline and 2,430,000 bytes. Quantization/full convergence CSVs must
agree, the executable hash must remain unchanged, and precision/pass settings
are checked explicitly. A solver PASS alone is not an objective PASS.

The diagnostic runner also accepts `--prolongation-bits 16` and
`--restriction-group-size 128`; these overrides affect compressed cases only
and are recorded explicitly in the summary. This experiment changes restriction
quantization and wire metadata, not CG. A successful outcome has not yet been
established. Diagnostic timings include observation and are not benchmark claims.

### Multi-step validation of ungrouped zero-preserving R12/P16

The ungrouped R12/P16 configuration passed the refinement-2, eight-rank,
one-step 10-cell test with 22 iterations and 2,213,992 total logical bytes,
versus 21 iterations and 2,430,000 bytes for original R16/P16. Its final raw
relative Euclidean and maximum absolute differences from direct were
9.208e-7 and 1e-7 (baseline 2.820e-6 and 3e-7). This is a single-problem result,
not a default change or general robustness/performance claim.

Validate the same configuration for ten steps on all three existing problems:

```bash
python3 tests/run_emi_mixed_precision_validation.py --steps 10 --ranks 8
```

No rebuild is needed after the successful grouped/control builds; this adds
Python validation infrastructure only. Defaults are 10cells at refinement 2,
weak_8_h16 at refinement 0, and strong_h16 at refinement 0. The block meshes
come from `../input_emi_mesh_blocks/weak_strong_scaling`; override with
`--mesh-dir` if necessary. To run only the first problem:

```bash
python3 tests/run_emi_mixed_precision_validation.py --cases 10cells --steps 10 --ranks 8
```

For each mesh the script runs direct (one rank), uncompressed CG, full original
R16/P16, and full zero-preserving R12/P16. References run only once per mesh.
The candidate always uses one restriction pass, no groups, the original
prolongation quantizer and unchanged CG. Physical parameters remain sigma_i=3,
sigma_e=20, C_m=1, R=0.1, dt=0.01; finalTime=steps*dt. All BDDC steps must
meet the true relative tolerance 1e-7 with at most 300 solver calls per step.
Direct residuals are checked against 1e-7 at every step as well.

The runner reuses `collect_emi_scaling_smoke.py` and the diagnostic runner's
coordinate/finite-error comparison helper. It does not repeat quantization-only
codec tests. CG observation and buffer capture are disabled; per-step true
residual reports and convergence CSVs remain enabled. Final VTK output is used
for raw accuracy comparisons, so intermediate-step solution differences are
not evaluated. Timing includes setup, residual checks and output and is not
a scaling benchmark.

Each mesh has its own matched 16/16 byte budget and final raw error thresholds.
The objective passes only if all solves complete, every step converges, the
candidate's total logical bytes are strictly lower, and both its final relative
Euclidean error and maximum absolute error versus direct are no worse than
16/16. The historical 2,430,000-byte target is NOT applied to different meshes
or ten-step runs. Mean-offset removal never affects the acceptance criterion.

Results, exact commands, executable hash, fixed parameters, logs, comparisons,
per-step iteration counts/residuals and byte metrics are saved in a fresh
`output/emi-mixed-validation-*` directory. `summary.json` and `summary.txt`
distinguish solver checks from the complete accuracy/communication objective.
A failed reference skips the remainder of that mesh; subsequent meshes are
still evaluated. A candidate failure or worse accuracy/traffic gives a nonzero
final exit status while retaining all outputs. No precision or tolerance
fallback is applied, and existing output directories are never reused.

The ten-step validation `emi-mixed-validation-20261101T085314Z-5184`
did not establish a generally successful R12/P16 configuration. The 10-cell
case converged and saved 8.04% of logical bytes but failed the raw final accuracy
comparison. The weak_8_h16 candidate failed at step 6: its best relative residual
was 4.6333e-7 at iteration 117 and its final residual was 127.678 at iteration
300. The strong_h16 candidate passed the objective with 3.27% fewer logical
bytes. Keep these outcomes separate; do not generalize the earlier one-step
success to all meshes or time steps.

Replay only the failed run with read-only CG observations, preserving its exact
stored command except for output directory and `--bddcCgDiagnostics 1`:

```bash
python3 tests/replay_emi_cg_diagnostics.py \
  --run output/emi-mixed-validation-20261101T085314Z-5184/weak_8_h16/full-zero-r12-p16
```

Do not rebuild for this replay: its executable hash must match the saved run.
The same solver failure is expected. The script checks byte-identical
convergence CSVs and the same success/failure outcome before summarizing the
CG observations. A successful diagnostic replay does NOT mean the solver
converged. New logs and `replay.json` are saved under a fresh
`output/emi-cg-replay-*` directory. No quantizer, metadata precision, physical
parameter, iteration limit, pass count or tolerance is changed.

### Experimental CG scalar from the original local residual

`--bddcCgOriginalResidual 1` changes the CG numerator from
`sum_i z_i^T rRes_i` to `sum_i z_i^T r_i`. In compressed CG, `r_i` is
already recomputed as `f_i-A_i*u_i` after each update. The default is zero,
preserving the previous behavior. The option is restricted to non-SDC BDDC CG
with a positive true-residual tolerance and no operator audit.

This is a numerical experiment: alpha and beta can change. It adds no exchange
or reduction and keeps the quantizer, payload layout, float32 scale and float32
minimum, one restriction pass, and CG recurrence formulas unchanged.
It does not make a nonlinear/history-dependent compressed preconditioner fixed
SPD or guarantee convergence. Original-local and assembled diagnostic rho
should be compared while also monitoring shared-copy continuity.

The weak-mesh replay showed that loss of accuracy and conjugacy preceded the
growth of shared-copy disagreement. At the best step-6 residual (iteration
117), the solution jump was about 3.1e-16; the much larger late jumps alone
therefore do not explain the onset of failure.

After rebuilding `emiModel` against the updated Kaskade headers, run a controlled experiment:

```bash
python3 tests/test_bddc_original_rho.py --kaskade-root /kaskade7
python3 tests/run_emi_mixed_precision_validation.py \
  --cases weak_8_h16 --steps 10 --ranks 8 \
  --candidate-original-rho --candidate-cg-diagnostics \
  --reference-campaign output/emi-mixed-validation-20261101T085314Z-5184
```

The direct, uncompressed CG and original 16/16 controls are rerun with the new
binary; only the zero-preserving R12/P16 candidate enables the new scalar and
observer. Saved commands must match except for output, observer and scalar
switches. Uncompressed and 16/16 convergence CSVs must be byte-identical to the
saved controls before running the candidate. Outputs go to a fresh directory.
Final raw L2/max error and total logical bytes are compared against these newly
run controls with the existing acceptance criteria; every step must meet 1e-7.
Candidate diagnostics add operator work and reductions, so these timings must
not be used to claim a performance improvement.

Local verification compiles the production `prolongate()` body in a small
vector/matrix harness, testing original/perturbed residuals, shared copies,
constraints, first/later CG calls and unchanged steepest descent behavior.
Mock runner tests check candidate-only flags, history parity and preservation
of parameters. Actual multi-step convergence remains a container integration
test, not a result established by these small tests.

For a matched corrected-rho comparison, use:

```bash
python3 tests/run_emi_mixed_precision_validation.py \
  --cases weak_8_h16 --steps 10 --ranks 8 \
  --matched-original-rho --matched-cg-diagnostics \
  --reference-campaign output/emi-mixed-validation-20261101T091721Z-5741
```

Both compressed runs now use `--bddcCgOriginalResidual 1` and the same
observer setting. Direct and uncompressed controls retain their original
settings. This runner extension requires no rebuild of a binary that already
supports `bddcCgOriginalResidual`. Existing output directories are preserved.
The uncompressed and R12/P16 convergence histories must reproduce the saved
run because their scalar settings are unchanged. The 16/16 history may change
because its scalar setting changes intentionally; this is recorded explicitly.

The prior candidate-only scalar experiment converged in all ten steps but
required 435 iterations versus the legacy 16/16 baseline's 365. Its logical
bytes were 12,601,784 versus 11,957,688, and both final raw accuracy measures
were worse. The matched experiment measures performance of the two compression
configurations with the same scalar formulation; it does not isolate bit depth
alone, since exact-zero restriction remains enabled only for R12/P16.
The acceptance criteria remain fewer total logical bytes and no worse final raw
L2/max errors against direct, with every step satisfying 1e-7. Equal observer
settings remove one timing asymmetry, but these remain instrumented validation
runs, not a timing benchmark.

To test the next intermediate precision, R14/P16, with all other candidate
settings held fixed:

```bash
python3 tests/run_emi_mixed_precision_validation.py \
  --cases weak_8_h16 --steps 10 --ranks 8 \
  --candidate-restriction-bits 14 \
  --matched-original-rho --matched-cg-diagnostics \
  --reference-campaign output/emi-mixed-validation-20261101T092402Z-5920
```

This is a Python runner extension; the current corrected-rho binary needs no
recompilation. The default candidate remains R12/P16. The new mode directory
and summary key are `full-zero-r14-p16`, and the selected restriction precision
is recorded in `environment.json` and the saved command. Prolongation stays at
16 bits, restriction stays zero-preserving with one pass and no groups, and
scale/minimum metadata remain float32. Both compressed runs use corrected rho
and the same diagnostics. The existing convergence and raw accuracy/byte
acceptance checks remain unchanged.

When referencing the prior R12 campaign, only the candidate restriction bit
depth is additionally allowed to differ. The uncompressed and corrected-rho
16/16 histories must reproduce the saved histories. Candidate history parity
is intentionally skipped when restriction precision differs. All previous
results are preserved; R14 convergence and efficiency remain to be measured.

### Potential reference and per-step mean-error investigation

Source inspection of `EMI_model.hh` shows an extracellular outer-boundary
term with hard-coded `gamma = 1e-5`: its residual is `-gamma*u_e*v`.
The semi-implicit Euler wrapper assembles `A = B - dt*L'`, so this contributes
`+dt*gamma` times a boundary mass matrix. Bulk terms depend on gradients and
inner-boundary terms on voltage differences. Consequently a common constant
shift is anchored by the extracellular boundary term, not by an exact mean-zero
constraint. The exposed `penalty` value is stored by EMIModel but is not used
in these weak-form terms. Changing it is not a remedy for this observation.
No boundary coefficient, constraint or physical parameter has been changed.

Provided the relevant connected component touches the extracellular boundary,
the common constant mode is therefore not an exact null mode. Its stiffness
must be measured from the assembled matrix; `dt*gamma = 1e-7` in these runs is
only the boundary-term multiplier, not an eigenvalue. The combination of weak
anchoring and an unweighted residual tolerance can permit a visible common
offset. This is an explanation to test, not proof that every observed error
is caused by that mode or that the offset is physically irrelevant.

Existing comparison outputs use nonconforming VTK samples written with
`setPrecision(7)` in `writeState`. Their sample arithmetic mean is neither a
unique-DOF mean nor a volume-weighted mean. In campaign
`emi-mixed-validation-20261101T094154Z-6214`, the constant sample offset
accounts for approximately 96.40% of squared R14 error on 10cells and 99.995%
on strong_h16, using `1-(centered_L2/raw_L2)^2`. These are final-state
exported-data diagnostics, not full-precision solution claims.
The direct path saved per-step VTK states; BDDC saved only initial/final states.
Existing files cannot reconstruct the BDDC mean-error time history.

`--potentialDiagnostics 1` adds read-only rank-zero CSV snapshots after every
accepted non-SDC update, plus the initial state and operator row/column sums.
It writes double coefficients at 17 significant digits, with correction, RHS
and independently evaluated residual. It never shifts the mean or edits a
solver vector. The disabled path does no diagnostic work. Enabled diagnostics
add matrix traversal, file I/O and a per-snapshot status broadcast, but no
compressed transfer payloads; do not use their timings for scaling claims.
Existing diagnostic files are not overwritten.

Rebuild EMI, then rerun the saved R14 campaign with:

```bash
python3 tests/run_emi_mixed_precision_validation.py \
  --cases 10cells strong_h16 --steps 10 --ranks 8 \
  --candidate-restriction-bits 14 \
  --matched-original-rho --matched-cg-diagnostics \
  --potential-diagnostics \
  --reference-campaign output/emi-mixed-validation-20261101T094154Z-6214
```

All BDDC convergence histories must match the saved histories, establishing
observer parity for the numerical trajectories. Existing raw VTK acceptance
criteria remain in force; the campaign can still report NOT MET. The runner
also creates `potential-analysis/per-step-errors.csv`,
`constant-mode.csv` and `interpretation.json`, even when the completed
campaign misses the accuracy/byte objective.

The additional report compares unique algebraic DOF states against direct at
every step. It records raw and centered error, signed mean error, the fraction
of squared error associated with the mean, correction mean differences and
their accumulation. These norms are explicitly distinct from VTK sample norms
and physical FE norms. Centered errors never replace raw acceptance tests.

Operator column sums provide `1^T*A*u` and an anchor-weighted mean-error proxy.
The report checks `1^T*b - 1^T*A*du - 1^T*r` for roundoff closure.
It also reports `1^T*b + 1^T*A*u_previous` and
`1^T*A*u_new + 1^T*r`: their agreement follows from the update identity,
and their magnitude tests the expected cancellation of bulk/interface terms
under the constant test function. The constant Rayleigh quotient is a probe,
not an eigenvalue estimate. Membrane-voltage and gradient accuracy are not
inferred from the centered potential error.

### Opt-in constant-mode correction of the converged linear solve

`--bddcConstantCorrection 1` applies one scalar energy correction after CG has
already met the unchanged true-residual tolerance. For the assembled correction
`x=du` and its current timestep RHS, it computes
`delta = (1^T*(b-A*x))/(1^T*A*1)`, then adds `delta` to every global correction
DOF. It does not subtract the arithmetic mean of the solution, reference a
direct solution, or change the boundary anchoring. The default is zero.

The proposal is computed on rank zero using long-double accumulation, broadcast
to all ranks, and applied to the replicated global correction. The full original
residual norm is then recomputed and broadcast. The timestep is advanced only
if this corrected iterate still meets the same tolerance. A nonpositive or
nonfinite constant stiffness, invalid delta, or failed residual recheck aborts
without advancing the state. No silent fallback, extra restriction pass,
precision escalation or relaxed tolerance is used. An energy correction can
increase the Euclidean residual, which is why the recheck is essential.

The correction is outside the CG recurrence: existing convergence and CG
diagnostic CSVs describe the pre-correction CG trajectory. Separate
`BDDC constant correction step ...` log entries contain delta, residual sum,
constant stiffness, and pre/post relative residuals. Global completion summaries
and full-precision potential snapshots describe the corrected accepted iterate.
Local solver iterates are discarded at the next RHS restart; they are not used
again after this final correction. Later timestep trajectories may change.

For the controlled R14/P16 versus R16/P16 comparison, both compressed runs enable
the same correction and original-residual rho; direct and uncompressed controls
retain their previous settings:

```bash
python3 tests/run_emi_mixed_precision_validation.py \
  --cases 10cells strong_h16 --steps 10 --ranks 8 \
  --candidate-restriction-bits 14 \
  --matched-original-rho --matched-cg-diagnostics \
  --matched-constant-correction --potential-diagnostics \
  --reference-campaign output/emi-mixed-validation-20261101T095708Z-6582
```

Rebuild EMI for this C++ change. First-step CG histories must match the saved
pre-correction histories; full-history parity is required only for modes whose
numerical settings are unchanged. Fresh output directories preserve all prior
results. Both raw VTK acceptance criteria and full-precision mean diagnostics
remain available. The saved `20261101T101620Z-6974` container campaign completed
all ten steps for both compressed configurations on `10cells` and `strong_h16`.
Its strict VTK accuracy objective passed only for `10cells`; convergence and
communication savings must be distinguished from that accuracy comparison.

The scalar proposal broadcast contains two doubles (validity and delta); the
extra true-residual measurement broadcasts two more. The runner adds
`32*(ranks-1)*steps` logical root-to-peer payload bytes to each corrected run's
compressed-transfer count. At eight ranks and ten steps this is 2,240 additional
bytes. Raw codec counters are preserved separately. This is logical delivery
accounting, not actual MPI wire traffic; other validation reductions and framing
remain outside the existing counters. Correction work includes a matrix
traversal for the proposal, a vector update and the true-residual recheck.

`tests/test_emi_constant_correction.cpp` checks the sign/formula, read-only
proposal, preservation of differences, weak grounding, null/nonfinite mode
rejection, and a case where the Euclidean residual increases. Runner tests
check matched flags, post-correction residual checks, payload accounting and
first-step parity. No quantizer or float32 scale/minimum metadata is modified.

### Single-pass low-bit restriction ladder using saved controls

The next controlled experiment lowers restriction precision to 12, 10 and 8
bits while holding prolongation at 16 bits. It reuses the direct, uncompressed
and 16/16 controls from the same executable and changes only restriction bits
and the output directory in each saved R14/P16 command. Corrected rho,
zero-preserving restriction, one pass, no groups, post-CG scalar correction,
physical parameters and the true relative tolerance `1e-7` remain fixed.
Scale and minimum metadata remain float32. This is mixed precision, not
uniform 8-bit communication.

Run in the container without rebuilding the executable used by this campaign:

```bash
python3 tests/run_emi_low_bit_one_pass.py \
  --reference-campaign output/emi-mixed-validation-20261101T101620Z-6974 \
  --cases 10cells strong_h16 \
  --bits 12 10 8
```

The runner verifies the executable hash, creates fresh output, runs six
candidates, and continues after a failed candidate. It does not rerun controls,
change precision automatically, add refinement passes or relax tolerance.
`PASS` requires every timestep to converge and finite final comparisons to
complete. Byte savings and solution errors are separate findings: a slightly
higher error than the 16-bit reference is no longer a convergence veto.
The former strict objective is still recorded in JSON, and the original mixed
validation runner and historical results retain their acceptance criteria.

Full-precision final DOF errors include raw and centered relative L2 error,
maximum error and mean offset. Centering is an analysis only; no additional
solution shift is applied by this runner. These diagnostics do not establish
accuracy of membrane voltages or gradients. Logical byte totals include the
existing constant-correction broadcasts; reused controls are not fresh timing
measurements. Successful 14-bit results do not establish success at 12 or 8
bits: that is what this experiment tests.

`python3 tests/test_emi_low_bit_one_pass.py` checks command isolation, control
reuse, continuation after failure, separate accuracy reporting, preservation of
reference results, and rejection of a changed executable using mock runs.

To isolate prolongation precision after the successful R12/P16 ladder, run:

```bash
python3 tests/run_emi_low_bit_one_pass.py \
  --reference-campaign output/emi-mixed-validation-20261101T101620Z-6974 \
  --comparison-campaign output/emi-low-bit-one-pass-20261101T103116Z-7366 \
  --cases 10cells strong_h16 \
  --bits 12 --prolongation-bits 12
```

This adds only two R12/P12 solves. Direct, uncompressed, R16/P16 and R12/P16
results are reused. The comparison campaign must use the same reference
environment and executable; its successful commands must match after changing
only prolongation precision and the output directory. The summary reports
total logical byte savings against both R16/P16 and saved R12/P16, together
with iteration counts and final full-precision errors. Failure is retained and
does not trigger a higher-precision retry. No compilation is needed.

Both directions now use 12 quantization bits, but their existing quantizer
policies remain unchanged: restriction is zero-preserving and prolongation
uses its original quantizer. Scale/minimum metadata remain float32, and exact
corner values retain their existing representation. This experiment tests
equal bit depths, not a change to the numerical algorithm or wire format.

### Optional flexible CG for one-pass R12/P12

`--bddcFlexibleCg 1` changes the search-direction recurrence to one-term
truncated flexible CG (Polak--Ribiere). It requires `--bddcCgOriginalResidual 1`,
non-SDC BDDC CG and a positive true residual tolerance. Standard CG remains the
default. The fixed-preconditioner condition estimate is rejected in this mode.

For original residual `r_k=b-A*x_k` and communicated/preconditioned correction
`z_k`, the new coefficient is
`beta_k = z_k^T*(r_k-r_(k-1)) / (r_(k-1)^T*z_(k-1))`.
The direction is `q_k=z_k+beta_k*q_(k-1)` and the existing step length remains
`alpha_k=(r_k^T*z_k)/(q_k^T*A*q_k)`. The first update is unchanged. With exact
line searches this coefficient makes adjacent directions A-conjugate for a
changing preconditioner. Negative beta is allowed; nonfinite beta is rejected.
The usual positive-rho/curvature and true-residual checks remain in place.
This is not full-history orthogonalization and does not guarantee convergence
with arbitrary nonlinear or history-dependent quantization.

The implementation stores one previous local residual per owned subdomain,
copies it once per update, and computes one additional local dot product on
later updates. The numerator occupies the otherwise unused second scalar in
the existing prolongation Allreduce. Thus it adds no reduction, MPI message,
payload bytes, matrix-vector product or restriction pass. Residual history is
cleared at each RHS restart. Quantizers, exact corners, float32 scale/minimum
metadata and post-CG constant correction are unchanged.

Rebuild in the container and run the focused validation:

```bash
(
set -e
python3 tests/test_bddc_flexible_cg.py
python3 tests/test_emi_flexible_cg_check.py
make -W EMI_model.cpp -o /kaskade7/libs/libkaskade.a \
  MPI=1 MPICXX=mpicxx \
  MPI_INCLUDE_FLAGS="-include /usr/lib64/mpi/gcc/openmpi4/include/mpi.h -DMUMPS_MPI_H"
python3 tests/run_emi_flexible_cg_check.py \
  --reference-campaign output/emi-mixed-validation-20261101T101620Z-6974 \
  --cases 10cells strong_h16
)
```

The driver preserves prior results in their existing directories. It runs four
configurations per problem: standard uncompressed CG, standard full R16/P16,
flexible full R16/P16, and flexible full zero-preserving-R12/original-P12.
All use the saved mesh/physics, ten steps, eight ranks, one restriction pass,
300-update limit and tolerance `1e-7`. Saved direct solutions are reused.
The rebuilt standard controls must reproduce saved convergence CSVs and logical
byte counts before flexible experiments proceed for that problem. A failed
flexible-16 run is recorded but does not prevent the flexible-12 test.

The summary reports bytes relative to both standard and flexible 16/16, iteration
counts and final full-precision DOF errors. Every timestep must converge; raw and
centered errors are retained separately with no hidden accuracy threshold.
Byte totals include the existing constant-correction broadcasts; other
diagnostic collectives and MPI framing remain outside these logical counters.
Use the saved phase timings only as validation timings, not performance claims.
The transfer policies still differ between restriction and prolongation even
when both alphabets are 12-bit. This experiment evaluates the established
configurations, not bit depth in isolation from quantizer choice.

The local C++ harness compiles production solver/subdomain methods and checks
a changing SPD preconditioner, fixed-PCG equivalence, adjacent conjugacy,
original-residual history, reset, negative beta, nonfinite rejection and exchange
counts. Driver tests check controls, parity gates, command isolation and failure
retention. These tests do not establish convergence on the EMI meshes; the
container run above is required. For the flexible-CG framework see
[Notay (2000), Flexible Conjugate Gradients](https://doi.org/10.1137/S1064827599362314).

Container checkpoint `output/emi-flexible-cg-20261101T105600Z-7824` completed
all four configurations on both problems, with ten steps, eight ranks, one
restriction pass and true tolerance `1e-7`. Standard-control history and byte
parity checks passed. The full-precision final relative L2 errors below are
against the saved direct solution; iteration counts and bytes cover all steps.

| Problem | Solver / restriction-prolongation bits | Iterations | Logical bytes | Final relative L2 |
| --- | --- | ---: | ---: | ---: |
| 10cells | Standard CG, 16/16 | 210 | 24,696,304 | 1.5588e-9 |
| 10cells | Flexible CG, 16/16 | 210 | 24,698,992 | 1.5923e-9 |
| 10cells | Flexible CG, 12/12 | 329 | 28,876,424 | 4.1044e-10 |
| strong_h16 | Standard CG, 16/16 | 396 | 134,795,008 | 2.5117e-9 |
| strong_h16 | Flexible CG, 16/16 | 396 | 135,192,416 | 2.5546e-9 |
| strong_h16 | Flexible CG, 12/12 | 717 | 199,114,888 | 1.5260e-9 |

This demonstrates single-pass 12/12 convergence for these two cases, but not
the communication objective: total bytes exceed standard 16/16 by 16.93% and
47.72%, respectively. Against flexible 16/16 the increases are 16.91% and
47.28%. The additional iterations outweigh the reduced payload per iteration.
These are validation results, not runtime or scaling claims. Keep the campaign
directories as research data; they are not part of the source-code commit.

To test flexible CG with **8-bit restriction and 8-bit prolongation**, reuse
the successful controls above with the same executable (no rebuild required):

```bash
python3 tests/run_emi_flexible_cg_check.py \
  --reference-campaign output/emi-mixed-validation-20261101T101620Z-6974 \
  --reuse-controls-from output/emi-flexible-cg-20261101T105600Z-7824 \
  --cases 10cells strong_h16 --bits 8
```

This launches only two new candidate runs. `--bits` sets both directional
precisions, defaults to 12, and does not change the solver, quantizer policies,
float32 metadata, one-pass schedule, saved timestep count or tolerance. Output
goes to a fresh directory. Reuse requires matching executable hash, reference
settings, passing controls and matching commands; reused rows are identified
explicitly. Their timings are historical, not new measurements. If the binary
has changed, omit `--reuse-controls-from` to validate the controls again.
An 8/8 failure is preserved without a higher-precision retry.

Further saved container checks, using the same controls and executable:

| Flexible CG bits | 10cells: first-step / total iterations | strong_h16: first-step / total iterations | Bytes saved vs flexible 16/16 (10cells / strong_h16) |
| --- | ---: | ---: | ---: |
| 16/16 | 21 / 210 | 40 / 396 | baseline |
| 14/14 | 26 / 225 | 42 / 421 | 7.21% / 4.59% |
| 13/13 | 36 / 283 | 68 / 562 | -9.05% / -21.75% |
| 12/12 | 36 / 329 | 83 / 717 | -16.91% / -47.28% |
| 8/8 | breakdown at call 7 of step 1 | breakdown at call 14 of step 1 | unavailable |

The additional campaigns are `emi-flexible-cg-20261101T111850Z-8325` (8/8),
`emi-flexible-cg-20261101T112121Z-8438` (14/14) and
`emi-flexible-cg-20261101T112507Z-8553` (13/13). Passing totals cover ten steps,
each meeting the true relative tolerance `1e-7`, with one restriction pass.
The first-step increase shows that the extra iteration cost at 12--13 bits
is already present before later timesteps. A missing parsed failure count is
not evidence of zero solver calls.

Final full-precision relative errors against direct for 14/14 are 8.2188e-10
and 2.1038e-9; for 13/13 they are 7.6280e-10 and 1.5068e-9, respectively.
Among the tested equal-bit candidates, 14/14 is the best for total logical
communication on these two problems. This does not establish runtime speedup,
global optimality or robustness on other meshes. Failed cases remain part of
the record.

Read-only diagnosis of the saved flexible 8/8 campaign identifies negative
`rho=r^T*z` as the rejection trigger, before the next beta/direction is formed:

| Problem | Rejected solver call | Accepted updates | Independently computed rho | True relative residual at failure |
| --- | ---: | ---: | ---: | ---: |
| 10cells | 7 | 6 | -2.6993465407952456e-5 | 0.434570 |
| strong_h16 | 14 | 13 | -3.9770263935861021e-5 | 0.540453 |

The malformed/interleaved MPI log caused the collector to report zero calls
for 10cells; its diagnostic CSV and original log both identify call 7. The
internal/global rho relative discrepancies at rejection are 7.54e-12 and
3.54e-14, respectively, so both calculations confirm the sign. Maximum recorded
adjacent-direction conjugacy defects on accepted updates are only 3.01e-12
and 7.40e-13. Shared correction copies agree to about 1e-16. Curvature is not
evaluated for the rejected directions; blank CSV entries must not be interpreted
as zero or negative curvature.

For the quadratic energy, the derivative along z is `-r^T*z`; these corrections
therefore increase energy for a sufficiently small positive step. Changing beta
or restarting alone does not repair this negative scalar. The evidence does
not isolate which transfer stage causes it. A controlled next diagnostic is
flexible R8/P16 and flexible R16/P8 against the failed R8/P8, retaining the
same quantizer policies, physics, one pass and tolerance. Such mixed-precision
runs would isolate sensitivity; they would not establish successful 8/8.
No solver or compression change was made for this diagnosis.

Run both directional isolation tests with the existing flexible-CG binary:

```bash
python3 tests/run_emi_flexible_cg_check.py \
  --reference-campaign output/emi-mixed-validation-20261101T101620Z-6974 \
  --reuse-controls-from output/emi-flexible-cg-20261101T105600Z-7824 \
  --cases 10cells strong_h16 \
  --directional-check
```

This launches four new runs: R8/P16 and R16/P8 on each problem. Both retain
zero-preserving restriction, the original prolongation quantizer, flexible CG,
one restriction pass, ten timesteps, eight ranks and tolerance `1e-7`. Relative
to failed R8/P8, only one directional precision changes in each candidate.
In particular, R16/P8 still uses **zero-preserving** restriction; it does not
switch to the original restriction quantizer used by the 16/16 controls.
Float32 scale/minimum metadata and post-CG correction are unchanged.

`--directional-check` and `--bits` are mutually exclusive. Controls are reused
only after the existing hash/settings checks; all new logs and commands go to a
fresh directory. Failure of one candidate does not prevent the other from
running. No compilation is needed because this extends only the Python driver.
Convergence would identify sensitivity to a transfer direction, not demonstrate
successful equal-bit 8/8 communication or prove a unique causal mechanism.

The directional campaign `emi-flexible-cg-20261101T113952Z-8669` passed all ten
steps for both mixed configurations. R8/P16 took 285 / 561 iterations and
24,935,216 / 155,751,168 logical bytes on 10cells / strong_h16. R16/P8 took
322 / 706 iterations and 28,336,640 / 198,009,200 bytes. None beats flexible
16/16 in total bytes. This supports investigating combined transfer error;
it does not prove that either quantizer alone is defective.

### Controlled exact-zero prolongation experiment

`--bddcProlongateZeroPreserving 1` is a new opt-in numerical experiment.
The original prolongation quantizer uses an offset grid, which need not contain
zero. The experimental mode applies the existing symmetric exact-zero grid to
the **weighted local correction values**, before optional integer graph lifting:
`Z=2^(b-1)`, `L=Z-1`, `q=Z+round(s*x)`, `x_hat=(q-Z)/s`, where the float32
scale is rounded down from `L/max(abs(x))`. For an all-zero buffer the scale is
one. The scale and reserved zero offset still occupy two float32 values.
No new per-value data, group descriptors, exchange, reduction or restriction
pass is added. Sender-local and receiver contributions use the same wire grid;
excluded corners remain exact. All ranks must enable the same mode.

This removes spurious corrections at exact zeros and makes the grid symmetric
under sign reversal. It does **not** guarantee positive `r^T*z` or convergence;
on asymmetric ranges it can coarsen the grid. It changes neither FCG nor the
restriction feedback. Encoded byte counts may increase or decrease because
the symbols change. Single-pass 8/8 success remains unverified until the PDE
experiment is run; the existing 14/14 results remain the validated choice.

Run in the container after syncing both repositories:

```bash
(
set -e
cd /BDDC/emiModelBddcSdc_new_scratch
mpicxx -std=c++20 -O2 -I/kaskade7 \
  tests/test_bddc_transfer_consistency.cpp -o /tmp/test_bddc_transfer_consistency
/tmp/test_bddc_transfer_consistency
make -W EMI_model.cpp -o /kaskade7/libs/libkaskade.a \
  MPI=1 MPICXX=mpicxx \
  MPI_INCLUDE_FLAGS="-include /usr/lib64/mpi/gcc/openmpi4/include/mpi.h -DMUMPS_MPI_H"
python3 tests/run_emi_flexible_cg_check.py \
  --reference-campaign output/emi-mixed-validation-20261101T101620Z-6974 \
  --cases 10cells strong_h16 \
  --prolongation-quantizer-check
)
```

This reruns uncompressed CG, standard 16/16, flexible 16/16, the original failed
flexible 8/8 and the new `flex-zero-r8-p8-zeroP` for each case. Only prolongation
quantization differs between the two 8/8 runs. Physics, ten steps, eight ranks,
`1e-7`, zero-preserving restriction, post-CG correction and one pass are fixed.
Standard control histories and bytes must reproduce the reference. Do not reuse
the old executable's controls after rebuilding. The expected failure of the old
8/8 run is retained and does not stop the new candidate; it makes the overall
driver exit nonzero even if the new candidate passes. Read each row separately.
All outputs go to a new directory. Assess both all-step convergence and total
logical bytes against flexible 16/16, with final full-precision error against
the saved direct solution. No precision fallback or relaxed tolerance is used.

Local verification covers the quantizer bound and zero/sign symmetry, option
forwarding, actual transfer/codec methods with minimal Dune container stubs,
and mocked driver continuation/parity checks. It is not an MPI/PDE validation.

### Existing steepest descent versus CG at 8/8

To test whether the existing BDDC steepest-descent path converges with 8-bit
compression, without changing any numerical implementation:

```bash
python3 tests/run_emi_sd_cg_comparison.py \
  --reference-campaign output/emi-mixed-validation-20261101T101620Z-6974 \
  --cases 10cells
```

This Python-only driver needs no rebuild of the current executable. It inherits
the saved problem (10cells with refinement 2, ten steps, eight ranks, unchanged
physical parameters), uses one restriction pass and tolerance `1e-7`, and writes
fresh commands, logs, metrics and full-precision comparisons. Six runs compare
uncompressed CG/SD, full-coded CG16/16, and full-coded CG8/8, FCG8/8 and SD8/8.
All three 8/8 modes use symmetric exact-zero restriction, original prolongation,
float32 metadata and the same coding options. The unsuccessful experimental
exact-zero prolongation mode is disabled. Failed 8-bit runs do not stop later
candidates. Failed uncompressed or 16-bit controls stop that case's experiments.

The current application supports post-solve constant correction and CG
diagnostics only with CG. Both are therefore disabled for **all** these runs,
including fresh controls; these results must not be directly substituted for
the older constant-corrected campaigns. Full-precision raw and centered final
coefficient errors versus the saved direct solution are both reported, without
a hidden accuracy threshold. The global true residual must pass at every step.
The SD implementation uses the restricted residual for its step numerator;
CG/FCG use the original-residual rho option. This tests the existing solver
paths and does not claim that only the beta recurrence differs.

The table compares **total logical bytes** against the fresh compressed CG16
baseline; an increased iteration count is acceptable if the total is smaller.
The default cap is 300 calls per step. `--iterations N` applies the same cap to
every run, without changing tolerance or precision. Calls on failed rows refer
to the failed step; missing failure counts are shown as unavailable, not zero.
Timings include setup and output and are retained in JSON, not used as a scaling
benchmark. The driver returns nonzero if any case fails, while keeping all data.
Use `--cases 10cells strong_h16` to include the second saved problem.

The saved `emi-sd-cg8-20261101T121744Z-9497` comparison passed uncompressed
CG/SD and full-coded CG16. All 8/8 modes failed in step 1: CG reached 300 calls
at relative residual 1.11486, FCG broke down at call 7 (0.43457), and SD broke
down at call 7 (1.06960). SD's step numerator was -4.25612e-5 while curvature
was positive (0.00613837). Its numerator uses the restricted residual, not an
independent true-residual dot product. No completed-solve byte saving is claimed.

### Read-only prolongation error isolation

The optional `--bddcProlongationDiagnostics 1` observes each proposed correction
before the solver accepts/rejects it, for CG, FCG or SD. The raw local/coarse
correction is already saved by BDDC. The observer reconstructs a separate
correction using **the actual transfer averaging weights** with unquantized
prolongation, including the interior Dirichlet/harmonic extension. It does not
replace the actual correction, call encoding routines, alter published residual
feedback, change alpha/beta, or advance the solution.

The resulting CSV has:

- `rho_before_p`: true global `(b-A*x)^T*z_shadow`, where only prolongation is exact.
- `rho_after_p`: the same product with the actual compressed correction.
- `rho_change_p`: `(b-A*x)^T*(z_actual-z_shadow)`, with a closure-error check.
- `rho_internal`: the solver's own numerator, which can differ from true rho in SD.
- Correction error norm/energy, restriction error, actual shared-copy jump and
  transfer-weight partition defect.

The shadow keeps the **actual quantized restriction and local/coarse-solve
results**. A negative shadow rho therefore means the problem already exists
without prolongation quantization at that iterate; it does not isolate an
implementation bug in restriction or prove how an exact-restriction trajectory
would evolve. Positive shadow rho and nonpositive actual rho identifies a sign
change caused by prolongation error at that iterate. Positive values in both
columns do not establish full convergence.

Replay the saved experiment after syncing and rebuilding both headers/application:

```bash
(
set -e
python3 tests/test_emi_prolongation_diagnostics.py --kaskade-root /kaskade7
make -W EMI_model.cpp -o /kaskade7/libs/libkaskade.a \
  MPI=1 MPICXX=mpicxx \
  MPI_INCLUDE_FLAGS="-include /usr/lib64/mpi/gcc/openmpi4/include/mpi.h -DMUMPS_MPI_H"
python3 tests/replay_emi_prolongation_diagnostics.py \
  --campaign output/emi-sd-cg8-20261101T121744Z-9497 \
  --case 10cells
)
```

The replay defaults to uncompressed CG, full CG16, full FCG8 and full SD8, each
with the observer off/on using the same rebuilt executable. Only the diagnostic
flag and fresh output directory change. It checks printed step/failure outcomes
against the saved run and compares all full-precision potential snapshots off/on.
This is outcome/snapshot parity, not a claim that SD exports every iteration.
Expected 8-bit solver failures are retained and are not failures of the replay.
The run stops if parity fails or the diagnostic is invalid. An uncompressed
relative shadow discrepancy above 1e-8 or weight-partition defect above 1e-10
invalidates interpretation; these are diagnostic sanity checks, not changed
solver tolerances. Every original result is preserved.

The observer adds two global sum reductions per observed proposal, local
Dirichlet solves, root matrix products, I/O, and output-status broadcasts.
These are diagnostic overhead, excluded from logical compression counters;
these runs must not be used for performance claims. Local tests exercise the
production shadow method/report on a known two-subdomain matrix, including a
positive-to-negative rho example, actual versus projection weights, inactive
DOFs and unchanged input vectors. A mocked replay verifies failure continuation
and detection of observer side effects. Full MPI/PDE validation remains required.


### Coarse-mean probe at the 8-bit failing correction

The replay `emi-prolongation-replay-20261101T123228Z-9906` preserved observer/source
parity. At step 1, call 7, FCG8 had true rho 5.86147e-6 with exact prolongation and
-2.69935e-5 with actual prolongation; SD8 changed from 5.97726e-5 to -4.69861e-5.
The corresponding relative correction errors were 0.003775 and 0.003878.
These are observed sign changes, not a proof that a replacement codec converges.

`--bddcProlongationConstraintCall 7` extends the observer at **step 1, call 7 only**.
It reads the actual scalar coarse constraint matrix rows (ones on disjoint
supports), verifies matching supports across owners, and rejects unsupported
weighted/overlapping rows. It does not use overlapping pairwise transfer lists
as coarse constraints. All numerical flags and float32 scale/minimum metadata
remain unchanged.

For each coarse group G it measures the interface error e between the actual
and exact-prolongation shadows, and its mean m_G. It lifts e restricted to G
through the existing interior Dirichlet solve, then measures r^T H(e_G).
Thus the group attribution includes interior effects, not only interface terms.
The sum is checked against the total rho change through an explicit
`unattributed_rho_error` column: a sizeable remainder requires investigation of
uncovered interfaces, shared-copy jumps, masking, or numerical cancellation.

The hypothetical corrected vector is z_actual - sum_G H(m_G 1_G). It removes
only coarse group means; within-group errors remain. A second shadow uses
float32-rounded means. Neither vector is passed to the solver. Positivity of
these shadow products is a necessary screening result, not convergence or
communication savings. The float32 test rounds the **globally known** mean;
it does not reproduce rounding of independently computed producer contributions.

Files `bddc-prolongation-constraints-{summary,groups,dofs}.csv` give the corrected
rho, remaining mean defects, individual group contributions, global DOF mapping,
and scalar-payload estimates. For a group with k owners the estimate is
k*(k-1)*4 bytes for one float32 contribution per directed owner delivery (8 bytes
for float64). It counts logical deliveries even within a rank and includes all
groups, even zero-error/exact corners. It assumes known topology and piggybacking;
headers, IDs, alignment, codebooks, possible additional messages/reductions and
any changed iteration count are excluded. This is a hypothetical payload model,
not measured traffic or an implemented one-pass correction protocol. A successful
future solver run must still beat matched CG16 total bytes.

After syncing, run inside the container:

```bash
(
set -e
cd /BDDC/emiModelBddcSdc_new_scratch
python3 tests/test_emi_prolongation_diagnostics.py --kaskade-root /kaskade7
python3 tests/test_emi_prolongation_replay.py
make -W EMI_model.cpp -o /kaskade7/libs/libkaskade.a \
  MPI=1 MPICXX=mpicxx \
  MPI_INCLUDE_FLAGS="-include /usr/lib64/mpi/gcc/openmpi4/include/mpi.h -DMUMPS_MPI_H"
python3 tests/replay_emi_prolongation_diagnostics.py \
  --campaign output/emi-sd-cg8-20261101T121744Z-9497 \
  --case 10cells --constraint-call 7
)
```

The off run disables both observers; the on run enables the prolongation observer
and this selected-call probe. Saved commands, PDE parameters, tolerance, ten
steps, iteration limits, 8/8 quantizers and solver updates are unchanged. Expected
8-bit failures remain expected. Missing target-call data or failed outcome/snapshot
parity invalidates the diagnostic. Fresh output directories preserve all earlier
results. The extra probe adds two reductions, root residual evaluation and two
harmonic solves per incident group plus one exact shadow solve per owned subdomain.
Do not use these heavily instrumented runs as timing benchmarks.

Local tests compile the production shadow/report with small known matrices;
they verify actual C-row extraction, complete interior rho attribution, restoration
of positive rho, float32 rounding, retention of within-group error, unsupported
constraint rejection and unchanged input state. Mocked replay tests check flag
isolation, target-call output validation and parity failure handling. A full
MPI application build and container replay are still required.

### Smaller prolongation groups: frozen-correction screening

The call-7 coarse-mean experiment (`emi-prolongation-replay-20261101T125121Z-10482`)
left FCG8 rho negative (-1.16919e-5), as well as SD8 (-2.14526e-5), despite
near-exact removal of group averages. The unexplained rho contributions were
about 1e-17. Float32 rounding of those mean corrections did not change this
conclusion. Mean preservation alone is therefore not implemented as a fix.

The optional `--bddcProlongationGroupingCall 7` instead screens quantization of
the **actual producer weighted buffers**, at step 1 call 7. It evaluates the
original affine quantizer with one grid per contiguous 8, 16, 32, 64, 128, 256 or
512 non-excluded values in producer-union order. These are quantization groups,
not coarse-constraint groups. Exact corners retain their original trailer values.
Each producer uses the same grid for its own contribution and all destinations;
quantization is not applied to an already-averaged global correction. Encoding
uses double parameters and decoding uses float32 scale/minimum, matching the
current affine quantizer. No exact-zero policy, clipping change, solver change,
additional restriction pass or higher-precision metadata is introduced.

The observer sums normalized producer contributions, performs the existing
interior harmonic extension on temporary vectors, and measures true r.z and the
error relative to the exact-prolongation shadow. All candidates see the same raw
correction from the original run; none changes the trajectory. Group size 0 at
the current bit depth must reproduce the actual correction within 1e-8 relative
(the same diagnostic sanity bound as the original shadow test). This checks the
producer model before any candidates are interpreted.

`bddc-prolongation-grouping.csv` records:

- 8-bit candidate grids, the original whole-buffer grid, and a whole-buffer
  16-bit shadow on the **same state and restriction result**.
- True rho, error norms, producer-group count and logical message count.
- Hypothetical packet bytes with the actual integer graph transform and current
  shared Huffman codebook. The table is read, never retrained or changed.
- Metadata, stream bytes, padding and exact-corner trailers separately.
- An alternative fixed-bit packed payload model, **not** the current quantization-
  only path, which still uses transmission-sized symbol slots.

For a grouped packet the proposed metadata model is uint32 run count followed by
(uint32 number of scalar symbols, float32 scale, float32 minimum) per contiguous
group run in packet order, padded to the transmission slot. This replaces the
legacy scale/minimum header and makes group membership self describing even
when a neighbor receives only part of a producer union. Huffman/tail stream-size
fields, separate stream byte rounding and final slot padding are counted.
Messages with encoding disabled by the existing threshold use transmission-sized
symbols in the model. Missing Huffman symbols invalidate the codec estimate;
they are never silently treated as a successful compressed payload.

These are **snapshot prolongation models**, not a new wire format implementation.
They exclude restriction traffic, MPI framing, training/reductions and effects of
changed iteration counts. The P16 shadow uses the current trajectory's codebook,
which need not be the one learned by a standalone 16-bit run. Beating this
same-state P16 model is only preliminary screening; total bytes must eventually
be compared with a fresh matched 16/16 solve.

After syncing both repositories, rebuild and replay:

```bash
(
set -e
cd /BDDC/emiModelBddcSdc_new_scratch
mpicxx -std=c++20 -O2 -I/kaskade7 \
  tests/test_bddc_transfer_consistency.cpp -o /tmp/test_bddc_transfer_consistency
/tmp/test_bddc_transfer_consistency
python3 tests/test_emi_prolongation_diagnostics.py --kaskade-root /kaskade7
make -W EMI_model.cpp -o /kaskade7/libs/libkaskade.a \
  MPI=1 MPICXX=mpicxx \
  MPI_INCLUDE_FLAGS="-include /usr/lib64/mpi/gcc/openmpi4/include/mpi.h -DMUMPS_MPI_H"
python3 tests/replay_emi_prolongation_diagnostics.py \
  --campaign output/emi-sd-cg8-20261101T121744Z-9497 \
  --case 10cells --grouping-call 7
)
```

Grouping replay defaults to CG16, FCG8 and SD8, each with diagnostics off/on.
All ten original steps, physics, tolerances and solver modes are retained. The
8-bit failures remain expected. The final intersection lists group sizes with
positive shadow rho and fewer modeled P bytes for **both** failing configurations;
an empty list means this screen found no shared candidate. Neither result proves
convergence. Old outputs remain untouched and observer/source outcome and
potential-snapshot parity are required. These heavily instrumented runs include
extra reductions and solves and must not be used as timing benchmarks.

Local validation compiles the production transfer routines against minimal Dune
stubs, compares whole-buffer reconstructions and packet byte counts, checks exact
corners/float32 singleton grids and overlapping three-owner interfaces, and
verifies unchanged transfer state/counters. The report harness checks normalized
producer assembly and harmonic extension on known matrices. Replay tests check
flag isolation, target rows, numerical parity and byte/convergence screening.
A full application build and MPI replay must still be performed in the container.

## Experimental one-step FGMRES comparison

`tests/run_emi_fgmres_check.py` tests whether right-preconditioned FGMRES(30)
can use the existing 8/8 codec when prolongation makes `r.z` negative. This
is a new optional outer solver, not a correction to the existing CG/FCG
recurrences. It does not guarantee 8-bit convergence or communication savings.
The solver stores each preconditioned Arnoldi vector, uses two-pass modified
Gram-Schmidt and Givens rotations, and accepts convergence only after explicitly
computing `||b-A*x||/||b|| <= 1e-7`. A small Hessenberg residual alone is insufficient.
The formulation follows flexible right GMRES (Saad, 1993); see also the
[PETSc FGMRES documentation](https://petsc.org/release/manualpages/KSP/KSPFGMRES/).

For this validation implementation, Arnoldi, the global matrix-vector products
and basis storage reside on rank zero. Each preconditioner application broadcasts
one normalized Arnoldi vector, applies distributed BDDC, and reduces weighted
local corrections to root. BDDC starts from zero on that RHS, including fresh
restriction feedback, while retaining factors, weights, codebooks and traffic
counters. Its output includes both the initial interior solve and the prolonged
correction. There is one restriction exchange per application and no CG update
or rho sign test. Thus this changes the outer algorithm and the lifetime of
feedback relative to the CG control; it is not merely a different beta formula.

The runner uses fresh direct, CG16, uncompressed FGMRES, FGMRES16 and FGMRES8
controls/candidate, in that order. All use the saved 10cells mesh and physics,
eight ranks for iterative solves, **one time step**, and tolerance `1e-7`.
FGMRES8 uses exact-zero restriction and original ungrouped affine prolongation.
Scale/minimum stay float32. No precision fallback, extra restriction pass,
constant correction or tolerance relaxation is applied. The CG16 first-step
residual summary must reproduce the saved source. A failed control stops the
experiment before the candidate. New output directories preserve earlier runs.

The report separates interface logical bytes from an outer-solver payload ledger.
The latter counts broadcasts/reductions as `(ranks-1)*payload` and allreduces
as `2*(ranks-1)*payload`, including command/status traffic. CG's existing global
solution assembly and true-residual broadcasts are counted on the same basis.
The sum is called **accounted logical bytes**, not total hardware network traffic:
coarse-grid traffic, codebook training, internal BDDC control reductions, MPI
framing and setup traffic remain excluded. Both interface-only and accounted
savings against fresh CG16 and FGMRES16 are saved. Positive interface-only savings
are insufficient if the added global-vector communication costs more.

Per-iteration true/estimated residuals are written to `bddc-fgmres-step1.csv`.
JSON retains full-precision raw/centered/max solution errors versus direct,
application/cycle/matvec counts, peak V/Z basis bytes on root, and validation
timings. Basis bytes exclude other work vectors and are not peak process RSS.
The explicit residual check adds a second global matvec per candidate. Root
orthogonalization and memory costs make this a validation experiment, not a
scalable production solver or timing benchmark. PASS means convergence and
finite comparisons; it does not impose an undisclosed solution-error threshold.

After syncing **both** repositories, run inside the container:

```bash
(
set -e
cd /BDDC/emiModelBddcSdc_new_scratch
python3 tests/test_bddc_preconditioner_apply.py --kaskade-root /kaskade7
mpicxx -std=c++20 -O2 -I. tests/test_emi_fgmres.cpp -o /tmp/test_emi_fgmres
/tmp/test_emi_fgmres
mpicxx -std=c++20 -O2 -DKASKADE_HAVE_MPI -I. \
  tests/test_emi_fgmres_mpi.cpp -o /tmp/test_emi_fgmres_mpi
mpirun --allow-run-as-root -np 2 /tmp/test_emi_fgmres_mpi
make -W EMI_model.cpp -o /kaskade7/libs/libkaskade.a \
  MPI=1 MPICXX=mpicxx \
  MPI_INCLUDE_FLAGS="-include /usr/lib64/mpi/gcc/openmpi4/include/mpi.h -DMUMPS_MPI_H"
python3 tests/run_emi_fgmres_check.py \
  --reference-campaign output/emi-sd-cg8-20261101T121744Z-9497
)
```

Local checks cover varying/negative preconditioners, restarting, true-residual
acceptance, breakdown and nonfinite handling; a production-method harness checks
the initial interior contribution and absence of CG updates. Driver tests check
fresh controls, immutable results and communication accounting. The MPI adapter
test exercises collective success and a nonfinite result on one rank. These
small tests do not replace the container build and EMI experiment above.

## Frozen first-input BDDC probe

`tests/run_emi_frozen_preconditioner_probe.py` compares uncompressed, R16/P16,
R8/P16, R16/P8 and R8/P8 on the **same first normalized RHS**. This is the first
Arnoldi input used by FGMRES. Each configuration starts in a separate MPI job,
uses one BDDC application with one restriction pass, writes the full-precision
`rhs`, `v`, assembled correction `z` and global matrix image `Az`, then exits
without advancing the physical state. It does not perform a convergence run or
print the usual completed-simulation marker. Normal solver behavior is unchanged
when `--bddcPreconditionerProbe` is zero (the default).

The driver reuses the saved CG16 command to retain mesh, physics, initial state,
eight ranks and tolerance 1e-7. Eight-bit restriction uses the existing exact-zero
quantizer; 16-bit restriction and all prolongations use the original quantizers.
No grouping, extra pass, precision fallback or constant correction is added.
Scale/minimum remain float32. Every job has independent fresh transfer feedback
and codebook state, so this is not a replay of later warmed-codebook applications.
The shared production BDDC application includes the initial interior solution.

`MEASURED` means one finite application and consistent input were verified, not
solver convergence. The driver rejects differing DOF counts, RHS values or input
vectors. Errors are `||z-z_ref||/||z_ref||` and
`||Az-Az_ref||/||Az_ref||`, with **uncompressed BDDC** as the reference; that
reference is not assumed to be the exact inverse. The angle between `v` and `Az`
and `sqrt(1-cosine^2)` describe the best residual obtainable using this one
direction. They are not a prediction of the full restarted solve. JSON also
records `v.z` and non-additivity of the two directional changes relative to
R16/P16. A large matrix-image error despite a small correction error motivates
investigating which components the quantizer damages; it does not establish a
specific cause or a successful fix by itself.

Run after syncing the updated EMI files (the existing library
`applyPreconditioner` API is reused):

```bash
(
set -e
cd /BDDC/emiModelBddcSdc_new_scratch
mpicxx -std=c++20 -O2 -DKASKADE_HAVE_MPI -I. \
  tests/test_emi_frozen_probe.cpp -o /tmp/test_emi_frozen_probe
mpirun --allow-run-as-root -np 2 /tmp/test_emi_frozen_probe /tmp/emi-frozen-probe-test.csv
make -W EMI_model.cpp -o /kaskade7/libs/libkaskade.a \
  MPI=1 MPICXX=mpicxx \
  MPI_INCLUDE_FLAGS="-include /usr/lib64/mpi/gcc/openmpi4/include/mpi.h -DMUMPS_MPI_H"
python3 tests/run_emi_frozen_preconditioner_probe.py \
  --reference-campaign output/emi-fgmres-20261101T133918Z-12380
)
```

New timestamped directories preserve all earlier campaigns. These five
applications are diagnostic work, not production communication or timing
measurements. The C++ adapter test checks normalization, assembly, `Az`, a
nonfinite worker result and a root-only output failure (the latter prints an
expected error before PASS). The Python tests check fixed-input comparisons,
directional settings, rejection of invalid data, and preservation of prior runs.

### Constant-split prolongation: controlled experiments

`--bddcProlongateConstantSplit 1` is an optional numerical representation, off by
default and permitted by EMI in frozen-probe mode or explicitly one-step FGMRES. For each
producer, it takes the arithmetic mean of the **unweighted, non-exact** interface
correction, rounds that mean to float32 as `c`, and quantizes `w*(u-c)` instead
of `w*u`. Receivers reconstruct `decoded + w*c` before accumulating the producer
contribution. The encoder subtracts the same rounded constant, so rounding `c`
does not silently discard its rounding remainder. This is neither a change to
the CG/FGMRES recurrence nor a post-solve mean correction. It changes the
quantization error and may change convergence; no improvement is guaranteed.

Each producer's averaging weights are cached from the existing setup exchange.
No new solver exchange or per-DOF weight payload is introduced. The cache adds
memory proportional to neighbor interface sizes. Computing the mean and restoring
the constant add linear work in the interface size. Each prolongation packet
adds an **8-byte prefix: uint32 `CSP1` format tag and float32 constant**. The
existing scale and minimum remain float32; exact corner trailers remain original
weighted values. The extra header is included in sent/received byte counters.
All participants must enable the same mode before weight exchange. This first
implementation requires scalar double transmission slots, active interfaces,
ungrouped affine prolongation, no zero-preserving prolongation and no legacy
ASCII dumps. Integer lifting, bit-length coding and Huffman still operate on
the quantized fluctuation codes. Grouped shadow diagnostics reject this mode.

The `--constant-split` probe runs fresh copies of the five existing configurations
and then one `r8-p8-split` candidate. The five controls must reproduce the saved
probe's RHS, normalized input, correction and matrix image exactly. The candidate
still uses R8 exact-zero restriction, P8 affine quantization, one restriction
pass and unchanged physical parameters. No physical time step is advanced.
The report measures actual interface packet bytes for this single application,
including the new header, metadata, padding, exact trailers and intra-rank
deliveries. It excludes coarse traffic, setup/codebooks, diagnostic reductions
and MPI framing. Single-application savings are not total-solve savings.

Sync both repositories and run in the container:

```bash
(
set -e
cd /BDDC/emiModelBddcSdc_new_scratch
mpicxx -std=c++20 -O2 -I/kaskade7 \
  tests/test_bddc_constant_split.cpp -o /tmp/test_bddc_constant_split
/tmp/test_bddc_constant_split
mpicxx -std=c++20 -O2 -I/kaskade7 \
  tests/test_bddc_transfer_consistency.cpp -o /tmp/test_bddc_transfer_consistency
/tmp/test_bddc_transfer_consistency
make -W EMI_model.cpp -o /kaskade7/libs/libkaskade.a \
  MPI=1 MPICXX=mpicxx \
  MPI_INCLUDE_FLAGS="-include /usr/lib64/mpi/gcc/openmpi4/include/mpi.h -DMUMPS_MPI_H"
python3 tests/run_emi_frozen_preconditioner_probe.py \
  --reference-campaign output/emi-fgmres-20261101T133918Z-12380 \
  --reference-probe output/emi-frozen-probe-20261101T134855Z-12782 \
  --constant-split
)
```

The new transfer test uses large non-float32-exact offsets, nonuniform weights,
three-owner interfaces, different local numbering, plain/full codec equivalence,
mixed/all/no exact trailers, malformed headers, missing setup weights and actual
packet-size accounting. Existing transfer tests cover the default modes. Local
transfer checks use minimal Dune stubs; they do not replace the container build
and MPI application probe. No long convergence run is launched automatically.

### One-step FGMRES8/8 with constant-split prolongation

After a successful frozen probe, `tests/run_emi_fgmres_check.py --constant-split`
runs fresh direct, CG16, uncompressed FGMRES and FGMRES16 controls followed by
`full-fgmres8-split`. It replaces the original failed FGMRES8 candidate in this
campaign; it does not rerun that 300-application test. The original runner mode
remains available without the flag. Both probe and split options are explicitly
disabled in every control. The split candidate is an actual time-step solve,
not a frozen probe.

The source must be a saved one-step FGMRES campaign with full-precision potential
snapshots. Each fresh control must reproduce its saved initial/final vectors and
residual summary; FGMRES controls must also reproduce their residual CSVs. A
changed or failing control stops the campaign before the candidate. The model
allows split prolongation only with FGMRES and either a frozen probe or an
explicit `maximumNumberOfTimeSteps` between 1 and 10. The runner defaults to
one step; `--steps 10 --constant-split` selects the longer validation below.
The runner retains restart 30, cap 300, tolerance `1e-7`, one restriction pass,
eight ranks, the same mesh/physics, and no post-solve constant correction.

The split codec and float32 metadata are unchanged from the frozen probe. Its
eight-byte headers are already included in interface byte counters; they must
not be added again. The report gives interface-only and accounted
(interface plus outer collective payload model) savings against both CG16 and
FGMRES16. The outer model has the same exclusions documented above, so it is
not total hardware network traffic. Raw and centered relative L2 errors, raw
maximum error and mean error against fresh direct are displayed. No mean is
removed from the actual solution. PASS checks convergence and finite errors;
there is no implicit physical-accuracy acceptance threshold. Failed runs retain
their residual CSV and print the best/last true residual without claiming byte
savings. Timings and root basis memory remain diagnostic measurements in JSON.

After syncing the updated EMI files, rebuild and run:

```bash
(
set -e
cd /BDDC/emiModelBddcSdc_new_scratch
make -W EMI_model.cpp -o /kaskade7/libs/libkaskade.a \
  MPI=1 MPICXX=mpicxx \
  MPI_INCLUDE_FLAGS="-include /usr/lib64/mpi/gcc/openmpi4/include/mpi.h -DMUMPS_MPI_H"
python3 tests/run_emi_fgmres_check.py \
  --reference-campaign output/emi-fgmres-20261101T133918Z-12380 \
  --constant-split
)
```

Runner tests cover unchanged physical inputs, registered command-line options,
split selection, fresh-control parity failures, partial failed residual traces,
missing codec markers, and cases where interface savings are outweighed by
outer communication. They do not replace the container solve above.

### Ten-step validation of the successful constant-split FGMRES8/8 case

Use the successful one-step split campaign as the reference. The ten-step runner
requires completed saved direct, CG16, uncompressed FGMRES, FGMRES16 and split8
runs. It reruns all five configurations with the same first-step inputs,
`dt`, physics, FGMRES(30), 300-application cap per step, one restriction pass,
and tolerance `1e-7`. Only final time and maximum step count increase. No codec,
solver recurrence, float32 metadata, or post-solve correction changes are made.
Codebooks/factors retain the application's existing lifetime; each BDDC
preconditioner application still starts with fresh RHS feedback.

Every new run must reproduce the saved first-step potential and residual;
FGMRES residual traces must also match exactly, including for split8. A control
failure stops the campaign before the candidate. Every one of the ten steps must
meet the true tolerance, and all ten accepted potential snapshots must exist.
`per-step-errors.csv` reports raw and centered relative L2, maximum absolute and
mean errors against the fresh direct solution at the same time, together with
iterations and true residual. Centering is diagnostic only. Each trajectory
advances its own state; it is not reset to direct between steps.

The final report compares cumulative interface and accounted bytes against
both fresh 16/16 controls. The outer ledger sums applications and overhead for
**each** step, rather than treating ten steps as one large solve. It includes
all constant-split headers in interface bytes exactly once. Failed runs never
receive a communication-saving claim. Solver-wrapper and application times are
retained as validation timings, not performance measurements. Previous campaigns
are preserved in place; the runner creates a new output directory.

After syncing the updated EMI files, run inside the container:

```bash
(
set -e
cd /BDDC/emiModelBddcSdc_new_scratch
make -W EMI_model.cpp -o /kaskade7/libs/libkaskade.a \
  MPI=1 MPICXX=mpicxx \
  MPI_INCLUDE_FLAGS="-include /usr/lib64/mpi/gcc/openmpi4/include/mpi.h -DMUMPS_MPI_H"
python3 tests/run_emi_fgmres_check.py \
  --reference-campaign output/emi-fgmres-20261101T145024Z-13570 \
  --constant-split --steps 10
)
```

This stage tests only 10cells. Weak/strong cases are not launched automatically.
The first-step result established convergence and lower interface bytes, but
not lower accounted total communication; the ten-step campaign must report both.

### Checkpoint: ten-step 10cells split8 result (2026-11-01)

Container campaign `output/emi-fgmres-20261101T145954Z-13879` passed all
ten steps at true relative tolerance `1e-7`, with one restriction pass and
unchanged float32 constant/scale/minimum metadata. First-step parity against
`emi-fgmres-20261101T145024Z-13570` passed for all five fresh runs.

| Mode | Applications (10 steps) | Interface bytes | Accounted bytes | Final raw relative L2 | Final centered relative L2 |
| --- | ---: | ---: | ---: | ---: | ---: |
| CG16/16 | 210 | 24,696,992 | 932,952,032 | 1.878023e-7 | 1.539172e-9 |
| Uncompressed FGMRES | 209 | 90,782,912 | 974,266,268 | 1.977538e-9 | 1.973368e-9 |
| FGMRES16/16 | 210 | 24,688,048 | 912,299,808 | 3.135092e-8 | 1.988078e-9 |
| FGMRES8/8 split | 227 | 13,686,424 | 971,481,052 | 2.335418e-5 | 1.393521e-9 |

Split8 applications per step were `[22,22,23,23,23,23,23,23,23,22]`;
all true residuals were at most `9.092860e-8`. Compared with FGMRES16,
interface bytes fell 44.56%, while accounted bytes rose 6.49%. The overall
communication objective is **not met**. Accounted bytes use the exclusions
above and are not total measured network traffic. The candidate's final
mean error was `1.956782e-6`, maximum absolute error `1.957663e-6`;
centering is only a diagnostic and does not change the accepted solution.
Solver-wrapper times were 45.217 s (CG16), 42.254 s (uncompressed FGMRES),
46.060 s (FGMRES16), and 47.955 s (split8). These single validation runs
include diagnostic/output costs and establish no timing speedup. Raw output
files remain outside source control; retain the named campaign separately.

### Continue split8 validation on the existing weak/strong meshes

After the ten-step 10cells campaign succeeds, use `--case weak_8_h16` and
`--case strong_h16` with `--constant-split --steps 10`. The runner validates
all saved 10cells control/candidate residual summaries and FGMRES traces as
preflight evidence, then runs five **fresh** configurations for each new mesh.
There is no residual/solution parity claim across different meshes. The default
10cells workflow retains its first-step parity checks.

Inputs come from the existing
`../input_emi_mesh_blocks/weak_strong_scaling/{case}_kaskade` files, using
`refine=0`, as in the earlier weak/strong campaigns. `--mesh-dir` can select
another location for those same files; this runner does not generate meshes.
The source physical parameters remain sigma_i=3, sigma_e=20, C_m=1, R=0.1,
dt=0.01. Restart 30, 300 applications per step, eight ranks, true tolerance
1e-7, one restriction pass, and the constant-split codec remain fixed.

Only the Python runner changes for this extension; an executable already built
for the successful ten-step 10cells campaign needs no rebuild. After syncing
the updated runner, use:

```bash
(
set -e
cd /BDDC/emiModelBddcSdc_new_scratch
for case_name in weak_8_h16 strong_h16; do
  python3 tests/run_emi_fgmres_check.py \
    --reference-campaign output/emi-fgmres-20261101T145954Z-13879 \
    --constant-split --steps 10 --case "$case_name"
done
)
```

The loop stops if a case fails; all completed and partial logs remain in its
new output directory. Each case has independent fresh direct/CG16/FGMRES16
references and its own summary and per-step error CSV. As before, a failing
control prevents running the candidate, and a failing candidate has no claimed
communication savings. No weak/strong success is assumed from 10cells.

### One-step 6/6 and 4/4 screening with the successful constant split

`--bits 6 4 --constant-split --steps 1 --case 10cells` changes only the common,
restriction and prolongation bit depths for the two new candidates. It retains
FGMRES(30), a 300-application limit, one restriction pass, true tolerance 1e-7,
exact-zero restriction, ungrouped affine constant-split prolongation, full
coding, float32 metadata and existing exact-corner handling. It adds no
precision fallback, refinement exchange or post-solve correction.

The runner first executes fresh direct, CG16, uncompressed FGMRES, FGMRES16
and split8 controls, requiring saved first-step parity for all five. A failed
control stops screening. Then split6 and split4 run independently; failure of
split6 does not prevent split4. A nonzero final exit status means at least one
run failed. Failures retain partial residual traces and have no claimed savings.
The report compares both interface and accounted bytes with uncompressed
FGMRES, CG16, FGMRES16 and split8, so the baseline and scope remain explicit.

After syncing the Python runner (no rebuild needed with the validated split8
executable), run:

```bash
cd /BDDC/emiModelBddcSdc_new_scratch
python3 tests/run_emi_fgmres_check.py \
  --reference-campaign output/emi-fgmres-20261101T145954Z-13879 \
  --constant-split --steps 1 --case 10cells --bits 6 4
```

Four-bit screening remains restricted to one step on 10cells. Six-bit ten-step
validation is available below and requires a successful one-step reference.

### Precision study: extended 4/4 cap and ten-step 6/6 validation

The saved 4/4 trace in `emi-fgmres-20261101T200034Z-14562` decreases from
1.0722646e-5 at application 200 to 1.4172615e-7 at application 300, with no
increases in that interval. Short plateaus occur within restart cycles; this
is evidence for a cap-only extension, not demonstrated convergence.

Sync the updated `tests/run_emi_fgmres_check.py` and the new
`tests/run_emi_precision_followup.py`, then run in the existing Kaskade container:

```bash
cd /BDDC/emiModelBddcSdc_new_scratch
python3 tests/run_emi_precision_followup.py \
  --reference-campaign output/emi-fgmres-20261101T145954Z-13879 \
  --screening-campaign output/emi-fgmres-20261101T200034Z-14562
```

This creates a new directory, runs 4/4 from the original initial state with cap
500, and requires exact CSV-prefix parity with the original 300-iteration run.
It runs cap 1,000 only if fresh controls pass, the 500-cap trace reproduces the
saved prefix, the failure is an iteration-limit exit, and the last 100 intervals
decrease beyond relative 1e-8 numerical noise with the final residual the best
in that interval (also allowing relative 1e-8 noise). This gate only schedules
another experiment; the solver tolerance remains
1e-7. Every extended run starts afresh, not from a checkpoint. Controls keep
their original 300 cap. A 1,000 run must reproduce the full 500 trace.

Six-bit ten-step validation then runs independently with fresh direct,
uncompressed FGMRES, CG16, FGMRES16 and split8 controls. The first six-bit step
must reproduce the successful one-step solution/residual trace. Every accepted
step is checked against 1e-7 and a fresh direct solution. No rebuild is needed
for these Python-only changes with the existing validated executable. The
numerical algorithms, physics, float32 metadata, restart=30, one restriction
pass and no-fallback policy are unchanged. Add `--dry-run` to save the plan and
print commands without executing the solver.

Standalone equivalents, if needed:

```bash
python3 tests/run_emi_fgmres_check.py \
  --reference-campaign output/emi-fgmres-20261101T145954Z-13879 \
  --screening-campaign output/emi-fgmres-20261101T200034Z-14562 \
  --constant-split --case 10cells --steps 1 --bits 4 --four-bit-limit 500

python3 tests/run_emi_fgmres_check.py \
  --reference-campaign output/emi-fgmres-20261101T145954Z-13879 \
  --screening-campaign output/emi-fgmres-20261101T200034Z-14562 \
  --constant-split --case 10cells --steps 10 --bits 6
```

`tests/analyze_emi_fgmres_precision.py` reads campaigns without altering them
and writes a new HTML report, CSV tables, SHA-256 input manifest, and PDF/SVG/PNG
figures. It compares accepted snapshots with both direct and uncompressed
FGMRES references at every step. It distinguishes interface bytes, accounted
outer-plus-interface payload, unavailable network/communication-time data,
failed runs, and one-step versus ten-step evidence. Requires NumPy and Matplotlib:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  python3 tests/analyze_emi_fgmres_precision.py
```

The default inputs include the November 1 uniform and mixed one-step screens,
the successful 4/4 cap extension, and three ten-step campaigns (including 6/6
for 10cells). To analyze newly completed runs, specify each campaign
with repeatable `--campaign PATH` arguments (these replace the default list).
The report does not estimate unmeasured codec error or MPI time, or certify
physical accuracy from centered coefficient errors. Mixed-direction experiments
are available as a subsequent one-step stage below.

### Follow-up results and controlled mixed-direction screen

The campaign `output/emi-precision-followup-20261101T202243Z-14983` establishes:

- `b4-cap500`: 4/4 converges at application 301, true relative residual
  5.778563937e-8. All 301 original observations (iterations 0 through 300)
  reproduce exactly. Only `--bddcIterations` (300 to 500) and the output
  directory changed. No 1,000-cap run was needed. This is one-step evidence.
- `b6-ten-steps`: 6/6 passes all ten steps with application counts
  `[31,32,32,32,32,33,32,33,32,33]`, totaling 322. Its 14,814,944 interface
  bytes exceed matched split8's 13,686,424 by 8.25%. It saves 83.68% versus
  uncompressed FGMRES's interface total, but 8/8 remains more efficient.
- The completed 4/4 run sends 9,570,528 interface bytes, 10.17% more than
  uncompressed FGMRES for one step. Convergence and communication efficiency
  are distinct outcomes.

The next experiment changes only **one directional precision argument** from
the split8 control. With the updated Python driver (no numerical-code rebuild):

```bash
cd /BDDC/emiModelBddcSdc_new_scratch
python3 tests/run_emi_fgmres_check.py \
  --reference-campaign output/emi-precision-followup-20261101T202243Z-14983/b6-ten-steps \
  --constant-split --case 10cells --steps 1 \
  --pairs 8/6 6/8
```

The runner executes fresh direct, CG16, uncompressed FGMRES, FGMRES16, and
split8 controls before R8/P6 and R6/P8. Controls must reproduce saved first-step
solutions and convergence traces. Mixed modes preserve the generic precision
default of 8; the explicit directional override changes only R or P. All
zero/split policies, full-codec stages, restart=30, cap=300, tolerance=1e-7,
float32 metadata and one restriction pass are identical to split8. Candidate
failures are retained and do not prevent the other candidate; control failures
stop the campaign. `--pairs` cannot be combined with `--bits` or a changed four-bit cap.
Ten-step mixed validation is restricted to screened R6/P8 as described below. R8/P4 and R4/P8 are supported for later investigation,
but are not part of the command above.

The precision-analysis script now includes the successful 4/4 extension and
6/6 ten-step results by default, checks the 4/4 prefix independently, and keeps
the original failed 300-cap campaign as historical evidence. Reported solution
errors compare both direct and uncompressed references. The three ten-step
cases have separate figures for all FGMRES residual histories.

### R6/P8 one-step result and ten-step validation

`emi-fgmres-20261101T203933Z-15690` passed both mixed one-step candidates.
R6/P8 kept 22 applications and used 1,150,088 interface bytes, saving 12.77%
versus split8 and 86.76% versus uncompressed FGMRES. R8/P6 needed 32
applications and 1,708,232 bytes, 29.56% more than split8. This directional
comparison applies to this input and these quantizer policies.
R6/P8 raw/centered errors versus direct were 3.996068e-4 / 3.748169e-9;
its larger mean offset is retained, not treated as a harmless gauge shift.
Its accounted total still exceeds FGMRES16 by 3.13%.

Run the next validation with the updated Python driver; no rebuild is needed:

```bash
python3 tests/run_emi_fgmres_check.py \
  --reference-campaign output/emi-precision-followup-20261101T202243Z-14983/b6-ten-steps \
  --screening-campaign output/emi-fgmres-20261101T203933Z-15690 \
  --constant-split --case 10cells --steps 10 \
  --pairs 6/8
```

The driver reruns controls and R6/P8, checks exact first-step trace/solution
parity against the successful screen, and requires every subsequent step to
meet 1e-7. One pass, restart 30, cap 300, float32 metadata and numerical
implementation are unchanged. Results go to a fresh directory. Ten-step
R6/P8 results are recorded below.

### R6/P8 validation on weak and strong meshes

`emi-fgmres-20261101T205726Z-16111` passed all ten 10cells steps for
R6/P8 with 229 applications (22, then nine steps of 23). The candidate sent
12,028,992 full-compressed interface bytes: 86.75% less than uncompressed
FGMRES, 51.28% less than compressed FGMRES16, and 12.11% less than split8.
Accounted interface-plus-outer bytes were 978,080,428, 0.68% more than split8.
Final direct-reference raw/centered relative errors were 5.289474e-5 and
1.601946e-9. Convergence does not establish physical accuracy or network savings.

The next stage keeps R6/P8, full compression, constant splitting, float32
metadata, restart 30, cap 300, one restriction pass and tolerance 1e-7 fixed.
Only the established mesh inputs/refinement change between cases. The runner
requires a completed ten-step 10cells R6/P8 reference and checks its settings,
all accepted residuals and all FGMRES traces before creating new results.
It records reference command/log hashes. No first-step equality is imposed
between different meshes. Fresh direct, CG16, uncompressed FGMRES, FGMRES16
and split8 controls precede the mixed candidate for each mesh.

After syncing the Python driver (no numerical rebuild), run:

```bash
cd /BDDC/emiModelBddcSdc_new_scratch
for case_name in weak_8_h16 strong_h16; do
  python3 tests/run_emi_fgmres_check.py \
    --reference-campaign output/emi-fgmres-20261101T205726Z-16111 \
    --constant-split --steps 10 --case "$case_name" \
    --pairs 6/8 || printf 'Case %s failed; inspect its retained logs.\n' "$case_name"
done
```

Each case writes a unique output directory. All ten steps must converge;
failed candidates are not assigned communication savings. A failed control
stops that case before the candidate. The shell loop still attempts the other
case. Prior results are untouched. The completed Linux validation results
are recorded below.

### Completed R6/P8 ten-step checkpoint

All 30 time steps passed the true relative residual tolerance 1e-7 across
10cells, weak_8_h16 and strong_h16. Each used one restriction pass, FGMRES(30),
full compression, constant-split prolongation and unchanged float32 metadata.
These PDE campaigns were run by the user in the Linux Kaskade container.

| Case | Campaign under output/ | Applications | Interface bytes | Interface saving vs uncompressed FGMRES | Interface saving vs split8 |
|---|---|---:|---:|---:|---:|
| 10cells | emi-fgmres-20261101T205726Z-16111 | 229 | 12,028,992 | 86.75% | 12.11% |
| weak_8_h16 | emi-fgmres-20261101T210711Z-16464 | 373 | 7,023,040 | 78.93% | 6.94% |
| strong_h16 | emi-fgmres-20261101T211157Z-16694 | 470 | 101,559,440 | 73.63% | 5.80% |

The totals cover ten steps, including all extra applications and packet metadata.
R6/P8 is mixed precision, not uniform six-bit communication. Accounted
interface-plus-outer totals were 978,080,428 / 141,010,796 / 1,364,359,440
bytes, respectively: 0.68% / 1.41% / 0.91% MORE than matched split8. Thus
interface savings do not imply lower accounted total or measured network traffic.

Final raw / mean-centered relative coefficient errors against direct:

- 10cells: 5.289474e-5 / 1.601946e-9.
- weak_8_h16: 1.603363e-4 / 3.351779e-9.
- strong_h16: 1.950395e-4 / 2.803261e-9.

Mean offsets dominate the raw errors; centering remains diagnostic only. These
results establish convergence on the tested cases, not physical-accuracy
certification or general convergence guarantees. Runtime measurements include
validation/output and are not replicated performance benchmarks. The analysis
script's default campaign selection retains the earlier precision-screening
cohort; its generated reports must not be mistaken for this newer checkpoint.

Local checkpoint validation: 22 experiment-driver tests and the follow-up
cap-decision test pass. The commit contains drivers, tests, analysis tooling
and this documentation; generated campaigns, binaries and logs are not included.
