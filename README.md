# EMI model with Kaskade7

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
