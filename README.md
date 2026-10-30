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
