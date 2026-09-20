# EMI, SDC, AA, and BDDC: Numerical Notes

This note describes the implementation in this repository. It records the
operator definitions and algorithmic assumptions that matter when changing
the discretization or comparing solver configurations. The main implementation
is in `EMI_model.hh`, `EMI_model.cpp`, and `EMI_bddc.hh`.

## EMI model and mesh contract

The unknown is a scalar potential `u` represented on a piecewise-continuous
space. Extracellular material tags share one potential component. Every
intracellular material tag receives its own disconnected component. This lets
the model represent intracellular/extracellular membrane interfaces while
keeping a single scalar unknown in the implementation.

The weak model includes:

- Bulk conductivity terms, using `sigma_i` or `sigma_e` according to the cell tag.
- Membrane capacitance `C_m` across intracellular/extracellular interfaces.
- The configured membrane current on those interfaces.
- Gap-junction conductance `R` between distinct intracellular tags.
- A penalty term on extracellular boundary faces, controlled by `penalty`.

The mesh is read from VTU cell data named `domain`. Every cell tag should be
classified in exactly one of the extracellular or intracellular lists. The
three text files (`extra_set`, `intra_set`, `excited`) each contain a count on
the first line followed by that many integer tags. The excited tags should be
intracellular tags. The initial potential is `0.5` on excited cells and `0`
elsewhere. `R_extra` is accepted and stored, but the current weak-form
assembly does not use it. Do not interpret changing this option as changing
the numerical model until that term is implemented.

## Time stepping and operator convention

The time-step count is

\[
N = \begin{cases}
\min(N_{\max},\lceil T/dt\rceil), & N_{\max}>0,\\
\lceil T/dt\rceil, & N_{\max}=0.
\end{cases}
\]

The final interval is clipped at `finalTime`. `maximumNumberOfTimeSteps` is a
cap, not a replacement for the requested final time. The scenario makefile
currently sets this cap to one for quick checks.

For SDC, the model assembles two operators separately:

- `M`: membrane-capacitance mass matrix (`tau = 0`, `Mass_stiff(1)`).
- `A`: spatial-operator Jacobian (`tau = 1`, `Mass_stiff(0)`).

Here `A` means the matrix produced by this implementation's signed weak-form
and `SemiImplicitEuler` conventions; do not reinterpret it as a positive
stiffness matrix without applying the corresponding sign change. The
collocation interval `i` uses

\[
J_i = M - \widehat S_{i-1,i} A,
\]

where `Shat` is the preconditioner/integration matrix selected by
`sdcSweepType` (Euler or LU). The interval RHS combines the mass difference
`M (u_{i-1} - u_i)`, the quadrature of residuals at the collocation nodes, and
terms involving corrections already computed in the current sweep. The
implementation's exact assembly is in `sdcIterationStepJacobi` and
`sdcIterationStepBddc`.

The Radau grid starts with `sdcStartCollocationPoints` collocation points
after the initial time point (the state vector therefore has one additional
entry). It can be refined between sweeps up to `sdcCollocationPoints`. Each SDC step runs at
least `minimumSdcSweeps` and at most `maximumSdcSweeps`. After the minimum, it
stops when either the correction norm is below `sdcTolerance`, or the estimated
remaining correction is below `sdcAbsoluteTolerance` with a contraction
estimate less than one. The contraction estimate is updated from consecutive
sweep correction norms and initialized by `sdcInitialContraction`.

`Mass_stiff` controls which Jacobian contributions the assembler sees:

| Value | Assembled Jacobian contributions |
|---|---|
| `1` | Membrane mass/capacitance only |
| `0` | Spatial operator only |
| greater than `1` | Combined implicit-Euler operator |

## Algebraic adaptivity

AA is applied between SDC sweeps. For an active DOF with largest correction
`duMax` over collocation nodes and estimated contraction `q`, the selection
estimate is `q * duMax / (1 - q)`. A DOF is retained when this exceeds
`algebraicAdaptivityTolerance`; if `q >= 1`, all current candidates are
retained. Selected DOFs are expanded by the union of the mass and spatial
operator matrix neighborhoods to avoid dropping directly coupled unknowns.

In ordinary SDC, the active set is a global DOF list. Its matrices are reduced
to active rows and columns using a global-to-compact index map. Residual
assembly uses the inverse DOF-to-cell map and includes every cell touching an
active residual row; filtering only to cells whose own DOFs are all active
would omit element contributions to active rows.

In SDC+BDDC, selection is converted to complete owner subdomains. If any
selected DOF or its matrix neighbor belongs to a subdomain, that whole
subdomain is activated for the next BDDC sweep. This preserves the
subdomain-local solve and interface-constraint structure. The same DOF-to-cell
map limits residual assembly to cells incident to the active subdomains.

The first-sweep residual reuse assumes the current EMI functional is
autonomous: initial collocation guesses are identical at every node, so their
residuals match. If explicit time-dependent forcing is added, this reuse must
be revisited.

## BDDC partition and local systems

Each distinct material tag defines one subdomain in the current driver; all
cells with the same tag are grouped together. Global DOFs on material
interfaces have local copies in multiple subdomains. Their multiplicity
weights sum to one and are used when scattering global residuals and combining
local corrections. The code stores the global-DOF-to-owner and
global-DOF-to-cell inverse maps for AA and residual assembly.

For BDDC without SDC, each subdomain receives the restricted combined
implicit-Euler matrix. For SDC+BDDC, the local mass and spatial operators are
extracted separately and reused; the interval-dependent local matrix is
`M_local - Shat[i-1][i] * A_local`. This is why SDC+BDDC does not need the
separate combined Euler matrix.

`bddcInterfaceTypes` is a bit mask: corner `1`, edge `2`, face `4`, and `7`
selects all three. `bddcIterations` and `bddcTolerance` control each BDDC
linear solve. `bddcVerbose` enables per-iteration residual output.

The `nThreads` option controls assembler and per-subdomain setup/extraction
work. The default is two; set it to the number of CPUs allocated to the job,
not necessarily the machine's total core count. The MPI communication layer is
not implemented in this EMI driver. Future MPI work is intended to modify the
original Kaskade BDDC communication path in `mg/bddc.hpp`.

## Compression

With `bddcCompression=0`, the standard `SpaceTransfer` is used. With
`bddcCompression=1`, the compressed transfer path configures DCT transformation,
quantization, bit-length encoding, and Huffman coding for restriction and
prolongation. `bddcCompressionBits` controls quantization precision. This
option concerns BDDC interface-transfer data; it does not compress the global
matrix or the VTK solution fields. The extended compression CSV/codebook/text
exchange outputs are not enabled by default in this project.

## Outputs and validation

All VTK files are written under the chosen `--dir`:

| Mode | Final VTK basename |
|---|---|
| EMI semi-implicit Euler | `emiLast.vtu` |
| SDC | `emiSDCLast.vtu` |
| SDC + BDDC | `emiSDCBDDCLast.vtu` |
| SDC + BDDC + compression | `emiSDCBDDCLastCompression.vtu` |
| BDDC | `emiBDDCLast.vtu` |
| BDDC + compression | `emiBDDCLastCompression.vtu` |

The initial field is written as `emiInitial.vtu` when VTK output is enabled.
Ordinary EMI also writes per-step fields; ordinary SDC writes per-step fields
when `vtk=1`. SDC+BDDC currently writes the initial and final fields.

For a small-mesh smoke test, run each target in `steps_to_run.mk` with the
default one-step cap and check that it exits successfully, reports finite
matrix diagnostics, and creates the expected final VTK file. For numerical
comparisons, use identical mesh, time-step, SDC, and BDDC tolerances:

- Compare SDC with SDC+BDDC without AA to measure the BDDC solve tolerance's effect.
- Compare AA on/off for the same solver and verify that both reach the SDC stopping criterion.
- Compare compressed and uncompressed BDDC solutions; report solution/error norms as well as transfer bytes.
- Increase `maximumNumberOfTimeSteps` (or set it to `0`) for a production-length run after the smoke test.

These are validation recommendations, not claims that all comparisons have
already been run. Compression quality and solver convergence should be measured
for each mesh and parameter set.

For SDC+BDDC, `--profile 1` reports cumulative wall time in selected phases:
local operator extraction, residual assembly, BDDC setup, BDDC iterations,
SDC updates, and AA selection. These categories do not cover every operation
(for example, VTK output), so they need not sum to total runtime. Use `--vtk 0`
for solver-focused timing and compare runs with the same mesh, thread count,
solver tolerances, and time steps. Profiling is disabled by default.
