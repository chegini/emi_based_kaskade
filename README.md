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

The scenario helper is available through `steps_to_run.mk`:

```bash
make -f steps_to_run.mk help
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
| EMI + SDC + BDDC + AA | `run-sdc-bddc-aa` | adds AA to SDC+BDDC |
| EMI + SDC + BDDC + compression | `run-sdc-bddc-compression` | compressed BDDC transfer |
| EMI + SDC + BDDC + AA + compression | `run-sdc-bddc-aa-compression` | compressed transfer plus AA |

For example:

```bash
make -f steps_to_run.mk run-sdc-bddc-aa-compression
```

The helper defaults to two threads and one time step for quick checks. Override
the thread count with `NTHREADS`, for example `make -f steps_to_run.mk
run-sdc-bddc NTHREADS=4`. The one-step limit is a test setting in
`steps_to_run.mk`; increase `TIME_ARGS` there for longer runs. For direct
executable runs, `--nThreads` has the same two-thread default and
`--maximumNumberOfTimeSteps 0` means use all steps implied by `finalTime` and
`dt`.

To run all nine configurations sequentially:

```bash
make -f steps_to_run.mk run-all
```

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

The compressed path configures DCT transform, quantization, bit-length coding,
and Huffman coding for both restriction and prolongation. If compression is
disabled, the standard `SpaceTransfer` implementation is used.

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
