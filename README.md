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

AA selects active subdomains for later SDC sweeps. If one active degree of
freedom belongs to a subdomain, the implementation activates the complete
local subdomain and neighboring owners required by the matrix neighborhood.
This preserves the subdomain-wise BDDC solve structure.
