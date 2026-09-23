# Reproducible Kaskade7 and MPI setup

This document describes how to build and run the EMI example from a fresh
login on a machine with Podman. The commands are divided into **host** and
**container** steps. The source trees are mounted into the container, so no
source copy is required.

## 1. Host prerequisites

On the host machine, verify Podman and the MPI launcher:

```bash
command -v podman
command -v mpirun
```

The MPI launcher is only needed for the MPI smoke test. A normal Kaskade
build does not require MPI.

Set paths appropriate for the machine:

```bash
KASKADE=/data/numerik/people/fchegini/kaskade7_new
BDDC=/data/numerik/people/fchegini/project/fatemeh/MicroCard/current_papers/BDDC_petsc_SDC_2023/gatevariables/20Nov2023/BDDC
```

The `BDDC` directory must contain `emiModelBddcSdc_new_scratch` and the
`input_emi_mesh` directory.

## 2. Load the base image

Check the local images:

```bash
podman images
```

The base image used by this project is:

```text
localhost/kaskade7-build-environment:latest
```

If it is distributed as an archive, load it on the host with:

```bash
gunzip -c /path/to/kaskade7-build-environment.tar.gz | podman load
```

## 3. Build the MPI-enabled image

The base image is an openSUSE Tumbleweed image. Create a derived image once
on each machine where MPI builds are needed:

```bash
mkdir -p ~/kaskade7-mpi-image
cd ~/kaskade7-mpi-image
```

Create `Containerfile` with this content:

```Dockerfile
FROM localhost/kaskade7-build-environment:latest

USER root

RUN zypper --non-interactive install --no-recommends \
    openmpi4 openmpi4-devel

ENV PATH=/usr/lib64/mpi/gcc/openmpi4/bin:$PATH
ENV LD_LIBRARY_PATH=/usr/lib64/mpi/gcc/openmpi4/lib64:/usr/lib64/mpi/gcc/openmpi4/lib:$LD_LIBRARY_PATH
```

Build it:

```bash
podman build -t localhost/kaskade7-mpi-build-environment:latest .
```

Some Podman setups do not preserve the MPI `PATH` from `Containerfile`.
The explicit environment variables in the next section handle that case.

## 4. Start the development container

Run this on the host, from any directory:

```bash
podman run -it --rm \
  -e PATH=/usr/lib64/mpi/gcc/openmpi4/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin \
  -e LD_LIBRARY_PATH=/usr/lib64/mpi/gcc/openmpi4/lib64:/usr/lib64/mpi/gcc/openmpi4/lib \
  -v "${KASKADE}:/kaskade7" \
  -v "${BDDC}:/BDDC" \
  -w /BDDC/emiModelBddcSdc_new_scratch \
  localhost/kaskade7-mpi-build-environment:latest
```

Inside the container, verify the MPI tools:

```bash
command -v mpicxx
command -v mpirun
mpicxx --version
mpirun --version
```

The expected paths are under:

```text
/usr/lib64/mpi/gcc/openmpi4/bin
```

If `vi` is unavailable in the image, edit files on the host through the
mounted source tree or use `sed` for small changes.

## 5. Normal non-MPI build and run

Inside the container:

```bash
make -f steps_to_run.mk clean
make -f steps_to_run.mk build
```

Run a quick EMI check:

```bash
make -f steps_to_run.mk run-emi
```

Available scenario targets are listed with:

```bash
make -f steps_to_run.mk help
```

The default build does not include MPI and remains usable on machines with no
MPI installation.

## 6. MPI-enabled build

Clean first because `make` does not know that a compiler or preprocessor flag
changed between builds:

```bash
make -f steps_to_run.mk clean
make -f steps_to_run.mk build-mpi MPICXX=mpicxx
```

The scratch `Makefile` deliberately puts the OpenMPI include directory before
`/usr/include/mumps`. MUMPS provides a compatibility file named `mpi.h`; if
it is found first, the MPI API is incomplete and compilation fails with
errors such as `MPI_Comm_size` or `MPI_Allreduce` not declared.

If the OpenMPI installation uses another path, override it when building:

```bash
make -f steps_to_run.mk build-mpi \
  MPICXX=mpicxx \
  MPI_INCLUDE_FLAGS=-I/path/to/openmpi/include
```

## 7. MPI smoke test

The current MPI support is only a startup and communication test. It does
not distribute BDDC subdomains yet. Each MPI rank runs the EMI model
independently after participating in an `MPI_Allreduce`.

Inside the container:

```bash
mkdir -p output

env -u DISPLAY -u XAUTHORITY \
  mpirun --allow-run-as-root -np 2 \
  ./emiModel \
  --mpi 1 \
  --bddc 0 \
  --vtk 0 \
  --finalTime 0.01 \
  --dt 0.01 \
  --maximumNumberOfTimeSteps 1 \
  --nThreads 2 \
  --dir output/mpi_smoke
```

A successful smoke test prints one line per rank, for example:

```text
MPI smoke test: rank 0/2, allreduce sum=3
MPI smoke test: rank 1/2, allreduce sum=3
```

The EMI model output may appear twice because both ranks currently execute
the model independently. That is expected for this smoke test.

The helper target is also available:

```bash
make -f steps_to_run.mk run-mpi-smoke MPI_RANKS=2 MPICXX=mpicxx
```

If running as root inside Podman, add `--allow-run-as-root` to the `mpirun`
command. If the helper target needs this permanently, set the MPI launcher
arguments in the local make configuration.

## 8. Returning to a normal build

MPI and non-MPI objects must not be mixed. Switch back explicitly:

```bash
make -f steps_to_run.mk clean
make -f steps_to_run.mk build
```

## 9. SSH forwarding warnings

Messages such as these concern SSH/X11 forwarding, not MPI or the EMI solver:

```text
Address already in use
Cannot assign requested address
Authorization required, but no authorization protocol specified
```

Use a different local forwarding port or keep the existing SSH tunnel. For
MPI command-line tests, unset display variables as shown above.
