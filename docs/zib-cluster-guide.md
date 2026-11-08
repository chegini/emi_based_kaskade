# EMI/BDDC on the ZIB cluster: image, build, validation and results

Updated 7 November 2026. Intended for a colleague starting with a Linux
workstation and a ZIB cluster account. Commands use Bash.

The workflow is:

1. Obtain matching source trees and meshes.
2. Obtain or build an MPI-enabled Podman image on a Linux workstation.
3. Save the image and convert it to an Apptainer SIF on the cluster.
4. Allocate compute nodes; build Kaskade and EMI inside that image.
5. Test MPI transport, then run the small case and Pepe.
6. Inspect and preserve the logs, numerical errors and configuration.

**If the tested SIF and compatible Kaskade library are already available,
start at section 4, then build EMI in section 5 and run section 6.**

## 1. Files and source versions

Use your own writable project directory on storage visible from every compute
node. The required layout is:

```text
workspace/
  kaskade7_new/
  BDDC/
    emiModelBddcSdc_new_scratch/
    input_emi_mesh/
    input_emi_mesh_blocks/           optional for other experiments
  containers/
```

Obtain the **matching working-tree snapshots** of Kaskade and EMI from Fatemeh,
including the transport probe and local-bypass changes. A checkout of `main`
or `master` alone is not sufficient evidence that these features are present.
At the time of writing, some transport/bypass changes are still uncommitted
or unpublished. No push is needed to run a supplied snapshot.

For a fresh image/toolchain, use a fresh Linux build tree without old object
files or libraries. Do not mix macOS objects, another image's libraries or
another MPI installation with this build. Keep existing experimental trees
and results intact.

Pepe requires these files under `BDDC/input_emi_mesh/`:

```text
pepe_combi_domi.vtu
pepe_combi_domi_extracellular.txt
pepe_combi_domi_intracellular.txt
pepe_combi_domi_excited.txt
```

The small validation case uses files in the EMI repository's `input/` folder:

```text
10Cells3d_10extra_mesh.vtu
10Cells3d_10extra_list_extracellular.txt
10Cells3d_10extra_list_intracellular.txt
10Cells3d_10extra_early_excited.txt
```

## 2. Obtain or build the Podman image — Linux workstation

Run this section on a Linux x86_64 workstation with Podman, sufficient disk
space and network access, not as a build job on the cluster login node.

### Recommended for reproducing the existing experiments

Use the supplied tested image archive or SIF. The image used in the recorded
experiments was:

```text
localhost/kaskade7-robin-cluster:20261105
```

If this image already exists in your workstation's Podman storage, save it:

```bash
emi_image=localhost/kaskade7-robin-cluster:20261105
emi_image_dir="$HOME/emi-container-export"
mkdir -p "$emi_image_dir"
podman image inspect "$emi_image"
podman save --format docker-archive \
  --output "$emi_image_dir/kaskade7-robin.tar" "$emi_image"
sha256sum "$emi_image_dir/kaskade7-robin.tar"
```

Use `podman save`, not `podman export`: this preserves the container image
and its metadata. A supplied archive can also be loaded with
`podman load --input /path/to/kaskade7-robin.tar`.

### Building a new image from the repository

This is a **new candidate environment**, not an exact recreation of the tested
20261105 image. The current repository recipe specifies openSUSE Leap 16.1,
GCC 16 and Dune 2.11; its external package repositories and upstream branches
are not fully pinned. Availability and compatibility must be checked by the
build and the validation sequence below. Keep the existing tested SIF.

Set `emi_kaskade` to the supplied source tree on the workstation:

```bash
emi_kaskade=/path/to/workspace/kaskade7_new
emi_recipe_dir=$(mktemp -d "$HOME/emi-image-build.XXXXXX")
cp -a "$emi_kaskade/InstallDependencies/." "$emi_recipe_dir/"

# Avoid building Dune specifically for the workstation CPU.
sed 's/-march=native/-march=x86-64 -mtune=generic/g' \
  "$emi_recipe_dir/Dockerfile.kaskade7-build-environment" \
  > "$emi_recipe_dir/Containerfile.base"

podman build \
  --file "$emi_recipe_dir/Containerfile.base" \
  --tag localhost/kaskade7-build-environment:emi-candidate \
  "$emi_recipe_dir"
```

Continue only after the base image builds successfully. Add Open MPI and
Python explicitly; the base recipe alone is not the MPI runtime used by the
cluster scripts:

```bash
cat > "$emi_recipe_dir/Containerfile.mpi" <<'EOF'
FROM localhost/kaskade7-build-environment:emi-candidate
RUN zypper --non-interactive refresh && \
    zypper --non-interactive install -y openmpi4 openmpi4-devel python3 && \
    zypper --non-interactive clean --all
ENV PATH=/usr/lib64/mpi/gcc/openmpi4/bin:/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin
ENV LD_LIBRARY_PATH=/usr/lib64/mpi/gcc/openmpi4/lib64:/usr/lib64/mpi/gcc/openmpi4/lib
RUN mpicxx --showme:version && python3 --version && \
    test -r /usr/lib64/mpi/gcc/openmpi4/include/mpi.h
EOF

podman build --file "$emi_recipe_dir/Containerfile.mpi" \
  --tag localhost/kaskade7-emi:mpi-candidate "$emi_recipe_dir"

emi_image=localhost/kaskade7-emi:mpi-candidate
emi_image_dir="$HOME/emi-container-export"
mkdir -p "$emi_image_dir"
podman save --format docker-archive \
  --output "$emi_image_dir/kaskade7-emi-candidate.tar" "$emi_image"
sha256sum "$emi_image_dir/kaskade7-emi-candidate.tar"
```

Preserve the recipe directory and image digest with the archive. If the
distribution's package names or MPI paths change, resolve that build failure
before proceeding; the project wrappers currently expect the Open MPI 4 paths
shown above. These fresh-image commands have been checked against the recipe,
but a fresh image has not been built as part of writing this guide.

## 3. Transfer the archive and create the SIF

Put the archive on cluster storage visible from every allocated node. If the
workstation already uses that shared filesystem, no copy is needed. Otherwise,
copy it from the workstation, replacing the username and destination:

```bash
scp "$emi_image_dir/kaskade7-robin.tar" \
  YOUR_USER@htc-login03:/path/to/workspace/containers/
```

For the fresh-image route, use `kaskade7-emi-candidate.tar` and a distinct SIF
name instead. Log in to the cluster and request a one-node preparation job:

```bash
ssh YOUR_USER@htc-login03
export LANG=C LC_ALL=C LC_CTYPE=C
salloc --partition=high-mem --nodes=1 --ntasks=1 \
  --cpus-per-task=2 --mem=24G --time=02:00:00 \
  --job-name=emi-prepare
```

The shell may still display `htc-login03`; this is normal for `salloc`.
Launch actual conversion/build work with `srun`:

```bash
emi_image_dir=/path/to/workspace/containers
emi_sif="$emi_image_dir/kaskade7-robin.sif"
srun --mpi=none --nodes=1 --ntasks=1 \
  apptainer build "$emi_sif" \
  "docker-archive:$emi_image_dir/kaskade7-robin.tar"
sha256sum "$emi_sif"
```

Do not overwrite a known-good SIF with a candidate image. If a SIF was supplied
already, skip conversion. A standard SIF is the image; the source, library
builds, meshes and results remain in the mounted shared directories.

## 4. Set cluster paths and allocate nodes

These are Fatemeh's current paths. A colleague should replace them with their
own supplied workspace and image paths:

```bash
export LANG=C LC_ALL=C LC_CTYPE=C
emi_kaskade=/data/numerik/people/fchegini/kaskade7_new
emi_bddc=/data/numerik/people/fchegini/project/fatemeh/MicroCard/current_papers/BDDC_petsc_SDC_2023/gatevariables/20Nov2023/BDDC
emi_repo="$emi_bddc/emiModelBddcSdc_new_scratch"
emi_sif=/data/numerik/people/fchegini/containers/kaskade-robin-20261105/kaskade7-robin.sif
export emi_kaskade emi_bddc emi_repo emi_sif
cd "$emi_repo"

test -s "$emi_sif"
test -s "$emi_bddc/input_emi_mesh/pepe_combi_domi.vtu"
test -f tests/run_emi_transport_probe.py
test -f tests/run_emi_mpi_routing_check.py
```

Each check should succeed before continuing. Re-establish these variables
after reconnecting. Record `git rev-parse HEAD` and `git status --short` in
both repositories; a commit hash alone does not describe a modified checkout.

If you are still in the preparation allocation, finish section 5 there, then
`exit` that allocation shell before requesting the two-node allocation. Do
not nest active allocations.

For the two-node tests:

```bash
salloc --partition=big --nodes=2 --ntasks=4 --ntasks-per-node=2 \
  --cpus-per-task=1 --mem=24G --time=02:00:00 \
  --job-name=emi-pepe-validation
echo "$SLURM_JOB_ID"
scontrol show hostnames "$SLURM_JOB_NODELIST"
srun --mpi=list
```

Expect two node names and an available PMIx plugin. The tested cluster exposes
`pmix_v5`. The observed `high-mem` and `small` partitions allow only one node;
use `big` for this two-node configuration. `--mem=24G` is per node. Check current
partition limits with `sinfo -o "%P %a %D %s %l"` if allocation is rejected.

## 5. Build on a compute node inside the SIF

### First-time Kaskade library build

Skip this only if the supplied Kaskade tree already has a library built with
this image and compatible source. The EMI build wrapper assumes that library
exists. A new image/toolchain needs a fresh Kaskade build tree.

In a fresh tree, install the container configuration if absent:

```bash
if [ ! -e "$emi_kaskade/Makefile.Local" ]; then
  cp "$emi_kaskade/Makefile.Local-container" "$emi_kaskade/Makefile.Local"
fi
```

An existing `Makefile.Local` must be inspected: it should use `/kaskade7` and
include the image's `/Makefile.Local`, not workstation-specific dependencies.
Do not replace a colleague's configuration blindly. Then build:

```bash
srun --mpi=none --nodes=1 --ntasks=1 --ntasks-per-node=1 \
  apptainer exec --cleanenv \
  --bind "$emi_kaskade:/kaskade7" --pwd /kaskade7 "$emi_sif" \
  /bin/bash --noprofile --norc -c '
    set -eu
    export LANG=C LC_ALL=C OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
    make -j1 kasklib
    test -s /kaskade7/libs/libkaskade.a
  '
```

### Build the EMI executable

Run after Kaskade is available, and repeat after relevant source/header changes:

```bash
srun --mpi=none --nodes=1 --ntasks=1 --ntasks-per-node=1 \
  apptainer exec --cleanenv \
  --bind "$emi_kaskade:/kaskade7" \
  --bind "$emi_bddc:/BDDC" \
  --pwd /BDDC/emiModelBddcSdc_new_scratch "$emi_sif" \
  /bin/bash tests/build_emi_mpi_routing_cluster.sh
```

Expect the transfer/unit checks to pass and `Cluster EMI build passed`.
The wrapper selects the image's `mpicxx` and supplies the MPI/MUMPS header
compatibility flags. Do not substitute the host PETSc/MPICH compiler or launcher.
Do not run `make clean` in the experimental EMI tree: that target removes
`output/` and other result files.

## 6. Run validation in the two-node allocation

Run these Python commands from the **host allocation shell**, with the working
directory set to `$emi_repo`. The scripts launch the numerical jobs with
Slurm and Apptainer. Do not wrap them in another `srun`, `mpirun` or container.
Stop at the first failed stage and inspect that stage's log.

### A. Transport probe

```bash
python3 tests/run_emi_transport_probe.py \
  --slurm-container "$emi_sif" --kaskade "$emi_kaskade"
```

Expect `build-and-mpi-info`, `one-node-no-tcp` and `two-nodes` to pass.
The one-node probe excludes TCP to exercise shared-memory transport; the
two-node probe additionally permits TCP. Both use Vader with its single-copy
mechanism disabled.

### B. Small application test

```bash
python3 tests/run_emi_mpi_routing_check.py \
  --case 10cells --ranks 4 --expected-nodes 2 \
  --slurm-container "$emi_sif" --kaskade "$emi_kaskade" \
  --slurm-btl self,vader,tcp --slurm-vader-single-copy none \
  --check-local-bypass
```

For a separate one-node allocation, use `--ranks 2 --expected-nodes 1`.
Do not tell the checker to expect two nodes when only one was allocated.

### C. Focused Pepe diagnostic

```bash
python3 tests/run_emi_mpi_routing_check.py \
  --case pepe --ranks 4 --expected-nodes 2 \
  --slurm-container "$emi_sif" --kaskade "$emi_kaskade" \
  --slurm-btl self,vader,tcp --slurm-vader-single-copy none \
  --setup-diagnostic-only
```

This runs one uncompressed CG solve with material-exchange progress logging.
It is a bounded setup/solve check, not the full routing comparison.

### D. Full Pepe routing and compression-bypass comparison

```bash
python3 tests/run_emi_mpi_routing_check.py \
  --case pepe --ranks 4 --expected-nodes 2 \
  --slurm-container "$emi_sif" --kaskade "$emi_kaskade" \
  --slurm-btl self,vader,tcp --slurm-vader-single-copy none \
  --check-local-bypass
```

This runs a direct reference and the CG, FCG and FGMRES comparisons, using
the runner's default 16-bit transform-first compressed configuration. It is
one physical step, with histograms off. Read the printed planned cases and
saved commands to verify the configuration.

| Case suffix | Meaning |
|---|---|
| `uncompressed-legacy` | Raw transfers using the original routing |
| `uncompressed-destination` | Raw transfers using destination routing |
| `compressed-split-legacy` | Compression and constant splitting, original routing |
| `compressed-split-destination` | Same compression, destination routing |
| `compressed-split-local-bypass` | Destination routing with same-node compression bypass |

Local bypass leaves local interfaces uncompressed. To preserve shared-copy
consistency, a DOF with multiple owners on the same node is kept exact on
all its copies, including remote copies. Thus some inter-node values also
remain exact. MPI still manages communication between processes on the same
node; this feature does not introduce application-managed shared-memory windows.

## 7. What has been validated, and what the results mean

As of this guide's date:

- The full local-bypass comparison passed for 10cells and Pepe on two nodes
  with `self,tcp`.
- The new shared-memory transport probes passed twice with
  `self,vader,tcp` and `--slurm-vader-single-copy none`.
- The focused Pepe uncompressed solve passed with those shared-memory settings.
- The full 16-case Pepe comparison also passed with these shared-memory
  settings: `output/emi-mpi-routing-pepe-20261107T120921Z-2561102`.
  See the [numerical and timing review](pepe-routing-review-20261107.md).

Vader remains enabled when single-copy is `none`; this selects its copying
fallback. It avoids the path associated with the earlier large-message stall
in this container configuration. The exact original cause is not established.
Inter-node communication here uses TCP; these are not InfiniBand/RDMA timings.

The runners print a new directory under `output/` and keep per-case logs and
summary data. For a run in progress, copy its printed directory into:

```bash
emi_run=/absolute/path/to/the/printed/output/directory
tail -n 60 "$emi_run/cg-uncompressed-legacy/run.log"
```

Check completion, true residuals, iterations, solution differences against
the direct and matched uncompressed references, and local/inter-node byte
counts where reported. `PASS` does not impose a scientific solution-error
threshold. Raw and centered errors answer different questions; retain both.
The current exported-vector errors are not automatically finite-element
Sobolev norms or membrane-voltage errors.

Logical bytes are application payload accounting, not measured network traffic.
Do not interpret same-node payloads as bytes crossing the network. Validation
timings include setup and checks; repeated controlled runs are needed for
performance claims. Preserve image checksum, both source snapshots, executable
identity, node layout, commands and logs with the results.

Longer time integration is a separate stage. Do not assume this one-step routing
runner supports arbitrary time-step counts. The existing experimental guards
restrict CG/FCG constant splitting to one step and FGMRES constant splitting
to at most ten steps. Establish the intended error tolerance before describing
a configuration as scientifically validated.

## 8. Run without keeping an interactive connection open

After image creation, builds and transport checks have succeeded, a batch job
can run the full Pepe validation. Save the following as `run-pepe.sbatch`,
replace the three paths, and submit from a directory where you can write logs:

```bash
#!/bin/bash
#SBATCH --partition=big
#SBATCH --nodes=2
#SBATCH --ntasks=4
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=02:00:00
#SBATCH --job-name=emi-pepe
#SBATCH --output=emi-pepe-%j.log
set -euo pipefail
export LANG=C LC_ALL=C LC_CTYPE=C
emi_kaskade=/path/to/workspace/kaskade7_new
emi_repo=/path/to/workspace/BDDC/emiModelBddcSdc_new_scratch
emi_sif=/path/to/workspace/containers/kaskade7-robin.sif
cd "$emi_repo"
python3 tests/run_emi_mpi_routing_check.py \
  --case pepe --ranks 4 --expected-nodes 2 \
  --slurm-container "$emi_sif" --kaskade "$emi_kaskade" \
  --slurm-btl self,vader,tcp --slurm-vader-single-copy none \
  --check-local-bypass
```

Submit and monitor on the login node:

```bash
sbatch run-pepe.sbatch
squeue -u "$USER"
```

This batch template does not build or alter the image. Do not rebuild the
shared executable while another job is using that working tree.

## 9. Troubleshooting and finishing

| Symptom | Action |
|---|---|
| `Node count specification invalid` | Check partition limits; use `big` for two nodes. |
| Expired/invalid Slurm job ID | Leave the old allocation shell and obtain a new allocation. Do not reuse the old job ID. |
| Large-message stall with default Vader | Use the tested probe and `--slurm-vader-single-copy none`. `--slurm-btl self,tcp` is the diagnostic fallback; it does not exercise same-node shared-memory transport. |
| PMIx `munge`/`shmem2` component messages | Preserve logs and check actual probe completion. Some recorded runs passed despite these messages. Do not change authentication settings merely to suppress them. |
| MPI startup fails under Apptainer | Use the supplied host runner. MPI execution preserves Slurm/PMIx environment; `--cleanenv` is used for builds, not those numerical launches. |
| Missing result files after failure | Inspect the first application/Slurm error. Missing summaries can be consequences of termination. |
| Exit 143 / `DUE TO TIME LIMIT` | Inspect job/step accounting and progress logs. Longer allocation time alone does not cure a transport stall; runner job steps may have their own limits. |
| Clock-skew build warning | Check shared-file and node clocks; retain build logs and verify the rebuild. Do not treat the warning alone as evidence of a solver error. |

Use `squeue -u "$USER"` to inspect jobs. For an interactive allocation, `exit`
the allocation shell when finished so the nodes are released. To cancel a
specific job, use `scancel JOB_ID` with the intended numeric job ID.

## References

The operational commands above were checked against the project scripts and
the supplied cluster logs. General container/MPI behavior is documented in:

- [Podman save](https://docs.podman.io/en/latest/markdown/podman-save.1.html).
- [Apptainer container builds](https://apptainer.org/docs/user/latest/build_a_container.html).
- [Apptainer MPI integration](https://apptainer.org/docs/user/latest/mpi.html).
- [Open MPI shared-memory transport](https://docs.open-mpi.org/en/main/tuning-apps/networking/shared-memory.html)
  (Open MPI 4 calls this BTL `vader`; newer versions use `sm`).

For the detailed project checks, also see `tests/emi_mpi_routing.txt`.
