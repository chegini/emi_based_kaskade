#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -ne 2 ]; then
  echo "Usage: RUN_NODE=htc-cmp024 $0 weak|strong 8|16|32|64" >&2
  exit 2
fi

mode="$1"
ranks="$2"
refine="${REFINE:-0}"

case "$mode:$ranks" in
  weak:8|weak:16|weak:32|weak:64)
    mesh="../input_emi_mesh_blocks/weak_strong_scaling/weak_${ranks}_h16_kaskade_mesh.vtu"
    extra="../input_emi_mesh_blocks/weak_strong_scaling/weak_${ranks}_h16_kaskade_list_extracellular.txt"
    intra="../input_emi_mesh_blocks/weak_strong_scaling/weak_${ranks}_h16_kaskade_list_intracellular.txt"
    excited="../input_emi_mesh_blocks/weak_strong_scaling/weak_${ranks}_h16_kaskade_early_excited.txt"
    ;;
  strong:8|strong:16|strong:32|strong:64)
    mesh="../input_emi_mesh_blocks/weak_strong_scaling/strong_h16_kaskade_mesh.vtu"
    extra="../input_emi_mesh_blocks/weak_strong_scaling/strong_h16_kaskade_list_extracellular.txt"
    intra="../input_emi_mesh_blocks/weak_strong_scaling/strong_h16_kaskade_list_intracellular.txt"
    excited="../input_emi_mesh_blocks/weak_strong_scaling/strong_h16_kaskade_early_excited.txt"
    ;;
  *)
    echo "Invalid case: use weak|strong and ranks 8|16|32|64" >&2
    exit 2
    ;;
esac

if [ "$ranks" -le 32 ]; then
  nodes=1
  tasks_per_node="$ranks"
else
  nodes=1
  tasks_per_node="$ranks"
fi

node_args=()
if [ -n "${RUN_NODE:-}" ]; then
  node_args+=(--nodelist="$RUN_NODE")
fi

tag="$mode-$ranks-r$refine"
  export_spec="ALL,MPI_RANKS=$ranks,THREADS=1,CPUS_PER_RANK=1,BITS=12,COMPRESSION=1,REFINE=$refine,FINAL_TIME=0.1,MAXIMUM_NUMBER_OF_TIME_STEPS=10,GRAPH_LIFTING=1,GRAPH_LIFTING_AFTER_QUANTIZATION=1,HUFFMAN=1,BITLENGTH=1,ZIGZAG=1,HISTOGRAM=0,MIN_INTERFACE_SIZE=0,THRESHOLD_INTERFACE_TYPES=4,INTERFACE_TYPES=7,EXCLUDE_CORNERS=1,BDD_ITERATIONS=300,MESH=$mesh,EXTRA=$extra,INTRA=$intra,EXCITED=$excited,CASE_TAG=$tag"

echo "===== submitting $tag: ranks=$ranks ====="
set +e
jid=$(sbatch --parsable --wait \
  --job-name="$tag" \
  --partition=big \
  --exclude=htc-cmp008 \
  "${node_args[@]}" \
  --nodes="$nodes" \
  --ntasks="$ranks" \
  --ntasks-per-node="$tasks_per_node" \
  --cpus-per-task=1 \
  --time=04:00:00 \
  --export="$export_spec" \
  run_emi_mpi_2nodes_8ranks.sbatch)
rc=$?
set -e

echo "Job ID: ${jid:-unknown}"
if [ -n "${jid:-}" ]; then
  sacct -n -X -j "$jid" \
    --format=JobID,JobName,State,ExitCode,Elapsed,MaxRSS,NodeList
fi
exit "$rc"
