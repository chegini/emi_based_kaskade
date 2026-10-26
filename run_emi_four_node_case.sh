#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -ne 1 ]; then
  echo "Usage: RUN_NODE=nodeA,nodeB,nodeC,nodeD $0 weak|strong" >&2
  exit 2
fi

mode="$1"
ranks=128
nodes=4
tasks_per_node=32
refine="${REFINE:-0}"
run_nodes="${RUN_NODE:-}"
node_args=()

if [ -n "$run_nodes" ]; then
  node_count=$(awk -F, '{print NF}' <<< "$run_nodes")
  if [ "$node_count" -ne 4 ]; then
    echo "RUN_NODE must contain exactly four comma-separated nodes." >&2
    echo "Example: RUN_NODE=htc-cmp103,htc-cmp105,htc-cmp106,htc-cmp108 $0 weak" >&2
    exit 2
  fi
  node_args+=(--nodelist="$run_nodes")
fi

case "$mode" in
  weak)
    mesh="../input_emi_mesh_blocks/weak_strong_scaling/weak_128_h16_kaskade_mesh.vtu"
    extra="../input_emi_mesh_blocks/weak_strong_scaling/weak_128_h16_kaskade_list_extracellular.txt"
    intra="../input_emi_mesh_blocks/weak_strong_scaling/weak_128_h16_kaskade_list_intracellular.txt"
    excited="../input_emi_mesh_blocks/weak_strong_scaling/weak_128_h16_kaskade_early_excited.txt"
    ;;
  strong)
    mesh="../input_emi_mesh_blocks/weak_strong_scaling/strong_h16_kaskade_mesh.vtu"
    extra="../input_emi_mesh_blocks/weak_strong_scaling/strong_h16_kaskade_list_extracellular.txt"
    intra="../input_emi_mesh_blocks/weak_strong_scaling/strong_h16_kaskade_list_intracellular.txt"
    excited="../input_emi_mesh_blocks/weak_strong_scaling/strong_h16_kaskade_early_excited.txt"
    ;;
  *)
    echo "Usage: RUN_NODE=nodeA,nodeB,nodeC,nodeD $0 weak|strong" >&2
    exit 2
    ;;
esac

tag="$mode-$ranks-r$refine"
export_spec="ALL,MPI_RANKS=$ranks,THREADS=1,CPUS_PER_RANK=1,BITS=12,COMPRESSION=1,REFINE=$refine,FINAL_TIME=0.1,MAXIMUM_NUMBER_OF_TIME_STEPS=10,GRAPH_LIFTING=1,GRAPH_LIFTING_AFTER_QUANTIZATION=1,HUFFMAN=1,BITLENGTH=1,ZIGZAG=1,HISTOGRAM=0,MIN_INTERFACE_SIZE=0,THRESHOLD_INTERFACE_TYPES=4,INTERFACE_TYPES=7,EXCLUDE_CORNERS=1,BDD_ITERATIONS=300,MESH=$mesh,EXTRA=$extra,INTRA=$intra,EXCITED=$excited,CASE_TAG=$tag"

echo "===== submitting $tag on $run_nodes ====="

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
