#!/usr/bin/env bash
set -euo pipefail

# Run this script from the ZIB login node in the EMI scratch directory.
LAUNCHER="run_emi_mpi_2nodes_8ranks.sbatch"

submit_case() {
  tag="$1"
  ranks="$2"
  mesh="$3"
  extra="$4"
  intra="$5"
  excited="$6"

  if [ "$ranks" -le 32 ]; then
    nodes=1
    tasks_per_node="$ranks"
  else
    nodes=2
    tasks_per_node=32
  fi

  export_spec="ALL,MPI_RANKS=$ranks,THREADS=1,CPUS_PER_RANK=1,BITS=12,COMPRESSION=1,REFINE=0,FINAL_TIME=0.1,MAXIMUM_NUMBER_OF_TIME_STEPS=10,GRAPH_LIFTING=1,GRAPH_LIFTING_AFTER_QUANTIZATION=1,HUFFMAN=1,BITLENGTH=1,ZIGZAG=1,HISTOGRAM=0,MIN_INTERFACE_SIZE=0,THRESHOLD_INTERFACE_TYPES=4,INTERFACE_TYPES=7,EXCLUDE_CORNERS=1,BDD_ITERATIONS=300,MESH=$mesh,EXTRA=$extra,INTRA=$intra,EXCITED=$excited,CASE_TAG=$tag"
  node_args=()
  if [ -n "${RUN_NODE:-}" ]; then
    node_args+=(--nodelist="$RUN_NODE")
  fi

  echo "===== submitting $tag: ranks=$ranks ====="

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
    "$LAUNCHER")

  sacct -n -X -j "$jid" \
    --format=JobID,State,ExitCode,Elapsed,MaxRSS,NodeList
}

SCALING_DIR="../input_emi_mesh_blocks/weak_strong_scaling"

# Weak scaling: one mesh block per MPI rank.
submit_case \
  weak-8 8 \
  "$SCALING_DIR/weak_8_h16_kaskade_mesh.vtu" \
  "$SCALING_DIR/weak_8_h16_kaskade_list_extracellular.txt" \
  "$SCALING_DIR/weak_8_h16_kaskade_list_intracellular.txt" \
  "$SCALING_DIR/weak_8_h16_kaskade_early_excited.txt"

submit_case \
  weak-16 16 \
  "$SCALING_DIR/weak_16_h16_kaskade_mesh.vtu" \
  "$SCALING_DIR/weak_16_h16_kaskade_list_extracellular.txt" \
  "$SCALING_DIR/weak_16_h16_kaskade_list_intracellular.txt" \
  "$SCALING_DIR/weak_16_h16_kaskade_early_excited.txt"

submit_case \
  weak-32 32 \
  "$SCALING_DIR/weak_32_h16_kaskade_mesh.vtu" \
  "$SCALING_DIR/weak_32_h16_kaskade_list_extracellular.txt" \
  "$SCALING_DIR/weak_32_h16_kaskade_list_intracellular.txt" \
  "$SCALING_DIR/weak_32_h16_kaskade_early_excited.txt"

submit_case \
  weak-64 64 \
  "$SCALING_DIR/weak_64_h16_kaskade_mesh.vtu" \
  "$SCALING_DIR/weak_64_h16_kaskade_list_extracellular.txt" \
  "$SCALING_DIR/weak_64_h16_kaskade_list_intracellular.txt" \
  "$SCALING_DIR/weak_64_h16_kaskade_early_excited.txt"

# Strong scaling: fixed global mesh, increasing MPI ranks.
STRONG_MESH="$SCALING_DIR/strong_h16_kaskade_mesh.vtu"
STRONG_EXTRA="$SCALING_DIR/strong_h16_kaskade_list_extracellular.txt"
STRONG_INTRA="$SCALING_DIR/strong_h16_kaskade_list_intracellular.txt"
STRONG_EXCITED="$SCALING_DIR/strong_h16_kaskade_early_excited.txt"

for ranks in 8 16 32 64; do
  submit_case \
    "strong-$ranks" "$ranks" \
    "$STRONG_MESH" \
    "$STRONG_EXTRA" \
    "$STRONG_INTRA" \
    "$STRONG_EXCITED"
done
