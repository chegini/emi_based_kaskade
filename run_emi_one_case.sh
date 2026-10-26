#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -ne 2 ]; then
  echo "Usage: RUN_NODE=htc-cmp024 $0 weak|strong|pepe|rizzo|robin|kermit|gonzo|animal 8|16|32|64" >&2
  exit 2
fi

mode="$1"
ranks="$2"
refine="${REFINE:-0}"
histogram="${HISTOGRAM:-0}"
compression="${COMPRESSION:-1}"
vtk="${VTK:-0}"

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
  pepe:8|pepe:16|pepe:32|pepe:64)
    mesh="../input_emi_mesh/pepe_combi_domi.vtu"
    extra="../input_emi_mesh/pepe_combi_domi_extracellular.txt"
    intra="../input_emi_mesh/pepe_combi_domi_intracellular.txt"
    excited="../input_emi_mesh/pepe_combi_domi_excited.txt"
    ;;
  rizzo:8|rizzo:16|rizzo:32|rizzo:64)
    mesh="../input_emi_mesh/rizzo_combi_domi.vtu"
    extra="../input_emi_mesh/rizzo_combi_domi_extracellular.txt"
    intra="../input_emi_mesh/rizzo_combi_domi_intracellular.txt"
    excited="../input_emi_mesh/rizzo_combi_domi_excited.txt"
    ;;
  robin:8|robin:16|robin:32|robin:64)
    mesh="../input_emi_mesh/robin_combi_domi.vtu"
    extra="../input_emi_mesh/robin_combi_domi_extracellular.txt"
    intra="../input_emi_mesh/robin_combi_domi_intracellular.txt"
    excited="../input_emi_mesh/robin_combi_domi_excited.txt"
    ;;
  kermit:8|kermit:16|kermit:32|kermit:64)
    mesh="../input_emi_mesh/kermit_mesh.vtu"
    extra="../input_emi_mesh/kermit_extracellular.txt"
    intra="../input_emi_mesh/kermit_intracellular.txt"
    excited="../input_emi_mesh/kermit_early_excited.txt"
    ;;
  gonzo:8|gonzo:16|gonzo:32|gonzo:64)
    mesh="../input_emi_mesh/gonzo_mesh.vtu"
    extra="../input_emi_mesh/gonzo_extracellular.txt"
    intra="../input_emi_mesh/gonzo_intracellular.txt"
    excited="../input_emi_mesh/gonzo_early_excited.txt"
    ;;
  animal:8|animal:16|animal:32|animal:64)
    mesh="../input_emi_mesh/animal_mesh.vtu"
    extra="../input_emi_mesh/animal_extracellular.txt"
    intra="../input_emi_mesh/animal_intracellular.txt"
    excited="../input_emi_mesh/animal_early_excited.txt"
    ;;
  *)
    echo "Invalid case: use weak|strong|pepe|rizzo|robin|kermit|gonzo|animal and ranks 8|16|32|64" >&2
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

tag="$mode-$ranks-r$refine-h$histogram"
run_name="${RUN_NAME:-$tag}"
bits="${BITS:-12}"
export_spec="ALL,MPI_RANKS=$ranks,THREADS=1,CPUS_PER_RANK=1,BITS=$bits,COMPRESSION=$compression,REFINE=$refine,FINAL_TIME=0.1,MAXIMUM_NUMBER_OF_TIME_STEPS=10,GRAPH_LIFTING=${GRAPH_LIFTING:-1},GRAPH_LIFTING_AFTER_QUANTIZATION=${GRAPH_LIFTING_AFTER_QUANTIZATION:-1},HUFFMAN=${HUFFMAN:-1},BITLENGTH=${BITLENGTH:-1},ZIGZAG=${ZIGZAG:-1},HISTOGRAM=$histogram,VTK=$vtk,RUN_NAME=$run_name,MIN_INTERFACE_SIZE=0,THRESHOLD_INTERFACE_TYPES=4,INTERFACE_TYPES=7,EXCLUDE_CORNERS=1,BDD_ITERATIONS=300,MESH=$mesh,EXTRA=$extra,INTRA=$intra,EXCITED=$excited,CASE_TAG=$tag"

echo "===== submitting $tag: ranks=$ranks ====="
echo "Bits: $bits, compression: $compression, refine: $refine, VTK: $vtk"
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
