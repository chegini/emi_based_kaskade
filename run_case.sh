#!/usr/bin/env bash
set -euo pipefail

BITS=${BITS:-32}
LIFTING_MODE=${LIFTING_MODE:-before_quantization}
ZIGZAG=${ZIGZAG:-1}
HUFFMAN=${HUFFMAN:-1}
BITLENGTH=${BITLENGTH:-1}
FINAL_TIME=${FINAL_TIME:-0.1}

case "$LIFTING_MODE" in
  none)
    LIFTING=0
    AFTER_QUANTIZATION=0
    LIFTING_NAME="no-graph-lifting"
    ;;
  before_quantization)
    LIFTING=1
    AFTER_QUANTIZATION=0
    LIFTING_NAME="graph-lifting-before-quantization"
    ;;
  after_quantization)
    LIFTING=0
    AFTER_QUANTIZATION=1
    LIFTING_NAME="graph-lifting-after-quantization"
    ;;
  *)
    printf 'Invalid LIFTING_MODE=%s\n' "$LIFTING_MODE" >&2
    printf 'Use: none, before_quantization, or after_quantization\n' >&2
    exit 2
    ;;
esac

MESH=${MESH:-../input_emi_mesh_blocks/SX_h16_original_mesh.vtu}
EXTRA=${EXTRA:-../input_emi_mesh_blocks/SX_h16_original_list_extracellular.txt}
INTRA=${INTRA:-../input_emi_mesh_blocks/SX_h16_original_list_intracellular.txt}
EXCITED=${EXCITED:-../input_emi_mesh_blocks/SX_h16_original_early_excited.txt}

CASE="${LIFTING_NAME}_bits-${BITS}_zigzag-${ZIGZAG}_huffman-${HUFFMAN}_bitlength-${BITLENGTH}"
RUN="output/session13/SX_h16/${CASE}"
LOG="logs/session13/SX_h16-${CASE}.log"

mkdir -p "$RUN" "$(dirname "$LOG")"

./emiModel \
  --mpi 0 \
  --input "$MESH" \
  --extra_set "$EXTRA" \
  --intra_set "$INTRA" \
  --excited "$EXCITED" \
  --bddc 1 \
  --bddcCompression 1 \
  --bddcCompressionBits "$BITS" \
  --bddcGraphLifting "$LIFTING" \
  --bddcGraphLiftingAfterQuantization "$AFTER_QUANTIZATION" \
  --bddcZigzag "$ZIGZAG" \
  --bddcHuffman "$HUFFMAN" \
  --bddcBitlength "$BITLENGTH" \
  --bddcCompressionReport 1 \
  --bddcTolerance 1e-7 \
  --bddcIterations 20 \
  --bddcVerbose 1 \
  --dt 0.01 \
  --finalTime "$FINAL_TIME" \
  --dir "$RUN" \
  2>&1 | tee "$LOG"
