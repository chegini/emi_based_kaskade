#!/usr/bin/env bash
set -euo pipefail

BITS=${BITS:-32}
MIN_INTERFACE_SIZE=${MIN_INTERFACE_SIZE:-0}
THRESHOLD_INTERFACE_TYPES=${THRESHOLD_INTERFACE_TYPES:-2}
LIFTING_MODE=${LIFTING_MODE:-2}
ZIGZAG=${ZIGZAG:-1}
HUFFMAN=${HUFFMAN:-1}
BITLENGTH=${BITLENGTH:-1}
FINAL_TIME=${FINAL_TIME:-0.1}
THREADS=${THREADS:-2}
HISTOGRAM=${HISTOGRAM:-0}
VTK=${VTK:-0}
EXCLUDE_CORNERS=${EXCLUDE_CORNERS:-0}

case "$EXCLUDE_CORNERS" in
  0|1) ;;
  *) printf 'Invalid EXCLUDE_CORNERS=%s; use 0 or 1\n' "$EXCLUDE_CORNERS" >&2; exit 2 ;;
esac

case "$LIFTING_MODE" in
  0)
    LIFTING=0
    AFTER_QUANTIZATION=0
    LIFTING_NAME="no-graph-lifting"
    ;;
  1)
    LIFTING=0
    AFTER_QUANTIZATION=1
    LIFTING_NAME="graph-lifting-after-quantization"
    ;;
  2)
    LIFTING=1
    AFTER_QUANTIZATION=0
    LIFTING_NAME="graph-lifting-before-quantization"
    ;;
  *)
    printf 'Invalid LIFTING_MODE=%s\n' "$LIFTING_MODE" >&2
    printf 'Use: 0 (none), 1 (after quantization), or 2 (before quantization)\n' >&2
    exit 2
    ;;
esac

MESH=${MESH:-../input_emi_mesh_blocks/SX_h16_original_mesh.vtu}
EXTRA=${EXTRA:-../input_emi_mesh_blocks/SX_h16_original_list_extracellular.txt}
INTRA=${INTRA:-../input_emi_mesh_blocks/SX_h16_original_list_intracellular.txt}
EXCITED=${EXCITED:-../input_emi_mesh_blocks/SX_h16_original_early_excited.txt}

CASE="${LIFTING_NAME}_bits-${BITS}_corners-${EXCLUDE_CORNERS}_zigzag-${ZIGZAG}_huffman-${HUFFMAN}_bitlength-${BITLENGTH}"
if [[ "$HISTOGRAM" == 1 ]]; then
  CASE="${CASE}_histogram"
fi
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
  --bddcCompressionMinInterfaceSize "$MIN_INTERFACE_SIZE" \
  --bddcCompressionThresholdInterfaceTypes "$THRESHOLD_INTERFACE_TYPES" \
  --bddcGraphLifting "$LIFTING" \
  --bddcGraphLiftingAfterQuantization "$AFTER_QUANTIZATION" \
  --bddcZigzag "$ZIGZAG" \
  --bddcHuffman "$HUFFMAN" \
  --bddcBitlength "$BITLENGTH" \
  --bddcExcludeCornersFromCompression "$EXCLUDE_CORNERS" \
  --bddcCompressionReport 1 \
  --bddcHistogram "$HISTOGRAM" \
  --bddcTolerance 1e-7 \
  --bddcIterations 20 \
  --bddcVerbose 1 \
  --nThreads "$THREADS" \
  --dt 0.01 \
  --finalTime "$FINAL_TIME" \
  --vtk "$VTK" \
  --dir "$RUN" \
  2>&1 | tee "$LOG"
