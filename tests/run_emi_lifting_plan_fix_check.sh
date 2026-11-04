#!/usr/bin/env bash
set -euo pipefail
cd /BDDC/emiModelBddcSdc_new_scratch
bash tests/run_emi_transfer_roundtrip_audit.sh \
  --reference-run output/emi-cg-split-check-pepe-quantization-only-20261104T171014Z-6196/cg16 \
  "$@"
