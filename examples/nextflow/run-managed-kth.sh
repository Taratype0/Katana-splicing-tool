#!/usr/bin/env bash
# One parameter file, native Slurm tasks, automatic successful-stage cleanup.
# Use the matching installed Katana package. This does not install environments.
set -euo pipefail
repo=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)
: "${SLURM_ACCOUNT:?Set SLURM_ACCOUNT to your own allocation}"
: "${PARAMS_FILE:?Set PARAMS_FILE to your edited params.complete.yaml}"
katana-flow run --pipeline "$repo/pipeline/main.nf" \
  -params-file "$PARAMS_FILE" -profile slurm,apptainer \
  -c "$repo/examples/nextflow/kth.config" \
  --slurm_account "$SLURM_ACCOUNT" "$@"
# Add -resume only to resume an existing managed run; --keep-work keeps caches.
