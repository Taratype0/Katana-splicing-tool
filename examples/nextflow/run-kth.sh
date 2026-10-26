#!/usr/bin/env bash
# Run from any directory. Edit the parameter YAML; this script submits no jobs
# until you actually run it. Load your site's Nextflow/Apptainer modules first.
set -euo pipefail
repo=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)
: "${SLURM_ACCOUNT:?Set SLURM_ACCOUNT to your own allocation}"
: "${PARAMS_FILE:?Set PARAMS_FILE to your edited params.complete.yaml}"
: "${WORK_DIR:?Set WORK_DIR to your study's Nextflow work directory}"
nextflow run "$repo/pipeline/main.nf" \
  -params-file "$PARAMS_FILE" \
  -profile slurm,apptainer \
  -c "$repo/examples/nextflow/kth.config" \
  --slurm_account "$SLURM_ACCOUNT" \
  -work-dir "$WORK_DIR" -resume "$@"
