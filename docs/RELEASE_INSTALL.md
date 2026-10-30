# Install Katana

Use the [independent runtime guide](INDEPENDENT_RUNTIME.md) for the current
release, downloadable Python wheel, Linux runtime, matching task images and
checksums. Older releases remain available under their versioned GitHub tags.

- **FASTQ or BAM:** `nextflow run tools/katana-controller/pipeline/main.nf` with
  sample/contrast sheets, references, resources, profile and output paths.
- **Existing upstream bundle:** the same command with `--mode downstream-only
  --bundle /data/katana_bundle`.
- **Imported Katana project:** `--mode downstream-only --project /data/project.json`.
- **Browser and selected targets:** install the wheel and follow the
  [browser/target guide](TARGET_COMMANDS.md).

All new pipeline stages are scheduled by Nextflow. There is no external Python
stage coordinator. Task images supply native Python/R tools; the standalone
Python package can also use the [Conda environment](RUNTIME_ENVIRONMENTS.md).

[Run guide](NEXTFLOW_SIMPLE.md) · [Parameters](NEXTFLOW_PARAMETERS.md) ·
[Samples and contrasts](NEXTFLOW_USAGE.md) · [Outputs and recovery](NATIVE_OUTPUTS.md)
