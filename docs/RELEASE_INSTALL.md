# Install Katana

Use the [independent runtime guide](INDEPENDENT_RUNTIME.md) for the current
release, downloadable Python wheel, Linux runtime, matching task images and
checksums. Older releases remain available under their versioned GitHub tags.

- **FASTQ or BAM:** `nextflow -bg run Taratype0/Katana-splicing-tool -r v2.0.0-beta.5.dev2026110301`
  with sample/contrast sheets, references, profile and output paths. Containers
  and supported databases download automatically; no manual JSON is required.
- **Existing upstream bundle:** the same command with `--mode downstream-only
  --bundle /data/katana_bundle`.
- **Imported Katana project:** `--mode downstream-only --project /data/project.json`.
- **Browser and selected targets:** install the wheel and follow the
  [browser/target guide](TARGET_COMMANDS.md).

All new pipeline stages are scheduled by Nextflow. There is no external Python
stage coordinator. Task images supply native Python/R tools; the standalone
Python package can also use the [Conda environment](RUNTIME_ENVIRONMENTS.md).

[Run guide](NEXTFLOW_SIMPLE.md) 路 [Parameters](NEXTFLOW_PARAMETERS.md) 路
[Samples and contrasts](NEXTFLOW_USAGE.md) 路 [Outputs and recovery](NATIVE_OUTPUTS.md)

The resource-policy update uses wheel/controller `2.0.0b5.dev2026110301`.
It intentionally reuses the validated `20261103` downstream image and database
catalog: Python/R scientific code is unchanged. The new Nextflow source carries
the updated task allocations; upgrading only an image does not update a workflow.
See [task allocations and upstream versions](RESOURCE_ALLOCATION.md).
