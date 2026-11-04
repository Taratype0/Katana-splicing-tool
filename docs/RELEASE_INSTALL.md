# Install Katana

Use the [independent runtime guide](INDEPENDENT_RUNTIME.md) for the current
release, downloadable Python wheel, Linux runtime, matching task images and
checksums. Older releases remain available under their versioned GitHub tags.

- **FASTQ or BAM:** `nextflow -bg run Taratype0/Katana-splicing-tool -r v2.0.0-beta.5.dev2026110404`
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

[Run guide](NEXTFLOW_SIMPLE.md) | [Parameters](NEXTFLOW_PARAMETERS.md) 路
[Samples and contrasts](NEXTFLOW_USAGE.md) | [Outputs and recovery](NATIVE_OUTPUTS.md)

The approved-report update uses wheel/controller `2.0.0b5.dev2026110404`.
The report stage uses its matching image. Analysis and retention keep the
validated `2.0.0b5.dev2026110403` image and existing database catalog, so a
report-only upgrade preserves their cache identities. Python/R scientific code
and task allocations are unchanged.
See [task allocations and upstream versions](RESOURCE_ALLOCATION.md).

The approved seven-chapter report uses the 4 November faceted GO/KEGG layout. The report stage is pinned separately to runtime 2.0.0b5.dev2026110404; other analysis stages retain 2.0.0b5.dev2026110403 for cache reuse. Native completed enrichment tables are redrawn without refitting. All network modules remain browsable; association and pathway significance are separate filters.

