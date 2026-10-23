# Integrated bulk analysis and report

The workflow connects existing upstream modules, Katana computations and the
approved seven-section report. It does not change the report layout, frozen
MotifTrace algorithms or existing study results. It is a development workflow;
clean-environment FASTQ/BAM production acceptance is still a separate release
requirement. It is not an approved nf-core pipeline.

## Execution

Install the development checkout on the execution machine using
`python -m pip install -e '.[flow,next]'`. Nextflow 25.10.4, a compatible Java,
the upstream tool environments, R packages, RSeQC, FIMO and configured resources
are required by their respective stages. Containers must actually contain those
tools; choosing a container profile does not build or certify an image.

Use a prepared project with completed upstream outputs to avoid repeating
sequencing processing:

```sh
nextflow run pipeline/main.nf -profile slurm \
  --mode downstream-only --project /data/study/project.json \
  --flow_config /data/study/flow.json --outdir /data/study/new-results -resume
```

Installed wheels contain the same workflow sources. Export once using
`katana-flow export-pipeline --outdir /data/tools/katana-pipeline`, then use
`nextflow run /data/tools/katana-pipeline/main.nf` with the same arguments.
An export never replaces an existing pipeline directory.

For a new bulk study:

```sh
nextflow run pipeline/main.nf -profile slurm \
  --mode complete --input /data/study/samples.csv \
  --contrasts /data/study/contrasts.csv \
  --fasta /references/genome.fa --gtf /references/annotation.gtf \
  --project_template /data/study/project-template.json \
  --flow_config /data/study/flow.json --outdir /data/study/results -resume
```

`--mode upstream-only` produces the verified native bundle and does not claim
to have completed Katana analysis or its report. Local execution uses
`-profile local`; Slurm account/partition settings belong in site configuration.
Do not use a local executor on a shared login node for heavy analyses.

## Calculation and result contracts

- Resource checks precede new upstream execution. A complete imported project
  is checked again before downstream tasks are scheduled.
- Only declared, selected contrast IDs are scheduled. The plan does not invent
  all pairwise comparisons or infer experimental purposes from sample names.
- Each task has its own manifest, output directory and catalog. Parallel jobs
  do not write concurrently to one DuckDB file.
- AS read support, landscape, DEG/DAS/DTU/DEXSeq integration, exact overlaps,
  four evidence-specific enrichments and RBP identity use the existing
  scientific workflows. Each comparison has independent MotifTrace and
  boundary/CLIP tasks. Backgrounds match detected/testable evidence.
- Gene-body QC is calculated from BAM and matched BED12; it is not inferred
  from FastQC. BED12 is generated from the matching GTF when no explicit BED
  is supplied. Coverage summaries wait for all declared sample measurements.
- Network discovery, module bootstrap, associations and annotation have
  explicit dependencies. Numerical linkage and leaf IDs are saved before
  plotting, using the unchanged fit. Edge-support bootstrap and module
  bootstrap remain distinct.
- Network configuration can be generated from the imported project's native
  rMATS sources and sample-vector identities. An explicit existing YAML may
  override configuration. An execution-only compatibility wrapper restores
  integral ordinal-rank indices for dynamicTreeCut; ranks, ordering, fitting
  parameters and frozen source files remain unchanged.
- Interaction tasks require explicit designs and real counts. Simple contrast
  differences are not relabeled as fitted interactions.
- MAEI prepares exact annotated SE inclusion/skipping transcript TPM and the
  official 5,440 kmer counts in the existing 14-window family, with one
  independent task per comparison/window. It then executes the official
  external MAPP model, Z-score and native mixture
  significance chain. It requires hash-bound source scripts and matching
  reference annotation, genomic FASTA and existing transcript quantification.
  An explicit already-prepared quantity/count input remains supported. It
  does not run MAPP's other eight
  subworkflows, create TPM from rMATS PSI, or claim calibrated MotifTrace/MAPP
  parallel evidence from an uncalibrated native fit. New-project preparation
  and native fitting have component tests; complete report integration and
  parallel calibration remain separate acceptance requirements.
- Completed task receipts bind parameters, method versions, input and output
  SHA256 and dependency identities. Unknown/failed outputs cannot be adopted
  as success. `-resume` reuses the exact Nextflow task cache; prior independent
  studies require an explicit verified adapter and are not automatically
  relabeled as this run.
- Collection waits for every required task. A partial/test collection cannot
  become a complete report. The report acceptance gate rejects missing fixed
  computations and missing chapter figures.

## Configuration

`project-template.json` supplies `metadata.species`, resource `assets`, optional
`sample_attributes`, and explicit JC/JCEC selection. Sample identities,
contrast direction and rMATS vector order come from the verified upstream
bundle. Resources retain their own species, build, biosample and source
provenance. Reference CLIP is not a measurement in the current experiment.

`flow.json` uses schema `katana.flow.config/v1`:

```json
{
  "schema": "katana.flow.config/v1",
  "contrast_ids": ["B_vs_A"],
  "event_set": "stable_significant",
  "parameters": {"differential.deg": {"fdr": 0.05, "min_abs_log2fc": 1.0}},
  "qc": {"reference_bed": "reference.bed", "rseqc_executable": "geneBody_coverage.py"},
  "network": {},
  "design_factors": ["condition"],
  "interactions": [],
  "mapp": {"adapter": "mapp.maei", "inputs": {}, "parameters": {}}
}
```

The empty MAPP resource inputs deliberately fail readiness. Supply the
official `native_root`, exact `source_sha256` bindings for its seven MAEI
files, and the official kmer `counter_sha256`. The workflow then prepares
included/total TPM and motif counts from the project's matching inputs.
Do not interpret this resource-free example as an executable complete study.
Declared multifactor studies need explicit `interactions` entries with `id`,
`design_path`, `reference_a` and `reference_b`. Previously approved thresholds
remain in [THRESHOLDS.md](THRESHOLDS.md). Scientific applicability is recorded
separately from a failed or unavailable computation.

Inspect without executing:

```sh
katana-flow preflight --project project.json --config flow.json
katana-flow plan --project project.json --config flow.json --out plan.json
```

The fixed report stays read-only. `katana-flow report` is the strict completion
entrypoint; `katana-report` may still render honestly labeled historical partial
collections. Neither command changes source statistics or invents missing data.

The calibrated completed MAPP/MotifTrace parallel import is a separate evidence
route. Its new-project calibration is not yet connected to the complete-mode
gate. Native MAEI fitting is available and tested, but does not remove that
remaining full-report acceptance requirement. A native fit is never relabeled
as calibrated parallel evidence.
