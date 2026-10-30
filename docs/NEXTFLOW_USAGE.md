# Katana Nextflow: usage

**Choose the entry that matches the files you already have.** FASTQ and genome
BAM can enter the complete bulk workflow. A verified upstream bundle can enter
downstream analysis without repeating read processing. Completed downstream
results can enter reporting directly. Selected genes/events use a separate
portable target workflow.

Katana is an independent Nextflow DSL2 workflow reusing nf-core modules. It is
maintained separately from official nf-core releases. The sections below follow the
input, contrasts, references, execution, profiles and resource organization of
the published nf-core guides, while documenting Katana's actual implementation.
Options from another pipeline are not automatically Katana options.

- [Quick start](QUICKSTART.md)
- [Complete parameter reference](NEXTFLOW_PARAMETERS.md)
- [Installation and release assets](RELEASE_INSTALL.md)
- [Accepted versions and tested boundaries](validation/entrypoints-runtime-20261025.md)

## Direct command-line example (no parameter YAML required)

The current source workflow generates project/settings JSON automatically. See the
[direct parameter guide](NEXTFLOW_SIMPLE.md) and [full YAML template](../examples/nextflow/params.complete.yaml).
The 28 October release includes this interface; older 25 October wheels use
the explicit configuration examples later in this guide.

`-params-file` is an optional, standard Nextflow way to load YAML/JSON values.
`study.yml` is not a required Katana file or a reserved filename. All native
parameters can instead be supplied explicitly, including `--outdir`.

From a checkout with a **matching installed package/container**, the complete
FASTQ or BAM command has this form. Replace the example paths, reference IDs,
allocation and images with real values; this is not a ready-to-run study command.
The sample sheet selects FASTQ versus BAM input.

```bash
nextflow run pipeline/main.nf \
  -profile slurm,apptainer \
  -c examples/nextflow/kth.config \
  --slurm_account YOUR_ALLOCATION \
  --mode complete \
  --input /data/project/samples.csv \
  --contrasts /data/project/contrasts.csv \
  --fasta /data/reference/genome.fa \
  --gtf /data/reference/annotation.gtf \
  --reference_build YOUR_BUILD \
  --annotation_release YOUR_RELEASE \
  --species 'Mus musculus' \
  --gene_sets /data/resources/pathways.gmt \
  --rbp_gene_sets /data/resources/rbp_gene_sets.gmt \
  --motif_database /data/resources/motifs.json \
  --clip_peaks /data/resources/matched_clip_peaks.tsv \
  --katana_container /data/images/katana-matching-version.sif \
  --dexseq_tools_container /data/images/dexseq-tools.sif \
  --outdir /data/project/results \
  -work-dir /data/project/work \
  -resume
```

The shell line continuations form one command. `--outdir` is the final output
root, `-work-dir` is the Nextflow work/cache directory, `-c` loads a Groovy site
configuration, and `-profile` selects runtime profiles. `-resume` reuses retained
matching work when present; it does not recover deleted work.

[The supplied kth.config](../examples/nextflow/kth.config) is a Katana Dardel
example with shared allocation and the Klemming bind. Load the site's runtime
modules before launch. Other clusters should replace it with their own site
configuration. The name `kth.config` is arbitrary, not a Nextflow option.

The official nf-core institutional profile is named
[`pdc_kth`](https://nf-co.re/configs/pdc_kth/), not `kth_config`. This Katana
checkout does not automatically import that institutional profile. Its current
`--project` option names a downstream project file, whereas the official profile
uses `--project` for the Slurm account; do not apply it unchanged. This example
uses `--slurm_account` instead and does not automatically route jobs to memory
nodes. Queue/account configuration is separate from biological parameters.

`project-template.json` and `flow.json` are now generated internally by the current
source workflow. The explicit JSON examples below remain available for advanced settings
and older releases; they are not additional requirements of the direct route.

## 1. Choose an entry

| Available input | Entry | What executes |
| --- | --- | --- |
| Bulk FASTQ | `--mode complete` | Alignment, quantification and native differential analyses, Katana downstream, seven-section report. |
| Bulk genome-aligned BAM | `--mode complete` with a BAM samplesheet | BAM/reference audits, quantification and differential analyses, Katana downstream and report. Existing alignments are retained. |
| FASTQ/BAM, requiring only upstream files | `--mode upstream-only` | Native input bundle; no Katana downstream or study report. |
| Verified Katana upstream bundle | `katana-flow import`, then `--mode downstream-only` | Downstream computations and report; no upstream read-processing rerun. |
| Existing Katana v2 project with native input assets | `--mode downstream-only` | New selected downstream computations, with readiness checks. |
| Completed Katana/nf-core-derived result files | Explicit project/source mapping, then reporting or downstream planning | Reuses mapped completed evidence; missing computations require adequate inputs and a new explicit plan. The mapping step connects each result to its samples and comparison. |
| Completed flow collection | `katana-flow report` | Strict report assembly; no statistical fits. |
| Completed report sources or portable report | `katana-report` | Read-only assembly/refresh, with honest missing-result states. |
| Selected genes/exact events | `katana-targets` or browser gene/event selection | Independent cards, requested Sashimi and portable target reports. |

The raw workflow accepts **bulk RNA-seq**. Declaring `scRNA`/`snRNA` does not
enable a raw single-cell workflow. Processed single-cell/nucleus inputs have
their own supported downstream adapters and biological-unit requirements.

## 2. Installation and directory layout

Install/load the matching release artifacts using [RELEASE_INSTALL.md](RELEASE_INSTALL.md).
The Nextflow execution host needs Nextflow **25.10.4 or newer**, compatible
Java and the selected execution engine. Native complete execution has been
accepted on Linux x86-64; other platforms and cluster profiles retain their
[documented qualification boundaries](validation/entrypoints-runtime-20261025.md).

An installed Katana wheel exports the workflow:

```sh
katana-flow export-pipeline --outdir tools/katana-pipeline
nextflow run tools/katana-pipeline/main.nf --help
```

An existing export directory is not overwritten. A source checkout can instead
use `pipeline/main.nf`. Record the checkout commit and runtime image identities;
a workflow exported by a wheel retains the source bundled in that wheel.

**Source support for `strandedness=auto`:** use the current `main` checkout's
`pipeline/main.nf`. The frozen
`v2.0.0-beta.5.dev20261025` wheel and exported workflow retain their original
explicit-strandedness requirement; they are not updated by these source changes.
The released Katana downstream container remains usable; the inference processes
use their own pinned RSeQC/gffread images. Obtain the current source workflow with:

```sh
git clone --branch main https://github.com/Taratype0/Katana-splicing-tool.git katana-source
nextflow run katana-source/pipeline/main.nf --help
```

For `auto`, replace the exported workflow path in the commands below with
`katana-source/pipeline/main.nf`, retaining the other required study inputs
and matching runtime configuration.

Example study layout:

```text
study/
  inputs/
    samples.csv
    contrasts.csv
    project-template.json
    flow.json
  reads/
  references/
  databases/
  tools/katana-pipeline/
  results/                    # created by execution
  work/                       # Nextflow cache products
```

Commands below run from `study/`. Keep project/configuration inputs outside
`results/` and `work/`: a project input directory must not contain those outputs.
All input files must be local or mounted on the execution machine; raw-sheet
HTTP/S3 URLs must first be staged explicitly. No workstation or cluster path
is built into these examples.

## 3. Samplesheet input

Use CSV, or a `.tsv` file with tab-separated fields and an exact header. One
row represents **one independent biological sample**, using unique IDs.
Copyable [FASTQ](../examples/nextflow/samples_fastq.csv),
[BAM](../examples/nextflow/samples_bam.csv),
[simple contrast](../examples/nextflow/contrasts.csv), and
[selected factorial contrast](../examples/nextflow/contrasts_factorial.csv)
templates are provided with [their instructions](../examples/nextflow/README.md).

### FASTQ

```csv
sample,condition,fastq_1,fastq_2,strandedness,read_length,subject
Control_1,Control,../reads/Control_1_R1.fastq.gz,../reads/Control_1_R2.fastq.gz,reverse,150,animal_01
Control_2,Control,../reads/Control_2_R1.fastq.gz,../reads/Control_2_R2.fastq.gz,reverse,150,animal_02
Control_3,Control,../reads/Control_3_R1.fastq.gz,../reads/Control_3_R2.fastq.gz,reverse,150,animal_03
Treated_1,Treated,../reads/Treated_1_R1.fastq.gz,../reads/Treated_1_R2.fastq.gz,reverse,150,animal_04
Treated_2,Treated,../reads/Treated_2_R1.fastq.gz,../reads/Treated_2_R2.fastq.gz,reverse,150,animal_05
Treated_3,Treated,../reads/Treated_3_R1.fastq.gz,../reads/Treated_3_R2.fastq.gz,reverse,150,animal_06
```

| Column | Requirement |
| --- | --- |
| `sample` | Required, unique; starts with a letter, followed by letters/digits/`_`/`.`/`-`. Use stable IDs, not display labels or paths. |
| `condition` | Required; starts with a letter, then letters/digits/underscores. Matches contrast levels exactly, including case. |
| `fastq_1` | Required for a FASTQ row. Relative paths resolve from the samplesheet directory. |
| `fastq_2` | Mate file for paired reads; blank for single-end. FASTQ pairing determines `single_end`. |
| `strandedness` | Required: `forward`, `reverse`, `unstranded`, or `auto` in the current `main` source workflow. Explicit values are preserved. |
| `read_length` | Positive nominal read length, or supply study-wide `--read_length`. All samples in this route must declare the same nominal length. |
| `subject` | Optional true biological-unit identity. If specified in a selected comparison, specify every participating subject. Do not reuse an identity for independent animals. |
| `source` | Optional; `fastq` for FASTQ rows, `genome_bam` for BAM rows. Empty/`auto` is resolved from files. |
| `assay` | Optional; current raw route accepts only `bulk_rna_seq`. |
| `single_end` | Optional FASTQ consistency check; required for BAM. Explicit boolean text. |
| `batch` | Optional recorded metadata. It **does not** enable batch-adjusted native models in the complete route. |

Each physical input path can belong to only one row. FASTQ inputs are staged as
canonical `.fastq.gz` symlinks for the retained upstream tools; the source must
already be gzip compressed. The staging probe checks only the first 100 records,
not whole-file integrity or every read-pair identity. Native complete
acceptance currently covers paired-end data; parser acceptance of single-end
metadata is not a claim of equally complete single-end qualification.

### Multiple lanes and technical replicates

For multiple lanes, repeat the biological sample's row and add a distinct
`lane` value for each row. Supply the original gzip FASTQs directly:

```csv
sample,lane,condition,fastq_1,fastq_2,strandedness,read_length,subject
Control_1,L1,Control,/data/c1_L1_R1.fq.gz,/data/c1_L1_R2.fq.gz,auto,150,mouse01
Control_1,L2,Control,/data/c1_L2_R1.fq.gz,/data/c1_L2_R2.fq.gz,auto,150,mouse01
```

Nextflow groups lanes into one biological sample and runs one `STAGE_FASTQ`
task per sample. That task concatenates gzip members in sorted lane order,
using exactly the same order for both mates, before QC and alignment. Single
lanes remain symlinks. No separate lane-preparation command is needed.
All biological metadata and library settings must agree across lanes; reused
files and repeated/missing lane identifiers are rejected. Source files are
unchanged. Lane paths and ordering are recorded in the staging audit. Its
first-100-record probe is not a full-file integrity or read-pair audit.
A lane must never become a new biological replicate. See the
[sample-identity examples](../examples/contracts/sample-identity-v1/README.md).

### Genome-aligned BAM

```csv
sample,condition,genome_bam,strandedness,single_end,read_length,subject
Control_1,Control,../alignments/Control_1.bam,reverse,false,150,animal_01
Control_2,Control,../alignments/Control_2.bam,reverse,false,150,animal_02
Control_3,Control,../alignments/Control_3.bam,reverse,false,150,animal_03
Treated_1,Treated,../alignments/Treated_1.bam,reverse,false,150,animal_04
Treated_2,Treated,../alignments/Treated_2.bam,reverse,false,150,animal_05
Treated_3,Treated,../alignments/Treated_3.bam,reverse,false,150,animal_06
```

Provide a genome-coordinate `.bam`, not a transcriptome BAM or a BAM with
removed sequences. The workflow audits header/reference contig compatibility,
available sequence MD5 values, primary-record sequences/CIGAR, mate flags and
read lengths before quantification, then creates its own sorted BAM and index.
Supplying the original index is not required by this raw sheet. The BAM audit
does not establish strandedness; `auto` adds the separate
inference step below. Upstream read-removal history cannot be reconstructed.

For variable observed lengths, `--allow_variable_read_length true` permits
positive primary-record lengths up to the declared nominal length. It does not
permit different nominal read lengths across the study. Within a selected
comparison, layout and strandedness must be consistent.

BAM QC describes the **supplied BAM records**. It cannot recover discarded
reads, original raw FASTQ QC or the original trimming/alignment history.

### Automatic strandedness in the source workflow

Set `strandedness` to `auto` for each FASTQ or genome-BAM sample whose library
direction is unknown. Declared `forward`/`reverse`/`unstranded` rows can coexist
with `auto` rows and keep their declared values. Only `auto` samples incur
inference.

FASTQ samples first obtain their STAR genome BAM; supplied BAMs are not
realigned. After the shared BAM/reference audit, sorting and indexing,
gffread 0.12.1 converts the matched **unfiltered** GTF to BED12. RSeQC 5.0.3
`infer_experiment.py` compares the alignments with these strand models. This
uses no preliminary Salmon quantification or additional alignment. It supports
the native paired- and single-end inference formats; complete-route runtime and
assay qualification still have their documented scope.

The defaults are `--strandedness_sample_size 200000` and
`--strandedness_mapq 30`, matching RSeQC. Sampling takes the first eligible
gene-overlapping alignment records in BAM order, rather than a random sample or
200,000 paired fragments; the tool may read more records to find them.

The sense and antisense fractions are normalized over their sum, excluding the
separately reported failed-to-determine fraction. A normalized fraction of at
least `--stranded_threshold 0.8` calls `forward` or `reverse`; otherwise an
absolute difference of at most `--unstranded_threshold 0.1` calls `unstranded`.
Other outcomes are **undetermined and stop execution**, with no unstranded
fallback. Missing/malformed evidence, no informative alignments, or an inferred
layout inconsistent with `single_end` also fail. All selected samples must be
resolved and every selected comparison must have consistent actual
strandedness before Salmon, DEXSeq or rMATS can start.

The requested and resolved values, native RSeQC text/stderr evidence, parameters
and inference audit are retained in pipeline QC outputs and bundle provenance.
Failures keep their diagnostic evidence. Inference describes the supplied BAM;
upstream strand filtering can bias it and cannot be reversed by this step.
See the [four parameter definitions](NEXTFLOW_PARAMETERS.md#3-raw-read-qc-alignment-and-quantification)
and [pinned module/tool provenance](../pipeline/modules/local/infer_strandedness.provenance.json).

## 4. Contrastsheet input and direction

Create only the comparisons you want. The workflow does not enumerate all
pairs of conditions.

```csv
contrast,treatment,control
Treated_vs_Control,Treated,Control
```

| Column | Meaning and validation |
| --- | --- |
| `contrast` | Unique output identifier; same safe-ID syntax as `sample`. |
| `treatment` | Target/numerator level from `samples.csv:condition`. |
| `control` | Reference/denominator level; must differ from treatment. |
| `description` | Optional human-readable explanation; no statistical effect. |
| `paired` | Optional `true`/`false` (also 1/0, yes/no); true selects matched-subject analysis. |
| `blocking` | Optional `subject` for a paired comparison; otherwise leave blank. |
| `design` / `formula` | Optional `~ condition` or `~ subject + condition`; must agree with `paired`. |

Use `contrast,treatment,control` as the first three columns for compatibility.
Unknown extra columns, duplicate contrast IDs, repeated directional treatment/
control pairs, missing levels and self-comparisons fail validation. This is not
the `id,variable,reference,target,blocking` sheet used by
`nf-core/differentialabundance`; map those fields explicitly before reuse.

For example, a study with four backgrounds and two treatments can request the
following **twelve selected questions**, without calculating every pair. These
are generic condition labels; include only rows relevant to your actual study:

```csv
contrast,treatment,control
KO_vs_CONTROL_NO_DP,KO_NO_DP,CONTROL_NO_DP
WTrescue_vs_KO_NO_DP,WTrescue_NO_DP,KO_NO_DP
MUTrescue_vs_KO_NO_DP,MUTrescue_NO_DP,KO_NO_DP
MUTrescue_vs_WTrescue_NO_DP,MUTrescue_NO_DP,WTrescue_NO_DP
KO_vs_CONTROL_DP,KO_DP,CONTROL_DP
WTrescue_vs_KO_DP,WTrescue_DP,KO_DP
MUTrescue_vs_KO_DP,MUTrescue_DP,KO_DP
MUTrescue_vs_WTrescue_DP,MUTrescue_DP,WTrescue_DP
CONTROL_DP_response,CONTROL_DP,CONTROL_NO_DP
KO_DP_response,KO_DP,KO_NO_DP
WTrescue_DP_response,WTrescue_DP,WTrescue_NO_DP
MUTrescue_DP_response,MUTrescue_DP,MUTrescue_NO_DP
```

Within-treatment rows quantify background differences; within-background rows
quantify the declared treatment response. Labels such as `WTrescue` do not make
the software infer restoration, causality or experimental success. A difference
between two response estimates is not automatically a fitted interaction.

At least **two independent biological samples per selected arm** are required
for upstream differential testing. The default complete downstream network
defaults to **six independent samples in total per selected comparison**, not
six per group. Set `network.parameters.minimum_network_samples` explicitly
(for example, `4` for two-plus-two) to use a different network minimum. Values
must be integers of at least three. Other analyses keep their own requirements;
this option is an eligibility setting, not a power estimate. See [small studies](QUICKSTART.md#small-studies).

The default model is `~ condition`. For matched samples, give each animal the
same `subject` value in both conditions, and set `paired=true` in the contrast:

```csv
contrast,treatment,control,paired
Treated_vs_Control,Treated,Control,true
```

Every selected subject must occur exactly once in each arm. The workflow aligns
both arms by subject identity, including rMATS BAM lists and SUPPA columns;
CSV row order does not establish pairing. Independent animals need different
subject IDs. Batch-adjusted and interaction formulas remain separate supported
analyses, not implicit extensions of this raw-input model.

DEG uses `~ subject + condition`; rMATS enables native `--paired-stats`;
SUPPA enables native `-pa`. DEXSeq exon usage and transcript usage include
`subject:exon` in both likelihood-ratio models, testing `condition:exon`.
DEXSeq's native fold-change summary remains a descriptive condition estimate;
its P/FDR are from the subject-adjusted test. A paired network retains the
observed-PSI discovery/kME algorithm and uses subject-blocked module association.
Its record distinguishes observations from independent subjects (for example,
six observations from three complete pairs). Paired module association needs
at least three complete pairs. See [paired implementation and tests](validation/paired-entry-20261028.md).

Effect directions retain the native methods:

- DESeq2: `log2(treatment/control)`.
- rMATS: treatment is group 1; delta PSI is group 1 minus group 2.
- SUPPA: its adapter supplies control first and treatment second.
- DEXSeq: the actual native fold-change column name remains authoritative;
  current exports record its control-over-treatment exception. Gene/stage-wise
  adjusted values do not carry an effect direction by themselves.

Do not flip a native statistic merely by renaming a comparison. An imported
table needs its actual fit direction, sample-vector order and testing scope.

## 5. References and downstream resources

Always supply a matched genomic FASTA and GTF. `--reference_build` and
`--annotation_release` record your declarations; `--genome` is a fallback build
label here, **not** an automatic iGenomes download selector. The workflow checks
contigs, coordinates and transcript/exon identities and records hashes.

Mitochondrial-encoded genes and annotated pseudogenes are excluded from the
analysis GTF **before native quantification/differential analysis**. The full
input GTF and exclusion audit are retained. Nuclear genes are not excluded
merely because a gene symbol starts with `MT`. Actual `gene_type` or
`gene_biotype` annotations are required for pseudogene filtering.

`project-template.json` supplies resource identities and presentation labels;
sample and contrast identities come from the genuine upstream bundle:

```json
{
  "name": "Example bulk study",
  "rmats_mode": "JCEC",
  "metadata": {"species": "Mus musculus", "annotation_release": "your-release"},
  "contrast_labels": {"Treated_vs_Control": "Treated versus control"},
  "sample_attributes": {},
  "assets": []
}
```

The JSON above shows the template structure, **not a resource-complete report
configuration**. Populate `assets` with the following real files:

| Asset kind | Use |
| --- | --- |
| `gene_sets_gmt` | Declared GO/KEGG or other functional gene sets, preserving their actual source/namespace. |
| `rbp_gene_sets` | RBP identity sets for DEG-based enrichment and network annotation. |
| `rbp_motif_database` | Species-matched motif models. |
| `fimo_runtime` | Verified FIMO worker/runtime descriptor; an executable alone is not this descriptor. |
| `rbp_binding_peaks` | Species/build-matched CLIP peaks and their experimental provenance. |

Each asset has `asset_id`, `kind`, `path`, optional `required` and `metadata`.
Paths resolve from the template. Versioned biological evidence requires
`metadata.source_url`, `version`, `sha256`, `species`, and matching
`genome_build` for coordinate-based evidence. Supply the actual expected
resource schema, not a renamed unrelated file. Resource readiness fails before
read processing if required files/runtime are missing. Conda/container software
does not supply missing experimental databases.

The selected JC/JCEC mode is explicit and recorded. The native bundle retains
both modes; the imported analysis project selects the requested one. Enrichment
backgrounds are the detected/testable genes for that evidence layer, intersected
with the mapped resource, not the entire reference genome. Foreground,
background and hit-gene tables are included in outputs.

## 6. Downstream analysis configuration

Save `inputs/flow.json`:

```json
{
  "schema": "katana.flow.config/v1",
  "contrast_ids": ["Treated_vs_Control"],
  "event_set": "stable_significant",
  "parameters": {
    "differential.deg": {"fdr": 0.05, "min_abs_log2fc": 1.0}
  },
  "analyses": {"coverage": false},
  "network": {
    "scope": "contrast",
    "bootstrap": {"enabled": false},
    "parameters": {"rmats_mode": "JCEC", "min_module_size": 20, "merge_cut_height": 0.25}
  },
  "design_factors": ["condition"],
  "interactions": [],
  "mapp": {"adapter": "mapp.maei", "inputs": {}, "parameters": {}}
}
```

This selects RS-DAS by default. The fixed per-sample read-support gate is
`IJC + SJC >= 10` in **every expected sample**. Ordinary DAS can be chosen with
`event_set: "significant"`; motif/CLIP tasks still explicitly require RS-DAS.
The report and supported downstream selections retain the same scope.
See [thresholds and their basis](THRESHOLDS.md).

Each selected comparison gets independent network discovery/followups and
MotifTrace tasks, using only its own samples. Bootstrap is off; Hub ranking uses
`|kME|`. MotifTrace analyzes the supplied/requested models without fixed
model-count, test-count or matrix-element caps, while the report shows up to
ten significant models and keeps complete families in outputs. No significant
model produces a concise non-significance notice, not a fabricated sequence
result.

MAPP uses the packaged hash-bound author implementation and prepares the
14-window MAEI family from actual transcript quantification and matching
annotation. Empty `mapp.inputs` is therefore permitted; absent necessary
biological data is not. It does not run all other MAPP subworkflows.

An explicit multifactor downstream project requires compatible interaction
design entries; independent contrast differences cannot replace a fitted
interaction. That downstream extension does not enable an unsupported upstream
`--design`. The [parameter guide](NEXTFLOW_PARAMETERS.md) separates these layers.

## 7. Validate before scheduling

Metadata and file existence checks, using the installed/exported script:

```sh
python tools/katana-pipeline/bin/katana_validate_inputs.py \
  --samplesheet inputs/samples.csv --contrasts inputs/contrasts.csv \
  --fasta references/genome.fa --gtf references/annotation.gtf \
  --outdir checks/raw-inputs

katana-flow preflight --project-template inputs/project-template.json \
  --config inputs/flow.json --out checks/template-preflight.json
```

The Python script's `--skip-file-check` is metadata-only and cannot certify BAM
content or reference compatibility. Template preflight checks required resources
before upstream execution; the genuine imported project is checked again after
upstream export. Neither preflight command executes scientific fits.

## 8. Run the complete workflow

Set `KATANA_IMAGE` and `DEXSEQ_IMAGE` to the matching loaded/built images from
the [release guide](RELEASE_INSTALL.md) and
[DEXSeq instructions](../containers/dexseq-tools/README.md).
These names are local variables, not claims of a public registry image.

```sh
nextflow run tools/katana-pipeline/main.nf -profile docker \
  --mode complete --input inputs/samples.csv --contrasts inputs/contrasts.csv \
  --fasta references/genome.fa --gtf references/annotation.gtf \
  --reference_build YOUR_BUILD --annotation_release YOUR_RELEASE \
  --project_template inputs/project-template.json --flow_config inputs/flow.json \
  --katana_container "$KATANA_IMAGE" --dexseq_tools_container "$DEXSEQ_IMAGE" \
  --outdir results -work-dir work
```

The BAM route uses the **same command with a BAM samplesheet**. A BAM-only run
does not build an unused STAR index or realign those BAMs. It still performs
audits, Salmon genome projection and the specified downstream calculations.

For upstream-only execution:

```sh
nextflow run tools/katana-pipeline/main.nf -profile docker \
  --mode upstream-only --input inputs/samples.csv --contrasts inputs/contrasts.csv \
  --fasta references/genome.fa --gtf references/annotation.gtf \
  --dexseq_tools_container "$DEXSEQ_IMAGE" --outdir upstream-results -work-dir upstream-work
```

It stops at `upstream-results/katana_bundle`, not a complete report.

## 9. Existing upstream/nf-core results

A **complete Katana native bundle** is directly importable:

```sh
katana-flow import --bundle upstream-results/katana_bundle \
  --project-template inputs/project-template.json --outdir imported-project

katana-flow preflight --project imported-project/project.json --config inputs/flow.json
katana-flow plan --project imported-project/project.json --config inputs/flow.json \
  --require-ready --out plans/downstream-plan.json

nextflow run tools/katana-pipeline/main.nf -profile docker \
  --mode downstream-only --project imported-project/project.json \
  --flow_config inputs/flow.json --katana_container "$KATANA_IMAGE" \
  --outdir downstream-results -work-dir downstream-work
```

The separate `plan` command is an inspection artifact; the command above plans
from `flow.json` inside Nextflow. To reuse an **already executed** sealed plan,
use the original path and the resume instructions below.

For external nf-core outputs, map files into an explicit compatible Katana v2
project/source manifest first. There is currently no universal command that
recognizes every arbitrary nf-core release folder:

| Upstream producer | Map the actual supplied evidence | What that producer alone does not establish |
| --- | --- | --- |
| `nf-core/rnaseq` | Genome BAMs/indexes, Salmon `quant.sf`, expression/QC tables, matched references and sample identities | rMATS/DAS, complete DTU/DEU, network, MotifTrace and all report computations. |
| `nf-core/rnasplice` | rMATS five event-class tables, JC/JCEC, native sample order, exon/transcript/other native results actually run | Katana-specific network and MotifTrace or additional missing evidence layers. |
| `nf-core/differentialabundance` | Actual differential result/matrix files, original model, direction, method and testing universe | Other methods or biological units absent from that run. |

The name `rnadifferential` is not a Katana mode; use the actual producer/version
in provenance. Do not relabel a gene-level q-value as an exon-level statistic,
infer rMATS vector order from file sorting, or overwrite native effects/FDR.
Read [EXISTING_RESULTS.md](EXISTING_RESULTS.md) for completed-source descriptors,
SHA-bound adoption, study identity checks, path mapping and computing only
missing work. Unsupported/missing native schemas need an explicit adapter;
report rendering is not that adapter or the missing analysis.

## 10. Report-only execution and output files

Strict reporting from a completed flow collection:

```sh
katana-flow report --collection completed-collection --outdir report-final
```

Read-only assembly from explicitly mapped completed source results:

```sh
katana-report --project results.json --outdir report-from-results \
  --compact --archive downloads/study-report.zip --result-store payload-store
```

`results.json` is the completed-result descriptor documented in
[EXISTING_RESULTS.md](EXISTING_RESULTS.md), not a CSV or an arbitrary directory.
Use `--no-archive` instead of `--archive` when no extra ZIP is wanted. See
[REPORTS.md](REPORTS.md) for refreshing an already portable report.

The final report retains the approved seven English sections, comparison
selector, shared study QC, RS-DAS/ordinary-DAS switch, Top 10, network graphs and
complete member/gene tables, significant motif/enrichment summaries and complete
data links. Full tables are paginated/searchable or downloadable. Target reports
remain independent; the fixed study format is not redesigned.

Important outputs:

| Output | Meaning |
| --- | --- |
| `results/katana_bundle/manifest.json` | Upstream native artifact ledger and file hashes. |
| `results/pipeline_info/` | Validated sheets/audits plus Nextflow trace, timeline, execution report and DAG. The Nextflow execution report is not the scientific Katana report. |
| Pipeline QC strandedness evidence and bundle provenance | `auto` requests, resolved values, native RSeQC text/stderr and inference decisions; declared manual settings remain explicit. |
| `results/bam/`, `quant/`, `dexseq_exon/`, `dtu/`, `DEG/` | Published upstream products. |
| `results/downstream/` | Isolated downstream task products and collection, according to published process filenames. |
| Final `REPORT.html`, `REPORT.zip` | Study report and portable download. Use the paths emitted by the actual command. |
| `flow.acceptance.json` | Strict completed-report acceptance, including verified calculation-task count and archive hash. |

A failed completion gate retains diagnostics and may name its ZIP
`REPORT.diagnostic.zip`; it does not issue successful `flow.acceptance.json`.
An empty Slurm queue or an HTML file alone is not completion evidence.

## 11. Deep dive: chosen genes/events, local or server

Select rows in the browser, choose **Add genes** or **Add events**, then export
the target request ZIP. For scripted selection:

```json
[
  {
    "source_artifact_id": "ACTUAL_ARTIFACT_ID",
    "row_ids": ["ACTUAL_SOURCE_ROW_ID"],
    "target_kind": "event",
    "filters": {"event_set": "stable_significant", "fdr": 0.05}
  }
]
```

Save it as `selections.json`, substituting real exported IDs. Then:

```sh
katana-targets export --project imported-project/project.json \
  --selections selections.json --outdir target-request --contrast-id Treated_vs_Control
katana-targets inspect --request target-request/TARGET_REQUEST.zip
katana-targets plan --request target-request/TARGET_REQUEST.zip --bindings bindings.json
katana-targets run --request target-request/TARGET_REQUEST.zip \
  --bindings bindings.json --outdir target-results
katana-targets results --outdir target-results
```

Use the actual archive path printed by `export`. `target_kind` can be `gene`
or `event`. Use `--module gene_card --module event_card` on export if no read
analysis is requested. A gene target includes all available measured events
for that gene; it is not silently restricted to significant events.

`bindings.json` maps exported resource IDs to accessible BAM/index/reference
files. Read IDs from `inspect`; example IDs below are placeholders:

```json
{
  "schema": "katana.target_bindings/v1",
  "resources": {
    "RESOURCE_ID_FOR_BAM": "alignments/sample.bam",
    "RESOURCE_ID_FOR_INDEX": "alignments/sample.bam.bai",
    "RESOURCE_ID_FOR_GTF": "references/annotation.gtf"
  }
}
```

Paths resolve from this binding file and must match the request's verification
records. BAMs remain where the computation runs; they are not bundled into the
portable request/report. Cards can reuse saved evidence without BAM access;
new Sashimi requires usable bound reads and indexes.

The output contains `INDEX.html`, `REPORT.zip`, per-target `REPORT.html`, cards,
sample plots, saved network/motif evidence, native Sashimi when requested and
available, complete TSV/Parquet/JSON and a hash-covered result manifest.
Deep dive reuses differential/network/motif results; it does not refit them.

Transfer the completed ZIP back to another machine, inspect offline, import
through the matching Katana browser or restore a resumable workspace:

```sh
katana-targets restore --results completed-targets.zip --outdir restored-targets
katana-targets results --outdir restored-targets
katana-open --no-open
```

An offline report can display evidence but cannot execute Python or read a
remote BAM. `katana-open` supplies the browser service; imported results do not
submit new work. Use a secure tunnel for a remote loopback service. See
[TARGET_COMMANDS.md](TARGET_COMMANDS.md) for full resource/cache contracts.

## 12. Profiles, Slurm and other execution methods

### Container profiles

- `-profile docker`: Linux Docker execution using caller UID/GID.
- `-profile slurm,apptainer` or `slurm,singularity`: native Slurm jobs using
  compute-node-visible SIF/images and paths; site qualification still applies.
- `-profile local`: local native environment already installed.
- `-profile conda`: enables upstream modules' declared Conda recipes; this does
  **not** by itself create the complete Katana Python/R/native environment.
  Use the locked [standalone environment](RUNTIME_ENVIRONMENTS.md).

For a Slurm site:

```sh
nextflow run tools/katana-pipeline/main.nf -profile slurm,apptainer \
  -params-file inputs/params.yaml -c site.config -work-dir shared-work
```

Put actual scientific/file settings in `params.yaml`. Set the site's account
and partition through `--slurm_account`/`--slurm_partition` or the params file.
Use `site.config` for process resources/executor settings, for example:

```groovy
process {
    withName: '.*:STAR_ALIGN' {
        cpus = 12
        memory = '72 GB'
        time = '16 h'
    }
}
```

Check actual fully qualified process names in the trace. A mismatched selector
does not change a resource request. `--max_cpus`, `--max_memory`, `--max_time`
are optional **per-job ceilings**, not global quotas or concurrency controls.
The Slurm profile sets `executor.queueSize = 0`; there is no application-imposed
four-task cap. Slurm still decides admission under account/site limits.
Independent per-comparison network/MotifTrace/MAPP tasks obey the actual DAG;
dependencies are not bypassed for parallelism.

### Standalone analysis, browser and Python

These alternatives use the same installed analysis code and evidence contracts:

```sh
# Discover and run one downstream workflow in a prepared project.
python -m katana_next.cli --project imported-project/project.json workflows
python -m katana_next.cli --project imported-project/project.json plan differential.deg \
  --parameters deg-parameters.json
python -m katana_next.cli --project imported-project/project.json run differential.deg \
  --parameters deg-parameters.json

# Open the browser service locally, or leave browser opening to the user.
katana-open
katana-open --no-open
```

The installed Python API is documented in [USAGE.md](USAGE.md) and
[TARGET_COMMANDS.md](TARGET_COMMANDS.md). A module run does not automatically
replace a complete-report flow. For target jobs, `katana-targets tasks --request
targets.zip --format tsv` lists real target IDs; each Slurm/PBS/Kubernetes job
can call `katana-targets run --target-id ACTUAL_TARGET_ID` with its own output
directory. Katana does not automatically submit to every cluster scheduler.

## 13. Resume, updates and storage

New runs use the single Nextflow entry above. Configuration, upstream models,
Katana analysis, report generation, project retention and delivery verification
are tasks in the same DAG. See [output directories](NATIVE_OUTPUTS.md).
Add `--clean_work` for native cleanup after successful delivery. Nextflow keeps
work referenced by other resumed runs; `nextflow clean <run-name> -n` previews
cleanup for a completed run. Failed work remains available for recovery.
`--storage_mode lean` reduces optional publication, independently of cleanup.

The earlier `katana-flow run` launcher remains for historical runs only. Its
[managed storage guide](STORAGE_MANAGEMENT.md) describes that compatibility route.
Keep an active historical run on its original software and paths to preserve
cache identity; do not replace it just to change output directory names.

Repeat an unchanged Nextflow command from the original launch directory with
`-resume`; retain `.nextflow/`, original work files and input identities.
An unchanged completed process is reused only if Nextflow's fingerprint and
files still match. Images, code, paths, references or parameter changes may
invalidate relevant caches. Do not modify cache records to force reuse.

For an already executed downstream plan:

```sh
katana-flow validate-plan --plan original-plan/plan.json \
  --project imported-project/project.json --out checks/resume-validation.json

nextflow run tools/katana-pipeline/main.nf -profile docker \
  --mode downstream-only --project imported-project/project.json \
  --flow_plan original-plan/plan.json --katana_container "$KATANA_IMAGE" \
  --outdir downstream-results -work-dir downstream-work -resume
```

Use the **exact original absolute plan path** from the real run. The relative
example is valid only when it resolves to that same path. `--flow_plan` bypasses
replanning; supplying `--flow_config` as well does not apply new settings.
For an approved presentation-only refresh of a compatible plan, change
`--flow_report_revision LABEL`; statistical task identities remain unchanged.

Archive and verify final products before removing old run directories. Deleting
`work/` prevents the corresponding Nextflow resume. GitHub source does not
replace raw data, BAMs needed for new Sashimi, biological metadata, reference
identities or completed scientific results. Portable reports retain complete
included output data but are not substitutes for every execution input.

## Documentation sources

Structure and terminology were checked against official
[nf-core/rnaseq usage](https://nf-co.re/rnaseq/usage),
[nf-core/rnasplice usage](https://nf-co.re/rnasplice/usage),
[rnasplice parameters](https://nf-co.re/rnasplice/parameters), and
[nf-core/differentialabundance usage](https://nf-co.re/differentialabundance/usage)
on 25 October 2026. Katana formats/defaults above are defined by this release's
own validator, workflow, importer, planner and CLI, not copied from those tools'
current options. See [NEXTFLOW_PARAMETERS.md](NEXTFLOW_PARAMETERS.md).

## Selecting analyses and recovering a run

The default runs the complete report. See [analysis and execution options](ANALYSIS_AND_EXECUTION_OPTIONS.md) for per-module opt-outs, resource-only resume, shared MAPP preparation and storage behaviour.
