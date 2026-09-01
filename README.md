# Katana RNA-seq + Splicing Candidate Screening Tool

Katana is a desktop GUI workflow for **comparison-first**, **direction-aware** RNA-seq and splicing follow-up analysis.

It is designed to work across projects, not only one specific dataset. Each project is scanned locally, each comparison keeps its own direction/configuration, and each analysis module runs independently instead of chaining downstream modules automatically.

For methods development, CoSplice provides a same-input head-to-head benchmark of its single
confidence-annotated WGCNA production partition against exact R WGCNA and classical baselines. The
benchmark evaluates valid partition recovery, matched-bootstrap stability, known-truth recovery,
independent-cohort preservation, locked module-score projection, paired phenotype association,
resistance to single-gene pseudo-modules, GO/KEGG annotation and cross-cohort gene/function
replication. Biological annotations validate and explain a module after discovery; they never define
modules or tune network parameters. A publication aggregation command builds a traceable
multi-dataset report with source fingerprints and explicit claim limits. An untouched fourth
TCGA SpliceSeq cohort is supported as a confirmatory target without validation-set event ranking or
parameter tuning; see
[`docs/COSPLICE_NETWORK.md`](docs/COSPLICE_NETWORK.md) and the
[`publication benchmark protocol`](docs/PUBLICATION_BENCHMARK.md).
Completed BAM-junction and rMATS caches can additionally be compared with
`katana-cosplice-input-benchmark`. It audits stable-feature reproducibility, sample structure,
leave-one-replicate-out condition effects, significant-SE junction support, and BAM-only discovery
without rescanning BAMs or conflating compound rMATS events with individual junctions.

## Core design principles

- **Comparison-first**: each comparison keeps its own pages, outputs, cache, and logs
- **Direction-aware**: canonical splicing and DEG directions are detected from project inputs/configuration, shown read-only in Pairing, and inherited downstream
- **Gene-centered candidate screening**: candidate ranking, selection, cards, and follow-up all start from per-comparison candidate genes
- **Independent module execution**: opening a page reads status/cache only; analysis runs only after the user clicks `Run`, `Generate`, `Apply`, or `Recalculate`
- **Project-specific annotation**: gene symbol mapping is loaded from the current project's own annotation TSV when available

---

## What Katana does

Katana currently supports:

- Project scanning and input discovery
- Pairing / comparison direction setup
- Cross-comparison ordering and visualization group setup
- Overview plots for significant splicing events
- One shared stable-AS filter using rMATS FDR/dPSI, per-replicate IJC+SJC support, and replicate direction consistency
- Transcript + splicing integration
- Per-comparison candidate gene ranking and selection
- Candidate cards + event follow-up
- Cross-comparison candidate matrix and AS dPSI pattern comparison
- Event-driven sashimi manifest generation and sashimi execution
- Jutils output browsing
- RBP support/enrichment over candidate genes
- Isoform candidate support from quant.sf / GTF / DTU
- Comparison-first DEG exploration with Volcano/MA plots
- GO:BP / GO:MF / GO:CC / KEGG enrichment from DEG, candidate-selection, or custom genes
- Standalone and embedded **Katana-CoSplice** network discovery for replicate-supported sample-level PSI programs across multi-condition cohorts

### Katana-CoSplice network discovery

`katana-cosplice` is also installed as a command-line program. It builds sample-level splice-usage matrices from either exact rMATS events or BAM-derived competing junction clusters. BAM mode extracts observed CIGAR `N` junctions, uses the current project's GTF for annotation, and keeps rMATS as an event/differential-evidence layer rather than limiting discovery to the rMATS catalogue. Both modes feed the same parameter-fingerprinted stable-AS QC, which applies condition coverage, PSI-range, read-support and replicate-consistency safeguards before network fitting.

The production network uses official R WGCNA for signed soft adjacency, exact TOM, Dynamic Tree Cut and eigengene-based module merging. With replicate IJC/SJC, CoSplice estimates beta-binomial PSI posteriors, rejects low-reliability events, and uses posterior-mean PSI for the fixed WGCNA partition. Repeated posterior draws quantify edge correlation SD, sign stability, co-assignment and fixed-module kME stability without replacing that partition. Core membership must be reproduced in two disjoint bootstrap halves. A matched Setd2 ablation rejected posterior-ensemble partitioning because its high ARI came with an over-dominant module. Event reliability is not multiplied directly into TOM adjacency because that can also collapse otherwise valid modules. Observed and covariate-residual networks, formal gained/lost/rewired tests, annotation-aware exact event alignment, and leakage-safe external preservation remain separate validation layers. PSI-only cohorts use unmodified `WGCNA::blockwiseModules`. Sparse kNN/Leiden and boundary-changing consensus remain historical benchmark modes, not competing production results.

CoSplice calculates module program scores, condition/covariate associations and gained/lost/rewired topology, then validates modules through bootstrap recovery, observed/residual preservation, RBP/motif support, DTU/isoform evidence, GO/KEGG enrichment and cross-cohort preservation. Biological annotations explain and validate modules but never define membership. Annotation-aware projection, Bayesian multi-cohort preservation and random-effects meta-analysis remain available downstream. Cached HTML/PDF reports and a module browser preserve provenance without recomputation. Every UI action and `stage --name` command runs exactly one stage; unchanged inputs return a parameter-fingerprinted cache hit. See [docs/COSPLICE_NETWORK.md](docs/COSPLICE_NETWORK.md) and [docs/COSPLICE_METHODS_3_5.md](docs/COSPLICE_METHODS_3_5.md).

**Benchmark scope clarification:** a cached final consensus is ranked only when its discovery
event universe exactly matches the benchmark universe. A full-universe consensus is never
silently subsetted into a top-N comparison; mismatched caches are reported as non-comparable.

For the exact WGCNA comparator, create/update the Conda environment and install the official R package once:

```powershell
conda env update -n katana-splicing-tool -f configs/environment.yml
conda run -n katana-splicing-tool Rscript scripts/install_wgcna.R
```

For a leakage-safe publication ablation, freeze discovery/tuning selection before the holdout is
loaded, then aggregate completed runs into one auditable panel:

```powershell
katana-cosplice-locked-ablation --help
katana-cosplice-locked-panel-report --root <completed-run-root>
```

The locked runner compares stability-only CoSplice, adaptive multi-objective CoSplice and exact R
WGCNA on the same event universe. The functional objective is activated only when predeclared
tuning-GMT enrichment passes its FDR gate; otherwise its weight is set to zero and the negative
result is retained. See
[`docs/COSPLICE_LOCKED_PANEL_VALIDATION.md`](docs/COSPLICE_LOCKED_PANEL_VALIDATION.md).

Optional / project-dependent:

- Jutils heatmap / PCA outputs
- IsoformSwitchAnalyzeR follow-up
- custom RBP gene set files

---

## System requirements

### Required

- Windows 10/11
- Python 3.11+ recommended
- Local checkout of this repository
- rMATS result tables
- DEG result tables

### Strongly recommended

- Annotation TSV with at least:
  - `gene_id`
  - `gene_name`
  - optionally `transcript_id`, `transcript_name`
- BAM + BAI files for sashimi
- GTF for downstream validation and optional isoform follow-up

### Optional external/runtime dependencies

- `Rscript` for optional R-based downstream modules
- `pysam` for `rmats2sashimiplot`
- Jutils runtime if you want to run Jutils from inside Katana
- IsoformSwitchAnalyzeR runtime if you want isoform follow-up

---

## Repository layout

```text
app/
  main.py
  window.py
  pages/

src/
  services/
  analysis/
  adapters/
  io/
  domain/

software/
  Jutils/
  SUPPA/
  rmats2sashimiplot/

configs/
docs/
tests/
```

---

## Installation

### Option 1: editable install

```powershell
cd E:\Project\katana
python -m venv .venv
.venv\Scripts\activate
pip install -e .
```

### Option 2: quick local run

```powershell
cd E:\Project\katana
python -m app.main
```

### Optional Windows helper

```powershell
.\scripts\setup_windows_dev.ps1
.\launch_katana.bat
```

---

## Launching the app

```powershell
cd E:\Project\katana
python -m app.main
```

Then:

1. Open a project directory
2. Confirm detected inputs on `1. Project`
3. Confirm directions on `2. Pairing`
4. Confirm comparison sets and visual groups
5. Run analysis modules one by one as needed

---

## Expected project inputs

Katana is built to scan a project root and discover usable inputs from that project.

### Minimal expected inputs

- rMATS comparison folders and event tables
- DEG tables

### Additional supported inputs

- DEXSeq gene / exon tables
- DTU gene / isoform tables
- Annotation TSV
- BAM lists or BAM roots
- GTF
- FASTA
- Jutils outputs
- isoform quantification files (`quant.sf`)

### Important behavior

- **Annotation is not hard-linked to one dataset**
- Katana will use the **current project's own annotation TSV**
- If a new project has a different annotation TSV, Katana should use that one instead
- Gene, AS-event and CoSplice-module selections are published through a project-local `AnalysisArtifactRegistry`
- Downstream modules read the selected artifact and its comparison scope; they do not silently rerun an upstream module
- Cached artifacts retain module ID, scope, comparison IDs, parameters, source artifact IDs and output paths
- The Workflow Dashboard can build a cache-only project report with an artifact inventory, module status and provenance JSON

---

## Workflow overview

### 1. Project

Detects:

- comparison inputs
- available software paths
- annotation TSV
- optional BAM / GTF / FASTA / Jutils / isoform resources

### 2. Pairing

Defines:

- canonical splicing-side direction
- canonical DEG direction
- comparison naming
- comparison order

Important:

- AS / rMATS direction is treated as canonical splicing direction
- DEXSeq / DTU inherit this splicing-side direction
- Pairing displays both canonical directions as a read-only audit
- legacy `reverse_splicing` / `reverse_deg` values remain readable for old projects and cached results, but are not editable in the Pairing UI

### 3. Comparison Sets

Defines comparison-pair relationships for cross-comparison analysis.

### 4. Visual Groups

Defines which comparisons are grouped together in downstream comparison views.

### 5. Per-comparison Analysis

#### 5.0 Analysis Thresholds

Threshold controls are kept on a dedicated page instead of occupying every analysis screen. Katana provides exactly three presets plus one manual mode:

- `Default practical`
- `Exploratory`
- `Stringent`
- `Custom`

Changing a preset or custom value invalidates parameter-mismatched caches. It never starts an analysis. The new values are used only after the user explicitly reruns the relevant module.

#### 5.1 Overview

- `5.1.1 Landscape`
- `5.1.2 Bar`
- `5.1.3 Pie`

These focus on significant splicing event composition and counts.

The downstream stable-AS gate is computed identically across overview, integration, candidate, and cross-comparison modules:

- event passes the active rMATS FDR and absolute standardized dPSI cutoff
- each group has the configured number of individually read-supported replicates
- each supported replicate uses `IJC + SJC` against the configured per-replicate minimum
- all cross-group replicate PSI pairs are used for direction-consistency scoring
- missing PSI or read support is reported as unavailable, not converted into false support

Every stable-AS event remains auditable in downstream tables through the group-wise total replicate counts, read-supported replicate counts, mean `IJC + SJC`, direction-consistency score, and support-availability status. Candidate ranking also records a per-gene stability summary for every retained significant AS event and the full stability metrics for its dominant event.

#### 5.2 Transcript + Splicing Integration

- `5.2.1 Transcript + Splicing Integration`
- `5.2.2 Event-type Subanalysis`

These views combine DEG and AS information using standardized directions.

#### 5.3 Candidate Gene Screening

- `5.3.1 Candidate Ranking Table`
- `5.3.2 Candidate Gene Selection`
- `5.3.3 Evidence Heatmap`
- `5.3.4 Candidate Cards + Event Follow-up`
- `5.3.5 Candidate Gene AS Heatmap`

Key behavior:

- each comparison has its own candidate table
- gene ranking is numeric
- candidate selection is independent per comparison tab
- blacklist is independent per comparison tab
- downstream heatmap/cards follow candidate selection and blacklist

#### 5.4 Cross-comparison Candidate Comparison

- `5.4.1 Candidate Matrix`: pairwise gene-level candidate membership for comparison A vs comparison B
- `5.4.2 Shared / Specific / Gained / Lost`: pairwise gene-level pattern explanation
- `5.4.3 Cross-comparison AS interaction modeling`: replicate-aware AS interaction/read-support follow-up
- `5.4.4 Venn / Overlap`: Venn-style intersection counts for gene candidates and AS events from cached 5.4.1/5.4.3 results
- `5.4.5 Selected-comparison Event Venn`: one direction-aware event Venn across user-selected comparisons
- `5.4.6 Selected-comparison Event Similarity`: a separate direction-aware event-similarity heatmap for the same selected comparisons

#### 5.5 Sashimi

- `5.5.2 Run Sashimi` lists available splicing events and lets you select specific events
- `5.5.3 Sashimi Preview` shows generated PNG previews and output files

Sashimi does **not** auto-run on page open.

#### 5.6 Native Splicing Pattern Explorer

- `5.6.1 Sample PSI Heatmaps`: absolute PSI, row-Z-score PSI and hierarchical event/sample pattern views
- `5.6.2 PCA & Replicate QC`: per-comparison sample PSI PCA and within-group versus between-group replicate similarity
- `5.6.3 Similarity & Overlap`: sample PSI Spearman similarity with the shared measured-event denominator, plus significant-event overlap and direction agreement
- `5.6.4 Native Output Browser`: cached TSV, PNG and `parameters.json` preview/export

This is Katana's native replacement for the external Jutils workflow. It preserves the useful Jutils concepts of sample-level PSI heatmaps, hierarchical clustering, PCA and event overlap, while adding explicit measured-event denominators, replicate QC, comparison-aware direction, independent module execution and reproducible cache parameters. Pattern clusters are derived from row-Z-scored sample/replicate PSI trajectories using weighted-linkage hierarchical clustering with city-block distance; they are not clusters of one or two comparison-level dPSI values. White heatmap cells mean PSI was unavailable for that event/sample, not zero.

Each run also writes `native_pattern_summary.tsv`, which records the event/sample/comparison denominators, pair count, median direction and weighted-cosine similarity, low-confidence pair count, active cutoffs and scientific scope. The UI reads this lightweight summary on reopen rather than scanning every large result table.

`5.3.5` is deliberately separate: it restricts the PSI heatmap to genes selected in `5.3.2`, whereas `5.6` is a global selected-experiment pattern/QC module. Opening a `5.6` page reads cache only. Rebuilding requires the explicit Native Run action and does not trigger candidate ranking, sashimi, RBP or isoform analysis. The native overview owns the Run controls; `5.6.1`-`5.6.4` are cache-only result views.

#### 5.7 DEG Explorer

Runs independently for one comparison and uses the read-only canonical DEG direction audited in Pairing.

Outputs per comparison:

- `deg_all.tsv`
- `deg_significant.tsv`
- `deg_up.tsv`
- `deg_down.tsv`
- `summary.tsv`
- `parameters.json`
- `volcano.png`
- `ma_plot.png`

Opening the page reads a matching cache only. `Run DEG Explorer` rebuilds only this module.

#### 5.8 GO / KEGG Enrichment

Runs g:Profiler over-representation analysis for one comparison using one of these explicit gene sources:

- significant DEG
- upregulated DEG
- downregulated DEG
- saved candidate gene selection
- custom gene list

It queries `GO:BP`, `GO:MF`, `GO:CC`, and `KEGG`, then writes the exact input gene list, enrichment TSV, PNG, and `parameters.json`. DEG/candidate inputs must already exist; this module never runs an upstream analysis automatically. Internet access is required for the g:Profiler request. Enrichment is support/prioritization evidence, not a causal claim.

### 6. Samples

Sample-level metadata review and support configuration.

### 7. RBP

RBP support/enrichment is an optional prioritization layer. It checks whether the current per-comparison candidate genes overlap built-in or custom RBP/splicing-factor gene sets.

Outputs:

- `rbp_enrichment_summary.tsv`
- `rbp_candidate_gene_matches.tsv`
- `07_rbp_support/cross_comparison/cross_rbp_enrichment_matrix.tsv`
- `07_rbp_support/cross_comparison/cross_rbp_matched_gene_matrix.tsv`
- `07_rbp_support/cross_comparison/plots/cross_rbp_enrichment_support.png`
- `07_rbp_support/cross_comparison/plots/cross_rbp_matched_gene_matrix.png`
- `07_rbp_support/cross_comparison/plots/pairwise_rbp_rank_support_<comparisonA>_vs_<comparisonB>.png`

The `7.2 Cross-comparison RBP support` page includes a `Plots / Preview` tab. Page refresh reads cached TSV/PNG outputs only; `Run Cross-comparison RBP support` regenerates this module without triggering 5.2 / 5.3 / 5.4 / 5.5.

This is not a causal RBP model; it is a transparent candidate-support and follow-up prioritization module.

### 8. Isoform

Isoform support has two parts:

- sample design / IsoformSwitchAnalyzeR setup
- candidate isoform support from selected candidate genes, `quant.sf`, GTF transcript annotation and optional DTU q values

Outputs:

- `candidate_isoform_gene_summary.tsv`
- `candidate_isoform_transcript_support.tsv`
- `08_isoform_candidate_support/cross_comparison/cross_isoform_gene_support_matrix.tsv`
- `08_isoform_candidate_support/cross_comparison/cross_isoform_transcript_support_matrix.tsv`
- `08_isoform_candidate_support/cross_comparison/plots/cross_isoform_gene_support_heatmap.png`
- `08_isoform_candidate_support/cross_comparison/plots/cross_isoform_dominant_isoform_switches.png`
- `08_isoform_candidate_support/cross_comparison/plots/pairwise_isoform_delta_tpm_<comparisonA>_vs_<comparisonB>.png`

The `8.2 Cross-comparison isoform support` page includes a `Plots / Preview` tab. Page refresh reads cached TSV/PNG outputs only; `Run Cross-comparison Isoform support` regenerates this module without triggering 5.2 / 5.3 / 5.4 / 5.5.

Isoform evidence is optional support. It does not replace DEG+rMATS candidate definition.

### 9. Input Check

Checks whether expected supporting inputs are present.

### 10. Settings

Allows path-level runtime configuration.

---

## Expected outputs

Katana requires an explicit project output directory. New projects prompt for
this directory during loading, and it can be changed later in `Settings`.
Katana rejects output locations inside the Katana source checkout so analysis
artifacts cannot pollute the tool repository. A typical selected output
directory is:

```text
katana_output/
```

For legacy projects that already store `katana_output` inside the source
checkout, use `Settings > Migrate Existing Outputs`. Katana copies the
recognized output tree to the selected external directory, updates the project
configuration only after the copy succeeds, and then removes the old output
tree. Unrecognized input/project folders are refused rather than moved.

Typical sections include:

```text
katana_output/
  comparison_config.yaml
  01_per_comparison_analysis/
  02_candidate_gene_screening/
  03_cross_comparison_candidate_comparison/
  04_program_heatmap/
  04a_candidate_gene_screening/
  05_deg_explorer/
  05b_functional_enrichment/
  06_native_splicing_patterns/
  06_sashimi/
  07_rbp_support/
  08_isoform_candidate_support/
```

Examples of important files:

- `gene_level_integrated_candidates.tsv`
- `tier1_candidates.tsv` ... `tier5_candidates.tsv`
- `cross_comparison_candidate_matrix.tsv`
- `heatmap_gene_list.tsv`
- generated sashimi event files / manifests / plot folders
- Native sample-PSI matrices, PCA/QC, similarity, overlap and pattern-cluster outputs
- RBP support/enrichment tables
- isoform candidate support tables

---

## Direction behavior

### Splicing-side canonical direction

Detected from the project input/configuration, displayed read-only in `Pairing`, then inherited by:

- rMATS
- DEXSeq
- DTU
- event-level dPSI interpretation
- sashimi-side comparison labels

### DEG canonical direction

Detected from the project input/configuration and displayed read-only in `Pairing`.

Downstream tables use `standardized_log2FC`; legacy projects keep their saved direction transforms for reproducibility.

---

## Sashimi behavior

Katana uses **event selection first**, not automatic bulk generation.

Katana now supports two sashimi runners:

- `Katana native`: dependency-light native sashimi-like plots from selected rMATS events, BAM/BAI files and optional GTF annotation.
- `rmats2sashimiplot legacy`: the original external runner, kept as an optional compatibility path.

Current intended flow:

1. Open `5.5.2 Run Sashimi`
2. Filter/search available splicing events
3. Multi-select events
4. Click `Generate selected manifest`
5. Choose `Katana native` or `rmats2sashimiplot legacy`
6. Click `Run selected sashimi`
7. Inspect:
   - `5.5.3 Sashimi Preview`

If sashimi fails, event-level error information should include:

- gene
- event_type
- event_id
- comparison
- error message
- output folder

Common failure cause:

- missing `bamnostic` and `samtools` for `Katana native`. Katana native uses `samtools` when available and falls back to the pure-Python `bamnostic` reader on Windows.
- missing `pysam` in the runtime used by `rmats2sashimiplot`

The standalone native runner can also be called directly:

```powershell
python scripts\katana_native_sashimi.py --manifest selected_manifest.tsv --output-root katana_output --gtf annotation.gtf
```

---

## Jutils behavior

Katana keeps two complementary paths:

- `5.6 Native Splicing Pattern Explorer` provides sample-level PSI heatmaps,
  PCA, replicate QC, pattern clustering, measured-event denominators and
  quantitative comparison similarity.
- The legacy compatibility page can still run bundled Jutils, browse generated
  outputs, preview tables/images/HTML and expose PDF files through open/download.
- Direction-aware event Venn analysis is provided in `5.4.5`; selected-event
  Native Sashimi is provided in `5.5`.

Note:

- Native 5.6 is functionally equivalent or stronger for pattern exploration,
  but it is not a pixel-identical reimplementation of every Jutils CLI option.
- A strict same-input Setd2 benchmark is documented in
  [docs/JUTILS_HEAD_TO_HEAD_VALIDATION_v1.0.1.md](docs/JUTILS_HEAD_TO_HEAD_VALIDATION_v1.0.1.md).
- many legacy Jutils heatmap/PCA outputs are PDFs
- embedded preview may not always be available
- Katana should still expose file paths and open/download actions clearly

---

## Testing

### Basic static check

```powershell
python -m py_compile app\window.py
python -m py_compile app\pages\results_page.py
python -m py_compile app\pages\sashimi_page.py
python -m py_compile app\pages\jutils_page.py
python -m py_compile app\pages\deg_page.py
python -m py_compile app\pages\enrichment_page.py
python -m py_compile src\services\project_service.py
pytest -q
```

### GUI smoke testing

Recommended manual checks:

1. Open a project
2. Confirm `Project / Pairing / Comparison Sets / Visual Groups`
3. Open:
   - `5.1.1 Landscape`
   - `5.2.1 Transcript + Splicing Integration`
   - `5.3.1 Candidate Ranking Table`
   - `5.3.2 Candidate Gene Selection`
   - `5.5.2 Run Sashimi`
   - `5.6.1 Sample PSI Heatmaps`
   - `5.6.2 PCA & Replicate QC`
   - `5.6.3 Similarity & Overlap`
   - `5.6.4 Native Output Browser`
   - `5.7 DEG Explorer`
   - `5.8 GO / KEGG Enrichment`
4. Confirm:
   - project load does not auto-run heavy downstream modules
   - each page only runs when its own button is clicked
   - candidate selection is independent per comparison
   - sashimi only runs on selected events

### Recommended release checklist

- project opens successfully
- no page-open crash
- no automatic chained execution
- candidate ranking loads
- candidate selection persists per comparison
- sashimi event selection works
- Native Splicing Pattern Explorer lists and previews cached TSV/PNG outputs
- optional pages are hidden when required inputs are absent

### Validation summary for the current development version

The current repository was smoke-tested on `2026-07-13` against a real local project root.

Validated behaviors:

- project load completes without automatically triggering heavy downstream analysis
- comparison-aware pages open without immediate crashes
- `5.1` overview pages open
- `5.2` transcript + splicing integration pages open
- `5.3.1` candidate ranking page opens
- `5.3.2` candidate gene selection page opens
- `5.5.2 / 5.5.3` sashimi pages open in distinct modes
- `5.6.1 / 5.6.2 / 5.6.3 / 5.6.4` Native Splicing Pattern Explorer pages open and detect cached PSI/PCA/QC/similarity outputs
- `7 RBP` and `8 Isoform` open as explicit downstream support modules
- `7.2` and `8.2` cross-comparison support modules generate TSV + PNG outputs and reload cached preview files without rerunning analysis
- project confirmation flags survive reload
- saved comparison IDs remain canonical when scanner aliases differ
- stable-AS calculations are shared across analysis layers and covered by regression tests
- candidate caches reject missing or mismatched parameter fingerprints
- `5.7` generates comparison-specific DEG TSV/PNG outputs and reloads them from cache without triggering another module
- `5.8` generates GO/KEGG TSV/PNG outputs from an existing gene list without triggering an upstream module

Genericity checks:

- no active source-level hard-coded `Ythdc1`, `WT`, `KO`, `DP`, or project-specific comparison labels remain in `app/` and `src/`
- annotation lookup is project-specific and loaded from the current project's discovered annotation TSV
- comparison order and direction are loaded from project state rather than hard-coded assumptions

Known runtime caveat:

- successful sashimi plotting still depends on the external `rmats2sashimiplot` runtime and its Python dependencies such as `pysam`

---

## Known limitations

- Some previews (especially PDF outputs) are not embedded and must be opened externally
- Sashimi execution requires the runtime dependencies expected by `rmats2sashimiplot`
- Some optional modules depend on project-specific resource availability
- GO/KEGG enrichment requires internet access to the g:Profiler API
- Large projects can still require manual module-by-module execution time

---

## Release note for current workflow

This version is intended as a **modular 1.0-style desktop workflow**, with:

- independent module execution
- comparison-aware outputs
- read-only canonical direction auditing with standardized downstream effects
- project-specific annotation loading
- event-driven sashimi follow-up

---

## GitHub repository short description

Suggested repository description:

`Comparison-first, direction-aware desktop workflow for RNA-seq + alternative splicing candidate screening, cross-comparison follow-up, and event-driven sashimi review.`

---

## Author

Prepared for release by **Tara**.
