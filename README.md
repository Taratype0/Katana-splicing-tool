# Katana RNA-seq + Splicing Candidate Screening Tool

Katana is a desktop GUI workflow for **comparison-first**, **direction-aware** RNA-seq and splicing follow-up analysis.

It is designed to work across projects, not only one specific dataset. Each project is scanned locally, each comparison keeps its own direction/configuration, and each analysis module runs independently instead of chaining downstream modules automatically.

## Core design principles

- **Comparison-first**: each comparison keeps its own pages, outputs, cache, and logs
- **Direction-aware**: splicing-side direction is inherited from Pairing; DEG final direction can be set manually
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
- Transcript + splicing integration
- Per-comparison candidate gene ranking and selection
- Candidate cards + event follow-up
- Cross-comparison candidate matrix and AS dPSI pattern comparison
- Event-driven sashimi manifest generation and sashimi execution
- Jutils output browsing
- RBP support/enrichment over candidate genes
- Isoform candidate support from quant.sf / GTF / DTU

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
- manual DEG final direction
- comparison naming
- comparison order

Important:

- AS / rMATS direction is treated as canonical splicing direction
- DEXSeq / DTU inherit this splicing-side direction
- DEG direction can be chosen manually when DESeq2 file naming is unreliable

### 3. Comparison Sets

Defines comparison-pair relationships for cross-comparison analysis.

### 4. Visual Groups

Defines which comparisons are grouped together in downstream comparison views.

### 5. Per-comparison Analysis

#### 5.1 Overview

- `5.1.1 Landscape`
- `5.1.2 Bar`
- `5.1.3 Pie`

These focus on significant splicing event composition and counts.

#### 5.2 Transcript + Splicing Integration

- `5.2.1 Transcript + Splicing Integration`
- `5.2.2 Event-type Subanalysis`

These views combine DEG and AS information using standardized directions.

#### 5.3 Candidate Gene Screening

- `5.3.1 Candidate Ranking Table`
- `5.3.2 Candidate Gene Selection`
- `5.3.3 Evidence Heatmap`
- `5.3.4 Candidate Cards + Event Follow-up`

Key behavior:

- each comparison has its own candidate table
- gene ranking is numeric
- candidate selection is independent per comparison tab
- blacklist is independent per comparison tab
- downstream heatmap/cards follow candidate selection and blacklist

#### 5.4 Cross-comparison Candidate Comparison

- `5.4.1 Candidate Matrix`: pairwise gene-level candidate membership for comparison A vs comparison B
- `5.4.2 Shared / Specific / Gained / Lost`: pairwise gene-level pattern explanation
- `5.4.3 Cross-comparison AS dPSI Pattern`: event-level shared significant AS dPSI pattern prioritization
- `5.4.4 Cross-comparison AS interaction modeling`: replicate-aware AS interaction/effect follow-up

#### 5.5 Sashimi

- `5.5.2 Run Sashimi` lists available splicing events and lets you select specific events
- `5.5.3 Sashimi Preview` shows generated PNG previews and output files

Sashimi does **not** auto-run on page open.

#### 5.6 Jutils

- `5.6.1 Jutils Heatmap`
- `5.6.2 Jutils PCA`
- `5.6.3 Jutils Output Browser`

If preview is not embedded (for example PDF output), Katana shows the output path and provides open/download actions.

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

Katana writes output under a project-specific output root, typically:

```text
katana_output/
```

Typical sections include:

```text
katana_output/
  comparison_config.yaml
  01_per_comparison_analysis/
  02_candidate_gene_screening/
  03_cross_comparison_candidate_comparison/
  04_program_heatmap/
  05_jutils/
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
- Jutils output files
- RBP support/enrichment tables
- isoform candidate support tables

---

## Direction behavior

### Splicing-side canonical direction

Defined in `Pairing`, then inherited by:

- rMATS
- DEXSeq
- DTU
- event-level dPSI interpretation
- sashimi-side comparison labels

### DEG final direction

Set manually in `Pairing` when needed.

Katana does **not** blindly trust DESeq2 file naming if the user overrides the final direction.

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

Katana can:

- run Jutils
- browse generated outputs
- preview tables/images/HTML
- guide the user to PDF outputs via open/download

Note:

- many Jutils heatmap/PCA outputs are PDFs
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
python -m py_compile src\services\project_service.py
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
   - `5.6.1 Jutils Heatmap`
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
- Jutils output browser can list generated files
- optional pages are hidden when required inputs are absent

### Validation summary for the current 1.0 candidate

The current repository was smoke-tested on `2026-05-07` against a real local project root.

Validated behaviors:

- project load completes without automatically triggering heavy downstream analysis
- comparison-aware pages open without immediate crashes
- `5.1` overview pages open
- `5.2` transcript + splicing integration pages open
- `5.3.1` candidate ranking page opens
- `5.3.2` candidate gene selection page opens
- `5.5.2 / 5.5.3` sashimi pages open in distinct modes
- `5.6.1 / 5.6.2 / 5.6.3` Jutils pages open and detect available output files
- `7 RBP` and `8 Isoform` open as explicit downstream support modules
- `7.2` and `8.2` cross-comparison support modules generate TSV + PNG outputs and reload cached preview files without rerunning analysis

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
- Large projects can still require manual module-by-module execution time

---

## Release note for current workflow

This version is intended as a **modular 1.0-style desktop workflow**, with:

- independent module execution
- comparison-aware outputs
- manual direction control where needed
- project-specific annotation loading
- event-driven sashimi follow-up

---

## GitHub repository short description

Suggested repository description:

`Comparison-first, direction-aware desktop workflow for RNA-seq + alternative splicing candidate screening, cross-comparison follow-up, and event-driven sashimi review.`

---

## Author

Prepared for release by **Tara**.
