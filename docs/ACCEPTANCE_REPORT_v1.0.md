# Katana v1.0 Acceptance Report

Date: 2026-05-07

Prepared for release by **Tara**.

## Scope

This acceptance pass focused on:

- project loading behavior
- independent module/page opening without chained execution
- key GUI page smoke tests
- cross-comparison (`5.4`) page validation
- cache-assisted downstream end-to-end regression from real upstream outputs
- genericity checks against project-specific hard-coding in active source code
- release documentation completeness

## Validation method

### 1. Static compilation

The following modules were compiled successfully with `python -m py_compile`:

- `app/window.py`
- `app/pages/results_page.py`
- `app/pages/candidate_gene_page.py`
- `app/pages/candidate_selection_page.py`
- `app/pages/sashimi_page.py`
- `app/pages/jutils_page.py`
- `src/services/project_service.py`
- `src/adapters/cards_adapter.py`

### 2. Real-project smoke test

The application was smoke-tested against a real local project root:

`E:\Project\RNASEQ\Result\Ythdc1_backup`

Validated behaviors:

- project load completed successfully
- project load did not automatically trigger heavy downstream analysis
- comparison selection state was normalized on load
- the following pages could be instantiated and refreshed without immediate crash:
  - `5.1.1 Landscape`
  - `5.1.2 Bar`
  - `5.1.3 Pie`
  - `5.2.1 Transcript + Splicing Integration`
  - `5.2.2.2 RI`
  - `5.3.1 Candidate Ranking`
  - `5.3.2 Candidate Gene Selection`
  - `5.5.2 Run Sashimi`
  - `5.5.3 Sashimi Preview`
  - `5.5.4 Failed Sashimi Jobs`
  - `5.6.1 Jutils Heatmap`
  - `5.6.2 Jutils PCA`
  - `5.6.3 Jutils Output Browser`
- main window initialization succeeded
- `7. Isoform` navigation was hidden when isoform inputs were absent

### 3. Cross-comparison (`5.4`) validation

With cached upstream results available, the following cross-comparison pages were opened and refreshed successfully without immediate crash:

- `5.4.1 Candidate Matrix`
- `5.4.2 Shared / Specific / Gained / Lost`
- `5.4.3 Direction Reversal`
- `5.4.6 Cross-comparison significant AS event dPSI pattern analysis`

Validation notes:

- `5.4.1 / 5.4.2 / 5.4.3` no longer behave as empty pages
- `5.4.6` pattern recalculation no longer fails with `NameError: np is not defined`
- for the currently configured comparison pairs in the validation project, `5.4.6` returned `0` matching rows under the default cutoffs, which is acceptable behavior when the UI reports `No events found`

### 4. Cache-assisted downstream end-to-end regression

Because a fresh uncached `5.2.1 Transcript + Splicing Integration` recomputation did not complete within long timeout windows in the current environment, downstream regression was validated from real cached upstream outputs already present in the project:

- `katana_output/04_tx_splicing_integration/data/tx_splicing_gene_table.tsv`
- `katana_output/04_tx_splicing_integration/data/tx_splicing_summary.tsv`

Using those real upstream outputs, the following downstream modules were exercised successfully:

- `5.3 Mechanism support preview`
- `5.3.1 Candidate ranking rebuild`
- `5.4 Cross-comparison candidate matrix`
- `5.4.6 Cross-comparison significant AS event dPSI pattern analysis`
- `5.5 Sashimi event listing and selected-event manifest generation`
- `5.6 Jutils output browser`

Observed regression results on the validation project:

- mechanism support rows: `9018`
- candidate ranking rows: `97429`
- cross-comparison candidate matrix rows: `3385`
- sashimi available events (selected comparison): `78894`
- selected-event sashimi manifest rows (2 selected events): `2`
- Jutils output browser rows: `52`

Interpretation:

- downstream modules can read real cached upstream outputs and continue independently
- downstream page opening does not require project-load-time chained execution
- sashimi and Jutils paths are live on the validation project
- this is a valid downstream regression pass, but it is **cache-assisted**, not a full uncached rerun from raw inputs

### 5. Genericity check

Active source under `app/` and `src/` was scanned for project-specific hard-coding.

Result:

- no active `Ythdc1`, `WT`, `KO`, `DP`, or similar project-specific labels remain in the main analysis/UI code paths
- the remaining hard-coded group palette in `src/adapters/cards_adapter.py` was replaced with dynamic group coloring
- annotation mapping remains project-specific and is loaded from the current project's discovered annotation TSV

## Acceptance result

### Accepted for 1.0 candidate release with known limitations

The current repository is acceptable as a **v1.0 release candidate** for GitHub publication, with the limitations listed below clearly documented.

## Known limitations

### 1. Sashimi runtime dependency

Sashimi event selection, manifest generation, and failed-job reporting are integrated, but successful plotting still depends on the local `rmats2sashimiplot` runtime and Python dependencies such as `pysam`.

### 2. Embedded preview limitations

Some generated outputs, especially PDFs, are not always embedded inline and may need to be opened externally.

### 3. Genericity validation depth

This release was validated through:

- source-level hard-coding scan
- project load and GUI smoke tests
- single-project real-data testing
- cross-comparison page validation
- cache-assisted downstream regression from real upstream outputs

It was **not** validated end-to-end across multiple unrelated external projects during this acceptance pass.

### 4. Fresh uncached `5.2.1` runtime

A fully fresh uncached rerun starting from `5.2.1 Transcript + Splicing Integration` remains too slow in the current environment for release-signoff level timing confidence.

As a result:

- downstream acceptance is confirmed from real cached upstream outputs
- full fresh end-to-end acceptance from raw inputs remains incomplete

## Release recommendation

Recommended release label:

- `Katana v1.0`

Recommended positioning:

- modular, comparison-first desktop workflow
- direction-aware RNA-seq + AS candidate screening
- per-comparison candidate follow-up
- event-driven sashimi review

## Related release documents

- `README.md`
- `docs/GITHUB_DESCRIPTION.md`
- `docs/RELEASE_NOTES_v1.0.md`
- `docs/GITHUB_RELEASE_BODY_v1.0.md`
