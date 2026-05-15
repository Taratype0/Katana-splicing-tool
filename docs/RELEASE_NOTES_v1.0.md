# Katana v1.0 Release Notes

## Overview

Katana v1.0 establishes a comparison-first, direction-aware desktop workflow for RNA-seq and alternative splicing candidate screening.

This release focuses on:

- per-comparison analysis instead of mixed all-project views
- explicit direction handling across DEG and splicing-side analyses
- independent module execution
- candidate-gene-centered follow-up
- event-driven sashimi selection
- project-specific annotation mapping

## Highlights

### Comparison-first workflow

- comparisons keep their own display names, outputs, cache, and logs
- comparison order is preserved through downstream visualization and selection
- analysis modules read upstream state without automatically chaining downstream execution

### Direction-aware analysis

- canonical splicing-side direction is configured in Pairing
- rMATS, DEXSeq, and DTU inherit the canonical splicing-side direction
- DEG final direction can be chosen manually when DESeq2 file naming is unreliable

### Candidate screening

- per-comparison candidate ranking table
- per-comparison candidate gene selection
- blacklist-aware downstream follow-up
- gene-centered candidate cards with event-level follow-up

### Cross-comparison follow-up

- candidate matrix
- shared / specific / gained / lost pattern summaries
- direction reversal and direction-vs-strength follow-up
- significant AS event dPSI pattern analysis across comparisons

### Sashimi workflow

- no automatic bulk sashimi generation on page open
- event list first, then user-driven multi-selection
- manifest generation and sashimi execution are restricted to selected events
- failed jobs are reported with event-level details

### Jutils workflow

- separate heatmap, PCA, and output-browser pages
- output files can be inspected even when embedded preview is unavailable

## Validation summary

Validated in the current release candidate:

- project load succeeds without automatic heavy downstream execution
- overview, transcript + splicing, candidate ranking, candidate selection, sashimi, and Jutils pages open without immediate crashes
- isoform navigation is hidden when isoform inputs are absent
- source-level project-specific hard-coding has been removed from active `app/` and `src/` code paths

## Known limitations

- PDF outputs are not always previewed inline
- sashimi execution still depends on the local `rmats2sashimiplot` runtime and Python dependencies such as `pysam`
- optional modules remain dependent on project-specific resources and tool availability

## Author

Prepared for release by **Tara**.
