# Approved Katana report specification

Specification: **1.0 — approved 22 October 2026**.

The user approved the bulk and single-cell dummy reports and explicitly asked
to record this design as settled. This document fixes report content and
behaviour. Approval of the design does not mean that the general report
generator or the full distributed package has passed final acceptance.

## Report contents

| Section | Required content |
| --- | --- |
| 1. Design, inputs and quality control | Declared samples and contrasts, data completeness, reference/annotation identities and available sequencing, read-length, coverage and 5′/3′ quality measurements. |
| 2. AS landscape | Per-comparison event-type bars and pies, raw versus read-support-only counts, PCA, PSI heatmaps, pattern clustering, within-group consistency, outliers and complete input/result links. Sample-similarity plots remain in outputs. |
| 3. Differential molecular results | DEG volcano/heatmap with gene labels; DAS event volcano and ΔPSI; DTU transcript usage; DEXSeq exon usage; four-layer gene coverage, intersections and classes. Each comparison has separate Top 10 genes/events/transcripts/exon bins with expandable original fields and full-result links. |
| 4. Cross-comparison overlap and interaction | Exact-event and applicable exon-usage comparisons, Venn/UpSet, directions and explicitly supported interaction tests. No automatic rescue/restoration claim. Each result links its complete table. |
| 5. GO, KEGG and pathway enrichment | GO BP/MF/CC, KEGG and supported ranked/joint pathway results by comparison and evidence set. Link exact input genes and the matching detected/testable background; do not print long gene lists in the report. Retain database identities, contributing genes and original P/FDR. |
| 6. Katana-Splice-Network | All modules and unassigned events; sizes, distributions, event/gene composition, members, structure, PSI patterns, activities, condition associations, hubs, continuous membership metrics, enrichment and contributing genes. Show available bootstrap, structure, sensitivity and differential-network results, and external preservation when computed. |
| 7. RBP and MotifTrace | Keep RBP identity/set enrichment, matched event-neighbourhood motif/CLIP evidence and RNA sequence consequences distinct. Include motif logos and target RNA context, multiple-model browsing, applicable gain/loss/retention, aggregate/window outputs and supported completed parallel results. |

The appendix records methods, actual parameters, exclusions, resources,
versions, source records, output locations and Katana citation information.
Do not invent a publication or DOI.

## Presentation and result access

- Reports are generated in English. The approved report palette is B blue/red:
  positive `#C63F4A`, negative `#2863A5`, other/unassigned `#CDD1D6`, and a
  blue–white–red diverging heatmap. Keep nucleotide identity colours separate.
- Full tables, matrices, sample-similarity plots and complete heatmap pages
  remain in the output directory with search, paging or expansion and download.
  Report previews never replace the complete results.
- Rank Top 10 within one assay, comparison and applicable cell-type scope:
  FDR ascending, absolute native effect descending, then stable object ID.
  Keep multiple events/transcripts/exon bins belonging to the same gene.
  Show fewer than ten when fewer qualify; do not fill with non-significant rows.
- Motif enrichment ordering retains each method's statistical family and its
  FDR/P ordering. Keep RBP-group statistics distinct from individual-model or
  site statistics; RNA matches are ordered by RNA position. No fabricated
  per-model/site P/FDR is permitted.

## DAS selection and scientific boundaries

- Default: **RS-DAS**, requiring native caller significance/effect criteria
  plus IJC + SJC ≥ 10 in every selected biological sample of both contrast arms.
- Ordinary DAS uses the same native significance/effect criteria without the
  read-support foreground gate. Retain original full tables, source P/FDR,
  effects and separate significance/read-support flags.
- One report-wide choice is accessible in report controls, section 3 and DAS
  Top 10 tables. Supporting DAS panels, exact-event overlaps, directions and
  four-layer summaries must use that same choice.
- Supporting downstream results must have matching selection provenance.
  Missing matching enrichment results are **not computed** for that view;
  never relabel RS results as ordinary DAS or reuse their fitted P/FDR for a
  changed foreground/background.
- Explicitly identify modules requiring RS. Completed imported collections
  retain their actual scope. Networks use their declared measurable PSI
  population, including eligible non-DAS events; this switch does not change it.
- New-project mitochondrial-encoded and pseudogene exclusions default on
  before differential inputs. Preserve raw inputs/results and describe the
  actual exclusions in methods/footnotes. Do not claim that an old fitted
  analysis was retrospectively fitted with those exclusions.
- Record runtime thresholds, methods and database versions. Previously agreed
  scientific defaults remain in [thresholds](THRESHOLDS.md); configuration and
  visualisation do not silently change frozen statistical algorithms.

## Single-cell and single-nucleus reports

Use the same seven-section structure with relevant cell/donor/cell-type QC,
embeddings, compositions and donor-aware result scopes. Biological replicates
are donors or the declared biological units, not individual cells. Keep each
cell type's results separate from cohort-level results. Only enable analyses
supported by the assay, inputs and implementation; label missing or unsupported
results explicitly. This does not expand the bulk raw-data upstream pipeline
into a single-cell FASTQ pipeline.

## Content outside the one-click report

- Selected gene/event dossiers, gene cards, Sashimi and focused read validation
  belong in interactive browsing or a separate user-selected analysis.
- No dedicated ORF/UTR/Pfam/NMD or predicted protein-consequence chapter.
- No generated experimental-purpose narrative, automatic rescue interpretation,
  integrated candidate priority or proposed experimental validation programme.
  Descriptive statistics and documented computational results remain available.

## Implementation and change rules

At approval, both dummy reports and their offline packages were verified. The
development package provides existing analyses, Python/CLI/desktop entrypoints
and selection-aware evidence exports. **Connecting this approved seven-section
design to arbitrary eligible completed projects remains implementation work.**

Next: implement that generic renderer, validate it using existing completed
real results, finish current-version distribution/headless installation checks,
and connect the downstream package to Nextflow. Completed upstream analyses
must be reused. Preserve the original application and frozen scientific code.

Routine packaging, adapters, tests and fixes that preserve these decisions can
proceed without renewed approval. Changes to approved report content, major
visual behaviour, scientific defaults or scope require a new explicit decision
and a new specification revision. Keep prior approved examples recoverable.
