# Approved Katana report specification

Specification: **1.1 — final report format approved 25 October 2026**.

The user approved the latest completed bulk report as the fixed format for
future one-click reports. It supersedes the earlier presentation examples;
the original bulk and single-cell design approval of 22 October 2026 remains
part of the history. Preserve this report's structure, visual style and
interaction when updating results or rendering another eligible project.
The machine-readable decision is in [FORMAT_APPROVAL.json](reports/FORMAT_APPROVAL.json).
Report-format approval and software/pipeline release acceptance are separate.

## Report contents

| Section | Required content |
| --- | --- |
| 1. Design, inputs and quality control | Declared samples and contrasts, data completeness, reference/annotation identities and available sequencing, read-length, coverage and 5′/3′ quality measurements. |
| 2. AS landscape | Per-comparison event-type bars and pies, raw versus read-support-only counts, PCA, PSI heatmaps, pattern clustering, within-group consistency, outliers and complete input/result links. Sample-similarity plots remain in outputs. |
| 3. Differential molecular results | DEG volcano/heatmap with gene labels; DAS event volcano and ΔPSI; DTU transcript usage; DEXSeq exon usage; four-layer gene coverage, intersections and classes. Each comparison has separate Top 10 genes/events/transcripts/exon bins with expandable original fields and full-result links. |
| 4. Cross-comparison overlap and interaction | Exact-event and applicable exon-usage comparisons, Venn/UpSet, directions and explicitly supported interaction tests. No automatic rescue/restoration claim. Each result links its complete table. |
| 5. GO, KEGG and pathway enrichment | GO BP/MF/CC, KEGG and supported ranked/joint pathway results by comparison and evidence set. Link exact input genes and the matching detected/testable background; do not print long gene lists in the report. Retain database identities, contributing genes and original P/FDR. |
| 6. Katana-Splice-Network | Direct interactive module network and complete searchable/downloadable gene lists; all modules and unassigned events, sizes, event/gene composition, members, PSI patterns, activities, condition associations, hubs, continuous membership metrics, significant enrichment and contributing genes. Graph display limits are explicit; complete membership is retained. Dendrograms and supporting diagnostics stay in outputs. |
| 7. RBP and MotifTrace | Keep RBP identity/set enrichment, matched event-neighbourhood motif/CLIP evidence and RNA sequence consequences distinct. Show logos and target RNA context for up to ten distinct motifs with recorded native test FDR <= 0.05; otherwise show a concise no-significance statement. Complete model, window and nonsignificant results remain in outputs; supported completed parallel methods retain their own correction families. |

The appendix records methods, actual parameters, exclusions, resources,
versions, source records, output locations and Katana citation information.
Do not invent a publication or DOI.

## Presentation and result access

- Reports are generated in English. The approved report palette is B blue/red:
  positive `#C63F4A`, negative `#2863A5`, other/unassigned `#CDD1D6`, and a
  blue–white–red diverging heatmap. Keep nucleotide identity colours separate.
- Open on one declared comparison, with a persistent comparison selector.
  Switching it updates all comparison-scoped chapters together. Shared study
  QC remains visible; grouped cross-comparison panels retain their declared
  axes and appear only for the selected comparison's group. Do not return to
  displaying every comparison together by default.
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

The generic Python/CLI report renderer and read-only completed-report refresh
are implemented. The accepted real bulk reference passed selected-comparison
browser checks, complete module-gene browsing, original-statistic and portable
file checks, offline complete-table reading, and archive inventory/CRC checks.
The accepted report and its complete archive remain immutable references.

Additional project adapters, assay-specific acceptance, current-version CI and
distribution checks, and complete pipeline-entry acceptance retain their own
evidence and status. This approval does not declare all of them complete or
turn single-cell design approval into real single-cell acceptance. Reuse
completed analyses and preserve the original application and scientific code.

Routine packaging, adapters, tests and fixes that preserve these decisions can
proceed without renewed approval. Changes to approved report content, major
visual behaviour, scientific defaults or scope require a new explicit decision
and a new specification revision. Keep prior approved examples recoverable.
Do not redesign this format during routine packaging or integration work.
