# Browser exploration and selected-target reports

Product decisions recorded on 23 October 2026. The selected-target workflow
below is an implementation specification, not a claim that its portable batch
command or complete browser workspace has already shipped.

## Browser first

The browser is the primary interface. A local or server-side Python service
performs computation; a separate desktop distribution is not a current release
target. Existing explicit desktop entry points remain compatibility routes.

Use a restrained consulting-report interface with dark blue, white and neutral
surfaces. Scientific data retain the approved report colours: positive
`#C63F4A`, negative `#2863A5`, unassigned `#CDD1D6`, and blue-white-red heatmaps.
This is a style direction, not an approved new interface mock-up or an official
third-party brand template.

Organise study browsing around the seven sections in
[the approved report specification](REPORT_SPECIFICATION.md), then provide a
separate **Gene / Event workspace**. Do not add selected-target dossiers to the
seven-section report body or silently change its approved format.

## Thresholds and recorded results

- Let users change DEG FDR and absolute log2 fold change, DAS FDR and absolute
  delta PSI, and the ordinary-DAS / RS-DAS view. Start with the agreed defaults.
- Refresh the corresponding displayed plots, tables and Top 10 from original
  effects and P/FDR. A display change does not refit upstream models.
- Save both the actual source-result parameters and the user's browsing
  settings. Do not describe an unexecuted parameter choice as a completed fit.
- A different enrichment foreground or background needs a correspondingly
  computed result. Mark the old enrichment as belonging to its original set;
  do not relabel its P/FDR as a result for the new selection.
- Subsetting an existing MotifTrace result for a target dossier preserves the
  original whole-family statistics. A newly requested event-set inference,
  different motif library, window or matching design is a separate calculation.
- Network discovery retains its own measurable PSI population. Moving a DAS
  threshold does not automatically refit the network.

## Select targets once

Users can add genes or precise events while browsing differential results,
intersections, enrichment contributors, network members or MotifTrace results.
Keep a persistent target basket across these pages.

Gene and event selection have different meanings:

- A gene target uses a resolved stable gene identifier, with its symbol as a
  display label. List its available events for the selected comparisons; show
  which events the user requested. Do not silently keep only three events or
  substitute the event with the largest delta PSI.
- An event target retains its exact event type, geometry, strand, reference,
  source identity and comparison-specific result row. A short rMATS numeric ID
  or a gene symbol alone is insufficient.
- Preserve the gene-to-event relationships, biological sample/donor identities,
  cell type where applicable, contrast arms/direction and JC/JCEC measurement.
- Selections may include non-significant results. Keep their original statistics
  and mark them accurately instead of removing them from the requested output.

## Portable target request

Export a small versioned request plus the necessary identity/evidence tables.
The request records:

| Field group | Required meaning |
| --- | --- |
| Request and study | Request ID, schema version, source project identity and manifest fingerprint, bulk/scRNA/snRNA mode. |
| Targets | Stable gene or exact event identity, source row identity and content fingerprint, gene/event relationships and expansion rule. |
| Design | Requested comparisons and direction, sample/biological-unit identities, donor and cell type where relevant. |
| References | Species, assembly, annotation and sequence identities, coordinate convention and strand. |
| Selection | Actual source parameters, browsing FDR/effect thresholds, ordinary/RS choice, count mode and exclusion policy. |
| Requested output | Gene cards, event cards, saved statistical/network/motif evidence, optional Sashimi and complete tables. |
| Resources | Logical asset identities and verification records, with paths bound separately on the execution machine. |
| Execution | Which completed results must be reused and which missing target computations the user has requested. |

A saved local `selection_id` is not itself a portable request. The destination
must rebind verified source contents and biological identities; it must not
require the original machine's artifact database UUIDs or match files by
basename alone. Keep local/server path bindings separate from the main request.
Do not include credentials in the request.

BAMs, their indexes and large references stay on the machine performing the
requested read analysis. Download target reports and interactive data, not BAMs
by default. A machine without BAM access can still explore saved evidence; it
must not fabricate a Sashimi result.

## Execution and output

The browser service and a future target CLI consume the same request and
planner. An offline report can export a request; actual Python/BAM computation
is performed by the local or server service/CLI.

Reuse completed analysis tables and existing verified target products. Compute
only the requested missing target views or read evidence. Independently track
target, comparison, cell type and module tasks. A failed target retains its
failure record while other requested targets remain available. Slurm execution
needs a real scheduler adapter; a local thread pool is not a Slurm executor.

Deliver a batch index, a separate report per selected gene/event, complete
machine-readable tables, PNG/SVG/PDF figures as supported, provenance and a ZIP.
Gene/event views combine the available DEG, DAS, DTU, DEXSeq, sample-level
measurements, network membership/activity and MotifTrace evidence. Sashimi keeps
the exact requested event and biological replicates. Distinguish complete,
failed, not requested, unavailable and non-estimable results.

## Existing components and remaining work

Existing components include table selection, selection validation/routing,
selected-gene PSI/expression views, exact-event native Sashimi, network and
MotifTrace browsers, and parameter/source-aware caching. The original richer
gene-card logic can be migrated without editing the original application.

Remaining work includes complete browsing-state capture, the cross-page target
basket, portable request export/import and resource rebinding, an integrated
per-target evidence report, local/Slurm task orchestration, and return-result
import. These require end-to-end acceptance before they are advertised as a
completed portable target workflow.

Routine implementation continues on the development branch. Formal releases
and additional distribution channels follow completed development and release
acceptance; an interface mock-up or passing unit tests alone are insufficient.
