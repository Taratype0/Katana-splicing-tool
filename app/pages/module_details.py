from __future__ import annotations


DETAILS: dict[str, str] = {
    "5.1.landscape": """
Purpose:
Landscape is the first global AS overview and the upstream stable-AS gate. It summarizes rMATS events across the selected comparisons before candidate ranking.

Inputs:
- Selected comparisons in the Pairing order.
- rMATS event summaries with standardized dPSI and FDR/q values.
- Current practical AS threshold used by the landscape cache.
- Stable-AS thresholds: minimum mean IJC+SJC reads per group, minimum replicate count per group, and minimum pairwise replicate PSI direction consistency.

Outputs:
- Significant AS event counts and stable AS event counts by comparison and event type.
- Preview table with comparison_id, comparison name, event_type, n_sig, n_stable and stable filter parameters.

Interpretation:
- This is event-count composition, not a candidate-gene score.
- A larger bar means more stable AS events of that type in that comparison.
- n_sig means rMATS FDR/dPSI significant. n_stable means n_sig plus read-support and replicate-direction consistency filters.
- Downstream 5.2 and 5.3 reruns inherit the stable AS definition for rMATS_significant evidence.
- Event colors are fixed across Landscape, Bar and Pie.

Run and cache behavior:
- Opening the page reads cached results only.
- Click Run Landscape Overview to apply stable-AS thresholds and rebuild only this 5.1.1 result.
""",
    "5.1.bar": """
Purpose:
Bar compares stable AS-event counts by event type across selected comparisons.

Inputs:
- Same cached rMATS landscape table as 5.1.1.
- Selected comparison order inherited from Pairing.
- Stable-AS thresholds inherited from 5.1 Landscape controls.

Outputs:
- Grouped event-type bar plot.
- Preview table with event counts.

Interpretation:
- This view is best for comparing SE/RI/A3SS/A5SS/MXE abundance between experiments.
- It is not gene-level ranking and does not use candidate selection.
- The plotted count is n_stable; the data table keeps n_sig as the pre-stability significant count.

Run and cache behavior:
- Opening the page reads cached results only.
- Click Run Bar to rebuild only this 5.1.2 result.
""",
    "5.1.pie": """
Purpose:
Pie shows the stable AS-event composition inside one selected comparison.

Inputs:
- One comparison selected from the comparison tabs.
- rMATS significant and stable event counts.

Outputs:
- Event-type pie chart for the selected comparison using n_stable.
- Preview table with event_type and stable event count.

Interpretation:
- Use this to inspect within-comparison AS composition.
- Use Landscape/Bar when the question is cross-comparison event-count comparison.
- If the pie is empty while n_sig exists, the stable read/replicate filters are stricter than the current event support.

Run and cache behavior:
- Opening the page reads cached results only.
- Click Run Pie to rebuild only this 5.1.3 result.
""",
    "5.2.integration": """
Purpose:
5.2.1 integrates transcript expression and splicing for one comparison.

Inputs:
- DEG table using standardized_log2FC.
- rMATS representative event table using standardized_dPSI.
- rMATS events must pass the stable-AS definition from 5.1 before they are treated as splicing_sig.
- Pairing/comparison direction; downstream pages should not guess direction again.

Outputs:
- Interactive scatter: x = standardized_log2FC, y = standardized_dPSI.
- Classes: DEG only, DAS only, DEG + DAS same direction, DEG + DAS opposite direction.
- Preview table with gene, FDR, dPSI/log2FC, flip flags and direction labels.

Interpretation:
- Same direction means transcript and representative AS event move in concordant standardized signs.
- Opposite direction means DEG and DAS signs disagree.
- This page is the main DEG + splicing relationship view; DEXSeq/DTU are support layers, not replacements for rMATS.
- rMATS/DAS evidence here means FDR/dPSI significant plus stable read support and replicate-direction consistency.

Run and cache behavior:
- Opening the page reads cached 5.2 results only.
- Click Run Transcript + Splicing to rebuild only 5.2.1/5.2 cached integration.
""",
    "5.2.event": """
Purpose:
Event-type subanalysis repeats the 5.2 transcript + splicing integration for one AS event type.

Inputs:
- Cached 5.2 integration table.
- Selected event type: SE, RI, A3SS, A5SS or MXE.

Outputs:
- Event-type-specific DEG + DAS scatter.
- Preview table with representative event_id, FDR, dPSI, log2FC and direction flags.

Interpretation:
- This is still DEG + DAS integration, only filtered to one event type.
- If a gene has multiple significant events of the same type, the representative event is selected by largest abs(standardized dPSI), with lower FDR as tie-breaker.

Run and cache behavior:
- Opening the page reads cached 5.2 results only.
- It does not run 5.3/5.5 or downstream modules.
""",
    "5.3.ranking": """
Purpose:
5.3.1 creates one unique gene-level candidate ranking per comparison.

Inputs:
- DEG and rMATS are the primary candidate definition.
- DEXSeq and DTU are support evidence when present and mappable.
- Gene annotation maps gene_id to gene_symbol; no project-specific hard-coding is used.
- rMATS_significant means rMATS FDR/dPSI significant and stable by the 5.1 read-support/replicate-consistency gate.

Outputs:
- gene_level_integrated_candidates.tsv and tier tables.
- Columns include rank, candidate_tier, evidence_class, evidence_count, significant AS event summary and support flags.

Interpretation:
- rank is numeric and within-comparison only.
- A gene should appear once per comparison in the ranking table.
- rMATS remains primary splicing evidence; DEXSeq/DTU can strengthen support but do not replace rMATS.
- If stable-AS thresholds are too strict, n_rMATS_significant_events may decrease even when raw rMATS FDR/dPSI events exist.

Run and cache behavior:
- Opening 5.3.1 reads cache only.
- Click Generate / Load Candidate Results to build or load the ranking.
""",
    "5.3.selection": """
Purpose:
5.3.2 defines the candidate genes that downstream visualizations will use.

Inputs:
- Cached 5.3.1 ranking for the selected comparison.
- Optional manual gene additions and per-comparison blacklist.

Outputs:
- Per-comparison selected candidate gene list.
- Per-comparison blacklist.

Interpretation:
- Selection is independent for each comparison tab.
- Default selection is top N genes by numeric rank.
- Blacklisted genes are removed from downstream heatmap, candidate cards and sashimi/event follow-up.

Run and cache behavior:
- Opening the page reads current selection state only.
- Loading ranking/Reset Top N affects only the current comparison.
""",
    "5.3.heatmap": """
Purpose:
5.3.3 visualizes multi-layer evidence for the selected candidate genes.

Inputs:
- Candidate gene selection from 5.3.2.
- Blacklist-filtered candidate table.

Outputs:
- Evidence heatmap and heatmap source gene list.
- Preview table for the heatmap input.

Interpretation:
- The heatmap does not randomly sample genes from a preview table.
- Gene source should be visible: candidate selection, top N, shortlist or custom list.

Run and cache behavior:
- Opening the page reads cached/available candidate selection only.
- It does not initialize selection or rebuild ranking automatically.
""",
    "5.3.cards": """
Purpose:
5.3.4 is gene-centered candidate/event follow-up.

Inputs:
- Candidate gene selection from 5.3.2.
- Significant rMATS AS events for the selected gene.
- DEG/transcriptome expression support when available.

Outputs:
- One gene-centered card with comparison pages.
- Significant AS event list per gene/comparison.
- Sashimi entry points for specific events.

Interpretation:
- Multiple significant events are shown as an event list; Katana should not silently hide them.
- Dominant event is selected by max abs(dPSI) among significant events unless otherwise stated.
- Sashimi can be launched from an event, but card preview itself does not run sashimi.

Run and cache behavior:
- Opening the page reads candidate selection and event cache only.
- It does not trigger sashimi or downstream modules.
""",
    "5.4.matrix": """
Purpose:
5.4.1 is a pairwise gene-level candidate membership matrix.

Inputs:
- Cached 5.3.1 candidate rankings for comparison A and comparison B.
- Selected comparison pair.
- The word candidate here means a gene-level candidate from 5.3.1, where DEG, rMATS, DEXSeq and DTU are summarized into one gene row.

Outputs:
- One row per unique gene candidate in A and/or B.
- is_candidate_A/B, tier_A/B, rank_A/B, evidence_class_A/B and direction_class_A/B.
- pair_gene_level_direction_changed plus separate DEG and dominant-AS direction-change flags.
- DEG/rMATS/DEXSeq/DTU columns are support evidence inherited from the 5.3.1 gene-level candidate table; this page does not recalculate those evidence layers.

Interpretation:
- pair_membership_count = 1 means A-only or B-only; 2 means shared in both sides of this pair.
- A and B mean the two comparisons selected by the comparison-set tab, not treatment/control groups inside one comparison.
- rank_A and rank_B are inherited from 5.3.1 and are not a new cross-comparison rank.
- best_pair_rank is just the lowest numeric rank observed across A and B.
- pair_gene_level_direction_changed is a gene-level descriptor. It is true when standardized_log2FC and/or the dominant rMATS dPSI changes sign between A and B.
- This is not strict same-event AS reversal. Use 5.4.3/5.4.4 for event-level AS comparisons.
- This is gene-level membership, not an event-level AS interaction test.
- Use 5.4.1 when you want the raw table: which genes are candidates in A, B, or both.

Run and cache behavior:
- Opening the page reads cached cross-comparison candidate matrix only.
- Use 5.4.3 for AS-event differences with replicate read support.
""",
    "5.4.pattern": """
Purpose:
5.4.2 explains pairwise gene-level candidate patterns.

Inputs:
- Same A-vs-B gene membership matrix as 5.4.1.
- Same gene-level candidate definition as 5.4.1: candidates are inherited from 5.3.1 and may carry DEG, rMATS, DEXSeq and DTU evidence classes.

Outputs:
- Shared, A-specific, B-specific, gene-level direction-changed, tier-changed and evidence-class-changed gene groups.
- A filtered table of the classified gene-level patterns.

Interpretation:
- shared_pair: candidate in both A and B.
- A-specific/lost in B: candidate in A only, relative to selected A -> B order.
- B-specific/gained in B: candidate in B only, relative to selected A -> B order.
- shared_DEG_direction_reversed: candidate in both A and B, and standardized_log2FC changes sign between A and B.
- shared_dominant_AS_direction_reversed: candidate in both A and B, and dominant rMATS standardized dPSI changes sign between A and B.
- shared_mixed_direction_changed: candidate in both A and B, and both DEG and dominant AS signs change.
- These are gene-level direction-change classes. They are not strict same-AS-event reversal calls.
- tier_changed/evidence_class_changed apply only to genes present in both sides.
- This page explains classes; 5.4.1 shows the raw matrix.
- If 5.4.1 and 5.4.2 look similar, that is because they use the same source genes. The difference is that 5.4.2 turns raw A/B membership into named biological follow-up classes.
- gained/lost are pair-order labels only. They are not automatically a time-course claim unless the comparison-set order is biologically ordered.
- Use 5.4.3/5.4.4 for strict event-level AS reversal and read-support interpretation.

Run and cache behavior:
- Opening reads cache only.
- It does not rebuild 5.3 rankings.
""",
    "5.4.dpsi": """
Purpose:
Legacy summary-only dPSI pattern analysis compares the same significant AS event between two comparisons using standardized dPSI only.

Inputs:
- rMATS significant events only.
- Event matching key: gene + event_type + event_id, with coordinate fallback.

Outputs:
- Same-direction large dPSI difference events.
- Opposite-direction events.
- All shared significant event pairwise dPSI delta table.
- Optional A-only/B-only summary table.

Interpretation:
- This is pattern prioritization, not a statistical interaction test.
- dPSI cutoff defaults to 0.10 = 10 percentage points PSI change, a practical AS filter, not a universal law.
- It only uses rMATS summary values: standardized dPSI and FDR.
- It does not use replicate-level IJC/SJC read counts.
- Result table choices:
  - same-direction large delta: both comparisons significant, same dPSI sign, but effect magnitude differs.
  - opposite-direction: both comparisons significant, dPSI sign reverses.
  - all shared significant: all event-level pair rows that pass significance in both comparisons.
  - A-only/B-only summary: significant in one comparison but not the other; optional follow-up, not the main shared-event analysis.
- This page is hidden from the main workflow because it is weaker than the read-count-aware 5.4.3 module.
- Use 5.4.3 when replicate-level IJC/SJC modeling is needed.

Run and cache behavior:
- Apply/Recalculate runs only this legacy dPSI summary if opened from a legacy/internal route.
- Opening reads cached tables only.
""",
    "5.4.interaction": """
Purpose:
5.4.3 tests/supports cross-comparison AS effect differences with replicate-level rMATS read counts.

Inputs:
- Shared or exclusive rMATS events.
- Replicate IJC/SJC counts and standardized dPSI.
- Selected comparison pair or all-pair batch mode.

Outputs:
- Unified event result table: one row per AS event in the selected comparison pair.
- The unified table contains shared modeled events, A-specific events, B-specific events, benchmark columns and model diagnostic columns together.
- All-pair summary table: one row per comparison pair for batch-level counts.
- Secondary split TSV files are still written for audit/export, but the main UI uses the unified table.

Interpretation:
- Quasi-binomial diagnostics estimate overdispersion using Pearson residual scale.
- Beta-binomial columns come from the formal likelihood-ratio interaction fit when available.
- A-specific means an event is significant in comparison A, then checked for counterpart read support in comparison B.
- B-specific means an event is significant in comparison B, then checked for counterpart read support in comparison A.
- This is the main AS-event cross-comparison module because it uses replicate read support, not only summary dPSI.
- Unified table columns:
  - unified_display_class separates Shared significant modeled, Opposite direction, Same direction + large dPSI delta, A-specific significant and B-specific significant.
  - support_class records counterpart read-support status for A/B-specific events.
  - benchmark_class records whether simple delta dPSI ranking and model-supported ranking agree.
  - model_diagnostic_class and beta_binomial_model_status record model reliability.
- The summary plot at the top is intentionally compact; the unified event table is the main output. Use Data/Export to inspect and filter rows.

Run and cache behavior:
- Apply/Recalculate runs only 5.4.3.
- It does not trigger 5.3, 5.5, 7 or 8.
""",
    "5.4.venn": """
Purpose:
5.4.4 is an event-level Venn-style overlap interpretation page.

Inputs:
- Selected comparison pair.
- Significant AS definition: rMATS_FDR <= 0.05 and user-set abs(standardized dPSI) cutoff.
- Counterpart read-support threshold: user-set minimum mean IJC+SJC reads.
- Optional A-only/B-only exclusive-side read-support filter.

Outputs:
- Significant AS event Venn: A-specific, shared significant modeled events, B-specific.
- Optional filtered Venn where the A/B-specific sides are restricted to:
  - counterpart comparison has enough reads but is not significant
  - counterpart comparison has low/absent reads
  - counterpart event is not found in the other comparison
- Summary table with union count and intersection fraction.

Interpretation:
- A and B are the two comparisons in the selected comparison pair, not treatment/control groups inside one comparison.
- AS event overlap answers: which significant AS events are specific to A, shared by A/B, or specific to B?
- A-specific AS events are significant in A and checked directly for B-side counterpart read support.
- B-specific AS events are significant in B and checked directly for A-side counterpart read support.
- Shared center always means significant in both comparisons.
- The read-support filter changes only the A-only and B-only sides. It does not redefine the shared-significant center.
- "read support" means mean replicate IJC+SJC junction reads for the counterpart event in the other comparison; it is not selecting individual reads.
- The Venn page is descriptive. It is not a statistical interaction test.
- Use 5.4.3 for the read-count-aware model and FDR; use 5.4.4 to quickly explain intersections.

Run and cache behavior:
- Opening reads cached event-overlap output only.
- Run 5.4.4 Event Venn rebuilds only the selected pair's event-level Venn/overlap result.
- It does not auto-run 5.3, 5.5, 7 or 8.
""",
    "5.4.venn_multi": """
Purpose:
5.4.5 and 5.4.6 are selected-comparison direction-aware AS-event overlap pages.
5.4.5 shows only the large Venn / exact-intersection view.
5.4.6 shows only the AS-event response similarity heatmap.

Inputs:
- Two or more comparisons auto-loaded from Pairing / selected comparisons.
- The UI shows horizontal comparison checkboxes that wrap to a new line when full. No manual comparison_id/display_name typing is required.
- Strict significant AS events only: rMATS_FDR <= 0.05 and user-set abs(standardized dPSI) cutoff.
- Event identity uses event_matching_key: gene/gene_id + event_type + event_id, with coordinate fallback when event_id is unavailable.

Outputs:
- selected_comparison_event_venn_events.tsv: one row per strict significant AS event observed in at least one selected comparison.
- selected_comparison_event_venn_pairwise_summary.tsv: pairwise direction-aware A-only/shared/B-only counts, identity overlap, opposite-direction overlap and Jaccard similarity.
- 5.4.5 Direction-aware Venn / exact intersections shows only the Venn/intersection view.
- 5.4.6 Direction-aware similarity heatmap shows only the pairwise direction-aware AS-event response Jaccard similarity heatmap.
- For two selected comparisons, the Venn page shows one large A/B direction-aware exact Venn plot.
- For three selected comparisons, the Venn page shows one large three-set Venn plot with all exact intersection regions.
- For four selected comparisons, the Venn page shows one four-set Venn-style exact membership plot with all 15 intersection counts.
- For more than four selected comparisons, the Venn page shows exact membership counts as a readable intersection summary, because true area-readable Venn plots are not reliable beyond four sets.

Interpretation:
- This is direction-aware event-level AS response overlap, not gene-level candidate overlap.
- The same AS event is matched by event_matching_key, not by gene alone.
- event_matching_key is gene/gene_id + event_type + event_id; if event_id is missing, rMATS genomic coordinates are used as fallback.
- This strict key keeps multiple AS events in the same gene separate and avoids merging nearby but distinct events.
- Venn membership is direction-aware: the plotted element is event_matching_key + standardized dPSI sign.
- A-only/B-only mean significant AS event responses specific to one selected comparison after considering direction.
- shared_A_and_B means the same event_matching_key is significant in both comparisons and has the same standardized dPSI sign.
- If the same event identity is significant in multiple comparisons but the dPSI signs are opposite, it is not counted as shared in the Venn; it is flagged in the table with has_opposite_direction and summarized as opposite_direction_shared.
- Direction-aware Jaccard = same-direction shared AS event responses / union AS event responses for that comparison pair.
- identity_jaccard remains in the summary table for reference, but the heatmap uses direction-aware Jaccard.
- Every signed AS event response is counted once according to its exact membership across selected comparisons.
- This page is descriptive overlap/similarity analysis. It is not an interaction-significance test.
- Use 5.4.3 for read-count-aware interaction modeling; use 5.4.4 for one selected A/B pair; use 5.4.5 when you want a large multi-comparison Venn; use 5.4.6 when you want the multi-comparison similarity heatmap.

Run and cache behavior:
- Opening either page reads only the in-memory/cache result.
- Click Run Direction-aware Venn / Similarity to calculate only the selected-comparison overlap cache used by 5.4.5 and 5.4.6.
- It does not rebuild candidate ranking, 5.4.3 interaction modeling, sashimi, RBP or isoform modules.
""",
    "5.5.sashimi": """
Purpose:
5.5 generates event-level sashimi visualization for selected AS events.

Inputs:
- Event catalog from rMATS/candidate event follow-up.
- BAM/BAI, GTF and output folder/tool settings.
- User-selected events only.

Outputs:
- Selected-event manifest.
- Native Katana or rmats2sashimiplot outputs.
- PNG previews when available.
- Event-level failure table if a selected event fails.

Interpretation:
- The page does not auto-generate all events.
- Select one or more events, then Generate manifest or Run selected sashimi.
- Failed jobs are per-event diagnostics, not a separate analysis result.

Run and cache behavior:
- Opening reads cached event catalog and generated PNGs only.
- Refresh event list only rebuilds the event catalog.
""",
    "5.6.jutils": """
Purpose:
5.6 is Katana's native sample-level splicing pattern explorer and replaces external Jutils in the primary workflow.

Inputs:
- rMATS IncLevel1/IncLevel2 sample PSI, IJC/SJC support and event metadata from the selected comparisons.
- Supervised significant-event selection or unsupervised top-variable PSI event selection.
- User-editable FDR, abs(dPSI) and Top N display limits.

Outputs:
- Absolute and row-Z-score sample PSI heatmaps with hierarchical event/sample ordering.
- Per-comparison and combined sample PSI PCA.
- Replicate/group QC comparing within-group and between-group sample similarity.
- Sample PSI Spearman similarity together with the shared measured-event denominator.
- Significant-event Jaccard, direction agreement, opposite-direction fraction and dPSI correlation between selected comparisons.
- Descriptive PSI-trajectory clusters, reusable TSV matrices, PNG previews and parameters.json.

Interpretation:
- Heatmap/PCA answer sample-pattern and replicate-QC questions; they are not candidate ranking or differential-splicing significance tests.
- Similarity must be interpreted with shared_measured_events. A high correlation based on few events is low-confidence.
- PSI-pattern clusters use sample/replicate trajectories and are capped at a small number of descriptive modules; they are not causal programs.
- Significant-event overlap and dPSI direction agreement are comparison-level summaries and remain separate from sample-level PSI correlation.
- Candidate-gene AS heatmap belongs to 5.3.5 and uses the genes selected in 5.3.2.

Run and cache behavior:
- Opening 5.6 pages reads existing TSV/PNG cache only.
- Run Native Splicing Pattern Analysis rebuilds only 5.6 outputs; it does not trigger 5.2/5.3/5.4/5.5.
""",
    "7.rbp": """
Purpose:
RBP support checks whether candidate genes overlap known or custom RBP/splicing-factor gene sets.
It is module 7 because it is a gene-set support layer, not transcript/isoform evidence.

Inputs:
- Candidate gene selection or cached top-ranked genes.
- Built-in or custom RBP/splicing-factor gene sets.

Outputs:
- Per-comparison enrichment/support tables.
- Cross-comparison RBP support matrices and PNG plots.

Interpretation:
- This is support prioritization only.
- It does not claim causal RBP regulation.
- Cross-comparison plots show whether RBP-related candidates are shared or comparison-specific.
- Difference from 8 Isoform: 7 uses RBP/splicing-factor gene lists; 8 uses quant.sf/GTF/DTU transcript evidence.

Run and cache behavior:
- 7.1 runs only per-comparison RBP support.
- 7.2 runs only cross-comparison RBP support.
- Opening the page reads cached TSV/PNG files only.
""",
    "8.isoform": """
Purpose:
Isoform support checks whether candidate genes have transcript-level usage/expression evidence.
It is module 8 because it uses quant.sf/GTF/DTU transcript support, not RBP/splicing-factor gene sets.

Inputs:
- quant.sf files, GTF annotation and optional DTU results.
- Candidate gene selection or cached top-ranked genes.

Outputs:
- Per-comparison isoform support tables.
- Cross-comparison isoform support matrices and PNG plots.

Interpretation:
- Isoform support is downstream validation/support, not the core candidate definition.
- DEG+rMATS remain the core candidate definition.
- Cross-comparison views prioritize genes with changing isoform support patterns.
- Difference from 7 RBP: 8 evaluates transcript/isoform abundance and DTU-style support; it does not test RBP gene-set overlap.

Run and cache behavior:
- 8.1 runs only per-comparison isoform support.
- 8.2 runs only cross-comparison isoform support.
- Opening the page reads cached TSV/PNG files only.
""",
    "9.cosplice": """
Purpose:
CoSplice discovers groups of alternative-splicing events whose sample-level PSI profiles vary together across a sufficiently large multi-condition cohort. It is intended for projects such as Setd2 with many total samples; it is not fit separately inside a single n=3 comparison.

Inputs:
- Matrix source is explicit and fingerprinted:
  * rMATS event matrix uses sample-level IncLevel, IJC and SJC for compound SE/RI/A3SS/A5SS/MXE events.
  * BAM competing-junction matrix reads CIGAR N junctions, groups junctions that share a donor or acceptor, and defines usage as junction reads / total reads in that competitor cluster. IJC is the selected junction count and SJC is the cluster total minus that junction.
- BAM mode requires BAM+BAI for every sample and a project-specific GTF. It retains annotation-independent observed junctions; rMATS is then used as an annotation/differential-evidence layer rather than the matrix universe.
- Sample metadata with sample_id and condition; optional batch/covariates for a residual network.
- Stable-AS thresholds for condition coverage, junction-read support, informative PSI range and replicate repeatability.

Interpretation boundary:
- A BAM junction is not the same statistical unit as one rMATS compound event, so raw event counts from the two modes must not be compared as if they were identical.
- BAM mode broadens discovery beyond the rMATS event catalogue, but stable-AS QC must still remove low-read, low-dynamic-range and replicate-inconsistent junctions.
- Network fitting requires at least 15 total samples. Smaller projects can validate matrix construction and stable-AS QC, but should not force a module network.

Independent stages:
1. Audit resolves sample identity and file coverage.
2. Matrix builds exact event-by-sample PSI/IJC/SJC matrices.
3. Stable-AS QC removes poorly measured, read-unsupported or irreproducible events.
4. Observed network captures total biological/technical co-variation; residual network removes explicitly configured covariates first.
5. Tight/diffuse programs retain both sparse high-specificity modules and broader WGCNA-style soft-TOM modules.
6. Association tests module program scores against sample factors.
7. Module validation measures bootstrap recovery, graph-parameter sensitivity and observed/residual preservation.
8. Report summarizes existing cached outputs.

Network interpretation:
- A module is a descriptive co-splicing program: events with correlated PSI trajectories across samples.
- The default implementation uses a signed bicor sparse-kNN graph and Louvain community discovery. It is a scalable co-splicing network, not a claim of exact R WGCNA or graphical-lasso implementation.
- Program score is the first principal component of standardized PSI for events in that module.
- Tight programs use the sparse CoSplice graph; diffuse programs use signed soft adjacency and TOM. They remain separate outputs rather than forcing one partition to serve both purposes.
- Hybrid confidence combines bootstrap co-assignment, absolute program membership and replicate repeatability. External preservation contributes only when a genuinely independent cohort is configured.
- Bootstrap support resamples replicates within condition and reports how often retained edges reappear.
- The module-resolution table reports whether module number and size are sensitive to community resolution.
- Module validation resamples samples within condition, refits the network, and reports event membership stability, module Jaccard and whole-partition ARI.
- Parameter sensitivity compares kNN and community-resolution settings. Observed/residual preservation reports whether modules survive removal of configured covariates.
- Robust/moderate/unstable labels describe internal bootstrap recovery, not external biological validation.
- Association and RBP/motif/DEG/DTU support prioritize mechanisms; they do not prove causal regulation.

Run and cache behavior:
- Opening or refreshing this page reads cosplice_config.yaml, state.json and cached PNG files only.
- Every button runs exactly one stage. Missing upstream output produces a Missing input error instead of automatically running another stage.
- The command-line package provides the same stages through katana-cosplice.
""",
}


def module_detail(key: str, **values: object) -> str:
    from katana_cosplice.branding import display_module_name

    text = DETAILS.get(key, "").strip()
    if not values:
        return display_module_name(text)
    safe_values = {name: "" if value is None else str(value) for name, value in values.items()}
    return display_module_name(text.format(**safe_values))
