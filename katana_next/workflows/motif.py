"""FASTA-to-PWM support for explicitly selected exact events."""
from pathlib import Path
import pandas as pd

from katana_cosplice.motif_scanner import IndexedFasta, _event_windows, load_motif_database, motif_database_audit, scan_event_motifs
from katana_next.science.bulk import EVENT_COORDINATES
from katana_next.workflows.base import WorkflowBase, WorkflowDefinition, ParameterDefinition, ArtifactWriter
from katana_next.workflows.sources import selected_source, selection_rows
from katana_next.workflows.bulk import _plot


class EventMotifWorkflow(WorkflowBase):
    version = "2.4.1"
    definition = WorkflowDefinition(
        workflow_id="prioritize.event_motifs", workspace_id="prioritize", label="Selected-event RBP motif scan",
        question="Which species-matched RBP motifs occur around the splice boundaries of selected events?",
        method="Strand-oriented FIMO site p-values and gene-independent GC/length-matched enrichment; optional explicitly exploratory relative-score scanning.",
        output_summary="Event-window extraction audit, motif catalog, exact selected-event hit table and RBP support plot.",
        analysis_unit="sequence window and PWM match, not regulatory causality",
        parameters=(ParameterDefinition("selection_id", "Saved AS event selection", "selection", "", "Exact DAS events with declared coordinates. Every selected event is audited."),
                    ParameterDefinition("mode", "Motif evidence mode", "select", "calibrated", "Calibrated requires a FIMO runtime asset. Exploratory is not a significance test.", ("calibrated", "exploratory")),
                    ParameterDefinition("motif_evidence", "Motif catalog evidence", "select", "all_evidence_separate", "Experimental, homology-inferred and model-predicted motifs are tested separately, never merged into stronger binding evidence.", ("direct_experiment", "all_evidence_separate")),
                    ParameterDefinition("site_pvalue", "FIMO site p-value threshold", "number", .0001, "Fixed site screening threshold for both foreground and controls. Global conservative site BH and paired enrichment BH are separate outputs.", minimum=.00000001, maximum=.01, step=.0001),
                    ParameterDefinition("window_bp", "Window on each side of boundary (bp)", "number", 100, "Merged genomic windows around each event exon boundary; minus-strand sequences are reverse-complemented.", minimum=10, maximum=1000, step=10),
                    ParameterDefinition("score_threshold", "Minimum relative PWM score", "number", .85, "Exploratory normalized score, NOT a p-value or FDR. Does not claim FIMO equivalence.", minimum=0, maximum=1, step=.01)),
        required_asset_kinds=("rbp_motif_database",),
        optional_asset_kinds=("fimo_runtime",),
        presentation_role="advanced",
        limitations="Sequence matches and matched enrichment prioritize RBPs, not causal binding. Enrichment uses one exact event per gene with unique matched control genes; controls are exact-type and within 5 percentage points GC and 20% length. Unmatched genes are reported, not silently counted as no hits. Exploratory mode provides no calibrated significance.",
    )

    def parent_artifacts(self, context, parameters):
        source = selected_source(context, parameters)
        return [source] if source else []

    def input_paths(self, context, parameters):
        paths = [*super().input_paths(context, parameters), *(a.path for a in self.parent_artifacts(context, parameters))]
        fasta = context.manifest.annotation.genome_fasta
        if fasta:
            paths.extend([fasta, Path(str(fasta) + ".fai")])
        if parameters.get("mode", "calibrated") == "calibrated":
            paths.append(Path(__file__).resolve().parents[2]/"scripts"/"fimo_worker.py")
        return paths

    def run(self, context, job, run_fingerprint, parameters):
        source = selected_source(context, parameters)
        if not source:
            raise ValueError("Select exact AS events first. Motif scan never selects events implicitly.")
        events = selection_rows(context, parameters, source, pd.read_parquet(source.path)).copy()
        if "comparison_id" not in events or events.comparison_id.nunique() != 1:
            raise ValueError("Select events from exactly one comparison. Pooled-condition motif inference is not supported.")
        if parameters.get("mode", "calibrated") == "calibrated":
            return self._run_calibrated(context, job, run_fingerprint, parameters, source, events)
        required = {"event_id", "chr", "strand"}
        if not required.issubset(events.columns):
            raise ValueError("Motif scan requires exact AS event coordinates, chromosome and strand.")
        fasta_path = context.manifest.annotation.genome_fasta
        if not fasta_path or not fasta_path.is_file():
            raise ValueError("Configure genome FASTA and its .fai index; transcript FASTA is not interchangeable.")
        fasta = IndexedFasta(fasta_path)
        events = events.drop_duplicates("event_id")
        events["event_key"] = events.event_id
        coordinate_cols = [c for c in EVENT_COORDINATES if c not in {"chr", "strand"} and c in events]
        events["event_coordinates"] = events[coordinate_cols].apply(lambda r: "|".join(str(int(v)) for v in r if pd.notna(v)), axis=1)
        audit = []
        for index, row in events.iterrows():
            chrom = str(row.chr)
            mapped = chrom if chrom in fasta.index else chrom.removeprefix("chr") if chrom.startswith("chr") else "chr" + chrom
            if mapped not in fasta.index:
                raise ValueError(f"Chromosome {chrom} is absent from the declared FASTA; do not reinterpret this as no motif hits.")
            if row.strand not in {"+", "-"}:
                raise ValueError(f"Invalid strand for event {row.event_id}")
            events.loc[index, "chr"] = mapped
            windows = _event_windows(row, int(parameters["window_bp"]))
            if not windows:
                raise ValueError(f"No splice boundary coordinates: {row.event_id}")
            for start, end in windows:
                sequence = fasta.fetch(mapped, start, end)
                audit.append({"row_id": f"{row.event_id}::{start}-{end}", "event_id": row.event_id,
                              "source_chromosome": chrom, "fasta_chromosome": mapped, "start_0based": start,
                              "end_exclusive": min(end, fasta.index[mapped][0]), "strand": row.strand,
                              "sequence_length": len(sequence), "ambiguous_bases": sum(b not in "ACGT" for b in sequence)})
        databases = [a for a in context.manifest.assets if a.kind == "rbp_motif_database"]
        hits, catalog = [], []
        for i, database in enumerate(databases):
            organism = database.metadata.get("species", context.manifest.metadata.get("species", ""))
            if not organism:
                raise ValueError("Declare motif database species explicitly; cross-species matching is not automatic.")
            if organism == "Homo_sapiens" and context.manifest.annotation.genome_build in {"GRCm39", "GRCm38", "mm10", "mm39"}:
                raise ValueError("Human motifs cannot be silently substituted for the mouse motif catalog.")
            motifs = load_motif_database(database.path, organism=organism)
            if not motifs:
                raise ValueError(f"No motifs match {organism} in {database.asset_id}")
            catalog.extend({"row_id": f"{database.asset_id}::{m.rbp}::{m.motif_id}", "RBP": m.rbp,
                            "motif_id": m.motif_id, "species": organism, "source": m.source} for m in motifs)
            job.report(.15 + .6 * i / len(databases), f"Scanning {len(events)} exact events against {len(motifs)} motifs")
            part = scan_event_motifs(events, fasta_path, database.path, window_bp=int(parameters["window_bp"]),
                                     score_threshold=float(parameters["score_threshold"]), organism=organism)
            part["database"] = database.asset_id
            hits.append(part)
        result = pd.concat(hits, ignore_index=True)
        result["row_id"] = result.database + "::" + result.event_key + "::" + result.RBP + "::" + result.motif_id
        result = result.merge(events[["event_key", "gene_id", "gene_symbol"]], on="event_key", how="left", validate="many_to_one")
        writer = ArtifactWriter(context, self.definition.workflow_id, run_fingerprint, parameters, self.version,
                                common_metadata={"method": self.definition.method, "limitations": self.definition.limitations,
                                                 "selected_event_count": len(events), "relative_score_not_significance": True,
                                                 "attract_catalog_filter": "species matched, non-mutated, Ensembl gene ID for mouse/human; structure-only identifiers excluded"})
        parents = (source.artifact_id,)
        writer.table("Motif extraction windows", "motif_windows.parquet", pd.DataFrame(audit), entity_type="splice_event", parents=parents)
        writer.table("Species-matched motif catalog", "motif_catalog.parquet", pd.DataFrame(catalog), entity_type="motif", parents=parents)
        writer.table("Selected-event motif matches", "event_motif_hits.parquet", result, entity_type="motif_hit", parents=parents)
        summary = result.groupby("RBP", as_index=False).agg(events_with_matches=("event_key", "nunique"), sites=("site_count", "sum"))
        summary["row_id"] = summary.RBP
        writer.table("RBP motif support summary", "motif_summary.parquet", summary, entity_type="rbp", parents=parents)
        writer.plot("RBP motif matches in selected events", "motif_support.vl.json", _plot(summary.to_dict("records"), "bar",
                    {"x": {"field": "events_with_matches", "type": "quantitative", "title": "Selected events with sequence matches"},
                     "y": {"field": "RBP", "type": "nominal", "sort": "-x"}, "color": {"value": "#0072B2"},
                     "tooltip": [{"field": "RBP"}, {"field": "events_with_matches"}, {"field": "sites"}]},
                    title="Putative RBP binding motifs", description="Exploratory sequence support; not binding validation or enrichment significance.", height=max(300, 22 * len(summary))), entity_type="motif_plot", parents=parents)
        return writer.finish()

    def _run_calibrated(self, context, job, run_fingerprint, parameters, source, events):
        from katana_next.science.motif_significance import extract_windows, match_background, calibrated_scan, matched_enrichment
        runtimes = [a for a in context.manifest.assets if a.kind == "fimo_runtime"]
        if len(runtimes) != 1:
            raise ValueError("Calibrated motif analysis requires one configured FIMO runtime resource. No silent fallback to relative scores.")
        fasta = context.manifest.annotation.genome_fasta
        if not fasta or not fasta.is_file():
            raise ValueError("Configure genome FASTA and .fai for calibrated motifs")
        comparison = events.comparison_id.iloc[0]
        universe = pd.read_parquet(source.path)
        universe = universe.loc[universe.comparison_id.eq(comparison)].drop_duplicates("event_id")
        from katana_next.science.motif_significance import measurable_event_background
        universe, eligibility = measurable_event_background(universe)
        if not {"event_id", "chr", "strand", "gene_id", "event_type"}.issubset(universe.columns):
            raise ValueError("Motif significance requires an exact AS-event source with gene IDs and coordinates")
        selected = set(events.event_id)
        if not selected <= set(universe.event_id):
            raise ValueError("Selected events are outside the measurable-event background")
        job.report(.05, "Extracting the declared tested-event background; selection is not the universe")
        def progress(done, total):
            job.check_cancelled()
            job.report(.05 + .25*done/max(1, total), f"Extracting motif background: {done}/{total} exact events")
        windows, features = extract_windows(universe, fasta, parameters["window_bp"], progress=progress)
        pairs = match_background(features, selected)
        gene_symbols = universe.drop_duplicates("gene_id").set_index("gene_id").gene_symbol
        for role in ("foreground", "background"):
            pairs[f"{role}_gene_symbol"] = pairs[f"{role}_gene"].map(gene_symbols)
        controls = set(pairs.loc[pairs.status.eq("matched"), "background_event"])
        if not controls:
            raise ValueError("No independent GC/length/type-matched control genes. Select events from an unfiltered tested-event catalog, or use explicitly exploratory mode.")
        windows = windows[windows.event_id.isin(selected | controls)].copy()
        windows["role"] = windows.event_id.map(lambda event: "foreground" if event in selected else "background")
        motifs = []
        audits = []
        for asset in context.manifest.assets:
            if asset.kind != "rbp_motif_database":
                continue
            organism = asset.metadata.get("species", context.manifest.metadata.get("species", ""))
            if not organism:
                raise ValueError("Explicit motif species is required")
            mouse = context.manifest.annotation.genome_build in {"GRCm39", "GRCm38", "mm10", "mm39"}
            if mouse and organism != "Mus_musculus":
                raise ValueError("Motif catalog species does not match the mouse genome")
            if context.manifest.annotation.genome_build in {"GRCh38", "GRCh37", "hg19", "hg38"} and organism != "Homo_sapiens":
                raise ValueError("Motif catalog species does not match the human genome")
            motifs.extend(load_motif_database(asset.path, organism=organism))
            audit = motif_database_audit(asset.path, organism=organism)
            audit["database"] = asset.asset_id
            audits.append(audit)
        if parameters.get("motif_evidence") == "direct_experiment":
            motifs = [m for m in motifs if m.evidence_class == "direct_experiment"]
        if not motifs:
            raise ValueError("No species-matched PWM motifs available")
        writer = ArtifactWriter(context, self.definition.workflow_id, run_fingerprint, parameters, self.version,
                                common_metadata={"method": self.definition.method, "limitations": self.definition.limitations,
                                                 "background_source": source.artifact_id, "selected_events": len(selected),
                                                 "matched_independent_genes": len(controls), "testing_unit": "one representative event per gene",
                                                 "PWM_probability_convention": "MEME probability matrix nsites=20 for smoothing, not measured observations; pseudocount=0.1"})
        job.report(.4, f"MEME FIMO: {len(motifs)} motifs, {len(windows)} strand-oriented windows")
        sites, catalog, diagnostics = calibrated_scan(motifs, windows, writer.directory/"fimo", runtimes[0].path, parameters["site_pvalue"])
        enrichment = matched_enrichment(sites, catalog, pairs)
        enrichment["comparison_id"] = comparison
        enrichment["FDR_family"] = "all tested RBP-by-motif-evidence classes within this comparison"
        pairs["comparison_id"] = comparison
        sites["comparison_id"] = comparison
        parents = (source.artifact_id,)
        writer.table("Selected motif input events and genes", "motif_input_events.parquet",
                     events.drop_duplicates("event_id"), entity_type="splice_event", parents=parents)
        for label, name, table, entity in (
            ("FIMO calibrated sites", "fimo_sites.parquet", sites, "motif_hit"),
            ("Matched-background motif enrichment", "motif_enrichment.parquet", enrichment, "rbp"),
            ("Foreground and control matching audit", "motif_matching.parquet", pairs, "splice_event"),
            ("Tested-event sequence background", "motif_background.parquet", features, "splice_event"),
            ("Measurable-event eligibility", "motif_eligibility.parquet", eligibility, "splice_event"),
            ("Species-matched motif catalog", "motif_catalog.parquet", catalog, "motif"),
            ("Native motif database availability audit", "motif_database_audit.parquet", pd.concat(audits, ignore_index=True), "motif"),
        ):
            table = table.copy()
            table["row_id"] = table.index.astype(str)
            writer.table(label, name, table, entity_type=entity, parents=parents)
        import json
        (writer.directory/"fimo"/"diagnostics.json").write_text(json.dumps(diagnostics, indent=2), encoding="utf-8")
        for path in (writer.directory/"fimo").iterdir():
            writer.existing_file(f"FIMO provenance: {path.name}", path, entity_type="motif_provenance", parents=parents)
        display = enrichment.sort_values(["enrichment_fdr", "enrichment_pvalue", "RBP"], na_position="last").head(30).copy()
        display["RBP_evidence"] = display.RBP + " / " + display.evidence_class
        writer.plot("RBP motif enrichment with matched background", "motif_enrichment.vl.json", _plot(display.to_dict("records"), "circle",
                    {"x": {"field": "hit_fraction_difference", "type": "quantitative", "title": "Foreground minus matched-background hit fraction"},
                     "y": {"field": "RBP_evidence", "type": "nominal", "sort": None, "title": "RBP / motif evidence"},
                     "color": {"field": "enrichment_fdr", "type": "quantitative"},
                     "size": {"field": "foreground_hit_genes", "type": "quantitative", "title": "Selected hit genes"},
                     "tooltip": [{"field": c} for c in ("RBP", "evidence_class", "matched_gene_pairs", "foreground_hit_genes", "background_hit_genes", "enrichment_pvalue", "enrichment_fdr")]},
                    title=f"{comparison}: RBP motif enrichment", description="Top 30 by FDR/P for display only; all tested RBPs are in the full table. Exact paired test, one event per gene; motif evidence classes are separate.", height=max(300, len(display)*24)), entity_type="motif_plot", parents=parents)
        return writer.finish()
