"""A saved single-cell junction selection goes directly to donor-resolved BAM plots."""
import hashlib
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib.figure import Figure
from matplotlib.patches import PathPatch, Rectangle
from matplotlib.path import Path as MplPath

from katana_next.science.single_cell import load_single_cell_metadata
from katana_next.science.single_cell_sashimi import junction_coordinates, regional_evidence
from katana_next.science.isoform_consequences import load_models, normalized_id
from katana_next.workflows.base import WorkflowBase, WorkflowDefinition, ParameterDefinition, ArtifactWriter
from katana_next.workflows.sources import selected_source, selection_rows, validate_source
from katana_next.workflows.selected_gene_views import save_figure
from katana_next.workflows.figure_exports import _RENDER_LOCK


class SingleCellSashimiWorkflow(WorkflowBase):
    version = "1.1.2"
    definition = WorkflowDefinition(
        workflow_id="single_cell.sashimi", workspace_id="single-cell", label="Selected-junction donor Sashimi",
        question="Do the selected junction choices have direct BAM support in each donor, condition and cell type?",
        method="Exact cluster junctions from the saved count run; indexed BAM fetch with cell barcodes and molecule deduplication; separate donor/condition/cell-type coverage tracks.",
        output_summary="PNG/SVG per selected cluster and donor page, molecule counts/usage, raw base coverage, BAM/barcode audit and coordinate-aligned GTF transcript tracks.",
        analysis_unit="donor x condition x cell type; raw coverage is read depth, arcs count junction molecules",
        required_asset_kinds=("single_cell_metadata",),
        parameters=(ParameterDefinition("selection_id", "Saved single-cell event selection", "selection", "", "Select junction rows from a completed single-cell count or differential table. All selected clusters are plotted."),
                    ParameterDefinition("flank", "Flanking bases", "number", 150, "Genomic bases outside the outermost junction.", minimum=20, maximum=2000, step=10),
                    ParameterDefinition("min_mapq", "Minimum mapping quality", "number", 20, "Applies to the local BAM evidence; lower coverage is not zero expression.", minimum=0, maximum=60, step=1)),
        limitations="Requires BAMs plus matching indexes and metadata. UMI/read-name counts and raw read coverage are different measures. No new differential significance is calculated; cells do not become biological replicates. Counts-only public datasets cannot validate this BAM bridge.",
    )

    def count_source(self, context, parameters):
        source = selected_source(context, parameters)
        queue, seen, matches = ([source] if source else []), set(), {}
        while queue:
            artifact = queue.pop()
            if artifact is None or artifact.artifact_id in seen:
                continue
            seen.add(artifact.artifact_id)
            validate_source(context, artifact)
            if artifact.analysis_id == "single_cell.quantify" and artifact.name == "Single-cell junction choices":
                matches[artifact.artifact_id] = artifact
            queue += [context.catalog.get_artifact(a) for a in artifact.parent_artifact_ids]
        if len(matches) > 1:
            raise ValueError("Selection combines multiple count runs. Choose one original source run.")
        return next(iter(matches.values()), None)

    def parent_artifacts(self, context, parameters):
        source = selected_source(context, parameters)
        counts = self.count_source(context, parameters)
        if not source or not source.analysis_id.startswith("single_cell.") or counts is None:
            raise ValueError("Save a selection from a completed single-cell count/differential run. Planning never runs quantification automatically.")
        return list({a.artifact_id: a for a in (source, counts) if a}.values())

    def input_paths(self, context, parameters):
        paths = super().input_paths(context, parameters) + [a.path for a in self.parent_artifacts(context, parameters)]
        if context.manifest.annotation.gtf:
            paths.append(context.manifest.annotation.gtf)
        metadata = load_single_cell_metadata(context)
        required = set(metadata.sample_id)
        if not required or not required.issubset({s.sample_id for s in context.manifest.samples}):
            raise ValueError("Single-cell metadata must map to declared project samples before BAM visualization")
        for sample in context.manifest.samples:
            if sample.sample_id not in required:
                continue
            bam = sample.assets.get("single_cell_bam") or sample.assets.get("bam")
            if bam is None or not bam.is_file():
                raise ValueError(f"Missing single-cell BAM for {sample.sample_id}. Counts-only inputs cannot reconstruct coverage.")
            indexes = [p for p in (Path(str(bam)+".bai"), bam.with_suffix(".bai"), Path(str(bam)+".csi")) if p.is_file()]
            if not indexes:
                raise ValueError(f"Missing BAM index for {sample.sample_id}: provide matching BAI or CSI")
            paths.extend([bam, *indexes])
        return list(dict.fromkeys(paths))

    def run(self, context, job, run_fingerprint, parameters):
        source = selected_source(context, parameters)
        counts_source = self.count_source(context, parameters)
        if not source or not source.analysis_id.startswith("single_cell.") or counts_source is None:
            raise ValueError("Save a selection from single-cell count/differential results and run single-cell quantification explicitly first.")
        selected = selection_rows(context, parameters, source, pd.read_parquet(source.path))
        if "event_id" not in selected:
            raise ValueError("Selected rows must contain exact junction event_id values.")
        counts = pd.read_parquet(counts_source.path)
        catalog = counts[["event_id", "cluster_id"]].drop_duplicates()
        if catalog.event_id.duplicated().any():
            raise ValueError("One junction maps to multiple clusters in the source counts.")
        if not set(selected.event_id).issubset(set(catalog.event_id)):
            raise ValueError("Selection and count catalog do not match; do not mix count runs.")
        cluster_ids = set(catalog[catalog.event_id.isin(selected.event_id)].cluster_id)
        if any(pd.isna(c) or not str(c).strip() for c in cluster_ids):
            raise ValueError("Cluster annotation is required to display competing junctions and their usage denominator.")
        metadata = load_single_cell_metadata(context)
        parents = tuple(a.artifact_id for a in self.parent_artifacts(context, parameters))
        writer = ArtifactWriter(context, self.definition.workflow_id, run_fingerprint, parameters, self.version,
                 common_metadata={"method": self.definition.method, "limitations": self.definition.limitations,
                                  "coverage_units": "raw aligned read depth (M/= /X only); not RPM/RPKM or UMI depth",
                                  "arc_units": "unique (sample, cell, UMI/read-name) per junction",
                                  "coordinates": "0-based half-open; XS strand matches exact junction identity, dot means unstranded",
                                  "usage_denominator": "sum of measured junction molecules within the complete source cluster"})
        gtf = context.manifest.annotation.gtf
        # Restrict model parsing to genes represented in the selected count catalog.
        selected_genes = set(counts.loc[counts.cluster_id.isin(cluster_ids), "gene_id"].dropna().map(normalized_id)) if "gene_id" in counts else set()
        models = {}
        if gtf and gtf.is_file() and selected_genes:
            from katana_next.science.isoform import transcript_gene_map
            mapping = transcript_gene_map(gtf, writer.directory / "transcript_gene_cache.parquet")
            requested = set(mapping.loc[mapping.gene_id.isin(selected_genes), "transcript_id"])
            models = load_models(gtf, requested, job.check_cancelled)
        all_counts, all_audits, all_coverage, annotation_rows = [], [], [], []
        for i, cluster in enumerate(sorted(cluster_ids)):
            job.check_cancelled()
            job.report(.10 + .75*i/max(1, len(cluster_ids)), f"Reading BAM support for cluster {i+1}/{len(cluster_ids)}")
            events = sorted(catalog.loc[catalog.cluster_id.eq(cluster), "event_id"].tolist())
            region, coverage, support, audit = regional_evidence(context.manifest.samples, metadata, events,
                            flank=int(parameters["flank"]), min_mapq=int(parameters["min_mapq"]), check_cancelled=job.check_cancelled)
            chrom, left, right, strand = region
            support["cluster_id"], audit["cluster_id"] = cluster, cluster
            all_counts.append(support)
            all_audits.append(audit)
            annotated = []
            for tx, model in models.items():
                ordered = sorted(model.exons)
                gaps = {(a[1], b[0]) for a, b in zip(ordered, ordered[1:])}
                if model.chrom == chrom and (strand == "." or model.strand == strand) and gaps & {(junction_coordinates(e)[1], junction_coordinates(e)[2]) for e in events}:
                    annotated.append((tx, model))
                    for n, (a, b) in enumerate(model.exons, 1):
                        annotation_rows.append(dict(cluster_id=cluster, transcript_id=tx, exon_number=n, chrom=chrom, start=a, end=b, strand=model.strand))
            for key, depth in coverage.items():
                frame = pd.DataFrame({"genomic_position_0based": np.arange(left, right), "read_depth": depth})
                frame["cluster_id"] = cluster
                for column, value in zip(("donor_id", "condition", "cell_type"), key):
                    frame[column] = value
                all_coverage.append(frame)
            tracks = list(coverage)
            ordered_events = sorted(events, key=lambda e: (junction_coordinates(e)[1]+junction_coordinates(e)[2])/2)
            event_pages = [ordered_events[n:n+6] for n in range(0, len(ordered_events), 6)]
            plot_pages = [(start, junction_page, visible_events) for start in range(0, len(tracks), 6)
                          for junction_page, visible_events in enumerate(event_pages)]
            for page_number, (start, junction_page, visible_events) in enumerate(plot_pages, 1):
                page = tracks[start:start+6]
                view_left = max(left, min(junction_coordinates(e)[1] for e in visible_events)-int(parameters["flank"]))
                view_right = min(right, max(junction_coordinates(e)[2] for e in visible_events)+int(parameters["flank"]))
                view = slice(view_left-left, view_right-left)
                visible_models = annotated[:6]
                with _RENDER_LOCK:
                    figure = Figure(figsize=(12, 2.1*len(page) + max(1.5, .35*len(visible_models))), layout="constrained")
                    axes = figure.subplots(len(page)+1, 1, sharex=True, squeeze=False).ravel()
                    ceiling = max(1, max((float(coverage[k].max()) for k in tracks), default=1))
                    max_count = max(1, int(support.junction_molecules.max()))
                    for axis, key in zip(axes, page):
                        depth = coverage[key]
                        axis.fill_between(np.arange(view_left, view_right), depth[view], color="#7398b2", alpha=.75)
                        axis.plot(np.arange(view_left, view_right), depth[view], color="#365e7b", lw=.6)
                        subset = support[(support.donor_id == key[0]) & (support.condition == key[1]) & (support.cell_type == key[2])].set_index("event_id")
                        for lane, event in enumerate(sorted(visible_events, key=lambda e: junction_coordinates(e)[2]-junction_coordinates(e)[1])):
                            _, a, b, _ = junction_coordinates(event)
                            count = int(subset.loc[event, "junction_molecules"])
                            if not count:
                                continue
                            ya = float(depth[max(0, a-left-1)])
                            yb = float(depth[min(len(depth)-1, b-left)])
                            apex = ceiling * (1.25 + .40*lane)
                            control = (apex - .125*(ya+yb))/.75
                            path = MplPath([(a, ya), (a, control), (b, control), (b, yb)], [1, 4, 4, 4])
                            color = "#D55E00" if event in set(selected.event_id) else "#0072B2"
                            axis.add_patch(PathPatch(path, fill=False, color=color, lw=.6+2.4*count/max_count))
                            axis.text((a+b)/2, apex+.05*ceiling, str(count), ha="center", va="bottom", fontsize=8, color=color)
                        axis.set_ylim(0, ceiling*(1.7+.4*len(visible_events)))
                        axis.set_yticks(np.linspace(0, ceiling, 3))
                        axis.set_ylabel("Read depth", fontsize=8)
                        axis.set_title(f"Donor {key[0]} | {key[1]} | {key[2]} | cells with regional reads {int(subset.cells_with_regional_reads.iloc[0])}/{int(subset.cells_in_metadata.iloc[0])}", loc="left", fontsize=9)
                        axis.spines[["top", "right"]].set_visible(False)
                    axis = axes[-1]
                    for n, (tx, model) in enumerate(visible_models):
                        exons = sorted(model.exons)
                        axis.plot([max(left, exons[0][0]), min(right, exons[-1][1])], [n, n], color="black", lw=.6)
                        for a, b in exons:
                            if min(b, right) > max(a, left):
                                axis.add_patch(Rectangle((max(a,left), n-.13), min(b,right)-max(a,left), .26, color="black"))
                        for upstream, downstream in zip(exons, exons[1:]):
                            a, b = max(left, upstream[1]), min(right, downstream[0])
                            if b <= a:
                                continue
                            positions = np.linspace(a, b, max(3, int(25*(b-a)/(right-left))))[1:-1]
                            axis.scatter(positions, np.full(len(positions), n), marker=">" if model.strand == "+" else "<", s=8, color="black")
                    if visible_models:
                        axis.set_yticks(range(len(visible_models)), [tx for tx, _ in visible_models], fontsize=6)
                        axis.set_ylim(-.5, len(visible_models)-.5)
                    else:
                        axis.text(.5, .5, "No matching annotated transcript; no isoform model inferred", ha="center", transform=axis.transAxes, fontsize=8)
                        axis.set_yticks([])
                    axis.set_xlabel(f"{chrom} genomic position (0-based) | strand {strand}")
                    axis.set_xlim(view_left, view_right)
                    axis.ticklabel_format(axis="x", style="plain", useOffset=False)
                    figure.suptitle(f"Single-cell junction evidence | {cluster} | junction page {junction_page+1}/{len(event_pages)}\nArcs: molecules; orange: selected; blue: competing. All {len(events)} cluster junctions retained in tables and usage denominator.", fontsize=10)
                    stem = hashlib.sha256(str(cluster).encode()).hexdigest()[:12]
                    save_figure(writer, figure, f"single_cell_sashimi_{stem}_{page_number}", f"Donor Sashimi: {cluster} | page {page_number}", parents)
        for name, filename, frame, entity in (
            ("Donor junction read support", "junction_support.parquet", pd.concat(all_counts, ignore_index=True), "donor_cell_type_event"),
            ("Regional BAM and barcode audit", "bam_audit.parquet", pd.concat(all_audits, ignore_index=True), "bam_audit"),
            ("Raw regional read depth", "coverage.parquet", pd.concat(all_coverage, ignore_index=True), "genomic_coverage"),
            ("Matched transcript exon coordinates", "annotation.parquet", pd.DataFrame(annotation_rows, columns=["cluster_id", "transcript_id", "exon_number", "chrom", "start", "end", "strand"]), "exon")):
            frame["row_id"] = [str(i) for i in range(len(frame))]
            writer.table(name, filename, frame, entity_type=entity, parents=parents)
        return writer.finish()
