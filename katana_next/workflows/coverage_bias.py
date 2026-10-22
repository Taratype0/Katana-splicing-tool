import pandas as pd

from katana_next.science.coverage_bias import read_gene_body_coverage, read_qualimap_summary, summarize_coverage
from katana_next.workflows.base import ArtifactWriter, WorkflowBase, WorkflowDefinition
from katana_next.workflows.figure_style import readable_plot as _plot


class CoverageBiasWorkflow(WorkflowBase):
    version = "1.1.1"
    definition = WorkflowDefinition(
        "quality.read_coverage_bias", "quality", "5-prime / 3-prime read coverage",
        "Are reads distributed across gene bodies, or concentrated at one end?",
        "Import RSeQC, long-form or explicitly sample-mapped Qualimap gene-body profiles; normalize each profile by its mean. Qualimap's native 100-bp bias statistic is not the 20-percent window ratio computed here.",
        "Sample-level coverage curves, 3-prime/5-prime ratios and metadata; no automatic sample exclusion.",
        "sequencing sample", required_asset_kinds=("gene_body_coverage",),
        optional_asset_kinds=("qualimap_rnaseq_metrics",),
        limitations="FastQC is not gene-body coverage. End-tagged assays have expected bias. Global coverage does not rule out local artifacts or biological alternative end usage.",
    )

    def run(self, context, job, run_fingerprint, parameters):
        assets = [a for a in context.manifest.assets if a.kind == "gene_body_coverage"]
        profiles = pd.concat([read_gene_body_coverage(a.path, a.metadata.get("sample_map"), sample_id=a.metadata.get("sample_id")) for a in assets], ignore_index=True)
        if profiles.duplicated(["sample_id", "percentile"]).any():
            raise ValueError("Coverage assets overlap; choose one profile per sample")
        metadata = pd.DataFrame([
            {"sample_id": s.sample_id, "biological_sample_id": s.biological_sample_id,
             "condition": s.condition, "batch": s.batch or "", "assay": s.attributes.get("assay", "unspecified")}
            for s in context.manifest.samples
        ])
        if metadata.empty:
            raise ValueError("Declare sample metadata before importing coverage QC")
        curves, metrics = summarize_coverage(profiles, metadata)
        oriented = curves.axis_orientation.eq("5-prime to 3-prime").all()
        writer = ArtifactWriter(context, self.definition.workflow_id, run_fingerprint, parameters, self.version,
            common_metadata={"method": self.definition.method, "limitations": self.definition.limitations,
                             "sample_exclusion": "none; diagnostic only"})
        table = writer.table("Gene-body coverage curves", "gene_body_coverage.parquet", curves, entity_type="coverage_profile")
        writer.table("Read coverage bias metrics", "coverage_bias_metrics.parquet", metrics, entity_type="sample_qc", parents=(table.artifact_id,))
        native_assets = [a for a in context.manifest.assets if a.kind == "qualimap_rnaseq_metrics"]
        if native_assets:
            native = pd.DataFrame([read_qualimap_summary(a.path, sample_id=a.metadata.get("sample_id")) for a in native_assets])
            if native.sample_id.duplicated().any() or not set(native.sample_id) <= set(profiles.sample_id):
                raise ValueError("Native metrics must map uniquely to the declared coverage samples")
            native = native.merge(metadata, on="sample_id", validate="one_to_one")
            writer.table("Native Qualimap terminal coverage metrics", "qualimap_terminal_metrics.parquet", native, entity_type="sample_qc", parents=(table.artifact_id,))
            writer.plot("Native Qualimap 5-prime / 3-prime balance", "qualimap_terminal_balance.vl.json", _plot(native.to_dict("records"),
                {"type": "circle", "size": 75}, {
                    "x": {"field": "condition", "type": "nominal", "title": "Declared condition"},
                    "y": {"field": "qualimap_five_to_three_bias", "type": "quantitative", "title": "Native Qualimap 5-prime / 3-prime bias", "scale": {"zero": False}},
                    "color": {"field": "condition", "type": "nominal"},
                    "tooltip": [{"field": c} for c in ("sample_id", "condition", "qualimap_five_prime_bias", "qualimap_three_prime_bias", "qualimap_five_to_three_bias", "bias_window_bp", "bias_top_transcripts")],
                }, title="Terminal balance must be interpreted alongside the full coverage shape", description="Preserved native Qualimap metrics. A ratio near one does not prove uniform coverage, RNA integrity or absence of local artifacts."), entity_type="coverage_plot")
        writer.plot("5-prime to 3-prime gene-body coverage", "coverage_bias.vl.json", _plot(curves.to_dict("records"), "line", {
            "x": {"field": "percentile", "type": "quantitative", "title": "Gene body: 5-prime (0) to 3-prime (100)" if oriented else "Source transcript position (0 to 100); native end-bias metrics are separate"},
            "y": {"field": "relative_coverage", "type": "quantitative", "title": "Coverage / sample mean"},
            "detail": {"field": "sample_id", "type": "nominal"},
            "color": {"field": "condition", "type": "nominal"},
            "tooltip": [{"field": c} for c in ("sample_id", "condition", "batch", "percentile", "coverage", "relative_coverage")],
        }, title="Read coverage along gene bodies", description=self.definition.limitations), entity_type="coverage_plot", parents=(table.artifact_id,))
        return writer.finish()
