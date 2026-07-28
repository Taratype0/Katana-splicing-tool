from __future__ import annotations

from pathlib import Path

import pandas as pd


class JutilsBenchmarkService:
    """Build the auditable native-versus-legacy Jutils comparison report."""

    def write(
        self,
        base: Path,
        output_root: Path,
        tables: dict[str, pd.DataFrame],
    ) -> list[Path]:
        legacy_root = output_root / "08_jutils"
        legacy_comparisons = self._legacy_comparisons(legacy_root)

        def shape(name: str) -> str:
            frame = tables.get(name, pd.DataFrame())
            return f"{frame.shape[0]}x{frame.shape[1]}" if not frame.empty else "missing"

        legacy_plots = legacy_root / "plots"
        legacy_pdf_count = len(list(legacy_plots.rglob("*.pdf"))) if legacy_plots.exists() else 0
        legacy_png_count = len(list(legacy_root.rglob("*.png")))
        native_plots = base / "plots"
        native_png_count = len(list(native_plots.rglob("*.png"))) if native_plots.exists() else 0
        native_tsv_count = len(tables)
        similarity_summary = self._similarity_summary(tables)

        rows = [
            {
                "benchmark_item": "validated_scope",
                "legacy_jutils": "Published/core scope: converted differential-splicing inputs, PSI heatmaps, PCA/MDS, Venn overlap and sashimi.",
                "native_jutils": "5.6 scope: sample-level PSI heatmaps, PCA, event-pattern clustering, replicate QC and quantitative comparison similarity. Katana 5.4.5 supplies direction-aware event Venn; Katana 5.5 supplies Native Sashimi.",
                "verdict": "Equivalent or stronger for splicing-pattern exploration when 5.6 is evaluated within the complete Katana workflow; not a claim of pixel-identical output or universal statistical superiority.",
            },
            {
                "benchmark_item": "input_comparisons",
                "legacy_jutils": ", ".join(legacy_comparisons) if legacy_comparisons else "not found",
                "native_jutils": f"comparison names inherited from selected rMATS inputs; dPSI matrix = {shape('dpsi_matrix')}",
                "verdict": "Compare within the same practical scope; legacy Jutils remains the compatibility reference.",
            },
            {
                "benchmark_item": "heatmap_pca_outputs",
                "legacy_jutils": f"PDF heatmap/PCA outputs: {legacy_pdf_count}; PNG preview outputs: {legacy_png_count}",
                "native_jutils": f"PNG plots: {native_png_count}; TSV result tables: {native_tsv_count}",
                "verdict": "Legacy is useful for classic Jutils PDFs; native is easier for Katana preview and downstream export.",
            },
            {
                "benchmark_item": "heatmap_result_semantics",
                "legacy_jutils": "per-comparison sample PSI heatmap/PCA using external Jutils filters; current Katana adapter calls heatmap/PCA with q<=0.05 and |dPSI|>=0.20",
                "native_jutils": "per-comparison absolute PSI heatmaps plus combined absolute/row-Z-score sample PSI matrices; selection uses the active FDR, |dPSI| and top-events parameters",
                "verdict": "Both workflows use sample-level PSI. Native additionally exposes reusable absolute and row-standardized matrices and does not treat sparse comparison-level dPSI vectors as samples.",
            },
            {
                "benchmark_item": "pca_result_semantics",
                "legacy_jutils": "per-comparison sample-level PCA PDFs from the external Jutils converted table",
                "native_jutils": f"per_comparison_sample_psi_pca.tsv = {shape('per_comparison_sample_psi_pca')}; sample_psi_pca.tsv = {shape('sample_psi_pca')}; replicate_group_separation.tsv = {shape('replicate_group_separation')}",
                "verdict": "Native keeps sample-level PSI PCA and adds an explicit within-group versus between-group replicate-similarity QC summary.",
            },
            {
                "benchmark_item": "venn_result_semantics",
                "legacy_jutils": "classic Venn diagram PDF for overlap of significant converted AS events across Jutils input comparisons",
                "native_jutils": "pairwise_as_similarity.tsv reports shared/union counts, Jaccard, signed Jaccard, direction agreement and confidence weighting; the interactive direction-aware Venn is provided by Katana 5.4.5",
                "verdict": "Katana retains a visual Venn in 5.4.5 and adds denominator-aware quantitative overlap in 5.6.3, so overlap is available in both visual and auditable table form.",
            },
            {
                "benchmark_item": "sample_level_PSI_matrix",
                "legacy_jutils": "implicit inside heatmap/PCA pipeline; not exposed as a reusable Katana TSV matrix",
                "native_jutils": f"sample_level_psi_matrix.tsv = {shape('sample_level_psi_matrix')}; sample_psi_pca.tsv = {shape('sample_psi_pca')}",
                "verdict": "Native is stronger for transparent reuse, QC and export.",
            },
            {
                "benchmark_item": "comparison_similarity_quantification",
                "legacy_jutils": "Venn/PCA visualization, but no explicit pairwise similarity table with confidence weighting",
                "native_jutils": f"pairwise_as_similarity.tsv = {shape('pairwise_as_similarity')}; {similarity_summary}",
                "verdict": "Native adds explicit AS similarity quantification and confidence weighting.",
            },
            {
                "benchmark_item": "low_overlap_warning",
                "legacy_jutils": "not explicit; high-looking correlation can be misleading if only a few events overlap",
                "native_jutils": "shared_event_support_fraction, similarity_evidence_weight, confidence_weighted_cosine_similarity, similarity_confidence",
                "verdict": "Native is safer for interpretation because low-overlap pairs are flagged instead of overinterpreted.",
            },
            {
                "benchmark_item": "event_program_summary",
                "legacy_jutils": "not exposed as a program summary table in Katana output",
                "native_jutils": f"psi_pattern_clusters.tsv = {shape('psi_pattern_clusters')}; psi_pattern_cluster_summary.tsv = {shape('psi_pattern_cluster_summary')}",
                "verdict": "Native clusters row-standardized sample/replicate PSI trajectories, not underdetermined two-value comparison dPSI vectors; clusters remain descriptive rather than causal.",
            },
            {
                "benchmark_item": "replicate_qc",
                "legacy_jutils": "sample separation can be inspected visually in PCA/heatmap",
                "native_jutils": f"replicate_psi_qc.tsv = {shape('replicate_psi_qc')}; replicate_group_separation.tsv = {shape('replicate_group_separation')}",
                "verdict": "Native exposes event-level PSI consistency and within-group versus between-group similarity; legacy relies more on visual inspection.",
            },
            {
                "benchmark_item": "dependency_and_reproducibility",
                "legacy_jutils": "external jutils.py compatibility path; output mostly PDFs and converted intermediate files",
                "native_jutils": "runs inside Katana service, writes parameters.json plus TSV/PNG outputs",
                "verdict": "Native is better integrated; legacy remains useful as compatibility/reference output.",
            },
            {
                "benchmark_item": "deliberate_behavior_differences",
                "legacy_jutils": "Supervised heatmaps can retain up to 1000 rows; --top applies only in unsupervised mode. MDS is available for non-Euclidean metrics.",
                "native_jutils": "Top N is applied per comparison in both modes for readable UI previews; PCA uses SVD; similarity adds signed and denominator-aware metrics.",
                "verdict": "These are intentional product differences. Native 5.6 is not a byte-for-byte reimplementation of every Jutils command-line option.",
            },
        ]
        benchmark = pd.DataFrame(rows)
        tsv_path = base / "native_vs_legacy_jutils_benchmark.tsv"
        md_path = base / "native_vs_legacy_jutils_benchmark.md"
        benchmark.to_csv(tsv_path, sep="\t", index=False)
        md_path.write_text(self._markdown_table(benchmark), encoding="utf-8")
        return [tsv_path, md_path]

    @staticmethod
    def _legacy_comparisons(legacy_root: Path) -> list[str]:
        manifest_path = legacy_root / "data" / "jutils_manifest.tsv"
        if not manifest_path.exists():
            return []
        try:
            manifest = pd.read_csv(manifest_path, sep="\t")
        except Exception:
            return []
        if "comparison_name" not in manifest.columns:
            return []
        return manifest["comparison_name"].dropna().astype(str).tolist()

    @staticmethod
    def _similarity_summary(tables: dict[str, pd.DataFrame]) -> str:
        similarity = tables.get("pairwise_as_similarity", pd.DataFrame())
        if similarity.empty or "confidence_weighted_cosine_similarity" not in similarity.columns:
            return "not available"
        sort_cols = ["confidence_weighted_cosine_similarity"]
        if "shared_significant_events" in similarity.columns:
            sort_cols.append("shared_significant_events")
        best = similarity.sort_values(sort_cols, ascending=False).iloc[0]

        def number(name: str) -> float:
            value = pd.to_numeric(pd.Series([best.get(name, 0)]), errors="coerce")
            return float(value.fillna(0).iloc[0])

        return (
            f"best pair {best.get('comparison_A_display_name', best.get('comparison_A', 'A'))} vs "
            f"{best.get('comparison_B_display_name', best.get('comparison_B', 'B'))}: "
            f"shared={int(number('shared_significant_events'))}, "
            f"support_fraction={number('shared_event_support_fraction'):.3f}, "
            f"confidence_weighted_cosine={number('confidence_weighted_cosine_similarity'):.3f}, "
            f"confidence={best.get('similarity_confidence', 'not_available')}"
        )

    @staticmethod
    def _markdown_table(frame: pd.DataFrame) -> str:
        def escape(value: object) -> str:
            return str(value).replace("|", "\\|").replace("\n", " ")

        headers = list(frame.columns)
        lines = [
            "# Native Jutils vs legacy Jutils benchmark",
            "",
            "This table compares cached legacy Jutils compatibility outputs with Katana Native Jutils outputs. It is a benchmark of practical output semantics, not a claim that one method is universally superior.",
            "",
            "| " + " | ".join(headers) + " |",
            "| " + " | ".join(["---"] * len(headers)) + " |",
        ]
        for row in frame.itertuples(index=False):
            lines.append("| " + " | ".join(escape(value) for value in row) + " |")
        lines.extend(
            [
                "",
                "Main conclusion: Katana preserves classic sample-PSI heatmap/PCA and visual Venn coverage, while its native explorer adds reusable matrices, measured-event denominators, explicit AS similarity, low-overlap warnings, replicate QC and PNG-first preview.",
            ]
        )
        return "\n".join(lines)
