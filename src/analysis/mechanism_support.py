from __future__ import annotations

from pathlib import Path

import pandas as pd

from src.analysis.common import clean_text, load_annotation_gene_symbol_lookup, normalize_gene_id, read_table
from src.domain.models.analysis import AnalysisRequest
from src.domain.models.project import ProjectConfig


class MechanismSupportAnalyzer:
    def run(self, project: ProjectConfig, request: AnalysisRequest) -> pd.DataFrame:
        rows: list[dict[str, object]] = []
        selected = set(request.comparison_ids) if request.comparison_ids else None
        annotation_lookup = load_annotation_gene_symbol_lookup(project.tool_paths.get("annotation_tsv"))

        for comparison in project.available_comparisons:
            if selected and comparison.comparison_id not in selected:
                continue
            deg = self._get_deg_genes(project, comparison.comparison_id, request)
            suppa = self._get_suppa_genes(project, comparison.comparison_id, request)
            dexseq = self._get_dexseq_summary(project, comparison.comparison_id, request)
            dtu = self._get_dtu_summary(project, comparison.comparison_id, request)
            symbol_lookup = self._build_gene_symbol_lookup(comparison, dexseq, dtu, annotation_lookup)
            union = deg | suppa | set(dexseq) | set(dtu)
            for gene_id in sorted(union):
                support_deg = int(gene_id in deg)
                support_suppa = int(gene_id in suppa)
                dexseq_info = dexseq.get(gene_id, {})
                dtu_info = dtu.get(gene_id, {})
                support_dexseq = int(bool(dexseq_info.get("DEXSeq_significant", False)))
                support_dtu = int(bool(dtu_info.get("DTU_significant", False)))
                rows.append(
                    {
                        "comparison_id": comparison.comparison_id,
                        "comparison_name": comparison.display_resolved_name,
                        "gene_id": gene_id,
                        "gene_symbol": (
                            annotation_lookup.get(gene_id)
                            or
                            dexseq_info.get("gene_symbol")
                            or dtu_info.get("gene_symbol")
                            or symbol_lookup.get(gene_id, gene_id)
                        ),
                        "support_DEG": support_deg,
                        "support_SUPPA": support_suppa,
                        "support_DEXSeq": support_dexseq,
                        "support_DTU": support_dtu,
                        "best_DEXSeq_qvalue": dexseq_info.get("best_DEXSeq_qvalue", pd.NA),
                        "n_DEXSeq_significant_exons": dexseq_info.get("n_DEXSeq_significant_exons", 0),
                        "best_DTU_qvalue": dtu_info.get("best_DTU_qvalue", pd.NA),
                        "n_DTU_significant_isoforms": dtu_info.get("n_DTU_significant_isoforms", 0),
                        "n_support_methods": support_deg + support_suppa + support_dexseq + support_dtu,
                    }
                )
        return pd.DataFrame(rows)

    def _build_gene_symbol_lookup(
        self,
        comparison,
        dexseq: dict[str, dict[str, object]],
        dtu: dict[str, dict[str, object]],
        annotation_lookup: dict[str, str],
    ) -> dict[str, str]:
        lookup: dict[str, str] = dict(annotation_lookup)
        if comparison.rmats_path is None:
            for source in (dexseq, dtu):
                for gene_id, payload in source.items():
                    gene_symbol = clean_text(payload.get("gene_symbol"))
                    if gene_id and gene_symbol and gene_id not in lookup:
                        lookup[gene_id] = gene_symbol
            return lookup
        rmats_post = comparison.rmats_path / "rmats_post"
        for event_type in ("SE", "RI", "A3SS", "A5SS", "MXE"):
            table = rmats_post / f"{event_type}.MATS.JC.txt"
            df = read_table(table)
            if df is None or df.empty or "GeneID" not in df.columns:
                continue
            symbols = df.get("geneSymbol", pd.Series(dtype=str)).astype(str)
            for gene_id, gene_symbol in zip(df["GeneID"].astype(str), symbols, strict=False):
                cleaned_id = normalize_gene_id(gene_id)
                cleaned_symbol = clean_text(gene_symbol)
                if cleaned_id and cleaned_symbol and cleaned_id not in lookup:
                    lookup[cleaned_id] = cleaned_symbol
        for source in (dexseq, dtu):
            for gene_id, payload in source.items():
                gene_symbol = clean_text(payload.get("gene_symbol"))
                if gene_id and gene_symbol and gene_id not in lookup:
                    lookup[gene_id] = gene_symbol
        return lookup

    def _get_deg_genes(self, project: ProjectConfig, comparison_id: str, request: AnalysisRequest) -> set[str]:
        comparison = next((item for item in project.available_comparisons if item.comparison_id == comparison_id), None)
        if comparison is None or comparison.deg_path is None:
            return set()
        df = read_table(comparison.deg_path)
        if df is None or "gene_id" not in df.columns or "padj" not in df.columns:
            return set()
        padj = pd.to_numeric(df["padj"], errors="coerce")
        out = set()
        for gene_id in df.loc[padj < request.thresholds.deg_padj, "gene_id"]:
            cleaned = normalize_gene_id(gene_id)
            if cleaned:
                out.add(cleaned)
        return out

    def _get_suppa_genes(self, project: ProjectConfig, comparison_id: str, request: AnalysisRequest) -> set[str]:
        if project.suppa_root is None:
            return set()
        comparison = next((item for item in project.available_comparisons if item.comparison_id == comparison_id), None)
        comparison_dir = self._find_named_dir(
            project.suppa_root,
            comparison.source_name("suppa") if comparison is not None else comparison_id,
        )
        if comparison_dir is None:
            return set()
        files = list(comparison_dir.rglob("*_local_diffsplice.dpsi"))
        if not files:
            return set()
        df = read_table(files[0])
        if df is None or df.empty:
            return set()

        out: set[str] = set()
        if df.shape[1] == 2:
            event_ids = df.index.astype(str)
            dpsi = pd.to_numeric(df.iloc[:, 0], errors="coerce").abs()
            pval = pd.to_numeric(df.iloc[:, 1], errors="coerce")
            keep = (pval < request.thresholds.splicing_fdr) & (dpsi >= request.thresholds.splicing_dpsi)
            for event_id in event_ids[keep.fillna(False)]:
                out.update(self._split_suppa_gene_field(event_id))
        elif df.shape[1] >= 3:
            dpsi = pd.to_numeric(df.iloc[:, 1], errors="coerce").abs()
            pval = pd.to_numeric(df.iloc[:, 2], errors="coerce")
            sub = df.loc[(pval < request.thresholds.splicing_fdr) & (dpsi >= request.thresholds.splicing_dpsi)]
            for event_id in sub.iloc[:, 0]:
                out.update(self._split_suppa_gene_field(event_id))
        return out

    def _first_symbol_value(self, subset: pd.DataFrame) -> str | None:
        for column in ("gene_symbol", "geneSymbol", "symbol", "gene"):
            if column in subset.columns:
                values = subset[column].dropna().astype(str)
                if not values.empty:
                    cleaned = clean_text(values.iloc[0])
                    if cleaned:
                        return cleaned
        return None

    def _get_dexseq_summary(self, project: ProjectConfig, comparison_id: str, request: AnalysisRequest) -> dict[str, dict[str, object]]:
        if project.dexseq_root is None:
            return {}
        comparison = next((item for item in project.available_comparisons if item.comparison_id == comparison_id), None)
        source_name = comparison.source_name("dexseq") if comparison is not None else comparison_id
        file = self._find_support_file(project.dexseq_root, source_name, [("perGeneQValue", ".csv")])
        df = read_table(file)
        if df is None or "groupID" not in df.columns or "padj" not in df.columns:
            return {}
        padj = pd.to_numeric(df["padj"], errors="coerce")
        summary: dict[str, dict[str, object]] = {}
        df = df.copy()
        df["normalized_gene_id"] = df["groupID"].map(normalize_gene_id)
        df["padj"] = padj
        if "featureID" in df.columns:
            df["featureID"] = df["featureID"].astype(str)
        grouped = df.dropna(subset=["normalized_gene_id"]).groupby("normalized_gene_id", dropna=False)
        for gene_id, subset in grouped:
            significant = subset.loc[subset["padj"] < request.thresholds.dexseq_qvalue].copy()
            if significant.empty:
                continue
            summary[str(gene_id)] = {
                "DEXSeq_significant": True,
                "best_DEXSeq_qvalue": pd.to_numeric(significant["padj"], errors="coerce").min(),
                "n_DEXSeq_significant_exons": int(len(significant)),
                "gene_symbol": self._first_symbol_value(subset),
            }
        return summary

    def _get_dtu_summary(self, project: ProjectConfig, comparison_id: str, request: AnalysisRequest) -> dict[str, dict[str, object]]:
        if project.dtu_root is None:
            return {}
        comparison = next((item for item in project.available_comparisons if item.comparison_id == comparison_id), None)
        source_name = comparison.source_name("dtu") if comparison is not None else comparison_id
        file = self._find_support_file(
            project.dtu_root,
            source_name,
            [
                ("perGeneQValue", ".tsv"),
                ("DEXSeqResults", ".tsv"),
                ("getAdjustedPValues", ".tsv"),
            ],
        )
        df = read_table(file)
        if df is None or "groupID" not in df.columns or "padj" not in df.columns:
            return {}
        df = df.copy()
        df["padj"] = pd.to_numeric(df["padj"], errors="coerce")
        df["normalized_gene_id"] = df["groupID"].map(normalize_gene_id)
        if "featureID" in df.columns:
            df["featureID"] = df["featureID"].astype(str)
        summary: dict[str, dict[str, object]] = {}
        grouped = df.dropna(subset=["normalized_gene_id"]).groupby("normalized_gene_id", dropna=False)
        for gene_id, subset in grouped:
            significant = subset.loc[subset["padj"] < request.thresholds.dtu_qvalue].copy()
            if significant.empty:
                continue
            summary[str(gene_id)] = {
                "DTU_significant": True,
                "best_DTU_qvalue": pd.to_numeric(significant["padj"], errors="coerce").min(),
                "n_DTU_significant_isoforms": int(len(significant)),
                "gene_symbol": self._first_symbol_value(subset),
            }
        return summary

    def _split_suppa_gene_field(self, event_id: object) -> list[str]:
        if pd.isna(event_id):
            return []
        gene_part = str(event_id).split(";")[0]
        return [cleaned for raw in gene_part.split("_and_") if (cleaned := clean_text(raw))]

    def _find_named_dir(self, root: Path, name: str) -> Path | None:
        for path in root.rglob("*"):
            if path.is_dir() and path.name == name:
                return path
        return None

    def _find_named_file(self, root: Path, filename: str) -> Path | None:
        matches = list(root.rglob(filename))
        return matches[0] if matches else None

    def _find_support_file(
        self,
        root: Path,
        source_name: str,
        templates: list[tuple[str, str]],
    ) -> Path | None:
        candidates = self._comparison_name_candidates(source_name)
        for candidate in candidates:
            for prefix, suffix in templates:
                direct = self._find_named_file(root, f"{prefix}.{candidate}{suffix}")
                if direct is not None:
                    return direct
                wildcard = sorted(root.rglob(f"{prefix}.*{candidate}*{suffix}"))
                if wildcard:
                    return wildcard[0]
        return None

    def _comparison_name_candidates(self, source_name: str) -> list[str]:
        values: list[str] = []
        raw = clean_text(source_name)
        if raw:
            values.append(raw)
            if "_vs_" in raw:
                values.append(raw.replace("_vs_", "-"))
            if "-vs-" in raw:
                values.append(raw.replace("-vs-", "-"))
            if " vs " in raw:
                values.append(raw.replace(" vs ", "-"))
            if "-" in raw and "_vs_" not in raw and "-vs-" not in raw and " vs " not in raw:
                values.append(raw.replace("-", "_vs_"))
        deduped: list[str] = []
        for value in values:
            if value and value not in deduped:
                deduped.append(value)
        return deduped
