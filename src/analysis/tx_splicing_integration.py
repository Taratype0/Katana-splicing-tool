from __future__ import annotations

from dataclasses import asdict
import json
import math
import re
from pathlib import Path

import pandas as pd

from src.analysis.common import EVENT_TYPES, clean_text, event_uid_from_row, normalize_gene_id, read_table
from src.analysis.stable_as import apply_stable_as_filter
from src.domain.models.analysis import AnalysisRequest
from src.domain.models.project import ProjectConfig


class TxSplicingIntegrator:
    @staticmethod
    def _safe_bool(value) -> bool:
        if value is None:
            return False
        try:
            if pd.isna(value):
                return False
        except (TypeError, ValueError):
            pass
        if isinstance(value, str):
            return value.strip().lower() in {"1", "true", "yes", "y", "t"}
        return bool(value)

    def run(self, project: ProjectConfig, request: AnalysisRequest) -> tuple[pd.DataFrame, pd.DataFrame]:
        gene_frames: list[pd.DataFrame] = []
        summary_rows: list[dict[str, object]] = []
        selected = set(request.comparison_ids) if request.comparison_ids else None

        for comparison in project.available_comparisons:
            if selected and comparison.comparison_id not in selected:
                continue
            gene_level = self._read_cached_gene_level(project, comparison, request)
            if gene_level is None:
                deg_df = self._read_deg(comparison, request)
                rmats_df = self._read_rmats_all(comparison, request)
                gene_level = self._build_gene_level_summary(
                    comparison.comparison_id,
                    comparison.display_resolved_name,
                    deg_df,
                    rmats_df,
                )
                self._write_cached_gene_level(project, comparison, request, gene_level)
            if gene_level.empty:
                continue
            gene_frames.append(gene_level)
            summary_rows.append(
                {
                    "comparison_id": comparison.comparison_id,
                    "comparison_name": comparison.display_resolved_name,
                    "n_gene_rows": len(gene_level),
                    "n_deg_sig": int(gene_level["deg_sig"].fillna(False).sum()),
                    "n_as_sig": int(gene_level["splicing_sig"].fillna(False).sum()),
                    "n_both": int((gene_level["combined_class"] == "both").sum()),
                    "n_de_only": int((gene_level["combined_class"] == "DE-only").sum()),
                    "n_as_only": int((gene_level["combined_class"] == "AS-only").sum()),
                    "n_neither": int((gene_level["combined_class"] == "neither").sum()),
                }
            )

        gene_table = pd.concat(gene_frames, axis=0, ignore_index=True) if gene_frames else pd.DataFrame()
        summary = pd.DataFrame(summary_rows)
        return gene_table, summary

    def _read_deg(self, comparison, request: AnalysisRequest) -> pd.DataFrame:
        if comparison.deg_path is None or not comparison.deg_path.exists():
            return pd.DataFrame(columns=["gene_id", "log2FC", "deg_padj", "deg_sig"])
        try:
            df = pd.read_csv(
                comparison.deg_path,
                sep="\t",
                low_memory=False,
                usecols=lambda column: column in {"gene_id", "log2FoldChange", "padj"},
            )
        except Exception:
            df = read_table(comparison.deg_path)
        if df is None or "gene_id" not in df.columns:
            return pd.DataFrame(columns=["gene_id", "log2FC", "deg_padj", "deg_sig"])
        df = df.copy()
        df["gene_id"] = df["gene_id"].map(normalize_gene_id)
        df["log2FC"] = pd.to_numeric(df["log2FoldChange"], errors="coerce")
        if comparison.reverse_deg:
            df["log2FC"] = -df["log2FC"]
        df["deg_padj"] = pd.to_numeric(df["padj"], errors="coerce")
        df["deg_sig"] = (df["deg_padj"] <= request.thresholds.deg_padj) & (
            df["log2FC"].abs() > request.thresholds.deg_log2fc
        )
        df = df.dropna(subset=["gene_id"]).copy()
        return df[["gene_id", "log2FC", "deg_padj", "deg_sig"]].copy()

    def _read_rmats_all(self, comparison, request: AnalysisRequest) -> pd.DataFrame:
        if comparison.rmats_path is None:
            return pd.DataFrame()
        rows = []
        rmats_dir = comparison.rmats_path / "rmats_post"
        for event_type in EVENT_TYPES:
            fp = rmats_dir / f"{event_type}.MATS.{request.rmats_mode}.txt"
            if not fp.exists():
                continue
            needed = {
                "GeneID",
                "geneSymbol",
                "FDR",
                "IncLevelDifference",
                "IncLevel1",
                "IncLevel2",
                "IJC_SAMPLE_1",
                "IJC_SAMPLE_2",
                "SJC_SAMPLE_1",
                "SJC_SAMPLE_2",
                *self._event_uid_columns(),
            }
            df = pd.read_csv(
                fp,
                sep="\t",
                low_memory=False,
                usecols=lambda column, needed=needed: column in needed,
            )
            if df.empty:
                continue
            df["comparison_id"] = comparison.comparison_id
            df["event_type"] = event_type
            df["gene_id"] = df["GeneID"].map(normalize_gene_id)
            df["geneSymbol"] = df.get("geneSymbol", pd.Series(dtype=str)).map(clean_text)
            df["event_FDR"] = pd.to_numeric(df["FDR"], errors="coerce")
            df["dPSI"] = pd.to_numeric(df["IncLevelDifference"], errors="coerce")
            if comparison.reverse_splicing:
                df["dPSI"] = -df["dPSI"]
            df["abs_dPSI"] = df["dPSI"].abs()
            evaluated = apply_stable_as_filter(
                df,
                fdr=df["event_FDR"],
                standardized_dpsi=df["dPSI"],
                reverse_splicing=comparison.reverse_splicing,
                thresholds=request.thresholds,
            )
            for column in evaluated.columns:
                df[column] = evaluated[column]
            df["splicing_sig"] = df["stable_as_event"]
            df["event_uid"] = pd.NA
            sig_mask = df["splicing_sig"].fillna(False)
            if sig_mask.any():
                df.loc[sig_mask, "event_uid"] = self._event_uid_series(df.loc[sig_mask], event_type)
            rows.append(
                df[
                    [
                        "comparison_id",
                        "gene_id",
                        "geneSymbol",
                        "event_type",
                        "event_uid",
                        "event_FDR",
                        "dPSI",
                        "abs_dPSI",
                        "rmats_sig_before_stability",
                        "stable_as_event",
                        "mean_junction_reads_group1",
                        "mean_junction_reads_group2",
                        "replicate_count_group1",
                        "replicate_count_group2",
                        "supported_replicate_count_group1",
                        "supported_replicate_count_group2",
                        "replicate_direction_consistency",
                        "stable_support_status",
                        "splicing_sig",
                    ]
                ].copy()
            )
        return pd.concat(rows, axis=0, ignore_index=True) if rows else pd.DataFrame()

    def _build_gene_level_summary(
        self,
        comparison_id: str,
        comparison_name: str,
        deg_df: pd.DataFrame,
        rmats_df: pd.DataFrame,
    ) -> pd.DataFrame:
        deg_working = deg_df.dropna(subset=["gene_id"]).copy() if deg_df is not None and not deg_df.empty else pd.DataFrame()
        if not deg_working.empty:
            deg_working["gene_id"] = deg_working["gene_id"].astype(str)
            deg_working = deg_working.drop_duplicates(subset=["gene_id"], keep="first")

        rmats_working = rmats_df.dropna(subset=["gene_id"]).copy() if rmats_df is not None and not rmats_df.empty else pd.DataFrame()
        if not rmats_working.empty:
            rmats_working["gene_id"] = rmats_working["gene_id"].astype(str)

        deg_genes = set(deg_working["gene_id"].astype(str).tolist()) if not deg_working.empty else set()
        splicing_genes = set(rmats_working["gene_id"].astype(str).tolist()) if not rmats_working.empty else set()
        all_genes = sorted(deg_genes | splicing_genes)
        if not all_genes:
            return pd.DataFrame()

        result = pd.DataFrame({"gene_id": all_genes})
        result["comparison_id"] = comparison_id
        result["comparison_name"] = comparison_name

        if not deg_working.empty:
            result = result.merge(
                deg_working[["gene_id", "log2FC", "deg_padj", "deg_sig"]],
                on="gene_id",
                how="left",
            )
        else:
            result["log2FC"] = math.nan
            result["deg_padj"] = math.nan
            result["deg_sig"] = False

        if not rmats_working.empty:
            symbol_lookup = (
                rmats_working.dropna(subset=["geneSymbol"])
                .drop_duplicates(subset=["gene_id"], keep="first")
                .set_index("gene_id")["geneSymbol"]
            )
            result["geneSymbol"] = result["gene_id"].map(symbol_lookup).fillna(result["gene_id"])

            rmats_sig = rmats_working.loc[rmats_working["splicing_sig"].fillna(False)].copy()
            if not rmats_sig.empty:
                sig_counts = rmats_sig.groupby("gene_id", sort=False).size().rename("n_sig_events_all")
                sig_types = rmats_sig.groupby("gene_id", sort=False)["event_type"].agg(
                    lambda values: ";".join(sorted(set(values.dropna().astype(str))))
                ).rename("sig_event_types")
                representative = (
                    rmats_sig.sort_values(["gene_id", "abs_dPSI", "event_FDR"], ascending=[True, False, True], na_position="last")
                    .drop_duplicates(subset=["gene_id"], keep="first")
                    .rename(
                        columns={
                            "event_uid": "representative_event_id",
                            "event_type": "representative_event_type",
                            "dPSI": "representative_dPSI",
                            "event_FDR": "representative_event_FDR",
                            "mean_junction_reads_group1": "representative_mean_junction_reads_group1",
                            "mean_junction_reads_group2": "representative_mean_junction_reads_group2",
                            "replicate_count_group1": "representative_replicate_count_group1",
                            "replicate_count_group2": "representative_replicate_count_group2",
                            "supported_replicate_count_group1": "representative_supported_replicate_count_group1",
                            "supported_replicate_count_group2": "representative_supported_replicate_count_group2",
                            "replicate_direction_consistency": "representative_replicate_direction_consistency",
                            "stable_support_status": "representative_stable_support_status",
                        }
                    )
                )
                result = result.merge(sig_counts, on="gene_id", how="left")
                result = result.merge(sig_types, on="gene_id", how="left")
                result = result.merge(
                    representative[
                        [
                            "gene_id",
                            "representative_event_id",
                            "representative_event_type",
                            "representative_dPSI",
                            "representative_event_FDR",
                            "representative_mean_junction_reads_group1",
                            "representative_mean_junction_reads_group2",
                            "representative_replicate_count_group1",
                            "representative_replicate_count_group2",
                            "representative_supported_replicate_count_group1",
                            "representative_supported_replicate_count_group2",
                            "representative_replicate_direction_consistency",
                            "representative_stable_support_status",
                        ]
                    ],
                    on="gene_id",
                    how="left",
                )
            else:
                result["n_sig_events_all"] = 0
                result["sig_event_types"] = ""
                result["representative_event_id"] = pd.NA
                result["representative_event_type"] = pd.NA
                result["representative_dPSI"] = math.nan
                result["representative_event_FDR"] = math.nan
                result["representative_mean_junction_reads_group1"] = math.nan
                result["representative_mean_junction_reads_group2"] = math.nan
                result["representative_replicate_count_group1"] = math.nan
                result["representative_replicate_count_group2"] = math.nan
                result["representative_supported_replicate_count_group1"] = math.nan
                result["representative_supported_replicate_count_group2"] = math.nan
                result["representative_replicate_direction_consistency"] = math.nan
                result["representative_stable_support_status"] = pd.NA
        else:
            result["geneSymbol"] = result["gene_id"]
            result["n_sig_events_all"] = 0
            result["sig_event_types"] = ""
            result["representative_event_id"] = pd.NA
            result["representative_event_type"] = pd.NA
            result["representative_dPSI"] = math.nan
            result["representative_event_FDR"] = math.nan
            result["representative_mean_junction_reads_group1"] = math.nan
            result["representative_mean_junction_reads_group2"] = math.nan
            result["representative_replicate_count_group1"] = math.nan
            result["representative_replicate_count_group2"] = math.nan
            result["representative_supported_replicate_count_group1"] = math.nan
            result["representative_supported_replicate_count_group2"] = math.nan
            result["representative_replicate_direction_consistency"] = math.nan
            result["representative_stable_support_status"] = pd.NA

        for column in ("n_sig_events_all",):
            result[column] = pd.to_numeric(result.get(column), errors="coerce").fillna(0).astype(int)
        result["sig_event_types"] = result.get("sig_event_types", pd.Series("", index=result.index)).fillna("").astype(str)
        result["deg_sig"] = result.get("deg_sig", pd.Series(False, index=result.index)).map(self._safe_bool)
        result["splicing_sig"] = result["n_sig_events_all"] > 0
        result["geneSymbol"] = result["geneSymbol"].fillna(result["gene_id"])

        classes = result.apply(
            lambda row: self._classify_direction(
                row.get("log2FC") if pd.notna(row.get("log2FC")) else 0.0,
                row.get("representative_dPSI") if pd.notna(row.get("representative_dPSI")) else 0.0,
                self._safe_bool(row.get("deg_sig")),
                self._safe_bool(row.get("splicing_sig")),
            ),
            axis=1,
        )
        result["direction_class"] = [item[0] for item in classes]
        result["agreement_class"] = [item[1] for item in classes]
        result["combined_class"] = result["direction_class"].map(self._combined_class)
        return result[
            [
                "comparison_id",
                "comparison_name",
                "gene_id",
                "geneSymbol",
                "log2FC",
                "deg_padj",
                "deg_sig",
                "representative_event_id",
                "representative_event_type",
                "representative_dPSI",
                "representative_event_FDR",
                "representative_mean_junction_reads_group1",
                "representative_mean_junction_reads_group2",
                "representative_replicate_count_group1",
                "representative_replicate_count_group2",
                "representative_supported_replicate_count_group1",
                "representative_supported_replicate_count_group2",
                "representative_replicate_direction_consistency",
                "representative_stable_support_status",
                "splicing_sig",
                "combined_class",
                "direction_class",
                "agreement_class",
                "n_sig_events_all",
                "sig_event_types",
            ]
        ].copy()

    @staticmethod
    def _event_uid_columns() -> list[str]:
        return [
            "chr",
            "strand",
            "exonStart_0base",
            "exonEnd",
            "upstreamES",
            "upstreamEE",
            "downstreamES",
            "downstreamEE",
            "riExonStart_0base",
            "riExonEnd",
            "longExonStart_0base",
            "longExonEnd",
            "shortES",
            "shortEE",
            "flankingES",
            "flankingEE",
            "1stExonStart_0base",
            "1stExonEnd",
            "2ndExonStart_0base",
            "2ndExonEnd",
        ]

    def _event_uid_series(self, frame: pd.DataFrame, event_type: str) -> pd.Series:
        return frame.apply(lambda row: event_uid_from_row(row, event_type), axis=1)

    def _read_cached_gene_level(self, project: ProjectConfig, comparison, request: AnalysisRequest) -> pd.DataFrame | None:
        cache_dir = self._comparison_cache_dir(project, comparison.comparison_id)
        if cache_dir is None:
            return None
        meta_path = cache_dir / "tx_splicing_cache_meta.json"
        table_path = cache_dir / "tx_splicing_gene_table.tsv"
        if not meta_path.exists() or not table_path.exists():
            return None
        try:
            meta = json.loads(meta_path.read_text(encoding="utf-8"))
            if meta.get("signature") != self._comparison_cache_signature(comparison, request):
                return None
            frame = pd.read_csv(table_path, sep="\t", low_memory=False)
        except Exception:
            return None
        return self._coerce_gene_level_types(frame)

    def _write_cached_gene_level(self, project: ProjectConfig, comparison, request: AnalysisRequest, frame: pd.DataFrame) -> None:
        if frame is None or frame.empty:
            return
        cache_dir = self._comparison_cache_dir(project, comparison.comparison_id)
        if cache_dir is None:
            return
        cache_dir.mkdir(parents=True, exist_ok=True)
        table_path = cache_dir / "tx_splicing_gene_table.tsv"
        meta_path = cache_dir / "tx_splicing_cache_meta.json"
        frame.to_csv(table_path, sep="\t", index=False)
        meta_path.write_text(
            json.dumps({"signature": self._comparison_cache_signature(comparison, request)}, indent=2, sort_keys=True),
            encoding="utf-8",
        )

    @staticmethod
    def _comparison_cache_dir(project: ProjectConfig, comparison_id: str) -> Path | None:
        if project.output_root is None:
            return None
        safe_id = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(comparison_id)).strip("_") or "comparison"
        return project.output_root / "04_tx_splicing_integration" / "cache" / safe_id

    def _comparison_cache_signature(self, comparison, request: AnalysisRequest) -> dict[str, object]:
        rmats_files: list[dict[str, object]] = []
        rmats_dir = comparison.rmats_path / "rmats_post" if comparison.rmats_path is not None else None
        if rmats_dir is not None:
            for event_type in EVENT_TYPES:
                path = rmats_dir / f"{event_type}.MATS.{request.rmats_mode}.txt"
                rmats_files.append(self._file_signature(path, label=event_type))
        return {
            "cache_version": 2,
            "comparison_id": comparison.comparison_id,
            "rmats_mode": request.rmats_mode,
            "thresholds": asdict(request.thresholds),
            "reverse_deg": bool(comparison.reverse_deg),
            "reverse_splicing": bool(comparison.reverse_splicing),
            "deg": self._file_signature(comparison.deg_path, label="DEG"),
            "rmats": rmats_files,
        }

    @staticmethod
    def _file_signature(path: Path | None, *, label: str) -> dict[str, object]:
        if path is None or not Path(path).exists():
            return {"label": label, "path": str(path or ""), "exists": False}
        resolved = Path(path)
        stat = resolved.stat()
        return {
            "label": label,
            "path": str(resolved),
            "exists": True,
            "size": stat.st_size,
            "mtime_ns": stat.st_mtime_ns,
        }

    @staticmethod
    def _coerce_gene_level_types(frame: pd.DataFrame) -> pd.DataFrame:
        if frame is None or frame.empty:
            return pd.DataFrame()
        working = frame.copy()
        for column in (
            "log2FC",
            "deg_padj",
            "representative_dPSI",
            "representative_event_FDR",
            "representative_mean_junction_reads_group1",
            "representative_mean_junction_reads_group2",
            "representative_replicate_count_group1",
            "representative_replicate_count_group2",
            "representative_supported_replicate_count_group1",
            "representative_supported_replicate_count_group2",
            "representative_replicate_direction_consistency",
        ):
            if column in working.columns:
                working[column] = pd.to_numeric(working[column], errors="coerce")
        for column in ("deg_sig", "splicing_sig"):
            if column in working.columns:
                working[column] = working[column].map(
                    lambda value: str(value).strip().lower() in {"true", "1", "yes"}
                    if pd.notna(value)
                    else False
                )
        if "n_sig_events_all" in working.columns:
            working["n_sig_events_all"] = pd.to_numeric(working["n_sig_events_all"], errors="coerce").fillna(0).astype(int)
        return working

    def _classify_direction(self, log2fc: float, dpsi: float, deg_sig: bool, as_sig: bool) -> tuple[str, str]:
        if deg_sig and as_sig:
            if log2fc > 0 and dpsi > 0:
                return "both_up_up", "concordant"
            if log2fc > 0 and dpsi < 0:
                return "both_up_down", "discordant"
            if log2fc < 0 and dpsi > 0:
                return "both_down_up", "discordant"
            if log2fc < 0 and dpsi < 0:
                return "both_down_down", "concordant"
            return "both", "single_layer"
        if deg_sig and not as_sig:
            return ("DE-only_up" if log2fc > 0 else "DE-only_down"), "single_layer"
        if as_sig and not deg_sig:
            return ("AS-only_posPSI" if dpsi > 0 else "AS-only_negPSI"), "single_layer"
        return "neither", "neither"

    def _combined_class(self, direction_class: str) -> str:
        if direction_class.startswith("DE-only"):
            return "DE-only"
        if direction_class.startswith("AS-only"):
            return "AS-only"
        if direction_class.startswith("both_") or direction_class == "both":
            return "both"
        return "neither"
