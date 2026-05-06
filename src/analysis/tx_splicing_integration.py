from __future__ import annotations

import math

import pandas as pd

from src.analysis.common import EVENT_TYPES, clean_text, event_uid_from_row, normalize_gene_id, read_table
from src.domain.models.analysis import AnalysisRequest
from src.domain.models.project import ProjectConfig


class TxSplicingIntegrator:
    def run(self, project: ProjectConfig, request: AnalysisRequest) -> tuple[pd.DataFrame, pd.DataFrame]:
        gene_frames: list[pd.DataFrame] = []
        summary_rows: list[dict[str, object]] = []
        selected = set(request.comparison_ids) if request.comparison_ids else None

        for comparison in project.available_comparisons:
            if selected and comparison.comparison_id not in selected:
                continue
            deg_df = self._read_deg(comparison, request)
            rmats_df = self._read_rmats_all(comparison, request)
            gene_level = self._build_gene_level_summary(
                comparison.comparison_id,
                comparison.display_resolved_name,
                deg_df,
                rmats_df,
            )
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
            df = pd.read_csv(fp, sep="\t", low_memory=False)
            if df.empty:
                continue
            df["comparison_id"] = comparison.comparison_id
            df["event_type"] = event_type
            df["gene_id"] = df["GeneID"].map(normalize_gene_id)
            df["geneSymbol"] = df.get("geneSymbol", pd.Series(dtype=str)).map(clean_text)
            df["event_uid"] = df.apply(lambda row: event_uid_from_row(row, event_type), axis=1)
            df["event_FDR"] = pd.to_numeric(df["FDR"], errors="coerce")
            df["dPSI"] = pd.to_numeric(df["IncLevelDifference"], errors="coerce")
            if comparison.reverse_splicing:
                df["dPSI"] = -df["dPSI"]
            df["abs_dPSI"] = df["dPSI"].abs()
            df["splicing_sig"] = (df["event_FDR"] < request.thresholds.splicing_fdr) & (
                df["abs_dPSI"] > request.thresholds.splicing_dpsi
            )
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
        rmats_sig = rmats_df[rmats_df["splicing_sig"]].copy() if not rmats_df.empty else pd.DataFrame()
        deg_genes = set(deg_df["gene_id"].astype(str).tolist()) if not deg_df.empty else set()
        splicing_genes = set(rmats_df["gene_id"].dropna().astype(str).tolist()) if not rmats_df.empty else set()
        all_genes = sorted(deg_genes | splicing_genes)
        deg_map = deg_df.set_index("gene_id")[["log2FC", "deg_padj", "deg_sig"]].to_dict("index") if not deg_df.empty else {}

        rows = []
        for gene_id in all_genes:
            deg_info = deg_map.get(gene_id, {})
            log2fc = deg_info.get("log2FC", math.nan)
            deg_padj = deg_info.get("deg_padj", math.nan)
            deg_sig = bool(deg_info.get("deg_sig", False))

            sub_all = rmats_df[rmats_df["gene_id"].astype(str) == gene_id].copy() if not rmats_df.empty else pd.DataFrame()
            sub_sig = rmats_sig[rmats_sig["gene_id"].astype(str) == gene_id].copy() if not rmats_sig.empty else pd.DataFrame()

            gene_symbol = None
            if not sub_all.empty:
                values = sub_all["geneSymbol"].dropna().astype(str)
                if not values.empty:
                    gene_symbol = values.iloc[0]
            if gene_symbol is None:
                gene_symbol = gene_id

            representative_event_id = None
            representative_event_type = None
            representative_dpsi = math.nan
            representative_event_fdr = math.nan
            splicing_sig = False
            n_sig_events_all = len(sub_sig)
            sig_event_types = ";".join(sorted(sub_sig["event_type"].astype(str).unique().tolist())) if not sub_sig.empty else ""

            if not sub_sig.empty:
                representative = sub_sig.sort_values(["abs_dPSI", "event_FDR"], ascending=[False, True]).iloc[0]
                representative_event_id = representative["event_uid"]
                representative_event_type = representative["event_type"]
                representative_dpsi = representative["dPSI"]
                representative_event_fdr = representative["event_FDR"]
                splicing_sig = True

            direction_class, agreement_class = self._classify_direction(
                log2fc if pd.notna(log2fc) else 0.0,
                representative_dpsi if pd.notna(representative_dpsi) else 0.0,
                deg_sig,
                splicing_sig,
            )

            rows.append(
                {
                    "comparison_id": comparison_id,
                    "comparison_name": comparison_name,
                    "gene_id": gene_id,
                    "geneSymbol": gene_symbol,
                    "log2FC": log2fc,
                    "deg_padj": deg_padj,
                    "deg_sig": deg_sig,
                    "representative_event_id": representative_event_id,
                    "representative_event_type": representative_event_type,
                    "representative_dPSI": representative_dpsi,
                    "representative_event_FDR": representative_event_fdr,
                    "splicing_sig": splicing_sig,
                    "combined_class": self._combined_class(direction_class),
                    "direction_class": direction_class,
                    "agreement_class": agreement_class,
                    "n_sig_events_all": n_sig_events_all,
                    "sig_event_types": sig_event_types,
                }
            )
        return pd.DataFrame(rows)

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
