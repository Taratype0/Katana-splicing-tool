from __future__ import annotations

import pandas as pd

from src.analysis.common import EVENT_TYPES, event_uid_from_row, normalize_gene_id
from src.domain.models.analysis import AnalysisRequest
from src.domain.models.project import ComparisonDefinition, ProjectConfig


class SplicingProgramComparator:
    def run(
        self,
        project: ProjectConfig,
        request: AnalysisRequest,
        comparison_a: str,
        comparison_b: str,
    ) -> tuple[pd.DataFrame, pd.DataFrame]:
        left = self._read_comparison(project, request, comparison_a).rename(
            columns={
                "gene_id": "gene_id_A",
                "gene_symbol": "gene_symbol_A",
                "event_type": "event_type_A",
                "fdr": "FDR_A",
                "dpsi": "dPSI_A",
                "abs_dpsi": "abs_dPSI_A",
                "is_sig": "is_sig_A",
            }
        )
        right = self._read_comparison(project, request, comparison_b).rename(
            columns={
                "gene_id": "gene_id_B",
                "gene_symbol": "gene_symbol_B",
                "event_type": "event_type_B",
                "fdr": "FDR_B",
                "dpsi": "dPSI_B",
                "abs_dpsi": "abs_dPSI_B",
                "is_sig": "is_sig_B",
            }
        )

        merged = left.merge(right, on="event_uid", how="outer")
        if merged.empty:
            return merged, pd.DataFrame(columns=["class", "n_events"])
        merged["gene_symbol"] = merged["gene_symbol_A"].combine_first(merged["gene_symbol_B"])
        merged["gene_id"] = merged["gene_id_A"].combine_first(merged["gene_id_B"])
        merged["event_type"] = merged["event_type_A"].combine_first(merged["event_type_B"])
        merged["is_sig_A"] = merged["is_sig_A"].fillna(False)
        merged["is_sig_B"] = merged["is_sig_B"].fillna(False)
        merged["class"] = merged.apply(
            lambda row: self._classify(
                row["is_sig_A"],
                row["is_sig_B"],
                row["dPSI_A"],
                row["dPSI_B"],
                request.thresholds.program_delta_dpsi,
            ),
            axis=1,
        )
        display_a = self._display_name(project, comparison_a)
        display_b = self._display_name(project, comparison_b)
        merged["class_label"] = merged["class"].map(
            {
                "A_only": f"Only in {display_a}",
                "B_only": f"Only in {display_b}",
                "shared_same_direction": "Shared same direction",
                "same_direction_large_delta": "Same direction, large dPSI gap",
                "opposite_direction": "Opposite direction",
                "shared_zero_direction": "Shared zero direction",
                "shared_unresolved": "Shared unresolved",
                "neither": "Not significant in either",
            }
        )
        merged["abs_delta_between_programs"] = (
            pd.to_numeric(merged["dPSI_A"], errors="coerce")
            - pd.to_numeric(merged["dPSI_B"], errors="coerce")
        ).abs()
        merged["comparison_A"] = comparison_a
        merged["comparison_B"] = comparison_b
        merged["comparison_A_name"] = display_a
        merged["comparison_B_name"] = display_b

        summary = merged["class_label"].value_counts(dropna=False).reset_index()
        summary.columns = ["class_label", "n_events"]
        summary["comparison_A"] = comparison_a
        summary["comparison_B"] = comparison_b
        summary["comparison_A_name"] = display_a
        summary["comparison_B_name"] = display_b
        return merged, summary

    def _display_name(self, project: ProjectConfig, comparison_id: str) -> str:
        comparison = next(
            (item for item in project.available_comparisons if item.comparison_id == comparison_id),
            None,
        )
        return comparison.resolved_name if comparison is not None else comparison_id

    def _read_comparison(
        self,
        project: ProjectConfig,
        request: AnalysisRequest,
        comparison_id: str,
    ) -> pd.DataFrame:
        comparison = next(
            (item for item in project.available_comparisons if item.comparison_id == comparison_id),
            None,
        )
        if comparison is None or comparison.rmats_path is None:
            return pd.DataFrame(columns=["event_uid"])

        rows: list[pd.DataFrame] = []
        rmats_post = comparison.rmats_path / "rmats_post"
        for event_type in EVENT_TYPES:
            fp = rmats_post / f"{event_type}.MATS.{request.rmats_mode}.txt"
            if not fp.exists():
                continue
            df = pd.read_csv(fp, sep="\t", low_memory=False)
            if df.empty:
                continue
            df["event_uid"] = df.apply(lambda row: event_uid_from_row(row, event_type), axis=1)
            df["event_type"] = event_type
            df["gene_id"] = df.get("GeneID", pd.Series(dtype=str)).map(normalize_gene_id)
            df["gene_symbol"] = df.get("geneSymbol", pd.Series(dtype=str)).astype(str)
            df["fdr"] = pd.to_numeric(df["FDR"], errors="coerce")
            df["dpsi"] = pd.to_numeric(df["IncLevelDifference"], errors="coerce")
            if comparison.reverse_splicing:
                df["dpsi"] = -df["dpsi"]
            df["abs_dpsi"] = df["dpsi"].abs()
            df["is_sig"] = (df["fdr"] < request.thresholds.splicing_fdr) & (
                df["abs_dpsi"] >= request.thresholds.splicing_dpsi
            )
            rows.append(
                df[
                    ["event_uid", "event_type", "gene_id", "gene_symbol", "fdr", "dpsi", "abs_dpsi", "is_sig"]
                ].copy()
            )

        if not rows:
            return pd.DataFrame(columns=["event_uid"])
        output = pd.concat(rows, axis=0, ignore_index=True)
        return output.sort_values(
            ["event_uid", "is_sig", "abs_dpsi"], ascending=[True, False, False]
        ).drop_duplicates("event_uid")

    def _classify(
        self,
        sig_a: bool,
        sig_b: bool,
        dpsi_a: float,
        dpsi_b: float,
        program_delta_cutoff: float,
    ) -> str:
        if sig_a and sig_b:
            if pd.notna(dpsi_a) and pd.notna(dpsi_b):
                if dpsi_a == 0 or dpsi_b == 0:
                    return "shared_zero_direction"
                if dpsi_a * dpsi_b > 0:
                    if abs(float(dpsi_a) - float(dpsi_b)) >= float(program_delta_cutoff):
                        return "same_direction_large_delta"
                    return "shared_same_direction"
                return "opposite_direction"
            return "shared_unresolved"
        if sig_a and not sig_b:
            return "A_only"
        if sig_b and not sig_a:
            return "B_only"
        return "neither"
