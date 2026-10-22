from __future__ import annotations

import json
import re
from datetime import datetime
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


_UNSAFE_PATH = re.compile(r"[^A-Za-z0-9_.-]+")


def _safe_output_name(value: str) -> str:
    token = _UNSAFE_PATH.sub("_", str(value or "").strip()).strip("._")
    return token or "comparison"


class CrossComparisonService:
    """Own the 5.4 output contract independently from project orchestration."""

    def __init__(self, output_root: str | Path) -> None:
        self.output_root = Path(output_root)

    def base_dir(self) -> Path:
        return self.output_root / "03_cross_comparison_candidate_comparison"

    @staticmethod
    def _safe_bool(value: object) -> bool:
        if pd.isna(value):
            return False
        if isinstance(value, str):
            return value.strip().casefold() in {"1", "true", "yes", "y"}
        return bool(value)

    def candidate_files(self) -> dict[str, Path]:
        base = self.base_dir()
        return {
            "matrix": base / "cross_comparison_candidate_matrix.tsv",
            "summary": base / "cross_comparison_summary.tsv",
        }

    def pattern_dir(self, pair_id: str | None) -> Path | None:
        return self.base_dir() / _safe_output_name(pair_id) if pair_id else None

    def pattern_files(self, pair_id: str | None) -> dict[str, Path]:
        base = self.pattern_dir(pair_id)
        if base is None:
            return {}
        return {
            "same_direction_large_delta": base / "same_direction_large_delta_significant_as_events.tsv",
            "opposite_direction": base / "opposite_direction_significant_as_events.tsv",
            "all_shared": base / "all_shared_significant_as_pairwise_dpsi_delta.tsv",
            "condition_specific": base / "optional_condition_specific_as_events.tsv",
            "summary": base / "summary_counts.tsv",
        }

    def interaction_dir(self, pair_id: str | None) -> Path | None:
        pattern = self.pattern_dir(pair_id)
        return pattern / "as_interaction_modeling" if pattern is not None else None

    def interaction_files(self, pair_id: str | None) -> dict[str, Path]:
        base = self.interaction_dir(pair_id)
        if base is None:
            return {}
        return {
            "unified_event_table": base / "unified_as_interaction_events.tsv",
            "shared_interaction": base / "shared_significant_as_interaction_model.tsv",
            "opposite_direction": base / "opposite_direction_interaction_events.tsv",
            "large_delta": base / "large_delta_interaction_events.tsv",
            "exclusive": base / "exclusive_significant_as_support.tsv",
            "benchmark": base / "delta_dpsi_vs_interaction_benchmark.tsv",
            "diagnostics": base / "interaction_diagnostics.tsv",
        }

    def interaction_batch_dir(self) -> Path:
        return self.base_dir() / "as_interaction_modeling_all_pairs"

    def interaction_batch_files(self) -> dict[str, Path]:
        base = self.interaction_batch_dir()
        return {
            "batch_summary": base / "all_pair_interaction_summary.tsv",
            "batch_top_interactions": base / "all_pair_top_interactions.tsv",
            "batch_exclusive": base / "all_pair_exclusive_support.tsv",
            "batch_benchmark": base / "all_pair_delta_dpsi_vs_interaction_benchmark.tsv",
            "batch_diagnostics": base / "all_pair_interaction_diagnostics.tsv",
            "parameters": base / "parameters.json",
        }

    @staticmethod
    def _write_json(path: Path, payload: dict[str, object]) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        temporary = path.with_suffix(path.suffix + ".tmp")
        temporary.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")
        temporary.replace(path)

    def write_pattern_outputs(
        self,
        *,
        pair: Any,
        tables: dict[str, pd.DataFrame],
        abs_cutoff: float,
        delta_cutoff: float,
    ) -> list[Path]:
        base = self.pattern_dir(pair.pair_id)
        if base is None:
            return []
        base.mkdir(parents=True, exist_ok=True)
        file_map = self.pattern_files(pair.pair_id)
        generated: list[Path] = []
        for key in ("same_direction_large_delta", "opposite_direction", "all_shared", "condition_specific"):
            path = file_map[key]
            tables.get(key, pd.DataFrame()).to_csv(path, sep="\t", index=False)
            generated.append(path)
        summary = pd.DataFrame(
            [
                {
                    "pair_id": pair.pair_id,
                    "pair_name": pair.resolved_name,
                    "comparison_A": pair.comparison_a,
                    "comparison_B": pair.comparison_b,
                    "n_same_direction_large_delta": len(tables.get("same_direction_large_delta", pd.DataFrame())),
                    "n_opposite_direction": len(tables.get("opposite_direction", pd.DataFrame())),
                    "n_all_shared": len(tables.get("all_shared", pd.DataFrame())),
                    "n_condition_specific": len(tables.get("condition_specific", pd.DataFrame())),
                    "abs_dPSI_significance_cutoff_used": float(abs_cutoff),
                    "large_delta_dPSI_cutoff_used": float(delta_cutoff),
                }
            ]
        )
        summary.to_csv(file_map["summary"], sep="\t", index=False)
        generated.append(file_map["summary"])
        return generated

    def write_interaction_outputs(
        self,
        *,
        pair: Any,
        tables: dict[str, pd.DataFrame],
        abs_cutoff: float,
        delta_cutoff: float,
        min_total_reads: int,
    ) -> list[Path]:
        base = self.interaction_dir(pair.pair_id)
        if base is None:
            return []
        base.mkdir(parents=True, exist_ok=True)
        file_map = self.interaction_files(pair.pair_id)
        generated: list[Path] = []
        for key in (
            "unified_event_table",
            "shared_interaction",
            "opposite_direction",
            "large_delta",
            "exclusive",
            "benchmark",
            "diagnostics",
        ):
            path = file_map[key]
            tables.get(key, pd.DataFrame()).to_csv(path, sep="\t", index=False)
            generated.append(path)
        parameters_path = base / "parameters.json"
        self._write_json(
            parameters_path,
            {
                "module": "5.4.3 Cross-comparison AS interaction modeling",
                "description": (
                    "Count-based binomial GLM interaction model for shared significant AS events using replicate-level "
                    "rMATS IJC/SJC counts. The interaction term tests whether the group effect differs between comparison "
                    "A and comparison B; Pearson overdispersion is used as a quasi-binomial-style SE correction. "
                    "Empirical-Bayes shrinkage stabilizes interaction log-odds for prioritization."
                ),
                "pair_id": pair.pair_id,
                "pair_name": pair.resolved_name,
                "comparison_A": pair.comparison_a,
                "comparison_B": pair.comparison_b,
                "abs_dPSI_significance_cutoff": float(abs_cutoff),
                "large_delta_dPSI_cutoff": float(delta_cutoff),
                "min_total_reads": int(min_total_reads),
            },
        )
        generated.append(parameters_path)
        return generated

    def write_interaction_batch_outputs(
        self,
        *,
        tables: dict[str, pd.DataFrame],
        abs_cutoff: float,
        delta_cutoff: float,
        min_total_reads: int,
    ) -> list[Path]:
        base = self.interaction_batch_dir()
        base.mkdir(parents=True, exist_ok=True)
        file_map = self.interaction_batch_files()
        generated: list[Path] = []
        for key in (
            "batch_summary",
            "batch_top_interactions",
            "batch_exclusive",
            "batch_benchmark",
            "batch_diagnostics",
        ):
            path = file_map[key]
            tables.get(key, pd.DataFrame()).to_csv(path, sep="\t", index=False)
            generated.append(path)
        self._write_json(
            file_map["parameters"],
            {
                "module": "5.4.3 Cross-comparison AS interaction modeling all-pair batch summary",
                "description": (
                    "Runs the replicate-aware AS interaction model independently for every selected comparison pair and "
                    "aggregates pair-level summaries, top events and A/B-specific read-support classifications."
                ),
                "abs_dPSI_significance_cutoff": float(abs_cutoff),
                "large_delta_dPSI_cutoff": float(delta_cutoff),
                "min_total_reads": int(min_total_reads),
                "generated_at": datetime.now().isoformat(timespec="seconds"),
            },
        )
        generated.append(file_map["parameters"])
        return generated

    @staticmethod
    def empty_pattern_tables() -> dict[str, pd.DataFrame]:
        columns = [
            "gene",
            "event_type",
            "event_id",
            "event_matching_key",
            "comparison_A",
            "comparison_B",
            "experiment_A_display_name",
            "experiment_B_display_name",
            "dPSI_A",
            "dPSI_B",
            "FDR_A",
            "FDR_B",
            "significant_A",
            "significant_B",
            "direction_A",
            "direction_B",
            "delta_dPSI",
            "abs_delta_dPSI",
            "same_direction",
            "opposite_direction",
            "same_direction_strength_delta",
            "event_class",
            "abs_dPSI_significance_cutoff_used",
            "large_delta_dPSI_cutoff_used",
        ]
        empty = pd.DataFrame(columns=columns)
        return {
            "same_direction_large_delta": empty.copy(),
            "opposite_direction": empty.copy(),
            "all_shared": empty.copy(),
            "condition_specific": empty.copy(),
        }

    def build_pattern_tables(
        self,
        *,
        pair: Any,
        comparison_a: pd.DataFrame,
        comparison_b: pd.DataFrame,
        abs_dpsi_cutoff: float,
        large_delta_dpsi_cutoff: float,
    ) -> dict[str, pd.DataFrame]:
        if comparison_a.empty and comparison_b.empty:
            return self.empty_pattern_tables()

        merged = comparison_a.merge(
            comparison_b,
            on="event_matching_key",
            how="outer",
            suffixes=("_A", "_B"),
        )
        if merged.empty:
            return self.empty_pattern_tables()

        def combined_text(column: str) -> pd.Series:
            return (
                merged.get(f"{column}_A", pd.Series(pd.NA, index=merged.index))
                .combine_first(merged.get(f"{column}_B", pd.Series(pd.NA, index=merged.index)))
                .fillna("")
                .astype(str)
            )

        merged["gene"] = combined_text("gene")
        merged["event_type"] = combined_text("event_type")
        merged["event_id"] = combined_text("event_id")
        merged["comparison_A"] = pair.comparison_a or ""
        merged["comparison_B"] = pair.comparison_b or ""
        merged["experiment_A_display_name"] = (
            comparison_a["comparison_display_name"].iloc[0]
            if not comparison_a.empty and "comparison_display_name" in comparison_a.columns
            else str(pair.comparison_a or "")
        )
        merged["experiment_B_display_name"] = (
            comparison_b["comparison_display_name"].iloc[0]
            if not comparison_b.empty and "comparison_display_name" in comparison_b.columns
            else str(pair.comparison_b or "")
        )
        for column in ("dPSI_A", "dPSI_B", "FDR_A", "FDR_B"):
            merged[column] = pd.to_numeric(merged.get(column), errors="coerce")
        for column in ("significant_A", "significant_B"):
            merged[column] = (
                merged.get(column, pd.Series(False, index=merged.index, dtype=bool))
                .map(self._safe_bool)
            )
        for column in ("direction_A", "direction_B"):
            merged[column] = merged.get(column, pd.Series("", index=merged.index)).fillna("").astype(str)

        merged["delta_dPSI"] = merged["dPSI_A"] - merged["dPSI_B"]
        merged["abs_delta_dPSI"] = merged["delta_dPSI"].abs()
        sign_a = np.sign(merged["dPSI_A"].fillna(0))
        sign_b = np.sign(merged["dPSI_B"].fillna(0))
        merged["same_direction"] = (sign_a == sign_b) & sign_a.ne(0) & sign_b.ne(0)
        merged["opposite_direction"] = sign_a.ne(sign_b) & sign_a.ne(0) & sign_b.ne(0)
        merged["same_direction_strength_delta"] = (
            merged["dPSI_A"].abs() - merged["dPSI_B"].abs()
        ).abs()
        merged["abs_dPSI_significance_cutoff_used"] = float(abs_dpsi_cutoff)
        merged["large_delta_dPSI_cutoff_used"] = float(large_delta_dpsi_cutoff)

        shared = merged.loc[merged["significant_A"] & merged["significant_B"]].copy()
        shared["event_class"] = "shared_significant"
        same_direction_large_delta = shared.loc[
            shared["same_direction"]
            & shared["same_direction_strength_delta"].ge(large_delta_dpsi_cutoff)
        ].copy()
        same_direction_large_delta["event_class"] = "same_direction_large_delta"
        opposite_direction = shared.loc[shared["opposite_direction"]].copy()
        opposite_direction["event_class"] = "opposite_direction"

        condition_specific = merged.loc[
            (merged["significant_A"] & ~merged["significant_B"])
            | (~merged["significant_A"] & merged["significant_B"])
        ].copy()
        condition_specific["event_class"] = np.where(
            condition_specific["significant_A"],
            "comparison_A_only",
            "comparison_B_only",
        )

        output_columns = self.empty_pattern_tables()["all_shared"].columns
        return {
            "same_direction_large_delta": same_direction_large_delta.reindex(
                columns=output_columns
            ).reset_index(drop=True),
            "opposite_direction": opposite_direction.reindex(columns=output_columns).reset_index(
                drop=True
            ),
            "all_shared": shared.reindex(columns=output_columns).reset_index(drop=True),
            "condition_specific": condition_specific.reindex(columns=output_columns).reset_index(
                drop=True
            ),
        }
