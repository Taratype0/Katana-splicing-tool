from __future__ import annotations

import json
import re
from datetime import datetime
from pathlib import Path
from typing import Any

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
