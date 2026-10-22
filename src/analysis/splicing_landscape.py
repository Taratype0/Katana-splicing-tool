from __future__ import annotations

import math

import pandas as pd

from src.analysis.stable_as import apply_stable_as_filter
from src.domain.models.analysis import AnalysisRequest
from src.domain.models.project import ProjectConfig

EVENT_TYPES = ["SE", "RI", "A3SS", "A5SS", "MXE"]


class SplicingLandscapeAnalyzer:

    def run(self, project: ProjectConfig, request: AnalysisRequest) -> pd.DataFrame:
        rows: list[dict[str, object]] = []
        mode = request.rmats_mode
        for comparison in project.available_comparisons:
            if request.comparison_ids and comparison.comparison_id not in request.comparison_ids:
                continue
            if comparison.rmats_path is None:
                continue
            rmats_post = comparison.rmats_path / "rmats_post"
            for event_type in EVENT_TYPES:
                fp = rmats_post / f"{event_type}.MATS.{mode}.txt"
                if not fp.exists():
                    rows.append(
                        {
                            "comparison_id": comparison.comparison_id,
                            "comparison_name": comparison.resolved_name,
                            "event_type": event_type,
                            "n_total": math.nan,
                            "n_sig": math.nan,
                            "n_stable": math.nan,
                            "stable_fraction_of_sig": math.nan,
                            "stable_min_junction_reads": request.thresholds.stable_min_junction_reads,
                            "stable_min_replicates_per_group": request.thresholds.stable_min_replicates_per_group,
                            "stable_min_pairwise_consistency": request.thresholds.stable_min_pairwise_consistency,
                        }
                    )
                    continue
                df = pd.read_csv(fp, sep="\t", low_memory=False)
                if df.empty:
                    continue
                fdr = pd.to_numeric(df["FDR"], errors="coerce")
                standardized_dpsi = pd.to_numeric(df["IncLevelDifference"], errors="coerce")
                if comparison.reverse_splicing:
                    standardized_dpsi = -standardized_dpsi
                evaluated = apply_stable_as_filter(
                    df,
                    fdr=fdr,
                    standardized_dpsi=standardized_dpsi,
                    reverse_splicing=comparison.reverse_splicing,
                    thresholds=request.thresholds,
                )
                n_sig = int(evaluated["rmats_sig_before_stability"].sum())
                n_stable = int(evaluated["stable_as_event"].sum())
                rows.append(
                    {
                        "comparison_id": comparison.comparison_id,
                        "comparison_name": comparison.resolved_name,
                        "event_type": event_type,
                        "n_total": int(len(df)),
                        "n_sig": n_sig,
                        "n_stable": n_stable,
                        "stable_fraction_of_sig": float(n_stable / n_sig) if n_sig else 0.0,
                        "stable_min_junction_reads": request.thresholds.stable_min_junction_reads,
                        "stable_min_replicates_per_group": request.thresholds.stable_min_replicates_per_group,
                        "stable_min_pairwise_consistency": request.thresholds.stable_min_pairwise_consistency,
                    }
                )
        return pd.DataFrame(rows)
