from __future__ import annotations

import math
from pathlib import Path

import pandas as pd

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
                        }
                    )
                    continue
                df = pd.read_csv(fp, sep="\t", low_memory=False)
                if df.empty:
                    continue
                fdr = pd.to_numeric(df["FDR"], errors="coerce")
                dpsi = pd.to_numeric(df["IncLevelDifference"], errors="coerce").abs()
                sig = (fdr < request.thresholds.splicing_fdr) & (dpsi >= request.thresholds.splicing_dpsi)
                rows.append(
                    {
                        "comparison_id": comparison.comparison_id,
                        "comparison_name": comparison.resolved_name,
                        "event_type": event_type,
                        "n_total": int(len(df)),
                        "n_sig": int(sig.fillna(False).sum()),
                    }
                )
        return pd.DataFrame(rows)
