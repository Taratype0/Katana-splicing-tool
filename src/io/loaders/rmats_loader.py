from __future__ import annotations

from pathlib import Path

import pandas as pd

from src.analysis.common import EVENT_TYPES, clean_text, event_uid_from_row, read_table
from src.domain.models.project import ComparisonDefinition


class RMATSLoader:
    REQUIRED_COLUMNS = {"GeneID", "geneSymbol", "FDR", "IncLevelDifference"}

    def load_event_table(
        self,
        comparison: ComparisonDefinition,
        mode: str,
        event_type: str,
    ) -> pd.DataFrame:
        if comparison.rmats_path is None:
            return pd.DataFrame()
        path = comparison.rmats_path / "rmats_post" / f"{event_type}.MATS.{mode}.txt"
        df = read_table(path)
        if df is None or df.empty:
            return pd.DataFrame()
        df = df.copy()
        df["comparison_id"] = comparison.comparison_id
        df["comparison_name"] = comparison.display_resolved_name
        df["event_type"] = event_type
        df["gene_id"] = df["GeneID"].map(clean_text)
        df["gene_symbol"] = df["geneSymbol"].map(clean_text)
        df["fdr"] = pd.to_numeric(df["FDR"], errors="coerce")
        df["dpsi"] = pd.to_numeric(df["IncLevelDifference"], errors="coerce")
        if comparison.reverse_splicing:
            df["dpsi"] = -df["dpsi"]
        df["abs_dpsi"] = df["dpsi"].abs()
        df["event_uid"] = df.apply(lambda row: event_uid_from_row(row, event_type), axis=1)
        return df

    def load_all_events(self, comparison: ComparisonDefinition, mode: str) -> pd.DataFrame:
        frames = [self.load_event_table(comparison, mode, event_type) for event_type in EVENT_TYPES]
        frames = [frame for frame in frames if not frame.empty]
        return pd.concat(frames, axis=0, ignore_index=True) if frames else pd.DataFrame()

    def validate(self, path: Path) -> list[str]:
        df = read_table(path)
        if df is None:
            return [f"Missing rMATS file: {path}"]
        missing = sorted(self.REQUIRED_COLUMNS - set(df.columns))
        if missing:
            return [f"rMATS file missing columns {missing}: {path.name}"]
        return []
