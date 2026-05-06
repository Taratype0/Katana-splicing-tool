from __future__ import annotations

from pathlib import Path

import pandas as pd

from src.analysis.common import read_table
from src.domain.models.project import ComparisonDefinition


class DEGLoader:
    REQUIRED_COLUMNS = {"gene_id", "log2FoldChange", "padj"}

    def load(self, comparison: ComparisonDefinition) -> pd.DataFrame:
        df = read_table(comparison.deg_path)
        if df is None or df.empty:
            return pd.DataFrame()
        df = df.copy()
        df["comparison_id"] = comparison.comparison_id
        df["comparison_name"] = comparison.display_resolved_name
        df["log2FC"] = pd.to_numeric(df["log2FoldChange"], errors="coerce")
        if comparison.reverse_deg:
            df["log2FC"] = -df["log2FC"]
        df["deg_padj"] = pd.to_numeric(df["padj"], errors="coerce")
        return df

    def validate(self, path: Path) -> list[str]:
        df = read_table(path)
        if df is None:
            return [f"Missing DEG file: {path}"]
        missing = sorted(self.REQUIRED_COLUMNS - set(df.columns))
        if missing:
            return [f"DEG file missing columns {missing}: {path.name}"]
        return []
