from __future__ import annotations

from pathlib import Path

from src.analysis.common import read_table


class CountsLoader:
    REQUIRED_COLUMNS = {"gene_id"}

    def validate(self, path: Path) -> list[str]:
        df = read_table(path)
        if df is None:
            return [f"Missing normalized counts file: {path}"]
        missing = sorted(self.REQUIRED_COLUMNS - set(df.columns))
        if missing:
            return [f"Normalized counts file missing columns {missing}: {path.name}"]
        return []

