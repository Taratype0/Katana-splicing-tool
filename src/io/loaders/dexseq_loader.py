from __future__ import annotations

from pathlib import Path

from src.analysis.common import read_table


class DEXSeqLoader:
    REQUIRED_COLUMNS = {"groupID", "padj"}

    def validate(self, path: Path) -> list[str]:
        df = read_table(path)
        if df is None:
            return [f"Missing DEXSeq file: {path}"]
        missing = sorted(self.REQUIRED_COLUMNS - set(df.columns))
        if missing:
            return [f"DEXSeq file missing columns {missing}: {path.name}"]
        return []

