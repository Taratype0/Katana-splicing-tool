from __future__ import annotations

from pathlib import Path

from src.analysis.common import read_table


class QuantLoader:
    REQUIRED_COLUMNS = {"Name", "Length", "EffectiveLength", "TPM", "NumReads"}

    def validate(self, path: Path) -> list[str]:
        df = read_table(path)
        if df is None:
            return [f"Missing quant.sf file: {path}"]
        missing = sorted(self.REQUIRED_COLUMNS - set(df.columns))
        if missing:
            return [f"quant.sf missing columns {missing}: {path.name}"]
        return []

