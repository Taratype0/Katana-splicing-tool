from __future__ import annotations

from pathlib import Path

from src.analysis.common import read_table


class DTULoader:
    RESULT_COLUMNS = {"groupID", "featureID", "padj"}
    STAGER_COLUMNS = {"geneID", "txID", "gene", "transcript"}

    def validate_result(self, path: Path) -> list[str]:
        df = read_table(path)
        if df is None:
            return [f"Missing DTU result file: {path}"]
        missing = sorted(self.RESULT_COLUMNS - set(df.columns))
        if missing:
            return [f"DTU result file missing columns {missing}: {path.name}"]
        return []

    def validate_stager(self, path: Path) -> list[str]:
        df = read_table(path)
        if df is None:
            return [f"Missing stageR file: {path}"]
        missing = sorted(self.STAGER_COLUMNS - set(df.columns))
        if missing:
            return [f"stageR file missing columns {missing}: {path.name}"]
        return []

