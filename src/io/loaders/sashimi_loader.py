from __future__ import annotations

from pathlib import Path

from src.analysis.common import read_table


class SashimiLoader:
    REQUIRED_COLUMNS = {
        "geneSymbol",
        "GeneID",
        "display_group",
        "class_interpretation",
        "event_type",
        "event_uid",
        "comparison_id",
        "label1",
        "label2",
        "b1_txt",
        "b2_txt",
        "event_file",
        "outdir",
    }

    def validate_manifest(self, path: Path) -> list[str]:
        df = read_table(path)
        if df is None:
            return [f"Missing rmats2sashimi manifest: {path}"]
        missing = sorted(self.REQUIRED_COLUMNS - set(df.columns))
        if missing:
            return [f"rmats2sashimi manifest missing columns {missing}: {path.name}"]
        return []

