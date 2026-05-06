from __future__ import annotations

from pathlib import Path

from src.analysis.common import read_table


class CardsLoader:
    SHORTLIST_COLUMNS = {"geneSymbol", "GeneID", "event_type", "comparison_A", "comparison_B"}
    EXPR_COLUMNS = {"gene_id", "group", "expr"}

    def validate_shortlist(self, path: Path) -> list[str]:
        df = read_table(path)
        if df is None:
            return [f"Missing shortlist file: {path}"]
        missing = sorted(self.SHORTLIST_COLUMNS - set(df.columns))
        if missing:
            return [f"Shortlist file missing columns {missing}: {path.name}"]
        return []

    def validate_expression_support(self, path: Path) -> list[str]:
        df = read_table(path)
        if df is None:
            return [f"Missing expression support file: {path}"]
        missing = sorted(self.EXPR_COLUMNS - set(df.columns))
        if missing:
            return [f"Expression support file missing columns {missing}: {path.name}"]
        return []

