from __future__ import annotations

from pathlib import Path

import pandas as pd

from src.analysis.common import clean_text, read_table


class SUPPALoader:
    def load(self, path: Path) -> pd.DataFrame:
        df = read_table(path)
        if df is None or df.empty:
            return pd.DataFrame()
        df = df.copy()
        if df.shape[1] == 2:
            df = df.reset_index().rename(columns={"index": "event_id", df.columns[0]: "dpsi", df.columns[1]: "pvalue"})
        else:
            df = df.rename(columns={df.columns[0]: "event_id", df.columns[1]: "dpsi", df.columns[2]: "pvalue"})
        df["event_id"] = df["event_id"].astype(str)
        df["dpsi"] = pd.to_numeric(df["dpsi"], errors="coerce")
        df["pvalue"] = pd.to_numeric(df["pvalue"], errors="coerce")
        df["gene_ids"] = df["event_id"].map(self._split_gene_field)
        return df

    def validate(self, path: Path) -> list[str]:
        df = read_table(path)
        if df is None:
            return [f"Missing SUPPA file: {path}"]
        if df.shape[1] < 2:
            return [f"SUPPA file has too few columns: {path.name}"]
        return []

    def _split_gene_field(self, event_id: str) -> list[str]:
        gene_part = str(event_id).split(";")[0]
        return [cleaned for raw in gene_part.split("_and_") if (cleaned := clean_text(raw))]

