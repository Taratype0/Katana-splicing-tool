from __future__ import annotations

from functools import lru_cache
from pathlib import Path
import re

import pandas as pd

EVENT_TYPES = ["SE", "RI", "A3SS", "A5SS", "MXE"]


def event_uid_from_row(row: pd.Series, event_type: str) -> str:
    cols = [
        "chr",
        "strand",
        "exonStart_0base",
        "exonEnd",
        "upstreamES",
        "upstreamEE",
        "downstreamES",
        "downstreamEE",
        "riExonStart_0base",
        "riExonEnd",
        "longExonStart_0base",
        "longExonEnd",
        "shortES",
        "shortEE",
        "flankingES",
        "flankingEE",
        "1stExonStart_0base",
        "1stExonEnd",
        "2ndExonStart_0base",
        "2ndExonEnd",
    ]
    vals = [event_type]
    for col in cols:
        if col in row.index:
            vals.append(str(row[col]))
    return "|".join(vals)


def clean_text(value: object) -> str | None:
    if pd.isna(value):
        return None
    text = str(value).strip().strip('"').strip("'")
    if not text or text.lower() == "nan":
        return None
    return text


def normalize_gene_id(value: object) -> str | None:
    text = clean_text(value)
    if text is None:
        return None
    match = re.match(r"^(?P<core>ENS[A-Z]*G[0-9]+)(?:\.\d+)?$", text, flags=re.IGNORECASE)
    if match:
        return match.group("core")
    if "." in text:
        left, right = text.rsplit(".", 1)
        if right.isdigit() and left:
            return left
    return text


def read_table(path: Path | None) -> pd.DataFrame | None:
    if path is None or not path.exists():
        return None
    if path.suffix.lower() == ".csv":
        try:
            return pd.read_csv(path, sep=",", low_memory=False)
        except Exception:
            return pd.read_csv(path, sep="\t", low_memory=False)
    return pd.read_csv(path, sep="\t", low_memory=False)


@lru_cache(maxsize=8)
def _annotation_gene_symbol_lookup_cached(annotation_path: str) -> dict[str, str]:
    path = Path(annotation_path)
    df = read_table(path)
    if df is None or df.empty:
        return {}

    gene_id_column = next(
        (column for column in ("gene_id", "GeneID", "groupID") if column in df.columns),
        None,
    )
    gene_name_column = next(
        (column for column in ("gene_name", "geneSymbol", "gene_symbol", "symbol", "gene") if column in df.columns),
        None,
    )
    if gene_id_column is None or gene_name_column is None:
        return {}

    working = df[[gene_id_column, gene_name_column]].copy()
    working["normalized_gene_id"] = working[gene_id_column].map(normalize_gene_id)
    working["resolved_gene_symbol"] = working[gene_name_column].map(clean_text)
    working = working.dropna(subset=["normalized_gene_id", "resolved_gene_symbol"]).copy()
    if working.empty:
        return {}
    working = working.drop_duplicates(subset=["normalized_gene_id"], keep="first")
    return dict(
        zip(
            working["normalized_gene_id"].astype(str),
            working["resolved_gene_symbol"].astype(str),
            strict=False,
        )
    )


def load_annotation_gene_symbol_lookup(annotation_path: Path | str | None) -> dict[str, str]:
    if annotation_path is None:
        return {}
    path = Path(annotation_path)
    if not path.exists():
        return {}
    return _annotation_gene_symbol_lookup_cached(str(path.resolve()))
