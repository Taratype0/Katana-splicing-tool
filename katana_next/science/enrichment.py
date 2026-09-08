from __future__ import annotations

import math
from pathlib import Path
from typing import Iterable

import pandas as pd

from katana_next.science.statistics import bh_adjust, hypergeometric_overlap


def read_gene_sets(path: Path) -> dict[str, set[str]]:
    sets: dict[str, set[str]] = {}
    if path.suffix.lower() == ".gmt":
        for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
            parts = [part.strip() for part in line.split("\t")]
            if len(parts) >= 3:
                sets[parts[0]] = {gene for gene in parts[2:] if gene}
        return sets
    frame = pd.read_csv(path, sep=None, engine="python", low_memory=False)
    set_column = next((column for column in ("set", "term", "pathway", "rbp") if column in frame.columns), None)
    gene_column = next((column for column in ("gene", "gene_symbol", "symbol", "target") if column in frame.columns), None)
    if set_column is None or gene_column is None:
        raise ValueError(f"Gene-set file requires set/term and gene/symbol columns: {path}")
    for term, group in frame.groupby(set_column):
        sets[str(term)] = set(group[gene_column].dropna().astype(str).str.strip())
    return sets


def overrepresentation(
    selected: Iterable[str],
    background: Iterable[str],
    gene_sets: dict[str, set[str]],
    *,
    minimum_overlap: int = 2,
) -> pd.DataFrame:
    selected_set = {str(value) for value in selected if str(value)}
    background_set = {str(value) for value in background if str(value)}
    selected_set &= background_set
    rows = []
    for term, members in gene_sets.items():
        tested_members = members & background_set
        overlap = selected_set & tested_members
        if len(overlap) < minimum_overlap:
            continue
        pvalue = hypergeometric_overlap(
            len(background_set), len(selected_set), len(tested_members), len(overlap)
        )
        expected = len(selected_set) * len(tested_members) / max(1, len(background_set))
        rows.append(
            {
                "row_id": term,
                "term": term,
                "selected_genes": len(selected_set),
                "background_genes": len(background_set),
                "term_genes_in_background": len(tested_members),
                "overlap_genes": len(overlap),
                "expected_overlap": expected,
                "fold_enrichment": len(overlap) / expected if expected else math.nan,
                "pvalue": pvalue,
                "genes": ", ".join(sorted(overlap)),
            }
        )
    result = pd.DataFrame(rows)
    if not result.empty:
        result["fdr"] = bh_adjust(result["pvalue"])
        result = result.sort_values(["fdr", "fold_enrichment"], ascending=[True, False])
    return result
