from __future__ import annotations

import math
import re
from pathlib import Path
from typing import Iterable

import pandas as pd

from katana_next.science.statistics import bh_adjust


def gene_token(value: object) -> str:
    if pd.isna(value):
        return ""
    text = str(value).strip().strip('"')
    if text.casefold() in {"", "nan", "none", "na"}:
        return ""
    return re.sub(r"^(ENS[A-Z]*G[0-9]+)\.[0-9]+$", r"\1", text)


def enrichment_background(frame: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Keep measured/testable rows without selecting on significance."""
    eligible = pd.Series(True, index=frame.index)
    reason = pd.Series("source_background", index=frame.index)
    available = [c for c in ("as_available", "deg_available", "dtu_available", "dexseq_available") if c in frame]
    if available:
        measured = frame[available].fillna(False).any(axis=1)
        eligible &= measured
        reason.loc[~measured] = "no_available_evidence"
    if {"model_status", "pvalue"} <= set(frame):
        testable = pd.to_numeric(frame.pvalue, errors="coerce").between(0, 1)
        eligible &= testable
        reason.loc[~testable] = "model_without_valid_test"
    audit = frame[[c for c in ("row_id", "comparison_id", "gene_id", "gene_symbol", "model_status", "pvalue") if c in frame]].copy()
    audit["background_eligible"] = eligible
    audit["background_reason"] = reason
    return frame.loc[eligible].copy(), audit


def canonical_gene_mapping(frame: pd.DataFrame) -> tuple[pd.Series, dict[str, str], pd.DataFrame]:
    """One biological gene per ID, with ambiguous symbol aliases kept out of tests."""
    ids = frame.get("gene_id", pd.Series("", index=frame.index)).map(gene_token)
    symbols = frame.get("gene_symbol", pd.Series("", index=frame.index)).map(gene_token)
    canonical = ids.where(ids.ne(""), symbols).rename("canonical_gene_id")
    aliases: dict[str, set[str]] = {}
    for key, identifier, symbol in zip(canonical, ids, symbols, strict=True):
        if key:
            for alias in {identifier, symbol, key} - {""}:
                aliases.setdefault(alias, set()).add(key)
    mapping = {alias: next(iter(keys)) for alias, keys in aliases.items() if len(keys) == 1}
    audit = pd.DataFrame([
        {"row_id": alias, "alias": alias, "canonical_gene_ids": ";".join(sorted(keys)),
         "mapping_status": "mapped" if len(keys) == 1 else "ambiguous"}
        for alias, keys in sorted(aliases.items())
    ], columns=["row_id", "alias", "canonical_gene_ids", "mapping_status"])
    return canonical, mapping, audit


def fill_gene_symbols(frame: pd.DataFrame, annotation: pd.DataFrame) -> pd.DataFrame:
    """Fill missing names from an unambiguous project GTF; preserve row identity."""
    result = frame.copy()
    ids = result.get("gene_id", pd.Series("", index=result.index)).map(gene_token)
    symbols = result.get("gene_symbol", pd.Series("", index=result.index)).map(gene_token)
    mapping = annotation[["gene_id", "gene_symbol"]].map(gene_token).drop_duplicates()
    mapping = mapping[mapping.gene_id.ne("") & mapping.gene_symbol.ne("")]
    mapping = mapping[~mapping.gene_id.duplicated(keep=False)].set_index("gene_id").gene_symbol
    missing = symbols.eq("") | symbols.eq(ids)
    resolved = ids.map(mapping)
    fill = missing & resolved.notna()
    result["gene_symbol"] = symbols.where(~fill, resolved)
    origin = pd.Series("source", index=result.index)
    origin.loc[missing] = "unresolved"
    origin.loc[fill] = "project_gtf"
    result["gene_symbol_source"] = origin
    return result


def map_gene_sets(gene_sets: dict[str, set[str]], aliases: dict[str, str]) -> tuple[dict[str, set[str]], pd.DataFrame]:
    mapped, audit = {}, []
    for term, members in gene_sets.items():
        tokens = {gene_token(member) for member in members} - {""}
        mapped[term] = {aliases[token] for token in tokens if token in aliases}
        audit.append({"row_id": term, "term": term, "input_identifiers": len(tokens),
                      "mapped_unique_genes": len(mapped[term]),
                      "unmapped_identifiers": ";".join(sorted(tokens - aliases.keys()))})
    return mapped, pd.DataFrame(audit)


def read_gene_sets(path: Path) -> dict[str, set[str]]:
    sets: dict[str, set[str]] = {}
    if path.suffix.lower() == ".gmt":
        for number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), start=1):
            if not line.strip():
                continue
            parts = [part.strip() for part in line.split("\t")]
            if len(parts) < 3 or not parts[0]:
                raise ValueError(f"Malformed GMT row {number}: {path}")
            if parts[0] in sets:
                raise ValueError(f"Duplicate GMT term {parts[0]}: {path}")
            sets[parts[0]] = {gene for gene in parts[2:] if gene}
        return sets
    frame = pd.read_csv(path, sep=None, engine="python")
    set_column = next((column for column in ("set", "term", "pathway", "rbp") if column in frame.columns), None)
    gene_column = next((column for column in ("gene", "gene_symbol", "symbol", "target") if column in frame.columns), None)
    if set_column is None or gene_column is None:
        raise ValueError(f"Gene-set file requires set/term and gene/symbol columns: {path}")
    for term, group in frame.groupby(set_column):
        sets[str(term)] = set(group[gene_column].dropna().astype(str).str.strip())
    return sets


def declared_gene_sets(manifest, kind: str) -> tuple[dict[str, set[str]], pd.DataFrame]:
    """Load every declared library with namespaced terms, never an arbitrary first file."""
    combined, rows = {}, []
    for asset in manifest.assets:
        if asset.kind != kind:
            continue
        if not asset.path.is_file():
            raise FileNotFoundError(f"Declared annotation library is missing: {asset.path}")
        species = asset.metadata.get("species")
        build = manifest.annotation.genome_build
        if species == "Mus_musculus" and build in {"hg38", "GRCh38", "hg19", "GRCh37"} or species == "Homo_sapiens" and build in {"mm39", "GRCm39", "mm10", "GRCm38"}:
            raise ValueError(f"Annotation species does not match project: {asset.asset_id}")
        sets = read_gene_sets(asset.path)
        if not sets:
            raise ValueError(f"Empty gene-set library: {asset.path}")
        combined.update({f"{asset.asset_id} / {term}": genes for term, genes in sets.items()})
        from katana_next.resources import sha256
        digest = sha256(asset.path)
        expected = asset.metadata.get("sha256")
        if expected and digest != expected:
            raise ValueError(f"Annotation content changed; register a new version: {asset.path}")
        rows.append({"row_id": asset.asset_id, "asset_id": asset.asset_id, "path": str(asset.path),
                     "species": species or "user declared / unspecified", "sha256": digest,
                     "terms": len(sets), "semantics": asset.metadata.get("semantics", "custom gene sets")})
    if not combined:
        raise ValueError(f"No {kind} library configured. Prepare or register a versioned resource bundle first.")
    return combined, pd.DataFrame(rows)


def overrepresentation(
    selected: Iterable[str],
    background: Iterable[str],
    gene_sets: dict[str, set[str]],
    *,
    minimum_overlap: int = 2,
    minimum_term_size: int = 1,
    maximum_term_size: int | None = None,
) -> pd.DataFrame:
    selected_set = {str(value) for value in selected if str(value)}
    background_set = {str(value) for value in background if str(value)}
    selected_set &= background_set
    rows = []
    for term, members in gene_sets.items():
        tested_members = members & background_set
        overlap = selected_set & tested_members
        if len(tested_members) < minimum_term_size or (maximum_term_size is not None and len(tested_members) > maximum_term_size):
            continue
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
                "genes": ", ".join(sorted(overlap)),
                "passes_display_overlap": len(overlap) >= minimum_overlap,
            }
        )
    result = pd.DataFrame(rows)
    if not result.empty:
        from scipy.stats import hypergeom
        result["pvalue"] = hypergeom.sf(result.overlap_genes.to_numpy() - 1, len(background_set),
                                        len(selected_set), result.term_genes_in_background.to_numpy())
        result["fdr"] = bh_adjust(result["pvalue"])
        result["tested_term_count"] = len(result)
        result = result.sort_values(["fdr", "fold_enrichment"], ascending=[True, False])
    return result


def library_overrepresentation(selected, tested, gene_sets, *, minimum_overlap=2):
    """GO/pathway ORA conditions on being annotated in that declared library."""
    libraries = {}
    for term, genes in gene_sets.items():
        libraries.setdefault(term.split(" / ")[0], {})[term] = genes
    rows = []
    for name, sets in libraries.items():
        background = set(tested) & set().union(*sets.values())
        if not background:
            raise ValueError(f"No tested genes map to annotation library {name}")
        part = overrepresentation(selected, background, sets, minimum_overlap=minimum_overlap)
        part["library"] = name
        part["tested_genes_before_annotation"] = len(tested)
        part["selected_genes_before_annotation"] = len(set(selected) & set(tested))
        part["background_rule"] = "tested and annotated genes in this library; unique canonical gene IDs"
        rows.append(part)
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()
