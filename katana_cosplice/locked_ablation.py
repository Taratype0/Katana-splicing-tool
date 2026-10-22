from __future__ import annotations

import hashlib
import json
from dataclasses import asdict, is_dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
from scipy.stats import hypergeom

from katana_cosplice.biological_interpretability import gene_from_event_key
from katana_cosplice.utils import bh_adjust


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def json_ready(value: object) -> object:
    if is_dataclass(value):
        return json_ready(asdict(value))
    if isinstance(value, Path):
        return str(value.resolve())
    if isinstance(value, dict):
        return {str(key): json_ready(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_ready(item) for item in value]
    if isinstance(value, np.generic):
        return value.item()
    return value


def write_json(path: Path, payload: dict[str, object]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(json_ready(payload), indent=2, sort_keys=True),
        encoding="utf-8",
    )
    return path


def file_record(path: Path, *, role: str) -> dict[str, object]:
    resolved = path.resolve()
    stat = resolved.stat()
    return {
        "role": role,
        "path": str(resolved),
        "size_bytes": stat.st_size,
        "sha256": sha256_file(resolved),
    }


def read_psi_matrix(path: Path) -> pd.DataFrame:
    matrix = pd.read_csv(path, sep="\t", index_col=0)
    matrix.index = matrix.index.astype(str)
    matrix.columns = matrix.columns.astype(str)
    matrix = matrix.apply(pd.to_numeric, errors="coerce")
    if not matrix.index.is_unique or not matrix.columns.is_unique:
        raise ValueError(f"PSI matrix labels must be unique: {path}")
    return matrix


def select_discovery_events(matrix: pd.DataFrame, max_events: int) -> pd.DataFrame:
    quality = matrix.notna().mean(axis=1).ge(0.75) & matrix.std(
        axis=1, skipna=True
    ).gt(1e-8)
    selected = matrix.loc[quality]
    if len(selected) > max_events:
        selected = selected.loc[
            selected.var(axis=1, skipna=True).nlargest(max_events).index
        ]
    selected = selected.apply(lambda row: row.fillna(row.mean()), axis=1)
    if selected.empty:
        raise ValueError("No discovery events passed missingness and variance QC.")
    return selected


def sample_frame(path: Path, columns: pd.Index, cohort: str) -> pd.DataFrame:
    source = pd.read_csv(path, sep="\t")
    sample_column = "sample_id" if "sample_id" in source else "sample"
    if sample_column not in source:
        raise ValueError(f"Sample metadata requires sample_id or sample: {path}")
    source[sample_column] = source[sample_column].astype(str)
    source = source.set_index(sample_column).reindex(columns.astype(str))
    if source.index.isna().any() or source.isna().all(axis=1).any():
        missing = source.index[source.isna().all(axis=1)].tolist()
        raise ValueError(f"Sample metadata does not cover matrix columns: {missing[:5]}")
    condition = (
        source["phenotype"].fillna(cohort).astype(str)
        if "phenotype" in source
        else pd.Series(cohort, index=source.index, dtype=object)
    )
    return pd.DataFrame(
        {
            "sample_id": columns.astype(str),
            "condition": condition.reindex(columns.astype(str)).to_numpy(),
            "cohort": cohort,
        }
    )


def combine_gmt_files(paths: Iterable[Path], output: Path) -> Path:
    seen: set[tuple[str, tuple[str, ...]]] = set()
    lines: list[str] = []
    for path in paths:
        with path.open(encoding="utf-8") as handle:
            for raw in handle:
                fields = raw.rstrip("\n").split("\t")
                if len(fields) < 3:
                    continue
                genes = tuple(sorted({gene.strip().casefold() for gene in fields[2:] if gene.strip()}))
                key = (fields[0].strip(), genes)
                if not key[0] or not genes or key in seen:
                    continue
                seen.add(key)
                lines.append("\t".join([key[0], fields[1].strip(), *genes]))
    if not lines:
        raise ValueError("No usable terms were found in the supplied GMT files.")
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return output


def load_gmt(path: Path) -> dict[str, set[str]]:
    terms: dict[str, set[str]] = {}
    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            genes = {gene.strip().casefold() for gene in fields[2:] if gene.strip()}
            if genes:
                terms[fields[0]] = genes
    if not terms:
        raise ValueError(f"No usable terms in GMT: {path}")
    return terms


def local_module_enrichment(
    cohort: str,
    method: str,
    assignments: pd.DataFrame,
    gene_sets: dict[str, set[str]],
) -> pd.DataFrame:
    module_column = next(
        name
        for name in ("consensus_module", "wgcna_module", "module")
        if name in assignments
    )
    frame = assignments[["event_key", module_column]].copy()
    frame["module"] = frame[module_column].astype(str)
    frame["gene"] = frame["event_key"].map(gene_from_event_key).str.casefold()
    frame = frame.loc[frame["module"].ne("grey") & frame["gene"].ne("")]
    universe = set(frame["gene"])
    rows: list[dict[str, object]] = []
    for module, group in frame.groupby("module", observed=True, sort=True):
        genes = set(group["gene"])
        for term, raw_members in gene_sets.items():
            members = raw_members & universe
            overlap = genes & members
            pvalue = (
                float(hypergeom.sf(len(overlap) - 1, len(universe), len(members), len(genes)))
                if overlap and universe and members
                else 1.0
            )
            rows.append(
                {
                    "cohort": cohort,
                    "method": method,
                    "module": str(module),
                    "term": term,
                    "module_genes": len(genes),
                    "term_genes_in_universe": len(members),
                    "overlap_count": len(overlap),
                    "overlap_genes": ";".join(sorted(overlap)),
                    "pvalue": pvalue,
                }
            )
    result = pd.DataFrame(rows)
    if not result.empty:
        result["FDR"] = result.groupby(
            ["cohort", "method", "module"], observed=True
        )["pvalue"].transform(bh_adjust)
    return result


def functional_replication(enrichment: pd.DataFrame, fdr: float = 0.05) -> pd.DataFrame:
    columns = [
        "method",
        "discovery_module",
        "best_holdout_module",
        "discovery_significant_terms",
        "holdout_significant_terms",
        "shared_significant_terms",
        "best_term_jaccard",
    ]
    if enrichment.empty:
        return pd.DataFrame(columns=columns)
    significant = enrichment.loc[pd.to_numeric(enrichment["FDR"], errors="coerce").le(fdr)]
    rows: list[dict[str, object]] = []
    for method, all_method_rows in enrichment.groupby("method", observed=True):
        cohorts = all_method_rows["cohort"].astype(str).unique().tolist()
        if "discovery" not in cohorts or "holdout" not in cohorts:
            continue
        method_rows = significant.loc[significant["method"].eq(method)]
        terms = {
            (str(cohort), str(module)): set(group["term"].astype(str))
            for (cohort, module), group in method_rows.groupby(
                ["cohort", "module"], observed=True
            )
        }
        discovery_modules = sorted(
            all_method_rows.loc[
                all_method_rows["cohort"].eq("discovery"), "module"
            ].astype(str).unique()
        )
        holdout_modules = sorted(
            all_method_rows.loc[
                all_method_rows["cohort"].eq("holdout"), "module"
            ].astype(str).unique()
        )
        for module in discovery_modules:
            source = terms.get(("discovery", module), set())
            best_module = ""
            best_target: set[str] = set()
            best_shared: set[str] = set()
            best_score = 0.0
            for target_module in holdout_modules:
                target = terms.get(("holdout", target_module), set())
                shared = source & target
                union = source | target
                score = len(shared) / len(union) if union else 0.0
                if score > best_score:
                    best_module = target_module
                    best_target = target
                    best_shared = shared
                    best_score = score
            rows.append(
                {
                    "method": str(method),
                    "discovery_module": module,
                    "best_holdout_module": best_module,
                    "discovery_significant_terms": len(source),
                    "holdout_significant_terms": len(best_target),
                    "shared_significant_terms": len(best_shared),
                    "best_term_jaccard": best_score,
                }
            )
    return pd.DataFrame(rows, columns=columns)
