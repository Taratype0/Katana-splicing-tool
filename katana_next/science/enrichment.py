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


DETECTED_BACKGROUND_SOURCES = ("detected_union", "das", "deg", "dtu", "dexseq")
_ASSAY_PREFIX = {"das": "as", "deg": "deg", "dtu": "dtu", "dexseq": "dexseq"}


def native_background_source_provenance(context, source) -> dict:
    """Trace a full native-result background without claiming external certification.

    The package can verify which files and complete local artifacts were read.
    It cannot discover an undisclosed user truncation of an imported file; that
    delivery-completeness limitation is recorded explicitly, never hidden behind
    a `certified=True` declaration. Explicit hit-only/subset inputs are refused.
    """
    from katana_next.infrastructure.fingerprints import fingerprint_file, fingerprint_value
    from katana_next.workflows.sources import validate_source
    native = {"differential.as": "Differential AS events",
              "differential.deg": "Differentially expressed genes",
              "differential.dtu": "DTU and exon-usage evidence"}
    known = {**native, "prioritize.evidence": "Candidate evidence matrix"}

    def check_scope(metadata):
        scope = str(metadata.get("result_scope", metadata.get("source_result_scope", "not_declared"))).casefold()
        if (metadata.get("significant_only") is True or metadata.get("filtered_to_significant") is True
                or metadata.get("full_results") is False
                or scope in {"significant_only", "hits_only", "significant", "selected", "subset", "top_n", "foreground_only"}):
            raise ValueError("Significant-only or subset inputs can be foregrounds only; provide the complete native result background.")
        return scope

    def artifact_record(artifact):
        validate_source(context, artifact)
        if artifact.analysis_id not in known or artifact.name != known[artifact.analysis_id]:
            raise ValueError("A complete native result source is required for background; custom/hit-only tables can be foregrounds only.")
        scope = check_scope(artifact.metadata)
        identity = fingerprint_file(artifact.path, force_full=True)
        observed = fingerprint_value({"run": artifact.run_fingerprint, "name": artifact.name,
                                      "path": artifact.path.name, "content": identity.model_dump(mode="json")})
        if artifact.metadata.get("artifact_contract_version") == "3.0-correctness" and observed != artifact.fingerprint:
            raise ValueError("Background source content changed after its recorded run completed.")
        return {"artifact_id": artifact.artifact_id, "analysis_id": artifact.analysis_id,
                "path": str(artifact.path), "sha256": identity.digest,
                "rows": artifact.row_count, "source_result_scope": scope}

    records = [artifact_record(source)]
    parents = []
    if source.analysis_id == "prioritize.evidence":
        parents = [context.catalog.get_artifact(identifier) for identifier in source.parent_artifact_ids]
        if not parents or any(p is None or p.analysis_id not in native for p in parents):
            raise ValueError("Candidate background has no traceable complete native parent results; rebuild from complete imported result files.")
        records.extend(artifact_record(parent) for parent in parents)
    analyses = {p.analysis_id for p in parents} if parents else {source.analysis_id}
    kinds = set()
    for analysis, values in (("differential.as", {"rmats_result"}), ("differential.deg", {"deg_result"}),
                             ("differential.dtu", {"dtu_result", "dexseq_result"})):
        if analysis in analyses:
            kinds.update(values)
    scopes = set(source.contrast_ids or [c.contrast_id for c in context.manifest.contrasts])
    inputs = []
    snapshots = {str(f.path.resolve()): f for f in context.snapshot.input_fingerprints}
    for asset in context.manifest.assets:
        if asset.kind not in kinds or asset.metadata.get("contrast_id") not in scopes:
            continue
        scope = check_scope(asset.metadata)
        recorded = snapshots.get(str(asset.path.resolve()))
        if recorded is None or fingerprint_file(asset.path, force_full=recorded.hash_method == "sha256-full") != recorded:
            raise ValueError("Native background input changed or was not fingerprinted in this project snapshot.")
        identity = fingerprint_file(asset.path, force_full=True)
        inputs.append({"asset_id": asset.asset_id, "kind": asset.kind, "contrast_id": asset.metadata.get("contrast_id"),
                       "path": str(asset.path), "sha256": identity.digest, "source_result_scope": scope})
    if not inputs:
        raise ValueError("No complete native input file is traceable for this background; a detection flag alone is insufficient.")
    return {"schema": "katana.native_result_background_provenance/v1", "artifacts": records,
            "native_inputs": inputs, "local_content_verification": "file hashes and completed artifact lineage checked",
            "delivery_completeness": "not_independently_verified",
            "completeness_assumption": "all rows in the registered native result inputs; explicit significant-only or subset inputs rejected",
            "external_execution_certification": "not_asserted"}


def _declared_flag(values: pd.Series, name: str) -> pd.Series:
    """Accept typed flags, never Python's truthiness of strings such as 'False'."""
    def flag(value):
        if pd.isna(value):
            return False
        if isinstance(value, str):
            if value.strip().casefold() in {"true", "1"}:
                return True
            if value.strip().casefold() in {"false", "0", ""}:
                return False
            raise ValueError(f"Invalid detected/tested flag in {name}: {value!r}")
        if value in (True, False, 0, 1):
            return bool(value)
        raise ValueError(f"Invalid detected/tested flag in {name}: {value!r}")
    return values.map(flag).astype(bool)


def detected_gene_background(frame: pd.DataFrame, source: str = "detected_union") -> tuple[pd.DataFrame, pd.DataFrame]:
    """Select only observed/testable genes, with an auditable assay scope.

    A GTF, genome-wide gene list, significance-only list or positive-count cutoff
    invented here is not a detection universe. Callers supply their native
    measured/tested flags. Non-significant eligible genes remain in background.
    Library mapping and canonical-gene deduplication happen later, unchanged.
    """
    if source not in DETECTED_BACKGROUND_SOURCES:
        raise ValueError(f"Unknown detected background source: {source}")
    if source == "detected_union":
        flags = []
        assay_masks = []
        for prefix in ("as", "deg", "dtu", "dexseq"):
            columns = [c for c in (f"{prefix}_detected", f"{prefix}_eligible", f"{prefix}_tested", f"{prefix}_available") if c in frame]
            if columns:
                flags.extend(columns)
                assay_masks.append(pd.concat([_declared_flag(frame[c], c) for c in columns], axis=1).all(axis=1))
        if assay_masks:
            # The union is a deliberately labelled multi-assay universe, not a
            # per-assay background and never a reference annotation universe.
            measured = pd.concat(assay_masks, axis=1).any(axis=1)
        elif any(c in frame for c in ("detected", "measurement_eligible", "tested")):
            flags = [c for c in ("detected", "measurement_eligible", "tested") if c in frame]
            measured = pd.concat([_declared_flag(frame[c], c) for c in flags], axis=1).all(axis=1)
        else:
            raise ValueError("Detected enrichment background requires native detected/tested/available evidence; a reference gene list is not a measured background.")
    else:
        prefix = _ASSAY_PREFIX[source]
        flags = [c for c in (f"{prefix}_detected", f"{prefix}_eligible", f"{prefix}_tested", f"{prefix}_available") if c in frame]
        if not flags:
            raise ValueError(f"No {source.upper()} detection/testability evidence is available; another assay or a genome-wide annotation cannot substitute.")
        measured = pd.concat([_declared_flag(frame[c], c) for c in flags], axis=1).all(axis=1)
    if "measurement_eligible" in frame and "measurement_eligible" not in flags:
        measured &= _declared_flag(frame.measurement_eligible, "measurement_eligible")
        flags.append("measurement_eligible")
    if {"model_status", "pvalue"} <= set(frame):
        measured &= pd.to_numeric(frame.pvalue, errors="coerce").between(0, 1)
        failed = frame.model_status.fillna("").astype(str).str.casefold().str.contains(
            r"failed|error|insufficient|not_tested|not_computed|unavailable|skipped", regex=True)
        measured &= ~failed
        if "pvalue" not in flags:
            flags.extend(["model_status", "pvalue"])
    columns = [c for c in ("row_id", "comparison_id", "cell_type", "assay", "gene_id", "gene_symbol") if c in frame]
    audit = frame[columns].copy()
    audit["background_eligible"] = measured
    audit["background_reason"] = measured.map({True: "native_detected_or_testable", False: "not_detected_or_not_testable"})
    audit["background_source"] = source
    audit["detection_evidence_columns"] = ";".join(flags)
    audit["detection_threshold_policy"] = "native caller eligibility; no additional count or significance threshold"
    return frame.loc[measured].copy(), audit


def attach_native_gene_measurement_flags(candidates: pd.DataFrame, events: pd.DataFrame,
                                        deg: pd.DataFrame, usage: pd.DataFrame) -> pd.DataFrame:
    """Carry caller eligibility into a gene table without changing its statistics.

    In particular, the presence of an rMATS row alone cannot override the native
    every-sample read-support flag. This adapter never applies a new count or
    significance threshold, and never interprets a reference annotation as data.
    """
    result = candidates.copy()
    keys = ["comparison_id", "gene_id"]
    target = pd.MultiIndex.from_frame(result[keys])

    def add(column, frame, eligible):
        measured = frame.loc[eligible, keys].drop_duplicates() if len(frame) else pd.DataFrame(columns=keys)
        result[column] = target.isin(pd.MultiIndex.from_frame(measured))

    if len(events):
        if "measurement_eligible" not in events:
            raise ValueError("Candidate background requires native AS measurement eligibility; rebuild the quality input explicitly.")
        eligible = _declared_flag(events.measurement_eligible, "measurement_eligible")
        if "fdr" in events:
            eligible &= pd.to_numeric(events.fdr, errors="coerce").between(0, 1)
        add("as_eligible", events, eligible)
    else:
        result["as_eligible"] = False
    if len(deg):
        add("deg_tested", deg, pd.to_numeric(deg.deg_fdr, errors="coerce").between(0, 1))
    else:
        result["deg_tested"] = False
    for prefix, kind in (("dtu", "dtu_result"), ("dexseq", "dexseq_result")):
        if len(usage):
            kinds = usage.get("source_kind", pd.Series("dtu_result", index=usage.index))
            add(f"{prefix}_tested", usage, kinds.eq(kind) & pd.to_numeric(usage.dtu_fdr, errors="coerce").between(0, 1))
        else:
            result[f"{prefix}_tested"] = False
    return result


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
