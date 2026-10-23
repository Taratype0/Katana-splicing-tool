"""Portable, completed-evidence target requests; no analysis or catalog on import.

The request's content-derived ID covers the metadata and complete member ledger.
This detects accidental alteration; it is not an authentication signature. A
publisher can additionally distribute the ZIP SHA256 through a trusted channel.
"""
from __future__ import annotations

from contextlib import contextmanager
from dataclasses import dataclass
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path, PurePosixPath, PureWindowsPath
import re
import stat
import tempfile
from typing import Any
from urllib.parse import urlsplit
from zipfile import ZipFile, ZIP_DEFLATED

import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq


SCHEMA = "katana.target_request/v1"
DEFAULT_MODULES = ("gene_card", "event_card", "sashimi")
MODULES = frozenset(DEFAULT_MODULES)
_MAX_MEMBERS = 10000
_MAX_BYTES = 8 * 1024**3
_GENE_COLUMNS = ("gene_id", "GeneID")
_EVENT_COLUMNS = ("event_id", "event_uid", "event_key")
_MODULE_COLUMNS = ("module_id", "module", "module_label")
_GEOMETRY_COLUMNS = (
    "chr", "chrom", "chromosome", "strand", "exonStart_0base", "exonEnd",
    "upstreamES", "upstreamEE", "downstreamES", "downstreamEE",
    "riExonStart_0base", "riExonEnd", "longExonStart_0base", "longExonEnd",
    "shortES", "shortEE", "flankingES", "flankingEE", "1stExonStart_0base",
    "1stExonEnd", "2ndExonStart_0base", "2ndExonEnd", "start", "end",
)


@dataclass(frozen=True)
class Bundle:
    request: dict
    root: Path


def _canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"),
                      ensure_ascii=True, allow_nan=False).encode("utf-8")


def _digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024**2), b""):
            h.update(block)
    return h.hexdigest()


def _request_id(request):
    return "targets-" + hashlib.sha256(_canonical({k: v for k, v in request.items()
                                                 if k != "request_id"})).hexdigest()


def _absolute(text):
    return (PureWindowsPath(text).is_absolute() or bool(PureWindowsPath(text).drive)
            or PurePosixPath(text).is_absolute() or text.startswith("~"))


def _portable(value, resource_paths=None):
    """Retain metadata without shipping machine paths or embedded credentials."""
    resource_paths = resource_paths or {}
    if isinstance(value, Path):
        value = str(value)
    if isinstance(value, dict):
        return {str(k): _portable(v, resource_paths) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_portable(v, resource_paths) for v in value]
    if isinstance(value, str):
        if _absolute(value):
            return "resource:" + resource_paths[value] if value in resource_paths else (
                "local-reference:" + hashlib.sha256(value.encode()).hexdigest())
        if "://" in value:
            parsed = urlsplit(value)
            if parsed.username or parsed.password or parsed.query or parsed.fragment:
                return "redacted-resource-uri:" + hashlib.sha256(value.encode()).hexdigest()
            return value
        # A method note can contain a path as part of a sentence. The original
        # source remains immutable; only the portable copy is sanitized.
        return re.sub(r"(?:[A-Za-z]:[\\/]|\\\\)[^\s\"'<>;,]+|(?<!\w)/(?:[\w.-]+/)+[^\s\"'<>;,]*",
                      lambda m: "local-reference:" + hashlib.sha256(m[0].encode()).hexdigest(), value)
    if value is None or isinstance(value, (bool, int)):
        return value
    if isinstance(value, float):
        if not __import__("math").isfinite(value):
            raise ValueError("Portable request metadata must contain finite JSON numbers")
        return value
    if hasattr(value, "item"):
        return _portable(value.item(), resource_paths)
    raise ValueError(f"Unsupported portable metadata type: {type(value).__name__}")


def _json(path):
    def pairs(items):
        values = {}
        for k, v in items:
            if k in values:
                raise ValueError(f"Duplicate JSON key: {k}")
            values[k] = v
        return values
    return json.loads(Path(path).read_text(encoding="utf-8-sig"), object_pairs_hook=pairs,
                      parse_constant=lambda v: (_ for _ in ()).throw(ValueError(f"Nonfinite JSON: {v}")))


def _relative(name):
    if not isinstance(name, str) or not name or "\\" in name or ":" in name or "\x00" in name:
        raise ValueError("Member paths must be portable relative paths")
    p = PurePosixPath(name)
    if p.is_absolute() or any(v in {"", ".", ".."} for v in name.split("/")) or name.startswith("~"):
        raise ValueError("Member path traversal or absolute path")
    reserved = {"CON", "PRN", "AUX", "NUL", *(f"COM{i}" for i in range(1, 10)), *(f"LPT{i}" for i in range(1, 10))}
    if any(part.endswith((".", " ")) or part.split(".")[0].upper() in reserved for part in p.parts):
        raise ValueError("Ambiguous or reserved cross-platform member path")
    if p.suffix.lower() in {".bam", ".bai", ".cram", ".crai", ".fa", ".fasta", ".gtf"}:
        raise ValueError("BAMs and large references must remain external resources")
    return p


def _member(root, name):
    p = _relative(name)
    root = Path(root).resolve()
    current = root
    for part in p.parts:
        current = current / part
        if current.is_symlink():
            raise ValueError("Bundle members cannot be symlinks")
    if not current.resolve().is_relative_to(root) or not current.is_file():
        raise ValueError(f"Missing or unsafe bundle member: {name}")
    if not stat.S_ISREG(current.stat().st_mode):
        raise ValueError("Bundle member must be a regular file")
    return current


def _column(table, candidates):
    return next((c for c in candidates if c in table.column_names), None)


def _in(table, column, values):
    if not column or not values:
        return pa.array([False] * table.num_rows)
    return pc.fill_null(pc.is_in(pc.cast(table[column], pa.string()),
                                value_set=pa.array(sorted(values), type=pa.string())), False)


def _scope(table, column, allowed):
    return pc.or_(_in(table, column, allowed | {""}), pc.is_null(table[column]))


def _target_mask(table, target, events):
    """Gene IDs are authoritative; shared loci/symbols do not alias genes."""
    gc, ec = _column(table, _GENE_COLUMNS), _column(table, _EVENT_COLUMNS)
    if gc:
        gene = _in(table, gc, {target["gene_id"]})
        if target["target_kind"] == "gene" or not ec:
            return gene
        return pc.and_(gene, _in(table, ec, events))
    return _in(table, ec, events)


def _value(row, columns):
    for col in columns:
        value = row.get(col)
        if value is not None and str(value).strip():
            return str(value)
    return ""


def _mode(artifact):
    for scope in (artifact.parameters, artifact.metadata):
        for key in ("count_mode", "rmats_mode", "selected_rmats_mode", "mode"):
            value = str(scope.get(key, "")).upper()
            if value in {"JC", "JCEC"}:
                return value
    return None


def _verified_source(context, artifact):
    from katana_next.infrastructure.fingerprints import fingerprint_file, fingerprint_value
    if artifact is None or artifact.kind != "table" or artifact.status != "complete":
        raise ValueError("Target source must be a completed table")
    if (artifact.project_id != context.manifest.project_id or not artifact.path.is_file() or
            artifact.metadata.get("project_snapshot_fingerprint") != context.snapshot.manifest_fingerprint):
        raise ValueError("Selected source is unavailable, foreign or stale for the current project")
    before = artifact.path.stat()
    digest = _digest(artifact.path)
    if artifact.metadata.get("artifact_contract_version") == "3.0-correctness":
        expected = fingerprint_value({"run": artifact.run_fingerprint, "name": artifact.name,
            "path": artifact.path.name,
            "content": fingerprint_file(artifact.path, force_full=True).model_dump(mode="json")})
        if expected != artifact.fingerprint:
            raise ValueError("Completed source content no longer matches its registered fingerprint")
        marker_path = (context.manifest.output_root / "artifacts" /
                       re.sub(r"[^\w.-]", "_", artifact.analysis_id) /
                       artifact.run_fingerprint[:16] / "parameters.json")
        marker = _json(marker_path)
        if (marker.get("completed") is not True or marker.get("workflow_id") != artifact.analysis_id
                or marker.get("workflow_version") != artifact.metadata.get("workflow_version")
                or marker.get("run_fingerprint") != artifact.run_fingerprint
                or fingerprint_value(marker.get("parameters")) != fingerprint_value(artifact.parameters)):
            raise ValueError("Source completion marker does not bind its original workflow, version and run")
        matching = [a for a in marker.get("generated_artifacts", [])
                    if a.get("artifact_id") == artifact.artifact_id]
        if len(matching) != 1 or matching[0].get("fingerprint") != artifact.fingerprint:
            raise ValueError("Source completion marker does not authenticate the selected artifact")
    table = pq.read_table(artifact.path)
    after = artifact.path.stat()
    if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns) or _digest(artifact.path) != digest:
        raise ValueError("Source changed during target export")
    return table, digest


def _resources(context, samples):
    resources, bindings, path_ids = [], {}, {}
    snapshots = {str(p.path): p for p in context.snapshot.input_fingerprints}

    def add(kind, path, sample_id=None):
        rid = "resource-" + hashlib.sha256(_canonical([kind, sample_id or "reference"])).hexdigest()[:24]
        fp = snapshots.get(str(path.resolve()))
        entry = {"resource_id": rid, "kind": kind, "sha256": None, "bytes": None,
                 "verification": {"state": "not_verified", "hash_method": "unavailable"}}
        if sample_id:
            entry["sample_id"] = sample_id
        if fp and fp.exists:
            entry["bytes"] = fp.size
            entry["verification"] = {"state": "source_snapshot", "hash_method": fp.hash_method,
                                     "digest": fp.digest}
            if fp.hash_method == "sha256-full":
                entry["sha256"] = fp.digest
        resources.append(entry)
        bindings[rid] = str(path.resolve())
        path_ids[str(path)] = rid
        path_ids[str(path.resolve())] = rid
        return rid

    sample_entries = []
    for sample in samples:
        row = {"sample_id": sample.sample_id, "biological_sample_id": sample.biological_sample_id,
               "condition": sample.condition, "attributes": _portable(sample.attributes)}
        if sample.batch is not None:
            row["batch"] = sample.batch
        for kind in ("bam", "bai", "cram", "crai"):
            if sample.assets.get(kind):
                row[kind + "_resource_id"] = add(kind, sample.assets[kind], sample.sample_id)
        sample_entries.append(row)
    references = {"genome_build": context.manifest.annotation.genome_build,
                  "coordinate_system": context.manifest.metadata.get("coordinate_system", "0-based-half-open")}
    for kind in ("gtf", "genome_fasta", "transcript_fasta"):
        path = getattr(context.manifest.annotation, kind)
        if path:
            references[kind + "_resource_id"] = add(kind, path)
    return resources, sample_entries, references, bindings, path_ids


def _sanitize_table(table, paths):
    def redact(value):
        if isinstance(value, str):
            converted = _portable(value, paths)
            return converted, int(converted != value)
        if isinstance(value, (list, tuple)):
            converted = [redact(v) for v in value]
            return [v for v, _ in converted], sum(n for _, n in converted)
        if isinstance(value, dict):
            converted = {k: redact(v) for k, v in value.items()}
            return {k: v for k, (v, _) in converted.items()}, sum(n for _, n in converted.values())
        return value, 0

    redactions = {}
    for i, field in enumerate(table.schema):
        typ = field.type.value_type if pa.types.is_dictionary(field.type) else field.type
        if not (pa.types.is_string(typ) or pa.types.is_large_string(typ) or pa.types.is_list(typ)
                or pa.types.is_large_list(typ) or pa.types.is_struct(typ) or pa.types.is_map(typ)):
            continue
        values = table.column(i).to_pylist()
        converted_counts = [redact(v) for v in values]
        converted = [v for v, _ in converted_counts]
        count = sum(n for _, n in converted_counts)
        if count:
            redactions[field.name] = count
            table = table.set_column(i, field, pa.array(converted, type=field.type))
    # Arrow schema metadata from external writers can also contain paths.
    metadata = table.schema.metadata
    if metadata:
        metadata = {k: _portable(v.decode("utf-8"), paths).encode("utf-8")
                    for k, v in metadata.items()}
        table = table.replace_schema_metadata(metadata)
    return table, redactions


def export_targets(context, selections, output_dir, *, modules=None, contrast_ids=None):
    """Export exact selections plus completed related evidence, never raw reads."""
    if not isinstance(selections, (list, tuple)) or not selections:
        raise ValueError("At least one exact source selection is required")
    from katana_next.infrastructure.manifest_store import ManifestStore
    _, current_manifest = ManifestStore().load(context.manifest_path)
    if current_manifest.model_dump(mode="json") != context.manifest.model_dump(mode="json"):
        raise ValueError("Source project manifest changed after opening this project")
    selected_modules = list(DEFAULT_MODULES if modules is None else modules)
    if not selected_modules or len(set(selected_modules)) != len(selected_modules) or set(selected_modules) - MODULES:
        raise ValueError("Choose unique supported target modules")
    known_contrasts = {c.contrast_id: c for c in context.manifest.contrasts}
    explicit = list(contrast_ids) if contrast_ids is not None else None
    if explicit is not None and (not explicit or len(set(explicit)) != len(explicit) or set(explicit) - known_contrasts.keys()):
        raise ValueError("Unknown or duplicate requested contrast IDs")
    tables, sources, source_hashes, selection_records, targets = {}, {}, {}, [], {}
    observed_contrasts = set()
    for spec in selections:
        if not isinstance(spec, dict) or spec.get("target_kind") not in {"gene", "event"}:
            raise ValueError("Each selection requires target_kind gene or event")
        row_ids = spec.get("row_ids")
        if (not isinstance(row_ids, (list, tuple)) or not row_ids or
                any(not isinstance(v, str) or not v for v in row_ids) or len(set(row_ids)) != len(row_ids)):
            raise ValueError("Selection row_ids must be nonempty, unique exact string identities")
        artifact = context.catalog.get_artifact(spec.get("source_artifact_id", ""))
        if artifact is None:
            raise ValueError("Unknown selection source artifact")
        if artifact.artifact_id not in tables:
            table, digest = _verified_source(context, artifact)
            tables[artifact.artifact_id], source_hashes[artifact.artifact_id] = table, digest
            sources[artifact.artifact_id] = artifact
        table = tables[artifact.artifact_id]
        if "row_id" not in table.column_names:
            raise ValueError("Selection source requires stable row_id")
        rows = table.filter(pc.is_in(pc.cast(table["row_id"], pa.string()), value_set=pa.array(row_ids)))
        found = rows["row_id"].to_pylist()
        if len(found) != len(row_ids) or set(found) != set(row_ids):
            raise ValueError("Selected rows are absent or duplicated in their exact source")
        filters = spec.get("filters", {})
        if not isinstance(filters, dict):
            raise ValueError("Browsing filters must be a JSON object")
        record = {"source_artifact_id": artifact.artifact_id, "analysis_id": artifact.analysis_id,
                  "row_ids": list(row_ids), "source_parameters": artifact.parameters,
                  "browse_filters": filters, "source_sha256": source_hashes[artifact.artifact_id],
                  "source_fingerprint": artifact.fingerprint, "count_mode": _mode(artifact)}
        selection_records.append(record)
        for row in rows.to_pylist():
            gene_id = _value(row, _GENE_COLUMNS)
            if not gene_id:
                raise ValueError("Exact target selection requires canonical gene_id; symbols alone are ambiguous")
            comparison = _value(row, ("comparison_id", "contrast_id"))
            if comparison:
                if comparison not in known_contrasts:
                    raise ValueError("Selected row refers to a foreign comparison")
                observed_contrasts.add(comparison)
                if explicit is not None and comparison not in explicit:
                    raise ValueError("Selected row lies outside the explicitly requested contrasts")
            key = (spec["target_kind"], gene_id)
            values = {"target_kind": spec["target_kind"], "gene_id": gene_id,
                      "gene_symbol": _value(row, ("gene_symbol", "geneSymbol")) or gene_id}
            if spec["target_kind"] == "event":
                event_id, event_type = _value(row, _EVENT_COLUMNS), _value(row, ("event_type",))
                geometry = {k: row[k] for k in _GEOMETRY_COLUMNS if row.get(k) is not None}
                if not event_id or not event_type or not _value(geometry, ("chr", "chrom", "chromosome")) or geometry.get("strand") not in {"+", "-"}:
                    raise ValueError("Exact event target requires event identity, type, chromosome and strand geometry")
                values.update(event_id=event_id, event_type=event_type, event_geometry=geometry)
                key = ("event", event_id)
                if key in targets and any(targets[key].get(k) != values[k] for k in ("gene_id", "event_type", "event_geometry")):
                    raise ValueError("The same event ID has conflicting gene or geometry identities")
            if key not in targets:
                tid = values["target_kind"] + "-" + hashlib.sha256(_canonical(values)).hexdigest()[:24]
                targets[key] = {"target_id": tid, **values, "comparison_ids": [], "source_selections": []}
            if record not in targets[key]["source_selections"]:
                targets[key]["source_selections"].append(record)
    requested = explicit or [c for c in known_contrasts if c in observed_contrasts] or list(known_contrasts)
    if not requested:
        raise ValueError("Target requests require declared comparisons")
    for target in targets.values():
        target["comparison_ids"] = requested.copy()
    contrasts = [known_contrasts[c].model_dump(mode="json") for c in requested]
    sample_ids = {s for c in contrasts for a in ("numerator", "denominator") for s in c[a]["sample_ids"]}
    samples = [s for s in context.manifest.samples if s.sample_id in sample_ids]
    resources, sample_rows, references, bindings, path_ids = _resources(context, samples)
    for target in targets.values():
        target["source_selections"] = _portable(target["source_selections"], path_ids)
    excluded = []
    for artifact in context.catalog.list_artifacts(context.manifest.project_id, kind="table"):
        if artifact.artifact_id in tables:
            continue
        try:
            table, digest = _verified_source(context, artifact)
        except (ValueError, OSError, KeyError) as exc:
            excluded.append({"source_artifact_id": artifact.artifact_id, "reason": _portable(str(exc))})
            continue
        tables[artifact.artifact_id], source_hashes[artifact.artifact_id] = table, digest
        sources[artifact.artifact_id] = artifact
    target_events = {t["target_id"]: ({t["event_id"]} if t["target_kind"] == "event" else set()) for t in targets.values()}
    target_modules = {t["target_id"]: set() for t in targets.values()}
    # Discover every measured event for gene targets before slicing a wide PSI
    # matrix or module-context table; no significance or Top-N gate is applied.
    for table in tables.values():
        gc, ec = _column(table, _GENE_COLUMNS), _column(table, _EVENT_COLUMNS)
        if gc and ec:
            for target in targets.values():
                if target["target_kind"] != "gene":
                    continue
                mask = _in(table, gc, {target["gene_id"]})
                cc = _column(table, ("comparison_id", "contrast_id"))
                if cc:
                    mask = pc.and_(mask, _scope(table, cc, set(requested)))
                target_events[target["target_id"]].update(str(v) for v in table[ec].filter(mask).to_pylist() if v is not None)
    for table in tables.values():
        gc, ec = _column(table, _GENE_COLUMNS), _column(table, _EVENT_COLUMNS)
        mc = _column(table, _MODULE_COLUMNS)
        if mc and (gc or ec):
            for target in targets.values():
                mask = _target_mask(table, target, target_events[target["target_id"]])
                target_modules[target["target_id"]].update(str(v) for v in table[mc].filter(mask).to_pylist() if v is not None)
    root = Path(output_dir).expanduser().resolve()
    archive = root.with_suffix(root.suffix + ".zip")
    template = root.with_name(root.name + ".local-bindings.template.json")
    if archive.exists() or template.exists():
        raise FileExistsError("Target archive or bindings template already exists; choose a new export directory")
    if root.exists() and (not root.is_dir() or any(root.iterdir())):
        raise FileExistsError("Target export requires a new or empty owned directory")
    root.mkdir(parents=True, exist_ok=True)
    claim = root / ".export.claim"
    with claim.open("x", encoding="ascii") as handle:
        handle.write(str(os.getpid()))
    evidence, members = [], {}
    for aid, table in tables.items():
        artifact = sources[aid]
        gc, ec, mc = (_column(table, c) for c in (_GENE_COLUMNS, _EVENT_COLUMNS, _MODULE_COLUMNS))
        keep, matches, target_masks = _in(table, None, set()), set(), {}
        for target in targets.values():
            tid = target["target_id"]
            related = _target_mask(table, target, target_events[tid])
            if not gc and not ec:
                if mc:
                    related = _in(table, mc, target_modules[tid])
                elif "sample_id" in table.column_names:
                    related = _in(table, "sample_id", sample_ids)
            if not pc.any(related).as_py():
                continue
            for column, known, chosen in (("comparison_id", set(known_contrasts), set(requested)),
                                           ("contrast_id", set(known_contrasts), set(requested)),
                                           ("sample_id", {s.sample_id for s in context.manifest.samples}, sample_ids)):
                if column not in table.column_names:
                    continue
                values = {str(v) for v in table[column].filter(related).to_pylist() if v is not None and v != ""}
                if values - known:
                    raise ValueError(f"Related evidence refers to a foreign {column}")
                related = pc.and_(related, _scope(table, column, chosen))
            if pc.any(related).as_py():
                matches.add(tid)
                keep = pc.or_(keep, related)
                target_masks[tid] = related
        if not matches:
            continue
        subset = table.filter(keep)
        export_positions = {int(source_row): position for position, source_row in
                            enumerate(pc.indices_nonzero(keep).to_pylist())}
        target_row_indices = {tid: [export_positions[int(source_row)] for source_row in
                                   pc.indices_nonzero(mask).to_pylist()]
                              for tid, mask in target_masks.items()}
        excluded_samples = set(s.sample_id for s in context.manifest.samples) - sample_ids
        if excluded_samples:
            subset = subset.select([c for c in subset.column_names if c not in excluded_samples])
        subset, redactions = _sanitize_table(subset, path_ids)
        eid = "evidence-" + hashlib.sha256(_canonical([aid, source_hashes[aid]])).hexdigest()[:24]
        rel = f"evidence/{eid}.parquet"
        path = root / rel
        path.parent.mkdir(exist_ok=True)
        pq.write_table(subset, path, compression="zstd")
        details = {"sha256": _digest(path), "bytes": path.stat().st_size}
        members[rel] = details
        metadata = _portable(artifact.metadata, path_ids)
        metadata.update({"source_artifact_id": aid, "source_table_rows": table.num_rows,
                         "exported_rows": subset.num_rows, "subset_without_statistical_refit": True,
                         "gene_event_expansion": "all measured rows; no significance or Top-N gate",
                         "count_mode": _mode(artifact), "path_redactions": redactions,
                         "target_row_indices": target_row_indices})
        metadata["observed_source_policy"] = _portable({
            "read_support_policy": artifact.metadata.get("read_support_policy"),
            "das_selection": artifact.metadata.get("das_selection"),
            "workflow_version": artifact.metadata.get("workflow_version")}, path_ids)
        source_policy = artifact.metadata.get("read_support_policy")
        metadata["source_policy_compatible_with_current"] = None
        if source_policy is not None:
            from katana_next.science.bulk import READ_SUPPORT_POLICY_VERSION
            metadata["source_policy_compatible_with_current"] = source_policy == READ_SUPPORT_POLICY_VERSION
        if artifact.metadata.get("artifact_contract_version") == "3.0-correctness":
            metadata["completion_verification"] = "verified_completion_marker"
        else:
            metadata["completion_verification"] = "legacy_declared"
            metadata["completion_caveat"] = ("Legacy source was declared complete in the current project; "
                "its exported bytes were SHA256 verified. Original execution terminal evidence was not recovered.")
        evidence.append({"evidence_id": eid, "analysis_id": artifact.analysis_id, "name": artifact.name,
                         "entity_type": artifact.entity_type, "path": rel, **details,
                         "source_sha256": source_hashes[aid], "source_fingerprint": artifact.fingerprint,
                         "parameters": _portable(artifact.parameters, path_ids), "metadata": metadata,
                         "target_ids": sorted(matches)})
    if not evidence:
        raise ValueError("No completed evidence for the selected targets")
    request = {"schema": SCHEMA, "created_at": datetime.now(timezone.utc).isoformat(),
               "project": {"project_id": context.manifest.project_id, "name": context.manifest.name,
                           "assay": context.manifest.metadata.get("assay_mode", "bulk_rna"),
                           "source_manifest_sha256": _digest(context.manifest_path),
                           "source_snapshot_fingerprint": context.snapshot.manifest_fingerprint,
                           "references": references, "samples": sample_rows, "contrasts": contrasts},
               "targets": list(targets.values()), "evidence": evidence, "resources": resources,
               "modules": selected_modules, "members": members,
               "execution": {"implicit_upstream_analysis": False, "statistics_refit": False,
                             "resource_binding_schema": "katana.target_bindings/v1",
                             "excluded_sources": excluded}}
    request["request_id"] = _request_id(request)
    request_path = root / "request.json"
    request_path.write_text(json.dumps(request, ensure_ascii=False, indent=2, allow_nan=False), encoding="utf-8")
    claim.unlink()
    validate_request(request, root)
    if archive.exists():
        raise FileExistsError("Target archive already exists; choose a new export directory")
    with ZipFile(archive, "x", compression=ZIP_DEFLATED) as z:
        for name in ["request.json", *sorted(members)]:
            z.write(root / name, name)
    with template.open("x", encoding="utf-8") as handle:
        json.dump({"schema": "katana.target_bindings/v1", "resources": bindings}, handle, indent=2)
    return {"request_id": request["request_id"], "request_path": str(request_path),
            "bundle_directory": str(root), "archive_path": str(archive),
            "local_bindings_template": str(template), "target_count": len(targets),
            "targets": request["targets"], "resources": resources}


def validate_request(request, root, *, verify_files=True):
    """Validate portable identities, member ledger and content before any run."""
    if not isinstance(request, dict) or request.get("schema") != SCHEMA:
        raise ValueError("Unsupported target request schema")
    for key in ("targets", "evidence", "resources", "modules"):
        if not isinstance(request.get(key), list):
            raise ValueError(f"Target request {key} must be a JSON array")
    if not isinstance(request.get("project"), dict) or not isinstance(request.get("members"), dict):
        raise ValueError("Target project and member ledger must be JSON objects")
    for key in ("targets", "evidence", "resources"):
        if any(not isinstance(item, dict) for item in request[key]):
            raise ValueError(f"Target request {key} entries must be JSON objects")
    if request.get("request_id") != _request_id(request):
        raise ValueError("Target request content identity mismatch")
    if not request.get("modules") or set(request["modules"]) - MODULES or len(set(request["modules"])) != len(request["modules"]):
        raise ValueError("Invalid target module selection")
    if _portable(request) != request:
        raise ValueError("Portable requests cannot contain local absolute paths or credentials")
    project = request.get("project", {})
    if not project.get("project_id") or not re.fullmatch("[a-f0-9]{64}", str(project.get("source_manifest_sha256", ""))) or not project.get("references", {}).get("genome_build") or not project.get("references", {}).get("coordinate_system"):
        raise ValueError("Missing project or reference identity")
    samples = project.get("samples", [])
    sample_ids = [s.get("sample_id") for s in samples]
    if not samples or any(not s.get("sample_id") or not s.get("biological_sample_id") or not s.get("condition") or not isinstance(s.get("attributes"), dict) for s in samples) or len(set(sample_ids)) != len(sample_ids):
        raise ValueError("Sample identities and biological units must be explicit and unique")
    contrasts = project.get("contrasts", [])
    cids = [c.get("contrast_id") for c in contrasts]
    if not contrasts or any(not c for c in cids) or len(set(cids)) != len(cids):
        raise ValueError("Contrasts must have unique identities")
    for contrast in contrasts:
        n, d = [contrast.get(a, {}).get("sample_ids", []) for a in ("numerator", "denominator")]
        if not n or not d or len(set(n)) != len(n) or len(set(d)) != len(d) or set(n) & set(d) or (set(n) | set(d)) - set(sample_ids):
            raise ValueError("Contrast arms must refer to declared disjoint biological samples")
    resources = request.get("resources", [])
    rids = [r.get("resource_id") for r in resources]
    if len(set(rids)) != len(rids) or any(not r for r in rids):
        raise ValueError("Resource identities must be unique")
    for resource in resources:
        if resource.get("sample_id") and resource["sample_id"] not in sample_ids:
            raise ValueError("Resource belongs to a foreign sample")
        sha = resource.get("sha256")
        if sha is not None and (not re.fullmatch("[a-f0-9]{64}", sha) or resource.get("verification", {}).get("hash_method") != "sha256-full"):
            raise ValueError("A sampled or unavailable resource fingerprint is not a full SHA256")
    for sample in samples:
        for key in ("bam_resource_id", "bai_resource_id", "cram_resource_id", "crai_resource_id"):
            if sample.get(key) and sample[key] not in rids:
                raise ValueError("Sample resource binding is unresolved")
    for key, value in project["references"].items():
        if key.endswith("_resource_id") and value not in rids:
            raise ValueError("Reference resource binding is unresolved")
    targets = request.get("targets", [])
    tids = [t.get("target_id") for t in targets]
    if not targets or len(set(tids)) != len(tids) or any(not t for t in tids):
        raise ValueError("Targets must have unique identities")
    for target in targets:
        if target.get("target_kind") not in {"gene", "event"} or not target.get("gene_id") or not target.get("comparison_ids") or set(target["comparison_ids"]) - set(cids):
            raise ValueError("Target gene identity or comparison scope is invalid")
        if target["target_kind"] == "event" and (not target.get("event_id") or not target.get("event_type") or not target.get("event_geometry")):
            raise ValueError("An exact event requires its native identity and geometry")
        if not target.get("source_selections") or any(not s.get("row_ids") or not s.get("source_sha256") or not isinstance(s.get("source_parameters"), dict) or not isinstance(s.get("browse_filters"), dict) for s in target["source_selections"]):
            raise ValueError("Target source and browsing snapshots must be recorded separately")
    members = request.get("members", {})
    for name, identity in members.items():
        if (not isinstance(identity, dict) or not isinstance(identity.get("bytes"), int) or isinstance(identity.get("bytes"), bool)
                or identity["bytes"] < 0 or not re.fullmatch("[a-f0-9]{64}", str(identity.get("sha256", "")))):
            raise ValueError("Every member requires a complete byte count and full SHA256")
    if not members or len(members) > _MAX_MEMBERS or sum(d.get("bytes", 0) for d in members.values()) > _MAX_BYTES:
        raise ValueError("Invalid or excessive bundle member ledger")
    names = [str(_relative(name)) for name in members]
    if "request.json" in members or len({n.casefold() for n in names}) != len(names):
        raise ValueError("Duplicate, aliased or self-referential bundle members")
    root = Path(root).resolve()
    if verify_files:
        actual = set()
        for path in root.rglob("*"):
            if path.is_symlink():
                raise ValueError("Bundle directories cannot contain symlinks")
            if path.is_file():
                actual.add(path.relative_to(root).as_posix())
        if actual != {"request.json", *members}:
            raise ValueError("Bundle member set does not match the complete ledger")
        on_disk = _json(_member(root, "request.json"))
        if on_disk != request:
            raise ValueError("Request file does not match the validated request")
        for name, identity in members.items():
            path = _member(root, name)
            if path.stat().st_size != identity.get("bytes") or _digest(path) != identity.get("sha256"):
                raise ValueError(f"Bundle member integrity mismatch: {name}")
    evidence = request.get("evidence", [])
    eids = [e.get("evidence_id") for e in evidence]
    if not evidence or len(set(eids)) != len(eids) or any(not e for e in eids):
        raise ValueError("Missing or duplicate completed evidence")
    evidence_paths = []
    for item in evidence:
        name = item.get("path", "")
        _relative(name)
        if name not in members or any(item.get(k) != members[name].get(k) for k in ("sha256", "bytes")):
            raise ValueError("Evidence does not match its member ledger")
        if not item.get("source_sha256") or not item.get("source_fingerprint") or not item.get("target_ids") or set(item["target_ids"]) - set(tids):
            raise ValueError("Evidence source or target identity is incomplete")
        evidence_paths.append(name)
        if verify_files:
            table = pq.read_table(_member(root, name))
            _, redactions = _sanitize_table(table, {})
            if redactions:
                raise ValueError("Portable evidence cannot contain machine absolute paths")
            row_map = item.get("metadata", {}).get("target_row_indices")
            if not isinstance(row_map, dict) or set(row_map) != set(item["target_ids"]):
                raise ValueError("Evidence requires precise per-target row membership")
            for indices in row_map.values():
                if (not isinstance(indices, list) or not indices or len(set(indices)) != len(indices)
                        or any(type(i) is not int or i < 0 or i >= table.num_rows for i in indices)):
                    raise ValueError("Invalid per-target evidence row indices")
            if set(i for indices in row_map.values() for i in indices) != set(range(table.num_rows)):
                raise ValueError("Unassigned unrelated evidence rows are not portable target inputs")
            if "sample_id" in table.column_names:
                present = {v for v in table["sample_id"].to_pylist() if v is not None and v != ""}
                if present - set(sample_ids):
                    raise ValueError("Evidence contains foreign sample identities")
            for column in ("comparison_id", "contrast_id"):
                if column in table.column_names:
                    present = {v for v in table[column].to_pylist() if v is not None and v != ""}
                    if present - set(cids):
                        raise ValueError("Evidence contains foreign contrast identities")
    if set(evidence_paths) != set(members) or len(set(evidence_paths)) != len(evidence_paths):
        raise ValueError("Every evidence member must be described exactly once")
    for target in targets:
        for selection in target["source_selections"]:
            matching = [e for e in evidence if
                        e.get("metadata", {}).get("source_artifact_id") == selection.get("source_artifact_id")
                        and target["target_id"] in e["target_ids"]
                        and e["source_sha256"] == selection["source_sha256"]
                        and e["source_fingerprint"] == selection["source_fingerprint"]]
            if len(matching) != 1:
                raise ValueError("Selected source must be authenticated by its exact portable evidence")
            if verify_files:
                source_table = pq.read_table(_member(root, matching[0]["path"]))
                if "row_id" not in source_table.column_names or not set(selection["row_ids"]).issubset(set(source_table["row_id"].to_pylist())):
                    raise ValueError("Selected exact rows are absent from their portable source evidence")
    return {"schema": SCHEMA, "valid": True, "request_id": request["request_id"],
            "target_count": len(targets), "member_count": len(members), "files_verified": bool(verify_files)}


@contextmanager
def open_bundle(path):
    """Open a moved ZIP/directory/request, verifying all members before yielding."""
    path = Path(path).expanduser().resolve()
    if path.is_file() and path.suffix.lower() == ".zip":
        with tempfile.TemporaryDirectory(prefix="katana-targets-") as temporary:
            root = Path(temporary)
            with ZipFile(path) as z:
                entries = z.infolist()
                if not entries or len(entries) > _MAX_MEMBERS + 1 or sum(e.file_size for e in entries) > _MAX_BYTES:
                    raise ValueError("Empty or excessive target archive")
                names = []
                for entry in entries:
                    _relative(entry.filename)
                    mode = entry.external_attr >> 16
                    if entry.is_dir() or entry.flag_bits & 1 or stat.S_ISLNK(mode) or (stat.S_IFMT(mode) not in {0, stat.S_IFREG}):
                        raise ValueError("Only unencrypted regular-file ZIP members are allowed")
                    names.append(entry.filename)
                if len({n.casefold() for n in names}) != len(names) or "request.json" not in names:
                    raise ValueError("Duplicate ZIP members or missing request")
                request = json.loads(z.read("request.json"))
                if set(names) != {"request.json", *request.get("members", {})}:
                    raise ValueError("Archive member set does not match request ledger")
                for name in names:
                    destination = root / name
                    destination.parent.mkdir(parents=True, exist_ok=True)
                    with z.open(name) as source, destination.open("xb") as output:
                        __import__("shutil").copyfileobj(source, output)
            validate_request(request, root)
            yield Bundle(request, root)
        return
    root = path if path.is_dir() else path.parent
    if path.is_file() and path.name != "request.json":
        raise ValueError("Open a target ZIP, directory or request.json")
    request = _json(_member(root, "request.json"))
    validate_request(request, root)
    yield Bundle(request, root)
