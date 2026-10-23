"""Portable selected-target execution, independent of an originating project DB.

This runner reuses saved statistical evidence. The only read calculation is
explicitly requested exact-event native Sashimi. It does not submit jobs.
"""
from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor, as_completed
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import re
import shutil
import uuid
import zipfile

import pandas as pd

from ._execution import (atomic_json, bounded, canonical, fingerprint, ledger,
                         owned_claim, portable, sha256, utc_now, verify_ledger,
                         write_json)
from .render import (html_page, render_index, render_sample_plots, render_table,
                     render_target, sample_measurements)
from .sashimi import GEOMETRY_COLUMNS, generate_sashimi, native_reader


RESULT_SCHEMA = "katana.target_results/v1"
TARGET_SCHEMA = "katana.target_product/v1"
GENE_COLUMNS = ("gene_id", "GeneID", "gene", "geneID", "gene_symbol", "geneSymbol")
EVENT_COLUMNS = ("event_id", "event_uid", "event_key", "event", "ID")
CONTRAST_COLUMNS = ("comparison_id", "contrast_id", "comparison", "contrast")


def _load_bindings(resources):
    if resources is None:
        return {}, None
    if isinstance(resources, (str, Path)):
        path = Path(resources).expanduser().resolve()
        document = json.loads(path.read_text(encoding="utf-8"))
        base = path.parent
    elif isinstance(resources, dict):
        document, base = resources, Path.cwd()
    else:
        raise TypeError("Resource bindings must be a JSON path or dictionary")
    if "schema" in document and document["schema"] != "katana.target_bindings/v1":
        raise ValueError("Unsupported resource bindings schema")
    values = document.get("resources", document)
    if not isinstance(values, dict):
        raise ValueError("Resource bindings must map logical resource IDs to paths")
    result = {}
    for rid, value in values.items():
        if rid in {"schema", "created_at"}:
            continue
        raw = value.get("path") if isinstance(value, dict) else value
        if raw is None or not isinstance(raw, (str, Path)):
            raise ValueError(f"No local resource path was bound for {rid}")
        path = Path(raw).expanduser()
        result[str(rid)] = (base / path).resolve() if not path.is_absolute() else path.resolve()
    return result, base


def _sampled_sha256(path, block_size=1024 * 1024):
    path = Path(path)
    size = path.stat().st_size
    digest = hashlib.sha256()
    digest.update(str(size).encode("ascii") + b"\0")
    with path.open("rb") as stream:
        for offset in sorted(set([0, max(0, size // 2 - block_size // 2), max(0, size - block_size)])):
            digest.update(str(offset).encode("ascii") + b"\0")
            stream.seek(offset)
            digest.update(stream.read(block_size))
    return digest.hexdigest()


def _legacy_sampled_sha256(path, block_size=1024 * 1024):
    """The existing Katana source-snapshot method, not a full file hash."""
    path = Path(path)
    size = path.stat().st_size
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for offset in (0, max(0, size // 2 - block_size // 2), max(0, size - block_size)):
            stream.seek(offset)
            digest.update(stream.read(block_size))
    digest.update(str(size).encode("ascii"))
    return digest.hexdigest()


def _bound_resources(request, bindings):
    unknown = sorted(set(bindings) - {r["resource_id"] for r in request.get("resources", [])})
    if unknown:
        raise ValueError("Bindings contain unknown logical resource IDs: " + ", ".join(unknown))
    result = {}
    for declaration in request.get("resources", []):
        rid = declaration["resource_id"]
        record = {"resource_id": rid, "kind": declaration.get("kind"), "available": False}
        result[rid] = record
        path = bindings.get(rid)
        if path is None:
            record["reason"] = "Resource has not been bound on this execution machine"
            continue
        if not path.is_file():
            record["reason"] = "Bound resource does not exist or is not a regular file"
            continue
        before = path.stat()
        if declaration.get("bytes") is not None and int(declaration["bytes"]) != before.st_size:
            record["reason"] = "Bound resource size does not match the exported identity"
            continue
        expected_sha = declaration.get("sha256")
        verification = declaration.get("verification", {})
        method = str(verification.get("hash_method", ""))
        if expected_sha:
            actual = sha256(path)
            if actual != expected_sha:
                record["reason"] = "Bound resource full SHA256 does not match the exported identity"
                continue
            identity = {"hash_method": "sha256-full", "digest": actual, "identity_verified": True, "full_content_verified": True}
        elif method.lower() in {"sha256", "sha256-full", "full_sha256"} and verification.get("digest"):
            actual = sha256(path)
            if actual != verification["digest"]:
                record["reason"] = "Bound resource full SHA256 does not match the exported verification"
                continue
            identity = {"hash_method": "sha256-full", "digest": actual, "identity_verified": True, "full_content_verified": True}
        else:
            # The source snapshot may intentionally sample a multi-GB BAM.
            # Verify exactly that declared method; never read a whole cohort
            # solely to render selected target views or call a sampled hash full.
            identity = {"hash_method": "sha256-sampled-v1", "digest": _sampled_sha256(path),
                        "identity_verified": False,
                        "full_content_verified": False,
                        "source_verification": portable(verification),
                        "note": "Sampled execution-machine fingerprint; full-file source equality was not verified"}
            if method.lower() in {"sha256-sampled-v1", "sampled-sha256-v1", "sha256-sampled-3x1mib"} and verification.get("digest"):
                sampled = (_legacy_sampled_sha256(path) if method.lower() == "sha256-sampled-3x1mib"
                           else _sampled_sha256(path))
                if sampled != verification["digest"]:
                    record["reason"] = "Bound resource sampled fingerprint does not match"
                    continue
                identity.update(hash_method=method, digest=sampled, sampled_source_match=True,
                                identity_verified=True,
                                note="Source sampled fingerprint matched; unsampled bytes were not verified")
            elif verification.get("digest") and method:
                # An unknown digest algorithm must not be silently certified.
                record["reason"] = f"Unsupported source resource verification method: {method}"
                continue
        after = path.stat()
        if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
            record["reason"] = "Bound resource changed while its content was verified"
            continue
        record.update(available=True, path=str(path), bytes=after.st_size, verification=identity)
    return result


def _resource_public(resources):
    return {rid: {k: v for k, v in record.items() if k != "path"} for rid, record in resources.items()}


def _method_identity():
    source = Path(__file__).parent
    files = {name: sha256(source / name) for name in ("runner.py", "render.py", "sashimi.py", "_execution.py")}
    import src.adapters.native_sashimi_adapter as adapter
    files["native_sashimi_adapter.py"] = sha256(Path(adapter.__file__))
    versions = {}
    for name in ("pandas", "numpy", "matplotlib", "pyarrow", "bamnostic", "pysam"):
        try:
            versions[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            continue
    return {"implementation": "katana.portable_targets/v1", "sources": files,
            "dependencies": versions, "bam_reader": native_reader()}


def _read_table(root, evidence):
    path = bounded(root, evidence["path"])
    if path.stat().st_size != int(evidence["bytes"]) or sha256(path) != evidence["sha256"]:
        raise ValueError(f"Evidence changed since export: {evidence['evidence_id']}")
    if path.suffix.lower() in {".parquet", ".pq"}:
        return pd.read_parquet(path)
    if path.suffix.lower() in {".json", ".jsonl"}:
        return pd.read_json(path, lines=path.suffix.lower() == ".jsonl")
    return pd.read_csv(path, sep="," if path.suffix.lower() == ".csv" else "\t", low_memory=False)


def _select_frame(frame, target):
    """Preserve exact event where present, otherwise gene/context associations."""
    selected = frame
    cids = set(map(str, target.get("comparison_ids", [])))
    for column in CONTRAST_COLUMNS:
        if cids and column in selected.columns:
            selected = selected.loc[selected[column].astype(str).isin(cids)]
            break
    if target["target_kind"] == "event":
        eid = str(target["event_id"])
        event_cols = [c for c in EVENT_COLUMNS if c in selected.columns]
        if event_cols:
            # Match the exporter's declared canonical identity, never another
            # row's historical numeric/UID alias when event_id is present.
            selected = selected.loc[selected[event_cols[0]].astype(str).eq(eid)]
            for name, value in target.get("event_geometry", {}).items():
                if name in selected.columns and value is not None:
                    if pd.api.types.is_numeric_dtype(selected[name]):
                        selected = selected.loc[pd.to_numeric(selected[name], errors="coerce").eq(float(value))]
                    else:
                        selected = selected.loc[selected[name].astype(str).eq(str(value))]
            if "event_type" in selected.columns and target.get("event_type"):
                selected = selected.loc[selected["event_type"].astype(str).eq(str(target["event_type"]))]
            return selected.copy()
    canonical = [c for c in ("gene_id", "GeneID", "geneID") if c in selected.columns]
    genes = canonical or [c for c in ("gene_symbol", "geneSymbol", "gene") if c in selected.columns]
    ids = {str(target["gene_id"])} if canonical else {str(v) for v in (target.get("gene_id"), target.get("gene_symbol")) if v is not None and str(v)}
    if genes and ids:
        mask = selected[genes[0]].astype(str).isin(ids)
        for c in genes[1:]:
            mask |= selected[c].astype(str).isin(ids)
        selected = selected.loc[mask]
    return selected.copy()


def _tables_for_target(request, root, target):
    result = []
    for evidence in request.get("evidence", []):
        if target["target_id"] not in evidence.get("target_ids", []):
            continue
        frame = _read_table(root, evidence)
        positions = evidence.get("metadata", {}).get("target_row_indices", {}).get(target["target_id"])
        if positions is not None:
            if (not isinstance(positions, list) or any(not isinstance(p, int) or isinstance(p, bool) or p < 0 or p >= len(frame) for p in positions)
                    or len(positions) != len(set(positions))):
                raise ValueError("Invalid authenticated target row positions")
            frame = frame.iloc[positions]
        result.append((evidence, _select_frame(frame, target)))
    if not result or not any(len(frame) for _, frame in result):
        raise ValueError("No complete selected evidence rows were exported for this target")
    return result


def _event_cards(request, target, tables):
    found = {}
    ordered_tables = sorted(tables, key=lambda item: len(set(item[1].columns) & set(GEOMETRY_COLUMNS)), reverse=True)
    for evidence, frame in ordered_tables:
        event_column = next((c for c in EVENT_COLUMNS if c in frame.columns), None)
        for _, row in frame.iterrows():
            event_type = row.get("event_type", evidence.get("metadata", {}).get("event_type", target.get("event_type")))
            if event_type not in {"SE", "RI", "A3SS", "A5SS", "MXE"}:
                continue
            if event_column and pd.isna(row[event_column]):
                continue
            eid = str(row[event_column]) if event_column else str(target.get("event_id", ""))
            if not eid:
                continue
            cid = next((str(row[c]) for c in CONTRAST_COLUMNS if c in frame.columns), None)
            if cid is None:
                cid = str(evidence.get("metadata", {}).get("comparison_id", ""))
            if not cid and len(target.get("comparison_ids", [])) == 1:
                cid = str(target["comparison_ids"][0])
            geometry = {c: portable(row[c]) for c in GEOMETRY_COLUMNS if c in frame.columns and pd.notna(row[c])}
            if target["target_kind"] == "event":
                geometry.update(target.get("event_geometry", {}))
            related = [(k, v) for k, v in found.items() if v["event_id"] == eid and v["comparison_id"] == cid
                       and v["event_type"] == event_type
                       and all(v["event_geometry"].get(c) == value for c, value in geometry.items())]
            key = related[0][0] if len(related) == 1 else fingerprint(
                {"event_id": eid, "comparison_id": cid, "event_type": event_type, "geometry": geometry})
            card = found.setdefault(key, {"event_id": eid, "event_type": event_type, "comparison_id": cid,
                                          "gene_id": target.get("gene_id"), "gene_symbol": target.get("gene_symbol"),
                                          "event_geometry": geometry, "evidence": []})
            card["evidence"].append({"evidence_id": evidence["evidence_id"], "analysis_id": evidence.get("analysis_id", ""),
                                     "source_sha256": evidence.get("source_sha256", evidence["sha256"]),
                                     "parameters": evidence.get("parameters", {}), "metadata": evidence.get("metadata", {}),
                                     "row": portable(row.to_dict())})
    if target["target_kind"] == "event" and not found:
        geometry = target.get("event_geometry", {})
        for cid in target.get("comparison_ids", [""]):
            key = fingerprint([target["event_id"], cid, geometry])
            found[key] = {"event_id": target["event_id"], "event_type": target.get("event_type"),
                          "comparison_id": cid, "gene_id": target.get("gene_id"),
                          "gene_symbol": target.get("gene_symbol"), "event_geometry": geometry, "evidence": []}
    return list(found.values())


def _sashimi_events(cards):
    events = []
    for card in cards:
        # Keep JC/JCEC and conflicting saved statistical families separate;
        # never choose the most significant or first convenient rMATS result.
        candidates = [e for e in card["evidence"] if any(k in e["row"] for k in ("IncLevel1", "IJC_SAMPLE_1"))]
        if not candidates:
            candidates = [{"row": {}, "metadata": {}, "parameters": {}, "evidence_id": "geometry-only"}]
        seen = set()
        for evidence in candidates:
            row = dict(evidence["row"])
            mode = str(evidence.get("metadata", {}).get("count_mode") or
                       evidence.get("parameters", {}).get("count_mode") or "unknown")
            family = fingerprint({"mode": mode, "source": evidence.get("source_sha256"),
                                  "native": {k: row.get(k) for k in ("FDR", "fdr", "PValue", "IncLevelDifference", "dPSI", "dpsi",
                                                                      "IncLevel1", "IncLevel2", "IJC_SAMPLE_1", "IJC_SAMPLE_2", "SJC_SAMPLE_1", "SJC_SAMPLE_2")}})
            if family in seen:
                continue
            seen.add(family)
            row.update(card["event_geometry"])
            row.update(event_id=card["event_id"], event_type=card["event_type"], comparison_id=card["comparison_id"],
                       source_evidence_id=evidence["evidence_id"], count_mode=mode, native_statistical_family=family)
            events.append(row)
    return events


def _target_identity(request, target, resources, method):
    samples = request["project"].get("samples", [])
    cids = set(target.get("comparison_ids", []))
    relevant = {str(s) for c in request["project"].get("contrasts", []) if c["contrast_id"] in cids
                for arm in ("numerator", "denominator") for s in c.get(arm, {}).get("sample_ids", [])}
    ids = {s.get(name) for s in samples if str(s["sample_id"]) in relevant for name in ("bam_resource_id", "bai_resource_id")}
    ids.add(request["project"].get("references", {}).get("gtf_resource_id"))
    bound = {rid: r for rid, r in _resource_public(resources).items() if rid in ids}
    return fingerprint({"project": request["project"], "target": target,
                        "evidence": [e for e in request.get("evidence", []) if target["target_id"] in e.get("target_ids", [])],
                        "modules": request.get("modules", []), "resources": bound, "method": method})


def _target_directory(root, target, identity):
    tid = target["target_id"]
    if not isinstance(tid, str) or not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]{0,127}", tid) or tid in {".", ".."}:
        raise ValueError("Target identity is not a safe portable path component")
    return root / "targets" / tid / identity[:20]


def _read_target_product(root, target, identity):
    directory = _target_directory(root, target, identity)
    candidates = [directory / "target-manifest.json"]
    if directory.parent.exists():
        candidates += sorted(directory.parent.glob(directory.name + "-retry-*/target-manifest.json"), reverse=True)
    available = []
    for candidate in candidates:
        if candidate.exists():
            document = json.loads(candidate.read_text(encoding="utf-8"))
            available.append((candidate, document))
    if not available:
        return None
    path, document = next(((p, d) for p, d in available if d.get("status") in {"complete", "partial"}), available[0])
    directory = path.parent
    if document.get("schema") != TARGET_SCHEMA or document.get("input_identity") != identity:
        raise ValueError("Target product ownership or input identity mismatch")
    verify_ledger(directory, document["members"])
    document["report"] = (directory / "REPORT.html").relative_to(root).as_posix()
    document["manifest"] = path.relative_to(root).as_posix()
    return document


def _target_execute(request, bundle_root, output_root, target, resources, method, progress=None):
    tid = target["target_id"]
    identity = _target_identity(request, target, resources, method)
    directory = _target_directory(output_root, target, identity)
    if progress:
        progress({"type": "target_started", "target_id": tid, "input_identity": identity})
    with owned_claim(output_root / ".claims" / (tid + ".json"), identity):
        existing = _read_target_product(output_root, target, identity)
        if existing and existing.get("status") in {"complete", "partial"}:
            if progress:
                progress({"type": "target_reused", "target_id": tid, "status": existing["status"]})
            return {**existing, "cache_hit": True}
        if directory.exists():
            # Failed attempts are kept; a fresh retry has its own products and
            # cannot erase diagnostics or overwrite a completed product.
            if existing and existing.get("status") == "failed":
                directory = directory.with_name(directory.name + "-retry-" + uuid.uuid4().hex[:8])
            else:
                raise ValueError(f"Refusing to overwrite unverified target products: {tid}")
        directory.mkdir(parents=True)
        record = {"schema": TARGET_SCHEMA, "request_id": request["request_id"],
                  "project_id": request["project"]["project_id"], "target_id": tid,
                  "target_kind": target["target_kind"], "gene_id": target.get("gene_id"),
                  "gene_symbol": target.get("gene_symbol"), "event_id": target.get("event_id"),
                  "comparison_ids": target.get("comparison_ids", []), "input_identity": identity,
                  "created_at": utc_now(), "method": method, "cache_hit": False,
                  "resource_verification": _resource_public(resources)}
        try:
            tables = _tables_for_target(request, bundle_root, target)
            cards = _event_cards(request, target, tables)
            gene_card = {"schema": "katana.gene_card/v1", "project_id": request["project"]["project_id"],
                         "gene_id": target.get("gene_id"), "gene_symbol": target.get("gene_symbol"),
                         "comparison_ids": target.get("comparison_ids", []), "expansion": "all_exported_measured_events",
                         "events": cards, "source_selections": target.get("source_selections", [])}
            write_json(directory / "gene-card.json", gene_card)
            write_json(directory / "event-cards.json", {"schema": "katana.event_cards/v1", "cards": cards})
            flat_cards = [{k: v for k, v in card.items() if k != "evidence"} | {"evidence_count": len(card["evidence"])} for card in cards]
            render_table(directory / "event-cards.html", "Complete exact-event cards", pd.DataFrame(flat_cards),
                         "Source effects and P/FDR are retained in complete evidence tables and event-cards.json.")
            evidence_outputs = []
            (directory / "tables").mkdir()
            for i, (item, frame) in enumerate(tables):
                stem = f"tables/evidence-{i:04d}"
                frame.to_parquet(directory / (stem + ".parquet"), index=True)
                frame.to_csv(directory / (stem + ".tsv"), sep="\t", index=False, float_format="%.17g")
                write_json(directory / (stem + ".json"), {"columns": list(frame.columns), "rows": frame.to_dict("records"),
                                                         "dtypes": {str(c): str(d) for c, d in frame.dtypes.items()}})
                render_table(directory / (stem + ".html"), item.get("name", item["evidence_id"]), frame,
                             "Original native effects, P/FDR, columns and source statistical-family identities; no refit.")
                evidence_outputs.append({"evidence_id": item["evidence_id"], "analysis_id": item.get("analysis_id", ""),
                                         "name": item.get("name", item["evidence_id"]), "rows": len(frame),
                                         "parquet": stem + ".parquet", "tsv": stem + ".tsv", "json": stem + ".json", "html": stem + ".html",
                                         "source_sha256": item.get("source_sha256", item["sha256"]),
                                         "source_fingerprint": item.get("source_fingerprint"), "parameters": item.get("parameters", {}),
                                         "metadata": item.get("metadata", {})})
            measurements = sample_measurements(request, target, tables)
            measurements.to_parquet(directory / "sample-measurements.parquet", index=False)
            measurements.to_csv(directory / "sample-measurements.tsv", sep="\t", index=False, float_format="%.17g")
            render_table(directory / "sample-measurements.html", "All saved biological sample measurements", measurements)
            plots = render_sample_plots(directory, measurements)
            sashimi = [{"status": "not_requested", "reason": "Native Sashimi was not requested", "files": []}]
            if "sashimi" in request.get("modules", []):
                sashimi = generate_sashimi(request, target, _sashimi_events(cards), resources, directory / "sashimi")
            status = ("failed" if any(r["status"] == "failed" for r in sashimi) else
                      "complete" if all(r["status"] in {"complete", "not_requested"} for r in sashimi) else "partial")
            record.update(status=status, event_count=len(cards), evidence=evidence_outputs,
                          samples=len(measurements), plots=plots, sashimi=sashimi)
            provenance = {"schema": "katana.target_provenance/v1", "request_id": request["request_id"],
                          "project": request["project"], "target": target, "method": method,
                          "evidence": evidence_outputs, "resources": _resource_public(resources),
                          "statistical_refit": False, "motif_rescan": False, "network_fit": False,
                          "sashimi": sashimi}
            write_json(directory / "target-provenance.json", provenance)
            write_json(directory / "interactive-data.json", {"schema": "katana.target_interactive/v1", "target": target,
                                                             "gene_card": gene_card, "event_cards": cards,
                                                             "evidence": evidence_outputs, "samples": measurements.to_dict("records"),
                                                             "plots": plots, "sashimi": sashimi})
            render_target(directory, request, target, evidence_outputs, cards, measurements, plots, sashimi, status)
        except Exception as exc:
            record.update(status="failed", reason=f"{type(exc).__name__}: {exc}")
            from html import escape
            (directory / "REPORT.html").write_text(html_page("Katana target failed", "<section><h2>Target diagnostics</h2><p>"
                                                              + escape(record["reason"]) + "</p></section>"), encoding="utf-8")
            write_json(directory / "failure.json", {"target_id": tid, "input_identity": identity, "reason": record["reason"]})
        record["members"] = ledger(directory, exclude=("target-manifest.json",))
        write_json(directory / "target-manifest.json", record)
        record["report"] = (directory / "REPORT.html").relative_to(output_root).as_posix()
        record["manifest"] = (directory / "target-manifest.json").relative_to(output_root).as_posix()
        if progress:
            progress({"type": "target_finished", "target_id": tid, "status": record["status"]})
        return record


def _capacity():
    try:
        return max(1, len(os.sched_getaffinity(0)))
    except (AttributeError, OSError):
        return max(1, os.cpu_count() or 1)


def _plan(request, root, resources, outdir=None):
    methods = _method_identity()
    tasks, errors, warnings = [], [], []
    for target in request["targets"]:
        tid = target["target_id"]
        try:
            tables = _tables_for_target(request, root, target)
            cards = _event_cards(request, target, tables)
            identity = _target_identity(request, target, resources, methods)
            reused = _read_target_product(Path(outdir), target, identity) if outdir else None
            tasks.append({"target_id": tid, "target_kind": target["target_kind"], "comparison_ids": target.get("comparison_ids", []),
                          "event_count": len(cards), "evidence_rows": sum(len(frame) for _, frame in tables),
                          "input_identity": identity, "reusable": bool(reused and reused["status"] in {"complete", "partial"}),
                          "modules": request.get("modules", []), "executor": "current_python_process"})
        except Exception as exc:
            errors.append({"target_id": tid, "message": f"{type(exc).__name__}: {exc}"})
    if "sashimi" in request.get("modules", []):
        for rid, item in resources.items():
            if item["kind"] in {"bam", "bai"} and not item["available"]:
                warnings.append({"resource_id": rid, "message": item["reason"], "affected_module": "sashimi"})
            elif item.get("available") and not item.get("verification", {}).get("identity_verified", False):
                warnings.append({"resource_id": rid, "message": "Execution-machine fingerprint exists; equality to source-machine resource is unverified"})
            elif item.get("available") and not item.get("verification", {}).get("full_content_verified", False):
                warnings.append({"resource_id": rid, "message": "Resource matches its declared sampled fingerprint; this is not full-file verification"})
        if not native_reader():
            warnings.append({"message": "Install bamnostic, pysam or samtools for native Sashimi"})
    return {"schema": "katana.target_plan/v1", "request_id": request["request_id"],
            "project_id": request["project"]["project_id"], "ready": not errors, "errors": errors, "warnings": warnings,
            "tasks": tasks, "resources": _resource_public(resources), "method": methods,
            "statistical_refit": False, "motif_rescan": False, "network_fit": False,
            "executor": "local Python, usable inside any scheduler allocation; no scheduler submission"}


def plan_targets(request_path, resources=None, *, outdir=None):
    from .bundle import open_bundle
    bindings, _ = _load_bindings(resources)
    with open_bundle(request_path) as bundle:
        verified = _bound_resources(bundle.request, bindings)
        return _plan(bundle.request, bundle.root, verified, outdir=outdir)


def _owner(root, request):
    root.mkdir(parents=True, exist_ok=True)
    owner = root / ".target-owner.json"
    expected = {"schema": "katana.target_output_owner/v1", "request_id": request["request_id"],
                "project_id": request["project"]["project_id"], "request_sha256": fingerprint(request)}
    def read_owner():
        # O_EXCL prevents two owners; allow a concurrently opened file to finish
        # its small JSON write without ever replacing an unknown owner.
        import time
        for attempt in range(20):
            try:
                return json.loads(owner.read_text(encoding="utf-8"))
            except json.JSONDecodeError:
                if attempt == 19:
                    raise
                time.sleep(.01)
    if owner.exists():
        if read_owner() != expected:
            raise ValueError("Output directory belongs to a different portable request")
        return
    if any(root.iterdir()):
        if owner.exists() and read_owner() == expected:
            return
        raise ValueError("Refusing an existing nonempty output directory without matching target ownership")
    try:
        fd = os.open(str(owner), os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
        with os.fdopen(fd, "w", encoding="utf-8") as stream:
            json.dump(expected, stream, sort_keys=True)
    except FileExistsError:
        if read_owner() != expected:
            raise ValueError("Concurrent target output owner does not match") from None


def _copy_inputs(request, source, destination):
    request_path = destination / "request.json"
    if request_path.exists():
        if json.loads(request_path.read_text(encoding="utf-8")) != request:
            raise ValueError("Existing copied request was modified")
    else:
        write_json(request_path, request)
    for relative, proof in request.get("members", {}).items():
        if relative == "request.json":
            continue
        target = bounded(destination, relative)
        origin = bounded(source, relative)
        if target.exists():
            if target.stat().st_size != proof["bytes"] or sha256(target) != proof["sha256"]:
                raise ValueError(f"Existing copied source evidence was modified: {relative}")
        else:
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(origin, target)


def _summarize(request, root, latest):
    current = {}
    old_path = root / "result-manifest.json"
    if old_path.exists():
        old = json.loads(old_path.read_text(encoding="utf-8"))
        if old.get("request_id") != request["request_id"]:
            raise ValueError("Foreign result manifest in target output directory")
        for item in old.get("targets", []):
            current[item["target_id"]] = item
    for item in latest:
        current[item["target_id"]] = item
    targets = []
    for target in request["targets"]:
        tid = target["target_id"]
        if tid in current:
            targets.append(current[tid])
        else:
            targets.append({"target_id": tid, "target_kind": target["target_kind"], "gene_id": target.get("gene_id"),
                            "gene_symbol": target.get("gene_symbol"), "event_id": target.get("event_id"), "status": "not_requested"})
    statuses = {t["status"] for t in targets}
    status = "complete" if statuses == {"complete"} else "failed" if statuses <= {"failed", "not_requested"} else "partial"
    return {"schema": RESULT_SCHEMA, "request_id": request["request_id"], "project_id": request["project"]["project_id"],
            "project_name": request["project"].get("name", ""), "status": status,
            "created_at": utc_now(), "report": "INDEX.html", "archive": "REPORT.zip",
            "targets": targets, "source_project_manifest_sha256": request["project"].get("source_manifest_sha256")}


def _archive(root, members):
    archive = root / "REPORT.zip"
    temporary = root / (".archive-" + uuid.uuid4().hex + ".tmp")
    try:
        with zipfile.ZipFile(temporary, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=6) as output:
            for relative in sorted(set(members) | {"result-manifest.json"}):
                output.write(bounded(root, relative), relative)
        os.replace(temporary, archive)
    finally:
        temporary.unlink(missing_ok=True)


def _accepted_batch_members(root, request, manifest):
    """Only accepted products enter a download, never unknown/orphan files."""
    names = {"request.json", "INDEX.html", *request["members"]}
    for target in manifest["targets"]:
        if not target.get("manifest"):
            continue
        receipt_path = bounded(root, target["manifest"])
        receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
        verify_ledger(receipt_path.parent, receipt.get("members"))
        names.add(target["manifest"])
        prefix = receipt_path.parent.relative_to(root)
        names.update((prefix / name).as_posix() for name in receipt["members"])
    result = {}
    for name in sorted(names):
        path = bounded(root, name)
        if path.suffix.lower() in {".bam", ".bai", ".cram", ".crai"}:
            raise ValueError("BAM/CRAM files are execution resources, not report download products")
        result[name] = {"sha256": sha256(path), "bytes": path.stat().st_size}
    return result


def run_targets(request_path, output_dir, resources=None, *, workers=None, target_ids=None, progress=None):
    from .bundle import open_bundle
    root = Path(output_dir).expanduser().resolve()
    bindings, _ = _load_bindings(resources)
    if workers is not None and (not isinstance(workers, int) or isinstance(workers, bool) or workers <= 0):
        raise ValueError("workers must be a positive integer")
    with open_bundle(request_path) as bundle:
        request = bundle.request
        requested_ids = set(target_ids) if target_ids is not None else {t["target_id"] for t in request["targets"]}
        available_ids = {t["target_id"] for t in request["targets"]}
        if not requested_ids or requested_ids - available_ids:
            raise ValueError("target_ids must be a nonempty subset of the portable request")
        resources_checked = _bound_resources(request, bindings)
        method = _method_identity()
        _owner(root, request)
        with owned_claim(root / ".claims" / "inputs.json", request["request_id"], wait_seconds=10):
            _copy_inputs(request, bundle.root, root)
        selected = [t for t in request["targets"] if t["target_id"] in requested_ids]
        latest = []
        count = min(len(selected), workers if workers is not None else _capacity())
        with ThreadPoolExecutor(max_workers=count) as pool:
            pending = {pool.submit(_target_execute, request, bundle.root, root, t,
                                   resources_checked, method, progress): t for t in selected}
            for future in as_completed(pending):
                target = pending[future]
                try:
                    latest.append(future.result())
                except Exception as exc:
                    # A foreign claim, product corruption, or unknown files are
                    # retained and reported. Other independent targets continue.
                    latest.append({"target_id": target["target_id"], "target_kind": target["target_kind"],
                                   "gene_id": target.get("gene_id"), "gene_symbol": target.get("gene_symbol"),
                                   "event_id": target.get("event_id"), "status": "failed",
                                   "reason": f"{type(exc).__name__}: {exc}", "cache_hit": False})
        # Scheduler-per-target callers share the directory safely. Each target
        # has a claim; the finalizer is exclusive and never removes another claim.
        with owned_claim(root / ".claims" / "batch-index.json", request["request_id"], wait_seconds=30):
            manifest = _summarize(request, root, latest)
            old_path = root / "result-manifest.json"
            all_reused = bool(latest) and all(t.get("cache_hit") for t in latest)
            if all_reused and old_path.exists():
                previous = json.loads(old_path.read_text(encoding="utf-8"))
                before = [(t["target_id"], t.get("input_identity"), t["status"]) for t in previous.get("targets", [])]
                after = [(t["target_id"], t.get("input_identity"), t["status"]) for t in manifest["targets"]]
                if before == after:
                    result = inspect_results(root)
                    result.update(status="cache_hit", result_status=previous["status"],
                                  cache={"hits": len(latest), "requested": len(latest)},
                                  requested_target_ids=sorted(requested_ids), execution_status=_execution_status(latest))
                    return result
            render_index(root, manifest)
            members = _accepted_batch_members(root, request, manifest)
            manifest["members"] = members
            atomic_json(old_path, manifest)
            _archive(root, members)
            result = inspect_results(root)
        result["cache"] = {"hits": sum(bool(t.get("cache_hit")) for t in latest), "requested": len(latest)}
        result["requested_target_ids"] = sorted(requested_ids)
        result["execution_status"] = _execution_status(latest)
        return result


def _execution_status(targets):
    statuses = {t["status"] for t in targets}
    return "complete" if statuses == {"complete"} else "failed" if statuses <= {"failed"} else "partial"


def inspect_results(output_dir):
    """Verify portable result members, nested target receipts and the ZIP."""
    path = Path(output_dir).expanduser().resolve()
    root = path.parent if path.is_file() else path
    manifest_path = path if path.is_file() else root / "result-manifest.json"
    if manifest_path.name != "result-manifest.json":
        raise ValueError("Expected result-manifest.json or its containing directory")
    document = json.loads(manifest_path.read_text(encoding="utf-8"))
    if document.get("schema") != RESULT_SCHEMA:
        raise ValueError("Unsupported target results schema")
    verify_ledger(root, document.get("members"))
    for target in document.get("targets", []):
        if not target.get("manifest"):
            continue
        receipt_path = bounded(root, target["manifest"])
        receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
        if (receipt.get("schema") != TARGET_SCHEMA or receipt.get("target_id") != target["target_id"]
                or receipt.get("request_id") != document["request_id"]
                or receipt.get("project_id") != document["project_id"]):
            raise ValueError("Nested target receipt identity mismatch")
        verify_ledger(receipt_path.parent, receipt.get("members"))
        if receipt.get("status") != target.get("status") or receipt.get("input_identity") != target.get("input_identity"):
            raise ValueError("Nested target status or input identity differs from the result manifest")
    statuses = {t["status"] for t in document.get("targets", [])}
    expected_status = "complete" if statuses == {"complete"} else "failed" if statuses <= {"failed", "not_requested"} else "partial"
    if document.get("status") != expected_status:
        raise ValueError("Batch result status differs from its actual target statuses")
    archive = bounded(root, document.get("archive", "REPORT.zip"))
    if archive.exists():
        with zipfile.ZipFile(archive) as stream:
            infos = stream.infolist()
            names = [i.filename for i in infos]
            if len(names) != len(set(names)):
                raise ValueError("Duplicate ZIP result members")
            expected = set(document["members"]) | {"result-manifest.json"}
            if set(names) != expected:
                raise ValueError("Result ZIP member set differs from the verified output")
            for info in infos:
                bounded(root, info.filename)
                if info.file_size > 0 and ((info.external_attr >> 16) & 0o170000) == 0o120000:
                    raise ValueError("Symlink ZIP result member")
                value = stream.read(info)
                if info.filename == "result-manifest.json":
                    if json.loads(value) != document:
                        raise ValueError("ZIP manifest differs from the accepted result manifest")
                elif (len(value) != document["members"][info.filename]["bytes"]
                      or hashlib.sha256(value).hexdigest() != document["members"][info.filename]["sha256"]):
                    raise ValueError(f"ZIP result member hash mismatch: {info.filename}")
            bad = stream.testzip()
            if bad:
                raise ValueError(f"ZIP CRC failure: {bad}")
    return {**document, "verified": True, "report": str(bounded(root, document.get("report", "INDEX.html"))),
            "archive": str(archive) if archive.exists() else None, "manifest": str(manifest_path),
            "archive_status": "verified" if archive.exists() else "unavailable"}
