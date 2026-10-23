"""Explicit transport of historical results, without scientific recomputation.

An adoption receipt proves an unchanged, bounded snapshot and its study binding.
It does not promote a file inventory to a native execution or numerical proof.
Original completion receipts and family manifests travel unchanged; the final
acceptance layer remains responsible for their method-specific interpretation.
"""
from __future__ import annotations

import copy
from datetime import datetime, timezone
from pathlib import Path
import shutil

from .runtime import canonical_hash, file_record, read_json, verify_files, write_json


SOURCE_SCHEMA = "katana.flow.historical_source/v1"
ADOPTION_SCHEMA = "katana.flow.historical_adoption/v1"


def _signed(value):
    result = dict(value)
    result["receipt_sha256"] = canonical_hash(result)
    return result


def _authenticate(value, schema):
    if (value.get("schema") != schema or
        value.get("receipt_sha256") != canonical_hash({k: v for k, v in value.items() if k != "receipt_sha256"})):
        raise ValueError("Historical snapshot content identity is invalid")


def _bounded(path, roots):
    path = Path(path).expanduser()
    if path.is_symlink():
        raise ValueError("Historical products must be regular files, not symlinks")
    path = path.resolve()
    owners = [root for root in roots if path.is_relative_to(root)]
    if not owners:
        raise ValueError(f"Historical result is outside the explicitly allowed roots: {path}")
    return path, max(owners, key=lambda root: len(root.parts))


def _path_members(record, *, include_exports=True):
    """Only documented result paths; labels and sequence strings stay opaque."""
    if record.get("path"):
        yield record["path"]
    meta = record.get("metadata", {})
    if include_exports:
        for value in meta.get("exports", {}).values():
            if isinstance(value, str):
                yield value
    for item in [*record.get("data_links", []), *meta.get("data_links", [])]:
        if item.get("path"):
            yield item["path"]


def _family_members(root):
    root = Path(root).resolve()
    marker = root / "COMPLETE.json"
    value = read_json(marker)
    members = value.get("members")
    if value.get("complete") is not True or not isinstance(members, dict) or not members:
        raise ValueError("Completed scientific family requires its full native member ledger")
    yield marker
    for name, expected in members.items():
        relative = Path(name)
        candidate = (root / relative).resolve()
        if relative.is_absolute() or not candidate.is_relative_to(root) or candidate == root:
            raise ValueError("Scientific family member escapes its completed directory")
        actual = file_record(candidate)
        wanted = expected if isinstance(expected, str) else expected.get("sha256")
        size = expected.get("bytes", expected.get("size_bytes", actual["size_bytes"])) if isinstance(expected, dict) else actual["size_bytes"]
        if wanted != actual["sha256"] or size != actual["size_bytes"]:
            raise ValueError("Native completed scientific family member changed")
        yield candidate


def _raw_artifacts(raw, descriptor, project):
    """Preserve complete catalog versions even when the report deduplicates."""
    records = raw.get("artifacts")
    if not isinstance(records, list):
        name = raw.get("artifact_catalog") or raw.get("catalog_export") or raw.get("native_artifact_catalog")
        catalog = Path(name) if name else project.parent / "native-artifact-catalog.json"
        catalog = catalog if catalog.is_absolute() else descriptor.parent / catalog
        import json
        export = json.loads(catalog.read_text(encoding="utf-8-sig"))
        records = export.get("artifacts") if isinstance(export, dict) else export
    if not isinstance(records, list):
        raise ValueError("Historical import requires an explicit artifact record list")
    result = []
    for item in records:
        if item.get("status", "complete") not in {"complete", "completed"}:
            continue
        row = copy.deepcopy(item)
        path = Path(row["path"])
        row["path"] = str((path if path.is_absolute() else project.parent / path).resolve())
        metadata = row.get("metadata", {})
        for key, value in metadata.get("exports", {}).items():
            if isinstance(value, str):
                path = Path(value)
                metadata["exports"][key] = str((path if path.is_absolute() else project.parent / path).resolve())
        result.append(row)
    return result


def seal_existing_results(descriptor_path, output_path, *, allowed_roots,
                          completion_evidence=(), family_roots=(), include_files=(),
                          purpose="report-results", include_exports=True):
    """Seal a read-only report/catalog snapshot with explicit local file bounds.

    The source project and original contrast definitions remain authoritative.
    All original completed artifact records are retained, including comparisons
    outside the report-only selection. ``completion_evidence`` is an explicit
    list of old receipt paths, not an assertion that every file is a valid fit.
    Additional parser-required tables can be supplied through ``include_files``.
    """
    from katana_next.reports.sources import load_source
    if purpose not in {"report-results", "analysis-inputs"} or (not include_exports and purpose != "analysis-inputs"):
        raise ValueError("Excluding human-readable exports requires an explicit analysis-inputs snapshot")
    descriptor_path = Path(descriptor_path).expanduser().resolve()
    output_path = Path(output_path).expanduser().resolve()
    if output_path.exists():
        raise FileExistsError("Existing historical snapshot must not be overwritten")
    roots = [Path(root).expanduser().resolve() for root in allowed_roots]
    if not roots or len(roots) != len(set(roots)) or any(not root.is_dir() for root in roots):
        raise ValueError("Historical import requires distinct existing allowed result roots")
    descriptor_identity = file_record(descriptor_path)
    raw = read_json(descriptor_path)
    if not raw.get("project_manifest"):
        raise ValueError("Historical import requires an explicit original project manifest")
    project_path = Path(raw["project_manifest"])
    project_path = (project_path if project_path.is_absolute() else descriptor_path.parent / project_path).resolve()
    project_identity = file_record(project_path)
    native = read_json(project_path)
    # Avoid reading a mutable database as if it were a sealed artifact export.
    if not (isinstance(raw.get("artifacts"), list) or raw.get("artifact_catalog") or
            raw.get("catalog_export") or raw.get("native_artifact_catalog") or
            (project_path.parent / "native-artifact-catalog.json").is_file()):
        raise ValueError("Export the historical artifact catalog before adopting results")
    loaded = load_source(descriptor_path)
    project_id = native.get("project_id")
    if not project_id or loaded["native_project_manifest_sha256"] != file_record(project_path)["sha256"]:
        raise ValueError("Historical source project identity is inconsistent")
    if loaded.get("supplemental"):
        raise ValueError("Seal supplemental studies separately with their explicit parent bindings")
    artifacts = loaded["artifacts"]
    all_artifacts = _raw_artifacts(raw, descriptor_path, project_path)
    if not artifacts:
        raise ValueError("Historical source has no actual completed artifacts")
    ids = [row["artifact_id"] for row in artifacts]
    if len(ids) != len(set(ids)):
        raise ValueError("Historical artifact IDs must be unique")
    known_contrasts = {row["contrast_id"] for row in native.get("contrasts", [])}
    paths = {descriptor_path, project_path, *map(Path, loaded["source_files"]), *map(Path, include_files), *map(Path, completion_evidence)}
    for record in [*all_artifacts, *artifacts, *loaded["extra_assets"]]:
        if record.get("project_id") != project_id or not set(record.get("contrast_ids", [])) <= known_contrasts:
            raise ValueError("Historical artifact belongs to another study or comparison")
        paths.update(map(Path, _path_members(record, include_exports=include_exports)))
    for view in loaded["extra_views"]:
        for section, keys in (("panels", ("path", "image", "svg")), ("tables", ("path", "source", "view")),
                              ("top10", ("full_path", "full_table"))):
            for item in view.get(section, []):
                paths.update(Path(item[key]) for key in keys if item.get(key))
                paths.update(Path(item["path"]) for item in item.get("data_links", []) if item.get("path"))
    for root in family_roots:
        paths.update(_family_members(root))
    # Keep producer completion markers and browser attachments as actual files,
    # without editing their old JSON or inferring a successful native fit.
    by_id = {row["artifact_id"]: row for row in [*all_artifacts, *loaded["extra_assets"]] if row.get("artifact_id")}
    from katana_next.reports.rich_results import _dependencies
    for path in list(paths):
        path, owner = _bounded(path, roots)
        for parent in (path.parent, *path.parents):
            if not parent.is_relative_to(owner):
                break
            marker = parent / "parameters.json"
            if marker.is_file():
                paths.add(marker)
                break
        if path.suffix.lower() in {".html", ".htm", ".js", ".css"}:
            closure = _dependencies(path, owner, by_id)
            paths.update(Path(row["path"]) for row in closure["local_dependencies"])
            if closure["portability_warnings"]:
                raise ValueError("Historical browser has unresolved local attachments: " + "; ".join(closure["portability_warnings"]))
    ledger = []
    for candidate in sorted(paths, key=str):
        path, owner = _bounded(candidate, roots)
        record = file_record(path)
        ledger.append({**record, "root_index": roots.index(owner), "relative_path": path.relative_to(owner).as_posix()})
    # Keep the view structure and all seven chapters. Only explicit paths are
    # resolved; no unavailable branch is removed and no statistics are changed.
    report = {"project_manifest": str(project_path), "artifacts": artifacts,
              "extra_assets": loaded["extra_assets"], "extra_views": loaded["extra_views"],
              "report_contrast_ids": loaded["report_selection"]["contrast_ids"],
              "report_contrast_groups": loaded["report_selection"]["groups"],
              "mode": loaded["mode"], "dummy": loaded["dummy"],
              "report_native_figures_in_main": loaded["report_native_figures_in_main"],
              "report_require_plot_data": loaded["report_require_plot_data"]}
    result = _signed({"schema": SOURCE_SCHEMA, "status": "sealed", "project_id": project_id,
        "source_project_manifest": project_identity, "source_descriptor": descriptor_identity,
        "allowed_roots": list(map(str, roots)), "files": ledger, "report_descriptor": report,
        "original_artifact_records": all_artifacts,
        "purpose": purpose, "include_exports": bool(include_exports),
        "excluded_export_scope": "None" if include_exports else "Original export pointers retained as provenance; only analysis input files are transported",
        "completion_evidence": [file_record(path) for path in completion_evidence],
        "family_roots": list(map(lambda root: str(Path(root).resolve()), family_roots)),
        "scientific_computation_performed": False,
        "scope": "Unchanged historical file snapshot; native method-specific completion is not inferred from inventory",
        "sealed_at": datetime.now(timezone.utc).isoformat()})
    verify_files([*result["files"], project_identity, descriptor_identity])
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("x", encoding="utf-8") as stream:
        import json
        json.dump(result, stream, indent=2, ensure_ascii=False)
    return output_path


def validate_existing_source(path):
    value = read_json(path)
    _authenticate(value, SOURCE_SCHEMA)
    if value.get("status") != "sealed" or value.get("scientific_computation_performed") is not False:
        raise ValueError("Historical import requires an explicitly sealed transport-only source")
    roots = list(map(Path, value["allowed_roots"]))
    for row in value["files"]:
        actual, owner = _bounded(row["path"], roots)
        if (row["root_index"] != roots.index(owner) or row["relative_path"] != actual.relative_to(owner).as_posix()):
            raise ValueError("Historical source file ledger changed its bounded ownership")
    verify_files(value["files"])
    verify_files([value["source_project_manifest"], value["source_descriptor"], *value["completion_evidence"]])
    for root in value["family_roots"]:
        list(_family_members(root))
    return value


def _rewrite(value, mapping):
    if isinstance(value, dict):
        return {key: _rewrite(child, mapping) for key, child in value.items()}
    if isinstance(value, list):
        return [_rewrite(child, mapping) for child in value]
    return mapping.get(value, value) if isinstance(value, str) else value


def adopt_existing_results(source_path, output_dir):
    """Copy an already sealed historical snapshot once, without any job calls."""
    source_path = Path(source_path).resolve()
    source = validate_existing_source(source_path)
    output = Path(output_dir).expanduser().resolve()
    receipt_path = output / "adoption.receipt.json"
    if receipt_path.exists():
        receipt = validate_adoption(output)
        if receipt["source_snapshot_sha256"] != file_record(source_path)["sha256"]:
            raise ValueError("Existing adoption belongs to another historical snapshot")
        return receipt_path
    if output.exists():
        raise FileExistsError("Preserve the unknown adoption attempt and choose a new output directory")
    if output == source_path.parent:
        raise ValueError("Historical snapshot cannot replace its sealed source directory")
    for root in map(Path, source["allowed_roots"]):
        if output == root or output.is_relative_to(root):
            raise ValueError("Historical snapshot output must be outside original source roots")
    output.mkdir(parents=True)
    write_json(output / "adoption.claim.json", {"source_snapshot": file_record(source_path)})
    mapping, products = {}, []
    for item in source["files"]:
        relative = Path("products") / str(item["root_index"]) / item["relative_path"]
        destination = output / relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(item["path"], destination)
        actual = file_record(destination)
        if (actual["sha256"], actual["size_bytes"]) != (item["sha256"], item["size_bytes"]):
            raise ValueError("Historical product changed while being adopted")
        mapping[item["path"]] = relative.as_posix()
        products.append({**actual, "path": relative.as_posix(), "source": item})
    descriptor = _rewrite(copy.deepcopy(source["report_descriptor"]), mapping)
    descriptor["historical_adoption"] = {"schema": ADOPTION_SCHEMA, "scientific_computation_performed": False,
                                         "source_project_manifest_sha256": source["source_project_manifest"]["sha256"]}
    descriptor_path = write_json(output / "report-project.json", descriptor)
    snapshot_path = output / "historical-source.json"
    shutil.copyfile(source_path, snapshot_path)
    for path in (descriptor_path, snapshot_path):
        products.append({**file_record(path), "path": path.relative_to(output).as_posix()})
    receipt = _signed({"schema": ADOPTION_SCHEMA, "status": "adopted", "project_id": source["project_id"],
        "source_snapshot_sha256": file_record(source_path)["sha256"],
        "source_project_manifest_sha256": source["source_project_manifest"]["sha256"],
        "outputs": products, "report_descriptor": descriptor_path.name,
        "scientific_computation_performed": False, "native_jobs_submitted": 0,
        "scope": "Transport and integrity verification only; original native evidence is preserved unchanged",
        "finished_at": datetime.now(timezone.utc).isoformat()})
    validate_existing_source(source_path)
    write_json(receipt_path, receipt)
    validate_adoption(output)
    return receipt_path


def validate_adoption(path):
    location = Path(path).resolve()
    receipt_path = location / "adoption.receipt.json" if location.is_dir() else location
    root = receipt_path.parent
    value = read_json(receipt_path)
    _authenticate(value, ADOPTION_SCHEMA)
    if value.get("status") != "adopted" or value.get("scientific_computation_performed") is not False or value.get("native_jobs_submitted") != 0:
        raise ValueError("Historical adoption cannot claim a new native scientific run")
    for row in value["outputs"]:
        relative = Path(row["path"])
        candidate = (root / relative).resolve()
        if relative.is_absolute() or not candidate.is_relative_to(root):
            raise ValueError("Historical adoption output escapes its owned directory")
        verify_files([{**row, "path": str(candidate)}])
    source = read_json(root / "historical-source.json")
    _authenticate(source, SOURCE_SCHEMA)
    if (file_record(root / "historical-source.json")["sha256"] != value["source_snapshot_sha256"] or
            source["source_project_manifest"]["sha256"] != value["source_project_manifest_sha256"] or
            source["project_id"] != value["project_id"]):
        raise ValueError("Adopted source and project identities differ")
    return value


def materialize_adopted_descriptor(path):
    root = Path(path).resolve()
    root = root if root.is_dir() else root.parent
    receipt = validate_adoption(root)
    mapping = {row["path"]: str((root / row["path"]).resolve()) for row in receipt["outputs"]}
    return _rewrite(read_json(root / receipt["report_descriptor"]), mapping)


def _study_identity(project, contrast_ids):
    contrasts = {row["contrast_id"]: row for row in project.get("contrasts", [])}
    if len(contrasts) != len(project.get("contrasts", [])) or not set(contrast_ids) <= set(contrasts):
        raise ValueError("Identity projection includes duplicate or unknown comparisons")
    return {"project_id": project["project_id"], "genome_build": project["annotation"]["genome_build"],
        "samples": [{key: row.get(key, {} if key == "attributes" else None)
                     for key in ("sample_id", "biological_sample_id", "condition", "batch", "attributes")}
                    for row in project.get("samples", [])],
        "contrasts": [contrasts[identity] for identity in contrast_ids]}


def create_identity_projection(adoption, compute_project, output_path, *, contrast_ids):
    """Certify study/design equality while recording both distinct file SHAs.

    This authorizes asset/path remapping only. It is not an annotation content
    certificate; actual reference and database file identities are bound in each
    computational task's own input ledger.
    """
    output_path = Path(output_path).resolve()
    if output_path.exists():
        raise FileExistsError("Identity projection certificates are immutable")
    descriptor = materialize_adopted_descriptor(adoption)
    receipt = validate_adoption(adoption)
    original = read_json(descriptor["project_manifest"])
    compute_project = Path(compute_project).resolve()
    overlay = read_json(compute_project)
    contrasts = list(contrast_ids)
    if not contrasts or len(contrasts) != len(set(contrasts)):
        raise ValueError("Identity projection requires unique selected comparisons")
    identity = _study_identity(original, contrasts)
    if identity != _study_identity(overlay, contrasts):
        raise ValueError("Compute overlay changes original biological samples, reference build or contrast arms")
    certificate = _signed({"schema": "katana.flow.study_identity_projection/v1", "status": "verified",
        "project_id": original["project_id"], "contrast_ids": contrasts,
        "historical_adoption_receipt_sha256": receipt["receipt_sha256"],
        "source_project_manifest_sha256": receipt["source_project_manifest_sha256"],
        "compute_project_manifest_sha256": file_record(compute_project)["sha256"],
        "study_identity": identity, "study_identity_sha256": canonical_hash(identity),
        "scope": "Study/design equality only; remapped resource bytes need separate task input verification"})
    write_json(output_path, certificate)
    return output_path


def validate_identity_projection(adoption, identity_projection, compute_project):
    """Read-only study-binding check, suitable for scheduling preflight."""
    receipt = validate_adoption(adoption)
    descriptor = materialize_adopted_descriptor(adoption)
    projection = read_json(identity_projection)
    _authenticate(projection, "katana.flow.study_identity_projection/v1")
    original = read_json(descriptor["project_manifest"])
    compute = read_json(compute_project)
    identity = _study_identity(compute, projection["contrast_ids"])
    if (projection.get("status") != "verified" or projection.get("project_id") != compute.get("project_id") or
        projection.get("historical_adoption_receipt_sha256") != receipt["receipt_sha256"] or
        projection.get("source_project_manifest_sha256") != receipt["source_project_manifest_sha256"] or
        projection.get("compute_project_manifest_sha256") != file_record(compute_project)["sha256"] or
        projection.get("study_identity_sha256") != canonical_hash(identity) or
        projection.get("study_identity") != identity or _study_identity(original, projection["contrast_ids"]) != identity):
        raise ValueError("Historical input identity projection does not match the actual compute project")
    return projection


def import_adopted_artifacts(adoption, identity_projection, compute_project, output_dir, *, artifact_ids=None):
    """Create task-owned historical inputs for later scientific computations.

    Returned artifacts retain their original scientific identities. Only the
    import operation runs now; no former native job or statistical fit is
    manufactured as part of the current flow lifecycle.
    """
    from katana_next.domain.models import AnalysisArtifact
    adoption = Path(adoption).resolve()
    adoption_root = adoption if adoption.is_dir() else adoption.parent
    projection = validate_identity_projection(adoption, identity_projection, compute_project)
    receipt = validate_adoption(adoption)
    descriptor = materialize_adopted_descriptor(adoption)
    source_snapshot = read_json(adoption_root / "historical-source.json")
    originals = source_snapshot["original_artifact_records"]
    known = {row["artifact_id"]: row for row in originals}
    if len(known) != len(originals):
        raise ValueError("Historical catalog contains conflicting original artifact identities")
    chosen = list(artifact_ids) if artifact_ids is not None else [row["artifact_id"] for row in descriptor["artifacts"]]
    if len(chosen) != len(set(chosen)) or not chosen or not set(chosen) <= set(known):
        raise ValueError("Adoption task requires distinct known historical artifact IDs")
    copied_by_source = {row["source"]["path"]: adoption_root / row["path"]
                        for row in receipt["outputs"] if row.get("source")}
    sources = set()
    markers = set()
    for identity_id in chosen:
        row = known[identity_id]
        if row.get("status", "complete") not in {"complete", "completed"}:
            raise ValueError("Historical import cannot turn an incomplete artifact into a completed result")
        sources.update(_path_members(row, include_exports=source_snapshot.get("include_exports", True)))
        for source_name, local in copied_by_source.items():
            if local.name != "parameters.json":
                continue
            marker = read_json(local)
            if (marker.get("completed") is True and marker.get("workflow_id") == row["analysis_id"] and
                marker.get("run_fingerprint") == row["run_fingerprint"] and
                identity_id in {item["artifact_id"] for item in marker.get("generated_artifacts", [])}):
                markers.add(source_name)
                break
        else:
            if row.get("metadata", {}).get("artifact_contract_version") == "3.0-correctness":
                raise ValueError("Historical scientific artifact has no authentic completed producer marker")
    output = Path(output_dir).resolve()
    if output.exists():
        raise FileExistsError("Preserve a previous historical task transport attempt")
    if output.is_relative_to(adoption_root):
        raise ValueError("Historical task cannot write inside its adopted source")
    output.mkdir(parents=True)
    mapping, outputs = {}, []
    for source_name in sorted(sources | markers):
        if source_name not in copied_by_source:
            raise ValueError("Historical scientific product is absent from the verified adoption ledger")
        local = copied_by_source[source_name]
        relative = local.relative_to(adoption_root)
        target = output / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(local, target)
        old = file_record(local)
        current = file_record(target)
        if (current["sha256"], current["size_bytes"]) != (old["sha256"], old["size_bytes"]):
            raise ValueError("Historical scientific product changed during task transport")
        mapping[source_name] = str(target)
        outputs.append(current)
    artifacts = []
    for identity_id in chosen:
        row = _rewrite(copy.deepcopy(known[identity_id]), mapping)
        row.setdefault("metadata", {})["historical_import"] = {
            "source_adoption_receipt_sha256": receipt["receipt_sha256"],
            "source_project_manifest_sha256": receipt["source_project_manifest_sha256"],
            "compute_project_manifest_sha256": projection["compute_project_manifest_sha256"],
            "scientific_computation_performed": False}
        artifacts.append(AnalysisArtifact.model_validate(row).model_dump(mode="json"))
    preserved = output / "evidence"
    preserved.mkdir()
    for original_path in (adoption_root / "adoption.receipt.json", adoption_root / "historical-source.json", Path(identity_projection)):
        target = preserved / original_path.name
        shutil.copyfile(original_path, target)
        outputs.append(file_record(target))
    validate_adoption(adoption)
    return {"artifacts": artifacts, "outputs": outputs,
        "producer_markers": [file_record(mapping[name]) for name in sorted(markers)],
        "inputs": [file_record(adoption_root / "adoption.receipt.json"), file_record(identity_projection)],
        "evidence": {"kind": "historical_import", "scientific_computation_performed": False,
            "native_jobs_submitted": 0, "source_adoption_receipt_sha256": receipt["receipt_sha256"],
            "source_project_manifest_sha256": receipt["source_project_manifest_sha256"],
            "compute_project_manifest_sha256": projection["compute_project_manifest_sha256"],
            "identity_projection_sha256": file_record(identity_projection)["sha256"],
            "imported_artifact_ids": chosen,
            "scope": "Existing scientific outputs copied unchanged; current task completes only transport and validation"}}


__all__ = ["seal_existing_results", "validate_existing_source", "adopt_existing_results",
           "validate_adoption", "materialize_adopted_descriptor", "create_identity_projection",
           "validate_identity_projection", "import_adopted_artifacts"]
