"""Collect independently completed task products for the fixed report renderer."""
from __future__ import annotations
from .model import receipt_plan_id

import copy
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import re
import shutil

from .runtime import (canonical_hash, file_hash, file_record, read_json,
                      validate_receipt, verify_files, write_json, _source_manifest)


COLLECTION_SCHEMA = "katana.flow.collection_receipt/v1"


def _safe(value):
    return re.sub(r"[^A-Za-z0-9._-]", "_", str(value))[:140] or "result"


def _rewrite(value, mapping):
    if isinstance(value, dict):
        return {key: _rewrite(child, mapping) for key, child in value.items()}
    if isinstance(value, list):
        return [_rewrite(child, mapping) for child in value]
    if isinstance(value, str):
        return mapping.get(value, value)
    return value


def validate_collection(path):
    path = Path(path)
    if path.is_dir():
        path /= "collection.receipt.json"
    value = read_json(path)
    if value.get("schema") != COLLECTION_SCHEMA or value.get("status") != "completed":
        raise ValueError("No completed flow collection is available")
    identity = {key: child for key, child in value.items() if key != "receipt_sha256"}
    if value.get("receipt_sha256") != canonical_hash(identity):
        raise ValueError("Flow collection receipt content identity is invalid")
    root = path.parent.resolve()
    ledger = []
    for record in value["outputs"]:
        relative = Path(record.get("relative_path", record["path"]))
        if relative.is_absolute() or ".." in relative.parts:
            raise ValueError("Collection output ledger must contain bounded relative paths")
        location = (root / relative).resolve()
        if not location.is_relative_to(root):
            raise ValueError("Collection output escapes its result directory")
        ledger.append({**record, "path": str(location)})
    verify_files(ledger)
    from .model import load_plan
    plan = load_plan(root / "plan.json")
    if plan["plan_id"] != value.get("plan_id"):
        raise ValueError("Collection and preserved execution plan identities differ")
    planned = {task["id"]: task for task in plan["tasks"]}
    recorded = {task["task_id"]: task for task in value.get("tasks", [])}
    proofs = {task["task_id"]: task for task in value.get("task_receipts", [])}
    if set(planned) != set(recorded) or set(planned) != set(proofs):
        raise ValueError("Collection does not preserve every planned task's completion proof")
    for task_id, spec in planned.items():
        relative = Path(proofs[task_id]["collection_path"])
        if relative.is_absolute() or ".." in relative.parts:
            raise ValueError("Task completion proof must stay inside the collection")
        proof = read_json(root / relative)
        unsigned = {key: item for key, item in proof.items() if key != "receipt_sha256"}
        if (proof.get("receipt_sha256") != canonical_hash(unsigned) or
            proof.get("receipt_sha256") != recorded[task_id].get("receipt_sha256") or
            proof.get("task_contract_sha256") != canonical_hash(spec) or
            proof.get("plan_id") != receipt_plan_id(plan, task_id) or
            proof.get("status") != recorded[task_id].get("status") or
            proof.get("status") not in {"completed", "reused", "not_applicable"} or
            proof.get("source_project_manifest_sha256") != value.get("source_project_manifest_sha256")):
            raise ValueError("Preserved task completion proof does not match the collected plan and products")
    return value


def materialize_report_descriptor(collection):
    """Resolve a portable collection at its present location, without editing it.

    Main artifact paths and static exports share the same verified relative
    ledger.  This keeps the existing report renderer's absolute-path interface
    independent of where the collection was downloaded or mounted.
    """
    root = Path(collection).expanduser().resolve()
    if root.is_file():
        root = root.parent
    receipt = validate_collection(root)
    descriptor = read_json(root / receipt.get("report_descriptor", "report-project.json"))
    catalog = read_json(root / receipt.get("artifact_catalog", "artifacts.json"))
    mapping = {item.get("relative_path", item["path"]): str((root / item.get("relative_path", item["path"])).resolve())
               for item in receipt["outputs"]}
    for item in receipt["outputs"]:
        for parent in Path(item.get("relative_path", item["path"])).parents:
            mapping[parent.as_posix()] = str((root / parent).resolve())
    descriptor = _rewrite(descriptor, mapping)
    records = catalog.get("artifacts", [])
    descriptor["artifacts"] = _rewrite(records, mapping)
    descriptor.pop("artifact_catalog", None)
    descriptor["project_manifest"] = str(root / receipt.get("project_manifest", "project.json"))
    return descriptor


def collect_results(plan_path, task_results, output_dir, *, project_path=None):
    """Copy verified completed products, keeping native tables and statistics.

    The report descriptor uses one read-only JSON catalogue; collection never
    opens or mutates the source project's database. Raw sequencing resources
    remain in the declared project and are not silently copied into a report.
    """
    plan_location = Path(plan_path).expanduser().resolve()
    from .model import load_plan
    plan = load_plan(plan_location)
    source_path = Path(project_path or plan["project_manifest"]).expanduser()
    source_path = (source_path if source_path.is_absolute() else plan_location.parent / source_path).resolve()
    source_sha = file_hash(source_path)
    if plan.get("project_manifest_sha256", plan.get("source_project_manifest_sha256")) != source_sha:
        raise ValueError("Source project changed after planning")
    tasks = {task["id"]: task for task in plan["tasks"]}
    receipts, roots = {}, {}
    for candidate in task_results:
        candidate = Path(candidate).expanduser().resolve()
        if candidate.is_dir():
            candidate /= "task.receipt.json"
        raw = validate_receipt(candidate, plan=plan)
        task_id = raw["task_id"]
        if task_id not in tasks or task_id in receipts:
            raise ValueError("Collection requires exactly one receipt per planned task")
        validate_receipt(candidate, plan=plan, task=tasks[task_id])
        if raw["source_project_manifest_sha256"] != source_sha:
            raise ValueError("Task products belong to a different source project revision")
        receipts[task_id] = (candidate, raw)
        roots[task_id] = candidate.parent
    if set(receipts) != set(tasks):
        raise ValueError(f"Planned tasks have no completed products: {sorted(set(tasks) - set(receipts))}")
    directory = Path(output_dir).expanduser().resolve()
    existing = directory / "collection.receipt.json"
    task_identity = {key: raw["receipt_sha256"] for key, (_path, raw) in sorted(receipts.items())}
    request_identity = canonical_hash({"plan_id": plan["plan_id"], "source": source_sha, "tasks": task_identity})
    if existing.exists():
        value = validate_collection(existing)
        if value.get("request_sha256") != request_identity:
            raise ValueError("Existing collection belongs to different completed tasks")
        return existing
    source = _source_manifest(source_path)
    for protected in (Path(__file__).resolve().parents[2], source.input_root, source.output_root, *roots.values()):
        if directory == protected.resolve() or directory.is_relative_to(protected.resolve()):
            raise ValueError("Collection must be separate from original inputs, source and task products")
    directory.mkdir(parents=True, exist_ok=True)
    if any(directory.iterdir()):
        raise FileExistsError("Collection directory contains an unknown attempt; preserve it and choose a new output directory")
    with (directory / "collection.claim.json").open("x", encoding="utf-8") as stream:
        json.dump({"request_sha256": request_identity}, stream)

    original_files = {}
    artifacts = {}
    generated = set()
    producers = {}
    for task_id, (_receipt_path, raw) in receipts.items():
        generated.update(a["artifact_id"] for a in raw["generated_artifacts"])
        producers.update({a["artifact_id"]: (tasks[task_id], raw) for a in raw["generated_artifacts"]})
        for item in raw["outputs"]:
            original_files[item["path"]] = item
        for artifact in raw["artifacts"]:
            artifact = copy.deepcopy(artifact)
            artifact.get("metadata", {}).pop("project_snapshot_fingerprint", None)
            artifact.get("metadata", {}).pop("flow_snapshot_binding", None)
            previous = artifacts.get(artifact["artifact_id"])
            if previous is not None and previous != artifact:
                raise ValueError("Collection contains conflicting versions of one artifact")
            artifacts[artifact["artifact_id"]] = artifact
    # Auxiliary joins are private execution views; the public collection keeps
    # the original scientific products and their exact comparison granularity.
    artifacts = {key: value for key, value in artifacts.items() if key in generated}
    if not artifacts:
        raise ValueError("Collection has no actual generated scientific products")
    mapping, copied, provenance = {}, [], []
    for source_name, item in sorted(original_files.items()):
        path = Path(source_name).resolve()
        owners = [(key, root) for key, root in roots.items() if path.is_relative_to(root)]
        if not owners:
            raise ValueError("Task output is outside every declared owned task result directory")
        owner, root = max(owners, key=lambda pair: len(pair[1].parts))
        relative = path.relative_to(root)
        if any(part in {"..", ".katana"} for part in relative.parts):
            raise ValueError("Internal catalogue files are not report products")
        destination = directory / "products" / _safe(owner) / relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        from .file_io import retain_owned_file
        retain_owned_file(path, destination)
        if file_hash(destination) != item["sha256"]:
            raise ValueError("Copied scientific product differs from its completed task")
        mapping[source_name] = destination.relative_to(directory).as_posix()
        parent = path.parent
        while parent.is_relative_to(root):
            mapping[str(parent)] = (Path("products") / _safe(owner) / parent.relative_to(root)).as_posix()
            if parent == root:
                break
            parent = parent.parent
        copied.append(destination)
        provenance.append({"source": item, "collection_path": destination.relative_to(directory).as_posix()})
    # Scientific JSON files can themselves be members of COMPLETE.json. Keep
    # every copied producer byte unchanged, including its native provenance.
    # Only the separate collection catalogue is rebased for the renderer.
    exported = []
    for artifact in artifacts.values():
        updated = _rewrite(artifact, mapping)
        producer, completion = producers[artifact["artifact_id"]]
        selection = completion.get("selection") or {}
        scopes = selection.get("contrast_ids", [])
        if producer.get("selection") and scopes:
            if updated.get("contrast_ids") and set(updated["contrast_ids"]) != set(scopes):
                raise ValueError("Generated artifact comparison labels differ from its validated selected source")
            updated["contrast_ids"] = list(scopes)
            updated.setdefault("metadata", {})["flow_contrast_binding"] = {
                "source": "validated saved selection", "contrast_ids": list(scopes),
                "source_artifact_id": selection["source_artifact_id"], "statistics_recomputed": False}
        if updated["path"] == artifact["path"]:
            raise ValueError("A completed scientific artifact is missing from the verified output ledger")
        updated.setdefault("metadata", {})["flow_collection"] = {
            "source_artifact_fingerprint": artifact["fingerprint"],
            "source_project_manifest_sha256": source_sha,
            "statistics_recomputed": False,
        }
        if (updated.get("analysis_id") == "programs.cosplice" and
            updated["metadata"].get("flow_stage") in {"network.fit", "network.associate", "network.annotate"} and
            updated["metadata"].get("flow_artifact_key") != "dendrogram"):
            updated["metadata"]["display_in_main"] = False
        exported.append(updated)
    project = source.model_dump(mode="json")
    project["output_root"] = "."
    project_file = write_json(directory / "project.json", project)
    catalog_file = write_json(directory / "artifacts.json", {"artifacts": exported})
    descriptor = {"project_manifest": "project.json", "artifact_catalog": "artifacts.json",
                  "report_require_plot_data": True,
                  "report_contrast_ids": plan.get("contrast_ids", [row["contrast_id"] for row in project.get("contrasts", [])]),
                  "flow_collection": {"plan_id": plan["plan_id"], "source_project_manifest_sha256": source_sha}}
    report = plan.get("report", {})
    for key in ("report_contrast_ids", "report_contrast_groups", "mode"):
        if key in report:
            descriptor[key] = report[key]
        elif key in plan:
            descriptor[key] = plan[key]
    descriptor["report_analyses"] = plan.get("configuration", {}).get("analyses", {})
    descriptor_file = write_json(directory / "report-project.json", descriptor)
    provenance_file = write_json(directory / "transport-provenance.json", {"products": provenance})
    plan_file = write_json(directory / "plan.json", plan)
    evidence = []
    for key, (path, _raw) in receipts.items():
        destination = directory / "evidence" / _safe(key) / "task.receipt.json"
        destination.parent.mkdir(parents=True, exist_ok=True)
        from .file_io import retain_owned_file
        retain_owned_file(path, destination)
        evidence.append(destination)
    from .presentation import build_presentation_assets
    rendered = build_presentation_assets(directory, plan, receipts, provenance)
    presentation_files = rendered.get("files", [])
    presentation_mapping = {str(path.resolve()): path.relative_to(directory).as_posix()
                            for path in [*copied, *evidence, *presentation_files]}
    for path in [*copied, *evidence, *presentation_files]:
        parent = path.parent
        while parent.is_relative_to(directory):
            presentation_mapping[str(parent.resolve())] = parent.relative_to(directory).as_posix()
            if parent == directory:
                break
            parent = parent.parent
    descriptor["extra_assets"] = _rewrite(rendered.get("extra_assets", []), presentation_mapping)
    descriptor["extra_views"] = _rewrite(rendered.get("extra_views", []), presentation_mapping)
    write_json(descriptor_file, descriptor)
    ledger = []
    for path in [*copied, *evidence, *presentation_files, project_file, catalog_file, descriptor_file, provenance_file, plan_file]:
        record = file_record(path)
        relative = path.relative_to(directory).as_posix()
        ledger.append({**record, "path": relative, "relative_path": relative})
    receipt = {"schema": COLLECTION_SCHEMA, "status": "completed", "plan_id": plan["plan_id"],
               "request_sha256": request_identity, "project_id": source.project_id,
               "source_project_manifest_sha256": source_sha,
               "tasks": [{"task_id": key, "status": raw["status"], "receipt_sha256": raw["receipt_sha256"]}
                         for key, (_path, raw) in sorted(receipts.items())],
               "task_receipts": [{"task_id": key, **file_record(path),
                                  "collection_path": (Path("evidence") / _safe(key) / "task.receipt.json").as_posix()}
                                 for key, (path, _raw) in sorted(receipts.items())],
               "artifact_count": len(exported), "outputs": ledger,
               "project_manifest": project_file.name, "artifact_catalog": catalog_file.name,
               "report_descriptor": descriptor_file.name,
               "completed_at": datetime.now(timezone.utc).isoformat(),
               "scope": "Verified task products collected; strict fixed-report acceptance is a separate gate"}
    receipt["receipt_sha256"] = canonical_hash(receipt)
    write_json(existing, receipt)
    validate_collection(existing)
    return existing


collect_tasks = collect_results

__all__ = ["collect_results", "collect_tasks", "validate_collection", "materialize_report_descriptor"]
