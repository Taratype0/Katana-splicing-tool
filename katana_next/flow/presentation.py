"""Approved report views from verified portable flow products.

This module never fits a network, rescans RNA or repeats association/enrichment
tests. Execution receipts and unchanged numerical families are prerequisites.
"""
from __future__ import annotations

import json
from pathlib import Path
import re

from .model import load_plan
from .runtime import canonical_hash, file_hash, file_record, read_json, write_json


NETWORK_PROOF_SCHEMA = "katana.flow.network_view_inputs/v1"
NETWORK_SCOPE = "flow_verified_network_complete"
MOTIF_PROOF_SCHEMA = "katana.flow.motif_view_inputs/v1"


def _safe(value):
    return re.sub(r"[^A-Za-z0-9._-]", "_", str(value))[:140] or "result"


def _inside(root, relative):
    relative = Path(relative)
    path = (root / relative).resolve()
    if relative.is_absolute() or not path.is_relative_to(root.resolve()):
        raise ValueError("Presentation proof path escapes its verified collection")
    return path


def _task_proof(root, record, plan):
    path = _inside(root, record["path"])
    if file_hash(path) != record["sha256"]:
        raise ValueError("Presentation task completion proof changed")
    value = read_json(path)
    if (value.get("schema") != "katana.flow.task_receipt/v1" or
        value.get("status") not in {"completed", "reused"} or
        value.get("plan_id") != plan["plan_id"] or
        value.get("receipt_sha256") != canonical_hash({k: v for k, v in value.items() if k != "receipt_sha256"})):
        raise ValueError("Presentation requires an actually completed flow task")
    task = next((item for item in plan["tasks"] if item["id"] == value.get("task_id")), None)
    if task is None or value.get("task_contract_sha256") != canonical_hash(task):
        raise ValueError("Presentation task does not match the sealed calculation plan")
    if not value.get("generated_artifacts") or not value.get("outputs"):
        raise ValueError("Presentation task has no actual calculated products")
    if value.get("source_project_manifest_sha256") != plan["project_manifest_sha256"]:
        raise ValueError("Presentation task belongs to another study revision")
    return task, value, path


def _verified_view_files(root, proof, task_values):
    verified = {}
    for item in proof["files"]:
        path = _inside(root, item["path"])
        actual = file_record(path)
        if (actual["sha256"], actual["size_bytes"]) != (item["sha256"], item["size_bytes"]):
            raise ValueError("Completed presentation numerical member changed")
        task = task_values[item["task_id"]]
        original = next((row for row in task["outputs"] if row["path"] == item["source_path"]), None)
        if original is None or (original["sha256"], original["size_bytes"]) != (item["sha256"], item["size_bytes"]):
            raise ValueError("Presentation input is not an unchanged product of its declared completed task")
        verified[path] = item
    return verified


def validate_network_view_proof(completion_path, input_paths, contrasts):
    """Authenticate local/cluster task receipts without inventing Slurm records."""
    path = Path(completion_path).resolve()
    proof = read_json(path)
    if proof.get("schema") != NETWORK_PROOF_SCHEMA or proof.get("status") != "completed":
        raise ValueError("Expected a complete portable network presentation proof")
    root = path.parent.parent
    plan = load_plan(root / "plan.json")
    if proof.get("plan_id") != plan["plan_id"] or proof.get("contrasts") != contrasts:
        raise ValueError("Network presentation has another plan or comparison population")
    required = {"network.fit", "network.bootstrap", "network.associate", "network.annotate"}
    task_values, records, paths = {}, {}, [path, root / "plan.json"]
    for record in proof["task_receipts"]:
        task, value, location = _task_proof(root, record, plan)
        adapter = task.get("adapter")
        if adapter not in required or adapter in records:
            raise ValueError("Network presentation requires one exact lifecycle for each required stage")
        records[adapter] = (task, value)
        task_values[task["id"]] = value
        paths.append(location)
    if set(records) != required:
        raise ValueError("Network fit, bootstrap, association and annotation must all actually complete")
    fit_task, fit = records["network.fit"]
    engine = fit.get("execution_evidence", {}).get("engine_root")
    expected = {"network.fit": "network_observed", "network.bootstrap": "consensus_modules",
                "network.associate": "associate_observed_declared_contrasts", "network.annotate": "annotate_modules"}
    for adapter, (task, value) in records.items():
        evidence = value.get("execution_evidence", {})
        if expected[adapter] not in evidence.get("executed_stages", []):
            raise ValueError("A network stage has no actual native calculation evidence")
        if adapter != "network.fit" and (fit_task["id"] not in task.get("dependencies", []) or
                                         evidence.get("source_engine_root") != engine):
            raise ValueError("Network postprocessing belongs to another fitted module population")
    verified = _verified_view_files(root, proof, task_values)
    if not set(map(lambda item: Path(item).resolve(), input_paths)) <= set(verified):
        raise ValueError("Every plotted network population must be in the exact completed task ledgers")
    precomputed = {name: _inside(root, relative) for name, relative in proof["precomputed"].items()}
    if not set(precomputed.values()) <= set(verified):
        raise ValueError("Complete association and enrichment families must be independently verified")
    for name in ("linkage", "bootstrap_members", "bootstrap_modules"):
        if _inside(root, proof[name]) not in verified:
            raise ValueError("Network tree and module bootstrap need their actual saved numerical products")
    precomputed["provenance"] = {"source_task_id": fit_task["id"], "flow_plan_id": plan["plan_id"]}
    return NETWORK_SCOPE, None, precomputed, [*paths, *verified]


def validate_motif_view_proof(completion_path, source, comparison_id):
    path = Path(completion_path).resolve()
    proof = read_json(path)
    if proof.get("schema") != MOTIF_PROOF_SCHEMA or proof.get("status") != "completed":
        raise ValueError("Expected a complete portable MotifTrace presentation proof")
    root = path.parent.parent
    plan = load_plan(root / "plan.json")
    if proof.get("plan_id") != plan["plan_id"] or proof.get("comparison_id") != comparison_id:
        raise ValueError("MotifTrace presentation has another plan or comparison")
    task, value, _location = _task_proof(root, proof["task_receipt"], plan)
    if task["workflow_id"] != "prioritize.motiftrace_consequences" or task.get("scope", {}).get("contrast_ids") != [comparison_id]:
        raise ValueError("MotifTrace presentation needs the exact unsplit comparison task")
    if value.get("job", {}).get("status") not in {"complete", "cache_hit"}:
        raise ValueError("MotifTrace scientific workflow has no completed actual job record")
    verified = _verified_view_files(root, proof, {task["id"]: value})
    marker = Path(source).resolve() / "COMPLETE.json"
    if marker not in verified:
        raise ValueError("MotifTrace completion marker was not produced by its completed comparison task")
    complete = read_json(marker)
    for name in complete.get("members", {}):
        if (marker.parent / name).resolve() not in verified:
            raise ValueError("MotifTrace presentation omitted a member of the complete scientific family")
    return proof


def build_presentation_assets(root, plan, receipts, transport):
    """Create the approved figures and browsers from completed task products.

    Called once while creating a new collection. It preserves all scientific
    bytes and imports already computed full test families without new fits.
    """
    root = Path(root).resolve()
    project = read_json(root / "project.json")
    project_sha = file_hash(root / "project.json")
    translated = {item["source"]["path"]: (root / item["collection_path"]).resolve() for item in transport}
    directory = root / "presentation"
    results = {"extra_assets": [], "extra_views": [], "files": []}

    def task_record(identity):
        path = root / "evidence" / _safe(identity) / "task.receipt.json"
        return {"path": path.relative_to(root).as_posix(), "sha256": file_hash(path)}

    def translated_path(path):
        original = str(Path(path).resolve())
        if original not in translated:
            raise ValueError("A report view input was not transported from a verified completed task: " + original)
        return translated[original]

    def proof_member(task_id, original):
        path = translated_path(original)
        return {**file_record(path), "path": path.relative_to(root).as_posix(),
                "source_path": str(Path(original).resolve()), "task_id": task_id}

    declared = {task["id"]: task for task in plan["tasks"]}
    selected = plan.get("contrast_ids", [c["contrast_id"] for c in project["contrasts"]])
    contrasts = [c for cid in selected for c in project["contrasts"] if c["contrast_id"] == cid]
    network = {task.get("adapter"): (task, receipts[task["id"]][1]) for task in plan["tasks"]
               if task.get("adapter", "").startswith("network.")}
    if network:
        required = {"network.fit", "network.bootstrap", "network.associate", "network.annotate"}
        if set(network) != required:
            raise ValueError("The fixed network report needs all four completed network stages")
        fit_task, fit = network["network.fit"]
        engine = Path(fit["execution_evidence"]["engine_root"])
        names = {"assignments_path": "04_network/observed/event_module_assignments.tsv",
                 "scores_path": "04_network/observed/module_program_scores.tsv",
                 "psi_path": "04_network/observed/selected_event_matrix.tsv",
                 "annotation_path": "02_matrix/event_metadata.tsv",
                 "samples_path": "02_matrix/sample_metadata_resolved.tsv",
                 "edges_path": "04_network/observed/network_edges.tsv"}
        native = {key: engine / relative for key, relative in names.items()}
        members = [proof_member(fit_task["id"], path) for path in native.values()]
        precomputed = {}
        for adapter in ("network.associate", "network.annotate"):
            task, value = network[adapter]
            for key, path in value["execution_evidence"].get("precomputed_views", {}).items():
                members.append(proof_member(task["id"], path))
                precomputed[key] = translated_path(path).relative_to(root).as_posix()
        if set(precomputed) != {"full_tests", "selected_tests", "trait_audit", "tested", "background", "input_genes"}:
            raise ValueError("Network views require the complete already computed contrast and ORA families")
        linkage = engine / "04_network/observed/event_linkage.npy"
        members.append(proof_member(fit_task["id"], linkage))
        bootstrap_task, bootstrap = network["network.bootstrap"]
        boot_root = Path(bootstrap["execution_evidence"]["engine_root"]) / "04c_consensus_modules"
        boot_members = boot_root / "consensus_event_module_assignments.tsv"
        boot_modules = boot_root / "consensus_module_summary.tsv"
        members += [proof_member(bootstrap_task["id"], path) for path in (boot_members, boot_modules)]
        proof = {"schema": NETWORK_PROOF_SCHEMA, "status": "completed", "plan_id": plan["plan_id"],
                 "contrasts": contrasts, "source_project_manifest_sha256": plan["project_manifest_sha256"],
                 "task_receipts": [task_record(task["id"]) for task, _value in network.values()],
                 "files": members, "precomputed": precomputed,
                 "linkage": translated_path(linkage).relative_to(root).as_posix(),
                 "bootstrap_members": translated_path(boot_members).relative_to(root).as_posix(),
                 "bootstrap_modules": translated_path(boot_modules).relative_to(root).as_posix()}
        proof_path = write_json(directory / "network-inputs.json", proof)
        from katana_next.reports.network_views import build_network_assets
        rendered = build_network_assets(**{key: translated_path(path) for key, path in native.items()},
                                       contrasts=contrasts, completion_manifest=proof_path,
                                       output_dir=directory / "network", project_id=project["project_id"])
        for artifact in rendered["artifacts"]:
            artifact.setdefault("metadata", {})["source_project_manifest_sha256"] = project_sha
        results["extra_assets"].extend(rendered["artifacts"])
        results["extra_views"].append({"project_id": project["project_id"], "source_project_manifest_sha256": project_sha,
                                       "status": "complete", "unavailable": rendered["unavailable"]})

    for identity, (_receipt_path, value) in receipts.items():
        task = declared[identity]
        if task["workflow_id"] != "prioritize.motiftrace_consequences" or value["status"] == "not_applicable":
            continue
        scopes = task.get("scope", {}).get("contrast_ids", [])
        if len(scopes) != 1:
            raise ValueError("MotifTrace report rendering requires an explicit single comparison scope")
        markers = [item["path"] for item in value["outputs"] if Path(item["path"]).name == "COMPLETE.json"]
        if len(markers) != 1:
            raise ValueError("MotifTrace task must supply exactly one complete unsplit scientific family")
        source_original = Path(markers[0]).parent
        marker = translated_path(markers[0])
        complete = read_json(marker)
        originals = [source_original / name for name in complete["members"]]
        proof = {"schema": MOTIF_PROOF_SCHEMA, "status": "completed", "plan_id": plan["plan_id"],
                 "comparison_id": scopes[0], "project_id": project["project_id"],
                 "source_project_manifest_sha256": plan["project_manifest_sha256"],
                 "render_project_manifest_sha256": project_sha,
                 "task_receipt": task_record(identity),
                 "complete_manifest": {"sha256": file_hash(marker), "bytes": marker.stat().st_size},
                 "files": [proof_member(identity, original) for original in [*originals, source_original / "COMPLETE.json"]]}
        proof_path = write_json(directory / ("motif-" + _safe(identity) + "-inputs.json"), proof)
        from katana_next.reports.motif_views import build_completed_motif_views
        rendered = build_completed_motif_views(marker.parent, directory / ("motif-" + _safe(identity)),
                    comparison_id=scopes[0], project_id=project["project_id"],
                    source_project_manifest_sha256=project_sha, completion_proof=proof_path, proof_kind="flow")
        results["extra_views"].append(rendered)
    if any(task.get("adapter") == "mapp.maei" for task in plan["tasks"]):
        from .mapp_views import build_mapp_views
        mapp = build_mapp_views(root, plan, receipts, transport)
        results["extra_views"].extend(mapp.get("extra_views", []))
    if any(task.get("adapter") == "qc.sequencing" for task in plan["tasks"]):
        from .qc_views import build_qc_views
        qc = build_qc_views(root, plan, receipts, transport)
        results["extra_views"].extend(qc.get("extra_views", []))
        results["files"].extend(qc.get("files", []))
    if directory.is_dir():
        results["files"] = list(dict.fromkeys([*results["files"], *[path for path in directory.rglob("*") if path.is_file()]]))
    return results
