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


def network_bootstrap_enabled(plan):
    """Only a sealed, explicit false can waive historical bootstrap products."""
    value = plan.get("configuration", {}).get("network", {}).get("bootstrap", {}).get("enabled", True)
    if not isinstance(value, bool):
        raise ValueError("The sealed network bootstrap policy must be boolean")
    return value


def _network_required_stages(plan):
    return {"network.fit", "network.associate", "network.annotate"} | (
        {"network.bootstrap"} if network_bootstrap_enabled(plan) else set())


def _network_presentation_policy(plan):
    enabled = network_bootstrap_enabled(plan)
    return {"enabled": enabled, "status": "completed" if enabled else "disabled_by_request",
            "hub_metric": plan.get("configuration", {}).get("network", {}).get("parameters", {}).get(
                "hub_metric", "abs_kME" if not enabled else "legacy_recorded_metrics")}


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


def _network_task_identity(task, project):
    """Re-derive a scoped task's samples from the immutable study declaration."""
    parameters = task.get("parameters", {})
    if "network_scope" not in parameters:
        return None  # Preserve historical single-cohort plans and their proofs.
    from types import SimpleNamespace
    from .network_inputs import network_scope_identity
    study = SimpleNamespace(
        samples=[SimpleNamespace(sample_id=row["sample_id"], biological_sample_id=row["biological_sample_id"])
                 for row in project["samples"]],
        contrasts=[SimpleNamespace(contrast_id=row["contrast_id"],
                    numerator=SimpleNamespace(sample_ids=row["numerator"]["sample_ids"]),
                    denominator=SimpleNamespace(sample_ids=row["denominator"]["sample_ids"]))
                   for row in project["contrasts"]])
    identity = network_scope_identity(study, parameters)
    if any(task.get("scope", {}).get(key) != value for key, value in identity.items()):
        raise ValueError("Network presentation task parameters and sample scope disagree")
    return identity


def _network_lifecycles(plan, receipts):
    """Group the explicitly selected stages by their actual fit dependency."""
    network = [task for task in plan["tasks"] if task.get("adapter", "").startswith("network.")]
    fits = {task["id"]: {"network.fit": (task, receipts[task["id"]][1])}
            for task in network if task["adapter"] == "network.fit"}
    for task in network:
        if task["adapter"] == "network.fit":
            continue
        parents = [identifier for identifier in task.get("dependencies", []) if identifier in fits]
        if len(parents) != 1 or task["adapter"] in fits[parents[0]]:
            raise ValueError("Each network followup needs one unique fitted lifecycle")
        fits[parents[0]][task["adapter"]] = task, receipts[task["id"]][1]
    required = _network_required_stages(plan)
    if any(set(stages) != required for stages in fits.values()):
        message = ("all four completed network stages for every fit" if network_bootstrap_enabled(plan)
                   else "every selected completed network stage for every fit")
        raise ValueError("The fixed network report needs " + message)
    return fits


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
    required = _network_required_stages(plan)
    policy = _network_presentation_policy(plan)
    if ("bootstrap" in proof and proof["bootstrap"] != policy) or (not policy["enabled"] and "bootstrap" not in proof):
        raise ValueError("Network presentation bootstrap policy differs from the sealed plan")
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
        raise ValueError("Every selected network fit, association, annotation and optional bootstrap must actually complete")
    fit_task, fit = records["network.fit"]
    identity = _network_task_identity(fit_task, read_json(root / "project.json"))
    if identity != proof.get("network_identity"):
        raise ValueError("Network presentation identity differs from its fitted sample population")
    engine = fit.get("execution_evidence", {}).get("engine_root")
    name = str(fit_task.get("parameters", {}).get("network_name", "observed"))
    expected = {"network.fit": "network_" + name, "network.bootstrap": "consensus_modules",
                "network.associate": "associate_" + name + "_declared_contrasts", "network.annotate": "annotate_modules"}
    for adapter, (task, value) in records.items():
        evidence = value.get("execution_evidence", {})
        if identity and (_network_task_identity(task, read_json(root / "project.json")) != identity or
                         evidence.get("network_identity") != identity):
            raise ValueError("Network presentation stages belong to different contrast/sample scopes")
        if expected[adapter] not in evidence.get("executed_stages", []):
            raise ValueError("A network stage has no actual native calculation evidence")
        if adapter != "network.fit" and (fit_task["id"] not in task.get("dependencies", []) or
                                         evidence.get("source_engine_root") != engine):
            raise ValueError("Network postprocessing belongs to another fitted module population")
    verified = _verified_view_files(root, proof, task_values)
    if identity:
        import pandas as pd
        marker = _inside(root, proof["scope_marker"])
        sample_metadata = _inside(root, proof["sample_metadata"])
        if marker not in verified or sample_metadata not in verified or read_json(marker) != identity:
            raise ValueError("Scoped network presentation needs the actual verified fit scope marker and samples")
        samples = pd.read_csv(sample_metadata, sep="\t", dtype=str)
        if (not samples.sample_id.is_unique or set(samples.sample_id) != set(identity["sample_ids"])):
            raise ValueError("Network report contains missing or foreign fitted samples")
        if ("biological_sample_id" not in samples or
            samples.set_index("sample_id").biological_sample_id.to_dict() !=
                dict(zip(identity["sample_ids"], identity["biological_sample_ids"]))):
            raise ValueError("Network report biological units differ from the fitted sample population")
        if [row["contrast_id"] for row in contrasts] != identity["contrast_ids"]:
            raise ValueError("Network report comparisons differ from this fitted scope")
    if not set(map(lambda item: Path(item).resolve(), input_paths)) <= set(verified):
        raise ValueError("Every plotted network population must be in the exact completed task ledgers")
    precomputed = {name: _inside(root, relative) for name, relative in proof["precomputed"].items()}
    if not set(precomputed.values()) <= set(verified):
        raise ValueError("Complete association and enrichment families must be independently verified")
    hub_keys = {"hub_events", "hub_genes", "hub_parameters"}
    if not policy["enabled"] and not hub_keys <= set(precomputed):
        raise ValueError("Absolute kME report ranking needs its actual completed Hub numerical family")
    if hub_keys <= set(precomputed):
        hub = read_json(precomputed["hub_parameters"])
        if (hub.get("schema") != "katana.network.hub_membership/v1" or
            hub.get("hub_metric") != policy["hub_metric"] or hub.get("network_refitted") is not False or
            hub.get("signed_kME_retained") is not True):
            raise ValueError("Completed Hub metadata differs from the declared observed kME policy")
        required_sources = {"assignments", "psi", "activity", "annotation"}
        if set(hub.get("sources", {})) != required_sources:
            raise ValueError("Observed Hub metrics lack their complete fitted source identities")
        source_names = {"assignments": f"04_network/{name}/event_module_assignments.tsv",
                        "psi": f"04_network/{name}/selected_event_matrix.tsv",
                        "activity": f"04_network/{name}/module_program_scores.tsv",
                        "annotation": "02_matrix/event_metadata.tsv"}
        annotation_engine = records["network.annotate"][1].get("execution_evidence", {}).get("engine_root")
        for key, source in hub["sources"].items():
            original = str((Path(engine) / source_names[key]).resolve())
            allowed_paths = {original}
            if annotation_engine:
                allowed_paths.add(str((Path(annotation_engine) / source_names[key]).resolve()))
            if (str(Path(source["path"]).resolve()) not in allowed_paths or
                not any(member["source_path"] == original and member["sha256"] == source["sha256"]
                        and member["task_id"] == fit_task["id"] for member in proof["files"])):
                raise ValueError("Observed Hub metrics belong to another fitted numerical population")
    numerical = ("linkage", "bootstrap_members", "bootstrap_modules") if policy["enabled"] else ("linkage",)
    if not policy["enabled"] and {"bootstrap_members", "bootstrap_modules"} & set(proof):
        raise ValueError("Disabled bootstrap cannot be represented by completed bootstrap products")
    for name in numerical:
        if _inside(root, proof[name]) not in verified:
            raise ValueError("Network tree and module bootstrap need their actual saved numerical products")
    precomputed["provenance"] = {"source_task_id": fit_task["id"], "flow_plan_id": plan["plan_id"],
                                 "network_identity": identity, "bootstrap": policy}
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
    lifecycles = _network_lifecycles(plan, receipts)
    for fit_id, network in lifecycles.items():
        fit_task, fit = network["network.fit"]
        identity = _network_task_identity(fit_task, project)
        scoped_contrasts = ([c for cid in identity["contrast_ids"] for c in contrasts if c["contrast_id"] == cid]
                            if identity else contrasts)
        if identity and [c["contrast_id"] for c in scoped_contrasts] != identity["contrast_ids"]:
            raise ValueError("A fitted network comparison is absent from the selected report")
        if not identity and len(lifecycles) != 1:
            raise ValueError("Multiple network fits need explicit contrast/cohort identities")
        suffix = "-" + _safe(fit_id) + "-" + canonical_hash(fit_id)[:8] if identity else ""
        display_contrast = identity["contrast_ids"][0] if identity and identity["network_scope"] == "contrast" else None
        name = str(fit_task.get("parameters", {}).get("network_name", "observed"))
        engine = Path(fit["execution_evidence"]["engine_root"])
        names = {"assignments_path": f"04_network/{name}/event_module_assignments.tsv",
                 "scores_path": f"04_network/{name}/module_program_scores.tsv",
                 "psi_path": f"04_network/{name}/selected_event_matrix.tsv",
                 "annotation_path": "02_matrix/event_metadata.tsv",
                 "samples_path": "02_matrix/sample_metadata_resolved.tsv",
                 "edges_path": f"04_network/{name}/network_edges.tsv"}
        native = {key: engine / relative for key, relative in names.items()}
        members = [proof_member(fit_task["id"], path) for path in native.values()]
        precomputed = {}
        for adapter in ("network.associate", "network.annotate"):
            task, value = network[adapter]
            for key, path in value["execution_evidence"].get("precomputed_views", {}).items():
                members.append(proof_member(task["id"], path))
                precomputed[key] = translated_path(path).relative_to(root).as_posix()
        base_keys = {"full_tests", "selected_tests", "trait_audit", "tested", "background", "input_genes"}
        hub_keys = {"hub_events", "hub_genes", "hub_parameters"}
        if (set(precomputed) not in (base_keys, base_keys | hub_keys) or
                (not network_bootstrap_enabled(plan) and not hub_keys <= set(precomputed))):
            raise ValueError("Network views require the complete already computed contrast and ORA families")
        linkage = engine / f"04_network/{name}/event_linkage.npy"
        members.append(proof_member(fit_task["id"], linkage))
        proof = {"schema": NETWORK_PROOF_SCHEMA, "status": "completed", "plan_id": plan["plan_id"],
                 "contrasts": scoped_contrasts, "source_project_manifest_sha256": plan["project_manifest_sha256"],
                 "task_receipts": [task_record(task["id"]) for task, _value in network.values()],
                 "files": members, "precomputed": precomputed,
                 "bootstrap": _network_presentation_policy(plan),
                 "linkage": translated_path(linkage).relative_to(root).as_posix()}
        if network_bootstrap_enabled(plan):
            bootstrap_task, bootstrap = network["network.bootstrap"]
            boot_root = Path(bootstrap["execution_evidence"]["engine_root"]) / "04c_consensus_modules"
            boot_members = boot_root / "consensus_event_module_assignments.tsv"
            boot_modules = boot_root / "consensus_module_summary.tsv"
            members += [proof_member(bootstrap_task["id"], path) for path in (boot_members, boot_modules)]
            proof.update(bootstrap_members=translated_path(boot_members).relative_to(root).as_posix(),
                         bootstrap_modules=translated_path(boot_modules).relative_to(root).as_posix())
        if identity:
            marker = engine / "network_scope.json"
            members.append(proof_member(fit_id, marker))
            proof.update(network_identity=identity,
                         scope_marker=translated_path(marker).relative_to(root).as_posix(),
                         sample_metadata=translated_path(native["samples_path"]).relative_to(root).as_posix())
        proof_path = write_json(directory / ("network" + suffix + "-inputs.json"), proof)
        from katana_next.reports.network_views import build_network_assets
        rendered = build_network_assets(**{key: translated_path(path) for key, path in native.items()},
                                       contrasts=scoped_contrasts, completion_manifest=proof_path,
                                       output_dir=directory / ("network" + suffix), project_id=project["project_id"])
        for artifact in rendered["artifacts"]:
            artifact.setdefault("metadata", {})["source_project_manifest_sha256"] = project_sha
            if identity:
                artifact["artifact_id"] += "-" + canonical_hash(identity)[:16]
                artifact["metadata"]["network_identity"] = identity
                artifact["metadata"]["interpretation_status"] = "minimum_met"
            if display_contrast:
                artifact["contrast_ids"] = [display_contrast]
                artifact["name"] += " | " + display_contrast
                artifact["metadata"]["caption"] += (f" Independent {display_contrast} network; "
                    f"{identity['n_samples']} declared biological samples. "
                    + ("Actual bootstrap diagnostics and continuous membership measures are retained."
                       if network_bootstrap_enabled(plan) else
                       "Hub ranking uses absolute observed kME; bootstrap was disabled by request and no resampling confidence is inferred."))
                for link in artifact["metadata"].get("data_links", []):
                    link["contrast_id"] = display_contrast
        for entry in rendered["unavailable"]:
            if display_contrast:
                entry["contrast_id"] = display_contrast
                entry["title"] += " | " + display_contrast
        # Persist the same identities returned to collection; an exported
        # registry must not revert to cohort-wide labels after relocation.
        write_json(directory / ("network" + suffix) / "artifact-records.json", rendered["artifacts"])
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
