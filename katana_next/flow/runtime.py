"""Isolated, blocking execution of the existing Katana analysis workflows.

Nextflow owns scheduling.  Every invocation owns a separate project catalogue;
parents are verified immutable files, never a shared writable DuckDB database.
"""
from __future__ import annotations

from datetime import datetime, timezone
from contextlib import contextmanager
import hashlib
import json
import os
from pathlib import Path
import shutil
from typing import Any
from threading import RLock
from uuid import NAMESPACE_URL, uuid5


RECEIPT_SCHEMA = "katana.flow.task_receipt/v1"
_MOTIF_DECLARATION_LOCK = RLock()


@contextmanager
def _declare_complete_motif_family(project, parameters):
    """Record the execution boundary before the unchanged kernel scores RNA."""
    from unittest.mock import patch
    import katana_next.workflows.motiftrace_consequences as workflow_module
    from katana_next.workflows.sources import selected_source
    original = workflow_module.analyze_consequences
    source = selected_source(project.context, parameters)
    policy = source.metadata.get("read_support_policy") if source is not None else None

    def run(*args, **kwargs):
        provenance = dict(kwargs.get("provenance") or {})
        provenance.update(statistical_family_split=False, read_support_policy=policy,
                          execution_boundary="one complete comparison; all predeclared models and tests",
                          flow_execution_adapter_sha256=file_hash(__file__))
        kwargs["provenance"] = provenance
        return original(*args, **kwargs)

    with _MOTIF_DECLARATION_LOCK, patch.object(workflow_module, "analyze_consequences", run):
        yield


def canonical_hash(value: Any) -> str:
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"),
                                     ensure_ascii=False, default=str).encode()).hexdigest()


def file_hash(path: str | Path) -> str:
    with Path(path).open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def read_json(path: str | Path) -> dict:
    value = json.loads(Path(path).read_text(encoding="utf-8-sig"))
    if not isinstance(value, dict):
        raise ValueError(f"Expected a JSON object: {path}")
    return value


def write_json(path: str | Path, value: Any) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, ensure_ascii=False, default=str),
                         encoding="utf-8")
    os.replace(temporary, path)
    return path


def file_record(path: str | Path) -> dict:
    path = Path(path).resolve()
    if not path.is_file():
        raise FileNotFoundError(path)
    before = path.stat()
    digest = file_hash(path)
    after = path.stat()
    if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
        raise ValueError(f"File changed while its content was verified: {path}")
    return {"path": str(path), "size_bytes": after.st_size, "sha256": digest}


def verify_files(records: list[dict]) -> None:
    for record in records:
        actual = file_record(record["path"])
        if actual["size_bytes"] != record.get("size_bytes") or actual["sha256"] != record.get("sha256"):
            raise ValueError(f"Verified content changed: {record['path']}")


def load_plan_task(plan_path: str | Path, task_id: str, project_path=None):
    path = Path(plan_path).expanduser().resolve()
    from .model import load_plan
    plan = load_plan(path)
    tasks = plan.get("tasks", [])
    if not isinstance(tasks, list) or len({t.get("id") for t in tasks}) != len(tasks):
        raise ValueError("Plan tasks require distinct IDs")
    matches = [t for t in tasks if t.get("id") == task_id]
    if len(matches) != 1:
        raise ValueError(f"Task is not declared in this plan: {task_id}")
    manifest = Path(project_path or plan["project_manifest"]).expanduser()
    manifest = (manifest if manifest.is_absolute() else path.parent / manifest).resolve()
    return path, plan, matches[0], manifest


def validate_receipt(path: str | Path, *, plan: dict | None = None,
                     task: dict | None = None, verify_inputs: bool = True) -> dict:
    path = Path(path)
    if path.is_dir():
        path = path / "task.receipt.json"
    value = read_json(path)
    if value.get("schema") != RECEIPT_SCHEMA or value.get("status") not in {"completed", "reused", "not_applicable"}:
        raise ValueError(f"Task has no successful or justified not-applicable receipt: {path}")
    signed = {k: v for k, v in value.items() if k != "receipt_sha256"}
    if value.get("receipt_sha256") != canonical_hash(signed):
        raise ValueError(f"Task receipt content identity is invalid: {path}")
    if plan is not None and value.get("plan_id") != plan.get("plan_id"):
        raise ValueError("Task receipt belongs to a different execution plan")
    if task is not None and (value.get("task_id") != task.get("id") or
                             value.get("task_contract_sha256") != canonical_hash(task)):
        raise ValueError("Task receipt belongs to a different task contract")
    if not value.get("generated_artifacts") or not value.get("outputs"):
        raise ValueError("A completed task requires real, verified products")
    verify_files(value["outputs"])
    if verify_inputs:
        verify_files(value.get("inputs", []))
    verify_files(value.get("producer_markers", []))
    for artifact in value.get("artifacts", []):
        if artifact.get("status") != "complete" or not Path(artifact.get("path", "")).is_file():
            raise ValueError("A task receipt contains an incomplete or unavailable artifact")
    return value


def _parent_receipts(paths, plan, task, manifest_sha):
    observed = {}
    expected = set(task.get("dependencies", []))
    declared = {row["id"]: row for row in plan["tasks"]}
    for path in paths:
        location = Path(path)
        if location.is_dir():
            location /= "task.receipt.json"
        header = read_json(location)
        identity = header.get("task_id")
        if identity not in declared or header.get("plan_id") != plan.get("plan_id"):
            raise ValueError("Dependency inputs include an unknown task or execution plan")
        # A level-wide Nextflow channel can supply earlier, unrelated receipts.
        # Only the declared parents are consumed or have their large ledgers
        # re-hashed. Their transitive products already travel in each receipt.
        if identity not in expected:
            continue
        value = validate_receipt(location, plan=plan, task=declared[identity])
        if value.get("source_project_manifest_sha256") != manifest_sha:
            raise ValueError("Dependency receipt belongs to a different source project revision")
        if value["task_id"] in observed:
            raise ValueError("A dependency receipt was supplied more than once")
        observed[value["task_id"]] = (location.resolve(), value)
    if set(observed) != expected:
        raise ValueError(f"Task dependencies differ from the plan; expected {sorted(expected)}, received {sorted(observed)}")
    return observed


def _source_manifest(path):
    from katana_next.infrastructure.manifest_store import ManifestStore
    return ManifestStore().load(path)[1]


def _import_parents(project, parents, manifest_sha):
    """Bind proved parent products to a private execution-only project copy.

    Only the snapshot token changes.  Scientific parameters, original producer
    identity, sample/contrast membership and files retain their original values.
    """
    from katana_next.domain.models import AnalysisArtifact
    artifacts, markers = {}, {}
    for _location, parent in parents.values():
        for marker in parent.get("producer_markers", []):
            markers[marker["path"]] = marker
        for raw in parent["artifacts"]:
            prior = artifacts.get(raw["artifact_id"])
            if prior is not None and prior != raw:
                # The private snapshot binding can differ between valid task
                # projections; the producer/file/science identities cannot.
                def identity(item):
                    item = json.loads(json.dumps(item))
                    meta = item.get("metadata", {})
                    meta.pop("project_snapshot_fingerprint", None)
                    meta.pop("flow_snapshot_binding", None)
                    return item
                if identity(prior) != identity(raw):
                    raise ValueError("Conflicting dependency records for one artifact identity")
            artifacts[raw["artifact_id"]] = raw
    for raw in artifacts.values():
        artifact = AnalysisArtifact.model_validate(raw)
        if artifact.project_id != project.project_id:
            raise ValueError("Parent artifact belongs to another project")
        meta = dict(artifact.metadata)
        meta["flow_snapshot_binding"] = {
            "source_project_manifest_sha256": manifest_sha,
            "producer_snapshot_fingerprint": meta.get("flow_snapshot_binding", {}).get(
                "producer_snapshot_fingerprint", meta.get("project_snapshot_fingerprint")),
            "source_artifact_fingerprint": artifact.fingerprint,
        }
        meta["project_snapshot_fingerprint"] = project.context.snapshot.manifest_fingerprint
        project.context.catalog.register_artifact(artifact.model_copy(update={"metadata": meta}))
        # Frozen source_complete intentionally checks a marker in this
        # catalogue's result store. Copy the authentic completed marker only.
        for marker in markers.values():
            original = read_json(marker["path"])
            if (original.get("workflow_id") == artifact.analysis_id and
                original.get("run_fingerprint") == artifact.run_fingerprint and
                artifact.artifact_id in {a["artifact_id"] for a in original.get("generated_artifacts", [])}):
                target = project.context.results.artifact_directory(artifact.analysis_id, artifact.run_fingerprint) / "parameters.json"
                if target.exists() and file_hash(target) != marker["sha256"]:
                    raise ValueError("Conflicting producer completion markers")
                if not target.exists():
                    shutil.copyfile(marker["path"], target)
                break
        else:
            if artifact.metadata.get("artifact_contract_version") == "3.0-correctness":
                raise ValueError("Dependency artifact has no verified completed producer marker")
    return artifacts


def _aggregate_parents(project, task):
    """Create an explicit, lossless view when a consumer needs all comparisons."""
    import pandas as pd
    from katana_next.workflows.base import ArtifactWriter
    groups = {}
    for artifact in project.artifacts():
        if artifact.kind != "table" or not artifact.analysis_id.startswith(("differential.", "quality.")):
            continue
        key = (artifact.analysis_id, artifact.name)
        groups.setdefault(key, []).append(artifact)
    created = []
    for (analysis, name), members in groups.items():
        if len(members) < 2:
            continue
        by_science = {}
        for member in members:
            params = {k: v for k, v in member.parameters.items()
                      if k not in {"contrast_ids", "source_artifact_ids", "quality_artifact_id", "selection_id"}}
            by_science.setdefault(canonical_hash(params), []).append(member)
        for compatible in by_science.values():
            if len(compatible) < 2:
                continue
            scopes = [set(a.contrast_ids or a.parameters.get("contrast_ids") or
                          [c.contrast_id for c in project.context.manifest.contrasts]) for a in compatible]
            if any(scopes[i] & scopes[j] for i in range(len(scopes)) for j in range(i)):
                # Alternate fits covering overlapping scopes must not be
                # silently treated as independent fragments of one result.
                continue
            frames = [pd.read_parquet(a.path) for a in compatible]
            columns = list(frames[0].columns)
            if any(list(frame.columns) != columns for frame in frames):
                raise ValueError("Cannot join parent tables with different column contracts")
            combined = pd.concat(frames, ignore_index=True)
            if "row_id" not in combined or combined["row_id"].astype(str).duplicated().any():
                raise ValueError("Parent table join would duplicate stable row identities")
            metadata = dict(compatible[0].metadata)
            for key in ("workflow_version", "read_support_policy", "das_selection"):
                if any(a.metadata.get(key) != metadata.get(key) for a in compatible):
                    raise ValueError("Cannot join parent tables with different method or read-support policies")
            fingerprint = canonical_hash({"operation": "lossless_contrast_concat/v1",
                                          "parents": [a.fingerprint for a in compatible]})
            parameters = dict(compatible[0].parameters)
            parameters["contrast_ids"] = [c.contrast_id for c in project.context.manifest.contrasts
                                           if any(c.contrast_id in scope for scope in scopes)]
            parameters.pop("source_artifact_ids", None)
            metadata["flow_lossless_parent_join"] = {"parent_artifact_ids": [a.artifact_id for a in compatible],
                                                    "statistics_recomputed": False}
            writer = ArtifactWriter(project.context, analysis, fingerprint, parameters,
                                    str(metadata.get("workflow_version", "")), common_metadata=metadata)
            result = writer.table(name, "joined.parquet", combined,
                                  entity_type=compatible[0].entity_type,
                                  parents=tuple(a.artifact_id for a in compatible))
            writer.finish()
            created.append(result)
    return created


def _selection(project, task, parameters):
    request = task.get("selection")
    if not request:
        return None
    import pandas as pd
    from katana_next.domain.models import SelectionSet
    from katana_next.workflows.sources import resolve_table
    from katana_next.science.das_selection import selected_event_mask
    scopes = list(request.get("contrast_ids") or parameters.get("contrast_ids") or [])
    query = {"contrast_ids": scopes} if scopes else {}
    source = resolve_table(project.context, request["analysis_id"], request.get("name"), query)
    if source is None:
        raise ValueError("Target selection requires a complete source table covering its comparisons")
    frame = pd.read_parquet(source.path)
    if scopes:
        if "comparison_id" not in frame:
            raise ValueError("Contrast selection requires explicit comparison_id values")
        frame = frame[frame["comparison_id"].astype(str).isin(scopes)]
    event_set = request.get("event_set", parameters.get("event_set", "stable_significant"))
    selected = frame[selected_event_mask(frame, event_set)]
    evidence = {"source_artifact_id": source.artifact_id, "source_rows": len(frame),
                "selected_rows": len(selected), "event_set": event_set, "contrast_ids": scopes}
    if selected.empty:
        return evidence
    ids = tuple(selected["row_id"].astype(str))
    identity = str(uuid5(NAMESPACE_URL, canonical_hash({"source": source.fingerprint, "ids": ids})))
    selection = SelectionSet(selection_id=identity, project_id=project.project_id,
                             name=f"Flow selection {task['id']}", entity_type=source.entity_type,
                             source_artifact_id=source.artifact_id, entity_ids=ids,
                             filters={"event_set": event_set, "contrast_ids": scopes})
    project.projects.save_selection(selection)
    parameters["selection_id"] = identity
    return evidence


def _stage_artifacts(project, task, directory, fingerprint):
    from .stages import execute_stage
    from katana_next.domain.models import AnalysisArtifact, ArtifactKind
    result = execute_stage(task, directory / "compute")
    if result.get("status") not in {"complete", "completed", "not_applicable"}:
        raise ValueError("Stage adapter did not complete its actual computation")
    outputs = result.get("outputs", [])
    if not outputs:
        raise ValueError("Stage adapter returned no actual products")
    artifacts = []
    analysis = task.get("workflow_id") or task.get("stage") or "flow.compute"
    for index, item in enumerate(outputs):
        verify_files([item])
        path = Path(item["path"]).resolve()
        if not path.is_relative_to(directory.resolve()):
            raise ValueError("Stage output must be owned by its isolated task directory")
        artifact = AnalysisArtifact(project_id=project.project_id,
            analysis_id=item.get("analysis_id", analysis),
            name=item.get("name", path.name), kind=ArtifactKind(item.get("kind", "file")),
            entity_type=item.get("entity_type", "stage_result"), run_fingerprint=fingerprint,
            fingerprint=canonical_hash(item), path=path,
            contrast_ids=tuple(task.get("scope", {}).get("contrast_ids", [])),
            parameters=task.get("parameters", {}),
            metadata={"flow_stage_evidence": result.get("evidence", {}),
                      "flow_stage": task.get("adapter"), "flow_artifact_key": item.get("artifact_key"),
                      "completion_verified": True, "project_snapshot_fingerprint": project.context.snapshot.manifest_fingerprint,
                      **item.get("metadata", {})})
        project.context.catalog.register_artifact(artifact)
        artifacts.append(artifact)
    result["status"] = "completed" if result["status"] == "complete" else result["status"]
    return artifacts, result


def _owned_parallel_inputs(manifest, task, directory):
    """Privately snapshot an explicitly requested completed import, unchanged.

    The frozen reader registers its completed source members as artifacts. Its
    manifest and relative member layout therefore need task-owned copies, while
    reuse must remain bound to the original declared resource fingerprints.
    """
    from types import SimpleNamespace
    from katana_next.workflows.motiftrace_parallel import MotifTraceParallelWorkflow
    workflow = MotifTraceParallelWorkflow()
    parameters = workflow.normalize_parameters(task.get("parameters") or {})
    context = SimpleNamespace(manifest=manifest)
    original_inputs = [file_record(path) for path in workflow.input_paths(context, parameters)]
    bundle, _ = workflow._load(context, parameters)
    asset = bundle.asset
    root = directory / "completed-parallel"
    root.mkdir()
    members = [bundle.manifest_path, *bundle.files.values()]
    owned = []
    for original in dict.fromkeys(members):
        target = root / original.relative_to(bundle.root)
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(original, target)
        source = next(item for item in original_inputs if item["path"] == str(original))
        if (file_hash(target), target.stat().st_size) != (source["sha256"], source["size_bytes"]):
            raise ValueError("Completed parallel source changed during its private snapshot")
        owned.append({"source": source, "path": str(target)})
    verify_files(original_inputs)
    copied_manifest = root / bundle.manifest_path.relative_to(bundle.root)
    copied_asset = asset.model_copy(update={"path": copied_manifest})
    assets = tuple(copied_asset if candidate.asset_id == asset.asset_id else candidate for candidate in manifest.assets)
    evidence = {"asset_id": asset.asset_id, "source_manifest": str(bundle.manifest_path),
                "source_manifest_sha256": file_hash(bundle.manifest_path),
                "private_manifest": str(copied_manifest), "members": owned,
                "statistics_recomputed": False, "source_modified": False}
    return manifest.model_copy(update={"assets": assets}), original_inputs, evidence


@contextmanager
def _explicit_parallel_scope(project, parameters):
    """Keep the validated single original contrast on this private importer.

    The generic service normally elides an all-project contrast selection. A
    completed parallel family requires its explicit identity even in a project
    containing only that one comparison. No other service instance is changed.
    """
    from unittest.mock import patch
    original = project.workflows._normalize_contrast_scope

    def retain(context, workflow, values):
        normalized = original(context, workflow, values)
        if workflow.definition.workflow_id == "prioritize.motiftrace_parallel" and "contrast_ids" not in normalized:
            requested = list(parameters.get("contrast_ids", []))
            if len(requested) != 1 or requested != [c.contrast_id for c in context.manifest.contrasts]:
                raise ValueError("Completed parallel import requires one validated original comparison")
            normalized["contrast_ids"] = requested
        return normalized

    with patch.object(project.workflows, "_normalize_contrast_scope", retain):
        yield


def execute_task(plan_path: str | Path, task_id: str, output_dir: str | Path,
                 dependency_receipts=(), *, project_path=None) -> Path:
    """Run one declared task to completion; a submitted/running job never passes.

    Existing valid products are returned untouched. Failed/unknown attempts are
    retained and require a new task directory, rather than deleting a claim.
    """
    from katana_next import __version__
    from katana_next.client import KatanaProject
    plan_location, plan, task, manifest_path = load_plan_task(plan_path, task_id, project_path)
    manifest_sha = file_hash(manifest_path)
    expected = plan.get("project_manifest_sha256", plan.get("source_project_manifest_sha256"))
    if expected and expected != manifest_sha:
        raise ValueError("Source project manifest changed after planning")
    source = _source_manifest(manifest_path)
    parents = _parent_receipts(dependency_receipts, plan, task, manifest_sha)
    if task.get("reuse_receipts"):
        raise ValueError("Cross-plan reuse receipts require an explicit verified import; reuse the existing same-plan task directory or Nextflow -resume instead of silently recomputing")
    base = {"plan_id": plan["plan_id"], "task": task, "source_project_manifest_sha256": manifest_sha,
            "package_version": __version__,
            "dependencies": {key: value[1]["receipt_sha256"] for key, value in sorted(parents.items())}}
    base_sha = canonical_hash(base)
    directory = Path(output_dir).expanduser().resolve()
    receipt_path = directory / "task.receipt.json"
    if receipt_path.exists():
        value = validate_receipt(receipt_path, plan=plan, task=task)
        if value.get("base_request_sha256") != base_sha:
            raise ValueError("Existing task products belong to changed inputs or execution parameters")
        return receipt_path
    for protected in (Path(__file__).resolve().parents[2], source.input_root, source.output_root):
        if directory == protected.resolve() or directory.is_relative_to(protected.resolve()):
            raise ValueError("An isolated flow task cannot write inside its source or existing result roots")
    directory.mkdir(parents=True, exist_ok=True)
    if any(directory.iterdir()):
        raise FileExistsError("Task directory already contains an unknown or failed attempt; preserve it and choose a new output directory")
    with (directory / "execution.claim.json").open("x", encoding="utf-8") as handle:
        json.dump({"pid": os.getpid(), "task_id": task_id, "base_request_sha256": base_sha}, handle)
    receipt = {"schema": RECEIPT_SCHEMA, "plan_id": plan["plan_id"], "task_id": task_id,
               "task_contract_sha256": canonical_hash(task), "workflow_id": task.get("workflow_id"),
               "stage": task.get("stage"), "source_project_manifest_sha256": manifest_sha,
               "source_project_manifest": str(manifest_path), "project_id": source.project_id,
               "base_request_sha256": base_sha, "package_version": __version__,
               "started_at": datetime.now(timezone.utc).isoformat(),
               "dependency_receipts": [{"task_id": key, **file_record(value[0])}
                                       for key, value in sorted(parents.items())]}
    try:
        from katana_next.domain.models import InputAsset
        derived_assets = []
        for parent_id, (_path, parent) in parents.items():
            asset = parent.get("execution_evidence", {}).get("project_asset")
            if asset:
                candidate = dict(asset)
                candidate.setdefault("asset_id", "flow-derived-" + parent_id)
                candidate.setdefault("required", True)
                candidate.setdefault("metadata", {})["flow_parent_receipt_sha256"] = parent["receipt_sha256"]
                derived_assets.append(InputAsset.model_validate(candidate))
        manifest = source.model_copy(update={"output_root": directory / "results",
                                             "assets": tuple(source.assets) + tuple(derived_assets)})
        original_import_inputs = []
        if task.get("workflow_id") == "prioritize.motiftrace_parallel":
            manifest, original_import_inputs, snapshot = _owned_parallel_inputs(manifest, task, directory)
            receipt["execution_evidence"] = {"owned_completed_parallel_snapshot": snapshot,
                                              "native_fitting_performed": False, "new_inference": False}
        private_manifest = write_json(directory / "project.json", manifest.model_dump(mode="json"))
        with KatanaProject(private_manifest) as project:
            _import_parents(project, parents, manifest_sha)
            _aggregate_parents(project, task)
            parameters = dict(task.get("parameters") or {})
            selection_evidence = _selection(project, task, parameters)
            receipt["selection"] = selection_evidence
            if task.get("stage") or task.get("adapter") or task.get("workflow_id", "").startswith("flow."):
                effective_task = json.loads(json.dumps(task))
                adapter = effective_task.get("adapter") or effective_task.get("workflow_id", "").removeprefix("flow.")
                if adapter.startswith("network.") or adapter == "qc.sequencing":
                    effective_task.setdefault("inputs", {})["project"] = str(private_manifest)
                for input_name, binding in effective_task.get("dependency_inputs", {}).items():
                    if not isinstance(binding, dict) or not binding.get("analysis_id"):
                        raise ValueError("Computed input bindings require an explicit producer workflow identity")
                    from katana_next.workflows.sources import resolve_table
                    requested = effective_task.get("scope", {}).get("contrast_ids", [])
                    if binding.get("kind") == "file":
                        candidates = [artifact for artifact in project.artifacts(binding["analysis_id"])
                                      if artifact.kind == "file" and artifact.name == binding.get("name")]
                        if len(candidates) != 1:
                            raise ValueError(f"Computed file input {input_name} needs one exact completed parent artifact")
                        source_artifact = candidates[0]
                    else:
                        source_artifact = resolve_table(project.context, binding["analysis_id"], binding.get("name"),
                                                        {"contrast_ids": requested} if requested else {})
                    if source_artifact is None:
                        raise ValueError(f"No verified completed parent table covers computed input {input_name}")
                    effective_task.setdefault("inputs", {})[input_name] = str(source_artifact.path)
                if adapter.startswith("network.") and adapter != "network.fit":
                    engines = [parent.get("execution_evidence", {}).get("engine_root")
                               for _path, parent in parents.values()
                               if parent.get("execution_evidence", {}).get("engine_root")]
                    if len(set(engines)) != 1:
                        raise ValueError("Network continuation needs one verified fit engine from its declared dependencies")
                    effective_task.setdefault("inputs", {})["engine_root"] = engines[0]
                    configurations = [parent.get("execution_evidence", {}).get("resolved_config")
                                      for _path, parent in parents.values()
                                      if parent.get("execution_evidence", {}).get("engine_root")]
                    if len(set(configurations)) != 1 or not configurations[0]:
                        raise ValueError("Network continuation needs its verified effective fit configuration")
                    effective_task["inputs"]["config"] = configurations[0]
                if adapter == "mapp.maei":
                    prepared = [parent.get("execution_evidence", {}).get("maei_inputs")
                                for _path, parent in parents.values()
                                if parent.get("execution_evidence", {}).get("maei_inputs")]
                    if prepared:
                        if len(prepared) != 1:
                            raise ValueError("Native MAPP fit needs one verified prepared quantity family")
                        effective_task.setdefault("inputs", {}).update(prepared[0])
                receipt["effective_task"] = effective_task
                paths = []
                def declared_paths(value):
                    if isinstance(value, dict):
                        for child in value.values():
                            declared_paths(child)
                    elif isinstance(value, list):
                        for child in value:
                            declared_paths(child)
                    elif isinstance(value, str) and Path(value).is_file():
                        paths.append(Path(value).resolve())
                declared_paths(effective_task.get("inputs", {}))
                if adapter.startswith("network.") and effective_task.get("inputs", {}).get("config"):
                    from katana_cosplice.config import load_config
                    from katana_next.workflows.programs import configuration_input_files
                    paths.extend(configuration_input_files(load_config(effective_task["inputs"]["config"])))
                receipt["inputs"] = [file_record(path) for path in dict.fromkeys(paths)]
                generated, stage = _stage_artifacts(project, effective_task, directory, canonical_hash(base))
                receipt["status"] = stage["status"]
                receipt["execution_evidence"] = stage.get("evidence", {})
            else:
                workflow = project.workflows.registry.get(task["workflow_id"])
                if task.get("workflow_version") and task["workflow_version"] != workflow.version:
                    raise ValueError("Workflow implementation version changed after planning")
                if task.get("selection") and not workflow.definition.supports_contrast_selection:
                    # Selection-only algorithms already obtain their exact
                    # comparison from the checked saved source rows. They do
                    # not expose the unrelated generic contrast-subset option.
                    declared = parameters.pop("contrast_ids", None)
                    selected_scopes = selection_evidence.get("contrast_ids", []) if selection_evidence else []
                    if declared and set(declared) != set(selected_scopes):
                        raise ValueError("Task comparison scope differs from its selected scientific source")
                normalized = workflow.normalize_parameters(parameters)
                if selection_evidence is not None and selection_evidence["selected_rows"] == 0:
                    # No target computation is scientifically applicable. Its
                    # audited foreground source is still a verified input;
                    # demanding a nonempty saved selection here would prevent
                    # the explicit zero-foreground result from being recorded.
                    paths = [project.context.catalog.get_artifact(selection_evidence["source_artifact_id"]).path]
                else:
                    paths = list(dict.fromkeys(Path(path).resolve() for path in workflow.input_paths(project.context, normalized)))
                receipt["inputs"] = list({record["path"]: record for record in
                                          [*original_import_inputs, *[file_record(path) for path in paths]]}.values())
                receipt["workflow_version"] = workflow.version
                if selection_evidence is not None and selection_evidence["selected_rows"] == 0:
                    from katana_next.workflows.base import ArtifactWriter
                    import pandas as pd
                    writer = ArtifactWriter(project.context, task["workflow_id"], canonical_hash(base), parameters, workflow.version)
                    audit = writer.table("No eligible foreground", "selection_audit.parquet",
                                         pd.DataFrame([selection_evidence]), entity_type="eligibility_audit",
                                         metadata={"not_applicable": True, "reason": "No events satisfy the declared foreground selection"})
                    writer.finish()
                    generated = [audit]
                    receipt["status"] = "not_applicable"
                    receipt["execution_evidence"] = {"reason": "No events satisfy the declared foreground selection", **selection_evidence}
                else:
                    from contextlib import nullcontext
                    manager = nullcontext()
                    if task["workflow_id"].startswith("programs."):
                        from .stages import persist_network_trees
                        manager = persist_network_trees()
                    elif task["workflow_id"] == "prioritize.motiftrace_consequences":
                        manager = _declare_complete_motif_family(project, parameters)
                    elif task["workflow_id"] == "prioritize.motiftrace_parallel":
                        manager = _explicit_parallel_scope(project, parameters)
                    with manager:
                        job = project.run(task["workflow_id"], parameters, timeout=None)
                    generated = [project.context.catalog.get_artifact(identifier) for identifier in job.artifact_ids]
                    if not generated or any(item is None for item in generated):
                        raise ValueError("The completed workflow has no registered products")
                    receipt["job"] = job.model_dump(mode="json")
                    receipt["status"] = "reused" if job.status == "cache_hit" else "completed"
                receipt["parameters"] = parameters
            verify_files(receipt["inputs"])
            all_artifacts = project.artifacts()
            receipt["artifacts"] = [a.model_dump(mode="json") for a in all_artifacts]
            receipt["generated_artifacts"] = [a.model_dump(mode="json") for a in generated]
            receipt["producer_markers"] = [file_record(path) for path in sorted((directory / "results").rglob("parameters.json"))]
            output_paths = set()
            imported = receipt.get("execution_evidence", {}).get("owned_completed_parallel_snapshot")
            if imported:
                output_paths.add(Path(imported["private_manifest"]).resolve())
            for artifact in all_artifacts:
                output_paths.add(artifact.path.resolve())
                for value in artifact.metadata.get("exports", {}).values():
                    if isinstance(value, str) and Path(value).is_file():
                        output_paths.add(Path(value).resolve())
            # Include auxiliary native products (FIMO/calibration details and
            # static figures), retaining the whole completed producer output.
            for root in (directory / "results" / "artifacts", directory / "compute"):
                if root.is_dir():
                    output_paths.update(path.resolve() for path in root.rglob("*")
                                        if path.is_file() and not path.is_symlink())
            receipt["outputs"] = [file_record(path) for path in sorted(output_paths)]
        if file_hash(manifest_path) != manifest_sha:
            raise ValueError("Source project manifest changed during task execution")
        receipt["finished_at"] = datetime.now(timezone.utc).isoformat()
        receipt["receipt_sha256"] = canonical_hash(receipt)
        write_json(receipt_path, receipt)
        validate_receipt(receipt_path, plan=plan, task=task)
        return receipt_path
    except BaseException as exc:
        receipt.update(status="failed", finished_at=datetime.now(timezone.utc).isoformat(),
                       error_type=type(exc).__name__, error=str(exc))
        receipt["receipt_sha256"] = canonical_hash({k: v for k, v in receipt.items() if k != "receipt_sha256"})
        write_json(directory / "failure.receipt.json", receipt)
        raise
