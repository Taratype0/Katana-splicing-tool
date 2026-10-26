"""Retain completed analysis inputs and products independently of Nextflow work.

Native tables, statistics and producer proofs are copied byte for byte.  A
separate transport catalogue binds their unchanged bytes to the retained paths;
it is not a new scientific fit or a claim that native completion markers moved.
"""
from __future__ import annotations

import copy
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import re
import shutil
from collections import Counter

from .collect import (materialize_report_descriptor, validate_collection,
                      _rewrite)
from .model import canonical_hash, file_sha256, read_json, write_json


SCHEMA = "katana.flow.retained_project/v1"


class _OwnedCopies:
    """Permit a move-by-link only from a registered, exclusively held stage.

    A path under work is insufficient: a site can stage an external user input
    by hardlink.  Account for every pre-existing alias before sharing its inode
    with retention, or fall back to an independent copy.
    """
    def __init__(self, scope, destination):
        from . import storage
        self.storage = storage
        self.scope, self.owner = storage._scope(scope)
        storage._assert_lock(self.scope, self.owner)
        if destination.resolve().is_relative_to(self.scope):
            raise ValueError("Retained project must be outside the managed cleanup scope")
        self.work = storage._bounded_path(self.scope / "work" / "downstream")
        self.registration = storage._read(self.scope / ".katana-stage-downstream.json")
        self._validate()
        rows = storage._inventory(self.work)
        self.original_aliases = Counter(tuple(row["node"]) for row in rows
                                        if not row["directory"] and not row["link"])
        self.retained_aliases = Counter()

    def _validate(self):
        storage = self.storage
        storage._assert_lock(self.scope, self.owner)
        registration = storage._read(self.scope / ".katana-stage-downstream.json")
        if (registration != self.registration or
                registration.get("scope_sha256") != self.owner["record_sha256"] or
                registration.get("stage") != "downstream" or
                registration.get("work_root") != str(self.work) or
                not self.work.is_dir() or storage._link(self.work.lstat()) or
                list(storage._node(self.work.lstat())) != registration.get("work_identity") or
                (self.scope / "cleanup-downstream.pending.json").exists() or
                (self.scope / "cleanup-downstream.receipt.json").exists()):
            raise ValueError("Owned retained source is not the live registered downstream stage")

    def allow(self, path):
        self._validate()
        resolved = path.resolve(strict=True)
        if not resolved.is_relative_to(self.work):
            return False
        # Reject any parent/file reparse escape at the actual link boundary.
        self.storage._bounded_path(resolved)
        info = resolved.lstat()
        if self.storage._link(info):
            return False
        node = self.storage._node(info)
        owned = self.original_aliases[node] + self.retained_aliases[node]
        return bool(self.original_aliases[node]) and info.st_nlink <= owned

    def record(self, destination):
        self.retained_aliases[self.storage._node(destination.lstat())] += 1


def _member(root, relative):
    relative = Path(relative)
    if relative.is_absolute() or ".." in relative.parts:
        raise ValueError("Retained project member must be a bounded relative path")
    result = root / relative
    if result.is_symlink() or not result.resolve().is_relative_to(root.resolve()):
        raise ValueError("Retained project member escapes its root")
    return result


def _record(root, path):
    return {"path": path.relative_to(root).as_posix(), "bytes": path.stat().st_size,
            "sha256": file_sha256(path)}


def _verified_copy(source, destination, copied_content=None, *, expected_sha256=None, owned_copies=None):
    """Deduplicate only new retained copies, never hardlink borrowed inputs."""
    before = (source.stat().st_size, source.stat().st_mtime_ns)
    expected = file_sha256(source)
    if expected_sha256 is not None and expected != expected_sha256:
        raise ValueError("A retained source changed after execution-bound input verification")
    if owned_copies is not None:
        owned_copies._validate()
    destination.parent.mkdir(parents=True, exist_ok=True)
    identity = (before[0], expected)
    prior = copied_content.get(identity) if copied_content is not None else None
    if prior is not None:
        if prior.is_symlink() or not prior.is_file() or file_sha256(prior) != expected:
            raise ValueError("An owned retained copy changed during content deduplication")
        try:
            os.link(prior, destination)
        except OSError:
            # Cross-filesystem or unsupported hardlinks never make borrowing
            # an input inode acceptable; fall back to an independent copy.
            shutil.copy2(source, destination)
    elif owned_copies is not None and owned_copies.allow(source):
        try:
            os.link(source, destination)
        except OSError:
            shutil.copy2(source, destination)
    else:
        shutil.copy2(source, destination)
    if ((source.stat().st_size, source.stat().st_mtime_ns) != before or
            file_sha256(source) != expected or file_sha256(destination) != expected):
        raise ValueError("A retained source changed during copying")
    if copied_content is not None:
        copied_content.setdefault(identity, destination)
    if owned_copies is not None:
        owned_copies.record(destination)


def _sidecars(path):
    candidates = [Path(str(path) + suffix) for suffix in (".fai", ".gzi", ".bai", ".crai")]
    if path.suffix.lower() == ".bam":
        candidates.append(path.with_suffix(".bai"))
    if path.suffix.lower() in {".fa", ".fasta", ".fna"}:
        candidates.append(path.with_suffix(".dict"))
    return [candidate for candidate in candidates if candidate.is_file()]


def _metadata_paths(metadata):
    """Only explicit supported resource overrides; never scan an environment."""
    settings = metadata.get("gene_filter", {})
    return [("gene_filter", key, settings[key]) for key in
            ("annotation_tsv", "gtf", "expression_matrix") if settings.get(key)]


def _registered_bundles(source, paths):
    """Retain only explicitly registered members of three existing formats."""
    kinds = {"rbp_sequence_results", "motiftrace_matched_windows", "motiftrace_parallel_results"}
    roots, metadata_paths, declared_hashes = [], [], {}
    for index, asset in enumerate(source.assets):
        supplied_root = asset.metadata.get("bundle_root")
        if asset.kind not in kinds:
            if supplied_root is not None:
                raise ValueError("Automatic retention does not support bundle_root for asset kind: " + asset.kind)
            continue
        manifest = asset.path.resolve(strict=True)
        if asset.kind != "motiftrace_parallel_results" and not supplied_root:
            raise ValueError("Registered completed bundle requires its explicit bundle_root")
        root = Path(supplied_root).expanduser() if supplied_root else manifest.parent
        root = (root if root.is_absolute() else source.input_root / root).resolve(strict=True)
        if not root.is_dir() or not manifest.is_relative_to(root):
            raise ValueError("Registered bundle manifest must be inside its declared root")
        if any(other != root and (root.is_relative_to(other) or other.is_relative_to(root)) for other in roots):
            raise ValueError("Overlapping completed bundle roots are not supported by automatic retention")
        roots.append(root)
        expected_manifest = asset.metadata.get("sha256")
        if not re.fullmatch(r"[0-9a-f]{64}", str(expected_manifest)) or file_sha256(manifest) != expected_manifest:
            raise ValueError("Registered bundle manifest SHA256 mismatch")
        value = read_json(manifest)
        if asset.kind == "rbp_sequence_results":
            if not isinstance(value, list):
                raise ValueError("Completed RBP member manifest must be a signed member list")
            members = value
        else:
            schema, key = (("motiftrace.matched_windows.v1", "members") if asset.kind == "motiftrace_matched_windows"
                           else ("motiftrace.completed_parallel.v1", "files"))
            if (not isinstance(value, dict) or value.get("schema") != schema or
                    value.get("state") != "complete" or not isinstance(value.get(key), dict)):
                raise ValueError("Completed bundle registration has the wrong schema or completion state")
            members = list(value[key].values())
        if not members:
            raise ValueError("Completed bundle has no registered members")
        prefix = Path("bundles") / hashlib.sha256(str(root).encode()).hexdigest()[:20]
        manifest_relative = (prefix / manifest.relative_to(root)).as_posix()
        prior = paths.get(manifest)
        if prior is not None and prior != manifest_relative:
            raise ValueError("Conflicting relative layout for one completed bundle member")
        paths[manifest] = manifest_relative
        declared_hashes[manifest] = expected_manifest
        seen = set()
        for member in members:
            if not isinstance(member, dict) or not isinstance(member.get("path"), str):
                raise ValueError("Completed bundle member requires a relative path and integrity record")
            name = member["path"]
            relative = PurePosixPath(name)
            if (not name or "\\" in name or ":" in name or relative.is_absolute() or ".." in relative.parts
                    or relative.as_posix() != name or name in seen or
                    not re.fullmatch(r"[0-9a-f]{64}", str(member.get("sha256"))) or
                    isinstance(member.get("bytes"), bool) or not isinstance(member.get("bytes"), int)
                    or member["bytes"] < 0):
                raise ValueError("Unsafe, duplicate or invalid registered bundle member")
            seen.add(name)
            path = root.joinpath(*relative.parts)
            resolved = path.resolve(strict=True)
            if path.is_symlink() or resolved != path or not resolved.is_relative_to(root) or not path.is_file():
                raise ValueError("Registered bundle member escapes its declared root")
            if path.stat().st_size != member["bytes"] or file_sha256(path) != member["sha256"]:
                raise ValueError("Registered bundle member SHA256/size mismatch: " + name)
            target = (prefix / Path(*relative.parts)).as_posix()
            prior = paths.get(resolved)
            if prior is not None and prior != target:
                raise ValueError("Conflicting relative layout for one completed bundle member")
            paths[resolved] = target
            prior_hash = declared_hashes.get(resolved)
            if prior_hash is not None and prior_hash != member["sha256"]:
                raise ValueError("Conflicting integrity for one completed bundle member")
            declared_hashes[resolved] = member["sha256"]
        if supplied_root is not None:
            metadata_paths.append([index, "bundle_root", (Path("inputs") / prefix).as_posix()])
    return metadata_paths, declared_hashes


def _project_payload(root):
    payload = read_json(root / "project-template.json")
    bindings = read_json(root / "resource-bindings.json")
    for group, key, relative in bindings["metadata_paths"]:
        payload["metadata"][group][key] = str(_member(root, relative).resolve())
    for index, key, relative in bindings.get("asset_metadata_paths", []):
        payload["assets"][index]["metadata"][key] = str(_member(root, relative).resolve())
    # This is a new filesystem identity, not a changed gene-filter policy.
    if payload["metadata"].get("gene_filter_identity"):
        from katana_next.infrastructure.manifest_store import ManifestStore
        from katana_next.domain.models import ProjectManifest
        from katana_next.science.gene_filter import gene_filter_details
        resolved = ManifestStore._resolve_paths(ProjectManifest.model_validate(payload), root)
        payload["metadata"]["gene_filter_identity"] = gene_filter_details(resolved)["identity"]
    return payload


def _bind_catalog(root):
    """Rebuild only path/snapshot-dependent catalogue and completion bindings."""
    from katana_next.domain.models import AnalysisArtifact
    from katana_next.infrastructure.fingerprints import fingerprint_file, fingerprint_value
    from katana_next.services.project_service import ProjectService

    payload = _project_payload(root)
    project_file = root / "project.json"
    if project_file.exists():
        current = read_json(project_file)
        if current != payload:
            # A relocated derived manifest can change explicit resource paths
            # and their path-sensitive filter identity, never a user's policy.
            def logical(value):
                value = copy.deepcopy(value)
                value.get("metadata", {}).pop("gene_filter_identity", None)
                for group, key, _relative in read_json(root / "resource-bindings.json")["metadata_paths"]:
                    value["metadata"][group].pop(key, None)
                for index, key, _relative in read_json(root / "resource-bindings.json").get("asset_metadata_paths", []):
                    value["assets"][index]["metadata"].pop(key, None)
                return value
            if logical(current) != logical(payload):
                raise ValueError("Retained project was edited; transport rebinding cannot overwrite user changes")
            write_json(project_file, payload)
    else:
        write_json(project_file, payload)
    context = ProjectService().open_project(root / "project.json")
    source_collection = root / "collection"
    validate_collection(source_collection)
    descriptor = materialize_report_descriptor(source_collection)
    source_catalog = read_json(source_collection / "artifacts.json")["artifacts"]
    originals = {artifact["artifact_id"]: artifact for artifact in source_catalog}
    native_markers = {}
    ledger = read_json(source_collection / "transport-provenance.json")["products"]
    for product in ledger:
        relative = product["collection_path"]
        if Path(relative).name != "parameters.json":
            continue
        path = _member(source_collection, relative)
        marker = read_json(path)
        if marker.get("completed") is True and marker.get("workflow_id") and marker.get("run_fingerprint"):
            key = (marker["workflow_id"], marker["run_fingerprint"])
            prior = native_markers.get(key)
            if prior and prior[0] != marker:
                raise ValueError("Conflicting original completion markers in retained collection")
            native_markers[key] = (marker, path)

    # Relative native product paths already resolve through the authenticated
    # collection descriptor; resource paths in descriptive metadata are rebound
    # separately. Original producer JSON files are never rewritten.
    resources = read_json(root / "resource-bindings.json")
    mapping = {old: str(_member(root, relative).resolve())
               for old, relative in resources["source_paths"].items()}
    artifacts = []
    for raw in descriptor["artifacts"]:
        raw = _rewrite(raw, mapping)
        artifact = AnalysisArtifact.model_validate(raw)
        if not artifact.path.is_file() or not artifact.path.resolve().is_relative_to(source_collection):
            raise ValueError("Retained catalogue artifact is not a verified collection member")
        original = originals[artifact.artifact_id]
        metadata = copy.deepcopy(artifact.metadata)
        metadata["storage_transport"] = {
            "schema": SCHEMA, "scientific_recalculation": False,
            "source_artifact_fingerprint": original["fingerprint"],
            "source_collection_receipt_sha256": read_json(source_collection / "collection.receipt.json")["receipt_sha256"],
        }
        metadata["project_snapshot_fingerprint"] = context.snapshot.manifest_fingerprint
        fingerprint = artifact.fingerprint
        if metadata.get("artifact_contract_version") == "3.0-correctness":
            key = (artifact.analysis_id, artifact.run_fingerprint)
            if key not in native_markers:
                raise ValueError("Retained artifact has no original completed producer proof")
            marker, marker_path = native_markers[key]
            matching = [item for item in marker.get("generated_artifacts", [])
                        if item.get("artifact_id") == artifact.artifact_id]
            if len(matching) != 1 or matching[0].get("fingerprint") != original["fingerprint"]:
                raise ValueError("Original completion marker does not authenticate the retained artifact")
            metadata["storage_transport"]["native_marker"] = marker_path.relative_to(root).as_posix()
            metadata["storage_transport"]["native_marker_sha256"] = file_sha256(marker_path)
            fingerprint = fingerprint_value({"run": artifact.run_fingerprint, "name": artifact.name,
                "path": artifact.path.name,
                "content": fingerprint_file(artifact.path, force_full=True).model_dump(mode="json")})
        updated = artifact.model_copy(update={"metadata": metadata, "fingerprint": fingerprint})
        context.catalog.register_artifact(updated)
        artifacts.append(updated)

    # Derived bindings live outside immutable scientific collection products.
    # They explicitly retain native marker identities, rather than pretending
    # that absolute-path fingerprints remain the same after relocation.
    for key, (native, native_path) in native_markers.items():
        members = [artifact for artifact in artifacts
                   if (artifact.analysis_id, artifact.run_fingerprint) == key]
        if not members:
            continue
        parameters = members[0].parameters
        if any(artifact.parameters != parameters for artifact in members):
            raise ValueError("One retained producer has conflicting scientific parameters")
        marker = copy.deepcopy(native)
        marker["parameters"] = parameters
        marker["generated_artifacts"] = [
            {"artifact_id": artifact.artifact_id, "name": artifact.name,
             "kind": artifact.kind.value, "path": str(artifact.path), "fingerprint": artifact.fingerprint}
            for artifact in members]
        marker["storage_transport"] = {"schema": SCHEMA, "scientific_recalculation": False,
            "native_marker": native_path.relative_to(root).as_posix(),
            "native_marker_sha256": file_sha256(native_path),
            "native_parameters_sha256": canonical_hash(native.get("parameters", {}))}
        write_json(context.results.artifact_directory(*key) / "parameters.json", marker)
    context.catalog.export_manifest(context.manifest.project_id)
    write_json(root / "catalog-binding.json", {
        "schema": "katana.flow.retained_catalog_binding/v1", "scientific_recalculation": False,
        "project_manifest_sha256": file_sha256(root / "project.json"),
        "project_snapshot_fingerprint": context.snapshot.manifest_fingerprint,
        "artifacts": [{"artifact_id": artifact.artifact_id, "fingerprint": artifact.fingerprint,
                       "source_fingerprint": originals[artifact.artifact_id]["fingerprint"]}
                      for artifact in artifacts]})
    return context


def retain_collection_project(collection, source_project, destination, *, owned_scope=None):
    """Create an independent analysis project; never modify or delete sources.

    Inputs are copied once per resolved source path (including explicit BAM
    indexes and existing FASTA indexes). No directory, software environment,
    unrelated run, or implicit raw-FASTQ tree is recursively imported. A managed
    caller may supply its currently held ``owned_scope`` to link confirmed
    stage-owned immutable bytes into retention before unlinking their work.
    Public callers omit it, and borrowed source inodes always remain isolated.
    """
    from katana_next.infrastructure.manifest_store import ManifestStore
    collection = Path(collection).expanduser().resolve()
    if collection.is_file():
        collection = collection.parent
    source_project = Path(source_project).expanduser().resolve()
    destination = Path(destination).expanduser().absolute()
    receipt = validate_collection(collection)
    if file_sha256(source_project) != receipt["source_project_manifest_sha256"]:
        raise ValueError("Retained project source differs from the collected scientific plan")
    _, source = ManifestStore().load(source_project)
    for protected in (collection, source.input_root, source.output_root, Path(__file__).resolve().parents[2]):
        if destination.resolve().is_relative_to(protected.resolve()):
            raise ValueError("Retained project must be outside source inputs, results, collection and software")
    if destination.exists() or destination.is_symlink():
        raise FileExistsError("Retained project destination already exists; preserve it")
    if not destination.parent.is_dir():
        raise ValueError("Retained project parent must already exist")
    owned_copies = _OwnedCopies(owned_scope, destination) if owned_scope is not None else None

    payload = source.model_dump(mode="json")
    paths = {}
    asset_metadata_paths, bundle_hashes = _registered_bundles(source, paths)
    missing_optional = []

    def retain(path, required=True):
        location = Path(path).expanduser().resolve()
        if location in paths:
            return paths[location]
        if location.exists() and not location.is_file():
            raise ValueError("Retained project accepts explicit files, not directory assets")
        parent_key = hashlib.sha256(str(location.parent).encode()).hexdigest()[:20]
        relative = (Path(parent_key) / location.name).as_posix()
        for existing, mapped in paths.items():
            if existing.parent == location.parent and Path(mapped).parts[0] == "bundles":
                relative = (Path(mapped).parent / location.name).as_posix()
                break
        if not location.is_file():
            if required:
                raise FileNotFoundError(f"Required retained input is missing: {location}")
            missing_optional.append({"path": str(location), "retained_path": "inputs/" + relative})
        paths[location] = relative
        return relative

    for asset in payload["assets"]:
        asset["path"] = retain(asset["path"], asset["required"])
    for sample in payload["samples"]:
        sample["assets"] = {key: retain(path) for key, path in sample["assets"].items()}
    for key in ("gtf", "genome_fasta", "transcript_fasta"):
        if payload["annotation"].get(key):
            payload["annotation"][key] = retain(payload["annotation"][key])
    metadata_paths = []
    for group, key, path in _metadata_paths(payload["metadata"]):
        location = Path(path)
        if not location.is_absolute():
            location = source.input_root / location
        relative = retain(location)
        metadata_paths.append([group, key, "inputs/" + relative])
    for path in list(paths):
        for sidecar in _sidecars(path):
            retain(sidecar)
    # A same-path edit after a calculation must not pair new resources with
    # old statistics. Preserve receipt bytes; compare explicit retained inputs
    # to their available original execution hashes before creating output.
    execution_inputs = {}
    for proof in receipt["task_receipts"]:
        value = read_json(_member(collection, proof["collection_path"]))
        for item in value.get("inputs", []):
            path = Path(item["path"]).resolve()
            prior = execution_inputs.get(path)
            if prior and (prior["sha256"], prior["size_bytes"]) != (item["sha256"], item["size_bytes"]):
                raise ValueError("Original task proofs contain conflicting input identities")
            execution_inputs[path] = item
    declared_hashes = {asset.path.resolve(): asset.metadata.get("upstream_sha256") or asset.metadata.get("sha256")
                       for asset in source.assets}
    declared_hashes.update(bundle_hashes)
    input_proofs = []
    verified_input_hashes = {}
    for path, relative in paths.items():
        if not path.is_file():
            continue
        expected = execution_inputs.get(path)
        digest = file_sha256(path)
        if expected and (path.stat().st_size != expected["size_bytes"] or digest != expected["sha256"]):
            raise ValueError("Retained input changed after scientific execution: " + str(path))
        if declared_hashes.get(path) and digest != declared_hashes[path]:
            raise ValueError("Retained input differs from its declared source SHA256: " + str(path))
        verified_input_hashes[path] = digest
        input_proofs.append({"retained_path": "inputs/" + relative, "sha256": digest,
            "verification": "original_task_input" if expected else
                            "declared_input_sha256" if declared_hashes.get(path) else "current_retained_copy"})
    payload["input_root"], payload["output_root"] = "inputs", "results"
    payload["metadata"]["storage_retention"] = {
        "schema": SCHEMA, "source_project_manifest_sha256": receipt["source_project_manifest_sha256"],
        "source_collection_receipt_sha256": receipt["receipt_sha256"], "scientific_recalculation": False,
    }
    destination.mkdir()
    write_json(destination / "retention.claim.json", {"schema": SCHEMA,
        "source_project_manifest_sha256": receipt["source_project_manifest_sha256"]})
    copied = []
    copied_content = {}
    for path, relative in sorted(paths.items(), key=lambda item: str(item[0])):
        if path.is_file():
            target = destination / "inputs" / relative
            _verified_copy(path, target, copied_content, expected_sha256=verified_input_hashes[path],
                           owned_copies=owned_copies)
            copied.append(target)
    # Every producer byte and proof remains unchanged in its own collection.
    collection_hashes = {item.get("relative_path", item["path"]): item["sha256"] for item in receipt["outputs"]}
    for path in sorted(collection.rglob("*")):
        if path.is_symlink() or not path.resolve().is_relative_to(collection):
            raise ValueError("Retained collection cannot contain symbolic links")
        if path.is_file():
            target = destination / "collection" / path.relative_to(collection)
            _verified_copy(path, target, copied_content,
                           expected_sha256=collection_hashes.get(path.relative_to(collection).as_posix()),
                           owned_copies=owned_copies)
            copied.append(target)
    template = write_json(destination / "project-template.json", payload)
    bindings = write_json(destination / "resource-bindings.json", {
        "source_paths": {str(path): "inputs/" + relative for path, relative in paths.items()},
        "metadata_paths": metadata_paths, "missing_optional_inputs": missing_optional,
        "asset_metadata_paths": asset_metadata_paths, "input_proofs": input_proofs})
    original = destination / "source-project.json"
    _verified_copy(source_project, original, copied_content,
                   expected_sha256=receipt["source_project_manifest_sha256"], owned_copies=owned_copies)
    copied.extend((template, bindings, original))
    data = {"schema": SCHEMA, "status": "completed", "scientific_recalculation": False,
            "source_project_manifest_sha256": receipt["source_project_manifest_sha256"],
            "source_collection_receipt_sha256": receipt["receipt_sha256"],
            "outputs": [_record(destination, path) for path in copied],
            "project_manifest": "project.json", "artifact_count": receipt["artifact_count"]}
    data["receipt_sha256"] = canonical_hash(data)
    write_json(destination / "retention.receipt.json", data)
    try:
        context = validate_retained_project(destination)
    except Exception:
        # An interrupted/failed retention attempt must never pass the cleanup
        # gate. Preserve its evidence and withhold the successful receipt.
        (destination / "retention.receipt.json").rename(destination / "retention.failed.json")
        raise
    return {"schema": SCHEMA, "status": "completed", "project_manifest": str(context.manifest_path),
            "receipt": str(destination / "retention.receipt.json"), "artifact_count": data["artifact_count"]}


def validate_retained_project(destination):
    """Revalidate bytes and bind the retained project at its current location.

    The derived catalogue may be recreated after moving to another filesystem.
    All original scientific products and their signed collection stay fixed.
    """
    root = Path(destination).expanduser().resolve()
    receipt = read_json(root / "retention.receipt.json")
    unsigned = {key: value for key, value in receipt.items() if key != "receipt_sha256"}
    if (receipt.get("schema") != SCHEMA or receipt.get("status") != "completed" or
            receipt.get("receipt_sha256") != canonical_hash(unsigned)):
        raise ValueError("No valid completed retained project receipt")
    seen = set()
    for item in receipt.get("outputs", []):
        if item["path"] in seen:
            raise ValueError("Duplicate retained project ledger member")
        seen.add(item["path"])
        path = _member(root, item["path"])
        if (not path.is_file() or path.stat().st_size != item["bytes"] or
                file_sha256(path) != item["sha256"]):
            raise ValueError("Retained project member is missing or changed: " + item["path"])
    mandatory = {"project-template.json", "resource-bindings.json", "source-project.json",
                 "collection/collection.receipt.json", "collection/artifacts.json"}
    if not mandatory.issubset(seen):
        raise ValueError("Retained project proof is incomplete")
    for proof in read_json(root / "resource-bindings.json")["input_proofs"]:
        if proof["retained_path"] not in seen or file_sha256(_member(root, proof["retained_path"])) != proof["sha256"]:
            raise ValueError("Retained resource no longer matches its execution/copy input proof")
    collection_receipt = validate_collection(root / "collection")
    if (collection_receipt["source_project_manifest_sha256"] != receipt["source_project_manifest_sha256"] or
            collection_receipt["receipt_sha256"] != receipt["source_collection_receipt_sha256"] or
            file_sha256(root / "source-project.json") != receipt["source_project_manifest_sha256"]):
        raise ValueError("Retained collection and original project identities differ")
    return _bind_catalog(root)


__all__ = ["retain_collection_project", "validate_retained_project"]
