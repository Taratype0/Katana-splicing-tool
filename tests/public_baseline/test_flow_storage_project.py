"""Retained results remain usable after deleting their generated source work."""
from __future__ import annotations

import json
import os
from pathlib import Path
import shutil

import pandas as pd
import pytest

from katana_next.flow.collect import collect_results, validate_collection
from katana_next.flow.model import file_sha256, seal_plan, write_json
from katana_next.flow.runtime import execute_task
from katana_next.flow.storage_project import retain_collection_project, validate_retained_project
from katana_next.full_fixture import create_full_fixture
from katana_next.services.project_service import ProjectService
from katana_next.targets.bundle import export_targets, open_bundle


def _collect_fixture(work, *, configure=None):
    work.mkdir(exist_ok=True)
    project = create_full_fixture(work / "project")
    payload = json.loads(project.read_text())
    from katana_next.sashimi_fixture import write_indexed_bam
    # Tiny real indexed BAMs exercise the existing native Sashimi reader.
    for sample in payload["samples"]:
        bam = Path(payload["input_root"]) / (sample["sample_id"] + ".bam")
        write_indexed_bam(bam)
        sample["assets"]["bam"] = str(bam)
        sample["assets"]["bai"] = str(bam) + ".bai"
    fasta = Path(payload["annotation"]["genome_fasta"])
    lines = fasta.read_bytes().splitlines(keepends=True)
    assert len(lines) == 2 and lines[0].startswith(b">")
    Path(str(fasta) + ".fai").write_text(
        f"{lines[0][1:].strip().decode()}\t{len(lines[1].strip())}\t{len(lines[0])}\t{len(lines[1].strip())}\t{len(lines[1])}\n")
    # Place one generated SE event inside the BAM fixture's 1-kb reference.
    # This changes only software fixture geometry before any analysis executes.
    for asset in payload["assets"]:
        path = Path(asset["path"])
        if path.name.startswith("SE.MATS."):
            frame = pd.read_csv(path, sep="\t")
            for key, value in {"upstreamES": 100, "upstreamEE": 120,
                               "exonStart_0base": 220, "exonEnd": 240,
                               "downstreamES": 320, "downstreamEE": 340}.items():
                frame.loc[0, key] = value
            frame.to_csv(path, sep="\t", index=False)
    if configure is not None:
        configure(payload, work)
    write_json(project, payload)
    plan = write_json(work / "plan.json", seal_plan({
        "project_manifest": str(project), "project_manifest_sha256": file_sha256(project),
        "tasks": [{"id": "as", "workflow_id": "differential.as", "dependencies": [],
                   "parameters": {}, "scope": {}}]}))
    result = execute_task(plan, "as", work / "as")
    collection_receipt = collect_results(plan, [result], work / "collection")
    return work, project, collection_receipt.parent


@pytest.fixture
def collected(tmp_path):
    return _collect_fixture(tmp_path / "owned-generated-work")


def test_source_work_removed_then_client_gene_event_export_and_cards(collected, tmp_path):
    work, source, collection = collected
    native = {path.relative_to(collection).as_posix(): file_sha256(path)
              for path in collection.rglob("*") if path.is_file()}
    retained = tmp_path / "retained"
    result = retain_collection_project(collection, source, retained)
    assert result["status"] == "completed"
    before = validate_retained_project(retained)
    artifact = next(a for a in before.catalog.list_artifacts(before.manifest.project_id)
                    if a.kind == "table" and a.name == "Differential AS events")
    frame = pd.read_parquet(artifact.path)
    row = frame[frame.event_type.eq("SE") & frame.exonStart_0base.eq(220)].iloc[0]
    # Recursive deletion is bounded to this test's own generated workspace.
    assert work.resolve().is_relative_to(tmp_path.resolve())
    shutil.rmtree(work)
    context = ProjectService().open_project(result["project_manifest"])
    assert all(path.is_file() for sample in context.manifest.samples for path in sample.assets.values())
    assert len(list((retained / "inputs").rglob("*.bam"))) == len(context.manifest.samples)
    bam_paths = [sample.assets["bam"] for sample in context.manifest.samples]
    # All generated BAMs have identical bytes. Deduplication stays entirely
    # within the new retention tree, rather than borrowing mutable source data.
    probe = tmp_path / "hardlink-probe"
    try:
        os.link(bam_paths[0], probe)
    except OSError:
        pass
    else:
        probe.unlink()
        assert os.path.samefile(bam_paths[0], bam_paths[1])
    assert Path(str(context.manifest.annotation.genome_fasta) + ".fai").is_file()
    for relative, digest in native.items():
        assert file_sha256(retained / "collection" / relative) == digest
    validate_collection(retained / "collection")
    for kind in ("gene", "event"):
        modules = ["gene_card", "event_card", "sashimi"] if kind == "event" else ["gene_card", "event_card"]
        exported = export_targets(context, [{"source_artifact_id": artifact.artifact_id,
            "row_ids": [str(row.row_id)], "target_kind": kind, "filters": {}}],
            tmp_path / (kind + "-request"), modules=modules)
        with open_bundle(exported["request_path"]) as bundle:
            assert bundle.request["targets"][0]["gene_id"] == str(row.gene_id)
            assert bundle.request["resources"]
        from katana_next.targets.runner import run_targets, inspect_results
        target_output = tmp_path / (kind + "-cards")
        ran = run_targets(exported["request_path"], target_output,
                          exported["local_bindings_template"], workers=1)
        assert inspect_results(target_output)["status"] == "complete"
        if kind == "event":
            assert ran["targets"][0]["sashimi"][0]["status"] == "complete"
            assert ran["targets"][0]["sashimi"][0]["reader"] == "bamnostic"
        assert ran


def test_retained_copy_can_rebind_at_new_location_without_original_work(collected, tmp_path):
    work, source, collection = collected
    retained = tmp_path / "first-retained"
    retain_collection_project(collection, source, retained)
    moved = tmp_path / "other-machine"
    shutil.copytree(retained, moved)
    assert work.resolve().is_relative_to(tmp_path.resolve())
    shutil.rmtree(work)
    assert retained.resolve().is_relative_to(tmp_path.resolve())
    shutil.rmtree(retained)
    context = validate_retained_project(moved)
    assert context.manifest_path == moved / "project.json"
    assert all(a.path.is_relative_to(moved) for a in context.catalog.list_artifacts(context.manifest.project_id))
    artifact = next(a for a in context.catalog.list_artifacts(context.manifest.project_id)
                    if a.kind == "table" and a.name == "Differential AS events")
    row = pd.read_parquet(artifact.path).iloc[0]
    result = export_targets(context, [{"source_artifact_id": artifact.artifact_id,
        "row_ids": [str(row.row_id)], "target_kind": "gene"}], tmp_path / "moved-request",
        modules=["gene_card"])
    assert Path(result["request_path"]).is_file()


def test_retention_rejects_wrong_source_existing_destination_and_corrupt_retention(collected, tmp_path):
    _, source, collection = collected
    wrong = tmp_path / "wrong.json"
    wrong.write_text(source.read_text() + "\n")
    with pytest.raises(ValueError, match="differs"):
        retain_collection_project(collection, wrong, tmp_path / "bad-source")
    occupied = tmp_path / "occupied"
    occupied.mkdir()
    with pytest.raises(FileExistsError):
        retain_collection_project(collection, source, occupied)
    retained = tmp_path / "retained"
    retain_collection_project(collection, source, retained)
    data = json.loads((retained / "retention.receipt.json").read_text())
    member = retained / next(item["path"] for item in data["outputs"] if item["path"].startswith("inputs/"))
    member.write_bytes(b"corrupt")
    with pytest.raises(ValueError, match="missing or changed"):
        validate_retained_project(retained)


def test_retention_cannot_live_in_owned_sources(collected):
    _, source, collection = collected
    with pytest.raises(ValueError, match="outside"):
        retain_collection_project(collection, source, collection / "unsafe")


def test_revalidation_does_not_overwrite_user_project_changes(collected, tmp_path):
    _, source, collection = collected
    retained = tmp_path / "retained"
    retain_collection_project(collection, source, retained)
    manifest = retained / "project.json"
    payload = json.loads(manifest.read_text())
    payload["metadata"]["user_note"] = "User's later project edit"
    write_json(manifest, payload)
    edited = manifest.read_bytes()
    with pytest.raises(ValueError, match="cannot overwrite user changes"):
        validate_retained_project(retained)
    assert manifest.read_bytes() == edited


def test_deduplication_never_borrows_source_inodes(collected, tmp_path):
    _, source, collection = collected
    source_payload = json.loads(source.read_text())
    source_bams = [Path(sample["assets"]["bam"]) for sample in source_payload["samples"]]
    retained = tmp_path / "retained"
    retain_collection_project(collection, source, retained)
    context = validate_retained_project(retained)
    retained_bams = [sample.assets["bam"] for sample in context.manifest.samples]
    assert source_bams[0].read_bytes() == source_bams[1].read_bytes()
    assert not os.path.samefile(source_bams[0], retained_bams[0])
    # Source edits are isolated from every retained equivalent copy.
    original = retained_bams[0].read_bytes()
    source_bams[0].write_bytes(b"source was later edited")
    assert all(path.read_bytes() == original for path in retained_bams)
    validate_retained_project(retained)


def test_retention_deduplication_falls_back_without_hardlink_support(tmp_path, monkeypatch):
    from katana_next.flow import storage_project
    sources = tmp_path / "borrowed-inputs"
    targets = tmp_path / "retained-inputs"
    sources.mkdir()
    first, second = sources / "one", sources / "two"
    first.write_bytes(b"identical verified payload")
    second.write_bytes(first.read_bytes())
    memo = {}
    storage_project._verified_copy(first, targets / "one", memo)
    def unavailable(*_args, **_kwargs):
        raise OSError("filesystem does not support hardlinks")
    monkeypatch.setattr(storage_project.os, "link", unavailable)
    storage_project._verified_copy(second, targets / "two", memo)
    assert (targets / "one").read_bytes() == (targets / "two").read_bytes()
    assert not os.path.samefile(targets / "one", targets / "two")
    assert not os.path.samefile(first, targets / "one")


def test_same_path_input_edit_after_scientific_execution_blocks_retention(collected, tmp_path):
    _, source, collection = collected
    payload = json.loads(source.read_text())
    resource = next(Path(asset["path"]) for asset in payload["assets"] if asset["kind"] == "rmats_result")
    resource.write_bytes(resource.read_bytes() + b"\n")
    destination = tmp_path / "rejected"
    with pytest.raises(ValueError, match="changed after scientific execution"):
        retain_collection_project(collection, source, destination)
    assert not destination.exists()


def test_managed_owned_link_handoff_then_actual_cleanup_keeps_sashimi(tmp_path):
    from katana_next.flow.storage import create_managed_scope, managed_run_lock, register_stage, cleanup_stage
    scope, final = tmp_path / "managed-run", tmp_path / "final"
    create_managed_scope(scope, protected_roots=[final])
    with managed_run_lock(scope):
        register_stage(scope, "downstream")
        work, source, collection = _collect_fixture(scope / "work" / "downstream")
        first_original_bam = Path(json.loads(source.read_text())["samples"][0]["assets"]["bam"])
        retain_collection_project(collection, source, final, owned_scope=scope)
        context = validate_retained_project(final)
        retained_bam = context.manifest.samples[0].assets["bam"]
        probe = tmp_path / "link-probe"
        try:
            os.link(first_original_bam, probe)
        except OSError:
            pass
        else:
            probe.unlink()
            assert os.path.samefile(first_original_bam, retained_bam)
        def accepted(_):
            validate_retained_project(final)
            return {"status": "completed"}
        cleaned = cleanup_stage(scope, "downstream", native_status={"engine": "nextflow", "state": "COMPLETED",
            "exit_code": 0, "engine_stopped": True}, checkpoint_validator=accepted)
        assert cleaned["status"] == "completed" and not work.exists()
    context = validate_retained_project(final)
    artifact = next(a for a in context.catalog.list_artifacts(context.manifest.project_id)
                    if a.kind == "table" and a.name == "Differential AS events")
    frame = pd.read_parquet(artifact.path)
    row = frame[frame.event_type.eq("SE") & frame.exonStart_0base.eq(220)].iloc[0]
    exported = export_targets(context, [{"source_artifact_id": artifact.artifact_id,
        "row_ids": [str(row.row_id)], "target_kind": "event"}], tmp_path / "owned-request",
        modules=["event_card", "sashimi"])
    from katana_next.targets.runner import run_targets
    ran = run_targets(exported["request_path"], tmp_path / "owned-cards",
                      exported["local_bindings_template"], workers=1)
    assert ran["status"] == "complete"
    assert ran["targets"][0]["sashimi"][0]["status"] == "complete"


def test_owned_link_policy_does_not_share_external_hardlinked_original(tmp_path):
    from katana_next.flow import storage_project
    from katana_next.flow.storage import create_managed_scope, managed_run_lock, register_stage
    scope, final = tmp_path / "managed-run", tmp_path / "final"
    outside = tmp_path / "user-original"
    outside.write_bytes(b"external original")
    create_managed_scope(scope, protected_roots=[outside, final])
    with managed_run_lock(scope):
        register_stage(scope, "downstream")
        staged = scope / "work" / "downstream" / "staged-input"
        try:
            os.link(outside, staged)
        except OSError:
            pytest.skip("Host filesystem does not support hardlinks")
        policy = storage_project._OwnedCopies(scope, final)
        destination = final / "file"
        storage_project._verified_copy(staged, destination, owned_copies=policy)
        assert not os.path.samefile(outside, destination)
        outside.write_bytes(b"externally changed")
        assert destination.read_bytes() == b"external original"


def test_owned_scope_requires_current_process_run_lock(collected, tmp_path):
    from katana_next.flow.storage import create_managed_scope, managed_run_lock, register_stage
    _, source, collection = collected
    scope = tmp_path / "scope"
    create_managed_scope(scope)
    with managed_run_lock(scope):
        register_stage(scope, "downstream")
    final = tmp_path / "final"
    with pytest.raises(RuntimeError, match="held exclusive"):
        retain_collection_project(collection, source, final, owned_scope=scope)
    assert not final.exists()


def _bundle_asset(payload, work, kind, *, root=None, name="registered.json", problem=None):
    root = root or work / ("signed-" + kind)
    root.mkdir(parents=True, exist_ok=True)
    member = root / "nested" / "values.tsv"
    member.parent.mkdir(exist_ok=True)
    member.write_bytes(b"model\tvalue\nM1\t0.03125\n")
    record = {"path": "nested/values.tsv", "bytes": member.stat().st_size, "sha256": file_sha256(member)}
    if problem == "escape":
        record["path"] = "../outside.tsv"
    manifest = root / name
    if kind == "rbp_sequence_results" or problem == "unknown":
        contents = [record]
    else:
        schema, key = (("motiftrace.matched_windows.v1", "members") if kind == "motiftrace_matched_windows"
                       else ("motiftrace.completed_parallel.v1", "files"))
        contents = {"schema": schema, "state": "complete", key: {"values": record}}
    write_json(manifest, contents)
    metadata = {"sha256": file_sha256(manifest)}
    if kind != "motiftrace_parallel_results":
        metadata["bundle_root"] = str(root)
    payload["assets"].append({"asset_id": "bundle-" + kind + "-" + name, "kind": kind,
                              "path": str(manifest), "required": True, "metadata": metadata})
    if problem == "changed":
        member.write_bytes(b"changed after signed registration")
    elif problem == "missing":
        member.unlink()
    return root, manifest, record


def test_registered_resource_bundles_keep_signed_member_layout_after_work_deletion_and_move(tmp_path):
    kinds = ("rbp_sequence_results", "motiftrace_matched_windows", "motiftrace_parallel_results")
    native = {}
    def configure(payload, work):
        for kind in kinds:
            _root, manifest, record = _bundle_asset(payload, work, kind)
            native[kind] = manifest.read_bytes(), record
    work, source, collection = _collect_fixture(tmp_path / "generated-work", configure=configure)
    retained = tmp_path / "retained"
    retain_collection_project(collection, source, retained)
    assert work.resolve().is_relative_to(tmp_path.resolve())
    shutil.rmtree(work)
    moved = tmp_path / "moved-retained"
    retained.rename(moved)
    context = validate_retained_project(moved)
    for asset in context.manifest.assets:
        if asset.kind not in kinds:
            continue
        original, record = native[asset.kind]
        assert asset.path.read_bytes() == original
        root = Path(asset.metadata["bundle_root"]) if "bundle_root" in asset.metadata else asset.path.parent
        assert root.is_relative_to(moved)
        assert file_sha256(root / record["path"]) == record["sha256"]
        assert (root / record["path"]).stat().st_size == record["bytes"]


@pytest.mark.parametrize("problem", ["escape", "changed", "missing", "unknown", "overlap"])
def test_invalid_or_unsupported_resource_bundle_refuses_before_creating_retention(tmp_path, problem):
    def configure(payload, work):
        kind = "unknown_result_bundle" if problem == "unknown" else "rbp_sequence_results"
        root, _manifest, _record = _bundle_asset(payload, work, kind, problem=problem)
        if problem == "overlap":
            _bundle_asset(payload, work, kind, root=root / "nested", name="second.json")
    _work, source, collection = _collect_fixture(tmp_path / "generated-work", configure=configure)
    destination = tmp_path / "refused"
    with pytest.raises((ValueError, OSError)):
        retain_collection_project(collection, source, destination)
    assert not destination.exists()
