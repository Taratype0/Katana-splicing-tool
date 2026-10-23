"""Read-only views from synthetic, content-bound completed-stage fixtures.

These tests execute the real view/association/ORA kernels, not native network
discovery. They verify transport, lineage and no second statistical fitting.
"""
import json
from pathlib import Path
import shutil

import numpy as np
import pandas as pd
import pytest

from katana_next.flow.model import seal_plan, write_json
from katana_next.flow.runtime import canonical_hash, file_record, file_hash
from katana_next.flow.presentation import build_presentation_assets, validate_network_view_proof
from test_report_network_views import _data
from test_report_completed_motif import completed


def _fixture_collection(root, project, tasks, products):
    """Make explicit synthetic lifecycle receipts; never claim native E2E."""
    root.mkdir()
    write_json(root / "project.json", project)
    plan = seal_plan({"project_manifest": str(root / "project.json"),
                      "project_manifest_sha256": file_hash(root / "project.json"),
                      "contrast_ids": [row["contrast_id"] for row in project["contrasts"]], "tasks": tasks})
    write_json(root / "plan.json", plan)
    receipts, transport = {}, []
    for task in plan["tasks"]:
        original_paths, evidence = products[task["id"]]
        outputs = [file_record(path) for path in original_paths]
        receipt = {"schema": "katana.flow.task_receipt/v1", "status": "completed",
                   "task_id": task["id"], "plan_id": plan["plan_id"],
                   "task_contract_sha256": canonical_hash(task),
                   "source_project_manifest_sha256": plan["project_manifest_sha256"],
                   "project_id": project["project_id"], "outputs": outputs,
                   "generated_artifacts": [{"artifact_id": "synthetic-view-input-" + task["id"]}],
                   "execution_evidence": evidence,
                   "job": {"status": "complete"}, "synthetic_lifecycle_fixture": True}
        receipt["receipt_sha256"] = canonical_hash(receipt)
        path = root / "evidence" / task["id"] / "task.receipt.json"
        path.parent.mkdir(parents=True)
        write_json(path, receipt)
        receipts[task["id"]] = path, receipt
        common = Path(original_paths[0]).parent
        while any(not Path(p).is_relative_to(common) for p in original_paths):
            common = common.parent
        for output in outputs:
            target = root / "products" / task["id"] / Path(output["path"]).relative_to(common)
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(output["path"], target)
            transport.append({"source": output, "collection_path": target.relative_to(root).as_posix()})
    return plan, receipts, transport


def _network_collection(tmp_path):
    from scipy.cluster.hierarchy import linkage
    from katana_next.reports.network_views import matched_module_associations
    from katana_next.workflows.module_functions import enrich_modules
    tmp_path.mkdir(exist_ok=True)
    source, samples, scores, contrasts = _data(tmp_path)
    fit = source / "fit-engine"
    native = fit / "04_network/observed"
    matrix = fit / "02_matrix"
    native.mkdir(parents=True); matrix.mkdir()
    for original, target in (("assignments.tsv", native / "event_module_assignments.tsv"),
                             ("scores.tsv", native / "module_program_scores.tsv"),
                             ("psi.tsv", native / "selected_event_matrix.tsv"),
                             ("annotation.tsv", matrix / "event_metadata.tsv"),
                             ("samples.tsv", matrix / "sample_metadata_resolved.tsv")):
        shutil.copyfile(source / original, target)
    pd.DataFrame({"event_a": ["e1"], "event_b": ["e2"], "adjacency": [.8]}).to_csv(native / "network_edges.tsv", sep="\t", index=False)
    np.save(native / "event_linkage.npy", linkage(pd.read_csv(source / "psi.tsv", sep="\t", index_col=0).values))
    boot = source / "bootstrap-engine"
    destination = boot / "04c_consensus_modules"
    destination.mkdir(parents=True)
    shutil.copyfile(source / "assignments.tsv", destination / "consensus_event_module_assignments.tsv")
    pd.DataFrame({"module": ["M01", "M02"], "confidence": [.9, .8]}).to_csv(destination / "consensus_module_summary.tsv", sep="\t", index=False)
    associate = source / "associate-engine"
    associate.mkdir()
    full, selected, audit = matched_module_associations(scores, samples, contrasts)
    full.to_csv(associate / "full.tsv", sep="\t", index=False)
    selected.to_csv(associate / "selected.tsv", sep="\t", index=False)
    write_json(associate / "audit.json", audit)
    annotation = source / "annotation-engine"
    annotation.mkdir()
    assigned = pd.read_csv(source / "assignments.tsv", sep="\t")
    annotated = pd.read_csv(source / "annotation.tsv", sep="\t")
    tested, background = enrich_modules(assigned, annotated, {"fixed_mouse": source / "sets.gmt"})
    tested.to_csv(annotation / "tested.tsv.gz", sep="\t", index=False)
    background.to_csv(annotation / "background.tsv", sep="\t", index=False)
    assigned.merge(annotated, on="event_key")[["module", "gene_symbol"]].drop_duplicates().to_csv(annotation / "genes.tsv", sep="\t", index=False)
    tasks = [{"id": "fit", "workflow_id": "flow.network.fit", "adapter": "network.fit", "dependencies": []}]
    tasks += [{"id": name, "workflow_id": "flow.network." + name, "adapter": "network." + name,
               "dependencies": ["fit"], "scope": {"contrast_ids": ["B_vs_A"]}}
              for name in ("bootstrap", "associate", "annotate")]
    products = {
        "fit": ([p for p in fit.rglob("*") if p.is_file()], {"engine_root": str(fit), "executed_stages": ["network_observed"]}),
        "bootstrap": ([p for p in boot.rglob("*") if p.is_file()], {"engine_root": str(boot), "source_engine_root": str(fit), "executed_stages": ["consensus_modules"]}),
        "associate": (list(associate.iterdir()), {"source_engine_root": str(fit), "executed_stages": ["associate_observed_declared_contrasts"],
                       "precomputed_views": {"full_tests": str(associate / "full.tsv"), "selected_tests": str(associate / "selected.tsv"), "trait_audit": str(associate / "audit.json")}}),
        "annotate": (list(annotation.iterdir()), {"source_engine_root": str(fit), "executed_stages": ["annotate_modules"],
                      "precomputed_views": {"tested": str(annotation / "tested.tsv.gz"), "background": str(annotation / "background.tsv"), "input_genes": str(annotation / "genes.tsv")}}),
    }
    return (*_fixture_collection(tmp_path / "collection", {"project_id": "study", "contrasts": contrasts}, tasks, products), source)


def test_generic_network_complete_views_never_repeat_completed_stats(tmp_path, monkeypatch):
    plan, receipts, transport, source = _network_collection(tmp_path)
    before = {str(path): file_hash(path) for path in source.rglob("*") if path.is_file()}
    def forbidden(*args, **kwargs):
        raise AssertionError("Report views must not repeat completed model or ORA fits")
    monkeypatch.setattr("katana_next.reports.network_views.matched_module_associations", forbidden)
    monkeypatch.setattr("katana_next.workflows.module_functions.enrich_modules", forbidden)
    root = tmp_path / "collection"
    result = build_presentation_assets(root, plan, receipts, transport)
    assert len(result["extra_assets"]) >= 10
    audit = json.loads((root / "presentation/network/network-report.receipt.json").read_text())
    assert audit["bootstrap_accepted"] is True
    assert audit["network_discovery_runs"] == audit["new_trait_fits"] == 0
    assert audit["new_module_ora"] is False
    assert not result["extra_views"][0]["unavailable"]
    assert "e5" in (root / "presentation/network/browser/ALL_MODULES.html").read_text()
    assert (root / "presentation/network/tables/matched_module_tests_full_family.tsv").read_bytes() == (source / "associate-engine/full.tsv").read_bytes()
    for item in result["extra_assets"]:
        if item["path"].endswith(".png"):
            assert any(link["role"] == "plot_data" for link in item["metadata"]["data_links"])
    assert before == {str(path): file_hash(path) for path in source.rglob("*") if path.is_file()}
    proof = json.loads((root / "presentation/network-inputs.json").read_text())
    assert "actual_sacct" not in proof
    changed = root / proof["files"][0]["path"]
    changed.write_bytes(b"corrupt")
    with pytest.raises(ValueError, match="member changed"):
        validate_network_view_proof(root / "presentation/network-inputs.json", [], proof["contrasts"])


def test_generic_motif_views_verify_complete_family_and_keep_all_models(completed, tmp_path):
    source, _old_proof = completed
    tasks = [{"id": "motif", "workflow_id": "prioritize.motiftrace_consequences", "dependencies": [],
              "scope": {"contrast_ids": ["B_vs_A"]}}]
    root = tmp_path / "collection"
    plan, receipts, transport = _fixture_collection(root, {"project_id": "study", "contrasts": [{"contrast_id": "B_vs_A"}]},
        tasks, {"motif": (list(source.iterdir()), {})})
    before = {path.name: file_hash(path) for path in source.iterdir()}
    result = build_presentation_assets(root, plan, receipts, transport)
    view = result["extra_views"][0]
    assert view["audit"]["proof_kind"] == "flow"
    assert view["audit"]["models"] == 3
    assert view["audit"]["new_scans"] == view["audit"]["new_statistical_tests"] == 0
    assert len(view["panels"]) == 8
    assert before == {path.name: file_hash(path) for path in source.iterdir()}


def test_real_motif_task_and_collection_preserve_the_sealed_family(tmp_path):
    from katana_next.full_fixture import create_full_fixture
    from katana_next.flow.runtime import execute_task, read_json
    from katana_next.flow.collect import collect_results, materialize_report_descriptor
    project = create_full_fixture(tmp_path / "source")
    native = read_json(project)
    cid = native["contrasts"][0]["contrast_id"]
    native["metadata"]["species"] = "Mus_musculus"
    native["assets"] = [asset for asset in native["assets"]
                        if asset["kind"] != "rmats_result" or
                        (asset["metadata"]["contrast_id"] == cid and asset["metadata"]["event_type"] == "SE")]
    write_json(project, native)
    fasta = project.parent / "inputs/fixture.genome.fa"
    sequence = fasta.read_text().splitlines()[1]
    length = len(sequence)
    rmats = next(asset for asset in native["assets"] if asset["kind"] == "rmats_result")
    chromosomes = sorted(pd.read_csv(rmats["path"], sep="\t")["chr"].unique())
    chunks, index, offset = [], [], 0
    for chromosome in chromosomes:
        header = ">" + chromosome + "\n"
        offset += len(header)
        index.append(f"{chromosome}\t{length}\t{offset}\t{length}\t{length + 1}\n")
        chunks.append(header + sequence + "\n")
        offset += length + 1
    fasta.write_text("".join(chunks), encoding="ascii", newline="\n")
    Path(str(fasta) + ".fai").write_text("".join(index))
    motif = tmp_path / "motif.tsv"
    motif.write_text("motif_id\tRBP\tspecies\tevidence\tsource\tmotif\n"
                     "M\tR\tMus_musculus\tuser_defined\tsynthetic test only\tAC\n")
    tasks = [{"id": "as", "workflow_id": "differential.as", "dependencies": [],
              "parameters": {"contrast_ids": [cid]}, "scope": {"contrast_ids": [cid]}},
             {"id": "motif", "workflow_id": "prioritize.motiftrace_consequences", "dependencies": ["as"],
              "parameters": {"contrast_ids": [cid], "custom_motif_path": str(motif), "cluster_permutations": 99},
              "scope": {"contrast_ids": [cid]},
              "selection": {"analysis_id": "differential.as", "name": "Differential AS events", "contrast_ids": [cid]}}]
    plan = seal_plan({"project_manifest": str(project), "project_manifest_sha256": file_hash(project),
                      "contrast_ids": [cid], "tasks": tasks})
    plan_path = write_json(tmp_path / "plan.json", plan)
    first = execute_task(plan_path, "as", tmp_path / "as")
    second = execute_task(plan_path, "motif", tmp_path / "motif", [first])
    task = read_json(second)
    marker = Path(next(record["path"] for record in task["outputs"] if Path(record["path"]).name == "COMPLETE.json"))
    frozen = read_json(marker)
    scientific_plan = read_json(marker.parent / "ANALYSIS_PLAN_BEFORE_SCORING.json")
    assert scientific_plan["provenance"]["statistical_family_split"] is False
    originals = {name: file_hash(marker.parent / name) for name in frozen["members"]}
    collection = collect_results(plan_path, [first, second], tmp_path / "collection")
    descriptor = materialize_report_descriptor(collection)
    motif_artifacts = [a for a in descriptor["artifacts"] if a["analysis_id"] == "prioritize.motiftrace_consequences"]
    assert motif_artifacts and all(a["contrast_ids"] == [cid] for a in motif_artifacts)
    from katana_next.reports.sources import load_source
    assert load_source(descriptor)["report_selection"]["contrast_ids"] == [cid]
    view = next(view for view in descriptor["extra_views"] if view.get("audit", {}).get("proof_kind") == "flow")
    assert view["audit"]["new_scans"] == view["audit"]["new_statistical_tests"] == 0
    assert originals == {name: file_hash(marker.parent / name) for name in frozen["members"]}
    copied_marker = next((tmp_path / "collection/products").rglob("COMPLETE.json"))
    assert copied_marker.read_bytes() == marker.read_bytes()
    assert all(file_hash(copied_marker.parent / name) == record["sha256"] for name, record in frozen["members"].items())
    from katana_next.reports import build_report
    report = build_report(descriptor, tmp_path / "rendered")
    assert (tmp_path / "rendered/REPORT.html").is_file()
    html = (tmp_path / "rendered/REPORT.html").read_text(encoding="utf-8")
    assert all(f'id="chapter-{chapter}"' in html for chapter in range(1, 8))
    assert report["contrasts"] == [cid]
