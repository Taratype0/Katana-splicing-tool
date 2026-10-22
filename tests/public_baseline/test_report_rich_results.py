"""Completed-result, scope and portable-copy contracts for rich report assets."""
from hashlib import sha256
import json
from pathlib import Path

import pytest

from katana_next.reports.rich_results import collect_rich_results, copy_rich_asset


def artifact(path, *, workflow="programs.cosplice", name=None, **changes):
    row = {"artifact_id": "artifact-1", "analysis_id": workflow,
           "path": str(path), "name": name or path.stem, "status": "complete",
           "run_fingerprint": "run-A", "contrast_ids": ["B_vs_A"],
           "parameters": {}, "metadata": {"workflow_version": "1.0.0"}}
    row.update(changes)
    return row


def available(rows):
    return [row for row in rows if row["status"] == "available"]


def test_missing_result_is_not_zero_or_borrowed(tmp_path):
    rows = collect_rich_results([], project_root=tmp_path)
    assert {row["chapter"] for row in rows} == {6, 7}
    assert all(row["status"] == "not_computed" and not row["path"] for row in rows)
    assert all(row["metadata"]["new_computation"] is False for row in rows)


def test_versioned_workflow_classification_not_filename_guess(tmp_path):
    image = tmp_path / "motiftrace_network.png"
    image.write_bytes(b"existing science figure")
    assert not available(collect_rich_results([artifact(image, workflow="unknown.research")]))
    assert not available(collect_rich_results([artifact(image, status="failed")]))
    assert not available(collect_rich_results([artifact(image, status="running")]))
    assert not available(collect_rich_results([artifact(image, status="cancelled")]))


def test_all_figures_tables_and_models_retained_with_scope(tmp_path):
    paths = [tmp_path / f"module_{index}.png" for index in range(15)]
    paths += [tmp_path / "all_members.tsv", tmp_path / "all_hub_metrics.tsv"]
    for path in paths:
        path.write_text("event\tmetric\nE1\t0.1\n", encoding="utf-8")
    artifacts = [artifact(path, artifact_id=f"id-{i}") for i, path in enumerate(paths)]
    rows = available(collect_rich_results(artifacts))
    assert len(rows) == 17
    assert all(row["chapter"] == 6 and row["contrast_id"] == "B_vs_A" for row in rows)
    assert all(row["metadata"]["all_rows_retained"] for row in rows)
    assert all(row["metadata"]["das_scope"]["mode"] == "measurement_population" for row in rows)
    assert all(row["metadata"]["new_computation"] is False for row in rows)


def test_same_file_different_contrast_and_run_preserved(tmp_path):
    table = tmp_path / "enrichment.tsv"
    table.write_text("RBP\tFDR\nX\t0.0314159265\n")
    first = artifact(table, workflow="prioritize.event_motifs")
    second = artifact(table, workflow="prioritize.event_motifs", contrast_ids=["C_vs_A"], run_fingerprint="run-B")
    rows = available(collect_rich_results([first, first, second]))
    assert len(rows) == 2
    assert {row["contrast_id"] for row in rows} == {"B_vs_A", "C_vs_A"}
    assert all(row["metadata"]["das_scope"]["mode"] == "rs_required" for row in rows)
    # Workflow requirements do not fabricate evidence that an old run used it.
    assert all(row["metadata"]["das_scope"]["recorded_event_set"] is None for row in rows)
    assert all(row["metadata"]["das_scope"]["scope_recorded"] is False for row in rows)


def test_fixed_import_scope_not_relabelled_as_das(tmp_path):
    table = tmp_path / "model_site_statistics.tsv"
    table.write_text("model\tq\nPWM1\tNA\n")
    row = available(collect_rich_results([artifact(table, workflow="prioritize.motiftrace_usage")]))[0]
    assert row["metadata"]["das_scope"]["mode"] == "fixed_completed_result"
    assert row["metadata"]["statistical_values_preserved"]


def test_rbp_identity_does_not_claim_motiftrace_or_clip(tmp_path):
    table = tmp_path / "rbp_members.tsv"
    table.write_text("gene_id\tRBP\nG1\tX\n")
    rows = collect_rich_results([artifact(table, workflow="prioritize.rbp_support")])
    assert len(available(rows)) == 1
    missing = {row["metadata"].get("component") for row in rows if row["status"] == "not_computed"}
    assert {"motif_enrichment", "motiftrace_sequence", "clip", "mapp"} <= missing
    assert "rbp_identity" not in missing


def test_only_completed_manifest_imports_and_all_exports(tmp_path):
    run = tmp_path / "run"
    run.mkdir()
    plots = run / "plots"
    plots.mkdir()
    table = run / "assignments.tsv"
    table.write_text("event\tmodule\nE1\tM1\n")
    image = plots / "distribution.png"
    image.write_bytes(b"plot")
    spec = run / "distribution.vl.json"
    spec.write_text('{"data":{"values":[]}}')
    marker = run / "parameters.json"
    payload = {"workflow_id": "programs.cosplice", "workflow_version": "1.0",
               "completed": True, "run_fingerprint": "from-manifest", "parameters": {},
               "generated_artifacts": [{"artifact_id": "assignment", "name": "All assignments",
                                        "path": str(table)}, {"artifact_id": "spec", "path": str(spec)}],
               "generated_plot_files": [str(image)], "generated_table_files": []}
    marker.write_text(json.dumps(payload))
    rows = available(collect_rich_results([], search_roots=[tmp_path]))
    assert {row["path"] for row in rows} == {str(table.resolve()), str(image.resolve()), str(spec.resolve())}
    assert len(rows) == 3
    assert next(row for row in rows if row["kind"] == "image")["title"] == "distribution"
    payload["completed"] = False
    marker.write_text(json.dumps(payload))
    assert not available(collect_rich_results([], search_roots=[tmp_path]))


def test_browser_dependencies_copied_relatively_originals_unchanged(tmp_path):
    run = tmp_path / "source"
    run.mkdir()
    (run / "plots").mkdir()
    (run / "assets").mkdir()
    image = run / "plots" / "all_modules.svg"
    image.write_text('<svg xmlns="http://www.w3.org/2000/svg"/>')
    style = run / "assets" / "style.css"
    style.write_text("body{color:blue}")
    table = run / "full_members.tsv"
    table.write_text("event\tmodule\nE1\tM1\nE2\tM2\n")
    browser = run / "BROWSER.html"
    browser.write_text('<link href="assets/style.css" rel="stylesheet"><img src="plots/all_modules.svg">'
                       '<a href="full_members.tsv">All members</a>')
    paths = [browser, image, style, table]
    before = {path: (path.stat().st_mtime_ns, sha256(path.read_bytes()).hexdigest()) for path in paths}
    row = next(row for row in available(collect_rich_results([artifact(browser)])) if row["kind"] == "html")
    assert len(row["metadata"]["local_dependencies"]) == 3
    assert not row["metadata"]["portability_warnings"]
    receipt = copy_rich_asset(row, tmp_path / "portable")
    assert receipt["originals_preserved"]
    assert Path(receipt["path"]).read_bytes() == browser.read_bytes()
    assert (tmp_path / "portable" / "plots" / "all_modules.svg").read_bytes() == image.read_bytes()
    assert {path: (path.stat().st_mtime_ns, sha256(path.read_bytes()).hexdigest()) for path in paths} == before


def test_app_api_links_rewrite_only_destination_and_keep_native_fdr(tmp_path):
    source = tmp_path / "source"
    source.mkdir()
    table = tmp_path / "original_model_statistics.tsv"
    table.write_text("model\tFDR\nM1\t0.0000123456789012345\n")
    browser = source / "motifs.html"
    url = "/api/projects/p1/artifacts/table-A/file?download=true"
    browser.write_text(f'<a href="{url}">All model results</a>')
    old_html, old_tsv = browser.read_bytes(), table.read_bytes()
    rows = collect_rich_results([artifact(browser, workflow="prioritize.rbp_sequence_results"),
                                 artifact(table, artifact_id="table-A", workflow="prioritize.rbp_sequence_results")])
    row = next(row for row in available(rows) if row["kind"] == "html")
    assert not row["metadata"]["portability_warnings"]
    receipt = copy_rich_asset(row, tmp_path / "portable")
    exported = Path(receipt["path"]).read_text()
    assert url not in exported
    assert '_related/table-A/original_model_statistics.tsv' in exported
    assert (tmp_path / "portable" / "_related" / "table-A" / table.name).read_bytes() == old_tsv
    assert browser.read_bytes() == old_html and table.read_bytes() == old_tsv
    assert receipt["link_rewrites"]


def test_unresolved_and_outside_browser_links_are_not_claimed_portable(tmp_path):
    source = tmp_path / "source"
    source.mkdir()
    secret = tmp_path / "outside.tsv"
    secret.write_text("outside")
    browser = source / "motifs.html"
    browser.write_text('<a href="../outside.tsv">outside</a><a href="missing.tsv">missing</a>'
                       '<a href="/api/projects/p1/artifacts/missing/file">app</a>')
    row = available(collect_rich_results([artifact(browser)]))[0]
    assert len(row["metadata"]["portability_warnings"]) == 3
    assert not row["metadata"]["local_dependencies"]


def test_copy_rejects_source_overlap_or_overwrite(tmp_path):
    source = tmp_path / "source"
    source.mkdir()
    table = source / "full.tsv"
    table.write_text("event\tPSI\nE1\t0.4\n")
    row = available(collect_rich_results([artifact(table)]))[0]
    with pytest.raises(ValueError, match="independent"):
        copy_rich_asset(row, source / "do-not-write")
    destination = tmp_path / "existing"
    destination.mkdir()
    with pytest.raises(ValueError, match="already exists"):
        copy_rich_asset(row, destination)
    assert not (source / "do-not-write").exists()


def test_plot_alternatives_retained_and_vega_file_data_copied(tmp_path):
    image = tmp_path / "distribution.png"
    image.write_bytes(b"existing PNG")
    svg = tmp_path / "distribution.svg"
    svg.write_text("<svg/>")
    data = tmp_path / "distribution.tsv"
    data.write_text("module\tn\nM1\t21\n")
    spec = tmp_path / "distribution.vl.json"
    spec.write_text(json.dumps({"data": {"url": data.name}}))
    row = artifact(spec, metadata={"exports": {"png": str(image), "svg": str(svg)}})
    rows = available(collect_rich_results([row]))
    assert len(rows) == 3
    assert next(r for r in rows if r["path"] == str(svg))["metadata"]["display_in_main"] is False
    assert next(r for r in rows if r["path"] == str(image))["metadata"]["display_in_main"] is True
    original = next(r for r in rows if r["kind"] == "vega")
    assert original["metadata"]["local_dependencies"][0]["path"] == str(data)


def test_app_link_nested_browser_dependencies_keep_relative_tree(tmp_path):
    root = tmp_path / "first"
    root.mkdir()
    second = tmp_path / "second"
    second.mkdir()
    (second / "assets").mkdir()
    style = second / "assets" / "browser.css"
    style.write_text("body{color:blue}")
    child = second / "child.html"
    child.write_text('<link href="assets/browser.css" rel="stylesheet">')
    browser = root / "main.html"
    browser.write_text('<a href="/api/projects/p1/artifacts/child/file">Child browser</a>')
    rows = collect_rich_results([artifact(browser), artifact(child, artifact_id="child")])
    row = next(row for row in available(rows) if row["path"] == str(browser))
    assert not row["metadata"]["portability_warnings"]
    receipt = copy_rich_asset(row, tmp_path / "portable")
    assert (tmp_path / "portable" / "_related" / "child" / "assets" / "browser.css").is_file()
    assert len(receipt["files"]) == 3


def test_cross_project_inputs_rejected_even_with_same_biology_label(tmp_path):
    first = tmp_path / "first.tsv"
    first.write_text("event\nE1\n")
    second = tmp_path / "second.tsv"
    second.write_text("event\nE2\n")
    with pytest.raises(ValueError, match="cohorts cannot be pooled"):
        collect_rich_results([artifact(first, project_id="study1"), artifact(second, project_id="study2")])


def test_registered_file_cannot_override_unfinished_producer_marker(tmp_path):
    table = tmp_path / "members.tsv"
    table.write_text("event\tmodule\nE1\tM1\n")
    marker = tmp_path / "parameters.json"
    value = {"workflow_id": "programs.cosplice", "completed": False,
             "generated_artifacts": [{"artifact_id": "artifact-1", "path": str(table)}]}
    marker.write_text(json.dumps(value))
    assert not available(collect_rich_results([artifact(table)]))
    value["completed"] = True
    marker.write_text(json.dumps(value))
    row = available(collect_rich_results([artifact(table)]))[0]
    assert row["metadata"]["completion_evidence"] == "completed_workflow_manifest"
    strict = artifact(table, metadata={"artifact_contract_version": "3.0-correctness"})
    value["generated_artifacts"] = []
    marker.write_text(json.dumps(value))
    assert not available(collect_rich_results([strict]))


def test_sequence_substep_does_not_claim_completed_motif_inference(tmp_path):
    image = tmp_path / "actual_sequence.png"
    image.write_bytes(b"existing sequence figure")
    rows = collect_rich_results([artifact(image, workflow="prioritize.motiftrace_consequences",
                                         metadata={"completion_scope": "sequence_consequences_only"})])
    assert len(available(rows)) == 1
    real = available(rows)[0]
    assert "inference_not_completed" in real["metadata"]["partial_status"]
    assert "not established" in real["caption"]
    missing = {row["metadata"].get("component") for row in rows if row["status"] == "not_computed"}
    assert "motif_enrichment" in missing
    assert "motiftrace_sequence" not in missing


def test_network_baseline_does_not_claim_bootstrap_confidence(tmp_path):
    image = tmp_path / "observed_modules.png"
    image.write_bytes(b"completed network baseline")
    base = artifact(image, metadata={"completion_scope": "observed_network_baseline_only"})
    rows = collect_rich_results([base])
    assert any(row["metadata"].get("component") == "network_bootstrap" and row["status"] == "not_computed" for row in rows)
    assert "own completed phase" in available(rows)[0]["caption"]
    boot = tmp_path / "bootstrap.tsv"
    boot.write_text("module\tconfidence\nM1\t0.82\n")
    rows = collect_rich_results([base, artifact(boot, artifact_id="bootstrap", metadata={"completion_scope": "bootstrap_confidence"})])
    assert not any(row["metadata"].get("component") == "network_bootstrap" for row in rows)


def test_explicit_independent_portable_root_retains_parent_table_links(tmp_path):
    root = tmp_path / "network-view"
    (root / "browser").mkdir(parents=True)
    (root / "tables").mkdir()
    table = root / "tables" / "all_members.tsv"
    table.write_text("event\tmodule\ne1\tM1\n")
    browser = root / "browser" / "ALL_MODULES.html"
    browser.write_text('<a href="../tables/all_members.tsv">Every member</a>')
    rows = collect_rich_results([artifact(browser, metadata={"portable_root": str(root), "display_in_main": False})])
    row = available(rows)[0]
    assert row["metadata"]["display_in_main"] is False
    assert not row["metadata"]["portability_warnings"]
    receipt = copy_rich_asset(row, tmp_path / "portable")
    assert (tmp_path / "portable" / "tables" / "all_members.tsv").read_bytes() == table.read_bytes()
    assert Path(receipt["path"]).name == "ALL_MODULES.html"


def test_declared_portable_root_cannot_relocate_unrelated_asset(tmp_path):
    source = tmp_path / "members.tsv"
    source.write_text("event\ne1\n")
    unrelated = tmp_path / "unrelated"
    unrelated.mkdir()
    with pytest.raises(ValueError, match="must contain"):
        collect_rich_results([artifact(source, metadata={"portable_root": str(unrelated)})])
