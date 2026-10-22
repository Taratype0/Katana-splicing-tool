"""Descriptive network defaults stay visible without reclassifying old results."""
from __future__ import annotations

import hashlib
import json
import posixpath
from types import SimpleNamespace
import zipfile

from bs4 import BeautifulSoup
import numpy as np
import pandas as pd
import pytest

from katana_next.full_fixture import create_full_fixture
from katana_next.services.project_service import ProjectService
from katana_next.services.module_dashboard import module_dashboard
from katana_next.workflows.base import ArtifactWriter
from katana_next.workflows.registry import WorkflowRegistry
from katana_next.workflows.reporting import ReproducibilityReportWorkflow, _network_label_policy_html
from katana_next.workflows.narrative_report import _append_network_module_methods
from katana_next.science.network_label_report import create_network_label_report


POLICY = {"policy_id": "network-descriptive-labels-opt-in-v1", "apply_threshold_labels": False, "hard_filter": False}


def test_report_html_pdf_and_portable_provenance_keep_label_off_distinct_from_unknown(tmp_path, monkeypatch):
    path = create_full_fixture(tmp_path / "project")
    context = ProjectService(tmp_path / "app").open_project(path)
    registry = WorkflowRegistry()
    writer = ArtifactWriter(context, "programs.cosplice", "recorded-network-labels-off", {}, registry.get("programs.cosplice").version,
                            common_metadata={"network_label_policy": dict(POLICY)})
    source = writer.table("Saved continuous metrics", "saved_metrics.parquet", pd.DataFrame([{"event": "e1", "kME": 0.9}]), entity_type="network")
    writer.finish()
    before = source.path.read_bytes()
    pdf_paragraphs = []
    from reportlab.platypus import SimpleDocTemplate
    build = SimpleDocTemplate.build

    def capture(document, items, *args, **kwargs):
        pdf_paragraphs.extend(item.getPlainText() for item in items if hasattr(item, "getPlainText"))
        return build(document, items, *args, **kwargs)

    monkeypatch.setattr(SimpleDocTemplate, "build", capture)
    ReproducibilityReportWorkflow().run(context, SimpleNamespace(report=lambda *_: None), "labels-report-new-run", {"latest_only": False, "include_result_files": True})
    outputs = context.catalog.list_artifacts(context.manifest.project_id, analysis_id="reports.bundle")
    html = next(a.path for a in outputs if a.path.suffix == ".html").read_text(encoding="utf-8")
    panel = BeautifulSoup(html, "html.parser").select_one(".network-label-policy")
    row = [cell.get_text() for cell in panel.select("tbody tr td")]
    assert row == [source.run_fingerprint, POLICY["policy_id"], "off", "no", "not recorded"]
    assert "without independent calibration" in panel.get_text()
    assert "not a calibrated probability" in panel.get_text()
    assert "FDR <=0.05" in panel.get_text()
    assert "threshold labels: off; hard filtering: no" in "\n".join(pdf_paragraphs)
    package_path = next(a.path for a in outputs if a.path.suffix == ".zip")
    with zipfile.ZipFile(package_path) as package:
        assert package.testzip() is None
        provenance = json.loads(package.read("provenance/provenance.json"))
        assert provenance["artifacts"][0]["metadata"]["network_label_policy"] == POLICY
    legacy = source.model_copy(update={"metadata": {}, "parameters": {"apply_threshold_labels": False}})
    historical = BeautifulSoup(_network_label_policy_html("programs.cosplice", [legacy]), "html.parser")
    assert [cell.get_text() for cell in historical.select("tbody tr td")][1:] == ["not recorded"] * 4
    assert source.path.read_bytes() == before


def test_one_click_appendix_keeps_original_structure_and_old_unknown_label_status(tmp_path):
    source = tmp_path / "approved-original.html"
    source.write_text("<html><body><main><section id='chapter7'><h2>Network results</h2><p>Original result.</p></section></main></body></html>", encoding="utf-8")
    original = source.read_bytes()
    output = tmp_path / "new-copy.html"
    output.write_bytes(original)
    record = {"status": "recorded", "network_module_parameters": {"min_module_size": 37, "merge_cut_height": 0.12},
              "network_label_policy": dict(POLICY), "network_run_fingerprint": "exact-linked-run"}
    appendix = _append_network_module_methods(output, record)
    soup = BeautifulSoup(output.read_text(encoding="utf-8"), "html.parser")
    assert soup.find(id="chapter7").get_text() == "Network resultsOriginal result."
    assert "Threshold labels: off. Hard filtering: no." in soup.select_one(".network-label-policy").get_text()
    assert json.loads(appendix.read_text())["network_label_policy"] == POLICY
    assert source.read_bytes() == original
    unknown = tmp_path / "historical-copy.html"
    unknown.write_bytes(original)
    _append_network_module_methods(unknown, {"status": "not_recorded", "reason": "No original policy record."})
    old_policy = BeautifulSoup(unknown.read_text(encoding="utf-8"), "html.parser").select_one(".network-label-policy").get_text()
    assert "Threshold labels: not recorded. Hard filtering: not recorded." in old_policy


def test_dashboard_keeps_every_metric_row_and_na_hub_labels_when_optional_labels_are_off(tmp_path):
    path = create_full_fixture(tmp_path / "project")
    context = ProjectService(tmp_path / "app").open_project(path)
    input_hashes = {p: hashlib.sha256(p.read_bytes()).hexdigest() for p in context.manifest.input_root.rglob("*") if p.is_file()}
    registry = WorkflowRegistry()
    events = [f"e{i:02}" for i in range(25)]  # More than any old implicit top-ten view.
    samples = [s.sample_id for s in context.manifest.samples]
    discovery = ArtifactWriter(context, "programs.cosplice", "discovery-labels-off", {}, registry.get("programs.cosplice").version,
                               common_metadata={"network_label_policy": dict(POLICY)})
    assignments = discovery.table("Fixed modules", "consensus_modules__assignments.parquet", pd.DataFrame({"event_key": events, "consensus_module": "M1"}), entity_type="module")
    rng = np.random.default_rng(21)
    discovery.table("Measured PSI", "matrix__psi.parquet", pd.DataFrame(rng.uniform(.1, .9, (25, len(samples))), columns=samples).assign(event_key=events), entity_type="matrix")
    discovery.finish()
    metrics = pd.DataFrame({"event_key": events, "gene": events, "module": "M1", "condition": "all",
                            "kWithin": np.arange(25, 0, -1, dtype=float), "kWithin_no_same_gene": 1.0,
                            "kME": .9, "rank_all_edges": np.arange(1, 26), "rank_no_same_gene": np.arange(1, 26),
                            "supported_degree": 24, "cross_gene_supported_degree": 24,
                            "hub_candidate": pd.array([pd.NA] * 25, dtype="boolean"),
                            "label_status": "not_applied", "estimation_status": "estimable"})
    hub = ArtifactWriter(context, "programs.hub_events", "hub-metrics-labels-off", {}, registry.get("programs.hub_events").version,
                         common_metadata={"source_discovery_run": assignments.run_fingerprint, "network_label_policy": dict(POLICY)})
    hub.table("Complete hub metrics", "hub_events__metrics.parquet", metrics, entity_type="hub")
    hub.finish()
    result = module_dashboard(context, assignments.artifact_id, "M1")
    assert result["hub_status"] == "labels_not_applied"
    assert result["hub_label_status"] == "not_applied"
    assert result["hub_events"] == []
    assert len(result["hub_comparison"]) == 25
    assert [row["rank_all_edges"] for row in result["hub_comparison"]] == list(range(1, 26))
    assert all(row["hub_candidate"] is None for row in result["hub_comparison"])
    assert result["hub_label_policy"] == POLICY
    assert input_hashes == {p: hashlib.sha256(p.read_bytes()).hexdigest() for p in input_hashes}


def test_standalone_report_companion_uses_current_tables_and_never_relabels_frozen_html(tmp_path):
    raw = tmp_path / "frozen-report.html"
    raw.write_text("<html><p>OLD CORE LABEL ONLY</p></html>", encoding="utf-8")
    current = tmp_path / "current-members.tsv"
    current.write_text("event\tkME\tmembership_role\ne1\t0.95\tnot_applied\ne2\t0.5\tnot_applied\n", encoding="utf-8")
    original = raw.read_bytes()
    out = tmp_path / "new-view" / "REPORT.html"
    assert create_network_label_report(raw, out, policy=POLICY, table_views={"members": current}) == out
    assert create_network_label_report(raw, out, policy=POLICY, table_views={"members": current}) == out
    report = BeautifulSoup(out.read_text(encoding="utf-8"), "html.parser")
    assert "Threshold labels: off" in report.get_text()
    assert "historical/internal reference" in report.get_text()
    assert "OLD CORE LABEL ONLY" not in report.get_text()
    assert (out.parent / "tables/members.tsv").read_bytes() == current.read_bytes()
    assert raw.read_bytes() == original
    provenance = json.loads((out.parent / "NETWORK_LABEL_REPORT_PROVENANCE.json").read_text())
    assert provenance["legacy_report_reference"]["sha256"] == hashlib.sha256(original).hexdigest()
    assert provenance["network_label_policy"] == POLICY
    before = out.read_bytes()
    current.write_text("changed-current-table\n", encoding="utf-8")
    with pytest.raises(FileExistsError, match="Refusing to replace"):
        create_network_label_report(raw, out, policy=POLICY, table_views={"members": current})
    assert out.read_bytes() == before
    with pytest.raises(ValueError, match="separate"):
        create_network_label_report(raw, raw, policy=POLICY, table_views={})


def test_portable_package_preserves_companion_report_attachment_links(tmp_path):
    manifest = create_full_fixture(tmp_path / "project")
    context = ProjectService(tmp_path / "app").open_project(manifest)
    writer = ArtifactWriter(context, "programs.cosplice", "current-companion-report", {}, WorkflowRegistry().get("programs.cosplice").version)
    raw = tmp_path / "raw-frozen-report.html"
    raw.write_text("<html>Historical core labels</html>", encoding="utf-8")
    view = tmp_path / "current-view.tsv"
    view.write_text("event\tlabel_status\ne1\tnot_applied\n", encoding="utf-8")
    companion = create_network_label_report(raw, writer.directory / "report__report.html", policy=POLICY, table_views={"members": view})
    writer.existing_file("Current label-view report", companion, entity_type="report", metadata={
        "preserve_relative_report_links": True,
        "network_label_policy": dict(POLICY),
        "exports": {"current_members": str(companion.parent / "tables/members.tsv"),
                    "label_provenance": str(companion.parent / "NETWORK_LABEL_REPORT_PROVENANCE.json")},
    })
    writer.finish()
    ReproducibilityReportWorkflow().run(context, SimpleNamespace(report=lambda *_: None), "portable-companion-bundle", {"latest_only": True, "include_result_files": True})
    output = next(a.path for a in context.catalog.list_artifacts(context.manifest.project_id, analysis_id="reports.bundle") if a.path.suffix == ".zip")
    with zipfile.ZipFile(output) as package:
        assert package.testzip() is None
        html_name = next(name for name in package.namelist() if name.endswith("/report__report.html"))
        report = BeautifulSoup(package.read(html_name), "html.parser")
        for link in report.find_all('a', href=True):
            attached = posixpath.normpath(posixpath.join(posixpath.dirname(html_name), link['href']))
            assert attached in package.namelist(), attached
        assert package.read(posixpath.join(posixpath.dirname(html_name), "tables/members.tsv")) == view.read_bytes()


def test_new_delivery_report_keeps_old_network_content_explicitly_scoped_as_reference(tmp_path):
    path = tmp_path / "generated-report.html"
    path.write_text("<html><body><main><section id='s7'><h2>7. Katana-Splice-Network</h2><p>Old core labels</p></section></main></body></html>", encoding="utf-8")
    _append_network_module_methods(path, {"status": "recorded", "network_module_parameters": {"min_module_size": 20, "merge_cut_height": .25}, "network_label_policy": POLICY}, mark_delivery_references=True)
    report = BeautifulSoup(path.read_text(encoding="utf-8"), "html.parser")
    section = report.find(id="s7")
    assert "Old core labels" in section.get_text()
    assert "historical/internal references" in section.select_one('.network-delivery-reference').get_text()
    assert "must not be relabelled" in section.get_text()
