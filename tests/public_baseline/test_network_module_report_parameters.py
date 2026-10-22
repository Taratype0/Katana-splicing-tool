"""The report must show executed network controls, never infer old defaults."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from types import SimpleNamespace
import zipfile

from bs4 import BeautifulSoup
import pandas as pd

from katana_next.domain.models import AnalysisArtifact, ArtifactKind
from katana_next.full_fixture import create_full_fixture
from katana_next.services.project_service import ProjectService
from katana_next.workflows.base import ArtifactWriter
from katana_next.workflows.programs import CoSpliceWorkflow
from katana_next.workflows.reporting import (
    ReproducibilityReportWorkflow,
    _network_module_parameter_records,
)


def _hashes(root: Path) -> dict[str, str]:
    return {
        str(path.relative_to(root)): hashlib.sha256(path.read_bytes()).hexdigest()
        for path in root.rglob("*") if path.is_file()
    }


def _record_network(context, run: str, minimum: int, merge: float):
    # Requested defaults intentionally differ: the report must use saved actuals.
    writer = ArtifactWriter(
        context, "programs.cosplice", run,
        {"min_module_size": 20, "merge_cut_height": 0.25},
        CoSpliceWorkflow.version,
        common_metadata={
            "question": "Which toy events form modules?",
            "network_module_parameters": {
                "min_module_size": minimum, "merge_cut_height": merge,
            },
        },
    )
    artifact = writer.table(
        "Saved toy network module summary", "module_summary.parquet",
        pd.DataFrame([{"row_id": "M1", "module": "M1", "events": 45}]),
        entity_type="network_module",
    )
    writer.finish()
    return artifact


def test_network_report_exports_actual_per_run_parameters_without_reanalysis(tmp_path, monkeypatch):
    manifest_path = create_full_fixture(tmp_path / "project")
    context = ProjectService(application_root=tmp_path / "state").open_project(manifest_path)
    source_hashes = _hashes(context.manifest.input_root)
    manifest_hash = hashlib.sha256(manifest_path.read_bytes()).hexdigest()
    first = _record_network(context, "run37-actual-network", 37, 0.12)
    second = _record_network(context, "run41-actual-network", 41, 0.17)
    old_path = context.manifest.output_root / "saved_historical_network.tsv"
    old_path.write_text("module\tevents\nold-module\t8\n", encoding="utf-8")
    historical = AnalysisArtifact(
        project_id=context.manifest.project_id, analysis_id="programs.cosplice",
        name="Historical network output", kind=ArtifactKind.TABLE,
        entity_type="network_module", fingerprint="historical-network-output",
        run_fingerprint="historical-network-run", path=old_path, row_count=1,
        parameters={"min_module_size": 999, "merge_cut_height": 0.99},
        metadata={"project_snapshot_fingerprint": context.snapshot.manifest_fingerprint},
    )
    context.catalog.register_artifact(historical)
    saved_outputs = {
        path: hashlib.sha256(path.read_bytes()).hexdigest()
        for path in [first.path, second.path, old_path]
    }

    # Capture semantic PDF paragraphs while still building the actual PDF file.
    from reportlab.platypus import SimpleDocTemplate
    original_build = SimpleDocTemplate.build
    pdf_paragraphs = []

    def build_and_capture(document, flowables, *args, **kwargs):
        pdf_paragraphs.extend(
            item.getPlainText() for item in flowables if hasattr(item, "getPlainText")
        )
        return original_build(document, flowables, *args, **kwargs)

    monkeypatch.setattr(SimpleDocTemplate, "build", build_and_capture)
    # Reporting is cache-only; attempting any scientific stage is an error.
    from katana_cosplice.pipeline import CoSplicePipeline
    monkeypatch.setattr(
        CoSplicePipeline, "run_stage",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(AssertionError("Unexpected reanalysis")),
    )
    workflow = ReproducibilityReportWorkflow()
    workflow.run(
        context, SimpleNamespace(report=lambda *_args: None), "report-actual-module-parameters",
        {"latest_only": False, "include_result_files": True},
    )
    artifacts = context.catalog.list_artifacts(context.manifest.project_id, analysis_id="reports.bundle")
    html_path = next(item.path for item in artifacts if item.path.suffix == ".html")
    html = BeautifulSoup(html_path.read_text(encoding="utf-8"), "html.parser")
    panel = html.select_one(".network-module-parameters")
    assert panel is not None
    rows = {
        cells[0].get_text(): [cell.get_text() for cell in cells[1:]]
        for row in panel.select("tbody tr") if (cells := row.find_all("td"))
    }
    assert rows == {
        "run37-actual-network": ["37", "0.12"],
        "run41-actual-network": ["41", "0.17"],
        "historical-network-run": ["not recorded", "not recorded"],
    }
    assert "999" not in panel.get_text()
    assert "0.25" not in panel.get_text()
    pdf_text = "\n".join(pdf_paragraphs)
    assert "Source run: run37-actual-network" in pdf_text
    assert "Minimum module size (events): 37" in pdf_text
    assert "Module merge cut height: 0.12" in pdf_text
    assert "Source run: run41-actual-network" in pdf_text
    assert "Minimum module size (events): 41" in pdf_text
    assert "Module merge cut height: 0.17" in pdf_text
    assert "Minimum module size (events): not recorded" in pdf_text
    pdf_path = next(item.path for item in artifacts if item.path.suffix == ".pdf")
    assert pdf_path.read_bytes().startswith(b"%PDF-")

    package_path = next(item.path for item in artifacts if item.path.suffix == ".zip")
    with zipfile.ZipFile(package_path) as package:
        assert package.testzip() is None
        provenance = json.loads(package.read("provenance/provenance.json"))
        by_run = {item["run_fingerprint"]: item for item in provenance["artifacts"]}
        assert by_run[first.run_fingerprint]["metadata"]["network_module_parameters"] == {
            "min_module_size": 37, "merge_cut_height": 0.12,
        }
        assert "network_module_parameters" not in by_run[historical.run_fingerprint]["metadata"]
        assert package.read("report/KATANA_EVIDENCE_REPORT.pdf") == pdf_path.read_bytes()
        assert package.read("report/KATANA_EVIDENCE_REPORT.html") == html_path.read_bytes()
        exported = [name for name in package.namelist() if name.startswith("artifacts/programs.cosplice/")]
        assert any(name.endswith("module_summary.parquet") for name in exported)
        assert any(name.endswith("saved_historical_network.tsv") for name in exported)
    assert source_hashes == _hashes(context.manifest.input_root)
    assert manifest_hash == hashlib.sha256(manifest_path.read_bytes()).hexdigest()
    assert saved_outputs == {path: hashlib.sha256(path.read_bytes()).hexdigest() for path in saved_outputs}
    assert "network_module_parameters" not in historical.metadata


def test_module_parameter_rows_preserve_actual_conflicts_and_skip_unrelated_analysis(tmp_path):
    manifest_path = create_full_fixture(tmp_path / "project")
    context = ProjectService(application_root=tmp_path / "state").open_project(manifest_path)
    artifact = _record_network(context, "same-run", 37, 0.12)
    conflicting = artifact.model_copy(update={
        "metadata": {"network_module_parameters": {"min_module_size": 41, "merge_cut_height": 0.17}},
    })
    rows = _network_module_parameter_records("programs.cosplice", [artifact, artifact, conflicting])
    assert len(rows) == 2  # Do not hide contradictory records from one source run.
    assert {(row["min_module_size"], row["merge_cut_height"]) for row in rows} == {
        ("37", "0.12"), ("41", "0.17"),
    }
    unrelated = artifact.model_copy(update={"analysis_id": "differential.deg", "metadata": {}})
    assert _network_module_parameter_records("differential.deg", [unrelated]) == []
