"""Portable target CLI and browser resources work outside the source checkout."""
from importlib import metadata
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import sysconfig
import zipfile

import pandas as pd


def _cli(root, *args):
    import katana_next.targets
    installed = os.environ.get("KATANA_WHEEL_ONLY") == "1"
    if installed:
        assert Path(katana_next.targets.__file__).resolve().is_relative_to(Path(sysconfig.get_path("purelib")))
    env = {**os.environ, "PYTHONDONTWRITEBYTECODE": "1", "MPLBACKEND": "Agg"}
    if not installed:
        env["PYTHONPATH"] = str(Path(katana_next.targets.__file__).resolve().parents[2])
    process = subprocess.run([sys.executable, *(["-I"] if installed else []), "-m", "katana_next.targets", *map(str, args)],
                             cwd=root, env=env, capture_output=True, text=True, timeout=120)
    assert process.returncode == 0, process.stderr
    return json.loads(process.stdout)


def test_installed_target_cli_relocates_without_the_original_project(tmp_path):
    from katana_next.full_fixture import create_full_fixture
    from katana_next.services.project_service import ProjectService
    from katana_next.workflows.base import ArtifactWriter
    from katana_next.science.bulk import READ_SUPPORT_POLICY_VERSION
    source = tmp_path / "source-machine"
    manifest = create_full_fixture(source)
    context = ProjectService().open_project(manifest)
    writer = ArtifactWriter(context, "differential.as", "portable-acceptance-source",
                            {"count_mode": "JC", "fdr": 0.05}, "3.0.0",
                            common_metadata={"read_support_policy": READ_SUPPORT_POLICY_VERSION})
    rows = pd.DataFrame([
        dict(row_id="selected", event_id="SE:event-1", event_type="SE", gene_id="G1", gene_symbol="Example",
             comparison_id="perturbation-vs-baseline", chr="chr1", strand="+", exonStart_0base=200, exonEnd=250,
             upstreamES=100, upstreamEE=150, downstreamES=300, downstreamEE=350, fdr=0.0123456789012345, dpsi=0.2),
        dict(row_id="not-significant", event_id="SE:event-2", event_type="SE", gene_id="G1", gene_symbol="Example",
             comparison_id="perturbation-vs-baseline", chr="chr1", strand="+", exonStart_0base=400, exonEnd=450,
             upstreamES=300, upstreamEE=350, downstreamES=500, downstreamEE=550, fdr=0.9, dpsi=0.001),
    ])
    artifact = writer.table("Complete measured AS", "events.parquet", rows, entity_type="splice_event")
    writer.finish()
    selection = tmp_path / "selection.json"
    selection.write_text(json.dumps([{"source_artifact_id": artifact.artifact_id, "row_ids": ["selected"],
                                     "target_kind": "gene", "filters": {"event_set": "stable_significant"}}]), encoding="utf-8")
    exported = _cli(tmp_path, "export", "--project", manifest, "--selections", selection,
                    "--module", "gene_card", "--module", "event_card", "--outdir", tmp_path / "export")
    destination = tmp_path / "other-machine"
    destination.mkdir()
    relocated = destination / "targets.zip"
    shutil.copyfile(exported["archive_path"], relocated)
    # This is generated fixture data owned by this test, never research data.
    shutil.rmtree(source)
    info = _cli(destination, "inspect", "--request", relocated)
    assert info["verified"] and info["targets"][0]["gene_id"] == "G1"
    tasks = _cli(destination, "tasks", "--request", relocated)
    assert len(tasks["tasks"]) == 1
    result = _cli(destination, "run", "--request", relocated, "--outdir", destination / "reports", "--workers", "1")
    assert result["status"] in {"complete", "cache_hit"}
    with zipfile.ZipFile(result["archive"]) as archive:
        assert archive.testzip() is None and "INDEX.html" in archive.namelist()
    products = {p: p.stat().st_mtime_ns for p in (destination / "reports").rglob("*") if p.is_file()}
    repeated = _cli(destination, "run", "--request", relocated, "--outdir", destination / "reports", "--workers", "1")
    assert repeated["status"] in {"cache_hit", "complete"}
    assert all(p.stat().st_mtime_ns == stamp for p, stamp in products.items())
    verified = _cli(destination, "results", "--outdir", destination / "reports")
    assert verified.get("status") != "failed"
    # A downloaded ZIP intentionally has no local owner or original exported
    # request ZIP. Restore must verify it and create both before cache reuse.
    returned = tmp_path / "return-machine"
    returned.mkdir()
    returned_zip = returned / "completed-targets.zip"
    shutil.copy2(result["archive"], returned_zip)
    restored = _cli(returned, "restore", "--results", returned_zip,
                    "--outdir", returned / "restored", "--project-id", context.manifest.project_id)
    assert restored["restored"] and restored["status"] == "complete"
    assert restored["restoration"]["original_completion_history_preserved"]
    assert Path(restored["request"]).name == "TARGET_REQUEST.zip"
    planned = _cli(returned, "plan", "--request", restored["request"], "--outdir", returned / "restored")
    assert planned["ready"] and all(task["reusable"] for task in planned["tasks"])
    restored_files = {p: (p.read_bytes(), p.stat().st_mtime_ns)
                      for p in (returned / "restored").rglob("*") if p.is_file()}
    reused = _cli(returned, "run", "--request", restored["request"],
                  "--outdir", returned / "restored", "--workers", "1")
    assert reused["status"] == "cache_hit" and reused["cache"]["hits"] == 1
    assert all((p.read_bytes(), p.stat().st_mtime_ns) == proof for p, proof in restored_files.items())


def test_browser_and_target_entrypoint_are_in_the_wheel():
    entries = metadata.distribution("katana-splicing-tool").entry_points
    if os.environ.get("KATANA_WHEEL_ONLY") == "1":
        entry = next(e for e in entries if e.name == "katana-targets")
        assert entry.value == "katana_next.targets.cli:main"
        from katana_next.web_assets import frontend_directory
        front = frontend_directory()
        assert front is not None and front.is_relative_to(Path(sysconfig.get_path("purelib")))
        assert list((front / "assets").glob("*.js"))
        assert "<div id=\"root\">" in (front / "index.html").read_text("utf-8")
