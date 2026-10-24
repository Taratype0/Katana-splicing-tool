"""Actual report archives accept immutable artifacts with pre-ZIP timestamps."""
import os
import zipfile

from katana_next.reports.engine import Builder


def test_real_report_zip_keeps_epoch_zero_artifact_bytes_and_source_mtime(tmp_path):
    # Portable, content-addressed collections may deliberately carry epoch 0.
    # The existing source and its copied report member must never be retimed.
    artifact = tmp_path / "existing-native-results.tsv"
    native_bytes = b"model_id\tP_value\tBonferroni_P\nACTG\t5e-310\t2.72e-306\n"
    artifact.write_bytes(native_bytes)
    os.utime(artifact, ns=(0, 0))
    before = artifact.read_bytes(), artifact.stat().st_mtime_ns
    assert before[1] == 0
    output = tmp_path / "report"
    output.mkdir()
    project = {"project_id": "epoch-report", "name": "Archive regression", "samples": [],
               "contrasts": [], "annotation": {}}
    builder = Builder({"project": project, "artifacts": [], "mode": "bulk"},
                      output, "stable_significant")
    record = builder.table({"path": str(artifact), "name": "Native result",
                            "analysis_id": "prioritize.mapp_native", "chapter": 7})
    copied = output / record["source"]
    assert copied.read_bytes() == native_bytes and copied.stat().st_mtime_ns == 0

    # Exercise the real seven-chapter renderer and its actual ZIP writer.
    builder.write()

    assert (artifact.read_bytes(), artifact.stat().st_mtime_ns) == before
    assert copied.read_bytes() == native_bytes and copied.stat().st_mtime_ns == 0
    with zipfile.ZipFile(output / "REPORT.zip") as archive:
        assert archive.testzip() is None
        assert archive.read(record["source"]) == native_bytes
        assert archive.getinfo(record["source"]).date_time == (1980, 1, 1, 0, 0, 0)
        assert archive.read("REPORT.html") == (output / "REPORT.html").read_bytes()
        assert archive.read("report-manifest.json") == (output / "report-manifest.json").read_bytes()
