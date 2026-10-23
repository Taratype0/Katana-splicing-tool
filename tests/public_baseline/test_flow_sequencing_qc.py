"""Existing native QC becomes verified portable chapter-one figures, no refit."""
import json
from pathlib import Path
import shutil
import zipfile

import pytest

from katana_next.flow.model import file_sha256, write_json
from katana_next.flow.runtime import canonical_hash
from katana_next.flow.sequencing import execute_stage, readiness
from katana_next.flow.qc_views import build_qc_views


FASTQC = """##FastQC\t0.12.1
>>Basic Statistics\tpass
#Measure\tValue
Filename\tS1_1.fastq.gz
Total Sequences\t1000
Sequence length\t150
%GC\t50
>>END_MODULE
>>Per base sequence quality\twarn
#Base\tMean\tMedian
1\t39.123456789\t40
10-14\t38.999999999\t40
>>END_MODULE
>>Sequence Length Distribution\tpass
#Length\tCount
150\t1000
>>END_MODULE
"""


def study(root):
    root.mkdir()
    assets = []
    for filename, phase in (("S1_1_fastqc.zip", "raw"), ("S1_1_val_1_fastqc.zip", "trimmed")):
        path = root / filename
        with zipfile.ZipFile(path, "w") as archive:
            archive.writestr("S1/fastqc_data.txt", FASTQC)
        assets.append({"asset_id": phase, "kind": "sequencing_qc", "path": filename,
                       "metadata": {"sample_id": "S1", "qc_phase": phase, "upstream_sha256": file_sha256(path)}})
    path = root / "S1_1.fastq.gz_trimming_report.txt"
    path.write_text("Total reads processed: 1,000\nReads with adapters: 250 (25.0%)\n"
                    "Reads written (passing filters): 900 (90.0%)\n")
    assets.append({"asset_id": "trim", "kind": "sequencing_qc", "path": path.name,
                   "metadata": {"sample_id": "S1", "upstream_sha256": file_sha256(path)}})
    path = root / "S1Log.final.out"
    path.write_text("Number of input reads | 1000\nUniquely mapped reads % | 79.66%\n"
                    "% of reads mapped to multiple loci | 10%\n% of reads mapped to too many loci | 1%\n"
                    "% of reads unmapped: too many mismatches | 1%\n% of reads unmapped: too short | 7%\n"
                    "% of reads unmapped: other | 1.34%\nAverage input read length | 297\n")
    assets.append({"asset_id": "star", "kind": "sequencing_qc", "path": path.name,
                   "metadata": {"sample_id": "S1", "upstream_sha256": file_sha256(path)}})
    project = {"schema_version": "2.0", "project_id": "qc-study", "name": "Native QC fixture",
               "input_root": ".", "output_root": "../unused-output", "annotation": {"genome_build": "test"},
               "samples": [{"sample_id": "S1", "biological_sample_id": "S1", "condition": "A"}],
               "contrasts": [], "assets": assets}
    return write_json(root / "project.json", project)


def test_sequencing_stage_and_views_keep_phases_native_values_and_transport(tmp_path):
    project = study(tmp_path / "original")
    spec = {"adapter": "qc.sequencing", "inputs": {"project": str(project)}}
    before = {path.name: (file_sha256(path), path.stat().st_mtime_ns)
              for path in project.parent.iterdir() if path.is_file()}
    assert readiness(spec) == []
    result = execute_stage(spec, tmp_path / "compute")
    assert result["status"] == "complete"
    assert result["evidence"]["upstream_jobs"] == 0
    assert result["evidence"]["raw_reads_or_bams_read"] is False
    assert result["evidence"]["native_qc_files"] == 4
    assert result["evidence"]["native_qc_figures"] >= 10
    assert {path.name: (file_sha256(path), path.stat().st_mtime_ns)
            for path in project.parent.iterdir() if path.is_file()} == before
    root = tmp_path / "relocated-collection"
    root.mkdir()
    shutil.copyfile(project, root / "project.json")
    transport = []
    for record in result["outputs"]:
        source = Path(record["path"])
        destination = root / "products" / source.relative_to(tmp_path / "compute")
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, destination)
        transport.append({"source": record, "collection_path": destination.relative_to(root).as_posix()})
    receipt = {"status": "completed", "outputs": result["outputs"], "execution_evidence": result["evidence"]}
    receipt["receipt_sha256"] = canonical_hash(receipt)
    plan = {"tasks": [{"id": "sequencing", "adapter": "qc.sequencing"}]}
    rendered = build_qc_views(root, plan, {"sequencing": (root / "receipt.json", receipt)}, transport)
    view = rendered["extra_views"][0]
    assert view["source_project_manifest_sha256"] == file_sha256(root / "project.json")
    assert len(view["panels"]) >= 10
    assert all(panel["chapter"] == 1 for panel in view["panels"])
    assert {panel["metadata"]["qc_phase"] for panel in view["panels"]} == {"raw", "trimmed", "trimming", "alignment"}
    assert not view.get("unavailable")
    assert all(Path(link["path"]).is_file() for panel in view["panels"] for link in panel["data_links"])
    # Rendering has not rounded or recomputed native quality observations.
    table = next(table for table in view["tables"] if table["name"] == "Raw | Complete original FastQC module table cells")
    assert "39.123456789" in Path(table["path"]).read_text()
    target = Path(view["panels"][0]["path"])
    target.write_bytes(target.read_bytes() + b"changed")
    with pytest.raises(ValueError, match="bytes changed"):
        build_qc_views(root, plan, {"sequencing": (root / "receipt.json", receipt)}, transport)


def test_changed_or_absent_native_sequencing_sources_do_not_pass_readiness(tmp_path):
    project = study(tmp_path / "original")
    spec = {"adapter": "qc.sequencing", "inputs": {"project": str(project)}}
    target = project.parent / "S1Log.final.out"
    target.write_text(target.read_text() + "changed\n")
    assert "SHA256 changed" in ";".join(readiness(spec))
    document = json.loads(project.read_text())
    document["assets"] = []
    write_json(project, document)
    assert "No completed native sequencing QC" in ";".join(readiness(spec))


def test_sequencing_qc_executes_through_task_and_relocatable_collection(tmp_path):
    from katana_next.flow.model import seal_plan
    from katana_next.flow.runtime import execute_task, validate_receipt
    from katana_next.flow.collect import collect_results, materialize_report_descriptor, validate_collection
    project = study(tmp_path / "original")
    original = project.read_bytes()
    plan = write_json(tmp_path / "plan.json", seal_plan({
        "project_manifest": str(project), "project_manifest_sha256": file_sha256(project),
        "report_profile": "test_subgraph", "tasks": [{"id": "sequencing", "workflow_id": "flow.qc.sequencing",
        "adapter": "qc.sequencing", "parameters": {}, "dependencies": [], "scope": {},
        "inputs": {"project": str(project)}}]}))
    receipt = execute_task(plan, "sequencing", tmp_path / "task")
    actual = validate_receipt(receipt)
    assert actual["status"] == "completed" and actual["execution_evidence"]["upstream_jobs"] == 0
    collection = collect_results(plan, [receipt], tmp_path / "collection")
    assert validate_collection(collection)["status"] == "completed"
    descriptor = materialize_report_descriptor(collection)
    views = descriptor["extra_views"]
    assert sum(len(view["panels"]) for view in views) >= 10
    moved = tmp_path / "different-machine-root"
    shutil.copytree(collection.parent, moved)
    relocated = materialize_report_descriptor(moved)
    assert all(Path(panel["path"]).is_relative_to(moved) for view in relocated["extra_views"] for panel in view["panels"])
    assert all(Path(link["path"]).is_file() for view in relocated["extra_views"] for panel in view["panels"] for link in panel["data_links"])
    assert project.read_bytes() == original


def test_each_provided_native_qc_phase_covers_eligible_declared_samples(tmp_path):
    project = study(tmp_path / "original")
    document = json.loads(project.read_text())
    document["samples"].append({"sample_id": "S2", "biological_sample_id": "S2", "condition": "A", "attributes": {"source": "fastq"}})
    write_json(project, document)
    spec = {"adapter": "qc.sequencing", "inputs": {"project": str(project)}}
    assert "missing declared eligible samples: S2" in ";".join(readiness(spec))
    document["samples"][-1]["attributes"]["source"] = "bam"
    write_json(project, document)
    assert readiness(spec) == []
    document["samples"][-1]["attributes"]["source"] = "genome_bam"
    write_json(project, document)
    assert readiness(spec) == []
