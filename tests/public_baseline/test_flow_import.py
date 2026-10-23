import json
import importlib.util

import pytest

from katana_next.flow.importer import import_bundle
from katana_next.flow.model import file_sha256, write_json


def bundle(tmp_path):
    root = tmp_path / "upstream"
    root.mkdir()
    native = root / "SE.MATS.JCEC.txt"
    native.write_text("ID\tFDR\tIncLevelDifference\n1\t0.01\t-0.2\n", encoding="utf-8")
    manifest = {"schema": "katana-pipeline-bundle/v1", "status": "complete", "reference": {"genome_build": "test-build"},
                "samples": [{"sample_id": "a", "condition": "A", "subject": "donor-a"},
                            {"sample_id": "b", "condition": "B", "subject": "donor-b"}],
                "comparisons": [{"comparison_id": "B_vs_A", "treatment": "B", "control": "A",
                                 "rmats_group1": "A", "rmats_group2": "B", "deg_target": "B", "deg_reference": "A"}],
                "rmats_sample_order": {"A": ["a"], "B": ["b"]},
                "artifacts": [{"kind": "rmats", "relative_path": native.name, "comparison_id": "B_vs_A",
                               "event_type": "SE", "mode": "JCEC", "status": "available",
                               "bytes": native.stat().st_size, "sha256": file_sha256(native)}]}
    write_json(root / "manifest.json", manifest)
    template = write_json(tmp_path / "template.json", {"metadata": {"species": "test"}})
    return root, manifest, template


def test_import_retains_native_direction_and_independent_units(tmp_path):
    root, manifest, template = bundle(tmp_path)
    original = file_sha256(root / "SE.MATS.JCEC.txt")
    project = import_bundle(root, template, tmp_path / "workspace")
    value = json.loads(project.read_text("utf-8"))
    asset = value["assets"][0]
    assert asset["metadata"]["reverse_direction"] is True
    assert asset["metadata"]["numerator_sample_ids"] == ["b"]
    assert value["samples"][0]["biological_sample_id"] == "donor-a"
    assert value["input_root"] == "bundle"
    assert file_sha256(project.parent / "bundle" / asset["path"]) == original
    assert file_sha256(root / "SE.MATS.JCEC.txt") == original
    with pytest.raises(ValueError, match="exists"):
        import_bundle(root, template, project.parent)


def test_import_rejects_changed_native_results_and_bad_order(tmp_path):
    root, manifest, template = bundle(tmp_path)
    (root / "SE.MATS.JCEC.txt").write_text("changed", encoding="utf-8")
    with pytest.raises(ValueError, match="SHA256"):
        import_bundle(root, template, tmp_path / "bad")
    assert not (tmp_path / "bad").exists()


def test_import_rejects_ambiguous_rmats_sample_vectors(tmp_path):
    root, manifest, template = bundle(tmp_path)
    manifest["rmats_sample_order"]["A"] = ["a", "a"]
    write_json(root / "manifest.json", manifest)
    with pytest.raises(ValueError, match="vector order"):
        import_bundle(root, template, tmp_path / "bad")


def test_native_dexseq_csv_preserves_exons_and_separate_gene_screen(tmp_path):
    import csv
    from katana_next.infrastructure.manifest_store import ManifestStore
    from katana_next.workflows.bulk import _audit_asset
    root, manifest, template = bundle(tmp_path)
    members = [
        ('dexseq', 'perGeneQValue.B_vs_A.csv', 'groupID,padj\ng1,0.02\n'),
        ('auxiliary', 'DEXSeqResults.B_vs_A.csv',
         'groupID,featureID,pvalue,padj,log2fold_A_B\ng1,E001,1e-310,2e-310,-0.7\ng1,E002,NA,NA,0.5\n')]
    before = {}
    for kind, name, content in members:
        path = root / 'dexseq_exon' / name
        path.parent.mkdir(exist_ok=True)
        path.write_text(content)
        before[path] = (path.read_bytes(), path.stat().st_mtime_ns)
        manifest['artifacts'].append({'kind': kind, 'relative_path': 'dexseq_exon/' + name,
            'comparison_id': 'B_vs_A', 'status': 'available', 'bytes': path.stat().st_size,
            'sha256': file_sha256(path), 'native_fold_change_direction': 'control_over_treatment_native_header_authoritative'})
    write_json(root / 'manifest.json', manifest)
    project = import_bundle(root, template, tmp_path / 'imported')
    _, imported = ManifestStore().load(project)
    exon = next(a for a in imported.assets if a.kind == 'dexseq_result')
    gene = next(a for a in imported.assets if a.kind == 'dexseq_gene_result')
    assert _audit_asset(exon.kind, exon.path)[0] == 'pass'
    assert gene.metadata['native_table_projection']['unit'] == 'gene'
    actual = list(csv.reader(exon.path.open(), delimiter='\t'))
    expected = list(csv.reader((root / 'dexseq_exon/DEXSeqResults.B_vs_A.csv').open()))
    assert actual == expected  # no float roundtrip, FDR underflow or NA substitution
    assert before == {p: (p.read_bytes(), p.stat().st_mtime_ns) for p in before}
    assert (project.parent / 'bundle/dexseq_exon/DEXSeqResults.B_vs_A.csv').read_bytes() == before[root / 'dexseq_exon/DEXSeqResults.B_vs_A.csv'][0]
    assert len(json.loads((project.parent / 'import.receipt.json').read_text())['format_projections']) == 1


def recorded_metadata(root, manifest, **overrides):
    record = {"reference_declarations": {"genome_build": "test-build", "annotation_release": "release-7",
               "validation": "user_declared_labels_not_independently_verified"},
              "gene_exclusions": {"exclude_pseudogenes": True, "exclude_mitochondrial": True,
                                   "application_stage": "before_quantification_counts_and_differential_models"}}
    record.update(overrides)
    path = write_json(root / "provenance/source_map.json", record)
    manifest["artifacts"].append({"kind": "provenance", "relative_path": "provenance/source_map.json",
        "status": "available", "bytes": path.stat().st_size, "sha256": file_sha256(path)})
    write_json(root / "manifest.json", manifest)
    return record


def test_legacy_verified_source_map_restores_reference_and_filter_provenance(tmp_path):
    root, manifest, template = bundle(tmp_path)
    manifest["reference"] = {}
    record = recorded_metadata(root, manifest)
    before = {str(p): (file_sha256(p), p.stat().st_mtime_ns) for p in root.rglob("*") if p.is_file()}
    result = json.loads(import_bundle(root, template, tmp_path / "imported").read_text())
    assert result["annotation"]["genome_build"] == "test-build"
    assert result["metadata"]["annotation_release"] == "release-7"
    assert result["metadata"]["upstream_metadata"]["gene_exclusions"] == record["gene_exclusions"]
    assert result["metadata"]["upstream_metadata"]["reference_declarations"]["validation"] == record["reference_declarations"]["validation"]
    assert result["assets"][0]["metadata"]["upstream_artifact"] == manifest["artifacts"][0]
    assert before == {str(p): (file_sha256(p), p.stat().st_mtime_ns) for p in root.rglob("*") if p.is_file()}


@pytest.mark.parametrize("mismatch", ["reference", "manifest", "template", "unverified"])
def test_conflicting_or_unverified_upstream_metadata_is_not_imported(tmp_path, mismatch):
    root, manifest, template = bundle(tmp_path)
    recorded_metadata(root, manifest)
    if mismatch == "reference":
        manifest["reference"]["genome_build"] = "other-build"
    elif mismatch == "manifest":
        manifest["gene_exclusions"] = {"exclude_pseudogenes": False}
    elif mismatch == "template":
        write_json(template, {"metadata": {"annotation_release": "other-release"}})
    else:
        (root / "provenance/source_map.json").write_text('{}')
    write_json(root / "manifest.json", manifest)
    with pytest.raises(ValueError, match="[Cc]onflict|SHA256"):
        import_bundle(root, template, tmp_path / "rejected")
    assert not (tmp_path / "rejected").exists()


def test_actual_exporter_importer_preserves_declarations_without_claiming_validation(tmp_path):
    from katana_next.flow.distribution import export_pipeline
    pipeline = export_pipeline(tmp_path / "pipeline")
    spec = importlib.util.spec_from_file_location("bundle_exporter", pipeline / "bin/katana_bundle.py")
    exporter = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(exporter)
    artifacts = []
    for event in exporter.EVENTS:
        for mode in exporter.MODES:
            name = f"{event}.MATS.{mode}.txt"
            path = tmp_path / name
            path.write_text("\t".join(sorted(exporter.RMATS_COLUMNS | exporter.COORDINATES[event])) + "\n")
            artifacts.append({"source": str(path), "relative_path": f"rmats/B_vs_A/rmats_post/{name}",
                              "kind": "rmats", "comparison_id": "B_vs_A"})
    for kind, name in (("deg", "B_vs_A.deseq2.results.tsv"), ("normalized_counts", "all.normalised_counts.tsv")):
        path = tmp_path / name
        columns = sorted(exporter.TABLE_COLUMNS[kind]) + (["a", "b"] if kind == "normalized_counts" else [])
        path.write_text("\t".join(columns) + "\n")
        artifacts.append({"source": str(path), "relative_path": f"expression/{name}", "kind": kind,
                          **({"comparison_id": "B_vs_A"} if kind == "deg" else {})})
    samples = tmp_path / "samples.csv"
    samples.write_text("sample,condition\na,A\nb,B\n")
    contrasts = tmp_path / "contrasts.csv"
    contrasts.write_text("contrast,treatment,control\nB_vs_A,B,A\n")
    declarations = {"genome_build": "declared-build", "annotation_release": "declared-release",
                    "validation": "user_declared_labels_not_independently_verified"}
    exclusions = {"exclude_pseudogenes": True, "exclude_mitochondrial": True}
    source = write_json(tmp_path / "source.json", {"artifacts": artifacts, "reference_declarations": declarations,
        "gene_exclusions": exclusions, "assay": "bulk", "assay_validation": "declared",
        "rmats_sample_order": {"A": ["a"], "B": ["b"]}})
    exported = tmp_path / "exported"
    manifest = exporter.export_bundle(source, exported, samples, contrasts)
    assert manifest["reference_declarations"] == declarations
    assert manifest["gene_exclusions"] == exclusions
    assert manifest["assay"] == "bulk"
    template = write_json(tmp_path / "template.json", {"metadata": {"species": "test"}})
    project = json.loads(import_bundle(exported, template, tmp_path / "imported").read_text())
    assert project["annotation"]["genome_build"] == "declared-build"
    assert project["metadata"]["upstream_metadata"]["reference_declarations"] == declarations
    assert project["metadata"]["upstream_metadata"]["gene_exclusions"] == exclusions
    assert project["metadata"]["annotation_release"] == "declared-release"
