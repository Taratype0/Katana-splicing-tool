"""Portable target input acceptance; synthetic data, no scientific execution."""
from __future__ import annotations

import copy
import hashlib
import json
from pathlib import Path
import shutil
import stat
from types import SimpleNamespace
from zipfile import ZipFile, ZIP_DEFLATED, ZipInfo

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from katana_next.domain.models import AnnotationReference, Contrast, ContrastArm, ProjectManifest, Sample
from katana_next.infrastructure.fingerprints import fingerprint_file
from katana_next.infrastructure.result_store import ResultStore
from katana_next.targets.bundle import export_targets, open_bundle, validate_request, _request_id
from katana_next.workflows.base import ArtifactWriter


class Catalog:
    def __init__(self):
        self.items = {}

    def register_artifact(self, artifact):
        self.items[artifact.artifact_id] = artifact
        return artifact

    def get_artifact(self, artifact_id):
        return self.items.get(artifact_id)

    def list_artifacts(self, project_id, *, kind=None, analysis_id=None):
        return [a for a in self.items.values() if a.project_id == project_id
                and (kind is None or a.kind == kind)
                and (analysis_id is None or a.analysis_id == analysis_id)]

    def export_manifest(self, _):
        pass


@pytest.fixture
def study(tmp_path):
    source = tmp_path / "source-machine"
    source.mkdir()
    gtf = source / "genes.gtf"
    gtf.write_text("chr1\ttest\texon\t101\t200\t.\t+\t.\tgene_id \"G1\";\n")
    samples = []
    for group in "ABC":
        for rep in (1, 2):
            sid = f"{group}{rep}"
            bam, bai = source / f"{sid}.bam", source / f"{sid}.bam.bai"
            bam.write_bytes(("synthetic resource " + sid).encode())
            bai.write_bytes(b"index fixture")
            samples.append(Sample(sample_id=sid, biological_sample_id=f"mouse-{sid}", condition=group,
                                  attributes={"donor_id": f"donor-{sid}", "cell_type": "declared-bulk"},
                                  assets={"bam": bam, "bai": bai}))
    contrasts = [Contrast(contrast_id=f"{g}_vs_A", label=f"{g} versus A",
                         numerator=ContrastArm(label=g, sample_ids=(f"{g}1", f"{g}2")),
                         denominator=ContrastArm(label="A", sample_ids=("A1", "A2"))) for g in "BC"]
    manifest = ProjectManifest(project_id="toy", name="Toy", input_root=source, output_root=source / "results",
                               annotation=AnnotationReference(genome_build="toy-v1", gtf=gtf),
                               samples=tuple(samples), contrasts=tuple(contrasts), metadata={"assay_mode": "bulk_rna"})
    manifest_path = source / "project.json"
    manifest_path.write_text(manifest.model_dump_json())
    paths = [gtf, *(p for s in samples for p in s.assets.values())]
    context = SimpleNamespace(manifest=manifest, manifest_path=manifest_path,
                              snapshot=SimpleNamespace(manifest_fingerprint="snapshot-v1",
                                                       input_fingerprints=tuple(fingerprint_file(p) for p in paths)),
                              catalog=Catalog(), results=ResultStore(manifest.output_root))

    def add(name, frame, *, analysis="fixture.events", parameters=None, metadata=None, finish=True):
        run = hashlib.sha256(name.encode()).hexdigest()
        writer = ArtifactWriter(context, analysis, run, parameters or {}, "1.0", common_metadata=metadata or {})
        artifact = writer.table(name, "table.parquet", frame, entity_type="splice_event" if "event_id" in frame else "gene")
        if finish:
            writer.finish()
        return artifact

    rows = []
    for cid in ("B_vs_A", "C_vs_A"):
        for i, (gene, event) in enumerate((("G1", "SE:e1"), ("G1", "SE:e2"), ("G2", "SE:e3"))):
            rows.append(dict(row_id=f"{cid}::{event}", gene_id=gene, gene_symbol=f"Symbol{gene}",
                             event_id=event, event_type="SE", comparison_id=cid, chr="chr1", strand="+",
                             exonStart_0base=100 + i * 300, exonEnd=200 + i * 300,
                             upstreamES=10 + i * 300, upstreamEE=80 + i * 300,
                             downstreamES=220 + i * 300, downstreamEE=280 + i * 300,
                             dpsi=[0.2, -0.01, 0.4][i], FDR=[1e-200, 0.9, 0.04][i],
                             significant=i != 1, stable_significant=i != 1,
                             source_path=str(source / "native" / "SE.MATS.JC.txt")))
    frame = pd.DataFrame(rows)
    frame["exonEnd"] = frame.exonEnd.astype("int32")
    events = add("All measured events", frame, parameters={"count_mode": "JC", "fdr": .05, "min_abs_dpsi": .1},
                 metadata={"test_family": "all tested events in one comparison", "read_support_policy": "declared-fixture"})
    genes = add("Gene results", pd.DataFrame([dict(row_id="B_vs_A::G1", comparison_id="B_vs_A", gene_id="G1", gene_symbol="SymbolG1", padj=.001, log2FoldChange=1.25),
                                              dict(row_id="B_vs_A::G2", comparison_id="B_vs_A", gene_id="G2", gene_symbol="SymbolG2", padj=.1, log2FoldChange=.1)]), analysis="fixture.deg")
    psi = pd.DataFrame({"row_id": ["SE:e1", "SE:e2", "SE:e3"], "event_id": ["SE:e1", "SE:e2", "SE:e3"],
                        **{s.sample_id: pd.Series([.1, .2, .3], dtype="float32") for s in samples}})
    add("PSI matrix", psi, analysis="fixture.psi")
    add("Membership", pd.DataFrame({"row_id": ["m1", "m2", "m3"], "event_id": ["SE:e1", "SE:e2", "SE:e3"],
                                   "gene_id": ["G1", "G1", "G2"], "module": ["M01", "M02", "M03"]}), analysis="fixture.network")
    add("Module traits", pd.DataFrame([{"row_id": f"{m}:{c}", "module": m, "comparison_id": c, "FDR": f}
                                      for m, f in (("M01", .01), ("M02", .9), ("M03", .03))
                                      for c in ("B_vs_A", "C_vs_A")]), analysis="fixture.network_traits")
    return SimpleNamespace(context=context, add=add, frame=frame, events=events, genes=genes, source=source,
                           select_gene={"source_artifact_id": genes.artifact_id, "row_ids": ["B_vs_A::G1"], "target_kind": "gene",
                                        "filters": {"browse_settings": {"fdr": .02, "event_set": "stable_significant"}}},
                           select_event={"source_artifact_id": events.artifact_id, "row_ids": ["C_vs_A::SE:e3"], "target_kind": "event", "filters": {}})


def _export(study, tmp_path, *, selections=None, contrasts=None, name="portable"):
    return export_targets(study.context, selections or [study.select_gene], tmp_path / name, contrast_ids=contrasts)


def _table(bundle, name):
    item = next(e for e in bundle.request["evidence"] if e["name"] == name)
    return pq.read_table(bundle.root / item["path"]), item


def _rewrite_request(root, change):
    path = Path(root) / "request.json"
    request = json.loads(path.read_text())
    change(request)
    request["request_id"] = _request_id(request)
    path.write_text(json.dumps(request))
    return request


def test_multi_source_exact_gene_event_export_moves_without_source_catalog(study, tmp_path):
    before = {p: (_digest(p), p.stat().st_mtime_ns) for p in study.source.rglob("*") if p.is_file()}
    result = _export(study, tmp_path, selections=[study.select_gene, study.select_event])
    assert result["target_count"] == 2
    with open_bundle(result["request_path"]) as bundle:
        assert bundle.request["modules"] == ["gene_card", "event_card", "sashimi"]
        assert {t["target_kind"] for t in bundle.request["targets"]} == {"gene", "event"}
        event = next(t for t in bundle.request["targets"] if t["target_kind"] == "event")
        assert event["event_id"] == "SE:e3" and event["event_geometry"]["exonStart_0base"] == 700
        table, _ = _table(bundle, "All measured events")
        assert table.num_rows == 6  # nonsignificant G1:e2 is not removed
        assert table.schema.field("exonEnd").type == pa.int32()
        assert table.to_pandas().FDR.tolist() == study.frame.FDR.tolist()
        assert set(table.to_pandas().source_path.str.startswith("local-reference:")) == {True}
        request_text = json.dumps(bundle.request)
        assert str(study.source) not in request_text
        assert bundle.request["project"]["samples"][0]["biological_sample_id"] == "mouse-A1"
    assert {p: (_digest(p), p.stat().st_mtime_ns) for p in before} == before
    destination = tmp_path / "other-machine"
    destination.mkdir()
    moved_zip = destination / "new-name.zip"
    shutil.copyfile(result["archive_path"], moved_zip)
    shutil.rmtree(study.source)
    with open_bundle(moved_zip) as moved:
        temporary = moved.root
        assert _table(moved, "All measured events")[0].num_rows == 6
        assert validate_request(moved.request, moved.root)["valid"]
    assert not temporary.exists()


def _digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def test_gene_expands_every_measured_event_and_selected_sample_module_context(study, tmp_path):
    result = _export(study, tmp_path, contrasts=["B_vs_A"])
    with open_bundle(result["bundle_directory"]) as bundle:
        events, _ = _table(bundle, "All measured events")
        assert events["event_id"].to_pylist() == ["SE:e1", "SE:e2"]
        assert events["stable_significant"].to_pylist() == [True, False]
        psi, _ = _table(bundle, "PSI matrix")
        assert psi.num_rows == 2 and "C1" not in psi.column_names and "B1" in psi.column_names
        traits, _ = _table(bundle, "Module traits")
        assert traits["module"].to_pylist() == ["M01", "M02"]
        assert {s["sample_id"] for s in bundle.request["project"]["samples"]} == {"A1", "A2", "B1", "B2"}
        assert len(bundle.request["project"]["contrasts"]) == 1


def test_exact_event_does_not_expand_to_other_events_of_same_gene(study, tmp_path):
    selection = {**study.select_event, "row_ids": ["B_vs_A::SE:e1"]}
    result = _export(study, tmp_path, selections=[selection])
    with open_bundle(result["archive_path"]) as bundle:
        table, _ = _table(bundle, "All measured events")
        assert table["row_id"].to_pylist() == ["B_vs_A::SE:e1"]
        assert _table(bundle, "PSI matrix")[0]["event_id"].to_pylist() == ["SE:e1"]


def test_source_thresholds_and_browse_thresholds_remain_separate(study, tmp_path):
    selection = {**study.select_event, "row_ids": ["B_vs_A::SE:e1"],
                 "filters": {"source_parameters": {"fdr": .05}, "browse_settings": {"fdr": .001, "event_set": "significant"},
                             "path": str(study.source)}}
    result = _export(study, tmp_path, selections=[selection])
    with open_bundle(result["archive_path"]) as bundle:
        s = bundle.request["targets"][0]["source_selections"][0]
        assert s["source_parameters"]["fdr"] == .05
        assert s["browse_filters"]["browse_settings"]["fdr"] == .001
        assert s["count_mode"] == "JC"
        assert str(study.source) not in json.dumps(s)
        table, item = _table(bundle, "All measured events")
        assert item["metadata"]["test_family"] == "all tested events in one comparison"
        assert table["FDR"].to_pylist() == [1e-200]


def test_external_resources_never_enter_zip_and_hash_method_is_honest(study, tmp_path):
    result = _export(study, tmp_path)
    with ZipFile(result["archive_path"]) as z:
        assert all(not n.endswith((".bam", ".bai", ".gtf")) for n in z.namelist())
        assert "local-bindings.template.json" not in z.namelist()
    bindings = json.loads(Path(result["local_bindings_template"]).read_text())
    assert bindings["schema"] == "katana.target_bindings/v1"
    assert any(str(study.source) in v for v in bindings["resources"].values())
    with open_bundle(result["request_path"]) as bundle:
        for resource in bundle.request["resources"]:
            assert resource["verification"]["hash_method"] == "sha256-full"
            assert len(resource["sha256"]) == 64
        rid = bundle.request["project"]["samples"][0]["bam_resource_id"]
        assert rid in bindings["resources"]


@pytest.mark.parametrize("bad", ["foreign", "stale", "failed", "incomplete", "modified"])
def test_rejects_invalid_selected_sources(study, tmp_path, bad):
    artifact = study.genes
    if bad == "foreign":
        artifact = artifact.model_copy(update={"project_id": "elsewhere"})
    elif bad == "stale":
        artifact = artifact.model_copy(update={"metadata": {**artifact.metadata, "project_snapshot_fingerprint": "older"}})
    elif bad == "failed":
        artifact = artifact.model_copy(update={"status": "failed"})
    elif bad == "incomplete":
        (artifact.path.parent / "parameters.json").unlink()
    else:
        changed = pd.read_parquet(artifact.path)
        changed.loc[0, "padj"] = .999
        changed.to_parquet(artifact.path)
    study.context.catalog.items[artifact.artifact_id] = artifact
    with pytest.raises((ValueError, FileNotFoundError)):
        _export(study, tmp_path)
    assert not (tmp_path / "portable").exists()


@pytest.mark.parametrize("change", [lambda s: s.update(row_ids=[]), lambda s: s.update(row_ids=["missing"]),
                                     lambda s: s.update(row_ids=["B_vs_A::G1", "B_vs_A::G1"]),
                                     lambda s: s.update(target_kind="module")])
def test_invalid_selection_is_not_inferred_or_truncated(study, tmp_path, change):
    selection = copy.deepcopy(study.select_gene)
    change(selection)
    with pytest.raises(ValueError):
        _export(study, tmp_path, selections=[selection])


def test_explicit_contrasts_cannot_silently_discard_selected_rows(study, tmp_path):
    with pytest.raises(ValueError, match="outside"):
        _export(study, tmp_path, selections=[study.select_event], contrasts=["B_vs_A"])


def test_request_and_payload_tampering_and_unlisted_files_rejected(study, tmp_path):
    result = _export(study, tmp_path)
    root = Path(result["bundle_directory"])
    request = json.loads(Path(result["request_path"]).read_text())
    tampered = copy.deepcopy(request)
    tampered["project"]["name"] = "changed"
    with pytest.raises(ValueError, match="identity mismatch"):
        validate_request(tampered, root)
    member = root / next(iter(request["members"]))
    data = member.read_bytes()
    member.write_bytes(data + b"tampered")
    with pytest.raises(ValueError, match="integrity"):
        with open_bundle(root):
            pass
    member.write_bytes(data)
    (root / "extra.txt").write_text("unlisted")
    with pytest.raises(ValueError, match="member set"):
        with open_bundle(root):
            pass


@pytest.mark.parametrize("change", [
    lambda r: r["project"]["samples"][0].pop("biological_sample_id"),
    lambda r: r["project"]["contrasts"][0]["numerator"].update(sample_ids=["foreign"]),
    lambda r: r["targets"][0].update(comparison_ids=["undeclared"]),
    lambda r: r["resources"][0]["verification"].update(hash_method="sha256-sampled-3x1MiB"),
    lambda r: r["evidence"][0].update(path="../outside.parquet"),
])
def test_resigned_request_still_requires_valid_identity_semantics(study, tmp_path, change):
    result = _export(study, tmp_path)
    request = _rewrite_request(result["bundle_directory"], change)
    with pytest.raises(ValueError):
        validate_request(request, result["bundle_directory"])


@pytest.mark.parametrize("evil", ["../outside", "/absolute", "C:/machine", "evidence/../alias", "evidence\\alias", "CON", "trailing."])
def test_unsafe_zip_member_is_rejected_before_extraction(study, tmp_path, evil):
    result = _export(study, tmp_path)
    unsafe = tmp_path / "unsafe.zip"
    shutil.copyfile(result["archive_path"], unsafe)
    with ZipFile(unsafe, "a") as z:
        z.writestr(evil, "unsafe")
    with pytest.raises(ValueError):
        with open_bundle(unsafe):
            pass
    assert not (tmp_path / "outside").exists()


def test_duplicate_and_symlink_archives_rejected(study, tmp_path):
    result = _export(study, tmp_path)
    duplicate = tmp_path / "duplicate.zip"
    shutil.copyfile(result["archive_path"], duplicate)
    with ZipFile(duplicate, "a") as z:
        with pytest.warns(UserWarning):
            z.writestr("request.json", "{}")
    with pytest.raises(ValueError, match="Duplicate"):
        with open_bundle(duplicate):
            pass
    symlink = tmp_path / "symlink.zip"
    shutil.copyfile(result["archive_path"], symlink)
    with ZipFile(symlink, "a") as z:
        entry = ZipInfo("link.parquet")
        entry.create_system = 3
        entry.external_attr = (stat.S_IFLNK | 0o777) << 16
        z.writestr(entry, "../../outside")
    with pytest.raises(ValueError, match="regular-file"):
        with open_bundle(symlink):
            pass


def test_jc_and_jcec_families_are_retained_without_relabel_or_refit(study, tmp_path):
    other = study.add("JCEC measured", study.frame.assign(FDR=lambda d: d.FDR * .5), parameters={"count_mode": "JCEC"})
    result = _export(study, tmp_path, selections=[study.select_event, {**study.select_event, "source_artifact_id": other.artifact_id}])
    with open_bundle(result["archive_path"]) as bundle:
        target = bundle.request["targets"][0]
        assert len(target["source_selections"]) == 2
        assert {s["count_mode"] for s in target["source_selections"]} == {"JC", "JCEC"}
        assert _table(bundle, "All measured events")[0]["FDR"].to_pylist() == [.04]
        assert _table(bundle, "JCEC measured")[0]["FDR"].to_pylist() == [.02]
        assert bundle.request["execution"]["statistics_refit"] is False


def test_existing_exports_are_never_overwritten(study, tmp_path):
    result = _export(study, tmp_path)
    before = Path(result["archive_path"]).read_bytes()
    with pytest.raises(FileExistsError):
        _export(study, tmp_path)
    assert Path(result["archive_path"]).read_bytes() == before


def test_multi_gene_context_keeps_precise_per_target_event_rows(study, tmp_path):
    second = {**study.select_gene, "row_ids": ["B_vs_A::G2"]}
    result = _export(study, tmp_path, selections=[study.select_gene, second])
    with open_bundle(result["archive_path"]) as bundle:
        table, item = _table(bundle, "PSI matrix")
        for target in bundle.request["targets"]:
            own = table.take(item["metadata"]["target_row_indices"][target["target_id"]])
            assert own["event_id"].to_pylist() == (["SE:e1", "SE:e2"] if target["gene_id"] == "G1" else ["SE:e3"])
    _rewrite_request(result["bundle_directory"], lambda r: r["evidence"][0]["metadata"]["target_row_indices"].update({r["targets"][0]["target_id"]: [999999]}))
    with pytest.raises(ValueError, match="row indices"):
        with open_bundle(result["request_path"]):
            pass


def test_changed_project_manifest_is_not_exported_as_original_snapshot(study, tmp_path):
    payload = json.loads(study.context.manifest_path.read_text())
    payload["annotation"]["genome_build"] = "other-build"
    study.context.manifest_path.write_text(json.dumps(payload))
    with pytest.raises(ValueError, match="manifest changed"):
        _export(study, tmp_path)


@pytest.mark.parametrize("field,value", [("comparison_id", "foreign-comparison"), ("sample_id", "foreign-sample")])
def test_related_foreign_evidence_is_rejected(study, tmp_path, field, value):
    study.add("Contaminated gene evidence", pd.DataFrame([{"row_id": "outside", "gene_id": "G1", field: value}]))
    with pytest.raises(ValueError, match="foreign"):
        _export(study, tmp_path)


def test_legacy_declared_source_is_compatible_but_not_native_verified(study, tmp_path):
    old = study.genes.model_copy(update={"metadata": {
        "project_snapshot_fingerprint": study.context.snapshot.manifest_fingerprint,
        "workflow_version": "historical"}})
    study.context.catalog.items[old.artifact_id] = old
    (old.path.parent / "parameters.json").unlink()
    result = _export(study, tmp_path)
    with open_bundle(result["archive_path"]) as bundle:
        _, old_record = _table(bundle, "Gene results")
        assert old_record["metadata"]["completion_verification"] == "legacy_declared"
        assert "not recovered" in old_record["metadata"]["completion_caveat"]
        _, new_record = _table(bundle, "All measured events")
        assert new_record["metadata"]["completion_verification"] == "verified_completion_marker"


def test_nested_path_values_redacted_while_nested_measurements_and_missing_preserved(study, tmp_path):
    frame = pd.DataFrame({"row_id": ["nested"], "gene_id": ["G1"],
                          "source_list": [[str(study.source / "raw.tsv"), "clean"]],
                          "measurements": [[.2, None, .9]]})
    study.add("Nested result", frame)
    result = _export(study, tmp_path)
    with open_bundle(result["archive_path"]) as bundle:
        table, record = _table(bundle, "Nested result")
        assert table["source_list"][0].as_py()[0].startswith("local-reference:")
        assert table["measurements"][0].as_py() == [.2, None, .9]
        assert record["metadata"]["path_redactions"] == {"source_list": 1}


def test_completed_frozen_method_version_is_not_forced_to_current_policy(study, tmp_path):
    original = study.events
    marker_path = original.path.parent / "parameters.json"
    marker = json.loads(marker_path.read_text())
    marker["workflow_version"] = "historical-frozen-0.1"
    marker_path.write_text(json.dumps(marker))
    metadata = {**original.metadata, "workflow_version": "historical-frozen-0.1", "read_support_policy": "historical-original"}
    preserved = original.model_copy(update={"metadata": metadata})
    study.context.catalog.items[original.artifact_id] = preserved
    result = _export(study, tmp_path, selections=[study.select_event])
    with open_bundle(result["archive_path"]) as bundle:
        _, item = _table(bundle, "All measured events")
        assert item["metadata"]["workflow_version"] == "historical-frozen-0.1"
        assert item["metadata"]["read_support_policy"] == "historical-original"
        assert item["metadata"]["completion_verification"] == "verified_completion_marker"


@pytest.mark.parametrize("kind", ["gene", "event"])
def test_shared_event_locus_or_gene_symbol_cannot_alias_canonical_gene(study, tmp_path, kind):
    overlap = study.frame.iloc[[0]].copy()
    overlap = pd.concat([overlap, overlap.assign(row_id="other-gene", gene_id="G2")], ignore_index=True)
    study.add("Shared physical locus annotations", overlap)
    selection = study.select_gene if kind == "gene" else {**study.select_event, "row_ids": ["B_vs_A::SE:e1"]}
    result = _export(study, tmp_path, selections=[selection])
    with open_bundle(result["archive_path"]) as bundle:
        table, _ = _table(bundle, "Shared physical locus annotations")
        assert table["gene_id"].to_pylist() == ["G1"]
        # An exact event also carries its gene-level expression context.
        genes, _ = _table(bundle, "Gene results")
        assert genes["gene_id"].to_pylist() == ["G1"]
