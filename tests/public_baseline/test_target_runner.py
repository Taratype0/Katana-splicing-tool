"""Portable target execution acceptance with real small indexed BAMs."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import shutil
import zipfile

import pandas as pd
import pytest

from katana_next.targets.bundle import _request_id
from katana_next.targets.runner import inspect_results, plan_targets, run_targets


def _write(path, value):
    path.write_text(json.dumps(value, indent=2), encoding="utf-8")


def _sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _fixture(root, *, sashimi=False, failed_target=False):
    source = root / "original-machine"
    source.mkdir()
    bundle = source / "portable"
    bundle.mkdir()
    (bundle / "evidence").mkdir()
    actual_resources = root / "execution-machine" / "local-data"
    actual_resources.mkdir(parents=True)
    resources, bindings, samples = [], {}, []
    for group in ("A", "B"):
        for rep in (1, 2):
            sid = f"{group}{rep}"
            sample = {"sample_id": sid, "biological_sample_id": f"mouse-{sid}", "condition": group,
                      "attributes": {"donor_id": f"donor-{sid}"}}
            if sashimi:
                from katana_next.sashimi_fixture import write_indexed_bam
                # Same basename across directories tests that track identities
                # come from declared biological samples rather than file stems.
                directory = actual_resources / sid
                directory.mkdir()
                bam = directory / "alignment.bam"
                write_indexed_bam(bam)
                bai = Path(str(bam) + ".bai")
                nonadjacent_index = directory / "separately-bound-index.bai"
                bai.rename(nonadjacent_index)
                for kind, path in (("bam", bam), ("bai", nonadjacent_index)):
                    rid = f"{sid}-{kind}"
                    sample[f"{kind}_resource_id"] = rid
                    resources.append({"resource_id": rid, "kind": kind, "sample_id": sid,
                                      "sha256": _sha(path), "bytes": path.stat().st_size,
                                      "verification": {"hash_method": "sha256-full", "digest": _sha(path)}})
                    bindings[rid] = str(path)
            samples.append(sample)
    targets = []
    for gene, kind, event in (("G1", "gene", None), ("G2", "event", "SE:e5")):
        target = {"target_id": kind + "-" + gene, "target_kind": kind, "gene_id": gene, "gene_symbol": "Symbol" + gene,
                  "comparison_ids": ["B_vs_A"], "source_selections": [{"row_ids": ["selected-" + gene],
                    "source_sha256": "b" * 64, "source_artifact_id": "native-artifact", "source_fingerprint": "frozen-source-native",
                    "source_parameters": {"fdr": .05, "count_mode": "JC"},
                    "browse_filters": {"fdr": .01, "event_set": "RS-DAS"}, "count_mode": "JC"}]}
        if event:
            target.update(event_id=event, event_type="SE", event_geometry={"chr": "chr1", "strand": "+",
                          "upstreamES": 100, "upstreamEE": 120, "exonStart_0base": 220, "exonEnd": 240,
                          "downstreamES": 320, "downstreamEE": 340})
        targets.append(target)
    if failed_target:
        targets.append({**targets[0], "target_id": "gene-missing", "gene_id": "NO_SUCH_GENE"})
    rows = []
    for i in range(1, 6):
        rows.append({"row_id": "selected-G1" if i == 1 else "selected-G2" if i == 5 else f"unselected-{i}",
                     "gene_id": "G1" if i < 5 else "G2", "gene_symbol": "SymbolG1" if i < 5 else "SymbolG2",
                     "event_id": f"SE:e{i}", "event_type": "SE", "comparison_id": "B_vs_A",
                     "chr": "chr1", "strand": "+", "upstreamES": 100, "upstreamEE": 120,
                     "exonStart_0base": 220, "exonEnd": 240, "downstreamES": 320, "downstreamEE": 340,
                     "FDR": 1e-200 if i == 1 else .9, "PValue": 3e-205 if i == 1 else .4,
                     "IncLevelDifference": .75, "IncLevel1": "0.9,0.8", "IncLevel2": "0.1,0.1",
                     "IJC_SAMPLE_1": "10,11", "SJC_SAMPLE_1": "1,1", "IJC_SAMPLE_2": "1,2", "SJC_SAMPLE_2": "10,10"})
    frame = pd.DataFrame(rows)
    frame["exonEnd"] = frame["exonEnd"].astype("int32")
    table = bundle / "evidence" / "native.parquet"
    frame.to_parquet(table, index=False)
    member = {"sha256": _sha(table), "bytes": table.stat().st_size}
    evidence = {"evidence_id": "completed-native", "analysis_id": "differential.as", "name": "Complete measured events",
                "entity_type": "splice_event", "path": "evidence/native.parquet", **member,
                "source_sha256": "b" * 64, "source_fingerprint": "frozen-source-native",
                "parameters": {"FDR": .05, "count_mode": "JC"}, "metadata": {"count_mode": "JC", "source_artifact_id": "native-artifact",
                    "target_row_indices": {t["target_id"]: list(range(4)) if t["target_kind"] == "gene" else [4] for t in targets}},
                "target_ids": [t["target_id"] for t in targets]}
    request = {"schema": "katana.target_request/v1", "created_at": "2026-10-23T00:00:00Z",
               "project": {"project_id": "portable-toy", "name": "Generated software acceptance fixture",
                           "assay": "bulk", "source_manifest_sha256": "c" * 64, "source_snapshot_fingerprint": "snapshot",
                           "references": {"genome_build": "synthetic", "coordinate_system": "0-based-half-open"},
                           "samples": samples, "contrasts": [{"contrast_id": "B_vs_A", "label": "B versus A",
                               "numerator": {"label": "B", "sample_ids": ["B2", "B1"]},
                               "denominator": {"label": "A", "sample_ids": ["A1", "A2"]}}]},
               "targets": targets, "evidence": [evidence], "resources": resources,
               "modules": ["gene_card", "event_card"] + (["sashimi"] if sashimi else []),
               "members": {"evidence/native.parquet": member},
               "execution": {"statistics_refit": False, "implicit_upstream_analysis": False}}
    request["request_id"] = _request_id(request)
    _write(bundle / "request.json", request)
    archive = root / "execution-machine" / "moved-targets.zip"
    with zipfile.ZipFile(archive, "w") as stream:
        for p in bundle.rglob("*"):
            if p.is_file():
                stream.write(p, p.relative_to(bundle).as_posix())
    # The destination has no original project DB or original bundle directory.
    shutil.rmtree(source)
    return archive, bindings, frame


def test_moved_bundle_all_gene_events_exact_event_native_statistics_cache(tmp_path):
    request, bindings, original = _fixture(tmp_path)
    output = tmp_path / "reports"
    planned = plan_targets(request, bindings)
    assert planned["ready"] and [t["event_count"] for t in planned["tasks"]] == [4, 1]
    result = run_targets(request, output, bindings, workers=2)
    assert result["status"] == "complete" and result["verified"]
    gene = next(t for t in result["targets"] if t["target_kind"] == "gene")
    event = next(t for t in result["targets"] if t["target_kind"] == "event")
    assert gene["event_count"] == 4 and event["event_count"] == 1
    saved = pd.read_parquet((output / gene["report"]).parent / gene["evidence"][0]["parquet"])
    pd.testing.assert_frame_equal(saved.reset_index(drop=True), original.iloc[:4].reset_index(drop=True))
    assert (saved["FDR"] == .9).sum() == 3 and saved["FDR"].iloc[0] == 1e-200
    samples = pd.read_parquet((output / event["report"]).parent / "sample-measurements.parquet")
    assert samples["sample_id"].tolist() == ["B2", "B1", "A1", "A2"]
    assert samples["value"].tolist() == [.9, .8, .1, .1]
    products = {p: p.stat().st_mtime_ns for p in output.rglob("*") if p.is_file()}
    repeated = run_targets(request, output, bindings, workers=2)
    assert repeated["status"] == "cache_hit" and repeated["cache"]["hits"] == 2
    assert all(p.stat().st_mtime_ns == value for p, value in products.items())
    assert inspect_results(output / "result-manifest.json")["verified"]
    with zipfile.ZipFile(result["archive"]) as stream:
        assert stream.testzip() is None and "INDEX.html" in stream.namelist()
        assert not any(name.endswith((".bam", ".bai")) for name in stream.namelist())


def test_subset_then_remaining_target_preserves_completed_products(tmp_path):
    request, bindings, _ = _fixture(tmp_path)
    output = tmp_path / "reports"
    first = run_targets(request, output, bindings, target_ids=["gene-G1"], workers=1)
    assert first["status"] == "partial"
    assert first["execution_status"] == "complete" and first["requested_target_ids"] == ["gene-G1"]
    target_root = (output / first["targets"][0]["report"]).parent
    products = {p: p.stat().st_mtime_ns for p in target_root.rglob("*") if p.is_file()}
    completed = run_targets(request, output, bindings, target_ids=["event-G2"], workers=1)
    assert completed["status"] == "complete"
    assert all(p.stat().st_mtime_ns == stamp for p, stamp in products.items())


def test_one_target_failure_does_not_abort_other_complete_targets(tmp_path, monkeypatch):
    request, bindings, _ = _fixture(tmp_path)
    from katana_next.targets import runner
    original = runner._tables_for_target
    def unavailable_one(request, root, target):
        if target["target_id"] == "event-G2":
            raise RuntimeError("Independent target read failed")
        return original(request, root, target)
    monkeypatch.setattr(runner, "_tables_for_target", unavailable_one)
    result = run_targets(request, tmp_path / "reports", bindings, workers=3)
    assert result["status"] == "partial"
    assert [t["status"] for t in result["targets"]] == ["complete", "failed"]
    assert "Independent target read failed" in result["targets"][1]["reason"]


def test_product_tampering_rejected_without_overwriting(tmp_path):
    request, bindings, _ = _fixture(tmp_path)
    output = tmp_path / "reports"
    result = run_targets(request, output, bindings, workers=1)
    report = output / result["targets"][0]["report"]
    report.write_text("modified independently", encoding="utf-8")
    with pytest.raises(ValueError, match="hash mismatch"):
        inspect_results(output)
    changed = run_targets(request, output, bindings, target_ids=["gene-G1"], workers=1)
    assert changed["status"] == "partial" and changed["targets"][0]["status"] == "failed"
    assert report.read_text(encoding="utf-8") == "modified independently"


def test_missing_and_corrupt_bam_are_unavailable_not_zero_read_tracks(tmp_path):
    request, bindings, _ = _fixture(tmp_path, sashimi=True)
    bad = Path(bindings["A1-bam"])
    bad.write_bytes(b"not the exported alignment")
    planned = plan_targets(request, bindings)
    assert planned["ready"]
    assert not planned["resources"]["A1-bam"]["available"]
    result = run_targets(request, tmp_path / "reports", bindings, target_ids=["event-G2"], workers=1)
    event = next(t for t in result["targets"] if t["target_id"] == "event-G2")
    assert event["status"] == "partial"
    assert event["sashimi"][0]["status"] == "unavailable"
    assert not event["sashimi"][0]["files"]


def test_actual_indexed_bam_separate_index_and_sample_order(tmp_path):
    request, bindings, original = _fixture(tmp_path, sashimi=True)
    before = {p: _sha(Path(p)) for p in bindings.values()}
    result = run_targets(request, tmp_path / "reports", bindings, target_ids=["event-G2"], workers=1)
    event = next(t for t in result["targets"] if t["target_id"] == "event-G2")
    assert event["status"] == "complete", event
    proof = event["sashimi"][0]
    assert proof["status"] == "complete" and proof["reader"] == "bamnostic"
    assert proof["sample_ids"] == [["B2", "B1"], ["A1", "A2"]]
    root = (Path(result["report"]).parent / event["report"]).parent
    counts_path = next(p for p in proof["files"] if p.endswith("native_sashimi_junction_counts.tsv"))
    counts = pd.read_csv(root / "sashimi" / counts_path, sep="\t")
    assert counts["sample"].drop_duplicates().tolist() == ["B2", "B1", "A1", "A2"]
    assert counts["junction_count_anchor8"].sum() > 0
    assert counts["FDR"].eq(.9).all()
    assert len([p for p in proof["files"] if p.endswith(".png")]) == 5
    assert all(Path(p).read_bytes() and _sha(Path(p)) == digest for p, digest in before.items())


def test_unknown_claim_not_removed_or_recomputed(tmp_path):
    request, bindings, _ = _fixture(tmp_path)
    output = tmp_path / "reports"
    run_targets(request, output, bindings, target_ids=["gene-G1"], workers=1)
    claim = output / ".claims" / "event-G2.json"
    claim.write_text('{"pid":98765,"identity":"unknown-live-execution"}', encoding="utf-8")
    original = claim.read_bytes()
    result = run_targets(request, output, bindings, target_ids=["event-G2"], workers=1)
    assert result["targets"][1]["status"] == "failed"
    assert claim.read_bytes() == original
    assert not (output / "targets" / "event-G2").exists()


def test_unowned_nonempty_output_is_never_overwritten(tmp_path):
    request, bindings, _ = _fixture(tmp_path)
    output = tmp_path / "reports"
    output.mkdir()
    stranger = output / "unrelated.txt"
    stranger.write_text("not Katana-owned", encoding="utf-8")
    with pytest.raises(ValueError, match="ownership"):
        run_targets(request, output, bindings)
    assert stranger.read_text(encoding="utf-8") == "not Katana-owned"


def test_actual_export_fixture_retains_every_analysis_and_target_scope(tmp_path):
    from target_runner_fixture import create_target_runner_fixture
    fixture = create_target_runner_fixture(tmp_path, sashimi=False)
    moved = tmp_path / "moved.zip"
    shutil.copyfile(fixture["export"]["archive_path"], moved)
    shutil.rmtree(fixture["source_dir"])
    result = run_targets(moved, tmp_path / "reports", fixture["bindings"], workers=2)
    assert result["status"] == "complete"
    assert len(result["targets"]) == 2
    gene = next(t for t in result["targets"] if t["target_kind"] == "gene")
    event = next(t for t in result["targets"] if t["target_kind"] == "event")
    assert gene["event_count"] == 4 and event["event_count"] == 1
    assert len(gene["evidence"]) == len(event["evidence"]) == 7
    for target in (gene, event):
        directory = Path(result["report"]).parent / Path(target["report"]).parent
        matrix = next(e for e in target["evidence"] if e["name"] == "Complete sample PSI matrix")
        frame = pd.read_parquet(directory / matrix["parquet"])
        assert set(frame["event_id"]) == ({f"SE:e{i}" for i in range(1, 5)} if target["target_kind"] == "gene" else {"SE:e5"})
        motifs = next(e for e in target["evidence"] if "motif sites" in e["name"])
        assert pd.read_parquet(directory / motifs["parquet"])["FDR"].eq(.03123456789012345).all()


def test_modes_reverse_direction_and_biological_units_remain_separate(tmp_path):
    from katana_next.targets.render import sample_measurements, render_sample_plots
    from katana_next.targets.runner import _sashimi_events
    request, _, _ = _fixture(tmp_path)
    from katana_next.targets.bundle import open_bundle
    with open_bundle(request) as bundle:
        doc = bundle.request
        target = doc["targets"][1]
    base = {"event_id": "SE:e5", "comparison_id": "B_vs_A", "IncLevel1": ".9,.8", "IncLevel2": ".1,.2", "reverse_direction": True}
    tables = []
    for mode, level, fdr in (("JC", ".1,.2", .01), ("JCEC", ".3,.4", .4)):
        item = {"evidence_id": mode, "source_fingerprint": "whole-family-" + mode, "metadata": {"count_mode": mode}}
        tables.append((item, pd.DataFrame([{**base, "IncLevel2": level, "FDR": fdr}])))
    samples = sample_measurements(doc, target, tables)
    assert len(samples) == 8
    assert samples.loc[samples["count_mode"] == "JC", "value"].tolist() == [.1, .2, .9, .8]
    assert samples.loc[samples["count_mode"] == "JCEC", "value"].tolist() == [.3, .4, .9, .8]
    assert samples["biological_sample_id"].str.startswith("mouse-").all()
    assert samples["donor_id"].str.startswith("donor-").all()
    plots = render_sample_plots(tmp_path, samples)
    assert len(plots) == 2 and {p["count_mode"] for p in plots} == {"JC", "JCEC"}
    card = {"event_id": "SE:e5", "event_type": "SE", "comparison_id": "B_vs_A", "event_geometry": target["event_geometry"],
            "evidence": [{**item, "row": frame.iloc[0].to_dict(), "source_sha256": mode} for mode, (item, frame) in zip(("JC", "JCEC"), tables)]}
    events = _sashimi_events([card])
    assert len(events) == 2 and {e["count_mode"] for e in events} == {"JC", "JCEC"}
    assert {e["FDR"] for e in events} == {.01, .4}
    assert len({e["native_statistical_family"] for e in events}) == 2


def test_transient_sashimi_failure_gets_independent_retry(tmp_path, monkeypatch):
    request, bindings, _ = _fixture(tmp_path, sashimi=True)
    from katana_next.targets import runner
    output = tmp_path / "reports"
    monkeypatch.setattr(runner, "generate_sashimi", lambda *a: [{"status": "failed", "reason": "Synthetic transient reader failure", "files": []}])
    first = run_targets(request, output, bindings, workers=1, target_ids=["event-G2"])
    failed = next(t for t in first["targets"] if t["target_id"] == "event-G2")
    assert failed["status"] == first["execution_status"] == "failed"
    failed_receipt = output / failed["manifest"]
    prior = failed_receipt.read_bytes()
    # An explicitly mocked successful read tests execution lifecycle, not BAM
    # correctness (the separate actual indexed-BAM test covers that).
    monkeypatch.setattr(runner, "generate_sashimi", lambda *a: [{"status": "complete", "files": []}])
    second = run_targets(request, output, bindings, workers=1, target_ids=["event-G2"])
    accepted = next(t for t in second["targets"] if t["target_id"] == "event-G2")
    assert accepted["status"] == "complete" and second["execution_status"] == "complete"
    assert accepted["manifest"] != failed["manifest"] and failed_receipt.read_bytes() == prior


def test_sampled_resource_identity_uses_declared_sampling_not_full_read(tmp_path, monkeypatch):
    from katana_next.targets import runner
    resource = tmp_path / "large-resource.bam"
    resource.write_bytes(b"a" * (4 * 1024 * 1024))
    declared = runner._legacy_sampled_sha256(resource)
    request = {"resources": [{"resource_id": "reads", "kind": "bam", "bytes": resource.stat().st_size,
                              "verification": {"hash_method": "sha256-sampled-3x1MiB", "digest": declared}}]}
    monkeypatch.setattr(runner, "sha256", lambda path: pytest.fail("Sampled resources must not trigger a full read"))
    result = runner._bound_resources(request, {"reads": resource})["reads"]
    assert result["available"]
    assert result["verification"]["hash_method"] == "sha256-sampled-3x1MiB"
    assert not result["verification"]["full_content_verified"]
    with resource.open("r+b") as stream:
        stream.write(b"changed sampled bytes")
    assert not runner._bound_resources(request, {"reads": resource})["reads"]["available"]


def test_reference_mismatch_never_produces_zero_coverage_success(tmp_path):
    from katana_next.targets.bundle import open_bundle
    from katana_next.targets.runner import _bound_resources
    from katana_next.targets.sashimi import generate_sashimi
    request, bindings, _ = _fixture(tmp_path, sashimi=True)
    with open_bundle(request) as bundle:
        doc = bundle.request
        target = dict(doc["targets"][1])
    target["event_geometry"] = {**target["event_geometry"], "chr": "absent-contig"}
    resources = _bound_resources(doc, {rid: Path(path) for rid, path in bindings.items()})
    records = generate_sashimi(doc, target, [{"event_id": "SE:e5", "comparison_id": "B_vs_A", "event_type": "SE"}],
                               resources, tmp_path / "sashimi")
    assert records[0]["status"] == "failed" and "absent" in records[0]["reason"]
    assert not records[0]["files"] and not list((tmp_path / "sashimi").rglob("*.png"))


def test_result_manifest_archive_escape_and_false_success_rejected(tmp_path):
    request, bindings, _ = _fixture(tmp_path)
    output = tmp_path / "reports"
    result = run_targets(request, output, bindings, target_ids=["gene-G1"], workers=1)
    path = Path(result["manifest"])
    original = json.loads(path.read_text("utf-8"))
    modified = {**original, "archive": "../foreign.zip"}
    _write(path, modified)
    with pytest.raises(ValueError, match="Unsafe portable"):
        inspect_results(output)
    _write(path, {**original, "status": "complete"})
    with pytest.raises(ValueError, match="actual target statuses"):
        inspect_results(output)


def test_independent_scheduler_targets_can_share_one_result_root(tmp_path):
    from concurrent.futures import ThreadPoolExecutor
    request, bindings, _ = _fixture(tmp_path)
    output = tmp_path / "reports"
    with ThreadPoolExecutor(max_workers=2) as pool:
        futures = [pool.submit(run_targets, request, output, bindings, target_ids=[target], workers=1)
                   for target in ("gene-G1", "event-G2")]
        results = [future.result() for future in futures]
    assert all(r["execution_status"] == "complete" for r in results)
    complete = inspect_results(output)
    assert complete["status"] == "complete"
    assert {t["target_id"] for t in complete["targets"]} == {"gene-G1", "event-G2"}


def test_canonical_event_and_gene_identifiers_do_not_match_foreign_aliases():
    from katana_next.targets.runner import _select_frame
    frame = pd.DataFrame([{"event_id": "canonical-one", "event_uid": "old-one", "gene_id": "G1", "gene_symbol": "Shared"},
                          {"event_id": "foreign", "event_uid": "canonical-one", "gene_id": "G2", "gene_symbol": "Shared"}])
    event = _select_frame(frame, {"target_kind": "event", "event_id": "canonical-one", "gene_id": "G1"})
    gene = _select_frame(frame, {"target_kind": "gene", "gene_id": "G1", "gene_symbol": "Shared"})
    assert event["event_id"].tolist() == gene["event_id"].tolist() == ["canonical-one"]


def test_single_cell_unfiltered_bam_is_not_labeled_cell_type_sashimi(tmp_path):
    from katana_next.targets.bundle import open_bundle
    from katana_next.targets.runner import _bound_resources
    from katana_next.targets.sashimi import generate_sashimi
    request, bindings, _ = _fixture(tmp_path, sashimi=True)
    with open_bundle(request) as bundle:
        doc = bundle.request
        target = doc["targets"][1]
    doc["project"]["assay"] = "scRNA"
    for sample in doc["project"]["samples"]:
        sample["attributes"]["cell_type"] = "Neuron"
    resources = _bound_resources(doc, {rid: Path(path) for rid, path in bindings.items()})
    records = generate_sashimi(doc, target, [{"event_id": "SE:e5", "comparison_id": "B_vs_A", "event_type": "SE"}],
                               resources, tmp_path / "sashimi")
    assert records[0]["status"] == "unavailable" and "cell-type-filtered" in records[0]["reason"]
    assert not list((tmp_path / "sashimi").rglob("*.png"))
