"""Adapter tests use real parsers/geometry; mocked tools are labelled explicitly."""
from dataclasses import replace
import hashlib
import json
import os
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import pytest

from katana_next.flow import stages


def test_annotation_bed_is_strand_oriented_and_preserves_source(tmp_path):
    gtf = tmp_path / "annotation.gtf"
    gtf.write_text('chr1\tx\texon\t11\t20\t.\t-\t.\tgene_id "g1"; transcript_id "t1";\n'
                   'chr1\tx\texon\t31\t45\t.\t-\t.\tgene_id "g1"; transcript_id "t1";\n')
    before = hashlib.sha256(gtf.read_bytes()).hexdigest()
    receipt = stages.execute_stage({"adapter": "reference.bed12", "inputs": {"gtf": str(gtf)}}, tmp_path / "result")
    bed = Path(receipt["outputs"][0]["path"])
    assert bed.read_text().strip().split("\t") == ["chr1", "10", "45", "t1", "0", "-", "10", "45", "0", "2", "10,15", "0,20"]
    assert receipt["status"] == "complete" and receipt["evidence"]["transcripts"] == 1
    assert hashlib.sha256(gtf.read_bytes()).hexdigest() == before
    assert receipt["outputs"][0]["sha256"] == hashlib.sha256(bed.read_bytes()).hexdigest()
    with pytest.raises(FileExistsError, match="Reconcile"):
        stages.execute_stage({"adapter": "reference.bed12", "inputs": {"gtf": str(gtf)}}, tmp_path / "result")


def test_overlapping_exons_are_not_silently_merged(tmp_path):
    gtf = tmp_path / "bad.gtf"
    gtf.write_text('chr1\tx\texon\t11\t30\t.\t+\t.\ttranscript_id "t1";\n'
                   'chr1\tx\texon\t25\t45\t.\t+\t.\ttranscript_id "t1";\n')
    with pytest.raises(ValueError, match="Overlapping exons"):
        stages.execute_stage({"adapter": "reference.bed12", "inputs": {"gtf": str(gtf)}}, tmp_path / "out")
    assert not (tmp_path / "out/stage.receipt.json").exists()


def test_iterative_plot_accepts_deep_unbalanced_original_linkage(tmp_path):
    n = 2500
    tree = np.array([[0, 1, .001, 2], *[[n+i-1, i+1, (i+1)/n, i+2] for i in range(1,n-1)]])
    original = tree.copy()
    output = stages._iterative_tree(tree, tmp_path / "deep.png", "observed")
    assert output.stat().st_size > 1000
    np.testing.assert_array_equal(tree, original)


def test_tree_is_persisted_before_graphics_failure_and_original_renderer_restored(tmp_path, monkeypatch):
    from katana_cosplice.network import NetworkResult
    import katana_cosplice.pipeline as pipeline
    tree = np.array([[0,1,.1,2], [3,2,.4,3]])
    matrix = pd.DataFrame([[.1,.2],[.2,.4],[.8,.7]], index=["e1","e2","e3"])
    result = NetworkResult(pd.DataFrame(), pd.DataFrame(), pd.DataFrame(),
                           pd.DataFrame(), pd.DataFrame(), matrix, tree)
    def broken(*_):
        raise RuntimeError("actual downstream graphics failure fixture")
    monkeypatch.setattr(pipeline, "render_network_plots", broken)
    with pytest.raises(RuntimeError, match="graphics failure"):
        with stages.persist_network_trees():
            pipeline.render_network_plots(result, pd.DataFrame(), tmp_path / "observed/plots", "observed")
    np.testing.assert_array_equal(np.load(tmp_path / "observed/event_linkage.npy", allow_pickle=False), tree)
    assert pd.read_csv(tmp_path / "observed/event_linkage_leaf_ids.tsv", sep="\t").event_key.tolist() == list(matrix.index)
    assert pipeline.render_network_plots is broken


def _bam(tmp_path):
    from katana_next.sashimi_fixture import write_indexed_bam
    bam = tmp_path / "input.bam"
    write_indexed_bam(bam)
    return bam


def test_qc_validates_real_bam_and_reads_mock_native_profile(tmp_path, monkeypatch):
    bam = _bam(tmp_path)
    bed = tmp_path / "ref.bed"
    bed.write_text("chr1\t0\t100\tt1\t0\t+\t0\t100\t0\t1\t100\t0\n")
    def mock_rseqc(command, out, name):
        # Transport fixture only: this is not reported as a native RSeQC test.
        (out / "gene_body.geneBodyCoverage.txt").write_text(
            "Percentile\t" + "\t".join(map(str, range(1,101))) + "\ninput\t" + "\t".join(["3"]*100) + "\n")
        (out / (name + ".command.json")).write_text(json.dumps({"fixture": "mock RSeQC", "argv": command}))
    monkeypatch.setattr(stages, "_run", mock_rseqc)
    result = stages.execute_stage({"adapter": "qc.genebody", "inputs": {
        "bam": str(bam), "reference_bed": str(bed)}, "parameters": {
        "rseqc_executable": sys.executable, "sample_id": "declared_mouse_1"}}, tmp_path / "qc")
    frame = pd.read_csv(tmp_path / "qc/gene_body_coverage.tsv", sep="\t")
    assert set(frame.sample_id) == {"declared_mouse_1"}
    assert len(frame) == 100 and frame.coverage.sum() == 300
    assert result["evidence"]["native_sample_label"] == "input"
    assert result["evidence"]["project_asset"]["kind"] == "gene_body_coverage"


def test_qc_rejects_reference_mismatch_before_tool_execution(tmp_path, monkeypatch):
    bam = _bam(tmp_path)
    bed = tmp_path / "ref.bed"
    bed.write_text("1\t0\t100\tt1\t0\t+\t0\t100\t0\t1\t100\t0\n")
    def forbidden(*_):
        pytest.fail("Reference mismatch should fail before launching RSeQC")
    monkeypatch.setattr(stages, "_run", forbidden)
    with pytest.raises(ValueError, match="reference chromosomes differ"):
        stages.execute_stage({"adapter": "qc.genebody", "inputs": {
            "bam": str(bam), "reference_bed": str(bed)}, "parameters": {
            "rseqc_executable": sys.executable, "sample_id": "mouse"}}, tmp_path / "qc")
    assert not (tmp_path / "qc/stage.receipt.json").exists()


def test_missing_native_mapp_is_a_preflight_failure_not_completed_import(tmp_path):
    issues = stages.readiness({"adapter": "mapp.maei", "inputs": {"native_root": str(tmp_path)},
                               "parameters": {"rscript": sys.executable}})
    assert any("official MAPP source" in issue for issue in issues)
    assert any("inclusion_table" in issue for issue in issues)
    assert any("source_sha256" in issue for issue in issues)
    assert not (tmp_path / "stage.receipt.json").exists()


def test_network_output_labels_match_the_existing_report_contract(tmp_path):
    output = tmp_path / "engine/04_network/observed/event_module_assignments.tsv"
    output.parent.mkdir(parents=True); output.write_text("event_key\tmodule\ne1\tblue\n")
    record = stages._ledger([output], "network.fit")[0]
    assert record["analysis_id"] == "programs.cosplice"
    assert record["artifact_key"] == "assignments"
    assert record["name"] == "4A. Observed network: assignments"


def _mapp_prepare_fixture(tmp_path):
    gtf = tmp_path / "annotation.gtf"
    gtf.write_text("".join(
        f'chr1\tfixture\texon\t{a+1}\t{b}\t.\t+\t.\tgene_id "G1"; transcript_id "{tx}";\n'
        for tx, blocks in (("INC.1", [(100,120),(220,340),(440,460)]),
                           ("EXC.1", [(100,120),(440,460)]))
        for a,b in blocks), encoding="utf-8")
    events = tmp_path / "events.tsv"
    pd.DataFrame([{"event_id": "SE:exact1", "event_type": "SE", "gene_id": "G1", "chr": "chr1",
                   "strand": "+", "upstreamEE": 120, "exonStart_0base": 220, "exonEnd": 340,
                   "downstreamES": 440, "comparison_id": "B_vs_A"},
                  {"event_id": "SE:no_skipped_tx", "event_type": "SE", "gene_id": "G1", "chr": "chr1",
                   "strand": "+", "upstreamEE": 120, "exonStart_0base": 220, "exonEnd": 340,
                   "downstreamES": 999, "comparison_id": "B_vs_A"}]).to_csv(events,sep="\t",index=False)
    fasta = tmp_path / "reference.fa"
    fasta.write_text(">chr1\n" + "ACGT"*250 + "\n",encoding="ascii",newline="\n")
    Path(str(fasta)+".fai").write_text("chr1\t1000\t6\t1000\t1001\n")
    quant = {}
    for sample, inclusion in (("A1",3),("A2",4),("B1",6),("B2",7)):
        path = tmp_path / (sample + ".quant.sf")
        pd.DataFrame({"Name": ["INC.1","EXC.1"], "TPM": [inclusion,10-inclusion]}).to_csv(path,sep="\t",index=False)
        quant[sample] = str(path)
    return {"adapter": "mapp.prepare", "inputs": {"events": str(events), "gtf": str(gtf),
              "genome_fasta": str(fasta), "quantification": quant},
            "parameters": {"sample_ids": list(quant), "comparison_id": "B_vs_A", "site": "3ss",
                           "window_start": -50, "window_end": 0}}


def test_mapp_bridge_keeps_exact_transcript_mapping_and_all_losses(tmp_path):
    from katana_next.flow.mapp_inputs import COUNTER, prepare_inputs
    spec = _mapp_prepare_fixture(tmp_path)
    native = tmp_path / "native"
    counter = native / COUNTER
    counter.parent.mkdir(parents=True); counter.write_text("# Transport fixture; does not impersonate the native algorithm\n")
    spec["inputs"].update(native_root=str(native), counter_sha256=stages._sha(counter))
    out = tmp_path / "prepare"; out.mkdir()
    def mock_counter(command, directory, name):
        # This test checks exact transcript identity, TPM and loss accounting.
        # The separate opt-in native test executes the actual official counter.
        pd.DataFrame(np.zeros((1,5440)),index=pd.Index(["SE:exact1"],name="coordinates"),
                     columns=[f"test-model-{i}" for i in range(5440)]).to_csv(directory / "NATIVE_SITECOUNTS.tsv",sep="\t")
        (directory / (name+".command.json")).write_text(json.dumps({"fixture": "mock transport"}))
    evidence, files = prepare_inputs(spec,out,run=mock_counter,sha=stages._sha,executable=stages._executable)
    inclusion = pd.read_csv(out / "INCLUSION_TABLE.tsv",sep="\t",index_col=0)
    assert inclusion.loc["SE:exact1","A1_included"] == 3
    assert inclusion.loc["SE:exact1","B2_included"] == 7
    assert inclusion.filter(like="_total").eq(10).all().all()
    audit = pd.read_csv(out / "all_event_eligibility.tsv",sep="\t").set_index("event_id")
    assert audit.loc["SE:no_skipped_tx","status"] == "no_exact_annotated_included_or_skipped_transcript"
    assert evidence["prepared_events"] == 1 and evidence["all_events"] == 2
    assert evidence["quantity_provenance"]["not_native_SUPPA_discovery"] is True
    assert (out / "sequences.fasta").read_text().splitlines()[1] == ("ACGT"*250)[170:220]


@pytest.mark.skipif(not os.environ.get("KATANA_TEST_MAPP_ROOT"), reason="Requires a separately installed official MAPP checkout")
def test_mapp_prepare_executes_the_actual_official_counter(tmp_path):
    from katana_next.flow.mapp_inputs import COUNTER
    spec = _mapp_prepare_fixture(tmp_path)
    native = Path(os.environ["KATANA_TEST_MAPP_ROOT"])
    spec["inputs"].update(native_root=str(native),counter_sha256=stages._sha(native / COUNTER))
    receipt = stages.execute_stage(spec,tmp_path / "actual-native-prepare")
    counts = pd.read_csv(tmp_path / "actual-native-prepare/NATIVE_SITECOUNTS.tsv",sep="\t",index_col=0)
    assert counts.shape == (1,5440)
    assert counts.loc["SE:exact1","u50to0.d0to0|ACGT"] == 12
    assert receipt["evidence"]["models"] == 5440
    assert receipt["status"] == "complete"


def test_preflight_deferred_paths_cannot_bypass_execution(tmp_path):
    spec = {"adapter": "reference.bed12", "inputs": {}, "readiness_deferred_inputs": ["gtf"]}
    assert stages.readiness(spec) == []
    with pytest.raises(ValueError,match="Resolve dependency"):
        stages.execute_stage(spec,tmp_path / "forbidden")
    assert not (tmp_path / "forbidden").exists()


def _network_interpretation_fixture(tmp_path):
    from katana_next.domain.models import ProjectManifest, Sample, Contrast, ContrastArm, InputAsset, AnnotationReference
    from katana_next.infrastructure.manifest_store import ManifestStore
    samples = [Sample(sample_id=f"{condition}{i}", biological_sample_id=f"{condition}{i}",condition=condition)
               for condition in ("A","B") for i in range(1,4)]
    gmt = tmp_path / "GO_BP.gmt"
    gmt.write_text("example_process\tfixture only\tG1\tG2\n")
    project = ProjectManifest(name="Generated network interpretation fixture",input_root=tmp_path,annotation=AnnotationReference(genome_build="test"),
        output_root=tmp_path / "source-results",samples=tuple(samples),assets=(InputAsset(asset_id="GO_BP",kind="gene_sets_gmt",path=gmt),),
        contrasts=(Contrast(contrast_id="B_vs_A",label="B versus A",
            numerator=ContrastArm(label="B",sample_ids=("B1","B2","B3")),
            denominator=ContrastArm(label="A",sample_ids=("A1","A2","A3"))),))
    manifest = ManifestStore().save(project,tmp_path / "project.json")
    engine = tmp_path / "engine"
    (engine / "04_network/observed").mkdir(parents=True)
    (engine / "02_matrix").mkdir()
    metadata = pd.DataFrame([{"sample_id":s.sample_id,"biological_sample_id":s.biological_sample_id,
                              "condition":s.condition} for s in samples])
    metadata.to_csv(engine / "02_matrix/sample_metadata_resolved.tsv",sep="\t",index=False)
    pd.DataFrame({"M1":[-.8,-.6,-.7,.6,.9,.7],"M2":[.1,-.2,.2,.2,.1,-.3]},
                 index=[s.sample_id for s in samples]).to_csv(engine / "04_network/observed/module_program_scores.tsv",sep="\t")
    pd.DataFrame({"event_key":["e1","e2","e3","e4"],"module":["M1","M1","M2","grey"]}).to_csv(
        engine / "04_network/observed/event_module_assignments.tsv",sep="\t",index=False)
    pd.DataFrame({"event_key":["e1","e2","e3","e4"],"gene_symbol":["G1","G2","G3","G4"]}).to_csv(
        engine / "02_matrix/event_metadata.tsv",sep="\t",index=False)
    return {"inputs":{"project":str(manifest)},"scope":{"contrast_ids":["B_vs_A"]}},engine


def test_network_interpretation_uses_exact_declared_arms_and_measured_background(tmp_path):
    spec,engine = _network_interpretation_fixture(tmp_path)
    associations,_ = stages._matched_network_associations(spec,engine,"observed")
    selected = pd.read_csv(associations["selected_tests"],sep="\t")
    full = pd.read_csv(associations["full_tests"],sep="\t")
    assert set(selected.contrast_id) == {"B_vs_A"} and len(selected) == 2
    assert len(full) > len(selected)
    assert set(full.correction_family_id) == {"B_vs_A::frozen_binary_trait_family"}
    audit = json.loads(Path(associations["trait_audit"]).read_text())
    assert audit[0]["sample_ids"]["numerator"] == ["B1","B2","B3"]
    enrichment,_ = stages._network_enrichment(spec,engine)
    background = pd.read_csv(enrichment["background"],sep="\t")
    assert set(background.gene_symbol) == {"G1","G2","G3","G4"}  # grey is measured background
    tested = pd.read_csv(enrichment["tested"],sep="\t")
    assert set(tested.module) == {"M1","M2"} and tested.background_genes.eq(4).all()


def test_network_continuation_requires_verified_dependency_engine_not_yaml_cache(tmp_path):
    from katana_cosplice.config import CoSpliceConfig
    from katana_next.workflows.programs import write_production_configuration
    samples = tmp_path / "samples.tsv"
    pd.DataFrame({"sample_id":["A","B","C"],"condition":["A","B","C"]}).to_csv(samples,sep="\t",index=False)
    old = tmp_path / "old-engine"; old.mkdir()
    source = tmp_path / "source.yaml"
    write_production_configuration(CoSpliceConfig(project_name="fixture",rmats_root=None,sample_metadata=samples,output_dir=old),source)
    with pytest.raises(ValueError,match="verified network.fit dependency"):
        stages.execute_stage({"adapter":"network.bootstrap","inputs":{"config":str(source)}},tmp_path / "bootstrap")
    assert list(old.iterdir()) == []


def test_auto_network_projection_preserves_native_orders_and_real_condition_labels(tmp_path):
    from katana_next.full_fixture import create_full_fixture
    from katana_next.flow.network_inputs import prepare_project_network,matrix_display_labels
    from katana_next.infrastructure.manifest_store import ManifestStore
    from katana_cosplice.pipeline import CoSplicePipeline
    original_manifest = create_full_fixture(tmp_path / "original")
    store = ManifestStore()
    _,project = store.load(original_manifest)
    # Exercise non-path-safe display conditions without changing sample IDs.
    rename={"Baseline":"WT/CONTROL", "Perturbation":"KO minus target", "Rescue":"Rescue"}
    project=project.model_copy(update={"samples":tuple(s.model_copy(update={"condition":rename[s.condition]}) for s in project.samples)})
    manifest=store.save(project,tmp_path / "display-project.json")
    original={a.path:stages._sha(a.path) for a in project.assets if a.kind=="rmats_result"}
    config,proof,_ = prepare_project_network(manifest,tmp_path / "prepared",{"rmats_mode":"JC"},stages._sha)
    pipeline=CoSplicePipeline(config)
    pipeline.run_stage("audit")
    with matrix_display_labels(proof["condition_aliases"]):
        pipeline.run_stage("matrix")
    observed=pd.read_csv(config.output_dir / "02_matrix/sample_metadata_resolved.tsv",sep="\t")
    assert set(observed.condition)==set(rename.values())
    assert set(observed.sample_id)=={s.sample_id for s in project.samples}
    mapping=json.loads(Path(proof["projection_mapping"]).read_text())
    assert len(mapping["projection"])==15 and not any(r["statistics_refitted"] for r in mapping["projection"])
    assert all(stages._sha(path)==digest for path,digest in original.items())


def test_raw_template_readiness_defers_only_genuine_future_inputs(tmp_path,monkeypatch):
    resource_template=tmp_path / "resource-template.json"
    resource_template.write_text(json.dumps({"metadata":{"rscript_path":sys.executable}}))
    seen=[]
    monkeypatch.setattr(stages,"_require_r_packages",lambda exe,packages,parameters:seen.append(packages))
    spec={"adapter":"network.fit","inputs":{"project":str(resource_template)},"readiness_deferred_inputs":["project"]}
    assert stages.readiness(spec)==[]
    assert seen==[("WGCNA",)]
    with pytest.raises(ValueError,match="Resolve dependency"):
        stages.execute_stage(spec,tmp_path / "blocked")


def test_dynamic_tree_ordinal_compat_changes_only_index_dtype_and_is_scoped():
    import importlib
    from scipy.stats import rankdata
    module=importlib.import_module("dynamicTreeCut.dynamicTreeCut")
    original=module.rankdata
    data=np.array([2.,2.,1.,3.,2.])
    with stages.persist_network_trees():
        ordinal=module.rankdata(data,method="ordinal")
        np.testing.assert_array_equal(ordinal,rankdata(data,method="ordinal"))
        assert ordinal.dtype.kind in {"i","u"}
        np.testing.assert_array_equal(module.rankdata(data,method="average"),rankdata(data,method="average"))
    assert module.rankdata is original
