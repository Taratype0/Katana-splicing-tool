"""Measured native MAEI inapplicability is distinct from missing computation."""
import json
from pathlib import Path
import sys

import pandas as pd
import pytest

from katana_next.flow import stages
from katana_next.flow.mapp_inputs import COUNTER


def fixture(tmp_path, *, mapped=False):
    gtf = tmp_path / "reference.gtf"
    gtf.write_text("".join(
        f'chr1\tfixture\texon\t{a+1}\t{b}\t.\t+\t.\tgene_id "G1"; transcript_id "{tx}";\n'
        for tx, blocks in (("INC.1", [(100,120),(220,340),(440,460)]),
                           ("EXC.1", [(100,120),(440,460)])) for a,b in blocks))
    events = tmp_path / "measured.tsv"
    pd.DataFrame([{"event_id":"SE:1", "event_type":"SE", "gene_id":"G1", "chr":"chr1",
                   "strand":"+", "upstreamEE":120, "exonStart_0base":220, "exonEnd":340,
                   "downstreamES":440 if mapped else 999, "comparison_id":"B_vs_A"}]).to_csv(events,sep="\t",index=False)
    fasta = tmp_path / "reference.fa"
    fasta.write_text(">chr1\n" + "ACGT"*250 + "\n",encoding="ascii",newline="\n")
    Path(str(fasta)+".fai").write_text("chr1\t1000\t6\t1000\t1001\n")
    quantities = {}
    for sample in ("A1","A2","B1","B2"):
        quant = tmp_path / (sample+".sf")
        pd.DataFrame({"Name":["INC.1","EXC.1"],"TPM":[3.,7.]}).to_csv(quant,sep="\t",index=False)
        quantities[sample] = str(quant)
    native = tmp_path / "native"
    counter = native / COUNTER
    counter.parent.mkdir(parents=True)
    counter.write_text("# Inapplicability fixture: no native scientific command may run.\n")
    return {"adapter":"mapp.prepare", "inputs":{"events":str(events), "gtf":str(gtf),
        "genome_fasta":str(fasta), "quantification":quantities, "native_root":str(native),
        "counter_sha256":stages._sha(counter)}, "parameters":{"python":sys.executable,
        "sample_ids":list(quantities), "comparison_id":"B_vs_A", "site":"3ss",
        "window_start":-50, "window_end":0}}


def forbidden(*_args, **_kwargs):
    pytest.fail("A verified zero-eligible universe must not manufacture a native fit")


def test_complete_zero_eligible_universe_propagates_as_proved_not_applicable(tmp_path, monkeypatch):
    monkeypatch.setattr(stages,"_run",forbidden)
    spec = fixture(tmp_path)
    prepared = stages.execute_stage(spec,tmp_path / "prepare")
    assert prepared["status"] == "not_applicable"
    proof = json.loads(Path(prepared["evidence"]["applicability_proof"]).read_text())
    assert proof["measured_events"] == 1 and proof["eligible_events"] == 0
    assert proof["complete_measured_universe"] is True and len(proof["source_inputs"]) == 9
    audit = pd.read_csv(proof["eligibility_table"],sep="\t")
    assert audit.status.tolist() == ["no_exact_annotated_included_or_skipped_transcript"]
    scope = {"contrast_ids":["B_vs_A"],"site":"3ss","window_start":-50,"window_end":0}
    fit = stages.execute_stage({"adapter":"mapp.maei", "scope":scope,
                                "inputs":prepared["evidence"]["maei_inputs"]},tmp_path / "fit")
    assert fit["status"] == "not_applicable" and fit["evidence"]["native_fitting_performed"] is False
    assert {Path(x["path"]).name for x in fit["outputs"]} == {"all_event_eligibility.tsv", "native_maei_applicability.json"}
    descendant = json.loads(Path(fit["evidence"]["applicability_proof"]).read_text())
    assert descendant["source_applicability_proof"]["sha256"] == stages._sha(Path(prepared["evidence"]["applicability_proof"]))
    assert Path(descendant["eligibility_table"]).parent == tmp_path / "fit"


def test_complete_non_se_universe_is_explicitly_not_applicable(tmp_path, monkeypatch):
    monkeypatch.setattr(stages,"_run",forbidden)
    spec = fixture(tmp_path)
    table = pd.read_csv(spec["inputs"]["events"],sep="\t")
    table["event_type"] = "RI"
    table.to_csv(spec["inputs"]["events"],sep="\t",index=False)
    result = stages.execute_stage(spec,tmp_path / "prepare")
    proof = json.loads(Path(result["evidence"]["applicability_proof"]).read_text())
    assert proof["source_event_rows"] == 1 and proof["comparison_event_rows"] == 1
    assert proof["measured_events"] == 0 and proof["loss_counts"] == {}
    assert len(pd.read_csv(proof["eligibility_table"],sep="\t")) == 0


def test_actual_private_task_receipts_preserve_native_fit_inapplicability(tmp_path, monkeypatch):
    from katana_next.domain.models import AnnotationReference, Contrast, ContrastArm, ProjectManifest, Sample
    from katana_next.flow.model import file_sha256, seal_plan, write_json
    from katana_next.flow.runtime import execute_task, validate_receipt
    from katana_next.infrastructure.manifest_store import ManifestStore
    monkeypatch.setattr(stages,"_run",forbidden)
    study = tmp_path / "study"
    study.mkdir()
    spec = fixture(study)
    manifest = ProjectManifest(name="Measured MAEI applicability fixture", input_root=study,
        output_root=study / "previous-results", annotation=AnnotationReference(genome_build="fixture",
            gtf=Path(spec["inputs"]["gtf"]), genome_fasta=Path(spec["inputs"]["genome_fasta"])),
        samples=tuple(Sample(sample_id=sid,biological_sample_id=sid,condition=sid[0],assets={"quant_sf":Path(path)})
                      for sid,path in spec["inputs"]["quantification"].items()),
        contrasts=(Contrast(contrast_id="B_vs_A",label="B versus A",
            numerator=ContrastArm(label="B",sample_ids=("B1","B2")),
            denominator=ContrastArm(label="A",sample_ids=("A1","A2"))),))
    source = ManifestStore().save(manifest,study / "project.json")
    before = source.read_bytes()
    scope = {"contrast_ids":["B_vs_A"],"site":"3ss","window_start":-50,"window_end":0}
    tasks = [{**spec,"id":"prepare", "workflow_id":"flow.mapp.prepare", "dependencies":[],"scope":scope},
             {"id":"fit","workflow_id":"flow.mapp.maei","adapter":"mapp.maei", "parameters":{},
              "dependencies":["prepare"], "scope":scope,"inputs":{}}]
    plan = write_json(tmp_path / "plan.json",seal_plan({"project_manifest":str(source),
        "project_manifest_sha256":file_sha256(source),"tasks":tasks}))
    prepared = execute_task(plan,"prepare",tmp_path / "prepare-task")
    fitted = execute_task(plan,"fit",tmp_path / "fit-task",[prepared])
    for receipt in (validate_receipt(prepared),validate_receipt(fitted)):
        assert receipt["status"] == "not_applicable"
        assert receipt["execution_evidence"]["native_fitting_performed"] is False
        assert receipt["generated_artifacts"] and receipt["outputs"]
        assert not any("AUTHOR_NATIVE_MIXTURE" in row["path"] for row in receipt["outputs"])
    assert source.read_bytes() == before


@pytest.mark.parametrize("mutation",["missing_schema","missing_quant","missing_comparison","missing_annotation","empty_fasta_index","missing_quantified_tx"])
def test_invalid_or_incomplete_inputs_cannot_be_reclassified_as_inapplicable(tmp_path, monkeypatch, mutation):
    monkeypatch.setattr(stages,"_run",forbidden)
    spec = fixture(tmp_path,mapped=mutation == "missing_quantified_tx")
    events = Path(spec["inputs"]["events"])
    if mutation == "missing_schema":
        pd.read_csv(events,sep="\t").drop(columns="event_id").to_csv(events,sep="\t",index=False)
    elif mutation == "missing_quant":
        Path(spec["inputs"]["quantification"]["A1"]).unlink()
    elif mutation == "missing_comparison":
        spec["parameters"]["comparison_id"] = "unmeasured_vs_A"
    elif mutation == "missing_annotation":
        Path(spec["inputs"]["gtf"]).write_text("# no declared exon annotation\n")
    elif mutation == "empty_fasta_index":
        Path(spec["inputs"]["genome_fasta"]+".fai").write_text("")
    else:
        pd.DataFrame({"Name":["EXC.1"],"TPM":[10.]}).to_csv(spec["inputs"]["quantification"]["A1"],sep="\t",index=False)
    with pytest.raises(ValueError):
        stages.execute_stage(spec,tmp_path / "bad")
    assert not (tmp_path / "bad/stage.receipt.json").exists()


@pytest.mark.parametrize("mutation",["audit_count","eligible","scope","source"])
def test_descendant_rejects_changed_or_mismatched_applicability_evidence(tmp_path, monkeypatch, mutation):
    monkeypatch.setattr(stages,"_run",forbidden)
    spec = fixture(tmp_path)
    result = stages.execute_stage(spec,tmp_path / "prepare")
    proof_path = Path(result["evidence"]["applicability_proof"])
    proof = json.loads(proof_path.read_text())
    scope = {"contrast_ids":["B_vs_A"],"site":"3ss","window_start":-50,"window_end":0}
    if mutation == "audit_count":
        proof["measured_events"] = 2
    elif mutation == "eligible":
        table = pd.read_csv(proof["eligibility_table"],sep="\t")
        table["status"] = "eligible"
        table.to_csv(proof["eligibility_table"],sep="\t",index=False)
        proof["eligibility_table_sha256"] = stages._sha(Path(proof["eligibility_table"]))
    elif mutation == "scope":
        scope["window_end"] = 25
    else:
        Path(spec["inputs"]["quantification"]["A1"]).write_text("changed\n")
    proof_path.write_text(json.dumps(proof))
    with pytest.raises(ValueError):
        stages.execute_stage({"adapter":"mapp.maei","scope":scope,
                              "inputs":{"applicability_proof":str(proof_path)}},tmp_path / "forbidden-fit")
    assert not (tmp_path / "forbidden-fit").exists()
