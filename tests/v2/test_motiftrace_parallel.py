import json
from types import SimpleNamespace

import pandas as pd
import pytest

from katana_next.science.motiftrace_parallel import (
    SCHEMA, MAPP_P, MT_P, ParallelModelBundle, digest, records, render_parallel,
)
from katana_next.workflows.motiftrace_parallel import MotifTraceParallelWorkflow
from katana_next.workflows.rbp_mechanism import RbpMechanismWorkflow


@pytest.fixture
def bundle(tmp_path):
    plan = dict(schema=SCHEMA, state="complete", species="human", genome_build="hg19", new_inference=False,
        family=dict(models=1, units_per_case=2, windows=2, MotifTrace_branches=4),
        methods={"MotifTrace": {}, "MAPP": {}}, limitations=[], files={},
        cases={"KD": dict(numerator=["k1", "k2"], denominator=["c1", "c2"], event_count=1, background_count=1)})
    def put(role, frame):
        name = role.replace(":", "_") + ".tsv"
        frame.to_csv(tmp_path / name, sep="\t", index=False)
        plan["files"][role] = dict(path=name, bytes=(tmp_path / name).stat().st_size, sha256=digest(tmp_path / name))
    put("models", pd.DataFrame([dict(model_id="m", literal_pattern="AGGAG")]))
    windows = pd.DataFrame([dict(model_id="m", pattern="AGGAG", case="KD", site="5ss", start=a, end=a+50,
                                 MAPP_fitted=f, **{MT_P: .06, MAPP_P: .001 if f else None}) for a,f in [(-50,True),(0,False)]])
    put("KD:windows", windows)
    put("KD:branches", pd.DataFrame([dict(model_id="m", case="KD", site="5ss", start=a, end=a+50, channel=c, score=s)
        for a in (-50,0) for c in ("PSI","fraction") for s in ("count","coverage")]))
    put("KD:usage", pd.DataFrame([dict(event_id="e", role="background")]))
    put("KD:design", pd.DataFrame([dict(sample="k1")]))
    put("window_sequences", pd.DataFrame([dict(case="KD", event_id="e", pattern="AGGAG", chr="chr1", strand="-",
                window_start0=100, window_end0=105, sequence_DNA="AGGAG", sequence_RNA="AGGAG")]))
    put("window_hits", pd.DataFrame([dict(case="KD", event_id="e", pattern_DNA="AGGAG", chr="chr1", strand="-",
              RNA_start0=0, RNA_end0=5, genomic_start0=100, genomic_end0=105, matched_RNA="AGGAG")]))
    put("window_availability", pd.DataFrame([dict(case="KD", event_id="e", pattern="AGGAG", PSI_KD_source_slot1=.2)]))
    (tmp_path / "calibration.json").write_text("{}")
    plan["files"]["KD:calibration"] = dict(path="calibration.json", bytes=2, sha256=digest(tmp_path / "calibration.json"))
    def seal():
        path = tmp_path / "manifest.json"
        path.write_text(json.dumps(plan))
        return SimpleNamespace(path=path, metadata=dict(sha256=digest(path)))
    return tmp_path, plan, put, seal


def load(seal):
    return ParallelModelBundle(seal(), species="human", genome_build="hg19")


def test_separate_P_and_NA_no_combined_discovery(bundle):
    _, _, _, seal = bundle
    p = load(seal).payload()["cases"]["KD"]
    assert p["windows"][0][MT_P] == .06
    assert p["windows"][0][MAPP_P] == .001
    assert p["windows"][1][MAPP_P] is None
    assert p["availability"][0]["PSI_KD_source_slot1"] == .2


@pytest.mark.parametrize("change,match", [
    (lambda p: p.update(state="running"), "completed"),
    (lambda p: p.update(new_inference=True), "inference"),
    (lambda p: p["files"]["models"].update(path="../models.tsv"), "Unsafe"),
    (lambda p: p["files"]["models"].update(path="C:/models.tsv"), "Unsafe"),
    (lambda p: p["family"].update(models=2), "model identities"),
    (lambda p: p["family"].update(units_per_case=1), "family"),
    (lambda p: p["cases"]["KD"].update(background_count=0), "background"),
])
def test_incomplete_or_changed_contract_rejected(bundle, change, match):
    _, plan, _, seal = bundle
    change(plan)
    with pytest.raises(ValueError, match=match):
        load(seal).payload()


def test_hash_species_and_sample_identity_rejected(bundle):
    root, _, _, seal = bundle
    asset = seal()
    with pytest.raises(ValueError, match="species/genome"):
        ParallelModelBundle(asset, species="mouse", genome_build="hg19")
    reader = load(seal)
    contrast = SimpleNamespace(contrast_id="KD", numerator=SimpleNamespace(sample_ids=["c1","c2"]), denominator=SimpleNamespace(sample_ids=["k1","k2"]))
    with pytest.raises(ValueError, match="identity mismatch"):
        reader.bind(contrast)
    (root / "models.tsv").write_text("changed")
    with pytest.raises(ValueError, match="Changed"):
        ParallelModelBundle(asset, species="human", genome_build="hg19")


def test_na_not_zero_and_motif_identity_strict(bundle):
    _, _, put, seal = bundle
    w = load(seal).table("KD:windows")
    w.loc[1, MAPP_P] = 1
    put("KD:windows", w)
    with pytest.raises(ValueError, match="non-applicability"):
        load(seal).payload()
    w.loc[1, MAPP_P] = None
    w.loc[0, "pattern"] = "GGAG"
    put("KD:windows", w)
    with pytest.raises(ValueError, match="sequence/model"):
        load(seal).payload()


def test_wrong_genomic_coordinates_and_ambiguous_join_rejected(bundle):
    _, _, put, seal = bundle
    hits = load(seal).table("window_hits")
    hits.loc[0, "genomic_start0"] = 99
    put("window_hits", hits)
    with pytest.raises(ValueError, match="coordinates differ"):
        load(seal).payload()
    hits.loc[0, "genomic_start0"] = 100
    put("window_hits", hits)
    seq = load(seal).table("window_sequences")
    put("window_sequences", pd.concat([seq, seq]))
    with pytest.raises(ValueError, match="Ambiguous"):
        load(seal).payload()


def test_full_family_branch_deletion_not_allowed(bundle):
    _, _, put, seal = bundle
    b = load(seal).table("KD:branches").iloc[:-1]
    put("KD:branches", b)
    with pytest.raises(ValueError, match="branches changed"):
        load(seal).payload()


def test_tiny_original_P_and_script_text_remain_data(tmp_path):
    assert records(pd.DataFrame([dict(P=1.234567890123456e-42)]))[0]["P"] == 1.234567890123456e-42
    out = tmp_path / "report.html"
    render_parallel({"hostile": "</script><script>alert(1)</script>", "P": 1.234567890123456e-42}, out)
    assert "</script><script>alert(1)" not in out.read_text(encoding="utf-8")
    assert "1.234567890123456e-42" in out.read_text(encoding="utf-8")


def test_product_registers_parallel_without_changing_default():
    from katana_next.workflows.registry import WorkflowRegistry
    assert isinstance(WorkflowRegistry().get("prioritize.motiftrace_parallel"), MotifTraceParallelWorkflow)
    workflow = RbpMechanismWorkflow()
    assert workflow.normalize_parameters({})["sequence_scope"] == "boundary_association"
    assert "completed_parallel_models" in next(p for p in workflow.definition.parameters if p.key == "sequence_scope").options
    assert workflow.parent_artifacts(None, {"sequence_scope": "completed_parallel_models"}) == []
