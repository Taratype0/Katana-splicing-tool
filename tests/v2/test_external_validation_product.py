from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from katana_next.client import KatanaProject
from katana_next.full_fixture import create_full_fixture
from katana_next.workflows.base import ArtifactWriter
from katana_next.workflows.external_validation import align_external, frozen_membership
from katana_next.workflows.registry import WorkflowRegistry


def matrices():
    rng = np.random.default_rng(9)
    reference = pd.DataFrame(rng.uniform(.1,.9,(12,8)), index=[f"e{i}" for i in range(12)], columns=[f"r{i}" for i in range(8)])
    target = pd.DataFrame(rng.uniform(.1,.9,(12,8)), index=reference.index, columns=[f"t{i}" for i in range(8)])
    samples = pd.DataFrame({"sample_id": target.columns, "biological_sample_id": target.columns})
    reference_samples = pd.DataFrame({"sample_id": reference.columns, "biological_sample_id": reference.columns})
    labels = pd.DataFrame({"event_key": reference.index,"consensus_module": ["M1"]*6+["M2"]*6})
    return reference,target,samples,reference_samples,labels


def test_external_alignment_excludes_unmeasurable_without_changing_membership():
    reference,target,samples,reference_samples,labels = matrices()
    target.iloc[0,:] = np.nan
    ref, external, audit = align_external(reference,target,samples,reference_samples,labels)
    assert len(ref) == len(external) == 11
    assert len(audit) == 12 and not audit.usable_in_both.iloc[0]
    assert set(audit.consensus_module) == {"M1","M2"}
    samples.loc[0,"biological_sample_id"] = reference_samples.biological_sample_id.iloc[0]
    with pytest.raises(ValueError,match="share biological"):
        align_external(reference,target,samples,reference_samples,labels)


@pytest.mark.parametrize("missing_event", [False, True])
def test_external_workflow_is_independent_and_resumes_cohort_checkpoints(tmp_path, monkeypatch, missing_event):
    import katana_next.workflows.external_validation as workflow
    manifest = create_full_fixture(tmp_path / "project")
    reference,target,samples,reference_samples,labels = matrices()
    if missing_event:
        target = target.drop(index="e0")
    target_path, sample_path = tmp_path/"external.tsv", tmp_path/"samples.tsv"
    target.to_csv(target_path,sep="\t")
    samples.to_csv(sample_path,sep="\t",index=False)
    cohort = SimpleNamespace(name="independent-cohort", psi_matrix=target_path, sample_metadata=sample_path)
    config = SimpleNamespace(external_cohorts=[cohort])
    monkeypatch.setattr(workflow,"external_config",lambda context: (SimpleNamespace(path=manifest),config,str(manifest)))
    calls = []
    def formal(ref, test, fixed, **kwargs):
        calls.append(kwargs)
        expected = labels[labels.event_key.isin(target.index)].reset_index(drop=True)
        pd.testing.assert_frame_equal(fixed,expected)
        assert ref.shape == test.shape == (len(target),8)
        return SimpleNamespace(wgcna=pd.DataFrame({"module":["M1","M2"],"Zsummary.pres":[3.,1.]}),
            netrep=pd.DataFrame({"module":["M1","M2"],"pvalue_avg.weight":[.01,.2],"pvalue_coherence":[.02,.4]}),
            audit=pd.DataFrame({"cohort":[cohort.name]}))
    monkeypatch.setattr(workflow,"run_formal_preservation",formal)
    with KatanaProject(manifest) as project:
        writer = ArtifactWriter(project.context,"programs.cosplice","frozen-discovery",{},WorkflowRegistry().get("programs.cosplice").version)
        matrix = writer.table("2. Build PSI matrix: psi","matrix__psi.parquet",reference.rename_axis("event_key").reset_index(),entity_type="matrix")
        writer.table("Samples","matrix__samples.parquet",reference_samples,entity_type="sample")
        writer.table("Events","matrix__events.parquet",pd.DataFrame({"event_key": reference.index}),entity_type="splice_event")
        writer.table("Assignments","network_observed__assignments.parquet",labels.rename(columns={"consensus_module":"module"}),entity_type="cosplice_module")
        writer.finish()
        parameters = {"matrix_artifact_id":matrix.artifact_id,"permutations":10,"soft_power":6}
        assert project.plan("programs.external_validation",parameters)["ready"]
        assert project.run("programs.external_validation",parameters).status == "complete"
        assert len(calls) == 1
        assert project.run("programs.external_validation",parameters).status == "cache_hit"
        # Explicit force re-registers outputs but keeps identical completed cohort computations.
        assert project.run("programs.external_validation",parameters,force=True).status == "complete"
        assert len(calls) == 1
        netrep = next(a for a in project.artifacts("programs.external_validation") if a.name == "External preservation netrep")
        results = pd.read_parquet(netrep.path)
        assert "batch_FDR_coherence" in results
        plot = next(a for a in project.artifacts("programs.external_validation") if a.name == "Official WGCNA module preservation")
        svg = open(plot.metadata["exports"]["svg"], encoding="utf-8").read()
        assert "M1" in svg and "M2" in svg and "NaN" not in svg
        assert {j.analysis_id for j in project.jobs.list()} == {"programs.external_validation"}


@pytest.mark.parametrize("side", ["reference", "external"])
def test_external_rejects_missing_or_duplicate_biological_identity(side):
    reference,target,samples,reference_samples,labels = matrices()
    metadata = reference_samples if side == "reference" else samples
    metadata.loc[0,"biological_sample_id"] = " "
    with pytest.raises(ValueError,match="declared biological"):
        align_external(reference,target,samples,reference_samples,labels)
    metadata.loc[0,"biological_sample_id"] = metadata.biological_sample_id.iloc[1]
    with pytest.raises(ValueError,match="aggregate technical"):
        align_external(reference,target,samples,reference_samples,labels)


def test_external_preserves_reasons_and_rejects_reference_out_of_range():
    reference,target,samples,reference_samples,labels = matrices()
    target = target.drop(index="e0")
    target.loc["e1"] = np.nan
    reference.loc["e2"] = .5
    _,_,audit = align_external(reference,target,samples,reference_samples,labels)
    assert len(audit) == 12
    reasons = audit.set_index("event_key").exclusion_reason
    assert reasons.e0 == "absent_external_event"
    assert "external_low_observation_count" in reasons.e1
    assert reasons.e2 == "reference_no_variation"
    reference.iloc[0,0] = 1.1
    with pytest.raises(ValueError,match="Reference PSI"):
        align_external(reference,target,samples,reference_samples,labels)


def test_refined_member_only_export_never_removes_frozen_null_background():
    reference,_,_,_,labels = matrices()
    reference.loc['raw-not-in-network'] = .5
    refined = labels.iloc[:3]
    ref, final = frozen_membership(reference, labels, refined)
    assert len(ref) == len(final) == 12
    assert 'raw-not-in-network' not in ref.index
    assert final.consensus_module.eq('grey').sum() == 9
    assert final.consensus_module.eq('M1').sum() == 3
    invalid = pd.concat([refined, pd.DataFrame({'event_key':['raw-not-in-network'], 'consensus_module':['M1']})])
    with pytest.raises(ValueError,match='unique subset'):
        frozen_membership(reference, labels, invalid)
