import sqlite3

import pandas as pd
import pytest

from katana_next.client import KatanaProject
from katana_next.domain.models import ProjectManifest, Sample, AnnotationReference, JobStatus
from katana_next.full_fixture import create_full_fixture
from katana_next.science.contrast_planner import plan_contrasts
from katana_next.science.attrition import cumulative_attrition


def factorial(tmp_path, *, paired=False, missing=False, confounded_batch=False):
    samples = []
    for genotype in ("WT", "MUT"):
        for dp in ("noDP", "DP"):
            if missing and (genotype,dp) == ("MUT","DP"):
                continue
            for replicate in range(3):
                identity = f"{genotype}-{dp}-{replicate}"
                samples.append(Sample(sample_id=identity, biological_sample_id=f"donor-{replicate}" if paired else identity,
                    condition=f"{genotype}:{dp}", batch=genotype if confounded_batch else f"b{replicate}",
                    attributes={"genotype": genotype, "DP": dp}))
    return ProjectManifest(name="Factorial fixture", input_root=tmp_path, output_root=tmp_path/"out",
                           annotation=AnnotationReference(genome_build="fixture"), samples=tuple(samples))


def test_factorial_simple_marginal_interaction_and_missing_cells(tmp_path):
    manifest = factorial(tmp_path)
    before = manifest.model_dump_json()
    plan = plan_contrasts(manifest, ["genotype","DP"], reference_levels={"genotype":"WT","DP":"noDP"})
    assert len(plan["contrasts"]) == 7
    assert any("genotype: MUT minus WT | DP=DP" == row["label"] for row in plan["contrasts"])
    assert all(c["estimable"] and c["status"] == "design_estimable" for c in plan["contrasts"])
    assert plan["residual_df"] == 8
    interaction = next(c for c in plan["contrasts"] if c["kind"] == "interaction")
    assert sorted(interaction["coefficients"].values()) == [-1,-1,1,1]
    assert manifest.model_dump_json() == before
    missing = plan_contrasts(factorial(tmp_path, missing=True), ["genotype","DP"])
    assert not next(c for c in missing["contrasts"] if c["kind"] == "interaction")["estimable"]
    assert any(c["estimable"] for c in missing["contrasts"])


def test_batch_confounding_paired_and_technical_units(tmp_path):
    confounded = plan_contrasts(factorial(tmp_path, confounded_batch=True), ["genotype","DP"], include_batch=True)
    assert any(not c["estimable"] for c in confounded["contrasts"] if c["label"].startswith("genotype:"))
    paired = factorial(tmp_path, paired=True)
    plan = plan_contrasts(paired, ["genotype","DP"])
    assert plan["paired_block"] and plan["residual_df"] == 6
    assert all(c["complete_paired_units"] == 3 for c in plan["contrasts"])
    duplicate = paired.samples[0].model_copy(update={"sample_id":"technical-duplicate"})
    plan2 = plan_contrasts(paired.model_copy(update={"samples":(*paired.samples,duplicate)}), ["genotype","DP"])
    assert plan2["residual_df"] == plan["residual_df"]
    assert plan2["technical_partitions_collapsed"] == 1
    with pytest.raises(ValueError, match="independent"):
        plan_contrasts(paired, ["genotype","DP"], pairing="unpaired")


def test_attrition_zero_stage_is_visible_and_reasons_conserve_rows():
    frame = pd.DataFrame({"comparison_id":["A","A","B"]})
    ledger, reason = cumulative_attrition(frame, [("input",[True]*3), ("reads",[True,False,False]), ("significance",[False]*3)])
    assert len(ledger) == 6
    assert ledger.loc[ledger.stage.eq("significance"),"events"].eq(0).all()
    assert ledger.excluded_at_stage.sum() == 3
    assert reason.tolist() == ["significance","reads","reads"]
    with pytest.raises(ValueError, match="reintroduces"):
        cumulative_attrition(frame, [("drop",[False]*3), ("bad",[True]*3)])


def test_job_history_survives_restart_and_retry_is_isolated(tmp_path):
    manifest = create_full_fixture(tmp_path/"fixture")
    with KatanaProject(manifest) as project:
        workflow = project.workflows.registry.get("differential.deg")
        original = workflow.run
        def fail(*_):
            raise RuntimeError("Injected failure for recovery acceptance")
        workflow.run = fail
        with pytest.raises(RuntimeError, match="Injected failure"):
            project.run("differential.deg")
        job_id = project.jobs.list()[0].job_id
        workflow.run = original
    with KatanaProject(manifest) as reopened:
        assert len(reopened.jobs.list()) == 1
        assert reopened.jobs.get(job_id).status == JobStatus.FAILED
        result = reopened.retry(job_id)
        assert result.status == JobStatus.COMPLETE
        assert {j.analysis_id for j in reopened.jobs.list()} == {"differential.deg"}
        assert len(reopened.jobs.list()) == 2
        changed = reopened.workflows.registry.get("differential.deg")
        changed.version = "new-method-version"
        with pytest.raises(ValueError, match="changed"):
            reopened.retry(job_id)


def test_dead_executor_is_marked_interrupted_without_automatic_run(tmp_path):
    from katana_next.infrastructure.job_journal import JobJournal
    from katana_next.domain.models import JobRecord
    journal = JobJournal(tmp_path)
    job = JobRecord(project_id="test", analysis_id="differential.deg", status=JobStatus.RUNNING, run_fingerprint="f")
    journal.save(job)
    with sqlite3.connect(journal.path) as db:
        db.execute("UPDATE jobs SET pid=2147483647")
    recovered = journal.read()[0]
    assert recovered.status == JobStatus.FAILED
    assert "Interrupted" in recovered.message
    other = JobJournal(tmp_path)
    try:
        journal.acquire_executor()
        with pytest.raises(ValueError, match="Another Katana"):
            other.acquire_executor()
    finally:
        journal.close(); other.close()


def test_reader_refreshes_foreign_progress_without_claiming_control(tmp_path):
    from katana_next.infrastructure.job_journal import JobJournal
    from katana_next.infrastructure.jobs import JobManager
    from katana_next.domain.models import JobRecord
    writer = JobJournal(tmp_path)
    record = JobRecord(project_id="reader-fixture", analysis_id="differential.deg",
                       status=JobStatus.RUNNING, run_fingerprint="f")
    writer.save(record)
    reader = JobManager()
    try:
        reader.bind_project(record.project_id, tmp_path)
        with pytest.raises(ValueError, match="another executor"):
            reader.cancel(record.job_id)
        writer.save(record.model_copy(update={"status":JobStatus.COMPLETE, "progress":1.0}))
        assert reader.get(record.job_id).status == JobStatus.COMPLETE
        assert reader.get(record.job_id).progress == 1.0
    finally:
        reader.shutdown()
        writer.close()


def test_lineage_retains_different_units_and_missing_ancestors(tmp_path):
    from katana_next.services.lineage import artifact_lineage
    from katana_next.workflows.base import ArtifactWriter
    with KatanaProject(create_full_fixture(tmp_path/"lineage")) as project:
        writer = ArtifactWriter(project.context,"acceptance.lineage","f",{},"test")
        parent = writer.table("Events","events.parquet",pd.DataFrame({"row_id":["e1","e2"]}),entity_type="splice_event")
        child = writer.table("Genes","genes.parquet",pd.DataFrame({"row_id":["g1"]}),entity_type="gene",parents=(parent.artifact_id,"missing-parent"))
        writer.finish()
        result = artifact_lineage(project.context,child.artifact_id)
        assert len(result["nodes"]) == 3
        assert {n["rows"] for n in result["nodes"]} == {None,1,2}
        assert len(result["edges"]) == 2
        assert project.jobs.list() == []


def test_module_member_gene_view_accepts_canonical_ids_without_symbols(tmp_path):
    from katana_next.workflows.base import ArtifactWriter
    from katana_next.workflows.registry import WorkflowRegistry
    from katana_next.domain.models import SelectionSet
    with KatanaProject(create_full_fixture(tmp_path/"module-path")) as project:
        project.run("landscape.bulk")
        project.run("differential.as")
        source = next(a for a in project.artifacts("differential.as") if a.name == "Differential AS events")
        events = pd.read_parquet(source.path)
        gene = events.gene_id.iloc[0]
        writer = ArtifactWriter(project.context,"programs.diagnostics","module-fixture",{},WorkflowRegistry().get("programs.diagnostics").version)
        members = writer.table("Annotated module members","members.parquet",pd.DataFrame({"row_id":["m1"],"gene_id":[gene],"module":["M1"]}),entity_type="module_member")
        writer.finish()
        selection = SelectionSet(project_id=project.project_id,name="M1 genes",entity_type="module_member",source_artifact_id=members.artifact_id,entity_ids=("m1",))
        project.projects.save_selection(selection)
        before = {j.job_id for j in project.jobs.list()}
        project.run("prioritize.selected_gene_views",{"selection_id":selection.selection_id})
        assert {j.analysis_id for j in project.jobs.list() if j.job_id not in before} == {"prioritize.selected_gene_views"}
        audit = next(a for a in project.artifacts("prioritize.selected_gene_views") if a.name == "Selected gene coverage audit")
        assert pd.read_parquet(audit.path).gene_id.tolist() == [gene]
