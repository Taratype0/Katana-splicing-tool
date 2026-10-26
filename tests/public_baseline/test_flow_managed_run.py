"""Storage orchestration tests; native/scientific work is explicitly generated.

Real Nextflow failure/resume is exercised by check_managed_storage.py; actual
retained-project cards and Sashimi are covered by test_flow_storage_project.
"""
import json
from pathlib import Path
import shutil
import sys
import zipfile

import pytest

from katana_next.flow import managed
from katana_next.flow.model import file_sha256, read_json, write_json


@pytest.fixture
def lifecycle(tmp_path, monkeypatch):
    from katana_next.flow import collect, importer, storage_project
    inputs = tmp_path / "borrowed"
    inputs.mkdir()
    for name in ("samples.csv", "contrasts.csv", "genome.fa", "annotation.gtf", "template.json", "flow.json", "project.json"):
        (inputs / name).write_text("generated input\n")
    write_json(inputs / 'flow.json', {})
    calls, fail = [], {"downstream": False, "acceptance": False}
    monkeypatch.setattr(managed, "_software_identity", lambda _: "generated-code-identity")
    monkeypatch.setattr(collect, "validate_collection", lambda _: {"plan_id": "generated-plan"})

    def native(command, work, scope, stage, *, resume=False):
        calls.append((stage, resume))
        work.mkdir(parents=True, exist_ok=True)
        (work / "native-intermediate").write_text("owned temporary data")
        logs = scope / "logs"
        logs.mkdir(exist_ok=True)
        (logs / (stage + ".log")).write_text("Generated lifecycle executor, not scientific acceptance\n")
        output = Path(command[command.index("--outdir") + 1])
        if stage == "upstream":
            bundle = output / "katana_bundle"
            bundle.mkdir(parents=True)
            data = bundle / "native.tsv"
            data.write_text("id\tvalue\na\t1\n")
            write_json(bundle / "manifest.json", {
                "schema": "katana-pipeline-bundle/v1", "status": "complete",
                "samples": [{"sample_id": "a", "condition": "A"}],
                "comparisons": [{"comparison_id": "B_vs_A"}],
                "artifacts": [{"relative_path": data.name, "status": "available",
                               "bytes": data.stat().st_size, "sha256": file_sha256(data)}]})
            (output / "models").mkdir()
            (output / "models" / "fit.rds").write_text("native statistics retained")
        else:
            if fail["downstream"]:
                raise RuntimeError("Generated native downstream failure")
            report = output / "report"
            report.mkdir(parents=True)
            (output / "collection").mkdir()
            (output / "collection" / "data.tsv").write_text("unchanged statistics")
            (report / "REPORT.html").write_text('<img src="plot.svg">')
            (report / "plot.svg").write_text('<svg xmlns="http://www.w3.org/2000/svg"/>')
            write_json(report / "report-manifest.json", {"generated_fixture": True})
            with zipfile.ZipFile(report / "REPORT.zip", "w") as archive:
                archive.write(report / "REPORT.html", "REPORT.html")
                archive.write(report / "plot.svg", "plot.svg")
            write_json(report / "flow.acceptance.json", {
                "schema": "katana.flow.report_acceptance/v1", "plan_id": "generated-plan",
                "status": "failed" if fail["acceptance"] else "completed",
                "archive_sha256": file_sha256(report / "REPORT.zip")})
            (output / "pipeline_info").mkdir()
            (output / "pipeline_info" / "trace.tsv").write_text("native stage diagnostics")
        return {"engine": "nextflow", "state": "COMPLETED", "exit_code": 0, "engine_stopped": True}

    def imported(bundle, template, out, **kwargs):
        out.mkdir(parents=True)
        write_json(out / "project.json", {"fixture": "generated imported inputs"})
        write_json(out / "import.receipt.json", {"status": "completed"})
        return out / "project.json"

    def retained(collection, source_project, destination, **kwargs):
        destination.mkdir()
        shutil.copytree(collection, destination / "collection")
        write_json(destination / "project.json", {"input": "independent.txt"})
        (destination / "independent.txt").write_text("retained evidence")
        write_json(destination / "retention.receipt.json", {"generated_fixture": True})
        return {"status": "completed"}

    def validate(destination):
        assert (destination / "independent.txt").read_text() == "retained evidence"
        return {"status": "completed"}

    monkeypatch.setattr(managed, "_native_phase", native)
    monkeypatch.setattr(importer, "import_bundle", imported)
    monkeypatch.setattr(storage_project, "retain_collection_project", retained)
    monkeypatch.setattr(storage_project, "validate_retained_project", validate)
    command = ["--outdir", str(tmp_path / "results"), "--input", str(inputs / "samples.csv"),
               "--contrasts", str(inputs / "contrasts.csv"), "--fasta", str(inputs / "genome.fa"),
               "--gtf", str(inputs / "annotation.gtf"), "--project_template", str(inputs / "template.json"),
               "--flow_config", str(inputs / "flow.json")]
    return command, inputs, calls, fail


def test_two_handoffs_delete_only_registered_work_and_keep_full_delivery(lifecycle):
    command, inputs, calls, _ = lifecycle
    result = managed.run(managed.parser().parse_args(command))
    out = Path(result["report"]).parent.parent
    scope = out.with_name(out.name + ".run")
    assert result["status"] == "completed" and result["native_resume_available"] is False
    assert calls == [("upstream", False), ("downstream", False)]
    assert not (scope / "work" / "upstream").exists()
    assert not (scope / "work" / "downstream").exists()
    assert (out / "upstream-support" / "models" / "fit.rds").read_text() == "native statistics retained"
    assert (out / "pipeline_info" / "trace.tsv").is_file()
    assert (out / "report" / "plot.svg").is_file()
    assert (scope / "logs" / "upstream.log").is_file()
    assert all(path.read_text() == "generated input\n" for path in inputs.iterdir() if path.name != 'flow.json')
    assert read_json(inputs / 'flow.json') == {}
    assert managed.run(managed.parser().parse_args([*command, "--resume"]))["status"] == "completed"
    assert len(calls) == 2


def test_failed_downstream_keeps_cache_and_resume_does_not_repeat_upstream(lifecycle):
    command, _, calls, fail = lifecycle
    fail["downstream"] = True
    with pytest.raises(RuntimeError, match="downstream failure"):
        managed.run(managed.parser().parse_args(command))
    out = Path(command[1])
    scope = out.with_name(out.name + ".run")
    assert not (scope / "work" / "upstream").exists()
    assert (scope / "work" / "downstream" / "native-intermediate").is_file()
    fail["downstream"] = False
    result = managed.run(managed.parser().parse_args([*command, "--resume"]))
    assert result["status"] == "completed"
    assert calls == [("upstream", False), ("downstream", False), ("downstream", True)]


def test_report_acceptance_failure_never_authorizes_final_cleanup(lifecycle):
    command, _, _, fail = lifecycle
    fail["acceptance"] = True
    with pytest.raises(ValueError, match="did not pass"):
        managed.run(managed.parser().parse_args(command))
    out = Path(command[1])
    assert (out.with_name(out.name + ".run") / "work" / "downstream").exists()
    assert not out.exists()


def test_keep_work_preserves_both_native_sessions(lifecycle):
    command, _, _, _ = lifecycle
    result = managed.run(managed.parser().parse_args([*command, "--keep-work"]))
    scope = Path(command[1]).with_name("results.run")
    assert result["native_resume_available"] is True
    assert all((scope / "work" / stage).is_dir() for stage in ("upstream", "downstream"))
    preview = managed.clean_completed_run(scope)
    assert preview["dry_run"] is True and (scope / "work" / "upstream").exists()
    cleaned = managed.clean_completed_run(scope, apply=True)
    assert cleaned["status"] == "completed"
    assert all(not (scope / "work" / stage).exists() for stage in ("upstream", "downstream"))


def test_crash_after_final_cleanup_recovers_without_recreating_deleted_work(lifecycle):
    command, _, calls, _ = lifecycle
    result = managed.run(managed.parser().parse_args(command))
    scope = Path(command[1]).with_name("results.run")
    state = read_json(scope / "managed-run.json")
    state["status"] = "failed_or_interrupted"
    write_json(scope / "managed-run.json", state)
    assert managed.run(managed.parser().parse_args([*command, "--resume"]))["status"] == "completed"
    assert len(calls) == 2 and not (scope / "work" / "downstream").exists()


def test_changed_config_or_report_asset_is_not_silently_reused(lifecycle):
    command, inputs, calls, _ = lifecycle
    managed.run(managed.parser().parse_args(command))
    flow = inputs / "flow.json"
    original = flow.read_text()
    flow.write_text("changed scientific configuration")
    with pytest.raises(ValueError, match="differs"):
        managed.run(managed.parser().parse_args([*command, "--resume"]))
    flow.write_text(original)
    (Path(command[1]) / "report" / "plot.svg").write_text("changed")
    with pytest.raises(ValueError, match="delivery member"):
        managed.run(managed.parser().parse_args([*command, "--resume"]))
    assert len(calls) == 2


def test_native_failure_without_terminal_evidence_cannot_be_blindly_resubmitted(tmp_path):
    scope = tmp_path / "scope"
    scope.mkdir()
    work = scope / "work" / "downstream"
    code = 'import sys; print("generated failure"); sys.exit(7)'
    with pytest.raises(RuntimeError, match="failed"):
        managed._native_phase([sys.executable, "-c", code], work, scope, "downstream")
    with pytest.raises(RuntimeError, match="termination is unverified"):
        managed._native_phase([sys.executable, "-c", code], work, scope, "downstream", resume=True)
    assert len(list((scope / "logs").glob("*.native.json"))) == 1


def test_default_runtime_cache_is_outside_reclaimable_stage_and_existing_setting_is_kept(tmp_path, monkeypatch):
    scope = tmp_path / "scope"
    scope.mkdir()
    monkeypatch.setenv("NXF_CONDA_CACHEDIR", str(tmp_path / "shared-conda"))
    code = 'import json,os,pathlib; pathlib.Path("env.json").write_text(json.dumps({k:v for k,v in os.environ.items() if k.startswith("NXF_")}))'
    work = scope / "work" / "upstream"
    managed._native_phase([sys.executable, "-c", code], work, scope, "upstream")
    observed = read_json(work / "launch" / "env.json")
    assert observed["NXF_CONDA_CACHEDIR"] == str(tmp_path / "shared-conda")
    for engine in ("SINGULARITY", "APPTAINER"):
        assert not Path(observed["NXF_" + engine + "_CACHEDIR"]).is_relative_to(work)


def test_tool_caches_inside_work_block_deletion(tmp_path):
    (tmp_path / "tasks" / "conda").mkdir(parents=True)
    with pytest.raises(ValueError, match="Tool cache"):
        managed._check_tool_cache_locations(tmp_path)


def test_nonstandard_environment_directory_is_not_deleted_as_an_intermediate(tmp_path):
    (tmp_path / "tasks" / "custom-env" / "conda-meta").mkdir(parents=True)
    with pytest.raises(ValueError, match="tool environment"):
        managed._check_tool_cache_locations(tmp_path)


def test_prelaunch_cache_validation_failure_is_known_not_started_and_can_resume(tmp_path, monkeypatch):
    scope = tmp_path / "scope"
    scope.mkdir()
    work = scope / "work" / "upstream"
    monkeypatch.setenv("NXF_CONDA_CACHEDIR", str(work / "tools"))
    command = [sys.executable, "-c", 'print("fixture native success")']
    with pytest.raises(ValueError, match="outside the reclaimable"):
        managed._native_phase(command, work, scope, "upstream")
    proof = read_json(scope / "logs" / "upstream-001.native.json")
    assert proof["state"] == "NOT_STARTED" and proof["engine_stopped"] is True
    monkeypatch.setenv("NXF_CONDA_CACHEDIR", str(tmp_path / "shared-conda"))
    assert managed._native_phase(command, work, scope, "upstream", resume=True)["exit_code"] == 0


def test_public_command_help_and_run_routing(capsys):
    from katana_next.flow.cli import main
    with pytest.raises(SystemExit) as stopped:
        main(["run", "--help"])
    assert stopped.value.code == 0
    assert "--keep-work" in capsys.readouterr().out
    for command in ("retain-project", "verify-retained", "clean-run"):
        with pytest.raises(SystemExit) as stopped:
            main([command, "--help"])
        assert stopped.value.code == 0


def test_managed_downstream_entry_has_no_upstream_submission(lifecycle):
    command, inputs, calls, _ = lifecycle
    result = managed.run(managed.parser().parse_args([
        "--mode", "downstream-only", "--project", str(inputs / "project.json"),
        "--flow_config", str(inputs / "flow.json"), "--outdir", command[1]]))
    assert result["status"] == "completed" and calls == [("downstream", False)]


def test_owned_bundle_links_never_link_borrowed_inputs_or_modify_projected_tables(tmp_path):
    from test_flow_import import bundle
    from katana_next.flow.importer import import_bundle
    source, _manifest, template = bundle(tmp_path)
    original = file_sha256(source / "SE.MATS.JCEC.txt")
    project = import_bundle(source, template, tmp_path / "linked", owned_link_root=tmp_path)
    assert file_sha256(source / "SE.MATS.JCEC.txt") == original
    assert project.is_file()
    with pytest.raises(ValueError, match="owned stage"):
        import_bundle(source, template, tmp_path / "outside" / "project", owned_link_root=source)
