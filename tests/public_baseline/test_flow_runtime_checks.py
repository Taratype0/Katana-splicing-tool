"""Runtime probes verify the selected task interpreter without analysis."""
import json
from pathlib import Path
import sys

import pytest

from katana_next.flow.model import seal_plan, write_json, file_sha256
from katana_next.flow import runtime_checks, stages


def planned(tmp_path, adapter=None, workflow="differential.as", parameters=None, assets=None, inputs=None):
    project = write_json(tmp_path / "project.json", {"assets": assets or []})
    task = {"id": "selected", "workflow_id": "flow." + adapter if adapter else workflow,
            "execution_group": "network_fit" if adapter == "network.fit" else "as",
            "parameters": parameters or {}, "dependencies": [], "inputs": inputs or {}}
    if adapter:
        task["adapter"] = adapter
    plan = seal_plan({"project_manifest": str(project), "project_manifest_sha256": file_sha256(project), "tasks": [task]})
    return write_json(tmp_path / "plan.json", plan), project


def test_failed_actual_runtime_import_blocks_probe_without_science(tmp_path, monkeypatch):
    plan, project = planned(tmp_path)
    def broken(*args, **kwargs):
        raise ValueError("selected environment has broken compiled dependency")
    monkeypatch.setattr(runtime_checks, "require_python_imports", broken)
    monkeypatch.setattr(stages, "execute_stage", lambda *_a, **_k: pytest.fail("No analysis allowed"))
    result = runtime_checks.check_task_runtime(plan, "selected", project_path=project, runtime_label="local-built-image")
    assert result["status"] == "failed" and "broken compiled" in str(result["errors"])
    assert result["scientific_computation"] is False
    assert result["runtime_label"] == "local-built-image"
    assert set(p.name for p in tmp_path.iterdir()) == {"project.json", "plan.json"}


def test_network_bootstrap_probes_declared_r_namespace_in_task_runtime(tmp_path, monkeypatch):
    plan, project = planned(tmp_path, "network.bootstrap", parameters={"rscript": "chosen-R"})
    seen = []
    monkeypatch.setattr(runtime_checks, "require_python_imports", lambda *a: None)
    monkeypatch.setattr(runtime_checks, "python_dependency_issues", lambda *a, **k: [])
    monkeypatch.setattr(stages, "_executable", lambda value, default: str(value or default))
    monkeypatch.setattr(stages, "_require_r_packages", lambda executable, modules, params: seen.append((executable, modules)))
    result = runtime_checks.check_task_runtime(plan, "selected", project_path=project)
    assert result["status"] == "ready" and seen == [("chosen-R", ("WGCNA",))]
    assert result["executables"]["Rscript"] == "chosen-R"


def test_fimo_probe_uses_declared_worker_not_host_python(tmp_path, monkeypatch):
    worker = write_json(tmp_path / "worker.json", {"python": "fimo-specific-python"})
    plan, project = planned(tmp_path, workflow="prioritize.rbp_mechanism",
                            assets=[{"kind": "fimo_runtime", "path": worker.name}])
    seen = []
    monkeypatch.setattr(runtime_checks, "require_python_imports", lambda executable, modules: seen.append((executable, modules)))
    monkeypatch.setattr(stages, "_executable", lambda value, default: str(value or default))
    result = runtime_checks.check_task_runtime(plan, "selected", project_path=project)
    assert result["status"] == "ready"
    assert seen[-1] == ("fimo-specific-python", ("pymemesuite.common", "pymemesuite.fimo"))
    assert seen[0][0] == sys.executable
    assert result["scientific_computation"] is False


def test_runtime_cli_records_failure_and_exits_nonzero(tmp_path, monkeypatch, capsys):
    from katana_next.flow.cli import main
    plan, project = planned(tmp_path)
    monkeypatch.setattr(runtime_checks, "check_task_runtime", lambda *a, **k: {"status": "failed", "errors": ["missing WGCNA"]})
    output = tmp_path / "runtime.json"
    code = main(["runtime-check", "--plan", str(plan), "--project", str(project), "--task-id", "selected", "--out", str(output)])
    assert code == 2 and json.loads(output.read_text())["status"] == "failed"
    assert json.loads(capsys.readouterr().out)["errors"] == ["missing WGCNA"]
    assert main(["runtime-check", "--plan", str(plan), "--project", str(project), "--task-id", "selected", "--out", str(output)]) == 1


def test_unknown_runtime_task_cannot_be_probed(tmp_path):
    plan, _ = planned(tmp_path)
    with pytest.raises(ValueError, match="absent from the sealed plan"):
        runtime_checks.check_task_runtime(plan, "not-declared")


@pytest.mark.parametrize("workflow, packages", [("differential.expression_response", ("DESeq2", "fgsea")),
                                                ("differential.replicate_interaction", ("DEXSeq", "BiocParallel"))])
def test_statistical_registry_workflows_probe_actual_R_namespace(tmp_path, monkeypatch, workflow, packages):
    plan, project = planned(tmp_path, workflow=workflow)
    seen = []
    monkeypatch.setattr(runtime_checks, "require_python_imports", lambda *a: None)
    monkeypatch.setattr(stages, "_executable", lambda value, default: str(value or default))
    monkeypatch.setattr(stages, "_require_r_packages", lambda executable, modules, params: seen.append((executable, modules)))
    result = runtime_checks.check_task_runtime(plan, "selected", project_path=project)
    assert result["status"] == "ready" and seen == [("Rscript", packages)]


def test_bundled_mapp_identity_is_checked_in_the_worker_without_fitting(tmp_path, monkeypatch):
    from katana_next.flow.bundled_runtime import mapp_input_defaults, MAPP_COMMIT
    inputs = mapp_input_defaults({})
    plan, project = planned(tmp_path, "mapp.prepare", inputs=inputs)
    before = file_sha256(plan)
    monkeypatch.setattr(runtime_checks, "require_python_imports", lambda *a: None)
    monkeypatch.setattr(runtime_checks, "python_dependency_issues", lambda *a, **k: [])
    result = runtime_checks.check_task_runtime(plan, "selected", project_path=project)
    assert result["status"] == "ready" and result["native_source"]["native_commit"] == MAPP_COMMIT
    assert "native_root" not in result["native_source"]
    assert file_sha256(plan) == before
    from katana_next.flow import bundled_runtime
    monkeypatch.setattr(bundled_runtime, "resolve_bundled_stage", lambda *a: (_ for _ in ()).throw(ValueError("worker source changed")))
    assert "worker source changed" in runtime_checks.check_task_runtime(plan, "selected", project_path=project)["errors"]


def test_raw_template_preflight_binds_real_bundled_mapp_sources(tmp_path, monkeypatch):
    from katana_next.flow.cli import template_preflight
    import shutil
    resource = tmp_path / "resource.txt"; resource.write_text("fixture")
    assets = [{"asset_id": kind, "kind": kind, "path": resource.name} for kind in
              ("gene_sets_gmt", "rbp_gene_sets", "rbp_motif_database", "fimo_runtime", "rbp_binding_peaks")]
    template = write_json(tmp_path / "template.json", {"metadata": {"species": "mus_musculus"}, "assets": assets})
    config = write_json(tmp_path / "config.json", {})
    from katana_next.flow import dependency_checks
    monkeypatch.setattr(dependency_checks, "python_dependency_issues", lambda *a, **k: [])
    monkeypatch.setattr(stages, "_executable", lambda value, default: str(value or default))
    monkeypatch.setattr(stages, "_require_r_packages", lambda *a: None)
    monkeypatch.setattr(shutil, "which", lambda name: str(name))
    monkeypatch.setattr(stages, "execute_stage", lambda *a, **k: pytest.fail("Preflight must not run analysis"))
    result = template_preflight(template, config)
    assert result["ready"] is True and result["errors"] == []
    assert json.loads(config.read_text(encoding="utf-8")) == {}
