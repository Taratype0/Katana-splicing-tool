"""Container boundaries reject host tools before imports, claims or calculation."""
import json
import os
from pathlib import Path
import sys

import pytest

from katana_next.flow import container_runtime as boundary, runtime_checks, stages
from katana_next.flow.model import seal_plan, write_json, file_sha256


@pytest.fixture
def image(tmp_path, monkeypatch):
    prefix = tmp_path / "image"
    (prefix / "bin").mkdir(parents=True)
    python = prefix / "bin/python.exe"
    python.write_bytes(b"image executable fixture")
    monkeypatch.setattr(sys, "prefix", str(prefix))
    monkeypatch.setattr(boundary, "_IMAGE_PREFIX", prefix)
    monkeypatch.setattr(sys, "executable", str(python))
    for key in ("PATH", "PYTHONNOUSERSITE", "PYTHONSAFEPATH"):
        monkeypatch.setenv(key, os.environ.get(key, ""))
    monkeypatch.setenv("KATANA_CONTAINER", "1")
    return prefix


def plan(tmp_path, parameters=None, metadata=None, assets=None, workflow="differential.as"):
    project = write_json(tmp_path / "project.json", {"metadata": metadata or {}, "assets": assets or []})
    task = {"id": "test", "workflow_id": workflow, "parameters": parameters or {}, "dependencies": []}
    value = seal_plan({"tasks": [task], "project_manifest": str(project),
                       "project_manifest_sha256": file_sha256(project)})
    return write_json(tmp_path / "plan.json", value), project, value, task


def test_container_clears_inherited_overrides_but_standalone_retains_them(image, monkeypatch):
    for key in ("PYTHONPATH", "PYTHONHOME", "R_LIBS_USER", "R_HOME", "LD_LIBRARY_PATH"):
        monkeypatch.setenv(key, "/mounted/host/environment")
    removed = boundary.sanitize_container_environment()
    assert {"PYTHONPATH", "PYTHONHOME", "R_LIBS_USER", "R_HOME", "LD_LIBRARY_PATH"} <= set(removed)
    assert all(key not in os.environ for key in removed)
    monkeypatch.setenv("KATANA_CONTAINER", "0")
    monkeypatch.setenv("R_LIBS_USER", "/explicit/standalone/library")
    assert boundary.sanitize_container_environment() == []
    assert os.environ["R_LIBS_USER"] == "/explicit/standalone/library"
    assert boundary.require_image_path("/explicit/standalone/Rscript", "Rscript") == "/explicit/standalone/Rscript"


@pytest.mark.parametrize("configuration", ["parameters", "metadata", "fimo", "expression-config"])
def test_mounted_runtime_rejected_before_imports(image, tmp_path, monkeypatch, configuration):
    outside = tmp_path / "old-conda" / "python.exe"
    outside.parent.mkdir(); outside.write_bytes(b"host fixture")
    monkeypatch.setattr(boundary.shutil, "which", lambda name: str(name))
    kwargs = {}
    if configuration == "parameters":
        kwargs["parameters"] = {"rscript": str(outside)}
    elif configuration == "metadata":
        kwargs["metadata"] = {"rscript_path": str(outside)}
    elif configuration == "fimo":
        worker = write_json(tmp_path / "fimo.json", {"python": str(outside)})
        kwargs.update(workflow="prioritize.event_motifs", assets=[{"kind": "fimo_runtime", "path": worker.name}])
    else:
        config = write_json(tmp_path / "expression.json", {"r_library": str(outside.parent), "gtf": "/data/genome.gtf"})
        kwargs.update(workflow="differential.expression_response", parameters={"config_path": str(config)})
    planned, project, _, _ = plan(tmp_path, **kwargs)
    monkeypatch.setattr(runtime_checks, "require_python_imports", lambda *_: pytest.fail("Host interpreter must not be probed"))
    result = runtime_checks.check_task_runtime(planned, "test", project_path=project)
    assert result["status"] == "failed" and "external host path" in str(result["errors"])
    assert not result["scientific_computation"]


def test_explicit_r_libraries_rejected_before_R_probe(image, tmp_path, monkeypatch):
    outside = tmp_path / "host-R-library"; outside.mkdir()
    monkeypatch.setattr(stages, "_r_packages", lambda *_: pytest.fail("Must not probe a host R library"))
    with pytest.raises(ValueError, match="external host path"):
        stages._require_r_packages(sys.executable, ("WGCNA",), {"r_libraries": [str(outside)]})
    (image / "library").mkdir()
    boundary.validate_runtime_settings({"r_library": str(image / "library"), "genome_fasta": "/mounted/data/genome.fa"})


def test_container_checks_fimo_below_imported_input_root(image, tmp_path, monkeypatch):
    monkeypatch.setattr(boundary.shutil, "which", lambda value: str(value))
    worker = write_json(tmp_path / "bundle/resources/fimo.json", {"python": sys.executable})
    project = {"input_root": "bundle", "assets": [{"kind": "fimo_runtime", "path": "resources/fimo.json"}]}
    task = {"workflow_id": "prioritize.rbp_mechanism"}
    boundary.validate_task_runtime_boundaries(task, project=project, project_path=tmp_path / "project.json")
    write_json(worker, {"python": str(tmp_path / "host/python")})
    with pytest.raises(ValueError, match="external host path"):
        boundary.validate_task_runtime_boundaries(task, project=project, project_path=tmp_path / "project.json")


def test_direct_task_rejects_host_runtime_before_creating_output(image, tmp_path, monkeypatch):
    from katana_next.flow import runtime
    planned, project, value, task = plan(tmp_path, parameters={"r_library": str(tmp_path / "host-library")})
    monkeypatch.setattr(runtime, "load_plan_task", lambda *args: (planned, value, task, project))
    output = tmp_path / "task-output"
    with pytest.raises(ValueError, match="external host path"):
        runtime.execute_task(planned, "test", output, project_path=project)
    assert not output.exists()


def test_image_symlink_cannot_escape_its_prefix(image, tmp_path):
    outside = tmp_path / "host.exe"; outside.write_bytes(b"external")
    link = image / "bin/alias.exe"
    try:
        link.symlink_to(outside)
    except OSError:
        pytest.skip("Host does not permit creating this symlink fixture")
    with pytest.raises(ValueError, match="external host path"):
        boundary.require_image_path(link, "Python")


def test_mounted_launcher_cannot_certify_its_own_prefix(image, tmp_path, monkeypatch):
    outside = tmp_path / "old-conda"
    monkeypatch.setattr(sys, "prefix", str(outside))
    monkeypatch.setattr(sys, "executable", str(outside / "bin/python"))
    with pytest.raises(ValueError, match="fixed image runtime"):
        boundary.require_image_path(sys.executable, "Python")


def test_singularity_clean_environment_scoped_to_katana_labels():
    root = Path(__file__).resolve().parents[2]
    config = (root / "pipeline/nextflow.config").read_text(encoding="utf-8")
    assert "withLabel: /katana_(control|downstream|report)/" in config
    assert "['singularity', 'apptainer'] ? '--cleanenv'" in config
    main = (root / "pipeline/main.nf").read_text(encoding="utf-8")
    assert "params.flow_executable in ['katana-flow', '/opt/conda/bin/katana-flow']" in main


@pytest.mark.parametrize("workflow_name", ["expression", "dexseq"])
@pytest.mark.parametrize("configured", [None, "Rscript"])
def test_statistical_workflow_uses_resolved_image_command_without_fitting(tmp_path, monkeypatch, workflow_name, configured):
    from types import SimpleNamespace
    import pandas as pd
    from katana_next.workflows import expression_response, dexseq_interaction
    module = expression_response if workflow_name == "expression" else dexseq_interaction
    workflow = module.ExpressionResponseWorkflow() if workflow_name == "expression" else module.DexseqInteractionWorkflow()
    context = SimpleNamespace(manifest=SimpleNamespace(metadata={"rscript_path": configured}))
    output = tmp_path / "output"
    library = tmp_path / "R-library"; library.mkdir()
    monkeypatch.setattr(module, "ArtifactWriter", lambda *a, **k: SimpleNamespace(directory=output))
    selected = []
    def which(name):
        selected.append(name)
        return sys.executable if name == "Rscript" else None
    monkeypatch.setattr(stages.shutil, "which", which)
    captured = {}
    class ReachedNativeCommand(Exception):
        pass
    def capture(command, **kwargs):
        captured.update(command=command, environment=kwargs["env"])
        raise ReachedNativeCommand
    monkeypatch.setattr(module.subprocess, "Popen", capture)
    if workflow_name == "expression":
        config = {key: key + ".tsv" for key in ("design", "tximport", "gene_names", "gmt")}
        config["r_library"] = str(library)
        monkeypatch.setattr(workflow, "_config", lambda *a: config)
        parameters = {}
    else:
        monkeypatch.setattr(workflow, "_inputs", lambda *a: (None, pd.DataFrame({"sample_id": ["s1"]})))
        parameters = {"reference_a": "A", "reference_b": "B", "min_reads": 10,
                      "min_replicates": 2, "r_library": str(library)}
    job = SimpleNamespace(report=lambda *a: None)
    with pytest.raises(ReachedNativeCommand):
        workflow.run(context, job, "fixture", parameters)
    assert selected == ["Rscript"]
    assert captured["command"][0] == str(Path(sys.executable).absolute())
    assert captured["environment"]["R_LIBS_USER"] == str(library.resolve())
