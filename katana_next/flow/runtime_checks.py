"""Non-computational probes in the exact runtime chosen for one flow task."""
from __future__ import annotations

import json
from pathlib import Path
import sys

from .dependency_checks import python_dependency_issues, require_python_imports
from .model import load_plan, canonical_hash


def check_task_runtime(plan_path, task_id, *, project_path=None, runtime_label=None):
    """Import software and query R namespaces; never open scientific matrices.

    This executes inside the task container immediately before its calculation.
    Earlier project preflight and later input/output acceptance remain separate.
    A declared external FIMO interpreter is tested rather than replaced by the
    current Python. No package is installed and no analysis is started.
    """
    from .stages import _executable, _require_r_packages
    plan = load_plan(plan_path)
    task = next((row for row in plan["tasks"] if row["id"] == task_id), None)
    if task is None:
        raise ValueError("Runtime probe task is absent from the sealed plan")
    task_identity = canonical_hash(task)
    adapter = task.get("adapter", "")
    parameters = task.get("parameters", {})
    errors, executables = [], {"python": sys.executable}
    project = {}
    if project_path is not None:
        project = json.loads(Path(project_path).read_text(encoding="utf-8-sig"))
    from .container_runtime import sanitize_container_environment, validate_task_runtime_boundaries, project_asset_path
    removed_environment = sanitize_container_environment()
    try:
        validate_task_runtime_boundaries(task, project=project, project_path=project_path)
    except (ValueError, OSError, KeyError) as exc:
        return {"schema": "katana.flow.runtime_probe/v1", "status": "failed",
                "plan_id": plan["plan_id"], "task_id": task_id,
                "task_contract_sha256": task_identity, "workflow_id": task["workflow_id"],
                "execution_group": task.get("execution_group"), "runtime_label": runtime_label,
                "executables": executables, "errors": [str(exc)],
                "removed_environment_overrides": removed_environment, "scientific_computation": False}
    native_source = None
    if task.get("inputs", {}).get("native_source_provider"):
        try:
            from .bundled_runtime import resolve_bundled_stage
            task = resolve_bundled_stage(task)
            native_source = {key: task["inputs"][key] for key in
                             ("native_source_provider", "native_commit", "counter_sha256", "source_sha256")}
        except (ValueError, OSError, KeyError) as exc:
            errors.append(str(exc))
    imports = ("katana_next.flow.runtime", "numpy", "pandas", "pyarrow", "scipy", "matplotlib")
    if not adapter:
        # Workflow execution loads the same registry; probe its transitive
        # imports here rather than discovering a missing Bio/compiled package
        # only after the task has opened its output catalogue.
        imports += ("katana_next.workflows.registry",)
    try:
        require_python_imports(sys.executable, imports)
    except (ValueError, OSError) as exc:
        errors.append(str(exc))
    errors.extend(python_dependency_issues(task, resolve_executable=_executable))
    try:
        if adapter in {"network.fit", "network.bootstrap", "mapp.maei"}:
            configured = parameters.get("rscript")
            if adapter.startswith("network."):
                configuration = task.get("inputs", {}).get("config")
                if not configured and configuration:
                    from katana_cosplice.config import load_config
                    configured = load_config(configuration).network.wgcna_rscript
                configured = configured or project.get("metadata", {}).get("rscript_path")
            rscript = _executable(configured, "Rscript")
            packages = ("Rcpp", "optparse") if adapter == "mapp.maei" else ("WGCNA",)
            _require_r_packages(rscript, packages, parameters)
            executables["Rscript"] = rscript
        if adapter == "qc.genebody":
            executables["geneBody_coverage.py"] = _executable(parameters.get("rseqc_executable"), "geneBody_coverage.py")
            executables["Rscript"] = _executable(parameters.get("rscript"), "Rscript")
        statistical_r = {"differential.expression_response": ("DESeq2", "fgsea"),
                         "differential.replicate_interaction": ("DEXSeq", "BiocParallel")}
        if task["workflow_id"] in statistical_r:
            rscript = _executable(project.get("metadata", {}).get("rscript_path"), "Rscript")
            library_parameters = parameters
            if task["workflow_id"] == "differential.expression_response" and parameters.get("config_path"):
                library_parameters = json.loads(Path(parameters["config_path"]).read_text(encoding="utf-8-sig"))
            _require_r_packages(rscript, statistical_r[task["workflow_id"]], library_parameters)
            executables["Rscript"] = rscript
        if task["workflow_id"] in {"prioritize.event_motifs", "prioritize.rbp_mechanism", "prioritize.se_rbp_impact"}:
            resources = [asset for asset in project.get("assets", []) if asset.get("kind") == "fimo_runtime"]
            if not resources:
                raise ValueError("Motif enrichment runtime probe requires the declared fimo_runtime resource")
            for index, asset in enumerate(resources):
                path = project_asset_path(asset, project, project_path)
                config = json.loads(path.read_text(encoding="utf-8-sig"))
                if config.get("wsl_distribution"):
                    raise ValueError("Pipeline task runtime requires a native FIMO worker configuration; a Windows WSL launcher is not portable inside its container")
                python = _executable(config.get("python"), sys.executable)
                require_python_imports(python, ("pymemesuite.common", "pymemesuite.fimo"))
                executables[f"fimo_worker_{index}"] = python
    except (ValueError, OSError, ImportError, KeyError) as exc:
        errors.append(str(exc))
    return {"schema": "katana.flow.runtime_probe/v1", "status": "failed" if errors else "ready",
            "plan_id": plan["plan_id"], "task_id": task_id,
            "task_contract_sha256": task_identity, "workflow_id": task["workflow_id"],
            "execution_group": task.get("execution_group"), "runtime_label": runtime_label,
            "executables": executables, "native_source": native_source, "errors": errors,
            "removed_environment_overrides": removed_environment,
            "scientific_computation": False, "scope": "Runtime imports/namespaces only; no scientific acceptance claim"}
