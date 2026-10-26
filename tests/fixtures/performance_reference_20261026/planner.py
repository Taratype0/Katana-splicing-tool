"""Build a fixed-report calculation graph from an explicit study configuration."""
from __future__ import annotations

from pathlib import Path
from copy import deepcopy
import re
import shutil

from .model import file_sha256, read_json, seal_plan


def _token(value):
    stem = re.sub(r"[^A-Za-z0-9_.-]", "_", str(value))[:70]
    from .model import canonical_hash
    return stem + "-" + canonical_hash(str(value))[:10]


def _configuration(config):
    if isinstance(config, dict):
        return config, Path.cwd()
    path = Path(config).expanduser().resolve()
    return read_json(path), path.parent


def create_plan(project, config):
    from katana_next.infrastructure.manifest_store import ManifestStore
    from katana_next.workflows.registry import WorkflowRegistry
    project_path, manifest = ManifestStore().load(project)
    settings, base = _configuration(config)
    settings = deepcopy(settings)
    if settings.get("schema", "katana.flow.config/v1") != "katana.flow.config/v1":
        raise ValueError("Unsupported flow configuration schema")
    declared = {c.contrast_id for c in manifest.contrasts}
    contrasts = settings.get("contrast_ids", [c.contrast_id for c in manifest.contrasts])
    if not contrasts or len(contrasts) != len(set(contrasts)) or set(contrasts) - declared:
        raise ValueError("Select unique declared contrast IDs; comparisons are never generated exhaustively")
    event_set = settings.get("event_set", "stable_significant")
    if event_set not in {"stable_significant", "significant"}:
        raise ValueError("Choose RS-DAS stable_significant or ordinary DAS significant")
    from .selection import analysis_selection
    enabled = analysis_selection(settings)
    settings["analyses"] = enabled
    registry, tasks = WorkflowRegistry(), []

    def local(value):
        path = Path(value).expanduser()
        return str((path if path.is_absolute() else base / path).resolve())

    def resource_inputs(values):
        path_keys = {"native_root", "inclusion_table", "sitecount_matrix", "design", "gtf", "genome_fasta", "events", "config", "reference_bed"}
        return {key: local(value) if key in path_keys and value else value for key, value in values.items()}

    def add(identity, workflow, dependencies=(), parameters=None, scope=None, **extras):
        params = dict(parameters or {})
        if workflow.startswith("flow."):
            method_version = "1"
        else:
            runner = registry.get(workflow)
            params = runner.normalize_parameters({**settings.get("parameters", {}).get(workflow, {}), **params})
            method_version = runner.version
        groups = {"project.audit": "audit", "quality.stable_as": "quality", "landscape.bulk": "landscape",
                  "differential.as": "as", "differential.deg": "deg", "differential.dtu": "dtu",
                  "prioritize.evidence": "evidence", "differential.cross_comparison": "overlap",
                  "prioritize.enrichment": "enrichment", "prioritize.rbp_support": "enrichment",
                  "prioritize.motiftrace_consequences": "motiftrace", "prioritize.rbp_mechanism": "motiftrace",
                  "prioritize.motiftrace_parallel": "motiftrace",
                  "quality.read_coverage_bias": "coverage", "differential.replicate_interaction": "interaction"}
        execution_group = groups.get(workflow, "reference" if workflow == "flow.reference.bed12" else
                                     "mapp_prepare" if workflow == "flow.mapp.prepare" else "qc" if workflow in {"flow.qc.genebody", "flow.qc.sequencing"} else
                                     "network_fit" if workflow == "flow.network.fit" else
                                     "network_followup" if workflow.startswith("flow.network.") else "mapp")
        task = {"id": identity, "workflow_id": workflow, "workflow_version": method_version, "execution_group": execution_group,
                "parameters": params, "dependencies": list(dependencies), "scope": scope or {"contrast_ids": contrasts},
                **extras}
        tasks.append(task)
        return identity

    add("audit", "project.audit")
    add("quality", "quality.stable_as", ["audit"], {"contrast_ids": contrasts})
    if enabled["landscape"]:
        add("landscape", "landscape.bulk", ["quality"], {"contrast_ids": contrasts})
    for suffix in ("as", "deg", "dtu"):
        add(suffix, "differential." + suffix, ["quality"], {"contrast_ids": contrasts})
    add("evidence", "prioritize.evidence", ["as", "deg", "dtu"], {"contrast_ids": contrasts, "event_set": event_set})
    if enabled["cross_comparison"] and len(contrasts) >= 2:
        add("overlap", "differential.cross_comparison", ["as"], {"contrast_ids": contrasts, "event_set": event_set})
    parallel_settings = settings.get("parameters", {}).get("prioritize.motiftrace_parallel", {})
    parallel_requested = bool(parallel_settings.get("parallel_bundle_asset_id")) or any(
        asset.kind == "motiftrace_parallel_results" for asset in manifest.assets)
    for cid in contrasts:
        token, scope = _token(cid), {"contrast_ids": [cid]}
        selection = {"analysis_id": "differential.as", "name": "Differential AS events",
                     "event_set": "stable_significant", "contrast_ids": [cid], "scope_policy": "RS-DAS required for motif/CLIP"}
        for foreground in (("deg", "das", "dtu", "dexseq") if enabled["enrichment"] else ()):
            add("enrichment-" + foreground + "-" + token, "prioritize.enrichment", ["evidence"],
                {"contrast_ids": [cid], "foreground_source": foreground, "background_source": "match_foreground", "event_set": event_set}, scope)
        if enabled["rbp_identity"]:
            add("rbp-identity-" + token, "prioritize.rbp_support", ["evidence"],
                {"contrast_ids": [cid], "foreground_source": "deg", "background_source": "match_foreground", "event_set": event_set}, scope)
        if enabled["motiftrace"]:
            add("motiftrace-" + token, "prioritize.motiftrace_consequences", ["as"],
                {"contrast_ids": [cid]}, scope, selection=selection)
        if enabled["motif_binding"]:
            add("motif-boundaries-" + token, "prioritize.rbp_mechanism", ["as"],
                {"contrast_ids": [cid], "binding_mode": "required"}, scope, selection=selection)
        if enabled["mapp"] and parallel_requested:
            add("motif-parallel-" + token, "prioritize.motiftrace_parallel", ["as"],
                {"contrast_ids": [cid]}, scope)

    qc = settings.get("qc", {})
    if enabled["sequencing_qc"]:
        add("sequencing-qc", "flow.qc.sequencing", ["audit"], {}, {"unit": "available_native_sequencing_QC"},
            adapter="qc.sequencing", inputs={"project": str(project_path)})
    if enabled["coverage"]:
        if not qc.get("reference_bed"):
            add("reference-bed", "flow.reference.bed12", ["audit"], {}, {"unit": "reference"},
                adapter="reference.bed12", inputs={"gtf": str(manifest.annotation.gtf or "")})
        coverage_tasks = []
        for sample in manifest.samples:
            inputs = {"bam": str(sample.assets.get("bam", sample.assets.get("genome_bam", ""))),
                      "reference_bed": local(qc["reference_bed"]) if qc.get("reference_bed") else ""}
            dependency_inputs = {} if qc.get("reference_bed") else {
                "reference_bed": {"analysis_id": "flow.reference.bed12", "name": "reference_transcripts.bed", "kind": "file"}}
            coverage_tasks.append(add("genebody-" + _token(sample.sample_id), "flow.qc.genebody", ["audit"] if qc.get("reference_bed") else ["reference-bed"],
                {k: v for k, v in qc.items() if k != "reference_bed"}, {"sample_id": sample.sample_id},
                adapter="qc.genebody", inputs=inputs, dependency_inputs=dependency_inputs))
        if coverage_tasks:
            add("coverage-summary", "quality.read_coverage_bias", coverage_tasks)

    if enabled["network"]:
        network = settings.get("network", {})
        if not isinstance(network, dict):
            raise ValueError("network must be a configuration object")
        bootstrap = network.get("bootstrap", {})
        if not isinstance(bootstrap, dict) or set(bootstrap) - {"enabled"}:
            raise ValueError("network.bootstrap must contain only the explicit enabled boolean")
        bootstrap_enabled = bootstrap.get("enabled", False)
        if not isinstance(bootstrap_enabled, bool):
            raise ValueError("network.bootstrap.enabled must be true or false")
        network["bootstrap"] = {"enabled": bootstrap_enabled}
        if not isinstance(network.get("parameters", {}), dict):
            raise ValueError("network.parameters must be a configuration object")
        network["parameters"] = dict(network.get("parameters", {}))
        network["parameters"].setdefault("hub_metric", "abs_kME")
        settings["network"] = network
        from .network_inputs import network_scope_identity, resolve_network_rmats_mode
        network_scope = network.get("scope", "contrast")
        if network_scope not in {"contrast", "cohort"}:
            raise ValueError("network.scope must be contrast or cohort")
        configurations = network.get("configs", {})
        if not isinstance(configurations, dict) or set(configurations) - set(contrasts):
            raise ValueError("network.configs must map selected comparison IDs to their scoped configurations")
        if network_scope == "cohort" and configurations:
            raise ValueError("Per-comparison network.configs cannot be used for a cohort discovery")
        for comparison_ids in ([[cid] for cid in contrasts] if network_scope == "contrast" else [contrasts]):
            params = dict(network.get("parameters", {}))
            if "bootstrap_enabled" in params and (not isinstance(params["bootstrap_enabled"], bool)
                                                   or params["bootstrap_enabled"] != bootstrap_enabled):
                raise ValueError("network.bootstrap.enabled conflicts with network.parameters.bootstrap_enabled")
            params["bootstrap_enabled"] = bootstrap_enabled
            params.setdefault("hub_metric", "abs_kME")
            if params["hub_metric"] not in {"abs_kME", "legacy_composite"}:
                raise ValueError("network.parameters.hub_metric must be abs_kME or legacy_composite")
            if bootstrap_enabled:
                iterations = params.get("bootstrap_iterations", 50)
                if isinstance(iterations, bool) or not isinstance(iterations, int) or iterations < 2:
                    raise ValueError("Enabled network bootstrap_iterations must be an integer of at least 2")
                params["bootstrap_iterations"] = iterations
            if "network_scope" in params and params["network_scope"] != network_scope:
                raise ValueError("network.scope conflicts with network.parameters.network_scope")
            params.update(network_scope=network_scope, contrast_ids=comparison_ids)
            identity = network_scope_identity(manifest, params, preserve_recorded_status=False)
            params.update(identity)
            params["rmats_mode"] = resolve_network_rmats_mode(manifest, params)
            suffix = "-" + _token(comparison_ids[0]) if network_scope == "contrast" else ""
            fit_id = "network-fit" + suffix
            inputs = {"project": str(project_path)}
            configuration = configurations.get(comparison_ids[0]) if network_scope == "contrast" else None
            configuration = configuration or network.get("config")
            if configuration:
                inputs["config"] = local(configuration)
            scope = {**identity, "unit": "measurable_contrast_PSI" if network_scope == "contrast" else "measurable_cohort_PSI"}
            add(fit_id, "flow.network.fit", ["audit"], params, scope,
                adapter="network.fit", inputs=inputs)
            for stage in (("bootstrap", "associate", "annotate") if bootstrap_enabled else ("associate", "annotate")):
                add("network-" + stage + suffix, "flow.network." + stage, [fit_id],
                    params, {**identity, "unit": "module"}, adapter="network." + stage, inputs=inputs)

    if enabled["interaction"]:
        interactions = settings.get("interactions", [])
        if not isinstance(interactions, list):
            raise ValueError("interactions must contain explicit factorial designs")
        for spec in interactions:
            parameters = {k: v for k, v in spec.items() if k != "id"}
            if parameters.get("design_path"):
                parameters["design_path"] = local(parameters["design_path"])
            add("interaction-" + _token(spec["id"]), "differential.replicate_interaction", ["audit"], parameters,
                {"interaction_id": spec["id"]})

    if enabled["mapp"]:
        mapp = settings.get("mapp", {})
        mapp_inputs = resource_inputs(mapp.get("inputs", {}))
        from .bundled_runtime import mapp_input_defaults
        mapp_inputs = mapp_input_defaults(mapp_inputs)
        if mapp_inputs.get("inclusion_table"):
            if len(contrasts) != 1:
                raise ValueError("Explicit prepared MAEI inputs require one exact selected comparison; use project-derived inputs for multiple comparisons")
            add("mapp", "flow.mapp.maei", ["audit"], mapp.get("parameters", {}),
                {"contrast_ids": contrasts}, adapter="mapp.maei", inputs=mapp_inputs)
        else:
            from .mapp_inputs import native_windows, window_id
            samples_by_id = {s.sample_id: s for s in manifest.samples}
            for contrast in manifest.contrasts:
                if contrast.contrast_id not in contrasts:
                    continue
                sample_ids = list(contrast.numerator.sample_ids) + list(contrast.denominator.sample_ids)
                prepare = add("mapp-prepare-" + _token(contrast.contrast_id), "flow.mapp.prepare", ["as"],
                    {**mapp.get("prepare_parameters", {}), "all_windows": True, "sample_ids": sample_ids,
                     "comparison_id": contrast.contrast_id}, {"contrast_ids": [contrast.contrast_id]},
                    adapter="mapp.prepare", inputs={**mapp_inputs, "gtf": str(manifest.annotation.gtf or ""),
                    "genome_fasta": str(manifest.annotation.genome_fasta or ""),
                    "quantification": {sid: str(samples_by_id[sid].assets.get("quant_sf", "")) for sid in sample_ids}},
                    dependency_inputs={"events": {"analysis_id": "differential.as", "name": "Differential AS events"}})
                for window in native_windows():
                    token = _token(contrast.contrast_id) + "-" + window["site"] + "-" + window_id(window["window_start"], window["window_end"])
                    add("mapp-fit-" + token, "flow.mapp.maei", [prepare], mapp.get("parameters", {}),
                        {"contrast_ids": [contrast.contrast_id], **window}, adapter="mapp.maei", inputs=mapp_inputs)
    # Bind the code actually imported by this interpreter, including source
    # checkouts that use dependencies from a separately installed distribution.
    from katana_next import __version__ as package_version
    return seal_plan({"project_manifest": str(project_path), "project_manifest_sha256": file_sha256(project_path),
                      "project_id": manifest.project_id, "package_version": package_version,
                      "report_specification": "report-approved-20261022-01", "report_profile": "complete",
                      "contrast_ids": contrasts, "event_set": event_set,
                      "applicability": {"cross_comparison": {
                          "status": "applicable" if len(contrasts) >= 2 else "not_applicable",
                          "reason": "At least two selected comparisons are required for cross-comparison analysis.",
                          "contrast_ids": list(contrasts)}},
                      "configuration": settings, "tasks": tasks})


def preflight_plan(plan):
    """Read-only early checks. Parent artifacts are produced later by this DAG."""
    from .model import validate_plan
    from katana_next.infrastructure.manifest_store import ManifestStore
    from katana_next.workflows.registry import WorkflowRegistry
    validate_plan(plan)
    path = Path(plan["project_manifest"])
    errors, warnings = [], []
    if file_sha256(path) != plan["project_manifest_sha256"]:
        errors.append({"task_id": "project", "error": "Project manifest changed after planning"})
    _, manifest = ManifestStore().load(path)
    from .selection import analysis_selection, resource_kinds
    enabled = analysis_selection(plan.get("configuration", {}))
    resources = resource_kinds(plan.get("configuration", {}))
    kinds = {a.kind for a in manifest.assets if a.path.is_file()}
    for a in manifest.assets:
        if a.required and (a.kind not in {"gene_sets_gmt", "rbp_gene_sets", "rbp_motif_database", "fimo_runtime", "rbp_binding_peaks"} or a.kind in resources) and not a.path.is_file():
            errors.append({"task_id": "project", "error": f"Missing required asset {a.asset_id}: {a.path}"})
    for kind in {"rmats_result", "deg_result", "dtu_result", "dexseq_result"} | resources:
        if kind not in kinds:
            errors.append({"task_id": "project", "error": f"Fixed report requires input asset kind: {kind}"})
    if (enabled["motiftrace"] or enabled["motif_binding"]) and not manifest.metadata.get("species"):
        errors.append({"task_id": "motiftrace", "error": "Explicit project species is required"})
    references = ["gtf"] + (["genome_fasta"] if any(enabled[k] for k in ("motiftrace", "motif_binding", "mapp")) else [])
    for key in references:
        value = getattr(manifest.annotation, key)
        if value is None or not value.is_file():
            errors.append({"task_id": "reference", "error": f"Missing {key}"})
    factors = plan["configuration"].get("design_factors", [])
    if "design_factors" not in plan["configuration"]:
        errors.append({"task_id": "design", "error": "Declare design_factors explicitly; use ['condition'] for a single-factor study"})
    if enabled["interaction"] and len(factors) >= 2 and not plan["configuration"].get("interactions"):
        errors.append({"task_id": "interaction", "error": "Declared multifactor study needs explicit interaction designs; simple contrasts cannot replace fitted interactions"})
    if len(factors) < 2:
        warnings.append("No multifactor interaction declared; report must state the single-factor scope")
    registry = WorkflowRegistry()
    for task in plan["tasks"]:
        if task.get("adapter"):
            from .stages import readiness
            try:
                if task.get("adapter") == "historical.adopt":
                    from .adoption import validate_identity_projection
                    inputs = task.get("inputs", {})
                    if not inputs.get("adoption") or not inputs.get("identity_projection"):
                        raise ValueError("Historical adoption requires explicit source and study-identity projection")
                    validate_identity_projection(inputs["adoption"], inputs["identity_projection"], path)
                    if task.get("dependencies"):
                        raise ValueError("Historical input transport must not masquerade as scientific continuation")
                    continue
                deferred = list(task.get("dependency_inputs", {}))
                by_id = {t["id"]: t for t in plan["tasks"]}
                if task.get("adapter") == "mapp.maei" and any(by_id[d].get("adapter") == "mapp.prepare" for d in task["dependencies"]):
                    deferred.extend(["inclusion_table", "sitecount_matrix", "design", "quantity_provenance"])
                issues = readiness({**task, "readiness_deferred_inputs": deferred})
                if isinstance(issues, dict):
                    issues = issues.get("errors", issues.get("issues", []))
                errors.extend({"task_id": task["id"], "error": str(issue)} for issue in issues)
            except (ValueError, OSError, ImportError, KeyError) as exc:
                errors.append({"task_id": task["id"], "error": str(exc)})
        else:
            runner = registry.get(task["workflow_id"])
            for kind in runner.definition.required_asset_kinds:
                if kind not in kinds and not (kind == "gene_body_coverage" and any(t.get("adapter") == "qc.genebody" for t in plan["tasks"])):
                    errors.append({"task_id": task["id"], "error": f"Missing input asset kind: {kind}"})
            if task["workflow_id"] == "prioritize.motiftrace_parallel":
                from types import SimpleNamespace
                try:
                    # This existing read-only adapter verifies the manifest and every
                    # member SHA, species/build, exact case, direction and sample arms.
                    runner.input_paths(SimpleNamespace(manifest=manifest), task["parameters"])
                except (ValueError, OSError, KeyError) as exc:
                    errors.append({"task_id": task["id"], "error": str(exc)})
    return {"schema": "katana.flow.preflight/v1", "plan_id": plan["plan_id"], "ready": not errors,
            "errors": errors, "warnings": warnings, "task_count": len(plan["tasks"]),
            "scope": "Read-only input and adapter checks. Native completion and numerical acceptance occur after execution."}
