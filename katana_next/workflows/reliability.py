from pathlib import Path

import pandas as pd

from katana_cosplice.reliability_audit import calibration_grid, fixed_module_audit, thinning_diagnostic
from katana_next.workflows.base import ArtifactWriter, ParameterDefinition, WorkflowBase, WorkflowDefinition
from katana_next.workflows.bulk import _plot
from katana_next.workflows.sources import source_complete, validate_source, resolve_table


def structural_keys(frame):
    from katana_cosplice.matrix import EVENT_COORDINATES
    result = pd.Series(pd.NA, index=frame.index, dtype="string")
    if not {"event_type", "chr", "strand"} <= set(frame):
        return result
    for kind, coordinates in EVENT_COORDINATES.items():
        selected = frame.event_type.eq(kind)
        if not selected.any() or not set(coordinates) <= set(frame):
            continue
        part = frame.loc[selected]
        nums = part[coordinates].apply(pd.to_numeric, errors="coerce")
        valid = nums.notna().all(axis=1) & (nums >= 0).all(axis=1) & (nums == nums.round()).all(axis=1)
        valid &= part.chr.notna() & part.strand.isin(["+", "-"])
        part = part.loc[valid]
        keys = kind + "|" + part.chr.astype(str) + "|" + part.strand.astype(str)
        for column in coordinates:
            keys += "|" + nums.loc[valid, column].astype("int64").astype(str)
        result.loc[part.index] = keys
    return result


def diagnostic_sources(context, parameters):
    sources = context.catalog.list_artifacts(context.manifest.project_id, analysis_id="programs.cosplice", kind="table")
    chosen = str(parameters.get("matrix_artifact_id", "")).strip()
    if chosen:
        matrix = context.catalog.get_artifact(chosen)
        validate_source(context, matrix)
        if matrix.analysis_id != "programs.cosplice" or matrix.path.stem != "matrix__psi":
            raise ValueError("Choose the original CoSplice sample-level PSI matrix artifact")
    else:
        matrix = next((a for a in sources if a.path.stem == "matrix__psi"
            and a.metadata.get("project_snapshot_fingerprint") == context.snapshot.manifest_fingerprint
            and source_complete(context, a)), None)
    if matrix is None:
        raise ValueError("Run CoSplice discovery explicitly first. Diagnostics never launch discovery.")
    result = {a.path.stem: a for a in sources if a.run_fingerprint == matrix.run_fingerprint}
    for key in ("matrix__psi", "matrix__events", "matrix__samples", "network_observed__assignments"):
        if key not in result:
            raise ValueError(f"The selected fixed CoSplice run lacks {key}")
        validate_source(context, result[key])
    for a in result.values():
        validate_source(context, a)
    return result


def read_matrix(artifact):
    frame = pd.read_parquet(artifact.path).drop(columns="row_id", errors="ignore")
    key = "event_key" if "event_key" in frame else frame.columns[0]
    if key != "event_key" and not str(key).startswith("Unnamed"):
        raise ValueError(f"Matrix event identity column is ambiguous: {artifact.name}")
    frame[key] = frame[key].astype(str)
    return frame.set_index(key).apply(pd.to_numeric, errors="raise")


class ModuleDiagnosticsWorkflow(WorkflowBase):
    version = "1.2.2"
    definition = WorkflowDefinition(
        "programs.diagnostics", "programs", "Inspect fixed modules",
        "Are module patterns driven by condition, shared event structure, or robust event coverage?",
        "All-event structural components and sampled pair-complete condition-adjusted correlations; fixed membership, not a new network or significance test.",
        "Module evidence cards, event/pair audits, coverage and condition-adjusted comparison plots.",
        "fixed module and event pair", parameters=(
            ParameterDefinition("matrix_artifact_id", "Source CoSplice run", "artifact", "", "Empty uses the latest compatible complete discovery run; every source is frozen to that run.", options=("programs.cosplice", "2. Build PSI matrix: psi")),
            ParameterDefinition("max_pairs_per_module", "Maximum diagnostic pairs per module", "number", 5000, "Uniform reproducible pair sampling, not event filtering. All events are retained.", minimum=1, maximum=100000),
            ParameterDefinition("include_batch", "Adjust for declared batch", "boolean", False, "Optional; rejects missing or confounded batch metadata."),
            ParameterDefinition("seed", "Diagnostic seed", "number", 1701, "Does not change production module membership.", minimum=0),
        ), limitations="Both events are regressed on the same measured samples; deficient designs or residual df <3 remain missing. Connected structural groups are not counts of independent events. Shared junctions flag potential dependence, not proof of artefact. Pooled residual correlation is not a condition-specific network, rewiring test or causal result. Donor/block identity is only used when explicitly declared.")

    def parent_artifacts(self, context, parameters):
        parents = list(diagnostic_sources(context, parameters).values())
        das = resolve_table(context, "differential.as", "Differential AS events", parameters)
        return parents + ([das] if das else [])

    def run(self, context, job, run_fingerprint, parameters):
        bound = parameters.get("source_artifact_ids", {}).get("programs.cosplice::2. Build PSI matrix: psi")
        sources = diagnostic_sources(context, {**parameters, "matrix_artifact_id": bound or parameters.get("matrix_artifact_id", "")})
        psi = read_matrix(sources["matrix__psi"])
        events = pd.read_parquet(sources["matrix__events"].path)
        samples = pd.read_parquet(sources["matrix__samples"].path)
        assignments = pd.read_parquet(sources.get("consensus_modules__assignments", sources["network_observed__assignments"]).path)
        if "module" not in assignments:
            column = "analysis_module" if "analysis_module" in assignments else "consensus_module"
            assignments = assignments.rename(columns={column: "module"})
        job.report(.1, "Inspecting the frozen event universe without changing module membership")
        result = fixed_module_audit(psi, samples, events, assignments,
            max_pairs=int(parameters["max_pairs_per_module"]), seed=int(parameters["seed"]), include_batch=bool(parameters["include_batch"]))
        writer = ArtifactWriter(context, self.definition.workflow_id, run_fingerprint, parameters, self.version,
            common_metadata={"method": self.definition.method, "limitations": self.definition.limitations,
                             "source_discovery_run": sources["matrix__psi"].run_fingerprint})
        parents = tuple(a.artifact_id for a in sources.values())
        das = resolve_table(context, "differential.as", "Differential AS events", parameters)
        if das:
            parents += (das.artifact_id,)
            membership = events.merge(assignments[["event_key", "module"]], on="event_key", how="left", validate="one_to_one")
            membership["structure_key"] = structural_keys(membership)
            bulk = pd.read_parquet(das.path)
            bulk["structure_key"] = structural_keys(bulk)
            mapped = bulk.dropna(subset=["structure_key"]).merge(
                membership.dropna(subset=["structure_key"])[["structure_key", "event_key", "module"]], on="structure_key", how="inner")
            mapped["row_id"] = [f"module-validation:{i}" for i in range(len(mapped))]
            writer.table("Module events ready for Native Sashimi", "validation_events.parquet", mapped,
                entity_type="splice_event", parents=parents,
                metadata={"matching": "Exact event-type/chromosome/strand/coordinates in the same project genome build; not gene-symbol or proximity matching"})
            writer.table("Module-to-DAS mapping audit", "event_mapping.parquet",
                membership[["event_key", "module", "structure_key"]].assign(
                    mapped_to_DAS=membership.event_key.isin(set(mapped.event_key))), entity_type="event_mapping", parents=parents)
        for column in ("consensus_membership_confidence", "fixed_module_kme_abs_mean"):
            if column in assignments:
                mapped_conf = assignments.set_index("event_key")[column]
                result["events"][column] = result["events"].event_key.map(mapped_conf)
        if "consensus_membership_confidence" in assignments:
            view = assignments.copy()
            view["confidence"] = pd.to_numeric(view.consensus_membership_confidence, errors="coerce")
            sensitivity = []
            assigned = view[~view.module.isin(["grey", "gray", "0", "unassigned", ""])]
            for cutoff in (None, 0., .5, .7, .9):
                kept = assigned if cutoff is None else assigned[assigned.confidence.ge(cutoff)]
                sensitivity.append({"diagnostic_cutoff": cutoff, "view": "All assigned" if cutoff is None else f"Confidence >= {cutoff:.1f}", "retained_events": len(kept),
                    "input_event_universe": len(psi), "coverage_fraction": len(kept)/max(1,len(psi)),
                    "represented_modules": kept.module.nunique(), "lost_modules": assigned.module.nunique()-kept.module.nunique(),
                    "missing_confidence_events": int(assigned.confidence.isna().sum()),
                    "meaning": "What-if display only; no cutoff selected and no production membership altered"})
            writer.table("Confidence versus event coverage sensitivity", "confidence_coverage.parquet", pd.DataFrame(sensitivity), entity_type="sensitivity", parents=parents)
            writer.plot("Confidence must be read alongside coverage", "confidence_coverage.vl.json",
                _plot(sensitivity, "bar", {"x": {"field": "view", "type": "ordinal", "sort": [r["view"] for r in sensitivity], "axis": {"labelAngle": -25}},
                    "y": {"field": "coverage_fraction", "type": "quantitative", "scale": {"domain": [0, 1]}, "title": "Assigned events / input events"},
                    "tooltip": [{"field": c} for c in sensitivity[0]]}, title="How much of the event universe remains?",
                    description="What-if display, not a new filter. Missing confidence is reported, not converted to zero. Check lost modules as well as retained events."),
                entity_type="module_diagnostics", parents=parents)
        for name, frame in result.items():
            if not len(frame.columns):
                frame = pd.DataFrame(columns=["module", "status"])
            frame = frame.copy()
            frame.insert(0, "row_id", [f"fixed-{name}:{i}" for i in range(len(frame))])
            writer.table(f"Fixed-module {name}", f"{name}.parquet", frame,
                entity_type="cosplice_module" if name == "modules" else "event_diagnostic" if name == "events" else "analysis_summary", parents=parents)
        modules = result["modules"]
        if len(modules):
            display = modules.melt(id_vars=["module", "n_events"], value_vars=["mean_observed_r_on_adjusted_pairs", "mean_condition_adjusted_r", "mean_adjusted_r_without_known_dependency"], var_name="measurement", value_name="mean_r")
            display["measurement"] = display.measurement.map({"mean_observed_r_on_adjusted_pairs": "Observed, matched pairs", "mean_condition_adjusted_r": "Condition-adjusted",
                "mean_adjusted_r_without_known_dependency": "Adjusted, dependency-screened"})
            display["module_label"] = display.module + " (n=" + display.n_events.astype(str) + ")"
            writer.plot("Observed versus condition-adjusted module coherence", "condition_coherence.vl.json",
                _plot(display.where(pd.notna(display), None).to_dict("records"), {"type": "point", "filled": True, "size": 85},
                    {"x": {"field": "mean_r", "type": "quantitative", "scale": {"domain": [-1, 1]}, "axis": {"tickCount": 5}, "title": "Mean pairwise Pearson r (not a significance test)"},
                     "y": {"field": "module_label", "type": "nominal", "title": "Frozen module (event count)"},
                     "yOffset": {"field": "measurement", "type": "nominal"},
                     "color": {"field": "measurement", "type": "nominal", "legend": {"orient": "bottom", "title": None, "labelLimit": 0}},
                     "tooltip": [{"field": c} for c in display.columns]},
                    title="Module correlation before and after condition adjustment",
                    description="Observed and adjusted means use the same evaluable pairs. Fixed modules; pair-complete pooled residual correlations, not condition-specific networks. The dependency-screened subset excludes shared junction/denominator and numerical complements; it is not proof of independence. Undefined correlations stay missing."),
                entity_type="module_diagnostics", parents=parents)
        job.report(.95, "Registering all module/event audits and explicit pair-sampling denominators")
        return writer.finish()


class ReadCalibrationWorkflow(WorkflowBase):
    version = "1.3.1"
    definition = WorkflowDefinition(
        "quality.read_calibration", "quality", "Calibrate read uncertainty",
        "Do uncertainty intervals cover the intended latent read fraction or length-normalized sample PSI?",
        "Known-truth simulations compare the effective-count approximation with a Jeffreys baseline; a separate length-aware calibration transforms sample read-fraction quantiles using the rMATS effective-length formula.",
        "Coverage, interval width, error and Monte Carlo uncertainty across depth/dispersion; optional count thinning and training-only gene-crossfit prediction of held biological replicates.",
        "simulated event and event/sample read fraction", parameters=(
            ParameterDefinition("simulated_events", "Independent simulated events per scenario", "number", 200, "Each event contributes one independent Monte Carlo unit; no adaptive tuning.", minimum=10, maximum=10000),
            ParameterDefinition("replicates", "Replicates per simulated condition", "number", 3, "Simulation design, not the number of actual biological replicates.", minimum=2, maximum=50),
            ParameterDefinition("seed", "Simulation seed", "number", 1701, "Frozen before comparing methods.", minimum=0),
            ParameterDefinition("include_count_thinning", "Also inspect existing CoSplice raw counts", "boolean", False, "Requires a complete cached discovery run. Does not rerun upstream analyses."),
            ParameterDefinition("include_heldout_replicates", "Predict held-out biological replicates", "boolean", False, "Uses completed Quality counts, never significant-event selection. Excludes held samples and test-fold genes from dispersion fitting."),
            ParameterDefinition("heldout_event_cap", "Validation event cap per comparison (0 = all)", "number", 0, "Optional explicit compute budget, seeded random events; retained/input denominators are exported. Never a biological threshold.", minimum=0),
            ParameterDefinition("predictive_draws", "Predictive Monte Carlo draws", "number", 2000, "Fixed-seed 95% predictive intervals for an observed held replicate, not latent biological truth.", minimum=100, maximum=10000),
            ParameterDefinition("prediction_model", "Held-repeat interval methods", "select", "legacy", "Hierarchical is an opt-in diagnostic candidate, not a validated replacement: real held-repeat testing found excessive width and worse interval scores. Compare coverage and width together; no production thresholds are changed.", options=("all", "legacy", "hierarchical")),
            ParameterDefinition("matrix_artifact_id", "Source CoSplice run", "artifact", "", "Select a complete fixed run for optional thinning.", options=("programs.cosplice", "2. Build PSI matrix: psi")),
        ), limitations="A simulation under specified models is not proof of universal calibration. Raw read fractions and length-normalized PSI are reported separately. Length-aware intervals target individual-sample measurement uncertainty, not a condition-mean beta-binomial model. Split reads are not independent biological replicates.")

    def parent_artifacts(self, context, parameters):
        parents = list(diagnostic_sources(context, parameters).values()) if parameters.get("include_count_thinning") else []
        if parameters.get("include_heldout_replicates"):
            source = resolve_table(context, "quality.stable_as", "Stable-AS event QC", parameters)
            if source is None:
                raise ValueError("Complete Quality explicitly before held-repeat prediction; no hidden upstream run")
            parents.append(source)
        return parents

    def run(self, context, job, run_fingerprint, parameters):
        job.report(.1, "Running predeclared known-truth scenarios; no real-data performance claim")
        result = calibration_grid(events=int(parameters["simulated_events"]), seed=int(parameters["seed"]), replicates=int(parameters["replicates"]))
        writer = ArtifactWriter(context, self.definition.workflow_id, run_fingerprint, parameters, self.version,
            common_metadata={"method": self.definition.method, "limitations": self.definition.limitations,
                "model_identity": "CoSplice method-of-moments overdispersion/effective-count approximation, not a full beta-binomial likelihood fit"})
        writer.table("Known-truth read-fraction calibration", "calibration.parquet", result, entity_type="calibration_scenario")
        from katana_next.science.psi_calibration import length_calibration_grid
        normalized = length_calibration_grid(events=int(parameters["simulated_events"]),
            seed=int(parameters["seed"]), replicates=int(parameters["replicates"]))
        writer.table("Length-normalized PSI calibration", "length_calibration.parquet", normalized,
            entity_type="calibration_scenario", metadata={"method": "Monotone rMATS effective-length transform of Jeffreys Beta quantiles; individual-sample count-conditional uncertainty", "source": "https://github.com/Xinglab/rmats-turbo/blob/master/rMATS_P/inclusion_level.py", "limitations": "Does not calibrate a condition-mean beta-binomial model; no biological superiority claim."})
        writer.plot("Length-normalized PSI interval coverage", "length_coverage.vl.json", _plot(normalized.to_dict("records"), {"type": "point", "filled": True},
            {"x": {"field": "depth", "type": "ordinal"}, "y": {"field": "observed_coverage", "type": "quantitative", "scale": {"domain": [0, 1]}},
             "color": {"field": "length_ratio", "type": "nominal"}, "tooltip": [{"field": c} for c in normalized.columns]},
            title="Individual-sample PSI: length-aware coverage", description="Known latent PSI simulations; examine coverage, width and error together."), entity_type="calibration_plot")
        writer.plot("Interval coverage by read depth", "coverage.vl.json", _plot(result.to_dict("records"), {"type": "line", "point": True},
            {"x": {"field": "depth", "type": "quantitative", "scale": {"type": "log"}},
             "y": {"field": "observed_coverage", "type": "quantitative", "scale": {"domain": [0, 1]}, "title": "Empirical coverage of nominal 95% intervals"},
             "color": {"field": "method", "type": "nominal"}, "strokeDash": {"field": "true_rho", "type": "nominal"},
             "tooltip": [{"field": c} for c in result.columns]}, title="Known-truth calibration, not biological validation",
             description="Uncertainty is measured across independent simulated events; compare coverage together with width and error."), entity_type="calibration_plot")
        if parameters.get("include_count_thinning"):
            job.check_cancelled()
            bound = parameters.get("source_artifact_ids", {}).get("programs.cosplice::2. Build PSI matrix: psi")
            sources = diagnostic_sources(context, {**parameters, "matrix_artifact_id": bound or parameters.get("matrix_artifact_id", "")})
            for key in ("matrix__ijc", "matrix__sjc"):
                if key not in sources:
                    raise ValueError("Read thinning requires original raw IJC/SJC matrices")
            thinning = thinning_diagnostic(read_matrix(sources["matrix__ijc"]), read_matrix(sources["matrix__sjc"]),
                pd.read_parquet(sources["matrix__samples"].path), seed=int(parameters["seed"]))
            writer.table("Disjoint-read reproducibility (not independent validation)", "thinning.parquet", thinning,
                entity_type="read_thinning", parents=tuple(a.artifact_id for a in sources.values()))
        if parameters.get("include_heldout_replicates"):
            from katana_next.science.replicate_prediction import heldout_arm, prediction_summary
            source = resolve_table(context,"quality.stable_as","Stable-AS event QC",parameters)
            all_events = pd.read_parquet(source.path)
            samples = {s.sample_id:s for s in context.manifest.samples}
            summaries, audits = [], []
            for ci, contrast in enumerate(context.manifest.contrasts):
                part = all_events[all_events.comparison_id.eq(contrast.contrast_id)].copy()
                if part.empty:
                    continue
                total = len(part)
                cap = int(parameters["heldout_event_cap"])
                if cap and len(part) > cap:
                    part = part.sort_values("event_id").sample(n=cap,random_state=int(parameters["seed"]))
                for arm in ("numerator","denominator"):
                    ids = tuple(getattr(contrast,arm).sample_ids)
                    # Mapping metadata can override raw vector order; require explicit agreement.
                    assets = [a for a in context.manifest.assets if a.kind == "rmats_result" and str(a.path) in set(part.source_path)]
                    if any(tuple(a.metadata.get(f"{arm}_sample_ids",ids)) != ids for a in assets):
                        raise ValueError("Held-repeat prediction requires source arm order identical to declared contrast order; revise mapping first")
                    reverse = part.reverse_direction.unique()
                    if len(reverse) != 1:
                        raise ValueError("Mixed source directions in one comparison")
                    index = 1 if (arm == "numerator") != bool(reverse[0]) else 2
                    def progress(sample, fold, begin, count):
                        job.check_cancelled()
                        if begin == 0:
                            job.report(.5+.4*ci/max(1,len(context.manifest.contrasts)), f"Held-repeat prediction: {contrast.label} | {sample} | gene fold {fold+1}/5")
                    detail, fits = heldout_arm(part,ids,[samples[s].biological_sample_id for s in ids],index=index,
                        seed=int(parameters["seed"]),draws=int(parameters["predictive_draws"]),report=progress,
                        model=parameters.get("prediction_model", "legacy"))
                    for frame in (detail,fits):
                        frame["comparison_id"],frame["arm"] = contrast.contrast_id,arm
                    writer.table(f"Held-repeat predictions: {contrast.label} | {arm}", f"heldout_{ci}_{arm}.parquet",detail,
                        entity_type="heldout_prediction",parents=(source.artifact_id,),metadata={"input_events":total,"evaluated_event_selection":len(part),
                        "limitations":"Gene-crossfit training-only plug-in dispersion, not full Bayesian inference; Monte Carlo intervals; no significance threshold is changed",
                        "hierarchical_mean_prior":"Beta(0.5,0.5)", "quadrature_nodes":256,
                        "point_estimator":"Pooled-count Jeffreys mean, not the hierarchical posterior mean"})
                    writer.table(f"Training dispersion audit: {contrast.label} | {arm}",f"heldout_fit_{ci}_{arm}.parquet",fits,
                        entity_type="calibration_fit",parents=(source.artifact_id,))
                    summaries.append(prediction_summary(detail))
                    audits.append({"comparison_id":contrast.contrast_id,"arm":arm,"input_events":total,"selected_events":len(part),"selection":"all" if len(part)==total else "explicit_seeded_random_compute_cap"})
            if summaries:
                summary=pd.concat(summaries,ignore_index=True)
                writer.table("Held biological replicate prediction summary","heldout_summary.parquet",summary,entity_type="calibration_summary",parents=(source.artifact_id,))
                writer.table("Held-repeat input coverage","heldout_selection.parquet",pd.DataFrame(audits),entity_type="analysis_coverage",parents=(source.artifact_id,))
                writer.plot("Held-repeat coverage and width","heldout_coverage.vl.json",_plot(summary.to_dict("records"),{"type":"point","filled":True},
                    {"x":{"field":"mean_interval_width","type":"quantitative","title":"Mean 95% predictive interval width"},
                     "y":{"field":"observed_predictive_coverage","type":"quantitative","scale":{"domain":[0,1]}},
                     "color":{"field":"method","type":"nominal"},"tooltip":[{"field":c} for c in summary.columns]},
                    title="Held biological repeats: coverage must be read together with width",description="Test-sample counts never fit dispersion. Observed replicate prediction is not latent PSI truth or a new significance test."),entity_type="calibration_plot",parents=(source.artifact_id,))
        return writer.finish()
