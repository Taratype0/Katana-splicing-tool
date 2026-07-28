from __future__ import annotations

from datetime import datetime, timezone
from html import escape
from pathlib import Path
from typing import Callable

import numpy as np
import pandas as pd

from katana_cosplice import __version__
from katana_cosplice.association import (
    associate_program_scores,
    factorial_program_association,
    render_association_plot,
    render_factorial_association_plot,
    summarize_program_scores,
)
from katana_cosplice.artifacts import ArtifactManifestStore
from katana_cosplice.audit import audit_inputs, resolve_rmats_root
from katana_cosplice.benchmark import render_benchmark_plots, render_benchmark_report, run_method_benchmark
from katana_cosplice.config import CoSpliceConfig
from katana_cosplice.factor_contrasts import analyze_factor_contrasts, render_factor_contrast_plots
from katana_cosplice.consensus import (
    build_bootstrap_coassignment_consensus,
    render_consensus_plots,
    save_coassignment_matrix,
)
from katana_cosplice.hybrid_programs import (
    add_external_preservation_confidence,
    build_hybrid_programs,
    render_hybrid_program_plots,
)
from katana_cosplice.matrix import build_sample_matrices
from katana_cosplice.meta_analysis import meta_analyze_modules, render_meta_analysis_plots
from katana_cosplice.module_annotation import annotate_modules, render_module_annotation_plot
from katana_cosplice.module_interpretation import interpret_modules, render_module_interpretation_plots
from katana_cosplice.module_models import fit_replicate_aware_models, render_module_statistics_plots
from katana_cosplice.network import fit_network, render_network_plots, residualize_matrix
from katana_cosplice.qc import stable_as_qc
from katana_cosplice.preservation import (
    PreservationResult,
    evaluate_external_module_preservation,
    render_preservation_plots,
)
from katana_cosplice.reporting import build_module_browser, build_validation_pdf, read_cached_table
from katana_cosplice.stage_registry import stage_spec
from katana_cosplice.utils import file_signature, fingerprint, write_json
from katana_cosplice.validation import render_validation_plots, render_validation_report, validate_network_modules


Progress = Callable[[str], None] | None


class CoSplicePipeline:
    def __init__(self, config: CoSpliceConfig, progress: Progress = None) -> None:
        self.config = config
        self.progress = progress
        self.output = config.output_dir
        self.output.mkdir(parents=True, exist_ok=True)
        self.state_path = self.output / "state.json"
        self.artifacts = ArtifactManifestStore(self.output, config.project_name)

    def _emit(self, message: str) -> None:
        if self.progress:
            self.progress(message)
        print(f"[katana-cosplice] {message}")

    def _state(self) -> dict[str, object]:
        if not self.state_path.exists():
            return {"project_name": self.config.project_name, "stages": {}}
        import json
        return json.loads(self.state_path.read_text(encoding="utf-8"))

    def _mark(self, stage: str, status: str, message: str, stage_fingerprint: str = "", outputs: list[Path] | None = None) -> None:
        stage_spec(stage)
        state = self._state()
        stages = state.setdefault("stages", {})
        previous = stages.get(stage, {})
        stages[stage] = {
            "status": status,
            "message": message,
            "fingerprint": stage_fingerprint,
            "updated_at": datetime.now(timezone.utc).isoformat(),
            "outputs": [str(path) for path in outputs] if outputs is not None else previous.get("outputs", []),
        }
        write_json(self.state_path, state)
        self.artifacts.update_stage(stage, status, message, stage_fingerprint, outputs)

    def status(self) -> dict[str, object]:
        return self._state()

    def run_stage(self, stage: str, force: bool = False) -> object:
        """Run exactly one named stage; missing dependencies are never auto-executed."""
        spec = stage_spec(stage)
        if not spec.implemented:
            raise NotImplementedError(f"CoSplice stage is not implemented: {stage}")
        actions: dict[str, Callable[[], object]] = {
            "audit": self.audit,
            "matrix": lambda: self.build_matrix(force=force),
            "stable_qc": lambda: self.stable_qc(force=force),
            "network_observed": lambda: self.fit_network("observed", force=force),
            "network_residual": lambda: self.fit_network("residual", force=force),
            "consensus_modules": lambda: self.consensus_modules(force=force),
            "hybrid_programs": lambda: self.hybrid_programs(force=force),
            "preserve_modules": lambda: self.preserve_modules(force=force),
            "associate_observed": lambda: self.associate("observed", force=force),
            "associate_residual": lambda: self.associate("residual", force=force),
            "validate_modules": lambda: self.validate_modules(force=force),
            "annotate_modules": lambda: self.annotate_modules(force=force),
            "benchmark_methods": lambda: self.benchmark_methods(force=force),
            "interpret_modules": lambda: self.interpret_modules(force=force),
            "fit_module_models": lambda: self.fit_module_models(force=force),
            "meta_analyze_modules": lambda: self.meta_analyze_modules(force=force),
            "factor_contrasts": lambda: self.factor_contrasts(force=force),
            "report": self.report,
        }
        return actions[stage]()

    def _fingerprint(self, stage: str, inputs: list[Path], settings: object) -> str:
        return fingerprint(
            {
                "schema": 1,
                "version": __version__,
                "stage": stage,
                "settings": settings,
                "inputs": [file_signature(path) for path in inputs],
            }
        )

    def audit(self) -> pd.DataFrame:
        stage = "audit"
        self._mark(stage, "running", "Auditing input files and sample identity.")
        try:
            result = audit_inputs(self.config)
            out_dir = self.output / "01_audit"
            out_dir.mkdir(parents=True, exist_ok=True)
            checks_path = out_dir / "input_audit.tsv"
            samples_path = out_dir / "sample_metadata_resolved.tsv"
            result.checks.to_csv(checks_path, sep="\t", index=False)
            result.samples.to_csv(samples_path, sep="\t", index=False)
            params = out_dir / "parameters.json"
            write_json(params, {"project": self.config.project_name, "passed": result.passed, "config": self.config.to_mapping()})
            status = "finished" if result.passed else "failed"
            self._mark(stage, status, "Input audit completed." if result.passed else "Input audit found blocking errors.", outputs=[checks_path, samples_path, params])
            return result.checks
        except Exception as exc:
            self._mark(stage, "failed", str(exc))
            raise

    def build_matrix(self, force: bool = False) -> dict[str, Path]:
        stage = "matrix"
        out_dir = self.output / "02_matrix"
        outputs = {
            "psi": out_dir / "sample_level_psi_matrix.tsv",
            "ijc": out_dir / "sample_level_ijc_matrix.tsv",
            "sjc": out_dir / "sample_level_sjc_matrix.tsv",
            "events": out_dir / "event_metadata.tsv",
            "samples": out_dir / "sample_metadata_resolved.tsv",
            "sources": out_dir / "rmats_source_selection.tsv",
            "integrity": out_dir / "matrix_integrity_summary.tsv",
            "parameters": out_dir / "parameters.json",
        }
        rmats_root = resolve_rmats_root(self.config.rmats_root)
        tables = list(rmats_root.glob(f"**/*.MATS.{self.config.rmats_mode}.txt"))
        bamlists = list((rmats_root / "bamlist").glob("*_bamlist.txt"))
        stage_fp = self._fingerprint(
            stage,
            [self.config.sample_metadata, *sorted(tables), *sorted(bamlists)],
            {"mode": self.config.rmats_mode, "conditions": self.config.selected_conditions},
        )
        state = self._state().get("stages", {}).get(stage, {})
        if not force and state.get("fingerprint") == stage_fp and all(path.exists() for path in outputs.values()):
            self._mark(stage, "cache_hit", "Loaded existing sample-level matrix cache.", stage_fp, list(outputs.values()))
            return outputs
        self._mark(stage, "running", "Building sample-level PSI/IJC/SJC matrices.", stage_fp)
        try:
            result = build_sample_matrices(self.config)
            out_dir.mkdir(parents=True, exist_ok=True)
            result.psi.to_csv(outputs["psi"], sep="\t", float_format="%.6g")
            result.ijc.to_csv(outputs["ijc"], sep="\t", float_format="%.6g")
            result.sjc.to_csv(outputs["sjc"], sep="\t", float_format="%.6g")
            result.events.to_csv(outputs["events"], sep="\t", index=False)
            result.samples.to_csv(outputs["samples"], sep="\t", index=False)
            result.sources.to_csv(outputs["sources"], sep="\t", index=False)
            integrity = pd.DataFrame(
                [
                    {"metric": "n_events", "value": result.psi.shape[0]},
                    {"metric": "n_samples", "value": result.psi.shape[1]},
                    {"metric": "n_conditions", "value": result.samples["condition"].nunique()},
                    {"metric": "psi_missing_fraction", "value": float(result.psi.isna().mean().mean())},
                    {"metric": "psi_min", "value": float(result.psi.min(skipna=True).min(skipna=True))},
                    {"metric": "psi_max", "value": float(result.psi.max(skipna=True).max(skipna=True))},
                    {"metric": "ijc_min", "value": float(result.ijc.min(skipna=True).min(skipna=True))},
                    {"metric": "sjc_min", "value": float(result.sjc.min(skipna=True).min(skipna=True))},
                    {"metric": "unique_event_keys", "value": bool(result.psi.index.is_unique)},
                    {"metric": "unique_sample_ids", "value": bool(result.psi.columns.is_unique)},
                ]
            )
            integrity.to_csv(outputs["integrity"], sep="\t", index=False)
            write_json(
                outputs["parameters"],
                {
                    "fingerprint": stage_fp,
                    "rmats_mode": self.config.rmats_mode,
                    "n_events": result.psi.shape[0],
                    "n_samples": result.psi.shape[1],
                    "selected_conditions": self.config.selected_conditions,
                },
            )
            self._mark(stage, "finished", f"Built {result.psi.shape[0]} events x {result.psi.shape[1]} samples.", stage_fp, list(outputs.values()))
            return outputs
        except Exception as exc:
            self._mark(stage, "failed", str(exc), stage_fp)
            raise

    def stable_qc(self, force: bool = False) -> dict[str, Path]:
        stage = "stable_qc"
        matrix_dir = self.output / "02_matrix"
        required = [matrix_dir / "sample_level_psi_matrix.tsv", matrix_dir / "sample_level_ijc_matrix.tsv", matrix_dir / "sample_level_sjc_matrix.tsv", matrix_dir / "sample_metadata_resolved.tsv"]
        if not all(path.exists() for path in required):
            raise FileNotFoundError("Missing input: sample-level matrix not built. Run build-matrix first.")
        out_dir = self.output / "03_stable_as_qc"
        outputs = {
            "metrics": out_dir / "stable_as_event_qc.tsv",
            "condition_metrics": out_dir / "stable_as_condition_qc.tsv",
            "summary": out_dir / "stable_as_qc_summary.tsv",
            "matrix": out_dir / "stable_psi_matrix.tsv",
            "parameters": out_dir / "parameters.json",
        }
        thresholds = self.config.stable_as.resolved()
        stage_fp = self._fingerprint(stage, required, thresholds)
        state = self._state().get("stages", {}).get(stage, {})
        if not force and state.get("fingerprint") == stage_fp and all(path.exists() for path in outputs.values()):
            self._mark(stage, "cache_hit", "Loaded existing stable-AS QC cache.", stage_fp, list(outputs.values()))
            return outputs
        self._mark(stage, "running", "Calculating read-supported replicate stability.", stage_fp)
        try:
            psi = pd.read_csv(required[0], sep="\t", index_col=0)
            ijc = pd.read_csv(required[1], sep="\t", index_col=0)
            sjc = pd.read_csv(required[2], sep="\t", index_col=0)
            samples = pd.read_csv(required[3], sep="\t")
            result = stable_as_qc(psi, ijc, sjc, samples, self.config.stable_as)
            out_dir.mkdir(parents=True, exist_ok=True)
            result.metrics.to_csv(outputs["metrics"], sep="\t", float_format="%.6g")
            result.condition_metrics.to_csv(outputs["condition_metrics"], sep="\t", index=False, float_format="%.6g")
            result.summary.to_csv(outputs["summary"], sep="\t", index=False, float_format="%.6g")
            result.stable_psi.to_csv(outputs["matrix"], sep="\t", float_format="%.6g")
            write_json(
                outputs["parameters"],
                {
                    "fingerprint": stage_fp,
                    "thresholds": thresholds,
                    "n_input_events": len(result.metrics),
                    "n_stable_events": len(result.stable_psi),
                    "condition_detail_scope": "stable_events_only",
                },
            )
            self._mark(stage, "finished", f"Stable-AS QC retained {len(result.stable_psi)} events.", stage_fp, list(outputs.values()))
            return outputs
        except Exception as exc:
            self._mark(stage, "failed", str(exc), stage_fp)
            raise

    def fit_network(self, network_name: str = "observed", force: bool = False) -> dict[str, Path]:
        if network_name not in {"observed", "residual"}:
            raise ValueError("network_name must be observed or residual")
        stage = f"network_{network_name}"
        qc_dir = self.output / "03_stable_as_qc"
        matrix_path = qc_dir / "stable_psi_matrix.tsv"
        metrics_path = qc_dir / "stable_as_event_qc.tsv"
        samples_path = self.output / "02_matrix" / "sample_metadata_resolved.tsv"
        required = [matrix_path, metrics_path, samples_path]
        if not all(path.exists() for path in required):
            raise FileNotFoundError("Missing input: stable-AS QC not run. Run qc first.")
        out_dir = self.output / "04_network" / network_name
        outputs = {
            "assignments": out_dir / "event_module_assignments.tsv",
            "edges": out_dir / "network_edges.tsv",
            "scores": out_dir / "module_program_scores.tsv",
            "diagnostics": out_dir / "network_diagnostics.tsv",
            "cut_sensitivity": out_dir / "module_resolution_sensitivity.tsv",
            "selected": out_dir / "selected_event_matrix.tsv",
            "parameters": out_dir / "parameters.json",
        }
        settings = self.config.network.__dict__ | {"network_name": network_name, "residual_covariates": self.config.residual_covariates}
        stage_fp = self._fingerprint(stage, required, settings)
        state = self._state().get("stages", {}).get(stage, {})
        if not force and state.get("fingerprint") == stage_fp and all(path.exists() for path in outputs.values()):
            self._mark(stage, "cache_hit", f"Loaded existing {network_name} network cache.", stage_fp, list(outputs.values()))
            return outputs
        self._mark(stage, "running", f"Fitting {network_name} co-splicing network.", stage_fp)
        try:
            matrix = pd.read_csv(matrix_path, sep="\t", index_col=0)
            metrics = pd.read_csv(metrics_path, sep="\t", index_col=0)
            samples = pd.read_csv(samples_path, sep="\t")
            if network_name == "residual":
                matrix = residualize_matrix(matrix, samples, self.config.residual_covariates)
            result = fit_network(matrix, metrics, self.config.network, samples=samples)
            out_dir.mkdir(parents=True, exist_ok=True)
            result.assignments.to_csv(outputs["assignments"], sep="\t", index=False)
            result.edges.to_csv(outputs["edges"], sep="\t", index=False, float_format="%.6g")
            result.program_scores.to_csv(outputs["scores"], sep="\t", float_format="%.6g")
            result.diagnostics.to_csv(outputs["diagnostics"], sep="\t", index=False)
            result.cut_height_sensitivity.to_csv(outputs["cut_sensitivity"], sep="\t", index=False)
            result.selected_matrix.to_csv(outputs["selected"], sep="\t", float_format="%.6g")
            plots = render_network_plots(result, samples, out_dir / "plots", network_name)
            write_json(outputs["parameters"], {"fingerprint": stage_fp, "settings": settings, "generated_plots": [str(path) for path in plots]})
            self._mark(stage, "finished", f"Finished {network_name} network with {result.assignments['module'].nunique()} labels.", stage_fp, [*outputs.values(), *plots])
            return outputs
        except Exception as exc:
            self._mark(stage, "failed", str(exc), stage_fp)
            raise

    def associate(self, network_name: str = "observed", force: bool = False) -> dict[str, Path]:
        stage = f"associate_{network_name}"
        scores_path = self.output / "04_network" / network_name / "module_program_scores.tsv"
        samples_path = self.output / "02_matrix" / "sample_metadata_resolved.tsv"
        if not scores_path.exists():
            raise FileNotFoundError(f"Missing input: {network_name} network not fit. Run fit-network --network {network_name} first.")
        out_dir = self.output / "05_association" / network_name
        outputs = {
            "table": out_dir / "module_trait_association.tsv",
            "plot": out_dir / "module_trait_association.png",
            "factorial": out_dir / "factorial_module_association.tsv",
            "coefficients": out_dir / "factorial_module_coefficients.tsv",
            "group_summary": out_dir / "module_program_score_group_summary.tsv",
            "design": out_dir / "association_design.tsv",
            "factorial_plot": out_dir / "factorial_module_association.png",
            "parameters": out_dir / "parameters.json",
        }
        factors = list(dict.fromkeys(self.config.association_factors))
        association_settings = {
            "factors": factors,
            "factorial_factors": self.config.factorial_factors,
            "batch_column": self.config.batch_column,
        }
        stage_fp = self._fingerprint(stage, [scores_path, samples_path], association_settings)
        state = self._state().get("stages", {}).get(stage, {})
        if not force and state.get("fingerprint") == stage_fp and outputs["table"].exists() and outputs["parameters"].exists():
            self._mark(stage, "cache_hit", f"Loaded existing {network_name} association cache.", stage_fp, list(outputs.values()))
            return outputs
        self._mark(stage, "running", f"Associating {network_name} program scores with sample factors.", stage_fp)
        try:
            scores = pd.read_csv(scores_path, sep="\t", index_col=0)
            samples = pd.read_csv(samples_path, sep="\t")
            table = associate_program_scores(scores, samples, factors)
            factorial, coefficients = factorial_program_association(
                scores,
                samples,
                self.config.factorial_factors,
                self.config.batch_column,
            )
            group_summary = summarize_program_scores(scores, samples, factors)
            design_columns = list(
                dict.fromkeys(
                    [
                        "sample_id",
                        *factors,
                        *self.config.factorial_factors,
                        *([self.config.batch_column] if self.config.batch_column else []),
                    ]
                )
            )
            design = samples.loc[:, [column for column in design_columns if column in samples.columns]].copy()
            out_dir.mkdir(parents=True, exist_ok=True)
            table.to_csv(outputs["table"], sep="\t", index=False, float_format="%.6g")
            plot = render_association_plot(table, outputs["plot"])
            factorial.to_csv(outputs["factorial"], sep="\t", index=False, float_format="%.6g")
            coefficients.to_csv(outputs["coefficients"], sep="\t", index=False, float_format="%.6g")
            group_summary.to_csv(outputs["group_summary"], sep="\t", index=False, float_format="%.6g")
            design.to_csv(outputs["design"], sep="\t", index=False)
            factorial_plot = render_factorial_association_plot(factorial, outputs["factorial_plot"])
            generated_plots = [path for path in (plot, factorial_plot) if path]
            write_json(outputs["parameters"], {"fingerprint": stage_fp, **association_settings, "generated_plots": [str(path) for path in generated_plots]})
            written = [
                outputs["table"],
                outputs["factorial"],
                outputs["coefficients"],
                outputs["group_summary"],
                outputs["design"],
                outputs["parameters"],
                *generated_plots,
            ]
            self._mark(stage, "finished", f"Finished {network_name} module-trait association.", stage_fp, written)
            return outputs
        except Exception as exc:
            self._mark(stage, "failed", str(exc), stage_fp)
            raise

    def consensus_modules(self, force: bool = False) -> dict[str, Path]:
        stage = "consensus_modules"
        matrix_path = self.output / "04_network" / "observed" / "selected_event_matrix.tsv"
        observed_path = self.output / "04_network" / "observed" / "event_module_assignments.tsv"
        samples_path = self.output / "02_matrix" / "sample_metadata_resolved.tsv"
        required = [matrix_path, observed_path, samples_path]
        if not all(path.exists() for path in required):
            raise FileNotFoundError(
                "Missing input: observed network cache is incomplete. Run network_observed first; "
                "consensus_modules will not run it automatically."
            )
        out_dir = self.output / "04c_consensus_modules"
        outputs = {
            "assignments": out_dir / "consensus_event_module_assignments.tsv",
            "modules": out_dir / "consensus_module_summary.tsv",
            "scores": out_dir / "consensus_module_program_scores.tsv",
            "bootstrap": out_dir / "coassignment_bootstrap_runs.tsv",
            "edges": out_dir / "bootstrap_coassignment_edges.tsv",
            "matrix": out_dir / "bootstrap_coassignment_matrix.npz",
            "diagnostics": out_dir / "consensus_diagnostics.tsv",
            "parameter_search": out_dir / "consensus_parameter_search.tsv",
            "parameters": out_dir / "parameters.json",
        }
        settings = {
            "network_partition": self.config.network.__dict__,
            "consensus": self.config.consensus.__dict__,
            "semantics": "outer_sample_stratified_bootstrap_coassignment",
        }
        stage_fp = self._fingerprint(stage, required, settings)
        state = self._state().get("stages", {}).get(stage, {})
        if not force and state.get("fingerprint") == stage_fp and all(path.exists() for path in outputs.values()):
            cached = [*outputs.values(), *sorted((out_dir / "plots").glob("*.png"))]
            self._mark(stage, "cache_hit", "Loaded existing bootstrap co-assignment consensus cache.", stage_fp, cached)
            return outputs
        self._mark(stage, "running", "Building bootstrap event co-assignment consensus modules.", stage_fp)
        try:
            matrix = pd.read_csv(matrix_path, sep="\t", index_col=0)
            observed = pd.read_csv(observed_path, sep="\t", low_memory=False)
            samples = pd.read_csv(samples_path, sep="\t", low_memory=False)
            result = build_bootstrap_coassignment_consensus(
                matrix,
                observed,
                samples,
                self.config.network,
                self.config.consensus,
                progress=self._emit,
            )
            out_dir.mkdir(parents=True, exist_ok=True)
            result.assignments.to_csv(outputs["assignments"], sep="\t", index=False, float_format="%.6g")
            result.module_summary.to_csv(outputs["modules"], sep="\t", index=False, float_format="%.6g")
            result.program_scores.to_csv(outputs["scores"], sep="\t", float_format="%.6g")
            result.bootstrap_runs.to_csv(outputs["bootstrap"], sep="\t", index=False, float_format="%.6g")
            result.coassignment_edges.to_csv(outputs["edges"], sep="\t", index=False, float_format="%.6g")
            result.diagnostics.to_csv(outputs["diagnostics"], sep="\t", index=False, float_format="%.6g")
            result.parameter_search.to_csv(
                outputs["parameter_search"], sep="\t", index=False, float_format="%.6g"
            )
            save_coassignment_matrix(result, outputs["matrix"])
            plots = render_consensus_plots(result, out_dir / "plots")
            write_json(
                outputs["parameters"],
                {
                    "fingerprint": stage_fp,
                    **settings,
                    "generated_plots": [str(path) for path in plots],
                    "selected_partition_parameters": result.parameter_search.loc[
                        result.parameter_search["selected"]
                    ].iloc[0].to_dict(),
                    "matrix_format": {
                        "path": str(outputs["matrix"]),
                        "arrays": ["event_keys", "coassignment"],
                    },
                    "interpretation_limit": (
                        "Consensus confidence is internal resampling evidence from the discovery cohort. "
                        "It improves module reproducibility but is not external preservation."
                    ),
                },
            )
            modules = int(result.assignments["consensus_module"].loc[
                result.assignments["consensus_module"].ne("grey")
            ].nunique())
            written = [*outputs.values(), *plots]
            self._mark(
                stage,
                "finished",
                f"Bootstrap co-assignment consensus finished with {modules} retained modules.",
                stage_fp,
                written,
            )
            return outputs
        except Exception as exc:
            self._mark(stage, "failed", str(exc), stage_fp)
            raise

    def hybrid_programs(self, force: bool = False) -> dict[str, Path]:
        """Retain tight sparse and diffuse TOM programs without rerunning upstream stages."""
        stage = "hybrid_programs"
        qc_dir = self.output / "03_stable_as_qc"
        samples_path = self.output / "02_matrix" / "sample_metadata_resolved.tsv"
        qc_path = qc_dir / "stable_as_event_qc.tsv"
        source = self.config.hybrid_programs.matrix_source
        residual_dir = self.output / "04_network" / "residual"
        residual_ready = all(
            (residual_dir / name).exists()
            for name in ("selected_event_matrix.tsv", "event_module_assignments.tsv")
        )
        if source == "residual" and not residual_ready:
            raise FileNotFoundError(
                "Missing input: residual tight network cache is incomplete. Run network_residual "
                "first; hybrid_programs will not run it automatically."
            )
        resolved_source = (
            "residual"
            if source == "residual" or (source == "auto" and residual_ready)
            else "observed"
        )
        network_dir = self.output / "04_network" / resolved_source
        matrix_path = network_dir / "selected_event_matrix.tsv"
        tight_path = network_dir / "event_module_assignments.tsv"
        required = [matrix_path, tight_path, qc_path, samples_path]
        if not all(path.exists() for path in required):
            raise FileNotFoundError(
                "Missing input: stable-AS QC or tight network cache is incomplete. Run stable_qc "
                f"and network_{resolved_source} first; hybrid_programs will not run them automatically."
            )
        coassignment_path = (
            self.output
            / "04c_consensus_modules"
            / "bootstrap_coassignment_matrix.npz"
        )
        optional_inputs = [coassignment_path] if coassignment_path.exists() else []
        cohort_inputs = [
            path
            for cohort in self.config.external_cohorts
            for path in (cohort.psi_matrix, cohort.sample_metadata)
        ]
        out_dir = self.output / "04d_hybrid_programs"
        outputs = {
            "events": out_dir / "hybrid_event_program_assignments.tsv",
            "programs": out_dir / "hybrid_program_summary.tsv",
            "scores": out_dir / "hybrid_program_scores.tsv",
            "bootstrap": out_dir / "diffuse_bootstrap_coassignment_runs.tsv",
            "overlap": out_dir / "tight_diffuse_program_overlap.tsv",
            "engine": out_dir / "diffuse_engine_audit.tsv",
            "preservation": out_dir / "external_program_preservation.tsv",
            "mapping": out_dir / "external_event_mapping.tsv",
            "null": out_dir / "external_preservation_null.tsv",
            "preservation_audit": out_dir / "external_preservation_audit.tsv",
            "parameters": out_dir / "parameters.json",
        }
        settings = {
            "matrix_source_requested": source,
            "matrix_source_resolved": resolved_source,
            "hybrid_programs": self.config.hybrid_programs.__dict__,
            "network": self.config.network.__dict__,
            "benchmark": self.config.benchmark.__dict__,
            "preservation": self.config.preservation.__dict__,
            "external_cohorts": [
                {
                    "name": cohort.name,
                    "psi_matrix": str(cohort.psi_matrix),
                    "sample_metadata": str(cohort.sample_metadata),
                }
                for cohort in self.config.external_cohorts
            ],
        }
        stage_fp = self._fingerprint(
            stage,
            [*required, *optional_inputs, *cohort_inputs],
            settings,
        )
        state = self._state().get("stages", {}).get(stage, {})
        if (
            not force
            and state.get("fingerprint") == stage_fp
            and all(path.exists() for path in outputs.values())
        ):
            cached = [*outputs.values(), *sorted((out_dir / "plots").glob("*.png"))]
            self._mark(
                stage,
                "cache_hit",
                "Loaded existing tight/diffuse program cache.",
                stage_fp,
                cached,
            )
            return outputs
        self._mark(
            stage,
            "running",
            f"Discovering tight and diffuse programs from the cached {resolved_source} matrix.",
            stage_fp,
        )
        try:
            matrix = pd.read_csv(matrix_path, sep="\t", index_col=0)
            qc = pd.read_csv(qc_path, sep="\t", index_col=0)
            samples = pd.read_csv(samples_path, sep="\t")
            tight = pd.read_csv(tight_path, sep="\t", low_memory=False)
            tight_coassignment: pd.Series | None = None
            if coassignment_path.exists():
                stored = np.load(coassignment_path, allow_pickle=True)
                coassignment_events = pd.Index(
                    stored["event_keys"].astype(str),
                    dtype=object,
                )
                coassignment = stored["coassignment"]
                positions = pd.Series(
                    np.arange(len(coassignment_events)),
                    index=coassignment_events,
                )
                tight_indexed = tight.set_index("event_key").reindex(matrix.index)
                tight_coassignment = pd.Series(np.nan, index=matrix.index, dtype=float)
                for _module, frame in tight_indexed.loc[
                    tight_indexed["module"].fillna("grey").ne("grey")
                ].groupby("module", observed=True):
                    available = frame.index.intersection(coassignment_events)
                    if len(available) < 2:
                        continue
                    indices = positions.reindex(available).to_numpy(dtype=int)
                    within = coassignment[np.ix_(indices, indices)].astype(float)
                    np.fill_diagonal(within, np.nan)
                    tight_coassignment.loc[available] = np.nanmean(within, axis=1)
            result = build_hybrid_programs(
                matrix,
                qc,
                samples,
                tight,
                self.config.network,
                self.config.benchmark,
                self.config.hybrid_programs,
                tight_coassignment=tight_coassignment,
                progress=self._emit,
            )

            preservation_results: list[PreservationResult] = []
            reference_resolved = matrix_path.resolve()
            reused_reference = [
                cohort.name
                for cohort in self.config.external_cohorts
                if cohort.psi_matrix.resolve() == reference_resolved
            ]
            if reused_reference:
                raise ValueError(
                    "External program preservation requires an independent target PSI "
                    "matrix. These cohorts reuse the discovery matrix: "
                    + ", ".join(reused_reference)
                )
            for program_type, column in (
                ("tight", "tight_program"),
                ("diffuse", "diffuse_program"),
            ):
                assignments = (
                    result.event_assignments.loc[
                        result.event_assignments[column].ne("grey"),
                        ["event_key", column],
                    ]
                    .rename(columns={column: "consensus_module"})
                    .copy()
                )
                preserved = evaluate_external_module_preservation(
                    matrix,
                    assignments,
                    self.config.external_cohorts,
                    self.config.preservation,
                    progress=self._emit,
                )
                for table in (
                    preserved.module_preservation,
                    preserved.event_mapping,
                    preserved.null_summary,
                    preserved.target_audit,
                ):
                    table.insert(0, "program_type", program_type)
                preservation_results.append(preserved)

            preservation = PreservationResult(
                pd.concat(
                    [item.module_preservation for item in preservation_results],
                    ignore_index=True,
                ),
                pd.concat(
                    [item.event_mapping for item in preservation_results],
                    ignore_index=True,
                ),
                pd.concat(
                    [item.null_summary for item in preservation_results],
                    ignore_index=True,
                ),
                pd.concat(
                    [item.target_audit for item in preservation_results],
                    ignore_index=True,
                ),
            )
            result.program_summary = add_external_preservation_confidence(
                result.program_summary,
                preservation.module_preservation,
                self.config.hybrid_programs,
            )
            out_dir.mkdir(parents=True, exist_ok=True)
            result.event_assignments.to_csv(
                outputs["events"], sep="\t", index=False, float_format="%.6g"
            )
            result.program_summary.to_csv(
                outputs["programs"], sep="\t", index=False, float_format="%.6g"
            )
            result.program_scores.to_csv(
                outputs["scores"], sep="\t", float_format="%.6g"
            )
            result.bootstrap_runs.to_csv(
                outputs["bootstrap"], sep="\t", index=False, float_format="%.6g"
            )
            result.overlap.to_csv(
                outputs["overlap"], sep="\t", index=False, float_format="%.6g"
            )
            result.engine_audit.to_csv(outputs["engine"], sep="\t", index=False)
            preservation.module_preservation.to_csv(
                outputs["preservation"], sep="\t", index=False, float_format="%.6g"
            )
            preservation.event_mapping.to_csv(outputs["mapping"], sep="\t", index=False)
            preservation.null_summary.to_csv(
                outputs["null"], sep="\t", index=False, float_format="%.6g"
            )
            preservation.target_audit.to_csv(
                outputs["preservation_audit"], sep="\t", index=False
            )
            plots = render_hybrid_program_plots(result, out_dir / "plots")
            plots.extend(render_preservation_plots(preservation, out_dir / "plots"))
            write_json(
                outputs["parameters"],
                {
                    "fingerprint": stage_fp,
                    **settings,
                    "generated_plots": [str(path) for path in plots],
                    "program_types": {
                        "tight": (
                            "Sparse CoSplice modules retained by the configured network "
                            "and membership filters."
                        ),
                        "diffuse": (
                            "WGCNA-style signed soft-adjacency/TOM modules that may retain "
                            "weak but coherent event correlations."
                        ),
                    },
                    "internal_confidence": (
                        "0.45 * bootstrap support/co-assignment + 0.35 * absolute program "
                        "membership + 0.20 * replicate repeatability"
                    ),
                    "overall_confidence": (
                        "Internal confidence when no independent cohort is configured; "
                        "otherwise 0.75 * internal + 0.25 * external preservation class score."
                    ),
                    "external_validation_limit": (
                        "External confidence is independent only when configured external "
                        "cohorts are genuinely independent and contain exact matching event keys."
                    ),
                },
            )
            written = [*outputs.values(), *plots]
            n_tight = int(
                result.program_summary["program_type"].eq("tight").sum()
            )
            n_diffuse = int(
                result.program_summary["program_type"].eq("diffuse").sum()
            )
            self._mark(
                stage,
                "finished",
                f"Hybrid programs finished: {n_tight} tight and {n_diffuse} diffuse.",
                stage_fp,
                written,
            )
            return outputs
        except Exception as exc:
            self._mark(stage, "failed", str(exc), stage_fp)
            raise

    def preserve_modules(self, force: bool = False) -> dict[str, Path]:
        stage = "preserve_modules"
        reference_path = self.output / "04_network" / "observed" / "selected_event_matrix.tsv"
        consensus_path = self.output / "04c_consensus_modules" / "consensus_event_module_assignments.tsv"
        required = [reference_path, consensus_path]
        if not all(path.exists() for path in required):
            raise FileNotFoundError(
                "Missing input: consensus module cache is incomplete. Run consensus_modules first; "
                "preserve_modules will not run it automatically."
            )
        cohort_inputs = [
            path
            for cohort in self.config.external_cohorts
            for path in (cohort.psi_matrix, cohort.sample_metadata)
        ]
        reference_resolved = reference_path.resolve()
        reused_reference = [
            cohort.name
            for cohort in self.config.external_cohorts
            if cohort.psi_matrix.resolve() == reference_resolved
        ]
        if reused_reference:
            raise ValueError(
                "External preservation requires an independent target PSI matrix. "
                "These cohorts reuse the discovery matrix: " + ", ".join(reused_reference)
            )
        out_dir = self.output / "06d_external_preservation"
        outputs = {
            "modules": out_dir / "external_module_preservation.tsv",
            "mapping": out_dir / "external_event_mapping.tsv",
            "null": out_dir / "external_preservation_null_summary.tsv",
            "audit": out_dir / "external_preservation_target_audit.tsv",
            "parameters": out_dir / "parameters.json",
        }
        settings = {
            "external_cohorts": [
                {
                    "name": cohort.name,
                    "psi_matrix": str(cohort.psi_matrix),
                    "sample_metadata": str(cohort.sample_metadata),
                }
                for cohort in self.config.external_cohorts
            ],
            "preservation": self.config.preservation.__dict__,
            "reference_modules": "bootstrap_coassignment_consensus",
            "event_mapping": "exact_event_key",
        }
        stage_fp = self._fingerprint(stage, [*required, *cohort_inputs], settings)
        state = self._state().get("stages", {}).get(stage, {})
        if not force and state.get("fingerprint") == stage_fp and all(path.exists() for path in outputs.values()):
            cached = [*outputs.values(), *sorted((out_dir / "plots").glob("*.png"))]
            self._mark(stage, "cache_hit", "Loaded existing external module-preservation cache.", stage_fp, cached)
            return outputs
        self._mark(stage, "running", "Projecting consensus modules into configured external cohorts.", stage_fp)
        try:
            reference = pd.read_csv(reference_path, sep="\t", index_col=0)
            consensus = pd.read_csv(consensus_path, sep="\t", low_memory=False)
            result = evaluate_external_module_preservation(
                reference,
                consensus,
                self.config.external_cohorts,
                self.config.preservation,
                progress=self._emit,
            )
            out_dir.mkdir(parents=True, exist_ok=True)
            result.module_preservation.to_csv(outputs["modules"], sep="\t", index=False, float_format="%.6g")
            result.event_mapping.to_csv(outputs["mapping"], sep="\t", index=False)
            result.null_summary.to_csv(outputs["null"], sep="\t", index=False, float_format="%.6g")
            result.target_audit.to_csv(outputs["audit"], sep="\t", index=False)
            plots = render_preservation_plots(result, out_dir / "plots")
            write_json(
                outputs["parameters"],
                {
                    "fingerprint": stage_fp,
                    **settings,
                    "generated_plots": [str(path) for path in plots],
                    "classification": {
                        "preserved": f"Zsummary_wgcna_like >= {self.config.preservation.preserved_zsummary:g}",
                        "partially_preserved": (
                            f"{self.config.preservation.partial_zsummary:g} <= Zsummary_wgcna_like < "
                            f"{self.config.preservation.preserved_zsummary:g}"
                        ),
                        "not_preserved": (
                            f"Zsummary_wgcna_like < {self.config.preservation.partial_zsummary:g}"
                        ),
                    },
                    "interpretation_limit": (
                        "Zsummary_wgcna_like is Katana's permutation-calibrated density/connectivity composite, "
                        "not the exact WGCNA modulePreservation implementation. External cohort evidence is "
                        "independent only when the configured target cohort is independent."
                    ),
                },
            )
            written = [*outputs.values(), *plots]
            completed = int(result.target_audit["status"].eq("finished").sum())
            if self.config.external_cohorts:
                message = (
                    f"External preservation finished for {completed}/{len(self.config.external_cohorts)} cohorts."
                )
            else:
                message = "No external cohorts configured; wrote empty auditable preservation outputs."
            self._mark(
                stage,
                "finished",
                message,
                stage_fp,
                written,
            )
            return outputs
        except Exception as exc:
            self._mark(stage, "failed", str(exc), stage_fp)
            raise

    def validate_modules(self, force: bool = False) -> dict[str, Path]:
        stage = "validate_modules"
        matrix_path = self.output / "04_network" / "observed" / "selected_event_matrix.tsv"
        observed_path = self.output / "04_network" / "observed" / "event_module_assignments.tsv"
        samples_path = self.output / "02_matrix" / "sample_metadata_resolved.tsv"
        required = [matrix_path, observed_path, samples_path]
        if not all(path.exists() for path in required):
            raise FileNotFoundError("Missing input: observed network not fit. Run 4A. Observed network first.")
        residual_path = self.output / "04_network" / "residual" / "event_module_assignments.tsv"
        inputs = [*required, *([residual_path] if residual_path.exists() else [])]
        out_dir = self.output / "06_module_validation"
        outputs = {
            "events": out_dir / "event_membership_stability.tsv",
            "modules": out_dir / "module_bootstrap_stability.tsv",
            "bootstrap": out_dir / "bootstrap_partition_stability.tsv",
            "parameters_table": out_dir / "network_parameter_sensitivity.tsv",
            "preservation_pairs": out_dir / "observed_residual_preservation_pairs.tsv",
            "preservation_summary": out_dir / "observed_residual_preservation_summary.tsv",
            "report": out_dir / "module_validation_report.html",
            "parameters": out_dir / "parameters.json",
        }
        settings = {
            "network": self.config.network.__dict__,
            "validation": self.config.validation.__dict__,
            "residual_network_available": residual_path.exists(),
        }
        stage_fp = self._fingerprint(stage, inputs, settings)
        state = self._state().get("stages", {}).get(stage, {})
        if not force and state.get("fingerprint") == stage_fp and all(path.exists() for path in outputs.values()):
            self._mark(stage, "cache_hit", "Loaded existing module-validation cache.", stage_fp, list(outputs.values()))
            return outputs
        self._mark(stage, "running", "Validating module bootstrap stability, parameter sensitivity and preservation.", stage_fp)
        try:
            matrix = pd.read_csv(matrix_path, sep="\t", index_col=0)
            observed = pd.read_csv(observed_path, sep="\t")
            samples = pd.read_csv(samples_path, sep="\t")
            residual = pd.read_csv(residual_path, sep="\t") if residual_path.exists() else None
            result = validate_network_modules(
                matrix,
                observed,
                samples,
                self.config.network,
                self.config.validation,
                residual_assignments=residual,
                progress=self._emit,
            )
            out_dir.mkdir(parents=True, exist_ok=True)
            result.event_stability.to_csv(outputs["events"], sep="\t", index=False, float_format="%.6g")
            result.module_stability.to_csv(outputs["modules"], sep="\t", index=False, float_format="%.6g")
            result.bootstrap_runs.to_csv(outputs["bootstrap"], sep="\t", index=False, float_format="%.6g")
            result.parameter_sensitivity.to_csv(outputs["parameters_table"], sep="\t", index=False, float_format="%.6g")
            result.preservation_pairs.to_csv(outputs["preservation_pairs"], sep="\t", index=False, float_format="%.6g")
            result.preservation_summary.to_csv(outputs["preservation_summary"], sep="\t", index=False, float_format="%.6g")
            plots = render_validation_plots(result, out_dir / "plots", self.config.validation)
            render_validation_report(result, outputs["report"])
            write_json(
                outputs["parameters"],
                {
                    "fingerprint": stage_fp,
                    **settings,
                    "generated_plots": [str(path) for path in plots],
                    "classification": {
                        "robust": f"validation_score >= {self.config.validation.robust_jaccard:g}",
                        "moderate": f"{self.config.validation.moderate_jaccard:g} <= validation_score < {self.config.validation.robust_jaccard:g}",
                        "unstable": f"validation_score < {self.config.validation.moderate_jaccard:g}",
                    },
                },
            )
            robust = int(result.module_stability["stability_class"].eq("robust").sum())
            moderate = int(result.module_stability["stability_class"].eq("moderate").sum())
            unstable = int(result.module_stability["stability_class"].eq("unstable").sum())
            written = [*outputs.values(), *plots]
            self._mark(
                stage,
                "finished",
                f"Module validation finished: {robust} robust, {moderate} moderate, {unstable} unstable.",
                stage_fp,
                written,
            )
            return outputs
        except Exception as exc:
            self._mark(stage, "failed", str(exc), stage_fp)
            raise

    def report(self) -> Path:
        stage = "report"
        report_path = self.output / "07_report" / "katana_cosplice_report.html"
        pdf_path = self.output / "07_report" / "katana_cosplice_validation_report.pdf"
        browser_path = self.output / "07_report" / "module_browser.html"
        self._mark(stage, "running", "Building cached-output report.")
        state = self._state()
        rows: list[str] = []
        for name, details in state.get("stages", {}).items():
            rows.append(
                f"<tr><td>{escape(name)}</td><td>{escape(str(details.get('status', '')))}</td>"
                f"<td>{escape(str(details.get('message', '')))}</td><td>{escape(str(details.get('updated_at', '')))}</td></tr>"
            )
        summary_blocks: list[str] = []
        qc_path = self.output / "03_stable_as_qc" / "parameters.json"
        if qc_path.exists():
            import json
            qc = json.loads(qc_path.read_text(encoding="utf-8"))
            summary_blocks.append(
                f"<section><h2>Stable-AS QC</h2><p><strong>{int(qc.get('n_stable_events', 0)):,}</strong> events passed "
                "the configured condition-coverage, junction-read, PSI dynamic-range and replicate-repeatability filters.</p></section>"
            )
        for network_name in ("observed", "residual"):
            diagnostics_path = self.output / "04_network" / network_name / "network_diagnostics.tsv"
            sensitivity_path = self.output / "04_network" / network_name / "module_resolution_sensitivity.tsv"
            if not diagnostics_path.exists():
                continue
            diagnostics = pd.read_csv(diagnostics_path, sep="\t")
            non_grey = diagnostics.loc[diagnostics["module"].ne("grey")]
            largest_fraction = float(diagnostics["largest_module_fraction"].iloc[0]) if "largest_module_fraction" in diagnostics else float("nan")
            warning = (
                "<p class='warning'>Dominant-module warning: the largest module contains at least 60% of selected events. "
                "Review the module-resolution sensitivity table and external validation before treating submodules as definitive.</p>"
                if largest_fraction >= 0.60
                else ""
            )
            sensitivity_link = f"<a href='../04_network/{network_name}/{sensitivity_path.name}'>module-resolution sensitivity TSV</a>" if sensitivity_path.exists() else ""
            summary_blocks.append(
                f"<section><h2>{network_name.title()} network</h2><p>{len(non_grey)} retained modules; "
                f"largest-module fraction {largest_fraction:.1%}. {sensitivity_link}</p>{warning}</section>"
            )
            factorial_path = self.output / "05_association" / network_name / "factorial_module_association.tsv"
            if factorial_path.exists():
                factorial = pd.read_csv(factorial_path, sep="\t")
                significant = factorial.loc[pd.to_numeric(factorial.get("FDR"), errors="coerce").le(0.05)]
                summary_blocks.append(
                    f"<section><h2>{network_name.title()} factorial association</h2>"
                    f"<p>{len(significant)} module-term associations pass FDR 0.05. Model factors: "
                    f"{escape(', '.join(self.config.factorial_factors) or 'not configured')}; batch covariate: "
                    f"{escape(self.config.batch_column or 'not configured')}.</p>"
                    "<p>These are within-cohort program associations, not independent validation or causal effects.</p></section>"
                )
        consensus_path = self.output / "04c_consensus_modules" / "consensus_module_summary.tsv"
        consensus_diagnostics_path = self.output / "04c_consensus_modules" / "consensus_diagnostics.tsv"
        if consensus_path.exists():
            consensus = pd.read_csv(consensus_path, sep="\t")
            retained = (
                consensus.loc[consensus["consensus_module"].ne("grey")]
                if "consensus_module" in consensus
                else consensus
            )
            diagnostics = (
                pd.read_csv(consensus_diagnostics_path, sep="\t")
                if consensus_diagnostics_path.exists()
                else pd.DataFrame()
            )
            observed_ari = "not available"
            if not diagnostics.empty and "adjusted_rand_index_observed_vs_consensus" in diagnostics:
                observed_ari = (
                    f"{float(diagnostics.iloc[0]['adjusted_rand_index_observed_vs_consensus']):.3f}"
                )
            summary_blocks.append(
                "<section><h2>Bootstrap co-assignment consensus</h2>"
                f"<p>{len(retained)} retained consensus modules; observed-versus-consensus ARI: "
                f"{escape(observed_ari)}.</p>"
                "<p>Modules are rebuilt from event-event bootstrap co-assignment probabilities. "
                "This measures internal assignment reproducibility; it is not external validation.</p></section>"
            )
        hybrid_path = (
            self.output
            / "04d_hybrid_programs"
            / "hybrid_program_summary.tsv"
        )
        if hybrid_path.exists():
            hybrid = pd.read_csv(hybrid_path, sep="\t")
            type_counts = hybrid["program_type"].value_counts()
            confidence_counts = hybrid["confidence_class"].value_counts()
            externally_evaluated = int(
                hybrid["external_validation_status"].eq("evaluated").sum()
            )
            summary_blocks.append(
                "<section><h2>Tight and diffuse co-splicing programs</h2>"
                f"<p>{int(type_counts.get('tight', 0))} tight sparse programs and "
                f"{int(type_counts.get('diffuse', 0))} diffuse soft-TOM programs were retained; "
                f"{int(confidence_counts.get('robust', 0))} robust, "
                f"{int(confidence_counts.get('moderate', 0))} moderate and "
                f"{int(confidence_counts.get('exploratory', 0))} exploratory confidence labels.</p>"
                f"<p>{externally_evaluated} program rows include independent-cohort preservation. "
                "All remaining confidence labels are internal-only and combine bootstrap "
                "co-assignment, program membership and replicate repeatability.</p></section>"
            )
        validation_path = self.output / "06_module_validation" / "module_bootstrap_stability.tsv"
        if validation_path.exists():
            validation = pd.read_csv(validation_path, sep="\t")
            counts = validation["stability_class"].value_counts()
            summary_blocks.append(
                "<section><h2>Module validation</h2>"
                f"<p>{int(counts.get('robust', 0))} robust, {int(counts.get('moderate', 0))} moderate and "
                f"{int(counts.get('unstable', 0))} unstable modules under bootstrap recovery. Parameter sensitivity is reported separately.</p>"
                "<p>This is internal robustness assessment, not independent external validation.</p></section>"
            )
        preservation_path = self.output / "06d_external_preservation" / "external_module_preservation.tsv"
        preservation_audit_path = self.output / "06d_external_preservation" / "external_preservation_target_audit.tsv"
        if preservation_path.exists():
            preservation = pd.read_csv(preservation_path, sep="\t")
            audit = (
                pd.read_csv(preservation_audit_path, sep="\t")
                if preservation_audit_path.exists()
                else pd.DataFrame()
            )
            classes = (
                preservation["preservation_class"].value_counts()
                if "preservation_class" in preservation
                else pd.Series(dtype="int64")
            )
            completed_targets = (
                int(audit["status"].eq("finished").sum()) if "status" in audit else 0
            )
            summary_blocks.append(
                "<section><h2>External module preservation</h2>"
                f"<p>{completed_targets} external cohorts completed; "
                f"{int(classes.get('preserved', 0))} preserved, "
                f"{int(classes.get('partially_preserved', 0))} partially preserved and "
                f"{int(classes.get('not_preserved', 0))} not preserved module-target results.</p>"
                "<p>Zsummary_wgcna_like is Katana's permutation-calibrated density/connectivity composite, "
                "not the exact WGCNA modulePreservation statistic. Independence requires a genuinely separate cohort.</p></section>"
            )
        annotation_path = self.output / "06b_module_annotation" / "module_evidence_summary.tsv"
        if annotation_path.exists():
            annotation = pd.read_csv(annotation_path, sep="\t")
            rbp_path = self.output / "06b_module_annotation" / "module_rbp_enrichment.tsv"
            rbp = pd.read_csv(rbp_path, sep="\t") if rbp_path.exists() else pd.DataFrame()
            rbp_significant = int(pd.to_numeric(rbp.get("FDR"), errors="coerce").le(0.05).sum()) if not rbp.empty else 0
            summary_blocks.append(
                "<section><h2>Module evidence annotation</h2>"
                f"<p>{len(annotation)} modules have cached DEG/DTU/RBP support summaries; "
                f"{rbp_significant} module–RBP-set enrichments pass FDR 0.05.</p>"
                "<p>These layers prioritize mechanisms but do not define modules or prove causal regulation.</p></section>"
            )
        benchmark_path = self.output / "06c_method_benchmark" / "method_benchmark_summary.tsv"
        if benchmark_path.exists():
            benchmark = pd.read_csv(benchmark_path, sep="\t")
            finished = benchmark.loc[benchmark["status"].eq("finished")]
            best = ""
            if not finished.empty and "mean_bootstrap_ARI" in finished:
                ranked = finished.dropna(subset=["mean_bootstrap_ARI"]).sort_values("mean_bootstrap_ARI", ascending=False)
                if not ranked.empty:
                    best = f" Highest internal bootstrap ARI: {escape(str(ranked.iloc[0]['method']))} ({float(ranked.iloc[0]['mean_bootstrap_ARI']):.3f})."
            summary_blocks.append(
                "<section><h2>Method benchmark</h2>"
                f"<p>{len(finished)} methods completed on one shared stable-AS event universe.{best}</p>"
                "<p>This is an internal same-cohort comparison and does not establish external superiority.</p></section>"
            )
        interpretation_path = self.output / "06e_module_interpretation" / "module_regulator_prioritization.tsv"
        interpretation_audit_path = self.output / "06e_module_interpretation" / "interpretation_source_audit.tsv"
        if interpretation_path.exists():
            regulators = pd.read_csv(interpretation_path, sep="\t")
            audit = pd.read_csv(interpretation_audit_path, sep="\t") if interpretation_audit_path.exists() else pd.DataFrame()
            loaded_sources = int(audit["status"].eq("loaded").sum()) if "status" in audit else 0
            summary_blocks.append(
                "<section><h2>Module functional interpretation</h2>"
                f"<p>{len(regulators)} module-RBP regulator candidates were prioritized from {loaded_sources} loaded support sources.</p>"
                "<p>Hub, motif, RBP-expression, DEG/DTU and GO/KEGG evidence prioritize candidate mechanisms; they do not establish causal RBP regulation.</p></section>"
            )
        motif_path = self.output / "06e_module_interpretation" / "event_rbp_motif_hits.tsv"
        isoform_path = self.output / "06e_module_interpretation" / "dominant_isoform_switches.tsv"
        isoform_enrichment_path = self.output / "06e_module_interpretation" / "module_isoform_switch_enrichment.tsv"
        motif = read_cached_table(motif_path)
        isoform = read_cached_table(isoform_path)
        isoform_enrichment = read_cached_table(isoform_enrichment_path)
        if not motif.empty or not isoform.empty:
            summary_blocks.append(
                "<section><h2>Sequence and isoform-switch follow-up</h2>"
                f"<p>{len(motif):,} event-motif hit rows and {len(isoform):,} descriptive dominant-isoform switches are cached; "
                f"{len(isoform_enrichment):,} module switch-enrichment rows were tested.</p>"
                "<p>Motifs are strand-oriented sequence matches. Isoform switches require dominant transcript change plus configured delta-IF/TPM filters; inferential support is inherited from DTU. Neither layer is causal evidence.</p></section>"
            )
        statistics_path = self.output / "06f_module_statistics" / "event_beta_binomial_condition_effect.tsv"
        benchmark_delta_path = self.output / "06f_module_statistics" / "beta_binomial_vs_simple_delta_dpsi.tsv"
        if statistics_path.exists():
            statistics = pd.read_csv(statistics_path, sep="\t")
            benchmark_delta = pd.read_csv(benchmark_delta_path, sep="\t") if benchmark_delta_path.exists() else pd.DataFrame()
            significant = int(pd.to_numeric(statistics.get("condition_effect_FDR"), errors="coerce").le(0.05).sum()) if not statistics.empty else 0
            overlap = float(benchmark_delta.iloc[0]["jaccard"]) if not benchmark_delta.empty else float("nan")
            summary_blocks.append(
                "<section><h2>Replicate-aware statistics</h2>"
                f"<p>{len(statistics)} beta-binomial event models completed; {significant} pass global FDR 0.05. "
                f"Top-list agreement with simple delta PSI (first configured top-N) is Jaccard {overlap:.3f}.</p>"
                "<p>This read-count likelihood test complements rather than replaces rMATS. Batch is included only when explicitly configured.</p></section>"
            )
        meta_path = self.output / "06g_cross_cohort_meta" / "random_effects_module_meta_analysis.tsv"
        meta_audit_path = self.output / "06g_cross_cohort_meta" / "meta_analysis_cohort_audit.tsv"
        consensus_meta_path = self.output / "06g_cross_cohort_meta" / "multi_cohort_module_consensus.tsv"
        if meta_path.exists():
            meta = pd.read_csv(meta_path, sep="\t")
            audit = pd.read_csv(meta_audit_path, sep="\t") if meta_audit_path.exists() else pd.DataFrame()
            completed = int(audit["status"].isin(["finished", "projected_effects_only"]).sum()) if "status" in audit else 0
            summary_blocks.append(
                "<section><h2>Cross-cohort module meta-analysis</h2>"
                f"<p>{completed} external cohorts contributed valid projected or independent-module results; {len(meta)} random-effects module contrasts were estimable.</p>"
                "<p>Events are aligned by exact canonical identity, same-build exact splice-site identity, or an explicit validated mapping. Gained/lost/rewired labels require independently discovered target modules and are never inferred from projection alone.</p></section>"
            )
        if consensus_meta_path.exists():
            consensus_meta = pd.read_csv(consensus_meta_path, sep="\t")
            consensus_states = consensus_meta.get("multi_cohort_state", pd.Series(dtype="object")).value_counts()
            summary_blocks.append(
                "<section><h2>Multi-cohort module consensus</h2>"
                f"<p>{int(consensus_states.get('preserved_across_cohorts', 0))} preserved across configured cohorts, "
                f"{int(consensus_states.get('partially_preserved', 0))} partially preserved and "
                f"{int(consensus_states.get('not_preserved', 0))} not preserved.</p>"
                "<p>Consensus states require FDR-supported event-level module overlap; strict rewiring additionally requires a negative permutation-supported connectivity change.</p></section>"
            )
        factor_path = self.output / "06h_factor_contrasts" / "factor_evidence_summary.tsv"
        if factor_path.exists():
            factor = pd.read_csv(factor_path, sep="\t")
            support_counts = factor.get("interpretation", pd.Series(dtype="object")).value_counts()
            summary_blocks.append(
                "<section><h2>Configured factor contrasts and phenocopy</h2>"
                f"<p>{len(factor)} predeclared context-perturbation summaries were tested; "
                f"{int(support_counts.get('partial concordant phenocopy support', 0))} met the configured partial-concordance rule.</p>"
                "<p>Results combine broad stable-AS signature correlation, significant-event overlap, direction agreement and "
                "difference-in-differences program tests. They prioritize phenocopy/occlusion follow-up and do not prove LOF, GOF or causality.</p></section>"
            )
        html = f"""<!doctype html>
<html><head><meta charset="utf-8"><title>{escape(self.config.project_name)} CoSplice report</title>
<style>body{{font-family:Arial,sans-serif;margin:36px;color:#17212b;max-width:1200px}}table{{border-collapse:collapse;width:100%}}td,th{{border:1px solid #cbd5df;padding:8px;text-align:left}}th{{background:#e9f1f2}}h1{{color:#087e8b}}section{{background:#f5f8f8;border-left:4px solid #087e8b;padding:10px 16px;margin:14px 0}}.warning{{color:#8b3a22;font-weight:700}}</style></head>
<body><h1>Katana-CoSplice: {escape(self.config.project_name)}</h1>
<p>Sample-level PSI co-splicing network. Modules are replicate-supported descriptive programs; module associations use the same discovery cohort and are prioritization evidence, not independent validation or causal regulation.</p>
<p><a href="katana_cosplice_validation_report.pdf">Open PDF validation report</a> | <a href="module_browser.html">Open module browser</a></p>
{''.join(summary_blocks)}
<h2>Stage provenance</h2>
<table><thead><tr><th>Stage</th><th>Status</th><th>Message</th><th>Updated</th></tr></thead><tbody>{''.join(rows)}</tbody></table>
</body></html>"""
        report_path.parent.mkdir(parents=True, exist_ok=True)
        report_path.write_text(html, encoding="utf-8")
        build_module_browser(self.output, self.config.project_name, browser_path)
        build_validation_pdf(self.output, self.config.project_name, state, pdf_path)
        self._mark(
            stage,
            "finished",
            "HTML, PDF and module browser generated from existing cached outputs.",
            outputs=[report_path, pdf_path, browser_path],
        )
        return report_path

    def annotate_modules(self, force: bool = False) -> dict[str, Path]:
        stage = "annotate_modules"
        assignments_path = self.output / "04_network" / "observed" / "event_module_assignments.tsv"
        event_metadata_path = self.output / "02_matrix" / "event_metadata.tsv"
        required = [assignments_path, event_metadata_path]
        if not all(path.exists() for path in required):
            raise FileNotFoundError("Missing input: observed network not available. Run network_observed first.")
        candidate_paths = (
            sorted(self.config.support_root.glob("*/gene_level_integrated_candidates.tsv"))
            if self.config.support_root is not None and self.config.support_root.exists()
            else []
        )
        optional_inputs = [*candidate_paths]
        if self.config.rbp_gene_set is not None:
            optional_inputs.append(self.config.rbp_gene_set)
        out_dir = self.output / "06b_module_annotation"
        outputs = {
            "events": out_dir / "event_module_annotation.tsv",
            "genes": out_dir / "module_gene_annotation.tsv",
            "summary": out_dir / "module_evidence_summary.tsv",
            "rbp": out_dir / "module_rbp_enrichment.tsv",
            "sources": out_dir / "support_source_audit.tsv",
            "parameters": out_dir / "parameters.json",
        }
        settings = {
            "network": "observed",
            "support_root": str(self.config.support_root or ""),
            "rbp_gene_set": str(self.config.rbp_gene_set or ""),
            "support_semantics": "annotation_only",
        }
        stage_fp = self._fingerprint(stage, [*required, *optional_inputs], settings)
        state = self._state().get("stages", {}).get(stage, {})
        if not force and state.get("fingerprint") == stage_fp and all(path.exists() for path in outputs.values()):
            cached = [*outputs.values(), *sorted((out_dir / "plots").glob("*.png"))]
            self._mark(stage, "cache_hit", "Loaded existing module evidence annotation cache.", stage_fp, cached)
            return outputs
        self._mark(stage, "running", "Annotating observed modules with cached DEG/DTU/RBP support.", stage_fp)
        try:
            assignments = pd.read_csv(assignments_path, sep="\t", low_memory=False)
            event_metadata = pd.read_csv(event_metadata_path, sep="\t", low_memory=False)
            result = annotate_modules(
                assignments,
                event_metadata,
                support_root=self.config.support_root,
                rbp_gene_set=self.config.rbp_gene_set,
            )
            out_dir.mkdir(parents=True, exist_ok=True)
            result.events.to_csv(outputs["events"], sep="\t", index=False, float_format="%.6g")
            result.genes.to_csv(outputs["genes"], sep="\t", index=False, float_format="%.6g")
            result.module_summary.to_csv(outputs["summary"], sep="\t", index=False, float_format="%.6g")
            result.rbp_enrichment.to_csv(outputs["rbp"], sep="\t", index=False, float_format="%.6g")
            result.support_sources.to_csv(outputs["sources"], sep="\t", index=False)
            plots = render_module_annotation_plot(result.module_summary, out_dir / "plots")
            write_json(
                outputs["parameters"],
                {
                    "fingerprint": stage_fp,
                    **settings,
                    "n_annotated_events": len(result.events),
                    "n_unique_module_genes": result.genes["gene_key"].nunique() if not result.genes.empty else 0,
                    "n_modules": result.module_summary["module"].nunique() if not result.module_summary.empty else 0,
                    "generated_plot_files": [str(path) for path in plots],
                    "interpretation_limit": (
                        "DEG, DTU and RBP overlap annotate/prioritize PSI-defined modules. "
                        "They do not establish causal RBP regulation or redefine module membership."
                    ),
                },
            )
            written = [*outputs.values(), *plots]
            self._mark(
                stage,
                "finished",
                f"Annotated {len(result.events)} events and {result.genes['gene_key'].nunique() if not result.genes.empty else 0} unique genes.",
                stage_fp,
                written,
            )
            return outputs
        except Exception as exc:
            self._mark(stage, "failed", str(exc), stage_fp)
            raise

    def interpret_modules(self, force: bool = False) -> dict[str, Path]:
        stage = "interpret_modules"
        consensus_dir = self.output / "04c_consensus_modules"
        observed_dir = self.output / "04_network" / "observed"
        consensus_sources = [
            consensus_dir / "consensus_event_module_assignments.tsv",
            consensus_dir / "consensus_module_program_scores.tsv",
        ]
        observed_sources = [
            observed_dir / "event_module_assignments.tsv",
            observed_dir / "module_program_scores.tsv",
        ]
        if all(path.exists() for path in consensus_sources):
            module_sources = consensus_sources
            module_source = "bootstrap_consensus"
        elif all(path.exists() for path in observed_sources):
            module_sources = observed_sources
            module_source = "observed_network_fallback"
        else:
            raise FileNotFoundError(
                "Missing input: no cached consensus or observed module result. Run network_observed first; "
                "consensus_modules is preferred for final interpretation."
            )
        required = [
            *module_sources,
            self.output / "03_stable_as_qc" / "stable_psi_matrix.tsv",
            self.output / "02_matrix" / "event_metadata.tsv",
        ]
        if not all(path.exists() for path in required):
            raise FileNotFoundError("Missing input: module, Stable-AS, or event metadata cache is incomplete.")
        optional_inputs: list[Path] = []
        for path in (
            self.config.rbp_gene_set,
            self.config.interpretation.rbp_expression_matrix,
            self.config.interpretation.rbp_motif_hits,
            self.config.interpretation.rbp_motif_database,
            self.config.interpretation.gene_sets_gmt,
            self.config.gtf,
            self.config.genome_fasta,
        ):
            if path is not None:
                optional_inputs.append(path)
        if self.config.quant_root is not None and self.config.quant_root.exists():
            optional_inputs.extend(sorted(self.config.quant_root.glob("**/quant.sf")))
        if self.config.support_root is not None and self.config.support_root.exists():
            optional_inputs.extend(sorted(self.config.support_root.glob("*/gene_level_integrated_candidates.tsv")))
        out_dir = self.output / "06e_module_interpretation"
        outputs = {
            "hub_events": out_dir / "module_hub_events.tsv",
            "hub_genes": out_dir / "module_hub_genes.tsv",
            "joint": out_dir / "module_joint_evidence.tsv",
            "motif": out_dir / "module_rbp_motif_enrichment.tsv",
            "motif_hits": out_dir / "event_rbp_motif_hits.tsv",
            "rbp_expression": out_dir / "rbp_module_expression_correlation.tsv",
            "rbp_expression_matrix": out_dir / "rbp_gene_expression_matrix.tsv",
            "gene_sets": out_dir / "module_go_kegg_enrichment.tsv",
            "isoform_switches": out_dir / "dominant_isoform_switches.tsv",
            "isoform_enrichment": out_dir / "module_isoform_switch_enrichment.tsv",
            "regulators": out_dir / "module_regulator_prioritization.tsv",
            "audit": out_dir / "interpretation_source_audit.tsv",
            "motif_parameters": out_dir / "motif_scan_parameters.json",
            "isoform_parameters": out_dir / "isoform_switch_parameters.json",
            "expression_parameters": out_dir / "rbp_expression_parameters.json",
            "parameters": out_dir / "parameters.json",
        }
        motif_inputs = [
            path for path in (
                self.output / "02_matrix" / "event_metadata.tsv",
                self.config.interpretation.rbp_motif_database,
                self.config.genome_fasta,
            ) if path is not None
        ]
        motif_settings = {
            "scanner_version": "2.0-merged-windows-preserve-motif-id",
            "database": str(self.config.interpretation.rbp_motif_database or ""),
            "organism": self.config.interpretation.rbp_motif_organism,
            "window_bp": self.config.interpretation.motif_window_bp,
            "score_threshold": self.config.interpretation.motif_score_threshold,
            "strand_oriented": True,
        }
        motif_fingerprint = self._fingerprint("interpret_modules_motif_scan", motif_inputs, motif_settings)
        motif_override: Path | None = None
        if outputs["motif_hits"].exists() and outputs["motif_parameters"].exists():
            import json
            cached_motif = json.loads(outputs["motif_parameters"].read_text(encoding="utf-8"))
            if cached_motif.get("fingerprint") == motif_fingerprint:
                motif_override = outputs["motif_hits"]
        sample_metadata_path = self.output / "02_matrix" / "sample_metadata_resolved.tsv"
        isoform_inputs = [
            path for path in (
                required[0], self.config.gtf, sample_metadata_path,
            ) if path is not None
        ]
        if self.config.quant_root is not None and self.config.quant_root.exists():
            isoform_inputs.extend(sorted(self.config.quant_root.glob("**/quant.sf")))
        isoform_settings = {
            "quant_root": str(self.config.quant_root or ""),
            "gtf": str(self.config.gtf or ""),
            "condition_column": self.config.statistics.condition_column,
            "min_gene_tpm": self.config.interpretation.isoform_min_gene_tpm,
            "min_abs_delta_if": self.config.interpretation.isoform_min_abs_delta_if,
        }
        isoform_fingerprint = self._fingerprint(
            "interpret_modules_isoform_switch", isoform_inputs, isoform_settings
        )
        isoform_override: Path | None = None
        if outputs["isoform_switches"].exists() and outputs["isoform_parameters"].exists():
            import json
            cached_isoform = json.loads(outputs["isoform_parameters"].read_text(encoding="utf-8"))
            if cached_isoform.get("fingerprint") == isoform_fingerprint:
                isoform_override = outputs["isoform_switches"]
        expression_inputs = [
            path for path in (self.config.gtf, self.config.interpretation.rbp_expression_matrix)
            if path is not None
        ]
        if self.config.interpretation.rbp_expression_matrix is None and self.config.quant_root is not None and self.config.quant_root.exists():
            expression_inputs.extend(sorted(self.config.quant_root.glob("**/quant.sf")))
        expression_settings = {
            "configured_matrix": str(self.config.interpretation.rbp_expression_matrix or ""),
            "quant_root": str(self.config.quant_root or ""),
            "gtf": str(self.config.gtf or ""),
            "aggregation": "transcript TPM summed to GTF gene_symbol",
        }
        expression_fingerprint = self._fingerprint(
            "interpret_modules_rbp_expression", expression_inputs, expression_settings
        )
        expression_override: Path | None = None
        if outputs["rbp_expression_matrix"].exists() and outputs["expression_parameters"].exists():
            import json
            cached_expression = json.loads(outputs["expression_parameters"].read_text(encoding="utf-8"))
            if cached_expression.get("fingerprint") == expression_fingerprint:
                expression_override = outputs["rbp_expression_matrix"]
        settings = {
            "interpretation": self.config.interpretation.__dict__,
            "support_root": str(self.config.support_root or ""),
            "rbp_gene_set": str(self.config.rbp_gene_set or ""),
            "module_source": module_source,
            "scientific_scope": "regulator prioritization; not causal inference",
        }
        stage_fp = self._fingerprint(stage, [*required, *optional_inputs], settings)
        state = self._state().get("stages", {}).get(stage, {})
        if not force and state.get("fingerprint") == stage_fp and all(path.exists() for path in outputs.values()):
            cached = [*outputs.values(), *sorted((out_dir / "plots").glob("*.png"))]
            self._mark(stage, "cache_hit", "Loaded existing module interpretation cache.", stage_fp, cached)
            return outputs
        self._mark(stage, "running", "Prioritizing module hubs, RBP support and functional gene sets.", stage_fp)
        try:
            result = interpret_modules(
                assignments=pd.read_csv(required[0], sep="\t", low_memory=False),
                program_scores=pd.read_csv(required[1], sep="\t", index_col=0),
                psi=pd.read_csv(required[2], sep="\t", index_col=0),
                event_metadata=pd.read_csv(required[3], sep="\t", low_memory=False),
                settings=self.config.interpretation,
                support_root=self.config.support_root,
                rbp_gene_set=self.config.rbp_gene_set,
                quant_root=self.config.quant_root,
                gtf=self.config.gtf,
                genome_fasta=self.config.genome_fasta,
                motif_hits_override=motif_override,
                isoform_switches_override=isoform_override,
                rbp_expression_override=expression_override,
                sample_metadata=(
                    pd.read_csv(sample_metadata_path, sep="\t", low_memory=False)
                    if sample_metadata_path.exists() else None
                ),
                condition_column=self.config.statistics.condition_column,
            )
            out_dir.mkdir(parents=True, exist_ok=True)
            tables = {
                "hub_events": result.hub_events,
                "hub_genes": result.hub_genes,
                "joint": result.joint_evidence,
                "motif": result.rbp_motif_enrichment,
                "motif_hits": result.rbp_motif_hits,
                "rbp_expression": result.rbp_expression_correlation,
                "gene_sets": result.gene_set_enrichment,
                "isoform_switches": result.isoform_switches,
                "isoform_enrichment": result.module_isoform_switch_enrichment,
                "regulators": result.regulator_ranking,
                "audit": result.source_audit,
            }
            for key, frame in tables.items():
                frame.to_csv(outputs[key], sep="\t", index=False, float_format="%.6g")
            result.rbp_expression_matrix.to_csv(
                outputs["rbp_expression_matrix"], sep="\t", index=True, float_format="%.6g"
            )
            write_json(outputs["motif_parameters"], {
                "fingerprint": motif_fingerprint,
                **motif_settings,
                "cache_status": "cache_hit" if motif_override is not None else "generated",
                "n_hit_rows": len(result.rbp_motif_hits),
            })
            write_json(outputs["isoform_parameters"], {
                "fingerprint": isoform_fingerprint,
                **isoform_settings,
                "cache_status": "cache_hit" if isoform_override is not None else "generated",
                "n_switch_rows": len(result.isoform_switches),
            })
            write_json(outputs["expression_parameters"], {
                "fingerprint": expression_fingerprint,
                **expression_settings,
                "cache_status": "cache_hit" if expression_override is not None else "generated",
                "n_gene_rows": len(result.rbp_expression_matrix),
            })
            plots = render_module_interpretation_plots(result, out_dir / "plots")
            write_json(outputs["parameters"], {
                "fingerprint": stage_fp, **settings,
                "generated_plot_files": [str(path) for path in plots],
                "method_basis": {
                    "hub": "WGCNA-style module eigengene membership (kME) plus consensus confidence/connectivity",
                    "motif": "rMAPS2-inspired positional motif-hit enrichment when a validated hit table is configured",
                    "enrichment": "one-sided hypergeometric test with Benjamini-Hochberg FDR",
                },
            })
            written = [*outputs.values(), *plots]
            self._mark(stage, "finished", f"Interpreted {result.hub_events['module'].nunique() if not result.hub_events.empty else 0} modules.", stage_fp, written)
            return outputs
        except Exception as exc:
            self._mark(stage, "failed", str(exc), stage_fp)
            raise

    def fit_module_models(self, force: bool = False) -> dict[str, Path]:
        stage = "fit_module_models"
        matrix_dir = self.output / "02_matrix"
        consensus_dir = self.output / "04c_consensus_modules"
        observed_dir = self.output / "04_network" / "observed"
        consensus_sources = [
            consensus_dir / "consensus_event_module_assignments.tsv",
            consensus_dir / "consensus_module_program_scores.tsv",
        ]
        observed_sources = [
            observed_dir / "event_module_assignments.tsv",
            observed_dir / "module_program_scores.tsv",
        ]
        if all(path.exists() for path in consensus_sources):
            module_sources = consensus_sources
            module_source = "bootstrap_consensus"
        elif all(path.exists() for path in observed_sources):
            module_sources = observed_sources
            module_source = "observed_network_fallback"
        else:
            raise FileNotFoundError(
                "Missing input: no cached consensus or observed module result. Run network_observed first; "
                "consensus_modules is preferred for final module statistics."
            )
        required = [
            matrix_dir / "sample_level_psi_matrix.tsv",
            matrix_dir / "sample_level_ijc_matrix.tsv",
            matrix_dir / "sample_level_sjc_matrix.tsv",
            matrix_dir / "sample_metadata_resolved.tsv",
            *module_sources,
        ]
        if not all(path.exists() for path in required):
            raise FileNotFoundError("Missing input: matrix or module cache is incomplete.")
        out_dir = self.output / "06f_module_statistics"
        outputs = {
            "events": out_dir / "event_beta_binomial_condition_effect.tsv",
            "modules": out_dir / "module_program_condition_effect.tsv",
            "variance": out_dir / "module_variance_partition.tsv",
            "benchmark": out_dir / "beta_binomial_vs_simple_delta_dpsi.tsv",
            "diagnostics": out_dir / "beta_binomial_model_diagnostics.tsv",
            "parameters": out_dir / "parameters.json",
        }
        settings = {
            "statistics": self.config.statistics.__dict__,
            "optional_batch_column": self.config.batch_column,
            "module_source": module_source,
            "excluded_covariates": "sex, age and cell type are not modeled by this stage",
        }
        stage_fp = self._fingerprint(stage, required, settings)
        state = self._state().get("stages", {}).get(stage, {})
        if not force and state.get("fingerprint") == stage_fp and all(path.exists() for path in outputs.values()):
            cached = [*outputs.values(), *sorted((out_dir / "plots").glob("*.png"))]
            self._mark(stage, "cache_hit", "Loaded existing replicate-aware statistics cache.", stage_fp, cached)
            return outputs
        self._mark(stage, "running", "Fitting beta-binomial event models and robust module program-score models.", stage_fp)
        try:
            result = fit_replicate_aware_models(
                psi=pd.read_csv(required[0], sep="\t", index_col=0),
                ijc=pd.read_csv(required[1], sep="\t", index_col=0),
                sjc=pd.read_csv(required[2], sep="\t", index_col=0),
                samples=pd.read_csv(required[3], sep="\t", low_memory=False),
                assignments=pd.read_csv(required[4], sep="\t", low_memory=False),
                program_scores=pd.read_csv(required[5], sep="\t", index_col=0),
                settings=self.config.statistics,
                batch_column=self.config.batch_column,
            )
            out_dir.mkdir(parents=True, exist_ok=True)
            for key, frame in {
                "events": result.event_models,
                "modules": result.module_models,
                "variance": result.variance_partition,
                "benchmark": result.simple_delta_benchmark,
                "diagnostics": result.model_diagnostics,
            }.items():
                frame.to_csv(outputs[key], sep="\t", index=False, float_format="%.6g")
            plots = render_module_statistics_plots(result, out_dir / "plots")
            write_json(outputs["parameters"], {
                "fingerprint": stage_fp, **settings,
                "generated_plot_files": [str(path) for path in plots],
                "model_scope": (
                    "The beta-binomial LRT is an additional replicate/read-count-aware condition test. "
                    "It does not replace rMATS FDR and is not a mixed model when no repeated donor/block identifier exists."
                ),
            })
            written = [*outputs.values(), *plots]
            self._mark(stage, "finished", f"Fitted {len(result.event_models)} event models and {len(result.module_models)} module models.", stage_fp, written)
            return outputs
        except Exception as exc:
            self._mark(stage, "failed", str(exc), stage_fp)
            raise

    def meta_analyze_modules(self, force: bool = False) -> dict[str, Path]:
        stage = "meta_analyze_modules"
        consensus_assignments = self.output / "04c_consensus_modules" / "consensus_event_module_assignments.tsv"
        observed_assignments = self.output / "04_network" / "observed" / "event_module_assignments.tsv"
        if consensus_assignments.exists():
            assignment_path = consensus_assignments
            module_source = "bootstrap_consensus"
        elif observed_assignments.exists():
            assignment_path = observed_assignments
            module_source = "observed_network_fallback"
        else:
            raise FileNotFoundError(
                "Missing input: no cached consensus or observed module assignments. Run network_observed first; "
                "consensus_modules is preferred for cross-cohort validation."
            )
        required = [
            self.output / "03_stable_as_qc" / "stable_psi_matrix.tsv",
            assignment_path,
            self.output / "02_matrix" / "sample_metadata_resolved.tsv",
        ]
        if not all(path.exists() for path in required):
            raise FileNotFoundError("Missing input: stable PSI, module assignments, or resolved sample metadata is incomplete.")
        cohort_inputs: list[Path] = []
        for cohort in self.config.external_cohorts:
            cohort_inputs.extend([cohort.psi_matrix, cohort.sample_metadata])
            if cohort.module_assignments is not None:
                cohort_inputs.append(cohort.module_assignments)
            if cohort.effect_table is not None:
                cohort_inputs.append(cohort.effect_table)
            for path in (cohort.event_metadata, cohort.event_mapping_table, cohort.preservation_table):
                if path is not None:
                    cohort_inputs.append(path)
        out_dir = self.output / "06g_cross_cohort_meta"
        outputs = {
            "alignment": out_dir / "exact_event_alignment.tsv",
            "matching": out_dir / "cross_cohort_module_matching.tsv",
            "states": out_dir / "cross_cohort_module_state.tsv",
            "effects": out_dir / "cohort_module_effects.tsv",
            "meta": out_dir / "random_effects_module_meta_analysis.tsv",
            "audit": out_dir / "meta_analysis_cohort_audit.tsv",
            "consensus": out_dir / "multi_cohort_module_consensus.tsv",
            "parameters": out_dir / "parameters.json",
        }
        settings = {
            "meta_analysis": self.config.meta_analysis.__dict__,
            "cohorts": [cohort.__dict__ for cohort in self.config.external_cohorts],
            "module_source": module_source,
            "event_identity": (
                "exact key plus same-build exact splice-site signature and optional explicit liftover mapping; "
                "no arbitrary coordinate tolerance"
            ),
        }
        stage_fp = self._fingerprint(stage, [*required, *cohort_inputs], settings)
        state = self._state().get("stages", {}).get(stage, {})
        if not force and state.get("fingerprint") == stage_fp and all(path.exists() for path in outputs.values()):
            cached = [*outputs.values(), *sorted((out_dir / "plots").glob("*.png"))]
            self._mark(stage, "cache_hit", "Loaded existing cross-cohort meta-analysis cache.", stage_fp, cached)
            return outputs
        self._mark(stage, "running", "Aligning exact AS events and meta-analyzing independent cohort module effects.", stage_fp)
        try:
            result = meta_analyze_modules(
                reference_psi=pd.read_csv(required[0], sep="\t", index_col=0),
                reference_assignments=pd.read_csv(required[1], sep="\t", low_memory=False),
                reference_metadata=pd.read_csv(required[2], sep="\t", low_memory=False),
                external_cohorts=self.config.external_cohorts,
                settings=self.config.meta_analysis,
                reference_name=self.config.project_name,
                reference_genome_build=self.config.genome_build,
            )
            out_dir.mkdir(parents=True, exist_ok=True)
            for key, frame in {
                "alignment": result.event_alignment,
                "matching": result.module_matching,
                "states": result.module_state,
                "effects": result.cohort_effects,
                "meta": result.random_effects_meta,
                "audit": result.cohort_audit,
                "consensus": result.multi_cohort_consensus,
            }.items():
                frame.to_csv(outputs[key], sep="\t", index=False, float_format="%.6g")
            plots = render_meta_analysis_plots(result, out_dir / "plots")
            write_json(outputs["parameters"], {
                "fingerprint": stage_fp, **settings,
                "generated_plot_files": [str(path) for path in plots],
                "scientific_boundary": (
                    "Independent target module assignments are matched on exact events and tested by hypergeometric overlap with BH-FDR. "
                    "Rewiring is not inferred from Jaccard alone. Without target assignments, structural state is not_assessed."
                ),
            })
            written = [*outputs.values(), *plots]
            self._mark(
                stage, "finished",
                f"Cross-cohort meta-analysis processed {len(self.config.external_cohorts)} configured cohorts using {module_source}.",
                stage_fp, written,
            )
            return outputs
        except Exception as exc:
            self._mark(stage, "failed", str(exc), stage_fp)
            raise

    def factor_contrasts(self, force: bool = False) -> dict[str, Path]:
        """Run only the configured factor/phenocopy contrasts from cached upstream artifacts."""
        stage = "factor_contrasts"
        settings = self.config.factor_contrasts
        if not settings.enabled:
            raise ValueError(
                "Configured factor contrasts are disabled. Set factor_contrasts.enabled=true and predeclare "
                "the primary, context and perturbation factors before running 6H."
            )
        score_path = (
            self.output
            / "04_network"
            / settings.network_source
            / "module_program_scores.tsv"
        )
        required = [
            self.output / "03_stable_as_qc" / "stable_psi_matrix.tsv",
            self.output / "02_matrix" / "sample_metadata_resolved.tsv",
            self.output / "02_matrix" / "event_metadata.tsv",
            score_path,
        ]
        missing = [str(path) for path in required if not path.exists()]
        if missing:
            raise FileNotFoundError(
                "Missing cached input for 6H factor contrasts: "
                + "; ".join(missing)
                + ". Run stable_qc and the configured network stage first. 6H will not run them automatically."
            )
        out_dir = self.output / "06h_factor_contrasts"
        outputs = {
            "events": out_dir / "event_factor_contrasts.tsv",
            "modules": out_dir / "module_factor_contrasts.tsv",
            "phenocopy": out_dir / "phenocopy_occlusion_summary.tsv",
            "evidence": out_dir / "factor_evidence_summary.tsv",
            "notes": out_dir / "presentation_notes.md",
            "parameters": out_dir / "parameters.json",
        }
        stage_settings = {
            "factor_contrasts": settings.__dict__,
            "network_source": settings.network_source,
            "scope": (
                "Predeclared stable-AS target-reference, perturbation and difference-in-differences contrasts; "
                "phenocopy support is not causal or functional-state proof."
            ),
        }
        stage_fp = self._fingerprint(stage, required, stage_settings)
        state = self._state().get("stages", {}).get(stage, {})
        if (
            not force
            and state.get("fingerprint") == stage_fp
            and all(path.exists() for path in outputs.values())
        ):
            cached = [*outputs.values(), *sorted((out_dir / "plots").glob("*.png"))]
            self._mark(stage, "cache_hit", "Loaded existing configured factor-contrast cache.", stage_fp, cached)
            return outputs

        self._mark(
            stage,
            "running",
            "Testing predeclared stable-AS and CoSplice-program factor contrasts; no upstream stage is being run.",
            stage_fp,
        )
        try:
            result = analyze_factor_contrasts(
                psi=pd.read_csv(required[0], sep="\t", index_col=0),
                samples=pd.read_csv(required[1], sep="\t", low_memory=False),
                event_metadata=pd.read_csv(required[2], sep="\t", low_memory=False),
                program_scores=pd.read_csv(required[3], sep="\t", index_col=0),
                settings=settings,
            )
            out_dir.mkdir(parents=True, exist_ok=True)
            result.event_contrasts.to_csv(outputs["events"], sep="\t", index=False, float_format="%.6g")
            result.module_contrasts.to_csv(outputs["modules"], sep="\t", index=False, float_format="%.6g")
            result.phenocopy_summary.to_csv(outputs["phenocopy"], sep="\t", index=False, float_format="%.6g")
            result.evidence_summary.to_csv(outputs["evidence"], sep="\t", index=False, float_format="%.6g")
            plots = render_factor_contrast_plots(result, out_dir / "plots", settings=settings)
            notes = [
                "# Configured factor contrasts and phenocopy evidence",
                "",
                "## Analysis scope",
                "",
                f"- Primary contrast: `{settings.target_level}` versus `{settings.reference_level}` within baseline `{settings.baseline_level}`.",
                f"- Context factors: `{', '.join(settings.context_factors) or 'none'}`.",
                f"- Perturbations: `{', '.join(settings.perturbation_levels)}` from factor `{settings.perturbation_factor}`.",
                f"- Program-score source: `{settings.network_source}` network.",
                "- Event tests use stable sample-level PSI and replicate-aware Welch contrasts.",
                "- Occlusion uses a difference-in-differences contrast: target perturbation effect minus reference perturbation effect.",
                "",
                "## Interpretation boundary",
                "",
                "Concordant event signatures, enriched overlap and perturbation attenuation provide phenocopy or occlusion support only. "
                "They do not prove loss of function, gain of function, pathway causality or direct regulation. Functional-state claims "
                "require orthogonal molecular evidence and preferably rescue/complementation.",
                "",
                "## Results",
                "",
            ]
            if result.evidence_summary.empty:
                notes.append("No configured contrast had enough samples to summarize.")
            else:
                notes.append(result.evidence_summary.to_markdown(index=False, floatfmt=".3g"))
            significant_interactions = result.module_contrasts.loc[
                result.module_contrasts["contrast_type"].eq("perturbation_occlusion_interaction")
                & result.module_contrasts["significant"].fillna(False)
            ].copy()
            notes.extend(
                [
                    "",
                    "## Significant program interactions",
                    "",
                    "A significant difference-in-differences result is not automatically occlusion. "
                    "`attenuated_same_direction` is compatible with occlusion; `opposite_response` is not.",
                    "",
                ]
            )
            if significant_interactions.empty:
                notes.append("No CoSplice program interaction passed the configured FDR cutoff.")
            else:
                interaction_columns = [
                    "module",
                    "context",
                    "perturbation",
                    "target_perturbation_effect",
                    "reference_perturbation_effect",
                    "effect",
                    "FDR",
                    "interaction_pattern",
                ]
                notes.append(
                    significant_interactions.loc[:, interaction_columns].to_markdown(index=False, floatfmt=".3g")
                )
            notes.extend(
                [
                    "",
                    "## Functional-state conclusion",
                    "",
                    "This RNA-level analysis alone is insufficient to classify the target as LOF or GOF. "
                    "A functional-state claim requires orthogonal target-activity evidence and rescue/complementation.",
                ]
            )
            outputs["notes"].write_text("\n".join(notes) + "\n", encoding="utf-8")
            write_json(
                outputs["parameters"],
                {
                    "fingerprint": stage_fp,
                    **stage_settings,
                    "generated_plot_files": [str(path) for path in plots],
                    "event_effect_definition": "mean PSI(target) - mean PSI(reference)",
                    "interaction_effect_definition": (
                        "[target perturbation - target baseline] - "
                        "[reference perturbation - reference baseline]"
                    ),
                    "scientific_limit": (
                        "Descriptive and replicate-aware prioritization; this stage is not a causal LOF/GOF classifier."
                    ),
                },
            )
            written = [*outputs.values(), *plots]
            self._mark(
                stage,
                "finished",
                f"Configured factor contrasts finished: {len(result.event_contrasts):,} event rows, "
                f"{len(result.module_contrasts):,} module rows and {len(result.phenocopy_summary)} phenocopy summaries.",
                stage_fp,
                written,
            )
            return outputs
        except Exception as exc:
            self._mark(stage, "failed", str(exc), stage_fp)
            raise

    def benchmark_methods(self, force: bool = False) -> dict[str, Path]:
        stage = "benchmark_methods"
        matrix_path = self.output / "03_stable_as_qc" / "stable_psi_matrix.tsv"
        qc_path = self.output / "03_stable_as_qc" / "stable_as_event_qc.tsv"
        samples_path = self.output / "02_matrix" / "sample_metadata_resolved.tsv"
        required = [matrix_path, qc_path, samples_path]
        if not all(path.exists() for path in required):
            raise FileNotFoundError("Missing input: Stable-AS QC cache is incomplete. Run matrix and stable_qc first.")
        out_dir = self.output / "06c_method_benchmark"
        outputs = {
            "universe": out_dir / "common_event_universe.tsv",
            "assignments": out_dir / "event_method_assignments.tsv",
            "summary": out_dir / "method_benchmark_summary.tsv",
            "bootstrap": out_dir / "bootstrap_partition_stability.tsv",
            "matching": out_dir / "module_matching.tsv",
            "audit": out_dir / "method_availability.tsv",
            "metric_ranking": out_dir / "method_metric_ranking.tsv",
            "interpretation": out_dir / "method_strengths_limitations.tsv",
            "report": out_dir / "method_benchmark_report.html",
            "parameters": out_dir / "parameters.json",
        }
        settings = {
            "network": self.config.network.__dict__,
            "benchmark": self.config.benchmark.__dict__,
            "comparison_semantics": "same_input_internal_head_to_head",
        }
        stage_fp = self._fingerprint(stage, required, settings)
        state = self._state().get("stages", {}).get(stage, {})
        if not force and state.get("fingerprint") == stage_fp and all(path.exists() for path in outputs.values()):
            cached = [*outputs.values(), *sorted((out_dir / "plots").glob("*.png"))]
            self._mark(stage, "cache_hit", "Loaded existing method benchmark cache.", stage_fp, cached)
            return outputs
        self._mark(stage, "running", "Benchmarking methods on one shared stable-AS event universe.", stage_fp)
        try:
            matrix = pd.read_csv(matrix_path, sep="\t", index_col=0)
            qc = pd.read_csv(qc_path, sep="\t", low_memory=False)
            samples = pd.read_csv(samples_path, sep="\t", low_memory=False)
            result = run_method_benchmark(
                matrix,
                qc,
                samples,
                self.config.network,
                self.config.benchmark,
                progress=self._emit,
            )
            out_dir.mkdir(parents=True, exist_ok=True)
            result.universe.to_csv(outputs["universe"], sep="\t", index=False, float_format="%.6g")
            result.assignments.to_csv(outputs["assignments"], sep="\t", index=False)
            result.method_summary.to_csv(outputs["summary"], sep="\t", index=False, float_format="%.6g")
            result.bootstrap_stability.to_csv(outputs["bootstrap"], sep="\t", index=False, float_format="%.6g")
            result.module_matching.to_csv(outputs["matching"], sep="\t", index=False, float_format="%.6g")
            result.method_audit.to_csv(outputs["audit"], sep="\t", index=False, float_format="%.6g")
            result.metric_ranking.to_csv(outputs["metric_ranking"], sep="\t", index=False, float_format="%.6g")
            result.method_interpretation.to_csv(outputs["interpretation"], sep="\t", index=False, float_format="%.6g")
            plots = render_benchmark_plots(result, out_dir / "plots")
            render_benchmark_report(result, outputs["report"])
            write_json(
                outputs["parameters"],
                {
                    "fingerprint": stage_fp,
                    **settings,
                    "n_common_events": len(result.universe),
                    "n_completed_methods": int(result.method_summary["status"].eq("finished").sum()),
                    "generated_plots": [str(path) for path in plots],
                    "interpretation_limit": (
                        "Internal same-cohort benchmark only. Higher internal stability or separation does not prove "
                        "external biological superiority; independent-cohort module preservation remains required."
                    ),
                },
            )
            written = [*outputs.values(), *plots]
            completed = int(result.method_summary["status"].eq("finished").sum())
            self._mark(stage, "finished", f"Method benchmark finished: {completed} methods completed.", stage_fp, written)
            return outputs
        except Exception as exc:
            self._mark(stage, "failed", str(exc), stage_fp)
            raise

    def run_all(self, include_residual: bool = True, force: bool = False) -> None:
        """Explicit full pipeline command; individual stage calls never use this method."""
        self.run_stage("audit")
        self.run_stage("matrix", force=force)
        self.run_stage("stable_qc", force=force)
        self.run_stage("network_observed", force=force)
        self.run_stage("consensus_modules", force=force)
        self.run_stage("associate_observed", force=force)
        if include_residual and self.config.residual_covariates:
            self.run_stage("network_residual", force=force)
            self.run_stage("associate_residual", force=force)
        self.run_stage("hybrid_programs", force=force)
        self.run_stage("validate_modules", force=force)
        self.run_stage("annotate_modules", force=force)
        self.run_stage("benchmark_methods", force=force)
        self.run_stage("preserve_modules", force=force)
        self.run_stage("interpret_modules", force=force)
        self.run_stage("fit_module_models", force=force)
        self.run_stage("meta_analyze_modules", force=force)
        if self.config.factor_contrasts.enabled:
            self.run_stage("factor_contrasts", force=force)
        self.run_stage("report")
