from __future__ import annotations

from pathlib import Path
from typing import Any
from dataclasses import asdict
from functools import wraps
import shutil
from threading import RLock

import pandas as pd
import yaml


def _read_engine_table(source: Path, separator: str) -> pd.DataFrame | None:
    try:
        return pd.read_csv(source, sep=separator, low_memory=False)
    except pd.errors.EmptyDataError:
        if source.name not in {
            "observed_residual_preservation_pairs.tsv",
            "observed_residual_preservation_summary.tsv",
        }:
            raise
        return None

from katana_cosplice.config import load_config
from katana_cosplice.pipeline import CoSplicePipeline
from katana_cosplice.qc import network_input_policy
from katana_cosplice.stage_registry import stage_spec
from katana_next.domain.models import ArtifactKind
from katana_next.infrastructure.jobs import JobContext
from katana_next.infrastructure.fingerprints import fingerprint_file, fingerprint_value
from katana_next.services.project_service import ProjectContext
from katana_next.workflows.base import (
    ArtifactWriter,
    ParameterDefinition,
    WorkflowBase,
    WorkflowDefinition,
)


DISCOVERY_STAGES = (
    "audit",
    "matrix",
    "stable_qc",
    "network_observed",
    "consensus_modules",
)
INTERPRETATION_STAGES = (
    "associate_observed",
    "validate_modules",
    "annotate_modules",
    "interpret_modules",
    "fit_module_models",
)
EXTENDED_STAGES = (
    "benchmark_methods",
    "preserve_modules",
    "meta_analyze_modules",
    "factor_contrasts",
    "differential_network",
    "multiview_fusion",
    "bayesian_preservation",
)


_ENGINE_LOCK = RLock()


def _serialized_engine(run):
    @wraps(run)
    def execute(self, context, job, *args, **kwargs):
        while not _ENGINE_LOCK.acquire(timeout=0.5):
            job.check_cancelled()
        try:
            return run(self, context, job, *args, **kwargs)
        finally:
            _ENGINE_LOCK.release()
    return execute


def configuration_input_files(config) -> list[Path]:
    paths: list[Path] = []
    def walk(value, key=""):
        if key == "output_dir":
            return
        if isinstance(value, dict):
            for name, item in value.items():
                walk(item, name)
        elif isinstance(value, (list, tuple)):
            for item in value:
                walk(item, key)
        elif isinstance(value, Path):
            if value.is_file():
                paths.append(value.resolve())
            elif value.is_dir():
                patterns = {"rmats_root": (f"*.MATS.{config.rmats_mode}.txt",),
                            "quant_root": ("quant.sf",), "bam_root": ("*.bam", "*.bai"),
                            "support_root": ("*.tsv", "*.csv", "*.txt")}.get(key, ())
                for pattern in patterns:
                    paths.extend(path.resolve() for path in value.rglob(pattern) if path.is_file())
    walk(asdict(config))
    return sorted(set(paths))


def configure_interpretation_resources(config, manifest, root: Path) -> None:
    """Bridge project resources to the standalone engine, never network membership."""
    from katana_next.resources import write_gmt, sha256
    from katana_next.science.enrichment import declared_gene_sets
    kinds = {a.kind for a in manifest.assets}
    if "gene_sets_gmt" in kinds:
        sets, audit = declared_gene_sets(manifest, "gene_sets_gmt")
        target = root / "declared_gene_sets.gmt"
        write_gmt(sets, target, "Katana declared resources; see resource_provenance.tsv")
        audit.to_csv(root / "resource_provenance.tsv", sep="\t", index=False)
        config.interpretation.gene_sets_gmt = target
    if "rbp_gene_sets" in kinds:
        sets, audit = declared_gene_sets(manifest, "rbp_gene_sets")
        rows = [{"gene_set": term, "gene_symbol": gene if not gene.startswith("ENS") else "",
                 "gene_id": gene if gene.startswith("ENS") else "", "source": "Declared project RBP resource"}
                for term, genes in sets.items() for gene in sorted(genes)]
        target = root / "declared_rbp_membership.tsv"
        pd.DataFrame(rows).to_csv(target, sep="\t", index=False)
        audit.to_csv(root / "rbp_resource_provenance.tsv", sep="\t", index=False)
        config.rbp_gene_set = target
    motifs = [a for a in manifest.assets if a.kind == "rbp_motif_database"]
    if len(motifs) > 1:
        raise ValueError("CoSplice engine accepts one motif database per run; select one explicitly in the project manifest.")
    if motifs:
        asset = motifs[0]
        if asset.metadata.get("sha256") and sha256(asset.path) != asset.metadata["sha256"]:
            raise ValueError("Motif database fingerprint changed")
        species = asset.metadata.get("species", "")
        if not species:
            raise ValueError("Declare the organism for the CoSplice motif database")
        config.interpretation.rbp_motif_database = asset.path
        config.interpretation.rbp_motif_organism = species


def configure_production_membership(config, parameters):
    policy = str(parameters.get("membership_policy", "backbone_confidence_only"))
    if policy not in {"backbone_confidence_only", "backbone_refined_consensus"}:
        raise ValueError("Unknown bootstrap membership policy")
    config.consensus.clustering_method = policy
    if policy == "backbone_confidence_only":
        # An old, larger consensus size cutoff must not silently delete modules
        # when the requested operation is confidence annotation only.
        config.consensus.min_module_size = config.network.min_module_size


def write_production_configuration(config, path: Path) -> Path:
    """Keep the active policy visible beside inert historical YAML settings."""
    payload = config.to_mapping()
    payload["network_input_policy"] = network_input_policy()
    payload["legacy_stable_as_settings_role"] = "ignored historical settings; active eligibility is network_input_policy"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(yaml.safe_dump(payload, sort_keys=False), encoding="utf-8")
    return path


class CoSpliceWorkflow(WorkflowBase):
    version = "3.0.0"
    definition = WorkflowDefinition(
        workflow_id="programs.cosplice",
        workspace_id="programs",
        label="cosplice programs",
        question="Which reproducible sample-level co-splicing programs exist and how do they associate with the declared design?",
        method="Signed blockwise WGCNA/TOM backbone with bootstrap membership confidence annotation; the default preserves backbone membership. Every cohort sample needs IJC+SJC >=10 and valid nonconstant PSI. DAS significance, phenotype and interpretation evidence never define input eligibility or membership.",
        output_summary="Modules, eigengene-like program scores, membership confidence, trait associations, validation and interpretation artifacts.",
        analysis_unit="exact splice event across independent biological samples",
        expected_artifacts=3,
        parameters=(
            ParameterDefinition(
                "scope",
                "Run scope",
                "select",
                "discovery",
                "Discovery builds the WGCNA backbone; Full adds associations, validation, interpretation and configured external analyses.",
                ("discovery", "interpretation", "full"),
            ),
            ParameterDefinition("force_engine", "Force cosplice rebuild", "boolean", False, "Rebuild selected engine stages even when their own fingerprints match."),
            ParameterDefinition("max_events", "Maximum events", "number", 0, "Zero uses every count-supported nonconstant event; a positive value is an explicitly recorded exploratory cap.", minimum=0, step=100),
            ParameterDefinition("bootstrap_iterations", "Bootstrap iterations", "number", 50, "Resamples complete donor/biological blocks within condition profiles. Without declared IDs, sample IDs are assumed independent; replicate numbers never imply pairing.", minimum=10, maximum=1000, step=10),
            ParameterDefinition("uncertainty_role", "Read uncertainty role", "select", "diagnostic_only",
                "Default keeps measured PSI and reports read uncertainty separately. Optional weighting is experimental; raw read-fraction posteriors never replace rMATS PSI.",
                ("diagnostic_only", "raw_reliability_weighted")),
            ParameterDefinition("membership_policy", "Bootstrap membership policy", "select", "backbone_confidence_only",
                "Default annotates confidence without rejecting backbone members. Hard filtering is an experimental sensitivity analysis, not a proven stability improvement.",
                ("backbone_confidence_only", "backbone_refined_consensus")),
        ),
        required_asset_kinds=("cosplice_config",),
        limitations="Formal condition-specific rewiring is disabled or labeled exploratory when independent samples per condition are insufficient.",
    )

    def input_paths(self, context: ProjectContext, parameters: dict[str, Any]) -> list[Path]:
        del parameters
        config_asset = self._config_asset(context)
        paths = [config_asset.path, *configuration_input_files(load_config(config_asset.path))]
        # Project references override standalone configuration in run(). Track
        # their contents too, including deletion rather than silently dropping it.
        annotation = context.manifest.annotation
        paths.extend(path for path in (
            annotation.gtf, annotation.genome_fasta, annotation.transcript_fasta
        ) if path is not None)
        for kind in {"rmats_result", "deg_result", "dtu_result", "dexseq_result", "gene_sets_gmt", "rbp_gene_sets", "rbp_motif_database"}:
            paths.extend(asset.path for asset in context.manifest.assets if asset.kind == kind and asset.path.exists())
        paths.extend(
            path
            for sample in context.manifest.samples
            for path in sample.assets.values()
            if path.exists()
        )
        return list(dict.fromkeys(path.resolve() for path in paths))

    @staticmethod
    def _config_asset(context: ProjectContext):
        asset = next(
            (
                item
                for item in context.manifest.assets
                if item.kind == "cosplice_config" and item.path.is_file()
            ),
            None,
        )
        if asset is None:
            raise ValueError("CoSplice requires a cosplice_config input asset.")
        return asset

    @_serialized_engine
    def run(self, context: ProjectContext, job: JobContext, run_fingerprint: str, parameters: dict[str, Any]) -> list[str]:
        config_asset = self._config_asset(context)
        engine_key = fingerprint_value({"project": context.snapshot.manifest_fingerprint,
            "inputs": [fingerprint_file(path) for path in self.input_paths(context, parameters)],
            "max_events": parameters.get("max_events", 0), "bootstrap_iterations": parameters.get("bootstrap_iterations", 50),
            "uncertainty_role": parameters.get("uncertainty_role", "diagnostic_only"),
            "membership_policy": parameters.get("membership_policy", "backbone_confidence_only"),
            "production_contract": "backbone-annotated-v9-all-sample-input", "network_input_policy": network_input_policy()})
        engine_root = context.manifest.output_root / "engines" / "cosplice" / engine_key[:24]
        engine_root.mkdir(parents=True, exist_ok=True)
        config = load_config(config_asset.path)
        config.output_dir = engine_root
        config.gene_exclusion = dict(context.manifest.metadata.get("gene_filter") or {})
        config.project_name = context.manifest.name
        config.genome_build = context.manifest.annotation.genome_build
        config.gtf = context.manifest.annotation.gtf or config.gtf
        config.genome_fasta = context.manifest.annotation.genome_fasta or config.genome_fasta
        config.transcript_fasta = context.manifest.annotation.transcript_fasta or config.transcript_fasta
        configure_interpretation_resources(config, context.manifest, engine_root)
        from katana_next.services.design_editor import sync_cosplice_design
        sync_cosplice_design(config, context.manifest, engine_root)
        config.network.module_detection = "wgcna_blockwise"
        config.network.engine = "tom_hierarchical"
        configure_production_membership(config, parameters)
        config.gene_prior.enabled = False
        config.multiobjective_selection.enabled = False
        configure_production_uncertainty(config, parameters)
        maximum = int(parameters.get("max_events", 0))
        config.network.max_events = maximum or None
        config.network.bootstrap_iterations = int(parameters.get("bootstrap_iterations", 50))
        config.consensus.bootstrap_iterations = int(parameters.get("bootstrap_iterations", 50))
        rscript = str(context.manifest.metadata.get("rscript_path", "")).strip()
        if rscript:
            config.network.wgcna_rscript = rscript
        resolved_config = engine_root / "katana_next_cosplice.yaml"
        write_production_configuration(config, resolved_config)
        pipeline = CoSplicePipeline(config)
        scope = str(parameters.get("scope", "discovery"))
        if scope == "discovery":
            stages = list(DISCOVERY_STAGES)
        elif scope == "interpretation":
            stages = list(INTERPRETATION_STAGES)
        elif scope == "full":
            stages = [*DISCOVERY_STAGES, *INTERPRETATION_STAGES]
            stages.extend(_enabled_extended_stages(config))
            stages.append("report")
        else:
            raise ValueError(f"Unknown CoSplice run scope: {scope}")
        if scope in {"discovery", "full"}:
            from katana_cosplice.audit import load_sample_metadata
            from katana_cosplice.resampling import bootstrap_plan
            design = load_sample_metadata(config.sample_metadata, config.selected_conditions)
            plan = bootstrap_plan(design, design.sample_id)
            plan.require_resampleable()
            info = plan.audit()
            job.report(.02, f"Bootstrap: {info['resampling_units']} {info['resampling_unit_field']} units; "
                f"{info['resampling_strata']} condition profiles; {info['fixed_singleton_strata']} fixed singleton strata; {info['resampling_scope']}")
        force = bool(parameters.get("force_engine", False))
        for index, stage in enumerate(stages, start=1):
            job.check_cancelled()
            job.report(
                0.04 + 0.76 * (index - 1) / max(1, len(stages)),
                f"CoSplice: {stage_spec(stage).title}",
            )
            pipeline.run_stage(stage, force=force)

        writer = ArtifactWriter(
            context,
            self.definition.workflow_id,
            run_fingerprint,
            parameters,
            self.version,
            common_metadata={
                "question": self.definition.question,
                "method": self.definition.method,
                "analysis_unit": self.definition.analysis_unit,
                "production_backbone": "signed WGCNA blockwiseModules",
                "network_measurement": "Observed sample-level PSI; approximate raw read fractions do not replace PSI",
                "network_input_policy": network_input_policy(),
                "read_uncertainty_role": parameters.get("uncertainty_role", "diagnostic_only"),
                "membership_semantics": "Default confidence annotation preserves the full-data WGCNA partition; the explicitly selected exploratory refinement can reject boundary events. Phenotype, RBP and functional terms never define membership.",
                "interpretation": "Module-trait associations are program-level effects. Individual event support remains available through membership and Native Sashimi.",
                "limitations": self.definition.limitations,
                "resolved_config": str(resolved_config),
                "executed_stages": stages,
            },
        )
        job.report(0.82, "Registering CoSplice engine outputs")
        registered = 0
        for stage in stages:
            for artifact_spec in stage_spec(stage).artifacts:
                source = engine_root / artifact_spec.relative_path
                if not source.is_file():
                    continue
                label = f"{stage_spec(stage).title}: {artifact_spec.key.replace('_', ' ')}"
                parents: tuple[str, ...] = ()
                suffix = source.suffix.lower()
                if suffix not in {".tsv", ".csv"}:
                    # Registered artifacts must not point into the mutable stage cache.
                    destination = writer.directory / f"{stage}__{artifact_spec.key}{suffix}"
                    shutil.copy2(source, destination)
                    source = destination
                if suffix in {".tsv", ".csv"}:
                    separator = "\t" if suffix == ".tsv" else ","
                    frame = _read_engine_table(source, separator)
                    if frame is None:
                        destination = writer.directory / f"{stage}__{artifact_spec.key}{suffix}"
                        shutil.copy2(source, destination)
                        writer.existing_file(
                            label + " (no result rows)", destination,
                            entity_type="cosplice_provenance",
                            metadata={"status": "not_assessed", "reason": "Engine emitted an empty optional observed-residual preservation table; no preservation conclusion available"},
                        )
                        registered += 1
                        continue
                    if "row_id" not in frame.columns:
                        frame.insert(0, "row_id", [f"{stage}:{artifact_spec.key}:{index}" for index in range(len(frame))])
                    writer.table(
                        label,
                        f"{stage}__{artifact_spec.key}.parquet",
                        frame,
                        entity_type=_entity_type(artifact_spec.key),
                        parents=parents,
                    )
                elif suffix in {".png", ".svg"}:
                    writer.existing_file(
                        label,
                        source,
                        entity_type="cosplice_visualization",
                        metadata={"mime_type": "image/png" if suffix == ".png" else "image/svg+xml"},
                    )
                elif suffix in {".html", ".pdf"}:
                    writer.existing_file(
                        label,
                        source,
                        kind=ArtifactKind.REPORT,
                        entity_type="cosplice_report",
                        metadata={"mime_type": "text/html" if suffix == ".html" else "application/pdf"},
                    )
                elif suffix in {".json", ".yaml", ".yml", ".npz"}:
                    writer.existing_file(label, source, entity_type="cosplice_provenance")
                registered += 1
        writer.existing_file(
            "Resolved CoSplice production configuration",
            write_production_configuration(config, writer.directory / "resolved_cosplice.yaml"),
            entity_type="analysis_configuration",
            metadata={"mime_type": "application/yaml"},
        )
        if registered == 0:
            raise RuntimeError("CoSplice completed without any registered output files.")
        job.report(0.97, f"CoSplice registered {registered} engine outputs")
        return writer.finish()


def configure_production_uncertainty(config, parameters):
    role = parameters.get("uncertainty_role", "diagnostic_only")
    if role not in {"diagnostic_only", "raw_reliability_weighted"}:
        raise ValueError("Production networks require observed PSI; unsupported read uncertainty role")
    config.uncertainty.network_mode = role
    config.uncertainty.enabled = True
    config.uncertainty.posterior_edge_diagnostics = False
    if role == "diagnostic_only":
        config.uncertainty.filter_low_reliability = False


def _enabled_extended_stages(config: Any) -> list[str]:
    stages = ["benchmark_methods", "preserve_modules", "meta_analyze_modules"]
    if config.factor_contrasts.enabled:
        stages.append("factor_contrasts")
    if config.differential_network.enabled:
        stages.append("differential_network")
    if config.multiview.enabled:
        stages.append("multiview_fusion")
    if config.external_cohorts and config.bayesian_preservation.enabled:
        stages.append("bayesian_preservation")
    return stages


def _entity_type(key: str) -> str:
    key = key.lower()
    if "event" in key or key in {"assignments", "selected"}:
        return "splice_event"
    if "module" in key or key in {"scores", "programs"}:
        return "cosplice_module"
    if "sample" in key:
        return "sample"
    if "gene" in key or "rbp" in key or "regulator" in key:
        return "gene"
    return "analysis_summary"
