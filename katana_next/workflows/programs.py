from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd

from katana_cosplice.config import load_config
from katana_cosplice.pipeline import CoSplicePipeline
from katana_cosplice.stage_registry import stage_spec
from katana_next.domain.models import ArtifactKind
from katana_next.infrastructure.jobs import JobContext
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


class CoSpliceWorkflow(WorkflowBase):
    version = "2.0.1"
    definition = WorkflowDefinition(
        workflow_id="programs.cosplice",
        workspace_id="programs",
        label="CoSplice production programs",
        question="Which reproducible sample-level co-splicing programs exist and how do they associate with the declared design?",
        method="Signed blockwise WGCNA/TOM backbone followed by bootstrap co-assignment rejection of unstable boundary events; phenotype and interpretation evidence never define membership.",
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
            ParameterDefinition("force_engine", "Force CoSplice rebuild", "boolean", False, "Rebuild selected engine stages even when their own fingerprints match."),
            ParameterDefinition("max_events", "Maximum events", "number", 0, "Zero uses every stable event; a positive value is an explicitly recorded exploratory cap.", minimum=0, step=100),
            ParameterDefinition("bootstrap_iterations", "Bootstrap iterations", "number", 50, "Sample-stratified co-assignment confidence iterations.", minimum=10, maximum=1000, step=10),
        ),
        required_asset_kinds=("cosplice_config",),
        limitations="Formal condition-specific rewiring is disabled or labeled exploratory when independent samples per condition are insufficient.",
    )

    def input_paths(self, context: ProjectContext, parameters: dict[str, Any]) -> list[Path]:
        del parameters
        config_asset = self._config_asset(context)
        paths = [config_asset.path]
        for kind in {"rmats_result", "deg_result", "dtu_result", "dexseq_result", "gene_sets_gmt", "rbp_motif_database"}:
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

    def run(self, context: ProjectContext, job: JobContext, run_fingerprint: str, parameters: dict[str, Any]) -> list[str]:
        config_asset = self._config_asset(context)
        engine_root = context.manifest.output_root / "engines" / "cosplice"
        engine_root.mkdir(parents=True, exist_ok=True)
        config = load_config(config_asset.path)
        config.output_dir = engine_root
        config.project_name = context.manifest.name
        config.genome_build = context.manifest.annotation.genome_build
        config.gtf = context.manifest.annotation.gtf or config.gtf
        config.genome_fasta = context.manifest.annotation.genome_fasta or config.genome_fasta
        config.transcript_fasta = context.manifest.annotation.transcript_fasta or config.transcript_fasta
        config.network.module_detection = "wgcna_blockwise"
        config.network.engine = "tom_hierarchical"
        config.consensus.clustering_method = "backbone_refined_consensus"
        maximum = int(parameters.get("max_events", 0))
        config.network.max_events = maximum or None
        config.network.bootstrap_iterations = int(parameters.get("bootstrap_iterations", 50))
        config.consensus.bootstrap_iterations = int(parameters.get("bootstrap_iterations", 50))
        rscript = str(context.manifest.metadata.get("rscript_path", "")).strip()
        if rscript:
            config.network.wgcna_rscript = rscript
        resolved_config = engine_root / "katana_next_cosplice.yaml"
        config.write(resolved_config)
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
                "membership_semantics": "Bootstrap co-assignment may reject unstable WGCNA boundary events to grey, but never repartitions stable modules or uses phenotype, RBP or functional terms to create membership.",
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
                if suffix in {".tsv", ".csv"}:
                    separator = "\t" if suffix == ".tsv" else ","
                    frame = pd.read_csv(source, sep=separator, low_memory=False)
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
            resolved_config,
            entity_type="analysis_configuration",
            metadata={"mime_type": "application/yaml"},
        )
        if registered == 0:
            raise RuntimeError("CoSplice completed without any registered output files.")
        job.report(0.97, f"CoSplice registered {registered} engine outputs")
        return writer.finish()


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
