from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from katana_next.domain.models import AnalysisArtifact
from katana_next.infrastructure.jobs import JobContext
from katana_next.science.single_cell import load_single_cell_metadata
from katana_next.science.single_cell_benchmark import benchmark_single_cell_methods
from katana_next.science.single_cell_landscape import build_single_cell_landscape
from katana_next.science.single_cell_programs import (
    build_single_cell_programs,
    load_declared_cell_embedding,
)
from katana_next.services.project_service import ProjectContext
from katana_next.workflows.base import (
    ArtifactWriter,
    ParameterDefinition,
    WorkflowBase,
    WorkflowDefinition,
)
from katana_next.workflows.bulk import (
    ORANGE,
    RED,
    TEAL,
    _base_config,
    _latest_table,
    _plot,
)


class SingleCellLandscapeWorkflow(WorkflowBase):
    version = "2.1.1"
    definition = WorkflowDefinition(
        workflow_id="single_cell.landscape",
        workspace_id="single-cell",
        label="Splicing landscape and similarity",
        question="Do donor-level splicing profiles separate by cell type or condition, and which units are similar?",
        method="Posterior PSI quality gates, logit-scale PCA/UMAP and pairwise Spearman similarity on donor-by-cell-type pseudobulks.",
        output_summary="Event measurement QC, selected PSI matrix, donor-level PCA/UMAP, clustered similarity and event heatmap.",
        analysis_unit="donor by cell type pseudobulk",
        expected_artifacts=8,
        parameters=(
            ParameterDefinition(
                "min_total_reads",
                "Minimum donor-event reads",
                "number",
                10,
                "An event is measured in a donor-cell-type unit only above this count.",
                minimum=1,
                step=1,
            ),
            ParameterDefinition(
                "min_unit_fraction",
                "Minimum measured-unit fraction",
                "number",
                0.70,
                "Fraction of donor-cell-type units required for landscape inclusion.",
                minimum=0.1,
                maximum=1,
                step=0.05,
            ),
            ParameterDefinition(
                "min_psi_variance",
                "Minimum PSI variance",
                "number",
                0.0025,
                "Removes invariant events before embedding.",
                minimum=0,
                maximum=0.25,
                step=0.0005,
            ),
            ParameterDefinition(
                "max_events",
                "Maximum landscape events",
                "number",
                500,
                "Explicit variance-and-reliability-ranked cap; never applied silently.",
                minimum=20,
                maximum=5000,
                step=20,
            ),
            ParameterDefinition(
                "heatmap_events",
                "Heatmap events",
                "number",
                60,
                "Most informative selected events shown in the readable heatmap.",
                minimum=10,
                maximum=200,
                step=10,
            ),
            ParameterDefinition(
                "random_seed",
                "Random seed",
                "number",
                17,
                "Fixed UMAP seed for reproducibility.",
                minimum=0,
                step=1,
            ),
        ),
        limitations="Embeddings and correlations are descriptive. Donor-cell-type combinations are not counted as additional independent donors for significance testing.",
    )

    def parent_artifacts(
        self, context: ProjectContext, parameters: dict[str, Any]
    ) -> list[AnalysisArtifact]:
        artifact = _latest_table(
            context, "single_cell.quantify", "Donor-cell-type pseudobulk", parameters
        )
        return [artifact] if artifact else []

    def input_paths(
        self, context: ProjectContext, parameters: dict[str, Any]
    ) -> list[Path]:
        return [
            artifact.path for artifact in self.parent_artifacts(context, parameters)
        ]

    def run(
        self,
        context: ProjectContext,
        job: JobContext,
        run_fingerprint: str,
        parameters: dict[str, Any],
    ) -> list[str]:
        parent = _latest_table(
            context, "single_cell.quantify", "Donor-cell-type pseudobulk", parameters
        )
        if parent is None:
            raise ValueError(
                "Run Junction usage and donor pseudobulk first. Landscape never starts it implicitly."
            )
        job.report(0.12, "Building the measured donor-by-cell-type PSI matrix")
        result = build_single_cell_landscape(
            pd.read_parquet(parent.path),
            min_total_reads=float(parameters.get("min_total_reads", 10)),
            min_unit_fraction=float(parameters.get("min_unit_fraction", 0.70)),
            min_psi_variance=float(parameters.get("min_psi_variance", 0.0025)),
            max_events=int(parameters.get("max_events", 500)),
            heatmap_events=int(parameters.get("heatmap_events", 60)),
            random_seed=int(parameters.get("random_seed", 17)),
        )
        common = {
            "question": self.definition.question,
            "method": self.definition.method,
            "analysis_unit": self.definition.analysis_unit,
            "interpretation": "Nearby units have similar measured splicing profiles; separation is a pattern, not a population-level test.",
            "limitations": self.definition.limitations,
        }
        writer = ArtifactWriter(
            context,
            self.definition.workflow_id,
            run_fingerprint,
            parameters,
            self.version,
            common_metadata=common,
        )
        qc = writer.table(
            "Single-cell landscape event QC",
            "single_cell_landscape_event_qc.parquet",
            result.event_qc,
            entity_type="splice_event",
            parents=(parent.artifact_id,),
        )
        matrix = writer.table(
            "Selected donor PSI matrix",
            "single_cell_selected_donor_psi.parquet",
            result.selected_psi,
            entity_type="donor_cell_type_event",
            parents=(qc.artifact_id,),
        )
        embedding = writer.table(
            "Donor splicing embedding",
            "single_cell_donor_embedding.parquet",
            result.embedding,
            entity_type="donor_cell_type",
            parents=(matrix.artifact_id,),
        )
        similarity = writer.table(
            "Donor splicing similarity",
            "single_cell_donor_similarity.parquet",
            result.similarity,
            entity_type="donor_pair",
            parents=(matrix.artifact_id,),
        )
        writer.table(
            "Single-cell landscape diagnostics",
            "single_cell_landscape_diagnostics.parquet",
            result.diagnostics,
            entity_type="analysis_diagnostic",
            parents=(matrix.artifact_id,),
        )
        writer.plot(
            "Donor splicing PCA and UMAP",
            "single_cell_donor_embedding.vl.json",
            _embedding_plot(result.embedding),
            entity_type="donor_cell_type",
            parents=(embedding.artifact_id,),
        )
        writer.plot(
            "Donor splicing similarity heatmap",
            "single_cell_donor_similarity.vl.json",
            _similarity_plot(result.similarity),
            entity_type="donor_pair",
            parents=(similarity.artifact_id,),
        )
        writer.plot(
            "Variable splicing event heatmap",
            "single_cell_event_heatmap.vl.json",
            _event_heatmap(result.heatmap),
            entity_type="donor_cell_type_event",
            parents=(matrix.artifact_id,),
        )
        job.report(0.96, "Single-cell landscape and similarity complete")
        return writer.finish()


class SingleCellProgramsWorkflow(WorkflowBase):
    version = "2.2.1"
    definition = WorkflowDefinition(
        workflow_id="single_cell.programs",
        workspace_id="single-cell",
        label="Cell-state and replicate co-splicing programs",
        question="Which splice events follow a declared cell-state manifold, and which donor-level programs are reproducible?",
        method="Exploratory graph autocorrelation with donor-stratified permutation, plus the shared official R WGCNA signed TOM adapter for donor programs and donor-block bootstrap confidence.",
        output_summary="State-associated events, program catalog, hub membership, cell/donor scores, network edges and confidence diagnostics.",
        analysis_unit="cell for exploratory state association; donor-condition within cell type for replicate-aware discovery",
        expected_artifacts=11,
        parameters=(
            ParameterDefinition(
                "min_cell_event_reads",
                "Minimum cell-event reads",
                "number",
                1,
                "Observed cell-event gate for state association.",
                minimum=1,
                step=1,
            ),
            ParameterDefinition(
                "min_cell_fraction",
                "Minimum measured-cell fraction",
                "number",
                0.08,
                "Sparse-data coverage gate for state association.",
                minimum=0.01,
                maximum=1,
                step=0.01,
            ),
            ParameterDefinition(
                "min_state_cells",
                "Minimum measured cells",
                "number",
                40,
                "Minimum embedded cells per tested event.",
                minimum=20,
                step=10,
            ),
            ParameterDefinition(
                "state_neighbors",
                "State-graph neighbors",
                "number",
                15,
                "kNN graph size on the declared expression/state embedding.",
                minimum=5,
                maximum=100,
                step=5,
            ),
            ParameterDefinition(
                "state_permutations",
                "Donor-stratified permutations",
                "number",
                100,
                "Permutations preserve donor composition while testing local state association.",
                minimum=20,
                maximum=2000,
                step=20,
            ),
            ParameterDefinition(
                "state_fdr",
                "State-association FDR",
                "number",
                0.10,
                "BH FDR across tested events.",
                minimum=0,
                maximum=1,
                step=0.01,
            ),
            ParameterDefinition(
                "min_independent_donors",
                "Minimum donors for replicate programs",
                "number",
                15,
                "Hard gate for WGCNA-style programs within each cell type.",
                minimum=6,
                step=1,
            ),
            ParameterDefinition(
                "min_unit_reads",
                "Minimum donor-event reads",
                "number",
                10,
                "Read gate for donor program discovery.",
                minimum=1,
                step=1,
            ),
            ParameterDefinition(
                "min_unit_fraction",
                "Minimum donor coverage",
                "number",
                0.70,
                "Required donor fraction per program event.",
                minimum=0.1,
                maximum=1,
                step=0.05,
            ),
            ParameterDefinition(
                "min_psi_variance",
                "Minimum PSI variance",
                "number",
                0.0025,
                "Removes invariant donor-level events.",
                minimum=0,
                maximum=0.25,
                step=0.0005,
            ),
            ParameterDefinition(
                "min_module_size",
                "Minimum module size",
                "number",
                5,
                "Small modules are returned to grey rather than forced.",
                minimum=3,
                step=1,
            ),
            ParameterDefinition(
                "max_events",
                "Maximum program events (0 = all eligible)",
                "number",
                0,
                "Zero retains all eligible events. A positive cap is exploratory; dense-network memory is checked before fitting.",
                minimum=0,
                step=20,
            ),
            ParameterDefinition(
                "bootstrap_iterations",
                "Donor bootstrap iterations",
                "number",
                20,
                "Condition-stratified bootstrap confidence for donor programs.",
                minimum=5,
                maximum=500,
                step=5,
            ),
            ParameterDefinition(
                "soft_power",
                "Signed soft power",
                "number",
                6,
                "WGCNA-style signed adjacency exponent.",
                minimum=1,
                maximum=20,
                step=1,
            ),
            ParameterDefinition(
                "random_seed",
                "Random seed",
                "number",
                17,
                "Reproducible graph permutation and bootstrap seed.",
                minimum=0,
                step=1,
            ),
        ),
        optional_asset_kinds=("single_cell_embedding",),
        limitations="Cell-state modules are exploratory and depend on an independently declared expression/state embedding. Confirmatory modules require enough independent donors and never use cells as pseudo-replicates.",
    )

    def parent_artifacts(
        self, context: ProjectContext, parameters: dict[str, Any]
    ) -> list[AnalysisArtifact]:
        names = ("Single-cell junction choices", "Donor-cell-type pseudobulk")
        return [
            artifact
            for name in names
            if (artifact := _latest_table(context, "single_cell.quantify", name, parameters
        ))
        ]

    def input_paths(
        self, context: ProjectContext, parameters: dict[str, Any]
    ) -> list[Path]:
        paths = [
            artifact.path for artifact in self.parent_artifacts(context, parameters)
        ]
        paths.extend(super().input_paths(context, parameters))
        return list(dict.fromkeys(path.resolve() for path in paths))

    def run(
        self,
        context: ProjectContext,
        job: JobContext,
        run_fingerprint: str,
        parameters: dict[str, Any],
    ) -> list[str]:
        raw = _latest_table(
            context, "single_cell.quantify", "Single-cell junction choices", parameters
        )
        pseudo = _latest_table(
            context, "single_cell.quantify", "Donor-cell-type pseudobulk", parameters
        )
        if raw is None or pseudo is None:
            raise ValueError(
                "Run Junction usage and donor pseudobulk first. Programs never start it implicitly."
            )
        metadata = load_single_cell_metadata(context)
        embedding, embedding_source = load_declared_cell_embedding(context, metadata)
        job.report(
            0.10, "Testing cell-state association with donor-stratified permutations"
        )
        result = build_single_cell_programs(
            metadata,
            pd.read_parquet(raw.path),
            pd.read_parquet(pseudo.path),
            embedding,
            embedding_source=embedding_source,
            min_cell_event_reads=float(parameters.get("min_cell_event_reads", 1)),
            min_cell_fraction=float(parameters.get("min_cell_fraction", 0.08)),
            min_state_cells=int(parameters.get("min_state_cells", 40)),
            state_neighbors=int(parameters.get("state_neighbors", 15)),
            state_permutations=int(parameters.get("state_permutations", 100)),
            state_fdr=float(parameters.get("state_fdr", 0.10)),
            min_independent_donors=int(parameters.get("min_independent_donors", 15)),
            min_unit_reads=float(parameters.get("min_unit_reads", 10)),
            min_unit_fraction=float(parameters.get("min_unit_fraction", 0.70)),
            min_psi_variance=float(parameters.get("min_psi_variance", 0.0025)),
            min_module_size=int(parameters.get("min_module_size", 5)),
            max_events=int(parameters.get("max_events", 0)),
            bootstrap_iterations=int(parameters.get("bootstrap_iterations", 20)),
            soft_power=int(parameters.get("soft_power", 6)),
            random_seed=int(parameters.get("random_seed", 17)),
        )
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
                "interpretation": "State association and donor program discovery are separate evidence tracks and are never collapsed into one p-value.",
                "limitations": self.definition.limitations,
            },
        )
        diagnostics = writer.table(
            "Single-cell program diagnostics",
            "single_cell_program_diagnostics.parquet",
            result.diagnostics,
            entity_type="analysis_diagnostic",
            parents=(raw.artifact_id, pseudo.artifact_id),
        )
        parents = (diagnostics.artifact_id,)
        catalog = writer.table(
            "Single-cell program event catalog",
            "single_cell_program_events.parquet",
            result.event_catalog,
            entity_type="cell_type_splice_event",
            parents=parents,
        )
        writer.plot(
            "Cell-state association map",
            "single_cell_state_association.vl.json",
            _state_association_plot(result.event_catalog),
            entity_type="cell_type_splice_event",
            parents=(catalog.artifact_id,),
        )
        writer.plot(
            "Program event membership",
            "single_cell_program_membership.vl.json",
            _program_membership_plot(result.event_catalog),
            entity_type="cell_type_splice_event",
            parents=(catalog.artifact_id,),
        )
        summary = writer.table(
            "Single-cell program summary",
            "single_cell_program_summary.parquet",
            result.program_summary,
            entity_type="splicing_program",
            parents=parents,
        )
        writer.plot(
            "Program condition associations",
            "single_cell_program_condition_associations.vl.json",
            _program_association_plot(result.program_summary),
            entity_type="splicing_program",
            parents=(summary.artifact_id,),
        )
        cell_score = writer.table(
            "Cell-state program scores",
            "single_cell_state_program_scores.parquet",
            result.cell_scores,
            entity_type="cell_splicing_program",
            parents=parents,
        )
        writer.plot(
            "Cell manifold by splicing program",
            "single_cell_state_program_embedding.vl.json",
            _state_program_embedding(result.cell_scores),
            entity_type="cell_splicing_program",
            parents=(cell_score.artifact_id,),
        )
        donor_score = writer.table(
            "Donor program scores",
            "single_cell_donor_program_scores.parquet",
            result.donor_scores,
            entity_type="donor_splicing_program",
            parents=parents,
        )
        writer.plot(
            "Donor program score profiles",
            "single_cell_donor_program_scores.vl.json",
            _donor_program_plot(result.donor_scores),
            entity_type="donor_splicing_program",
            parents=(donor_score.artifact_id,),
        )
        writer.table(
            "Single-cell program network edges",
            "single_cell_program_edges.parquet",
            result.network_edges,
            entity_type="splice_event_pair",
            parents=parents,
        )
        job.report(0.96, "Cell-state and replicate co-splicing programs complete")
        return writer.finish()


class SingleCellBenchmarkWorkflow(WorkflowBase):
    version = "2.2.2"
    definition = WorkflowDefinition(
        workflow_id="single_cell.benchmark",
        workspace_id="single-cell",
        label="Single-cell method benchmark",
        question="How do quantification and program results compare on exactly shared cells, events and donors?",
        method="Reference-PSI RMSE/Spearman and cell-type recovery, plus matched bootstrap comparison of uncertainty-weighted TOM, unweighted TOM and official R WGCNA when available.",
        output_summary="Method availability, accuracy metrics, module assignments, partition metrics and matched-bootstrap stability.",
        analysis_unit="exact shared cell-event pairs and exact shared donor-event matrices",
        expected_artifacts=8,
        parameters=(
            ParameterDefinition(
                "min_cell_event_reads",
                "Minimum cell-event reads",
                "number",
                1,
                "Observed pair gate for raw/posterior accuracy.",
                minimum=1,
                step=1,
            ),
            ParameterDefinition(
                "min_donors",
                "Minimum donors",
                "number",
                15,
                "Program head-to-head gate within each cell type.",
                minimum=6,
                step=1,
            ),
            ParameterDefinition(
                "min_unit_reads",
                "Minimum donor-event reads",
                "number",
                10,
                "Common event-universe read gate.",
                minimum=1,
                step=1,
            ),
            ParameterDefinition(
                "min_unit_fraction",
                "Minimum donor coverage",
                "number",
                0.70,
                "Common event-universe coverage gate.",
                minimum=0.1,
                maximum=1,
                step=0.05,
            ),
            ParameterDefinition(
                "min_psi_variance",
                "Minimum PSI variance",
                "number",
                0.0025,
                "Common event-universe variance gate.",
                minimum=0,
                maximum=0.25,
                step=0.0005,
            ),
            ParameterDefinition(
                "min_module_size",
                "Minimum module size",
                "number",
                5,
                "Identical module-size gate across methods.",
                minimum=3,
                step=1,
            ),
            ParameterDefinition(
                "max_events",
                "Maximum benchmark events",
                "number",
                300,
                "Explicit identical event cap across methods.",
                minimum=20,
                maximum=3000,
                step=20,
            ),
            ParameterDefinition(
                "neighbors",
                "Embedding neighbors",
                "number",
                15,
                "kNN size for the external-manifold estimator.",
                minimum=5,
                maximum=100,
                step=5,
            ),
            ParameterDefinition(
                "bootstrap_iterations",
                "Matched bootstrap iterations",
                "number",
                10,
                "The same donor resamples are reused for every method.",
                minimum=3,
                maximum=100,
                step=1,
            ),
            ParameterDefinition(
                "soft_power",
                "Signed soft power",
                "number",
                6,
                "Shared adjacency power for the weighted and unweighted official R TOM adapter.",
                minimum=1,
                maximum=20,
                step=1,
            ),
            ParameterDefinition(
                "random_seed",
                "Random seed",
                "number",
                17,
                "Shared seed for clustering and matched bootstrap.",
                minimum=0,
                step=1,
            ),
        ),
        optional_asset_kinds=(
            "single_cell_embedding",
            "single_cell_reference_psi",
            "single_cell_reference_modules",
            "single_cell_competitor_psi",
        ),
        limitations="Imported competitors are scored but not claimed as executed by Katana. Without an external reference, PSI accuracy cannot be ranked. A fixture benchmark validates software, not publication-level superiority.",
    )

    def parent_artifacts(
        self, context: ProjectContext, parameters: dict[str, Any]
    ) -> list[AnalysisArtifact]:
        names = ("Single-cell junction choices", "Donor-cell-type pseudobulk")
        return [
            artifact
            for name in names
            if (artifact := _latest_table(context, "single_cell.quantify", name, parameters
        ))
        ]

    def input_paths(
        self, context: ProjectContext, parameters: dict[str, Any]
    ) -> list[Path]:
        paths = [
            artifact.path for artifact in self.parent_artifacts(context, parameters)
        ]
        paths.extend(super().input_paths(context, parameters))
        return list(dict.fromkeys(path.resolve() for path in paths))

    def run(
        self,
        context: ProjectContext,
        job: JobContext,
        run_fingerprint: str,
        parameters: dict[str, Any],
    ) -> list[str]:
        raw = _latest_table(
            context, "single_cell.quantify", "Single-cell junction choices", parameters
        )
        pseudo = _latest_table(
            context, "single_cell.quantify", "Donor-cell-type pseudobulk", parameters
        )
        if raw is None or pseudo is None:
            raise ValueError(
                "Run Junction usage and donor pseudobulk first. Benchmark never starts it implicitly."
            )
        metadata = load_single_cell_metadata(context)
        embedding, _source = load_declared_cell_embedding(context, metadata)
        reference = _optional_asset_table(context, "single_cell_reference_psi")
        reference_modules = _optional_asset_table(
            context, "single_cell_reference_modules"
        )
        competitor = _optional_asset_table(context, "single_cell_competitor_psi")
        job.report(0.08, "Benchmarking PSI estimators on exact shared cell-event pairs")
        result = benchmark_single_cell_methods(
            metadata,
            pd.read_parquet(raw.path),
            pd.read_parquet(pseudo.path),
            embedding=embedding,
            reference_psi=reference,
            reference_modules=reference_modules,
            competitor_psi=competitor,
            min_cell_event_reads=float(parameters.get("min_cell_event_reads", 1)),
            min_donors=int(parameters.get("min_donors", 15)),
            min_unit_reads=float(parameters.get("min_unit_reads", 10)),
            min_unit_fraction=float(parameters.get("min_unit_fraction", 0.70)),
            min_psi_variance=float(parameters.get("min_psi_variance", 0.0025)),
            min_module_size=int(parameters.get("min_module_size", 5)),
            max_events=int(parameters.get("max_events", 300)),
            neighbors=int(parameters.get("neighbors", 15)),
            bootstrap_iterations=int(parameters.get("bootstrap_iterations", 10)),
            soft_power=int(parameters.get("soft_power", 6)),
            random_seed=int(parameters.get("random_seed", 17)),
        )
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
                "interpretation": "Read each metric separately; Katana does not manufacture a composite winner.",
                "limitations": self.definition.limitations,
            },
        )
        audit = writer.table(
            "Single-cell benchmark method audit",
            "single_cell_benchmark_method_audit.parquet",
            result.method_audit,
            entity_type="method_audit",
            parents=(raw.artifact_id, pseudo.artifact_id),
        )
        quant = writer.table(
            "Single-cell PSI benchmark",
            "single_cell_psi_benchmark.parquet",
            result.quantification_summary,
            entity_type="method_performance",
            parents=(audit.artifact_id,),
        )
        writer.plot(
            "Single-cell PSI accuracy",
            "single_cell_psi_accuracy.vl.json",
            _psi_benchmark_plot(result.quantification_summary),
            entity_type="method_performance",
            parents=(quant.artifact_id,),
        )
        program = writer.table(
            "Single-cell program benchmark",
            "single_cell_program_benchmark.parquet",
            result.program_summary,
            entity_type="method_performance",
            parents=(audit.artifact_id,),
        )
        writer.plot(
            "Single-cell program head-to-head",
            "single_cell_program_benchmark.vl.json",
            _program_benchmark_plot(result.program_summary),
            entity_type="method_performance",
            parents=(program.artifact_id,),
        )
        writer.table(
            "Single-cell benchmark module assignments",
            "single_cell_benchmark_assignments.parquet",
            result.program_assignments,
            entity_type="cell_type_splice_event",
            parents=(audit.artifact_id,),
        )
        bootstrap = writer.table(
            "Single-cell matched bootstrap stability",
            "single_cell_benchmark_bootstrap.parquet",
            result.bootstrap_stability,
            entity_type="method_bootstrap",
            parents=(audit.artifact_id,),
        )
        writer.plot(
            "Matched donor-bootstrap stability",
            "single_cell_bootstrap_stability.vl.json",
            _bootstrap_plot(result.bootstrap_stability),
            entity_type="method_bootstrap",
            parents=(bootstrap.artifact_id,),
        )
        job.report(0.96, "Single-cell method benchmark complete")
        return writer.finish()


def _optional_asset_table(context: ProjectContext, kind: str) -> pd.DataFrame | None:
    asset = next(
        (
            item
            for item in context.manifest.assets
            if item.kind == kind and item.path.is_file()
        ),
        None,
    )
    if asset is None:
        return None
    return pd.read_csv(asset.path, sep=None, engine="python")


def _empty_plot(title: str, message: str) -> dict[str, Any]:
    return {
        "$schema": "https://vega.github.io/schema/vega-lite/v5.json",
        "description": message,
        "autosize": {"type": "fit", "contains": "padding"},
        "width": "container",
        "height": 240,
        "title": title,
        "data": {"values": [{"message": message}]},
        "mark": {"type": "text", "fontSize": 16, "color": "#64748b"},
        "encoding": {"text": {"field": "message"}},
        "config": _base_config(),
    }


def _embedding_plot(frame: pd.DataFrame) -> dict[str, Any]:
    if frame.empty or not {"pca_1", "pca_2", "umap_1", "umap_2"}.issubset(frame):
        return _empty_plot(
            "Donor splicing PCA and UMAP",
            "No eligible donor-level PSI matrix was available.",
        )
    base = frame.to_dict("records")
    panels = []
    for prefix, title in (("pca", "PCA"), ("umap", "UMAP")):
        panels.append(
            {
                "title": title,
                "data": {"values": base},
                "mark": {
                    "type": "point",
                    "filled": True,
                    "size": 105,
                    "opacity": 0.86,
                    "stroke": "white",
                },
                "encoding": {
                    "x": {
                        "field": f"{prefix}_1",
                        "type": "quantitative",
                        "axis": {"title": f"{title} 1"},
                    },
                    "y": {
                        "field": f"{prefix}_2",
                        "type": "quantitative",
                        "axis": {"title": f"{title} 2"},
                    },
                    "color": {
                        "field": "condition",
                        "type": "nominal",
                        "scale": {"scheme": "tableau10"},
                    },
                    "shape": {"field": "cell_type", "type": "nominal"},
                    "tooltip": [
                        {"field": "donor_id"},
                        {"field": "cell_type"},
                        {"field": "condition"},
                        {"field": "embedding_method"},
                    ],
                },
            }
        )
    return {
        "$schema": "https://vega.github.io/schema/vega-lite/v5.json",
        "description": "Donor-cell-type splicing embedding.",
        "autosize": {"type": "fit", "contains": "padding"},
        "hconcat": panels,
        "config": _base_config(),
    }


def _similarity_plot(frame: pd.DataFrame) -> dict[str, Any]:
    if frame.empty or "spearman_similarity" not in frame:
        return _empty_plot(
            "Donor splicing similarity",
            "No eligible donor-level similarity matrix was available.",
        )
    return _plot(
        frame.to_dict("records"),
        {"type": "rect"},
        {
            "x": {
                "field": "unit_b",
                "type": "nominal",
                "sort": {"field": "unit_b_order", "op": "min"},
                "axis": {"title": None, "labelAngle": -90, "labelLimit": 130},
            },
            "y": {
                "field": "unit_a",
                "type": "nominal",
                "sort": {"field": "unit_a_order", "op": "min"},
                "axis": {"title": None, "labelLimit": 150},
            },
            "color": {
                "field": "spearman_similarity",
                "type": "quantitative",
                "scale": {"scheme": "redblue", "domainMid": 0},
                "legend": {"title": "Spearman r"},
            },
            "tooltip": [
                {"field": "unit_a"},
                {"field": "unit_b"},
                {"field": "spearman_similarity", "format": ".2f"},
            ],
        },
        title="Similarity of measured donor-level splicing profiles",
        description="Pairwise Spearman correlation on selected posterior PSI events.",
        height=520,
    )


def _event_heatmap(frame: pd.DataFrame) -> dict[str, Any]:
    if frame.empty or "standardized_psi" not in frame:
        return _empty_plot(
            "Variable splicing events", "No events passed the declared landscape gates."
        )
    return _plot(
        frame.to_dict("records"),
        {"type": "rect"},
        {
            "x": {
                "field": "unit_id",
                "type": "nominal",
                "axis": {
                    "title": "Donor :: cell type",
                    "labelAngle": -90,
                    "labelLimit": 120,
                },
            },
            "y": {
                "field": "event_label",
                "type": "nominal",
                "axis": {"title": None, "labelLimit": 150},
            },
            "color": {
                "field": "standardized_psi",
                "type": "quantitative",
                "scale": {"scheme": "redblue", "domainMid": 0},
                "legend": {"title": "Standardized PSI"},
            },
            "tooltip": [
                {"field": "event_id"},
                {"field": "gene_symbol"},
                {"field": "donor_id"},
                {"field": "cell_type"},
                {"field": "condition"},
                {"field": "standardized_psi", "format": ".2f"},
            ],
        },
        title="Most informative measured splicing events",
        description="Rows are reliability-ranked events; colors are event-wise standardized posterior PSI.",
        height=max(420, min(900, 13 * frame["event_id"].nunique())),
    )


def _state_association_plot(frame: pd.DataFrame) -> dict[str, Any]:
    if frame.empty or "program_track" not in frame:
        return _empty_plot(
            "Cell-state association", "No cell-state association result was available."
        )
    values = frame[frame["program_track"].eq("cell_state_exploratory")].copy()
    values["minus_log10_fdr"] = -np.log10(
        pd.to_numeric(values["fdr"], errors="coerce").clip(lower=1e-300)
    )
    return _plot(
        values.to_dict("records"),
        {
            "type": "point",
            "filled": True,
            "size": 75,
            "opacity": 0.82,
            "stroke": "white",
        },
        {
            "x": {
                "field": "state_autocorrelation",
                "type": "quantitative",
                "axis": {"title": "State-graph autocorrelation"},
            },
            "y": {
                "field": "minus_log10_fdr",
                "type": "quantitative",
                "axis": {"title": "-log10 FDR"},
            },
            "color": {
                "field": "state_associated",
                "type": "nominal",
                "scale": {"domain": [True, False], "range": [RED, "#94a3b8"]},
                "legend": {"title": "State associated"},
            },
            "shape": {"field": "module", "type": "nominal"},
            "tooltip": [
                {"field": "gene_symbol"},
                {"field": "event_id"},
                {"field": "state_autocorrelation", "format": ".3f"},
                {"field": "fdr", "format": ".2e"},
                {"field": "measured_cells"},
                {"field": "measured_donors"},
                {"field": "module"},
            ],
        },
        title="Cell-state-associated splicing events",
        description="P-values come from donor-stratified permutation on a declared external manifold.",
        height=430,
    )


def _program_membership_plot(frame: pd.DataFrame) -> dict[str, Any]:
    if frame.empty or not {"module", "module_kme"}.issubset(frame):
        return _empty_plot(
            "Program event membership", "No retained splicing program was available."
        )
    values = frame[frame["module"].ne("grey")].copy()
    values["event_label"] = values["gene_symbol"].where(
        values["gene_symbol"].ne(""), values["event_id"]
    )
    values = (
        values.sort_values(
            "module_kme", key=lambda column: column.abs(), ascending=False
        )
        .groupby(["program_track", "cell_type_scope", "module"])
        .head(20)
    )
    return _plot(
        values.to_dict("records"),
        {"type": "bar", "cornerRadiusEnd": 3},
        {
            "x": {
                "field": "module_kme",
                "type": "quantitative",
                "axis": {"title": "Module membership (kME)"},
            },
            "y": {
                "field": "event_label",
                "type": "nominal",
                "sort": "-x",
                "axis": {"title": None, "labelLimit": 150},
            },
            "color": {
                "field": "module",
                "type": "nominal",
                "scale": {"scheme": "tableau10"},
            },
            "facet": {"field": "cell_type_scope", "type": "nominal", "columns": 2},
            "tooltip": [
                {"field": "program_track"},
                {"field": "cell_type_scope"},
                {"field": "module"},
                {"field": "event_id"},
                {"field": "module_kme", "format": ".3f"},
                {"field": "bootstrap_membership_confidence", "format": ".2f"},
            ],
        },
        title="Hub events and module membership",
        description="Top absolute-kME events per program; grey events are excluded.",
        height=520,
    )


def _program_association_plot(frame: pd.DataFrame) -> dict[str, Any]:
    required = {
        "program_track",
        "cell_type_scope",
        "module",
        "max_condition_mean_difference",
        "condition_omnibus_fdr",
    }
    if frame.empty or not required.issubset(frame):
        return _empty_plot(
            "Program-condition associations",
            "No donor-level program association test was available.",
        )
    values = frame.copy()
    values["minus_log10_fdr"] = -np.log10(
        pd.to_numeric(values["condition_omnibus_fdr"], errors="coerce").clip(
            lower=1e-300
        )
    )
    values["program_label"] = (
        values["cell_type_scope"].astype(str) + " :: " + values["module"].astype(str)
    )
    return _plot(
        values.to_dict("records"),
        {"type": "point", "filled": True, "size": 135, "stroke": "white"},
        {
            "x": {
                "field": "max_condition_mean_difference",
                "type": "quantitative",
                "axis": {"title": "Largest donor-level condition mean difference"},
            },
            "y": {
                "field": "minus_log10_fdr",
                "type": "quantitative",
                "axis": {"title": "-log10 omnibus FDR"},
            },
            "color": {
                "field": "program_track",
                "type": "nominal",
                "scale": {"range": [TEAL, ORANGE]},
            },
            "shape": {"field": "cell_type_scope", "type": "nominal"},
            "tooltip": [
                {"field": "program_label"},
                {"field": "program_track"},
                {"field": "highest_score_condition"},
                {"field": "lowest_score_condition"},
                {"field": "max_condition_mean_difference", "format": ".3f"},
                {"field": "condition_omnibus_fdr", "format": ".2e"},
                {"field": "condition_test_status"},
                {"field": "condition_test_batch_adjusted"},
            ],
        },
        title="Condition association of donor-level splicing programs",
        description="Nested-model F test compares batch-only with batch-plus-condition models; each point is one retained program.",
        height=420,
    )


def _state_program_embedding(frame: pd.DataFrame) -> dict[str, Any]:
    if frame.empty or not {"module", "state_1", "state_2"}.issubset(frame):
        return _empty_plot(
            "Cell-state program activity",
            "No exploratory cell-state program was retained.",
        )
    modules = frame["module"].drop_duplicates().head(6)
    values = frame[frame["module"].isin(modules)]
    return _plot(
        values.to_dict("records"),
        {"type": "circle", "size": 34, "opacity": 0.78},
        {
            "x": {
                "field": "state_1",
                "type": "quantitative",
                "axis": {"title": "Declared state 1"},
            },
            "y": {
                "field": "state_2",
                "type": "quantitative",
                "axis": {"title": "Declared state 2"},
            },
            "color": {
                "field": "program_score",
                "type": "quantitative",
                "scale": {"scheme": "redblue", "domainMid": 0},
                "legend": {"title": "Program score"},
            },
            "facet": {"field": "module", "type": "nominal", "columns": 3},
            "tooltip": [
                {"field": "cell_id"},
                {"field": "donor_id"},
                {"field": "cell_type"},
                {"field": "condition"},
                {"field": "module"},
                {"field": "program_score", "format": ".2f"},
            ],
        },
        title="Splicing program activity on the declared cell-state manifold",
        description="Exploratory cell score; donor remains the inferential unit.",
        height=440,
    )


def _donor_program_plot(frame: pd.DataFrame) -> dict[str, Any]:
    if frame.empty or "program_score" not in frame:
        return _empty_plot(
            "Donor program scores", "No donor-level splicing program was retained."
        )
    return _plot(
        frame.to_dict("records"),
        {
            "type": "point",
            "filled": True,
            "size": 72,
            "opacity": 0.82,
            "stroke": "white",
        },
        {
            "x": {
                "field": "condition",
                "type": "nominal",
                "axis": {"title": None, "labelAngle": -35},
            },
            "y": {
                "field": "program_score",
                "type": "quantitative",
                "axis": {"title": "Program score"},
            },
            "color": {
                "field": "condition",
                "type": "nominal",
                "scale": {"scheme": "tableau10"},
            },
            "shape": {"field": "program_track", "type": "nominal"},
            "column": {
                "field": "module",
                "type": "nominal",
                "header": {"title": "Program"},
            },
            "tooltip": [
                {"field": "program_track"},
                {"field": "cell_type_scope"},
                {"field": "module"},
                {"field": "donor_id"},
                {"field": "condition"},
                {"field": "program_score", "format": ".3f"},
            ],
        },
        title="Program activity across biological donors",
        description="Each point is one donor; cells do not increase the replicate count.",
        height=390,
    )


def _psi_benchmark_plot(frame: pd.DataFrame) -> dict[str, Any]:
    if frame.empty or not {"rmse", "spearman"}.issubset(frame):
        return _empty_plot(
            "Single-cell PSI benchmark",
            "Declare reference PSI to compute accuracy metrics.",
        )
    return _plot(
        frame.to_dict("records"),
        {"type": "point", "filled": True, "size": 150, "stroke": "white"},
        {
            "x": {
                "field": "rmse",
                "type": "quantitative",
                "axis": {"title": "PSI RMSE (lower is better)"},
            },
            "y": {
                "field": "spearman",
                "type": "quantitative",
                "axis": {"title": "PSI Spearman r (higher is better)"},
            },
            "size": {
                "field": "coverage_fraction",
                "type": "quantitative",
                "legend": {"title": "Reference coverage"},
            },
            "color": {
                "field": "method",
                "type": "nominal",
                "scale": {"scheme": "tableau10"},
            },
            "tooltip": [
                {"field": "method"},
                {"field": "rmse", "format": ".4f"},
                {"field": "spearman", "format": ".3f"},
                {"field": "coverage_fraction", "format": ".2%"},
                {"field": "cell_type_ari", "format": ".3f"},
                {"field": "cell_type_nmi", "format": ".3f"},
                {"field": "uses_external_embedding"},
            ],
        },
        title="Single-cell PSI recovery on exact shared pairs",
        description="No winner is declared when reference truth is absent.",
        height=420,
    )


def _program_benchmark_plot(frame: pd.DataFrame) -> dict[str, Any]:
    if frame.empty or not {"cell_type", "method"}.issubset(frame):
        return _empty_plot(
            "Single-cell program benchmark",
            "No cell type passed the independent-donor benchmark gate.",
        )
    metrics = [
        metric
        for metric in (
            "reference_ari",
            "silhouette_cosine",
            "mean_bootstrap_ari",
            "grey_fraction",
        )
        if metric in frame
    ]
    if not metrics:
        return _empty_plot(
            "Single-cell program benchmark",
            "No comparable program metrics were available.",
        )
    values = frame.melt(
        id_vars=["cell_type", "method"],
        value_vars=metrics,
        var_name="metric",
        value_name="value",
    )
    return _plot(
        values.to_dict("records"),
        {"type": "bar", "cornerRadiusEnd": 3},
        {
            "x": {
                "field": "method",
                "type": "nominal",
                "axis": {"title": None, "labelAngle": -35, "labelLimit": 190},
            },
            "y": {
                "field": "value",
                "type": "quantitative",
                "axis": {"title": "Metric value"},
            },
            "color": {
                "field": "method",
                "type": "nominal",
                "scale": {"scheme": "tableau10"},
                "legend": None,
            },
            "column": {"field": "metric", "type": "nominal", "header": {"title": None}},
            "row": {
                "field": "cell_type",
                "type": "nominal",
                "header": {"title": "Cell type"},
            },
            "tooltip": [
                {"field": "cell_type"},
                {"field": "method"},
                {"field": "metric"},
                {"field": "value", "format": ".3f"},
            ],
        },
        title="Matched program-method comparison",
        description="All methods use the same donors and events; metrics are reported separately.",
        height=260,
    )


def _bootstrap_plot(frame: pd.DataFrame) -> dict[str, Any]:
    if frame.empty or "status" not in frame:
        return _empty_plot(
            "Matched donor bootstrap",
            "No comparable donor-level partition was available.",
        )
    complete = frame[frame["status"].eq("complete")]
    return _plot(
        complete.to_dict("records"),
        {"type": "boxplot", "extent": "min-max", "size": 24},
        {
            "x": {
                "field": "method",
                "type": "nominal",
                "axis": {"title": None, "labelAngle": -35, "labelLimit": 190},
            },
            "y": {
                "field": "adjusted_rand_index",
                "type": "quantitative",
                "scale": {"domain": [0, 1]},
                "axis": {"title": "Bootstrap ARI"},
            },
            "color": {
                "field": "method",
                "type": "nominal",
                "scale": {"scheme": "tableau10"},
                "legend": None,
            },
            "column": {
                "field": "cell_type",
                "type": "nominal",
                "header": {"title": "Cell type"},
            },
            "tooltip": [
                {"field": "method"},
                {"field": "cell_type"},
                {"field": "iteration"},
                {"field": "adjusted_rand_index", "format": ".3f"},
            ],
        },
        title="Matched donor-bootstrap partition stability",
        description="The same stratified donor resamples are used for every available method.",
        height=360,
    )
