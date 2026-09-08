from __future__ import annotations

import json
from pathlib import Path
from zipfile import ZipFile

import numpy as np
import pandas as pd
import pytest
import katana_cosplice.cli as cosplice_cli
from sklearn.metrics import adjusted_rand_score

from katana_cosplice.audit import audit_inputs, load_sample_metadata
from katana_cosplice.association import factorial_program_association
from katana_cosplice.benchmark import run_method_benchmark
from katana_cosplice.config import (
    BenchmarkConfig,
    ConsensusConfig,
    CoSpliceConfig,
    ExternalCohortConfig,
    HybridProgramConfig,
    NetworkConfig,
    PreservationConfig,
    StableASConfig,
    ValidationConfig,
    ModuleInterpretationConfig,
    ModuleStatisticsConfig,
    MetaAnalysisConfig,
    FactorContrastConfig,
)
from katana_cosplice.consensus import (
    _backbone_confidence_only_membership,
    _integrate_kme_membership_evidence,
    build_bootstrap_coassignment_consensus,
)
from katana_cosplice.hybrid_programs import build_hybrid_programs
from katana_cosplice.matrix import EVENT_COORDINATES, build_sample_matrices
from katana_cosplice.module_annotation import annotate_modules
from katana_cosplice.module_interpretation import interpret_modules
from katana_cosplice.module_models import fit_replicate_aware_models
from katana_cosplice.meta_analysis import meta_analyze_modules
from katana_cosplice.factor_contrasts import analyze_factor_contrasts, render_factor_contrast_plots
from katana_cosplice.event_harmonization import annotation_aware_event_alignment
from katana_cosplice.isoform_switch import quantify_isoform_switches
from katana_cosplice.motif_scanner import _event_windows, load_motif_database, scan_event_motifs
from katana_cosplice.network import (
    _adaptive_membership_cutoffs,
    _effective_consensus_bootstrap_iterations,
    _sparse_knn_edges,
    fit_network,
    signed_soft_adjacency,
    topological_overlap_matrix,
)
from katana_cosplice.pipeline import CoSplicePipeline
from katana_cosplice.preservation import evaluate_external_module_preservation
from katana_cosplice.project import CoSpliceProject
from katana_cosplice.qc import stable_as_qc
from katana_cosplice.stage_registry import stage_spec, stage_specs
from katana_cosplice.validation import validate_network_modules


def _synthetic_project(tmp_path: Path) -> CoSpliceConfig:
    rmats = tmp_path / "rmats"
    (rmats / "bamlist").mkdir(parents=True)
    (rmats / "bamlist" / "A_bamlist.txt").write_text("A_1_sorted.bam,A_2_sorted.bam,A_3_sorted.bam", encoding="utf-8")
    (rmats / "bamlist" / "B_bamlist.txt").write_text("B_1_sorted.bam,B_2_sorted.bam,B_3_sorted.bam", encoding="utf-8")
    post = rmats / "A-B" / "rmats_post"
    post.mkdir(parents=True)
    for event_index, (event_type, coordinates) in enumerate(EVENT_COORDINATES.items(), start=1):
        row: dict[str, object] = {
            "GeneID": f'"ENSMUSG{event_index:011d}.1"',
            "geneSymbol": f'"Gene{event_index}"',
            "chr": "chr1",
            "strand": "+",
            "IncLevel1": "0.20,0.22,0.21",
            "IncLevel2": "0.70,0.72,0.71",
            "IJC_SAMPLE_1": "30,31,32",
            "SJC_SAMPLE_1": "50,50,49",
            "IJC_SAMPLE_2": "70,72,71",
            "SJC_SAMPLE_2": "20,20,21",
        }
        for offset, column in enumerate(coordinates):
            row[column] = 1000 * event_index + offset * 10
        pd.DataFrame([row]).to_csv(post / f"{event_type}.MATS.JCEC.txt", sep="\t", index=False)
    metadata = pd.DataFrame(
        {
            "sample": ["A_1", "A_2", "A_3", "B_1", "B_2", "B_3"],
            "condition": ["A"] * 3 + ["B"] * 3,
        }
    )
    metadata_path = tmp_path / "samples.csv"
    metadata.to_csv(metadata_path, index=False)
    return CoSpliceConfig("synthetic", rmats, metadata_path, tmp_path / "output")


def test_audit_and_matrix_preserve_sample_order(tmp_path: Path) -> None:
    config = _synthetic_project(tmp_path)
    audit = audit_inputs(config)
    assert audit.passed

    result = build_sample_matrices(config)
    assert result.psi.shape == (5, 6)
    assert result.psi.columns.tolist() == ["A_1", "A_2", "A_3", "B_1", "B_2", "B_3"]
    assert np.isclose(result.psi.loc[result.events.iloc[0]["event_key"], "A_1"], 0.20)
    assert result.events["event_key"].is_unique


def test_qc_and_network_build_modules() -> None:
    rng = np.random.default_rng(7)
    samples = pd.DataFrame(
        {
            "sample_id": [f"S{i:02d}" for i in range(18)],
            "condition": np.repeat(["A", "B", "C", "D", "E", "F"], 3),
        }
    )
    latent_a = np.repeat(np.linspace(0.2, 0.8, 6), 3)
    latent_b = np.repeat(np.linspace(0.8, 0.2, 6), 3)
    matrix = np.vstack(
        [latent_a + rng.normal(0, 0.015, 18) for _ in range(25)]
        + [latent_b + rng.normal(0, 0.015, 18) for _ in range(25)]
    )
    psi = pd.DataFrame(np.clip(matrix, 0.01, 0.99), index=[f"E{i}" for i in range(50)], columns=samples["sample_id"])
    ijc = pd.DataFrame(80.0, index=psi.index, columns=psi.columns)
    sjc = pd.DataFrame(50.0, index=psi.index, columns=psi.columns)
    qc = stable_as_qc(
        psi,
        ijc,
        sjc,
        samples,
        StableASConfig(min_repeatability=0.1, min_psi_variance=0.0001),
    )
    assert len(qc.stable_psi) == 50
    network = fit_network(
        qc.stable_psi,
        qc.metrics,
        NetworkConfig(min_module_size=5, module_cut_height=0.4, edge_threshold=0.01, bootstrap_iterations=5),
        samples=samples,
    )
    assert not network.assignments.empty
    assert not network.program_scores.empty
    assert network.assignments["event_key"].is_unique
    assert "bootstrap_support" in network.edges
    assert "graph_bootstrap_support" in network.edges
    assert not network.cut_height_sensitivity.empty
    assert len(network.selected_matrix) == 50
    assert network.diagnostics["selection_rule"].eq("all_qc_passed").all()
    assert {
        "module_before_refinement",
        "module_membership",
        "best_module_membership",
        "second_best_module_membership",
        "module_membership_margin",
        "module_refinement_action",
    }.issubset(network.assignments.columns)
    assert network.diagnostics["module_refinement"].all()

    declared = NetworkConfig(consensus_bootstrap_iterations=20)
    assert _effective_consensus_bootstrap_iterations(1_000, declared) == 20
    assert _effective_consensus_bootstrap_iterations(7_000, declared) == 20

    capped = fit_network(
        qc.stable_psi,
        qc.metrics,
        NetworkConfig(min_module_size=5, max_events=40, edge_threshold=0.01, bootstrap_iterations=0),
        samples=samples,
    )
    assert len(capped.selected_matrix) == 40
    assert capped.diagnostics["selection_rule"].eq("top_40_by_stability_score").all()

    validation = validate_network_modules(
        network.selected_matrix,
        network.assignments,
        samples,
        NetworkConfig(min_module_size=5, edge_threshold=0.01, bootstrap_iterations=0),
        ValidationConfig(
            bootstrap_iterations=2,
            knn_neighbors=[25, 50],
            community_resolutions=[1.0],
        ),
        residual_assignments=network.assignments,
    )
    assert len(validation.bootstrap_runs) == 2
    assert validation.event_stability["event_key"].is_unique
    assert not validation.module_stability.empty
    assert not validation.parameter_sensitivity.empty
    assert validation.preservation_summary["best_jaccard"].eq(1.0).all()


def test_qc_requires_condition_coverage_and_dynamic_psi_range() -> None:
    samples = pd.DataFrame(
        {
            "sample_id": ["A1", "A2", "A3", "B1", "B2", "B3"],
            "condition": ["A", "A", "A", "B", "B", "B"],
        }
    )
    psi = pd.DataFrame(
        [
            [0.20, 0.21, 0.19, 0.70, 0.71, 0.69],
            [0.48, 0.49, 0.50, 0.51, 0.50, 0.49],
            [0.20, np.nan, np.nan, 0.70, 0.71, 0.69],
        ],
        index=["dynamic", "flat", "undercovered"],
        columns=samples["sample_id"],
    )
    reads = pd.DataFrame(50.0, index=psi.index, columns=psi.columns).where(psi.notna())
    result = stable_as_qc(
        psi,
        reads,
        reads,
        samples,
        StableASConfig(min_repeatability=0.0, min_psi_variance=0.0),
    )
    assert bool(result.metrics.loc["dynamic", "stable_as"])
    assert not bool(result.metrics.loc["flat", "passes_dynamic_range"])
    assert not bool(result.metrics.loc["undercovered", "passes_condition_support"])
    assert len(result.condition_metrics) == int(result.metrics["stable_as"].sum()) * 2
    assert set(result.condition_metrics["event_key"]) == set(result.stable_psi.index)
    assert set(result.summary["rule"]) >= {"passes_read_support", "passes_repeatability", "stable_as"}


def test_metadata_rejects_duplicate_normalized_sample_ids(tmp_path: Path) -> None:
    path = tmp_path / "samples.tsv"
    pd.DataFrame(
        {
            "sample_id": ["A_1.bam", "A_1_sorted.bam"],
            "condition": ["A", "A"],
        }
    ).to_csv(path, sep="\t", index=False)
    with pytest.raises(ValueError, match="duplicate normalized sample_id"):
        load_sample_metadata(path)


def test_matrix_rejects_bamlist_condition_mismatch(tmp_path: Path) -> None:
    config = _synthetic_project(tmp_path)
    (config.rmats_root / "bamlist" / "A_bamlist.txt").write_text(
        "A_1_sorted.bam,A_2_sorted.bam,B_1_sorted.bam",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="does not match metadata condition A"):
        build_sample_matrices(config)


def test_qc_rejects_out_of_range_psi() -> None:
    samples = pd.DataFrame({"sample_id": ["A1", "A2", "B1", "B2"], "condition": ["A", "A", "B", "B"]})
    psi = pd.DataFrame([[0.2, 0.3, 1.2, 0.4]], index=["event"], columns=samples["sample_id"])
    reads = pd.DataFrame(20.0, index=psi.index, columns=psi.columns)
    with pytest.raises(ValueError, match=r"within \[0, 1\]"):
        stable_as_qc(psi, reads, reads, samples, StableASConfig())


def test_pipeline_does_not_run_missing_upstream_implicitly(tmp_path: Path) -> None:
    pipeline = CoSplicePipeline(_synthetic_project(tmp_path))
    with pytest.raises(FileNotFoundError, match="Run build-matrix first"):
        pipeline.stable_qc()
    assert not (tmp_path / "output" / "03_stable_as_qc" / "stable_psi_matrix.tsv").exists()
    with pytest.raises(FileNotFoundError, match="Observed network first"):
        pipeline.validate_modules()
    with pytest.raises(FileNotFoundError, match="network_observed first"):
        pipeline.consensus_modules()


def test_stage_registry_is_ordered_and_dependencies_are_explicit() -> None:
    specs = stage_specs(visible_only=True)
    assert [spec.order for spec in specs] == sorted(spec.order for spec in specs)
    assert len({spec.key for spec in specs}) == len(specs)
    assert stage_spec("stable_qc").dependencies == ("matrix",)
    assert stage_spec("associate_observed").dependencies == ("network_observed",)
    assert stage_spec("consensus_modules").dependencies == ("network_observed",)
    assert stage_spec("preserve_modules").dependencies == ("consensus_modules",)
    assert stage_spec("annotate_modules").dependencies == ("network_observed",)
    assert stage_spec("benchmark_methods").dependencies == ("stable_qc",)
    assert stage_spec("interpret_modules").dependencies == ("network_observed",)
    assert stage_spec("interpret_modules").optional_dependencies == ("consensus_modules",)
    assert stage_spec("fit_module_models").dependencies == ("matrix", "network_observed")
    assert stage_spec("fit_module_models").optional_dependencies == ("consensus_modules",)
    assert stage_spec("meta_analyze_modules").dependencies == ("consensus_modules",)
    assert stage_spec("factor_contrasts").dependencies == ("stable_qc",)
    assert all(spec.explicit_only for spec in specs)


def test_configured_factor_contrasts_are_directional_and_rendered(tmp_path: Path) -> None:
    rows: list[dict[str, str]] = []
    for context in ("noDP", "DP"):
        for genotype in ("WT", "MUT"):
            for treatment in ("Baseline", "CLK3"):
                for replicate in range(1, 4):
                    rows.append(
                        {
                            "sample_id": f"{genotype}_{context}_{treatment}_{replicate}",
                            "genotype": genotype,
                            "DP_status": context,
                            "treatment": treatment,
                        }
                    )
    samples = pd.DataFrame(rows)
    columns = samples["sample_id"].tolist()
    psi = pd.DataFrame(index=["event_1", "event_2", "event_3"], columns=columns, dtype=float)
    scores = pd.DataFrame(index=columns, columns=["M01", "M02"], dtype=float)
    for row in samples.itertuples(index=False):
        context_shift = 0.03 if row.DP_status == "DP" else 0.0
        genotype_shift = 0.12 if row.genotype == "MUT" else 0.0
        drug_shift = 0.22 if row.treatment == "CLK3" else 0.0
        attenuation = -0.14 if row.genotype == "MUT" and row.treatment == "CLK3" else 0.0
        jitter = (int(row.sample_id.rsplit("_", 1)[1]) - 2) * 0.005
        psi.loc["event_1", row.sample_id] = 0.25 + context_shift + genotype_shift + drug_shift + attenuation + jitter
        psi.loc["event_2", row.sample_id] = 0.70 - genotype_shift - drug_shift - attenuation + jitter
        psi.loc["event_3", row.sample_id] = 0.45 + jitter
        scores.loc[row.sample_id, "M01"] = genotype_shift + drug_shift + attenuation + jitter
        scores.loc[row.sample_id, "M02"] = jitter
    event_metadata = pd.DataFrame(
        {
            "event_key": psi.index,
            "gene_id": ["G1", "G2", "G3"],
            "gene_symbol": ["Gene1", "Gene2", "Gene3"],
            "event_type": ["SE", "RI", "A3SS"],
        }
    )
    settings = FactorContrastConfig(
        enabled=True,
        primary_factor="genotype",
        target_level="MUT",
        reference_level="WT",
        context_factors=["DP_status"],
        perturbation_factor="treatment",
        baseline_level="Baseline",
        perturbation_levels=["CLK3"],
        network_source="residual",
        fdr_cutoff=0.05,
        min_abs_delta_psi=0.10,
        min_samples_per_group=3,
    )
    result = analyze_factor_contrasts(psi, scores, samples, event_metadata, settings)
    direct = result.event_contrasts.loc[
        result.event_contrasts["contrast_type"].eq("primary_baseline")
        & result.event_contrasts["event_key"].eq("event_1")
    ]
    assert len(direct) == 2
    assert np.allclose(direct["effect"], 0.12)
    interaction = result.event_contrasts.loc[
        result.event_contrasts["contrast_type"].eq("perturbation_occlusion_interaction")
        & result.event_contrasts["event_key"].eq("event_1")
    ]
    assert len(interaction) == 2
    assert np.allclose(interaction["effect"], -0.14)
    assert set(result.phenocopy_summary["context"]) == {"DP_status=noDP", "DP_status=DP"}
    plots = render_factor_contrast_plots(result, tmp_path / "plots")
    assert {path.name for path in plots} == {
        "factor_phenocopy_dashboard.png",
        "factor_phenocopy_signature_scatter.png",
        "module_perturbation_occlusion_heatmap.png",
    }


def test_bootstrap_coassignment_builds_confident_consensus_modules() -> None:
    rng = np.random.default_rng(101)
    samples = pd.DataFrame(
        {
            "sample_id": [f"S{i:02d}" for i in range(18)],
            "condition": np.repeat(list("ABCDEF"), 3),
        }
    )
    pattern_a = np.repeat(np.linspace(0.15, 0.85, 6), 3)
    pattern_b = np.repeat(np.linspace(0.85, 0.15, 6), 3)
    values = np.vstack(
        [pattern_a + rng.normal(0, 0.01, 18) for _ in range(12)]
        + [pattern_b + rng.normal(0, 0.01, 18) for _ in range(12)]
    )
    matrix = pd.DataFrame(
        np.clip(values, 0.01, 0.99),
        index=[f"event_{index:02d}" for index in range(24)],
        columns=samples["sample_id"],
    )
    observed = pd.DataFrame(
        {
            "event_key": matrix.index,
            "module": ["M01"] * 12 + ["M02"] * 12,
        }
    )
    result = build_bootstrap_coassignment_consensus(
        matrix,
        observed,
        samples,
        NetworkConfig(
            min_module_size=5,
            knn_neighbors=8,
            consensus_knn_neighbors=[6, 8],
            consensus_bootstrap_iterations=0,
            community_consensus_runs=1,
            bootstrap_iterations=0,
            edge_threshold=0.01,
        ),
        ConsensusConfig(
            bootstrap_iterations=4,
            min_coassignment=0.50,
            min_module_size=5,
            dense_max_events=100,
            output_edge_threshold=0.25,
        ),
    )
    assert result.assignments["event_key"].is_unique
    assert result.coassignment_matrix.shape == (24, 24)
    assert np.allclose(np.diag(result.coassignment_matrix), 1.0)
    assert result.assignments["consensus_module"].ne("grey").any()
    assert not result.module_summary.empty
    assert len(result.bootstrap_runs) == 4
    assert "consensus_membership_confidence" in result.assignments
    assert "membership_role" in result.assignments
    assert set(result.assignments["membership_role"]).issubset(
        {"core", "peripheral", "uncertain", "grey"}
    )
    assert {
        "candidate_graph_module",
        "analysis_module",
        "analysis_assignment_scope",
        "partition_source",
        "confidence_action",
        "membership_evidence",
        "full_data_kME",
        "bootstrap_mean_abs_kME",
        "bootstrap_q10_abs_kME",
        "bootstrap_kME_sign_consistency",
        "split_half_membership_role_a",
        "split_half_membership_role_b",
        "split_half_core_confirmed",
    }.issubset(result.assignments)
    assert set(result.assignments["analysis_assignment_scope"]).issubset(
        {
            "wgcna_backbone_confident",
            "wgcna_backbone_low_confidence",
            "unassigned",
        }
    )
    assert set(result.assignments["partition_source"]) == {
        "full_data_wgcna_backbone"
    }
    assert set(result.assignments["confidence_action"]) == {"annotation_only"}
    assert adjusted_rand_score(
        observed["module"], result.assignments["consensus_module"]
    ) == 1.0
    assert result.parameter_search["selected"].sum() == 1
    assert {
        "dynamic_deep_split",
        "selection_score",
        "robust_bootstrap_stability",
    }.issubset(result.parameter_search)
    assert set(result.parameter_search["method"]) == {
        "tom_backbone_confidence_only"
    }
    assert set(result.parameter_search["partition_policy"]) == {
        "fixed_wgcna_backbone"
    }
    assert set(result.parameter_search["constraints_role"]) == {
        "diagnostic_only"
    }
    assert set(result.parameter_search["selection_method"]) == {
        "fixed_wgcna_backbone_no_selection"
    }
    assert "backbone_cut_height" in result.parameter_search
    assert not bool(result.diagnostics.loc[0, "partition_changed_from_backbone"])
    assert result.selection_audit["audit_type"].eq(
        "tom_backbone_selection"
    ).any()


def test_tom_backbone_strengthens_shared_neighbors() -> None:
    correlation = np.asarray(
        [
            [1.0, 0.9, 0.9, 0.0],
            [0.9, 1.0, 0.9, 0.0],
            [0.9, 0.9, 1.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ]
    )
    settings = NetworkConfig(soft_power=2, min_module_size=2)
    adjacency = signed_soft_adjacency(correlation, settings)
    tom = topological_overlap_matrix(adjacency)

    assert NetworkConfig().engine == "tom_hierarchical"
    assert NetworkConfig().module_detection == "wgcna_blockwise"
    assert ConsensusConfig().clustering_method == "backbone_confidence_only"
    assert np.allclose(tom, tom.T)
    assert np.allclose(np.diag(tom), 1.0)
    assert tom[0, 1] > tom[0, 3]


def test_confidence_only_keeps_low_confidence_backbone_events_assigned() -> None:
    coassignment = np.asarray(
        [
            [1.0, 0.60, 0.40, 0.40],
            [0.60, 1.0, 0.40, 0.40],
            [0.40, 0.40, 1.0, 0.60],
            [0.40, 0.40, 0.60, 1.0],
        ]
    )
    labels, roles, confidence, second, margin, candidates = (
        _backbone_confidence_only_membership(
            coassignment,
            np.asarray(["M01", "M01", "M02", "M02"], dtype=object),
            2,
            ConsensusConfig(
                core_membership_min=0.90,
                peripheral_membership_min=0.80,
                peripheral_membership_margin=0.30,
            ),
        )
    )

    assert set(labels) == {"C01", "C02"}
    assert np.all(roles == "uncertain")
    assert np.allclose(confidence, 0.60)
    assert np.allclose(second, 0.40)
    assert np.allclose(margin, 0.20)
    assert np.array_equal(labels, candidates)


def test_bootstrap_kme_can_rescue_stable_diffuse_backbone_membership() -> None:
    labels = np.asarray(["C01", "C01", "C02", "C02", "grey"], dtype=object)
    coassignment_roles = np.asarray(
        ["uncertain", "core", "uncertain", "uncertain", "grey"],
        dtype=object,
    )
    margins = np.asarray([0.12, 0.20, -0.05, 0.10, 0.0])
    kme = pd.DataFrame(
        {
            "bootstrap_mean_abs_kME": [0.82, 0.84, 0.90, 0.58, np.nan],
            "bootstrap_q10_abs_kME": [0.66, 0.68, 0.72, 0.38, np.nan],
            "bootstrap_kME_sign_consistency": [0.90, 0.95, 0.95, 0.80, np.nan],
        }
    )

    roles, evidence = _integrate_kme_membership_evidence(
        labels,
        coassignment_roles,
        margins,
        kme,
        ConsensusConfig(),
    )

    assert roles.tolist() == ["core", "core", "uncertain", "peripheral", "grey"]
    assert evidence.tolist() == [
        "bootstrap_kME",
        "coassignment_and_bootstrap_kME",
        "insufficient_confidence",
        "bootstrap_kME",
        "unassigned",
    ]


def test_hybrid_programs_retain_tight_and_diffuse_assignments() -> None:
    rng = np.random.default_rng(303)
    samples = pd.DataFrame(
        {
            "sample_id": [f"S{index:02d}" for index in range(24)],
            "condition": np.repeat(list("ABCDEFGH"), 3),
        }
    )
    diffuse_a = np.repeat(np.linspace(0.15, 0.85, 8), 3)
    diffuse_b = np.repeat(np.linspace(0.85, 0.15, 8), 3)
    matrix = pd.DataFrame(
        np.clip(
            np.vstack(
                [diffuse_a + rng.normal(0, 0.025, 24) for _ in range(12)]
                + [diffuse_b + rng.normal(0, 0.025, 24) for _ in range(12)]
            ),
            0.01,
            0.99,
        ),
        index=[f"event_{index:02d}" for index in range(24)],
        columns=samples["sample_id"],
    )
    tight = pd.DataFrame(
        {
            "event_key": matrix.index,
            "module": ["M01"] * 6 + ["grey"] * 6 + ["M02"] * 6 + ["grey"] * 6,
            "module_membership": [0.92] * 6 + [np.nan] * 6 + [0.90] * 6 + [np.nan] * 6,
            "mean_incident_bootstrap_support": [0.88] * 6
            + [np.nan] * 6
            + [0.86] * 6
            + [np.nan] * 6,
        }
    )
    qc = pd.DataFrame(
        {"replicate_repeatability": np.full(len(matrix), 0.82)},
        index=matrix.index,
    )
    result = build_hybrid_programs(
        matrix,
        qc,
        samples,
        tight,
        NetworkConfig(
            method="pearson",
            min_module_size=4,
            soft_power=4,
        ),
        BenchmarkConfig(max_events=30, bootstrap_iterations=1),
        HybridProgramConfig(
            diffuse_engine="python_soft_tom",
            bootstrap_iterations=2,
            min_module_size=4,
            python_tom_cut_height=0.75,
            python_max_events=100,
        ),
    )

    assert result.event_assignments["event_key"].is_unique
    assert {"tight_program", "diffuse_program", "program_scope"}.issubset(
        result.event_assignments
    )
    assert result.event_assignments["tight_program"].ne("grey").sum() == 12
    assert result.event_assignments["diffuse_program"].ne("grey").any()
    assert set(result.program_summary["program_type"]) == {"tight", "diffuse"}
    assert result.program_summary["overall_confidence"].between(0, 1).all()
    assert len(result.bootstrap_runs) == (
        2
        * result.event_assignments.loc[
            result.event_assignments["diffuse_program"].ne("grey"),
            "diffuse_program",
        ].nunique()
    )


def test_external_preservation_uses_independent_target_and_exact_events(tmp_path: Path) -> None:
    rng = np.random.default_rng(202)
    event_keys = [f"event_{index:02d}" for index in range(24)]
    reference_samples = [f"R{index:02d}" for index in range(18)]
    target_samples = [f"T{index:02d}" for index in range(15)]
    ref_a = np.repeat(np.linspace(0.15, 0.85, 6), 3)
    ref_b = np.repeat(np.linspace(0.85, 0.15, 6), 3)
    target_a = np.repeat(np.linspace(0.18, 0.82, 5), 3)
    target_b = np.repeat(np.linspace(0.82, 0.18, 5), 3)
    reference = pd.DataFrame(
        np.vstack(
            [ref_a + rng.normal(0, 0.01, 18) for _ in range(12)]
            + [ref_b + rng.normal(0, 0.01, 18) for _ in range(12)]
        ),
        index=event_keys,
        columns=reference_samples,
    )
    target = pd.DataFrame(
        np.vstack(
            [target_a + rng.normal(0, 0.015, 15) for _ in range(12)]
            + [target_b + rng.normal(0, 0.015, 15) for _ in range(12)]
        ),
        index=event_keys,
        columns=target_samples,
    )
    target_path = tmp_path / "target_psi.tsv"
    metadata_path = tmp_path / "target_samples.tsv"
    target.to_csv(target_path, sep="\t")
    pd.DataFrame(
        {"sample_id": target_samples, "condition": np.repeat(list("ABCDE"), 3)}
    ).to_csv(metadata_path, sep="\t", index=False)
    assignments = pd.DataFrame(
        {
            "event_key": event_keys,
            "consensus_module": ["C01"] * 12 + ["C02"] * 12,
        }
    )
    result = evaluate_external_module_preservation(
        reference,
        assignments,
        [ExternalCohortConfig("independent", target_path, metadata_path)],
        PreservationConfig(
            permutations=10,
            min_shared_events=5,
            dense_max_events=100,
            preserved_zsummary=2.0,
            partial_zsummary=0.0,
            random_seed=9,
        ),
    )
    assert set(result.module_preservation["consensus_module"]) == {"C01", "C02"}
    assert result.event_mapping["event_key"].nunique() == 24
    assert result.target_audit.iloc[0]["status"] == "finished"
    assert {
        "Zdensity",
        "Zconnectivity",
        "Zsummary_wgcna_like",
        "membership_profile_correlation",
        "program_score_preservation",
        "preservation_class",
    }.issubset(result.module_preservation.columns)


def test_module_annotation_adds_support_without_changing_membership(tmp_path: Path) -> None:
    support_root = tmp_path / "support" / "A-vs-B"
    support_root.mkdir(parents=True)
    pd.DataFrame(
        {
            "comparison_id": ["A-vs-B", "A-vs-B"],
            "gene_id": ["ENSMUSG1.1", "ENSMUSG2.2"],
            "gene_symbol": ["Srsf1", "Gene2"],
            "DEG_significant": [True, False],
            "rMATS_significant": [True, True],
            "DEXSeq_significant": [False, True],
            "DTU_significant": [True, False],
            "DE_FDR": [0.01, 0.8],
            "best_rMATS_FDR": [0.001, 0.02],
            "best_DEXSeq_qvalue": [0.8, 0.01],
            "best_DTU_qvalue": [0.02, 0.9],
        }
    ).to_csv(support_root / "gene_level_integrated_candidates.tsv", sep="\t", index=False)
    assignments = pd.DataFrame(
        {
            "event_key": ["E1", "E2", "E3"],
            "module": ["M01", "M01", "M02"],
            "network_connectivity": [0.8, 0.7, 0.5],
        }
    )
    events = pd.DataFrame(
        {
            "event_key": ["E1", "E2", "E3"],
            "gene_id": ["ENSMUSG1", "ENSMUSG2", "ENSMUSG3"],
            "gene_symbol": ["Srsf1", "Gene2", "Gene3"],
            "event_type": ["SE", "RI", "A3SS"],
        }
    )

    result = annotate_modules(assignments, events, support_root=tmp_path / "support")

    assert result.events[["event_key", "module"]].to_records(index=False).tolist() == assignments[["event_key", "module"]].to_records(index=False).tolist()
    assert result.genes["gene_key"].nunique() == 3
    srsf1 = result.genes.loc[result.genes["gene_symbol"].eq("Srsf1")].iloc[0]
    assert bool(srsf1["DEG_significant"])
    assert bool(srsf1["DTU_significant"])
    assert bool(srsf1["RBP_gene"])
    assert set(result.module_summary["module"]) == {"M01", "M02"}
    assert not result.rbp_enrichment.empty


def test_module_interpretation_uses_configured_rbp_sources_without_causal_claims(tmp_path: Path) -> None:
    samples = [f"S{index}" for index in range(6)]
    assignments = pd.DataFrame(
        {
            "event_key": ["E1", "E2", "E3", "E4"],
            "consensus_module": ["C01", "C01", "C02", "C02"],
            "consensus_membership_confidence": [0.95, 0.90, 0.92, 0.88],
        }
    )
    psi = pd.DataFrame(
        [
            [0.1, 0.2, 0.3, 0.7, 0.8, 0.9],
            [0.15, 0.25, 0.35, 0.65, 0.75, 0.85],
            [0.9, 0.8, 0.7, 0.3, 0.2, 0.1],
            [0.85, 0.75, 0.65, 0.35, 0.25, 0.15],
        ],
        index=["E1", "E2", "E3", "E4"],
        columns=samples,
    )
    scores = pd.DataFrame({"C01": np.linspace(-1, 1, 6), "C02": np.linspace(1, -1, 6)}, index=samples)
    events = pd.DataFrame(
        {
            "event_key": ["E1", "E2", "E3", "E4"],
            "gene_id": ["G1", "G2", "G3", "G4"],
            "gene_symbol": ["Srsf1", "Gene2", "Gene3", "Gene4"],
            "event_type": ["SE", "RI", "A3SS", "A5SS"],
        }
    )
    motif_path = tmp_path / "motif_hits.tsv"
    pd.DataFrame({"event_key": ["E1", "E2"], "RBP": ["Srsf1", "Srsf1"]}).to_csv(motif_path, sep="\t", index=False)
    expression_path = tmp_path / "rbp_expression.tsv"
    pd.DataFrame([np.linspace(1, 10, 6)], index=["Srsf1"], columns=samples).to_csv(expression_path, sep="\t")
    gmt_path = tmp_path / "genes.gmt"
    gmt_path.write_text("GO_TEST\tdescription\tSrsf1\tGene2\n", encoding="utf-8")

    result = interpret_modules(
        assignments,
        psi,
        scores,
        events,
        ModuleInterpretationConfig(
            rbp_expression_matrix=expression_path,
            rbp_motif_hits=motif_path,
            gene_sets_gmt=gmt_path,
            min_module_genes=2,
        ),
    )

    assert result.hub_events["event_key"].nunique() == 4
    assert not result.rbp_motif_enrichment.empty
    assert not result.rbp_expression_correlation.empty
    assert not result.gene_set_enrichment.empty
    assert result.regulator_ranking["interpretation"].str.contains("not causal").all()


def test_beta_binomial_statistics_and_simple_delta_head_to_head() -> None:
    samples = pd.DataFrame(
        {"sample_id": ["A1", "A2", "A3", "B1", "B2", "B3"], "condition": ["A"] * 3 + ["B"] * 3}
    )
    psi = pd.DataFrame(
        [[0.10, 0.12, 0.11, 0.80, 0.82, 0.81], [0.40, 0.42, 0.41, 0.44, 0.43, 0.45]],
        index=["strong", "weak"], columns=samples["sample_id"],
    )
    totals = pd.DataFrame(100.0, index=psi.index, columns=psi.columns)
    ijc = (psi * totals).round()
    sjc = totals - ijc
    assignments = pd.DataFrame({"event_key": psi.index, "consensus_module": ["C01", "C01"]})
    scores = pd.DataFrame({"C01": [-1.1, -0.9, -1.0, 1.0, 1.1, 0.9]}, index=samples["sample_id"])

    result = fit_replicate_aware_models(
        psi, ijc, sjc, samples, assignments, scores,
        ModuleStatisticsConfig(min_total_reads=10, optimizer_max_iterations=100),
    )

    assert set(result.event_models["event_key"]) == {"strong", "weak"}
    strong = result.event_models.set_index("event_key").loc["strong"]
    assert strong["simple_max_delta_PSI"] > 0.6
    assert strong["condition_effect_FDR"] < 0.05
    assert not result.module_models.empty
    assert not result.simple_delta_benchmark.empty


def test_beta_binomial_rejects_condition_batch_confounding() -> None:
    samples = pd.DataFrame(
        {
            "sample_id": ["A1", "A2", "A3", "B1", "B2", "B3"],
            "condition": ["A"] * 3 + ["B"] * 3,
            "batch": ["batch_A"] * 3 + ["batch_B"] * 3,
        }
    )
    psi = pd.DataFrame([[0.1, 0.2, 0.15, 0.7, 0.8, 0.75]], index=["E1"], columns=samples["sample_id"])
    ijc = (psi * 100).round()
    sjc = 100 - ijc
    assignments = pd.DataFrame({"event_key": ["E1"], "module": ["M01"]})
    scores = pd.DataFrame({"M01": [-1.0, -0.8, -0.9, 0.8, 1.0, 0.9]}, index=samples["sample_id"])

    with pytest.raises(ValueError, match="confounded"):
        fit_replicate_aware_models(
            psi,
            ijc,
            sjc,
            samples,
            assignments,
            scores,
            ModuleStatisticsConfig(min_total_reads=10),
            batch_column="batch",
        )


def test_cross_cohort_meta_uses_exact_events_and_independent_assignments(tmp_path: Path) -> None:
    samples = ["A1", "A2", "A3", "B1", "B2", "B3"]
    events = [f"E{index}" for index in range(8)]
    reference = pd.DataFrame(
        np.vstack([
            np.linspace(0.1, 0.8, 6) + index * 0.002 if index < 4 else np.linspace(0.8, 0.1, 6) + index * 0.002
            for index in range(8)
        ]), index=events, columns=samples,
    )
    target = reference + 0.01
    metadata = pd.DataFrame({"sample_id": samples, "condition": ["A"] * 3 + ["B"] * 3})
    assignments = pd.DataFrame({"event_key": events, "consensus_module": ["C01"] * 4 + ["C02"] * 4})
    target_assignments = pd.DataFrame({"event_key": events, "consensus_module": ["T01"] * 4 + ["T02"] * 4})
    psi_path = tmp_path / "target.tsv"
    metadata_path = tmp_path / "metadata.tsv"
    assignments_path = tmp_path / "assignments.tsv"
    target.to_csv(psi_path, sep="\t")
    metadata.to_csv(metadata_path, sep="\t", index=False)
    target_assignments.to_csv(assignments_path, sep="\t", index=False)

    result = meta_analyze_modules(
        reference,
        assignments,
        metadata,
        [ExternalCohortConfig("external", psi_path, metadata_path, module_assignments=assignments_path)],
        MetaAnalysisConfig(min_shared_events=3, min_cohorts=2),
        reference_name="discovery",
    )

    assert result.event_alignment["exact_match"].sum() == 8
    assert set(result.module_state["module_state"]) == {"consensus"}
    assert not result.random_effects_meta.empty
    assert result.random_effects_meta["method"].str.contains("DerSimonian-Laird").all()
    assert result.module_matching["significant_exact_event_overlap"].all()


def test_cross_cohort_meta_reports_not_configured_instead_of_inventing_validation() -> None:
    samples = ["A1", "A2", "B1", "B2"]
    events = [f"E{index}" for index in range(4)]
    psi = pd.DataFrame(np.arange(16, dtype=float).reshape(4, 4) / 20.0, index=events, columns=samples)
    metadata = pd.DataFrame({"sample_id": samples, "condition": ["A", "A", "B", "B"]})
    assignments = pd.DataFrame({"event_key": events, "consensus_module": ["C01"] * 4})

    result = meta_analyze_modules(
        psi,
        assignments,
        metadata,
        [],
        MetaAnalysisConfig(min_shared_events=3, min_cohorts=2),
    )

    assert result.cohort_audit.iloc[0]["status"] == "not_configured"
    assert result.module_state.empty


def test_method_benchmark_uses_one_shared_event_universe() -> None:
    rng = np.random.default_rng(19)
    sample_ids = [f"S{index:02d}" for index in range(18)]
    conditions = np.repeat(["A", "B", "C", "D", "E", "F"], 3)
    latent_a = np.repeat(np.linspace(0.15, 0.85, 6), 3)
    latent_b = np.repeat(np.linspace(0.85, 0.15, 6), 3)
    values = np.vstack(
        [latent_a + rng.normal(0, 0.02, 18) for _ in range(30)]
        + [latent_b + rng.normal(0, 0.02, 18) for _ in range(30)]
    )
    matrix = pd.DataFrame(np.clip(values, 0.01, 0.99), index=[f"E{index:03d}" for index in range(60)], columns=sample_ids)
    qc = pd.DataFrame(
        {
            "event_key": matrix.index,
            "stability_score": np.linspace(1.0, 0.1, len(matrix)),
        }
    )
    samples = pd.DataFrame({"sample_id": sample_ids, "condition": conditions})
    settings = BenchmarkConfig(
        methods=["katana_sparse", "pearson_hierarchical", "graphical_lasso"],
        max_events=50,
        bootstrap_iterations=2,
        graphical_lasso_alpha=0.15,
        neutral_knn_neighbors=8,
    )
    result = run_method_benchmark(
        matrix,
        qc,
        samples,
        NetworkConfig(min_module_size=5, knn_neighbors=8, edge_threshold=0.05, bootstrap_iterations=0),
        settings,
    )

    assert len(result.universe) == 50
    assert result.assignments["event_key"].is_unique
    assert set(result.assignments.columns) >= {
        "module_katana_sparse",
        "module_pearson_hierarchical",
        "module_graphical_lasso",
    }
    assert set(result.method_summary["status"]).issubset({"finished", "degenerate"})
    assert result.method_summary.set_index("method").loc["katana_sparse", "status"] == "finished"
    assert set(result.metric_ranking.columns) >= {
        "method",
        "metric",
        "metric_rank",
        "normalized_descriptive_score",
    }
    assert result.method_interpretation["method"].is_unique
    assert "claim_scope" in result.method_interpretation
    assert len(result.bootstrap_stability) == 2 * 3
    assert not result.module_matching.empty


def test_method_benchmark_accepts_read_only_dataframe_storage() -> None:
    rng = np.random.default_rng(190)
    sample_ids = [f"S{index:02d}" for index in range(12)]
    values = rng.uniform(0.05, 0.95, size=(24, 12))
    values.setflags(write=False)
    matrix = pd.DataFrame(
        values,
        index=[f"E{index:03d}" for index in range(24)],
        columns=sample_ids,
        copy=False,
    )
    qc = pd.DataFrame(
        {
            "event_key": matrix.index,
            "stability_score": np.linspace(1.0, 0.1, len(matrix)),
        }
    )
    samples = pd.DataFrame(
        {
            "sample_id": sample_ids,
            "condition": np.repeat(["A", "B", "C", "D"], 3),
        }
    )

    result = run_method_benchmark(
        matrix,
        qc,
        samples,
        NetworkConfig(
            min_module_size=4,
            knn_neighbors=5,
            edge_threshold=0.05,
            bootstrap_iterations=0,
        ),
        BenchmarkConfig(
            methods=["katana_sparse"],
            max_events=24,
            bootstrap_iterations=0,
            neutral_knn_neighbors=5,
        ),
    )

    assert result.method_summary.iloc[0]["status"] in {"finished", "degenerate"}


def test_method_benchmark_includes_cached_final_consensus_without_rebuilding() -> None:
    rng = np.random.default_rng(191)
    sample_ids = [f"S{index:02d}" for index in range(18)]
    conditions = np.repeat(["A", "B", "C", "D", "E", "F"], 3)
    latent_a = np.repeat(np.linspace(0.15, 0.85, 6), 3)
    latent_b = np.repeat(np.linspace(0.85, 0.15, 6), 3)
    matrix = pd.DataFrame(
        np.clip(
            np.vstack(
                [latent_a + rng.normal(0, 0.02, 18) for _ in range(20)]
                + [latent_b + rng.normal(0, 0.02, 18) for _ in range(20)]
            ),
            0.01,
            0.99,
        ),
        index=[f"E{index:03d}" for index in range(40)],
        columns=sample_ids,
    )
    qc = pd.DataFrame(
        {
            "event_key": matrix.index,
            "stability_score": np.linspace(1.0, 0.1, len(matrix)),
        }
    )
    samples = pd.DataFrame({"sample_id": sample_ids, "condition": conditions})
    consensus = pd.DataFrame(
        {
            "event_key": matrix.index,
            "consensus_module": ["C01"] * 20 + ["C02"] * 20,
        }
    )
    result = run_method_benchmark(
        matrix,
        qc,
        samples,
        NetworkConfig(
            min_module_size=5,
            knn_neighbors=8,
            consensus_knn_neighbors=[6, 8],
            consensus_bootstrap_iterations=0,
            community_consensus_runs=1,
            edge_threshold=0.05,
            bootstrap_iterations=0,
        ),
        BenchmarkConfig(
            methods=["katana_sparse"],
            max_events=40,
            bootstrap_iterations=2,
            neutral_knn_neighbors=8,
        ),
        consensus_assignments=consensus,
    )

    assert "module_katana_consensus" in result.assignments
    consensus_summary = result.method_summary.set_index("method").loc["katana_consensus"]
    assert consensus_summary["status"] == "finished"
    assert np.isfinite(consensus_summary["mean_bootstrap_ARI"])
    assert set(result.bootstrap_stability["method"]) == {
        "katana_sparse",
        "katana_consensus",
    }
    audit = result.method_audit.set_index("method").loc["katana_consensus"]
    assert "loaded from cache" in audit["message"]


def test_method_benchmark_rejects_consensus_from_different_event_universe() -> None:
    sample_ids = [f"S{index}" for index in range(12)]
    matrix = pd.DataFrame(
        np.random.default_rng(17).normal(size=(20, 12)),
        index=[f"E{index:03d}" for index in range(20)],
        columns=sample_ids,
    )
    qc = pd.DataFrame(
        {"event_key": matrix.index, "stability_score": np.arange(20, 0, -1)}
    )
    samples = pd.DataFrame(
        {"sample_id": sample_ids, "condition": ["A"] * 6 + ["B"] * 6}
    )
    consensus = pd.DataFrame(
        {
            "event_key": [*matrix.index, "EXTRA_EVENT"],
            "consensus_module": ["C01"] * 10 + ["C02"] * 11,
        }
    )
    result = run_method_benchmark(
        matrix,
        qc,
        samples,
        NetworkConfig(
            min_module_size=3,
            knn_neighbors=5,
            consensus_knn_neighbors=[5],
            consensus_bootstrap_iterations=0,
            community_consensus_runs=1,
            bootstrap_iterations=0,
        ),
        BenchmarkConfig(
            methods=["katana_sparse"],
            max_events=20,
            bootstrap_iterations=1,
            neutral_knn_neighbors=5,
        ),
        consensus_assignments=consensus,
    )

    consensus_summary = result.method_summary.set_index("method").loc[
        "katana_consensus"
    ]
    assert consensus_summary["status"] == "unavailable"
    audit = result.method_audit.set_index("method").loc["katana_consensus"]
    assert "different event universe" in audit["message"]


def test_shared_neighbor_smoothing_is_sparse_and_auditable() -> None:
    rng = np.random.default_rng(31)
    latent = np.linspace(-1.0, 1.0, 18)
    values = np.vstack(
        [latent + rng.normal(0, 0.08, 18) for _ in range(12)]
        + [-latent + rng.normal(0, 0.08, 18) for _ in range(12)]
    )
    matrix = pd.DataFrame(values, index=[f"E{index:02d}" for index in range(24)])
    base = NetworkConfig(
        min_module_size=4,
        knn_neighbors=5,
        edge_threshold=0.01,
        edge_smoothing="none",
        bootstrap_iterations=0,
    )
    smooth = NetworkConfig(
        min_module_size=4,
        knn_neighbors=5,
        edge_threshold=0.01,
        edge_smoothing="shared_neighbor",
        shared_neighbor_weight=0.50,
        shared_neighbor_top_k=5,
        bootstrap_iterations=0,
    )
    base_edges = _sparse_knn_edges(matrix, base)
    smooth_edges = _sparse_knn_edges(matrix, smooth)

    assert len(smooth_edges[0]) >= len(base_edges[0])
    assert np.isfinite(smooth_edges[3]).all()
    assert (smooth_edges[4] > 0).any()
    assert len(smooth_edges[0]) <= len(base_edges[0]) + len(matrix) * smooth.shared_neighbor_top_k


def test_adaptive_membership_cutoff_uses_module_distribution_and_bounds() -> None:
    statistics = pd.DataFrame(
        {
            "assigned_module": ["M01"] * 5 + ["grey"],
            "module_membership": [0.30, 0.32, 0.34, 0.36, 0.38, np.nan],
            "best_module": ["M01"] * 6,
            "best_module_membership": [0.30, 0.32, 0.34, 0.36, 0.38, 0.01],
        }
    )
    statistics.attrs["n_samples"] = 100
    adaptive = NetworkConfig(
        min_module_size=5,
        min_module_membership=0.45,
        adaptive_membership=True,
        adaptive_membership_floor=0.08,
    )
    fixed = NetworkConfig(
        min_module_size=5,
        min_module_membership=0.45,
        adaptive_membership=False,
    )

    adaptive_cutoffs = _adaptive_membership_cutoffs(statistics, adaptive)
    fixed_cutoffs = _adaptive_membership_cutoffs(statistics, fixed)
    assert np.all((adaptive_cutoffs >= 0.08) & (adaptive_cutoffs <= 0.45))
    assert adaptive_cutoffs[0] > 0.08
    assert adaptive_cutoffs[-1] == 0.45
    assert np.allclose(fixed_cutoffs, 0.45)


def test_pipeline_single_stage_writes_artifact_manifest_without_chaining(tmp_path: Path) -> None:
    pipeline = CoSplicePipeline(_synthetic_project(tmp_path))
    pipeline.run_stage("audit")

    manifest_path = tmp_path / "output" / "artifact_manifest.json"
    assert manifest_path.exists()
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    assert manifest["stages"]["audit"]["status"] == "finished"
    assert manifest["stages"]["audit"]["explicit_only"] is True
    assert "matrix" not in manifest["stages"]
    assert not (tmp_path / "output" / "02_matrix" / "sample_level_psi_matrix.tsv").exists()


def test_registered_cli_command_dispatches_its_stage(monkeypatch: pytest.MonkeyPatch) -> None:
    calls: list[tuple[str, bool]] = []

    class FakePipeline:
        def __init__(self, _config: object) -> None:
            pass

        def run_stage(self, stage: str, force: bool = False) -> None:
            calls.append((stage, force))

    monkeypatch.setattr(cosplice_cli, "load_config", lambda _path: object())
    monkeypatch.setattr(cosplice_cli, "CoSplicePipeline", FakePipeline)

    assert cosplice_cli.main(["meta-analyze-modules", "--config", "unused.yaml", "--force"]) == 0
    assert calls == [("meta_analyze_modules", True)]


def test_project_facade_exposes_contracts_without_running_analysis(tmp_path: Path) -> None:
    project = CoSpliceProject(_synthetic_project(tmp_path))
    assert project.state()["stages"] == {}
    assert project.artifact_manifest()["stages"] == {}
    assert project.stages()[0]["key"] == "audit"
    assert not project.config.output_dir.exists()


def test_factorial_association_uses_configured_factors_and_batch() -> None:
    samples = pd.DataFrame(
        {
            "sample_id": [f"S{i:02d}" for i in range(24)],
            "genotype": np.tile(np.repeat(["WT", "MUT"], 6), 2),
            "DP_status": np.repeat(["DP", "noDP"], 12),
            "treatment": np.tile(np.repeat(["Baseline", "Drug"], 3), 4),
            "batch": np.tile([1, 2, 3], 8),
        }
    )
    genotype_effect = samples["genotype"].map({"WT": -1.0, "MUT": 1.0}).to_numpy()
    scores = pd.DataFrame(
        {
            "M01": genotype_effect + np.tile([-0.05, 0.0, 0.05], 8),
            "M02": np.linspace(-0.2, 0.2, 24),
        },
        index=samples["sample_id"],
    )

    anova, coefficients = factorial_program_association(
        scores,
        samples,
        ["genotype", "DP_status", "treatment"],
        "batch",
    )

    assert not anova.empty
    assert not coefficients.empty
    assert anova["model_term"].str.contains('Q("genotype")', regex=False).any()
    assert anova["model_term"].str.contains('Q("batch")', regex=False).any()
    assert {"partial_eta_squared", "FDR", "model_r_squared"}.issubset(anova.columns)


def test_motif_scanner_reads_indexed_fasta_and_reports_strand_oriented_hits(tmp_path: Path) -> None:
    fasta = tmp_path / "genome.fa"
    sequence = "C" * 20 + "AAAUUUAAA".replace("U", "T") + "C" * 71
    fasta.write_text(f">chr1\n{sequence}\n", encoding="ascii")
    Path(str(fasta) + ".fai").write_text(f"chr1\t{len(sequence)}\t6\t{len(sequence)}\t{len(sequence) + 1}\n", encoding="ascii")
    motifs = tmp_path / "motifs.tsv"
    pd.DataFrame({"RBP": ["RbpA"], "motif_id": ["M1"], "motif": ["AAA"], "source": ["test"]}).to_csv(
        motifs, sep="\t", index=False
    )
    events = pd.DataFrame({
        "event_key": ["Gene1|SE|chr1|+|20|29"],
        "chr": ["chr1"], "strand": ["+"], "event_coordinates": ["20|29"],
    })

    result = scan_event_motifs(events, fasta, motifs, window_bp=12, score_threshold=0.8)

    assert len(result) == 1
    assert result.iloc[0]["RBP"] == "RbpA"
    assert result.iloc[0]["site_count"] >= 2
    assert bool(result.iloc[0]["strand_oriented"])


def test_motif_scanner_reads_official_attract_zip_shape_and_filters_organism(tmp_path: Path) -> None:
    archive_path = tmp_path / "ATtRACT.zip"
    metadata = (
        "Gene_name\tMatrix_id\tOrganism\tDatabase\n"
        "RbpMouse\tM1\tMus_musculus\tATtRACT\n"
        "RbpHuman\tM2\tHomo_sapiens\tATtRACT\n"
    )
    pwm = ">M1 3\n1 0 0 0\n1 0 0 0\n1 0 0 0\n>M2 3\n0 1 0 0\n0 1 0 0\n0 1 0 0\n"
    with ZipFile(archive_path, "w") as archive:
        archive.writestr("ATtRACT_db.txt", metadata)
        archive.writestr("pwm.txt", pwm)

    motifs = load_motif_database(archive_path, organism="Mus_musculus")

    assert len(motifs) == 1
    assert motifs[0].rbp == "RbpMouse"
    assert motifs[0].pwm is not None and motifs[0].pwm.shape == (3, 4)


def test_motif_scanner_reads_cisbp_style_long_pwm(tmp_path: Path) -> None:
    database = tmp_path / "cisbp_pwm.tsv"
    pd.DataFrame({
        "RBP_Name": ["RbpA"] * 3,
        "Motif_ID": ["M001"] * 3,
        "position": [1, 2, 3],
        "A": [9, 8, 7],
        "C": [1, 1, 1],
        "G": [0, 1, 1],
        "T": [0, 0, 1],
        "Database": ["CisBP-RNA"] * 3,
    }).to_csv(database, sep="\t", index=False)

    motifs = load_motif_database(database)

    assert len(motifs) == 1
    assert motifs[0].rbp == "RbpA"
    assert motifs[0].motif_id == "M001"
    assert motifs[0].source == "CisBP-RNA"
    assert motifs[0].pwm is not None
    np.testing.assert_allclose(motifs[0].pwm.sum(axis=1), np.ones(3))


def test_motif_event_windows_merge_overlapping_splice_site_windows() -> None:
    windows = _event_windows(pd.Series({"event_coordinates": "100|105|300"}), window_bp=10)

    assert windows == [(90, 115), (290, 310)]


def test_isoform_switch_uses_quant_gtf_and_inherits_dtu_support(tmp_path: Path) -> None:
    gtf = tmp_path / "annotation.gtf"
    gtf.write_text(
        'chr1\ttest\ttranscript\t1\t100\t.\t+\t.\tgene_id "ENSG1.1"; gene_name "Gene1"; transcript_id "TX1.1"; transcript_name "Tx1";\n'
        'chr1\ttest\ttranscript\t1\t100\t.\t+\t.\tgene_id "ENSG1.1"; gene_name "Gene1"; transcript_id "TX2.1"; transcript_name "Tx2";\n',
        encoding="utf-8",
    )
    quant_root = tmp_path / "quant"
    for sample, values in {"A1": (80, 20), "A2": (75, 25), "B1": (20, 80), "B2": (25, 75)}.items():
        directory = quant_root / sample
        directory.mkdir(parents=True)
        pd.DataFrame({
            "Name": ["TX1.1", "TX2.1"], "Length": [100, 100], "EffectiveLength": [90, 90],
            "TPM": values, "NumReads": values,
        }).to_csv(directory / "quant.sf", sep="\t", index=False)
    metadata = pd.DataFrame({"sample_id": ["A1", "A2", "B1", "B2"], "condition": ["A", "A", "B", "B"]})
    module_genes = pd.DataFrame({"module": ["M01"], "gene_id": ["ENSG1"], "DTU_significant": [True], "best_DTU_qvalue": [0.01]})

    switches, enrichment = quantify_isoform_switches(
        quant_root, gtf, metadata, module_genes, condition_column="condition",
        min_gene_tpm=1.0, min_abs_delta_if=0.1,
    )

    assert len(switches) == 1
    assert switches.iloc[0]["dominant_transcript_A"] == "TX1"
    assert switches.iloc[0]["dominant_transcript_B"] == "TX2"
    assert bool(switches.iloc[0]["DTU_significant"])
    assert enrichment.iloc[0]["overlap_switch_genes"] == 1
    assert set(enrichment["switch_scope"]) == {"descriptive_switch", "DTU_supported_switch"}


def test_beta_binomial_random_intercept_uses_repeated_donor_blocks() -> None:
    samples = pd.DataFrame({
        "sample_id": [f"D{donor}_{condition}" for donor in range(1, 5) for condition in ("A", "B")],
        "condition": [condition for _ in range(1, 5) for condition in ("A", "B")],
        "donor": [f"D{donor}" for donor in range(1, 5) for _ in ("A", "B")],
    })
    psi = pd.DataFrame(
        [[0.15, 0.72, 0.18, 0.75, 0.12, 0.68, 0.20, 0.78]],
        index=["E1"], columns=samples["sample_id"],
    )
    totals = pd.DataFrame(120.0, index=psi.index, columns=psi.columns)
    ijc = (psi * totals).round()
    sjc = totals - ijc
    assignments = pd.DataFrame({"event_key": ["E1"], "module": ["M01"]})
    scores = pd.DataFrame({"M01": [-1.0, 0.8, -0.9, 0.9, -1.1, 0.7, -0.8, 1.0]}, index=samples["sample_id"])

    result = fit_replicate_aware_models(
        psi, ijc, sjc, samples, assignments, scores,
        ModuleStatisticsConfig(
            min_total_reads=10, min_samples_per_condition=2,
            random_effect_column="donor", quadrature_points=5, optimizer_max_iterations=100,
        ),
    )

    assert len(result.event_models) == 1
    row = result.event_models.iloc[0]
    assert row["random_effect_column"] == "donor"
    assert np.isfinite(row["random_intercept_sd"])


def test_cross_cohort_meta_summarizes_multiple_external_cohorts(tmp_path: Path) -> None:
    samples = ["A1", "A2", "A3", "B1", "B2", "B3"]
    events = [f"E{index}" for index in range(8)]
    reference = pd.DataFrame(
        np.vstack([
            np.linspace(0.1, 0.8, 6) + index * 0.001 if index < 4
            else np.linspace(0.8, 0.1, 6) + index * 0.001
            for index in range(8)
        ]), index=events, columns=samples,
    )
    metadata = pd.DataFrame({"sample_id": samples, "condition": ["A"] * 3 + ["B"] * 3})
    assignments = pd.DataFrame({"event_key": events, "consensus_module": ["C01"] * 4 + ["C02"] * 4})
    cohorts = []
    for number, offset in ((1, 0.01), (2, -0.01)):
        psi_path = tmp_path / f"target_{number}.tsv"
        metadata_path = tmp_path / f"metadata_{number}.tsv"
        assignments_path = tmp_path / f"assignments_{number}.tsv"
        (reference + offset).to_csv(psi_path, sep="\t")
        metadata.to_csv(metadata_path, sep="\t", index=False)
        assignments.rename(columns={"consensus_module": "module"}).to_csv(assignments_path, sep="\t", index=False)
        cohorts.append(ExternalCohortConfig(
            f"external_{number}", psi_path, metadata_path, module_assignments=assignments_path
        ))

    result = meta_analyze_modules(
        reference, assignments, metadata, cohorts,
        MetaAnalysisConfig(min_shared_events=3, min_cohorts=2, min_supported_cohort_fraction=0.5),
        reference_name="discovery",
    )

    assert set(result.multi_cohort_consensus["reference_module"]) == {"C01", "C02"}
    assert (result.multi_cohort_consensus["n_external_cohorts_assessed"] == 2).all()
    assert (result.multi_cohort_consensus["multi_cohort_state"] == "preserved_across_cohorts").all()


def test_annotation_aware_alignment_matches_same_build_splice_sites() -> None:
    reference = pd.Index(["OldGene|SE|chr1|+|100|200|300"])
    target = pd.Index(["NewGene|SE|chr1|+|100|200|300"])

    result = annotation_aware_event_alignment(
        reference, target, "external", reference_genome_build="GRCm39", target_genome_build="GRCm39"
    )

    assert result.shared_events == reference.tolist()
    assert result.reference_to_target[reference[0]] == target[0]
    assert result.mapping.iloc[0]["matching_rule"].startswith("same-build exact splice-site")
    assert result.mapping.iloc[0]["alignment_confidence"] == "same_build_exact_splice_sites"


def test_annotation_aware_alignment_rejects_ambiguous_same_build_signature() -> None:
    reference = pd.Index(
        [
            "GeneA|SE|chr1|+|100|200|300",
            "GeneB|SE|chr1|+|100|200|300",
        ]
    )
    target = pd.Index(["GeneC|SE|chr1|+|100|200|300"])
    result = annotation_aware_event_alignment(
        reference,
        target,
        "external",
        reference_genome_build="GRCm39",
        target_genome_build="GRCm39",
    )
    assert result.shared_events == []
    assert result.mapping["alignment_confidence"].eq("unmatched").all()


def test_annotation_aware_alignment_requires_explicit_mapping_across_builds() -> None:
    reference = pd.Index(["OldGene|SE|chr1|+|100|200|300"])
    target = pd.Index(["NewGene|SE|chr1|+|110|210|310"])

    unmatched = annotation_aware_event_alignment(
        reference,
        target,
        "external",
        reference_genome_build="GRCm39",
        target_genome_build="GRCm38",
    )
    assert unmatched.shared_events == []

    mapping = pd.DataFrame({
        "reference_event_key": [reference[0]],
        "target_event_key": [target[0]],
    })
    matched = annotation_aware_event_alignment(
        reference,
        target,
        "external",
        reference_genome_build="GRCm39",
        target_genome_build="GRCm38",
        mapping_table=mapping,
    )
    assert matched.shared_events == reference.tolist()
    assert matched.mapping.iloc[0]["matching_rule"] == "explicit validated event/liftover mapping"


def test_cross_cohort_rewiring_requires_negative_significant_connectivity(tmp_path: Path) -> None:
    samples = ["A1", "A2", "A3", "B1", "B2", "B3"]
    events = [f"E{index}" for index in range(8)]
    reference = pd.DataFrame(
        np.vstack([
            np.linspace(0.1, 0.8, 6) + index * 0.001 if index < 4
            else np.linspace(0.8, 0.1, 6) + index * 0.001
            for index in range(8)
        ]),
        index=events,
        columns=samples,
    )
    metadata = pd.DataFrame({"sample_id": samples, "condition": ["A"] * 3 + ["B"] * 3})
    assignments = pd.DataFrame({"event_key": events, "consensus_module": ["C01"] * 4 + ["C02"] * 4})
    psi_path = tmp_path / "target.tsv"
    metadata_path = tmp_path / "metadata.tsv"
    assignments_path = tmp_path / "assignments.tsv"
    preservation_path = tmp_path / "preservation.tsv"
    reference.to_csv(psi_path, sep="\t")
    metadata.to_csv(metadata_path, sep="\t", index=False)
    assignments.rename(columns={"consensus_module": "module"}).to_csv(assignments_path, sep="\t", index=False)
    pd.DataFrame({
        "consensus_module": ["C01", "C02"],
        "connectivity_correlation": [-0.7, -0.7],
        "connectivity_empirical_pvalue": [0.001, 0.001],
        "connectivity_empirical_FDR": [0.01, 0.20],
    }).to_csv(preservation_path, sep="\t", index=False)
    cohort = ExternalCohortConfig(
        "external",
        psi_path,
        metadata_path,
        module_assignments=assignments_path,
        preservation_table=preservation_path,
    )

    result = meta_analyze_modules(
        reference,
        assignments,
        metadata,
        [cohort],
        MetaAnalysisConfig(min_shared_events=3, min_cohorts=1, min_supported_cohort_fraction=0.5),
    )

    states = result.module_state.set_index("reference_module")["module_state"].to_dict()
    assert states["C01"] == "rewired_supported_by_connectivity"
    assert states["C02"] == "consensus"
    consensus = result.multi_cohort_consensus.set_index("reference_module")["multi_cohort_state"].to_dict()
    assert consensus["C01"] == "rewired_across_cohorts"
    assert consensus["C02"] == "preserved_across_cohorts"
