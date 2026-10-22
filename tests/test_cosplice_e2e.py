from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from katana_cosplice.config import (
    ConsensusConfig,
    CoSpliceConfig,
    ExternalCohortConfig,
    GenePriorConfig,
    HybridProgramConfig,
    MultiObjectiveSelectionConfig,
    NetworkConfig,
    StableASConfig,
    ValidationConfig,
)
from katana_cosplice.matrix import EVENT_COORDINATES, EVENT_TYPES
from katana_cosplice.pipeline import CoSplicePipeline
from katana_cosplice.stage_registry import stage_specs


def _write_pair_tables(root: Path, left: str, right: str, condition_index: dict[str, int]) -> None:
    post = root / f"{left}-{right}" / "rmats_post"
    post.mkdir(parents=True)
    jitter = np.asarray([-0.01, 0.0, 0.01])
    for type_index, event_type in enumerate(EVENT_TYPES):
        rows: list[dict[str, object]] = []
        for local_index in range(12):
            event_index = type_index * 12 + local_index
            increasing = event_index < 30

            def replicate_values(condition: str) -> np.ndarray:
                mean = 0.15 + 0.14 * condition_index[condition]
                if not increasing:
                    mean = 1.0 - mean
                return np.clip(mean + jitter, 0.02, 0.98)

            left_psi = replicate_values(left)
            right_psi = replicate_values(right)
            row: dict[str, object] = {
                "GeneID": f'"ENSMUSG{event_index:011d}.1"',
                "geneSymbol": f'"Gene{event_index:03d}"',
                "chr": "chr1",
                "strand": "+" if event_index % 2 == 0 else "-",
                "IncLevel1": ",".join(f"{value:.3f}" for value in left_psi),
                "IncLevel2": ",".join(f"{value:.3f}" for value in right_psi),
                "IJC_SAMPLE_1": ",".join(str(int(round(value * 120))) for value in left_psi),
                "SJC_SAMPLE_1": ",".join(str(int(round((1.0 - value) * 120))) for value in left_psi),
                "IJC_SAMPLE_2": ",".join(str(int(round(value * 120))) for value in right_psi),
                "SJC_SAMPLE_2": ",".join(str(int(round((1.0 - value) * 120))) for value in right_psi),
            }
            for offset, column in enumerate(EVENT_COORDINATES[event_type]):
                row[column] = 100_000 + event_index * 1_000 + offset * 100
            rows.append(row)
        pd.DataFrame(rows).to_csv(post / f"{event_type}.MATS.JCEC.txt", sep="\t", index=False)


def _full_project(tmp_path: Path) -> CoSpliceConfig:
    rmats = tmp_path / "rmats"
    bamlist = rmats / "bamlist"
    bamlist.mkdir(parents=True)
    conditions = list("ABCDEF")
    condition_index = {condition: index for index, condition in enumerate(conditions)}
    metadata_rows: list[dict[str, object]] = []
    for condition in conditions:
        sample_ids = [f"{condition}_{replicate}" for replicate in range(1, 4)]
        (bamlist / f"{condition}_bamlist.txt").write_text(
            ",".join(f"{sample}_sorted.bam" for sample in sample_ids),
            encoding="utf-8",
        )
        metadata_rows.extend(
            {
                "sample_id": sample,
                "condition": condition,
                "batch": str(replicate),
            }
            for replicate, sample in enumerate(sample_ids, start=1)
        )
    for left, right in (("A", "B"), ("C", "D"), ("E", "F")):
        _write_pair_tables(rmats, left, right, condition_index)
    metadata_path = tmp_path / "samples.tsv"
    pd.DataFrame(metadata_rows).to_csv(metadata_path, sep="\t", index=False)
    return CoSpliceConfig(
        project_name="synthetic_full",
        rmats_root=rmats,
        sample_metadata=metadata_path,
        output_dir=tmp_path / "output",
        residual_covariates=["batch"],
        association_factors=["condition", "batch"],
        factorial_factors=["condition"],
        batch_column="batch",
        stable_as=StableASConfig(),
        network=NetworkConfig(
            module_detection="dynamic_tree_cut",
            min_module_size=8,
            knn_neighbors=10,
            bootstrap_iterations=2,
            max_output_edges=5_000,
            random_seed=11,
        ),
        validation=ValidationConfig(
            bootstrap_iterations=2,
            knn_neighbors=[10],
            community_resolutions=[1.0],
            random_seed=13,
        ),
        consensus=ConsensusConfig(
            bootstrap_iterations=3,
            min_coassignment=0.50,
            min_module_size=8,
            dense_max_events=200,
            random_seed=17,
        ),
        hybrid_programs=HybridProgramConfig(
            diffuse_engine="python_soft_tom",
            bootstrap_iterations=2,
            min_module_size=8,
            python_tom_cut_height=0.90,
            python_max_events=200,
            random_seed=19,
        ),
    )


def _write_gene_prior_exon_counts(
    path: Path,
    metadata_path: Path,
    n_genes: int = 60,
) -> None:
    samples = pd.read_csv(metadata_path, sep="\t")
    condition_codes = {
        condition: index
        for index, condition in enumerate(
            samples["condition"].drop_duplicates().astype(str)
        )
    }
    sample_signal = np.asarray(
        [
            condition_codes[str(condition)] + int(str(sample_id).rsplit("_", 1)[-1]) / 10
            for sample_id, condition in zip(
                samples["sample_id"],
                samples["condition"],
                strict=True,
            )
        ],
        dtype=float,
    )
    rows: dict[str, np.ndarray] = {}
    for gene_index in range(n_genes):
        gene_id = f"ENSMUSG{gene_index:011d}.1"
        program = gene_index % 3
        for exon_index in range(3):
            scale = 8.0 + program * 2.0 + exon_index
            phase = (gene_index % 5) * 0.15 + exon_index * 0.35
            rows[f"{gene_id}_exon{exon_index + 1}"] = (
                100.0
                + gene_index
                + scale * sample_signal
                + 4.0 * np.sin(sample_signal + phase)
            )
    pd.DataFrame(
        rows,
        index=samples["sample_id"].astype(str),
    ).T.to_csv(path, sep="\t")


def test_gene_prior_runs_through_network_consensus_and_cache(tmp_path: Path) -> None:
    config = _full_project(tmp_path)
    exon_counts = tmp_path / "exon_counts.tsv"
    _write_gene_prior_exon_counts(exon_counts, config.sample_metadata)
    config.gene_prior = GenePriorConfig(
        enabled=True,
        exon_count_matrix=exon_counts,
        weight=0.10,
        scaffold_knn_neighbors=8,
        scaffold_min_module_size=4,
        random_seed=23,
    )
    config.validate()
    pipeline = CoSplicePipeline(config)

    pipeline.run_stage("audit")
    pipeline.run_stage("matrix")
    pipeline.run_stage("stable_qc")
    pipeline.run_stage("network_observed")

    network_dir = config.output_dir / "04_network" / "observed"
    assert (network_dir / "gene_mantel_prior_matrix.npz").exists()
    assert (network_dir / "gene_mantel_prior_gene_audit.tsv").exists()
    assert (network_dir / "gene_mantel_prior_event_audit.tsv").exists()
    edges = pd.read_csv(network_dir / "network_edges.tsv", sep="\t")
    assert {
        "base_AS_adjacency",
        "gene_mantel_prior",
        "gene_prior_adjacency_boost",
    }.issubset(edges)
    assert edges["gene_prior_adjacency_boost"].ge(0).all()

    pipeline.run_stage("consensus_modules")
    consensus_dir = config.output_dir / "04c_consensus_modules"
    assert (consensus_dir / "consensus_event_module_assignments.tsv").exists()
    assert (consensus_dir / "bootstrap_coassignment_matrix.npz").exists()

    pipeline.run_stage("network_observed")
    assert pipeline.status()["stages"]["network_observed"]["status"] == "cache_hit"
    pipeline.run_stage("consensus_modules")
    assert pipeline.status()["stages"]["consensus_modules"]["status"] == "cache_hit"


def test_multiobjective_consensus_pipeline_uses_tuning_not_holdout(
    tmp_path: Path,
) -> None:
    config = _full_project(tmp_path)
    tuning_matrix = tmp_path / "tuning_psi.tsv"
    tuning_samples = tmp_path / "tuning_samples.tsv"
    tuning_gmt = tmp_path / "tuning_terms.gmt"
    config.external_cohorts = [
        ExternalCohortConfig("tuning", tuning_matrix, tuning_samples),
        ExternalCohortConfig(
            "untouched_holdout",
            tmp_path / "missing_holdout.tsv",
            tmp_path / "missing_holdout_samples.tsv",
        ),
    ]
    config.multiobjective_selection = MultiObjectiveSelectionConfig(
        enabled=True,
        tuning_cohorts=["tuning"],
        holdout_cohorts=["untouched_holdout"],
        functional_gene_sets_gmt=tuning_gmt,
        min_shared_events_per_module=3,
        min_module_genes=2,
    )
    config.validate()
    pipeline = CoSplicePipeline(config)
    for stage in ("audit", "matrix", "stable_qc", "network_observed"):
        pipeline.run_stage(stage)

    reference = pd.read_csv(
        config.output_dir
        / "04_network"
        / "observed"
        / "selected_event_matrix.tsv",
        sep="\t",
        index_col=0,
    )
    rng = np.random.default_rng(991)
    target = reference.copy()
    target.columns = [f"T{index:02d}" for index in range(target.shape[1])]
    target += rng.normal(0, 0.002, target.shape)
    target.to_csv(tuning_matrix, sep="\t")
    pd.DataFrame(
        {
            "sample_id": target.columns,
            "condition": ["A"] * (len(target.columns) // 2)
            + ["B"] * (len(target.columns) - len(target.columns) // 2),
        }
    ).to_csv(tuning_samples, sep="\t", index=False)
    event_metadata = pd.read_csv(
        config.output_dir / "02_matrix" / "event_metadata.tsv",
        sep="\t",
    )
    genes = (
        event_metadata["gene_symbol"]
        .dropna()
        .astype(str)
        .drop_duplicates()
        .tolist()
    )
    midpoint = max(2, len(genes) // 2)
    tuning_gmt.write_text(
        "TERM_A\tdescription\t"
        + "\t".join(genes[:midpoint])
        + "\nTERM_B\tdescription\t"
        + "\t".join(genes[midpoint:])
        + "\n",
        encoding="utf-8",
    )

    pipeline.run_stage("consensus_modules")
    out_dir = config.output_dir / "04c_consensus_modules"
    search = pd.read_csv(
        out_dir / "consensus_parameter_search.tsv", sep="\t"
    )
    audit = pd.read_csv(
        out_dir / "multiobjective_selection_audit.tsv", sep="\t"
    )
    assert {
        "objective_internal_stability",
        "objective_external_preservation",
        "objective_functional_coherence",
        "multiobjective_selection_score",
        "pareto_optimal",
    }.issubset(search)
    assert search["selected"].sum() == 1
    assert not audit.empty
    assert set(
        audit.loc[
            audit["audit_type"].eq("external_preservation"), "cohort"
        ]
    ) == {"tuning"}
    assert (out_dir / "plots" / "multiobjective_candidate_tradeoff.png").exists()

    pipeline.run_stage("consensus_modules")
    assert (
        pipeline.status()["stages"]["consensus_modules"]["status"]
        == "cache_hit"
    )


def test_full_explicit_pipeline_produces_contract_artifacts(tmp_path: Path) -> None:
    config = _full_project(tmp_path)
    config.validate()
    pipeline = CoSplicePipeline(config)
    pipeline.run_all(include_residual=True)

    state = pipeline.status()["stages"]
    expected = {
        spec.key
        for spec in stage_specs(visible_only=True)
        if (
            spec.key != "factor_contrasts"
            or config.factor_contrasts.enabled
        )
        and (
            spec.key != "multiview_fusion"
            or config.multiview.enabled
        )
    }
    assert set(state) == expected
    assert all(details["status"] == "finished" for details in state.values())

    psi = pd.read_csv(config.output_dir / "02_matrix" / "sample_level_psi_matrix.tsv", sep="\t", index_col=0)
    stable = pd.read_csv(config.output_dir / "03_stable_as_qc" / "stable_psi_matrix.tsv", sep="\t", index_col=0)
    condition_qc = pd.read_csv(config.output_dir / "03_stable_as_qc" / "stable_as_condition_qc.tsv", sep="\t")
    assignments = pd.read_csv(config.output_dir / "04_network" / "observed" / "event_module_assignments.tsv", sep="\t")
    association = pd.read_csv(config.output_dir / "05_association" / "observed" / "module_trait_association.tsv", sep="\t")
    association_design = pd.read_csv(config.output_dir / "05_association" / "observed" / "association_design.tsv", sep="\t")
    validation = pd.read_csv(config.output_dir / "06_module_validation" / "module_bootstrap_stability.tsv", sep="\t")
    consensus = pd.read_csv(
        config.output_dir / "04c_consensus_modules" / "consensus_event_module_assignments.tsv",
        sep="\t",
    )

    assert psi.shape == (60, 18)
    assert stable.shape == (60, 18)
    assert len(condition_qc) == 60 * 6
    assert assignments["event_key"].is_unique
    assert assignments["module"].ne("grey").any()
    assert not association.empty
    assert association_design.columns.is_unique
    assert list(association_design.columns).count("batch") == 1
    assert not validation.empty
    assert consensus["event_key"].is_unique
    assert (config.output_dir / "04c_consensus_modules" / "bootstrap_coassignment_matrix.npz").exists()
    assert (config.output_dir / "04c_consensus_modules" / "consensus_parameter_search.tsv").exists()
    assert not (config.output_dir / "04d_hybrid_programs").exists()
    assert (config.output_dir / "07_report" / "katana_cosplice_report.html").exists()

    manifest = json.loads((config.output_dir / "artifact_manifest.json").read_text(encoding="utf-8"))
    assert set(manifest["stages"]) == set(state)
    assert all(stage["explicit_only"] for stage in manifest["stages"].values())
    assert all(
        artifact["exists"]
        for stage in manifest["stages"].values()
        for artifact in stage["generated_artifacts"]
    )

    pipeline.run_stage("matrix")
    assert pipeline.status()["stages"]["matrix"]["status"] == "cache_hit"
    assert pipeline.status()["stages"]["stable_qc"]["status"] == "finished"

    # The historical hybrid stage remains explicitly callable for benchmark
    # compatibility, but it is no longer part of the production run-all path.
    pipeline.run_stage("hybrid_programs")
    assert pipeline.status()["stages"]["hybrid_programs"]["status"] == "finished"
    assert pipeline.status()["stages"]["network_observed"]["status"] == "finished"

    hybrid_events = pd.read_csv(
        config.output_dir
        / "04d_hybrid_programs"
        / "hybrid_event_program_assignments.tsv",
        sep="\t",
    )
    hybrid_programs = pd.read_csv(
        config.output_dir / "04d_hybrid_programs" / "hybrid_program_summary.tsv",
        sep="\t",
    )
    assert hybrid_events["event_key"].is_unique
    assert {"tight_program", "diffuse_program", "program_scope"}.issubset(
        hybrid_events
    )
    assert set(hybrid_programs["program_type"]) == {"tight", "diffuse"}
    assert hybrid_programs["overall_confidence"].between(0, 1).all()
    assert (
        config.output_dir
        / "04d_hybrid_programs"
        / "plots"
        / "tight_diffuse_program_confidence.png"
    ).exists()

    pipeline.run_stage("hybrid_programs")
    assert pipeline.status()["stages"]["hybrid_programs"]["status"] == "cache_hit"
