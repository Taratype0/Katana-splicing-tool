from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from katana_cosplice.config import GenePriorConfig
from katana_cosplice.gene_prior import (
    build_gene_mantel_prior,
    load_gene_mantel_prior,
    save_gene_mantel_prior,
)
from katana_cosplice.network import _reinforce_with_gene_prior


def _write_exon_counts(path: Path) -> tuple[pd.Index, pd.Index]:
    samples = pd.Index([f"S{index}" for index in range(1, 7)])
    counts = pd.DataFrame(
        [
            [10, 12, 15, 25, 28, 30],
            [30, 28, 25, 15, 12, 10],
            [8, 10, 12, 20, 23, 25],
            [25, 23, 20, 12, 10, 8],
        ],
        index=[
            "ENSG000001_exon1",
            "ENSG000001_exon2",
            "ENSG000002_exon1",
            "ENSG000002_exon2",
        ],
        columns=samples,
    )
    counts.to_csv(path, sep="\t")
    events = pd.Index(
        [
            "ENSG000001|SE|event1",
            "ENSG000001|SE|event2",
            "ENSG000002|RI|event3",
        ]
    )
    return events, samples


def test_gene_prior_build_save_load_and_bounded_reinforcement(tmp_path: Path) -> None:
    counts_path = tmp_path / "exon_counts.tsv"
    events, samples = _write_exon_counts(counts_path)
    settings = GenePriorConfig(enabled=True, exon_count_matrix=counts_path, weight=0.2)
    prior = build_gene_mantel_prior(events, samples, settings)

    assert list(prior.gene_adjacency.index) == ["ENSG000001", "ENSG000002"]
    assert prior.event_audit["gene_prior_available"].all()

    matrix_path = tmp_path / "prior.npz"
    genes_path = tmp_path / "genes.tsv"
    event_path = tmp_path / "events.tsv"
    save_gene_mantel_prior(prior, matrix_path, genes_path, event_path)
    loaded = load_gene_mantel_prior(matrix_path, events)

    base = np.asarray([0.20, 0.40, 0.0])
    rows = np.asarray([0, 0, 1])
    cols = np.asarray([1, 2, 2])
    reinforced, values = _reinforce_with_gene_prior(
        base,
        events,
        rows,
        cols,
        loaded,
        settings.weight,
    )
    assert values[0] == 1.0
    assert reinforced[0] > base[0]
    assert reinforced[2] == 0.0
    assert np.all(reinforced <= base + settings.weight * base * (1.0 - base) + 1e-12)


def test_gene_prior_uses_explicit_feature_mapping(tmp_path: Path) -> None:
    samples = pd.Index([f"S{index}" for index in range(1, 7)])
    counts = pd.DataFrame(
        np.arange(24, dtype=float).reshape(4, 6) + 1,
        index=["feature_a", "feature_b", "feature_c", "feature_d"],
        columns=samples,
    )
    counts_path = tmp_path / "counts.tsv"
    counts.to_csv(counts_path, sep="\t")
    mapping_path = tmp_path / "mapping.tsv"
    pd.DataFrame(
        {
            "feature_id": counts.index,
            "gene_id": ["GENE1", "GENE1", "GENE2", "GENE2"],
        }
    ).to_csv(mapping_path, sep="\t", index=False)
    prior = build_gene_mantel_prior(
        pd.Index(["GENE1|SE|a", "GENE2|SE|b"]),
        samples,
        GenePriorConfig(
            enabled=True,
            exon_count_matrix=counts_path,
            feature_gene_map=mapping_path,
        ),
    )
    assert set(prior.gene_adjacency.index) == {"GENE1", "GENE2"}
