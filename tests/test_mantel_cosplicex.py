from __future__ import annotations

import numpy as np
import pandas as pd

from katana_cosplice.mantel_cosplicex import (
    canberra_distance_vectors,
    exon_gene_id,
    mantel_similarity,
    select_gene_profiles,
    unsigned_adjacency,
)


def test_exon_gene_id_uses_last_underscore() -> None:
    assert exon_gene_id("ENSMUSG000001_12345") == "ENSMUSG000001"


def test_select_gene_profiles_keeps_unique_multi_exon_genes() -> None:
    counts = pd.DataFrame(
        {
            "S1": [10, 5, 3, 4, 8],
            "S2": [12, 4, 4, 3, 9],
            "S3": [8, 7, 5, 2, 10],
            "S4": [9, 6, 2, 5, 11],
        },
        index=["G1_1", "G1_2", "G2_1", "G2_2", "G3_1"],
    )
    profiles, audit = select_gene_profiles(counts, min_exons=2)

    assert list(profiles) == audit["gene_id"].tolist()
    assert set(profiles) == {"G1", "G2"}
    assert audit["gene_id"].is_unique


def test_mantel_similarity_matches_correlated_distance_profiles() -> None:
    counts = pd.DataFrame(
        {
            "S1": [10, 1, 20, 2],
            "S2": [8, 3, 16, 6],
            "S3": [4, 7, 8, 14],
            "S4": [1, 10, 2, 20],
        },
        index=["G1_1", "G1_2", "G2_1", "G2_2"],
    )
    profiles, _ = select_gene_profiles(counts)
    vectors = canberra_distance_vectors(profiles)
    similarity = mantel_similarity(vectors)
    adjacency = unsigned_adjacency(similarity, soft_power=6)

    assert np.isclose(similarity.loc["G1", "G2"], 1.0)
    assert np.isclose(adjacency.loc["G1", "G2"], 1.0)
    assert np.allclose(np.diag(adjacency), 1.0)
