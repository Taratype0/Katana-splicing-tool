import numpy as np
import pandas as pd
import pytest

from katana_next.science.embedding_comparison import compare_embeddings


def fixture():
    values = np.random.default_rng(17).normal(size=(40, 2))
    return pd.DataFrame({"cell_id": [f"c{i}" for i in range(40)], "umap_1": values[:, 0],
        "umap_2": values[:, 1], "cell_type": ["A"] * 20 + ["B"] * 20})


def test_comparison_aligns_cell_ids_and_is_rotation_invariant():
    a = fixture()
    b = a.copy()
    b["umap_1"], b["umap_2"] = a.umap_2 * 3 + 10, -a.umap_1 * 3 - 8
    summary, rows = compare_embeddings(a, b.sample(frac=1, random_state=2), k=5)
    assert summary["mean_neighbor_overlap_fraction"] == 1
    assert summary["mean_neighbor_jaccard"] == 1
    assert summary["pairwise_distance_spearman"] == pytest.approx(1)
    assert rows.shared_neighbors.eq(5).all()
    assert summary["chance_overlap_permutation_mean"] == pytest.approx(5 / 39, abs=.02)


@pytest.mark.parametrize("problem", ["duplicate", "metadata", "nan", "too_few"])
def test_comparison_rejects_invalid_pairs(problem):
    a, b = fixture(), fixture()
    if problem == "duplicate":
        b.loc[0, "cell_id"] = b.cell_id.iloc[1]
    elif problem == "metadata":
        b.loc[0, "cell_type"] = "different"
    elif problem == "nan":
        b.loc[0, "umap_1"] = np.nan
    else:
        b = b.iloc[:3]
    with pytest.raises(ValueError):
        compare_embeddings(a, b, k=5)


def test_pre_umap_comparison_uses_all_components_not_plot_axes():
    a = fixture().drop(columns=["umap_1", "umap_2"])
    values = np.random.default_rng(7).normal(size=(40, 4))
    for i in range(4):
        a[f"pc_{i + 1}"] = values[:, i]
    summary, _ = compare_embeddings(a, a, k=5)
    assert summary["mean_neighbor_overlap_fraction"] == 1


def test_paired_plot_rejects_reversed_data_sources(tmp_path):
    import importlib.util
    from pathlib import Path
    path = Path(__file__).resolve().parents[2] / "scripts" / "demo_paired_cell_embeddings.py"
    spec = importlib.util.spec_from_file_location("paired_demo", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    a, b = fixture(), fixture()
    a["embedding_source"], b["embedding_source"] = "splicing", "expression"
    with pytest.raises(ValueError, match="titles"):
        module.paired_plot(a, b, "cell_type", tmp_path / "must_not_write")
    assert not (tmp_path / "must_not_write.png").exists()
