from __future__ import annotations

import importlib.util
from pathlib import Path

import pandas as pd


def _module():
    script = (
        Path(__file__).parents[1]
        / "scripts"
        / "benchmark_three_public_spliceseq.py"
    )
    spec = importlib.util.spec_from_file_location(
        "benchmark_three_public_spliceseq",
        script,
    )
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_selection_summary_returns_one_recommended_method() -> None:
    module = _module()
    internal = pd.DataFrame(
        [
            {
                "method": "CoSplice",
                "valid_partition": True,
                "min_bootstrap_ARI": 0.7,
                "silhouette_cosine": 0.5,
            },
            {
                "method": "Exact WGCNA",
                "valid_partition": False,
                "min_bootstrap_ARI": 0.9,
                "silhouette_cosine": 0.6,
            },
        ]
    )
    external = pd.DataFrame(
        [
            {"method": "CoSplice", "best_module_Zsummary": 2.5},
            {"method": "Exact WGCNA", "best_module_Zsummary": 3.0},
        ]
    )
    result = module._selection_summary(internal, external)
    assert result["recommended_production_method"].sum() == 1
    assert result.iloc[0]["method"] == "CoSplice"


def test_external_matrix_can_keep_full_eligible_event_catalogue(
    tmp_path: Path,
) -> None:
    module = _module()
    path = tmp_path / "psi.tsv"
    pd.DataFrame(
        {
            "s1": [0.1, 0.2, 0.3],
            "s2": [0.2, 0.4, 0.3],
            "s3": [0.3, 0.6, 0.3],
            "s4": [0.4, 0.8, 0.3],
        },
        index=["event_a", "event_b", "constant"],
    ).to_csv(path, sep="\t")

    discovery = module._load_matrix(path, 1)
    external = module._load_matrix(path, None)

    assert len(discovery) == 1
    assert external.index.tolist() == ["event_a", "event_b"]
