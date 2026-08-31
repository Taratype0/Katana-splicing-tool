from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd


def _module():
    script = (
        Path(__file__).parents[1]
        / "scripts"
        / "benchmark_three_independent_public_datasets.py"
    )
    spec = importlib.util.spec_from_file_location(
        "benchmark_three_independent_public_datasets",
        script,
    )
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_select_event_universe_uses_one_shared_variance_ranking() -> None:
    module = _module()
    matrix = pd.DataFrame(
        {
            "S1": [0.1, 0.2, np.nan, 0.5],
            "S2": [0.9, 0.2, 0.4, 0.5],
            "S3": [0.2, 0.2, np.nan, 0.5],
            "S4": [0.8, 0.2, 0.7, 0.5],
        },
        index=["variable_a", "constant", "variable_b", "constant_b"],
    )
    selected = module._select_event_universe(matrix, 2)
    assert list(selected.index) == ["variable_a", "variable_b"]
    assert selected.notna().all().all()
