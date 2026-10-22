from __future__ import annotations

import importlib.util
from pathlib import Path

import pandas as pd


def _load_benchmark_module():
    script = Path(__file__).parents[1] / "scripts" / "benchmark_public_spliceseq.py"
    spec = importlib.util.spec_from_file_location("benchmark_public_spliceseq", script)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_public_benchmark_outputs_use_configured_cohort_names(tmp_path: Path) -> None:
    module = _load_benchmark_module()
    internal = pd.DataFrame(
        [
            {"method": "CoSplice", "mean_bootstrap_ARI": 0.8, "min_bootstrap_ARI": 0.7},
            {"method": "Exact WGCNA", "mean_bootstrap_ARI": 0.9, "min_bootstrap_ARI": 0.8},
        ]
    )
    external = pd.DataFrame(
        [
            {"method": "CoSplice", "median_Zsummary": 2.1},
            {"method": "Exact WGCNA", "median_Zsummary": 0.6},
        ]
    )
    plot = tmp_path / "benchmark.png"
    module._dashboard(
        internal,
        external,
        plot,
        discovery_name="Discovery-Cohort",
        external_name="Validation-Cohort",
    )
    report = module._write_report(
        tmp_path,
        internal,
        external,
        discovery_name="Discovery-Cohort",
        external_name="Validation-Cohort",
        parameters={
            "n_discovery_events": 150,
            "n_discovery_samples": 30,
            "n_external_samples": 24,
            "n_exact_shared_events": 140,
            "bootstraps": 20,
            "permutations_per_module": 1000,
        },
    )

    assert plot.exists() and plot.stat().st_size > 0
    text = report.read_text(encoding="utf-8")
    assert "Discovery-Cohort" in text
    assert "Validation-Cohort" in text
    assert "not the exact `WGCNA::modulePreservation` statistic" in text
