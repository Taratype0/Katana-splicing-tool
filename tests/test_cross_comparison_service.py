from dataclasses import dataclass

import pandas as pd

from src.services.cross_comparison_service import CrossComparisonService
from src.services.cross_interaction_mixin import CrossInteractionMixin


@dataclass
class _Pair:
    pair_id: str = "a_vs_b"
    resolved_name: str = "A compared with B"
    comparison_a: str = "a"
    comparison_b: str = "b"


def test_cross_comparison_service_writes_pattern_contract(tmp_path) -> None:
    service = CrossComparisonService(tmp_path)
    tables = {
        "same_direction_large_delta": pd.DataFrame({"event": ["e1"]}),
        "opposite_direction": pd.DataFrame({"event": ["e2"]}),
        "all_shared": pd.DataFrame({"event": ["e1", "e2"]}),
        "condition_specific": pd.DataFrame({"event": ["e3"]}),
    }
    generated = service.write_pattern_outputs(
        pair=_Pair(),
        tables=tables,
        abs_cutoff=0.1,
        delta_cutoff=0.2,
    )
    assert len(generated) == 5
    assert all(path.exists() for path in generated)
    summary = pd.read_csv(service.pattern_files("a_vs_b")["summary"], sep="\t")
    assert summary.iloc[0]["n_all_shared"] == 2
    assert summary.iloc[0]["comparison_A"] == "a"


def test_cross_comparison_service_writes_interaction_parameters_atomically(tmp_path) -> None:
    service = CrossComparisonService(tmp_path)
    generated = service.write_interaction_outputs(
        pair=_Pair(),
        tables={"unified_event_table": pd.DataFrame({"event": ["e1"]})},
        abs_cutoff=0.1,
        delta_cutoff=0.15,
        min_total_reads=10,
    )
    parameters = service.interaction_dir("a_vs_b") / "parameters.json"
    assert parameters in generated
    assert parameters.exists()
    assert not parameters.with_suffix(".json.tmp").exists()


def test_cross_comparison_service_builds_directional_pattern_tables(tmp_path) -> None:
    service = CrossComparisonService(tmp_path)
    columns = {
        "comparison_display_name": ["A", "A", "A"],
        "gene": ["G1", "G2", "G3"],
        "event_type": ["SE", "RI", "A3SS"],
        "event_id": ["e1", "e2", "e3"],
        "event_matching_key": ["G1|SE|e1", "G2|RI|e2", "G3|A3SS|e3"],
        "dPSI": [0.3, 0.2, -0.25],
        "FDR": [0.01, 0.02, 0.03],
        "significant": [True, "true", True],
        "direction": ["up", "up", "down"],
    }
    comparison_a = pd.DataFrame(columns)
    comparison_b = pd.DataFrame(
        {
            **columns,
            "comparison_display_name": ["B", "B", "B"],
            "dPSI": [0.1, -0.3, -0.1],
            "significant": [True, True, False],
            "direction": ["up", "down", "down"],
        }
    )

    tables = service.build_pattern_tables(
        pair=_Pair(),
        comparison_a=comparison_a,
        comparison_b=comparison_b,
        abs_dpsi_cutoff=0.1,
        large_delta_dpsi_cutoff=0.15,
    )

    assert set(tables["all_shared"]["event_id"]) == {"e1", "e2"}
    assert tables["same_direction_large_delta"]["event_id"].tolist() == ["e1"]
    assert tables["opposite_direction"]["event_id"].tolist() == ["e2"]
    assert tables["condition_specific"]["event_id"].tolist() == ["e3"]
    assert tables["condition_specific"]["event_class"].tolist() == ["comparison_A_only"]


def test_cross_interaction_mixin_builds_empty_single_and_batch_contracts() -> None:
    single = CrossInteractionMixin._empty_cross_as_interaction_tables()
    batch = CrossInteractionMixin._empty_cross_as_interaction_batch_tables()

    assert {"shared_interaction", "exclusive", "benchmark", "diagnostics"} <= set(single)
    assert list(batch["batch_top_interactions"].columns) == list(single["shared_interaction"].columns)
    assert list(batch["batch_exclusive"].columns) == list(single["exclusive"].columns)
    assert list(batch["batch_benchmark"].columns) == list(single["benchmark"].columns)
    assert list(batch["batch_diagnostics"].columns) == list(single["diagnostics"].columns)
