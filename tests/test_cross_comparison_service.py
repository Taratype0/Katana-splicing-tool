from dataclasses import dataclass

import pandas as pd

from src.services.cross_comparison_service import CrossComparisonService


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
