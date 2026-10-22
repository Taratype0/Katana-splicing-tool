from pathlib import Path

import pandas as pd

from katana_cosplice.locked_ablation import (
    combine_gmt_files,
    functional_replication,
    load_gmt,
    local_module_enrichment,
    select_discovery_events,
)


def test_discovery_selection_is_variance_ranked_and_imputed() -> None:
    matrix = pd.DataFrame(
        {
            "s1": [0.1, 0.5, 0.2],
            "s2": [0.9, 0.5, None],
            "s3": [0.2, 0.5, 0.4],
            "s4": [0.8, 0.5, 0.6],
        },
        index=["A|ES|1", "B|ES|2", "C|ES|3"],
    )
    selected = select_discovery_events(matrix, 1)
    assert selected.index.tolist() == ["A|ES|1"]
    assert not selected.isna().any().any()


def test_local_functional_replication_uses_independent_modules(tmp_path: Path) -> None:
    first = tmp_path / "go.gmt"
    second = tmp_path / "kegg.gmt"
    first.write_text("splice\tdesc\ta\tb\n", encoding="utf-8")
    second.write_text("other\tdesc\tc\td\n", encoding="utf-8")
    merged = combine_gmt_files([first, second], tmp_path / "merged.gmt")
    terms = load_gmt(merged)
    discovery = pd.DataFrame(
        {"event_key": ["A|ES|1", "B|ES|2", "C|ES|3"], "consensus_module": ["C01", "C01", "grey"]}
    )
    holdout = pd.DataFrame(
        {"event_key": ["A|ES|1", "B|ES|2", "D|ES|4"], "consensus_module": ["C02", "C02", "grey"]}
    )
    enrichment = pd.concat(
        [
            local_module_enrichment("discovery", "method", discovery, terms),
            local_module_enrichment("holdout", "method", holdout, terms),
        ],
        ignore_index=True,
    )
    replication = functional_replication(enrichment, fdr=1.0)
    assert replication.loc[0, "best_holdout_module"] == "C02"
    assert replication.loc[0, "best_term_jaccard"] == 1.0
