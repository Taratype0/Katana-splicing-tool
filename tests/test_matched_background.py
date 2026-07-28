from __future__ import annotations

import pandas as pd

from katana_cosplice.matched_background import (
    build_matched_background_mut_coverage,
    build_significant_event_venn_summary,
)


def _events(rows: list[tuple[str, float, float]]) -> pd.DataFrame:
    return pd.DataFrame(rows, columns=["event_key", "FDR", "dPSI"])


def test_mut_signature_is_the_fixed_coverage_denominator() -> None:
    result = build_matched_background_mut_coverage(
        _events(
            [
                ("e1", 0.001, 0.20),
                ("e2", 0.002, -0.30),
                ("e3", 0.003, 0.15),
                ("not_sig", 0.50, 0.50),
            ]
        ),
        {
            "CLK3": _events([("e1", 0.001, 0.25), ("e2", 0.001, 0.31), ("x", 0.001, 0.50)]),
            "Pla": _events([("e1", 0.001, 0.18), ("e3", 0.20, 0.20), ("y", 0.001, -0.50)]),
        },
        context="DP",
    )

    summary = result.inhibitor_summary.set_index("inhibitor")
    assert summary.loc["CLK3", "n_MUT_significant"] == 3
    assert summary.loc["CLK3", "n_MUT_events_covered"] == 2
    assert summary.loc["CLK3", "n_same_direction"] == 1
    assert summary.loc["Pla", "n_MUT_events_covered"] == 1
    assert summary.loc["CLK3+Pla", "n_MUT_events_covered"] == 1

    classes = result.event_catalog.set_index("event_key")["coverage_class"].to_dict()
    assert classes == {
        "e1": "covered_by_both",
        "e2": "CLK3_only",
        "e3": "MUT_only",
    }


def test_inhibitor_dpsi_must_already_use_the_mut_canonical_direction() -> None:
    result = build_matched_background_mut_coverage(
        _events([("same", 0.001, 0.20), ("opposite", 0.001, -0.20)]),
        {"CLK3": _events([("same", 0.001, 0.30), ("opposite", 0.001, 0.30)])},
        context="noDP",
    )

    catalog = result.event_catalog.set_index("event_key")
    assert bool(catalog.loc["same", "CLK3_same_direction"])
    assert bool(catalog.loc["opposite", "CLK3_opposite_direction"])


def test_venn_summary_uses_all_three_complete_significant_event_sets() -> None:
    result = build_significant_event_venn_summary(
        {
            "MUT": _events(
                [
                    ("mut_only", 0.01, 0.2),
                    ("mut_clk", 0.01, 0.2),
                    ("all_three", 0.01, 0.2),
                    ("small_effect", 0.01, 0.05),
                ]
            ),
            "CLK3": _events(
                [
                    ("clk_only", 0.01, 0.2),
                    ("mut_clk", 0.01, -0.2),
                    ("all_three", 0.01, 0.2),
                ]
            ),
            "Pla": _events(
                [
                    ("pla_only", 0.01, 0.2),
                    ("all_three", 0.01, -0.2),
                ]
            ),
        },
        context="DP",
    )

    counts = result.region_summary.set_index("region")["n_events"].to_dict()
    assert counts["MUT_only"] == 1
    assert counts["CLK3_only"] == 1
    assert counts["PlaB_only"] == 1
    assert counts["MUT_CLK3_only"] == 1
    assert counts["MUT_CLK3_PlaB"] == 1
    assert len(result.membership_catalog) == 5
