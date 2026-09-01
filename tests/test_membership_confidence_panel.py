from __future__ import annotations

import pandas as pd

from katana_cosplice.membership_validation import (
    ARM_CORE,
    ARM_FULL,
    ARM_MATCHED,
    build_membership_arms,
)


def test_membership_arms_match_core_counts_per_module() -> None:
    assignments = pd.DataFrame(
        {
            "event_key": ["A", "B", "C", "D", "E", "F"],
            "consensus_module": ["M01", "M01", "M01", "M02", "M02", "grey"],
            "membership_role": ["core", "peripheral", "core", "peripheral", "core", "uncertain"],
            "full_data_kME": [0.70, 0.95, 0.80, 0.90, 0.60, 1.0],
        }
    )

    arms, audit = build_membership_arms(assignments)

    assert set(arms) == {ARM_FULL, ARM_MATCHED, ARM_CORE}
    assert set(arms[ARM_FULL]["event_key"]) == {"A", "B", "C", "D", "E"}
    assert set(arms[ARM_CORE]["event_key"]) == {"A", "C", "E"}
    assert set(arms[ARM_MATCHED]["event_key"]) == {"B", "C", "D"}
    assert audit[ARM_MATCHED].tolist() == audit[ARM_CORE].tolist()


def test_membership_arms_exclude_grey_even_when_core_labeled() -> None:
    assignments = pd.DataFrame(
        {
            "event_key": ["A", "B"],
            "consensus_module": ["M01", "grey"],
            "membership_role": ["core", "core"],
            "full_data_kME": [0.8, 0.9],
        }
    )

    arms, _audit = build_membership_arms(assignments)

    assert arms[ARM_CORE]["event_key"].tolist() == ["A"]
    assert arms[ARM_MATCHED]["event_key"].tolist() == ["A"]
