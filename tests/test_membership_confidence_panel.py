from __future__ import annotations

import pandas as pd

from katana_cosplice.membership_validation import (
    ARM_CONFIDENCE,
    ARM_CORE,
    ARM_FULL,
    ARM_MATCHED,
    add_continuous_membership_confidence,
    build_membership_arms,
)


def _with_confidence_evidence(frame: pd.DataFrame) -> pd.DataFrame:
    frame = frame.copy()
    frame["consensus_membership_confidence"] = 0.8
    frame["consensus_confidence_margin"] = 0.6
    frame["bootstrap_q10_abs_kME"] = frame["full_data_kME"].abs().clip(0, 1)
    frame["bootstrap_kME_sign_consistency"] = 0.9
    frame["split_half_membership_role_a"] = frame["membership_role"]
    frame["split_half_membership_role_b"] = frame["membership_role"]
    return frame


def test_membership_arms_match_core_counts_per_module() -> None:
    assignments = _with_confidence_evidence(pd.DataFrame(
        {
            "event_key": ["A", "B", "C", "D", "E", "F"],
            "consensus_module": ["M01", "M01", "M01", "M02", "M02", "grey"],
            "membership_role": ["core", "peripheral", "core", "peripheral", "core", "uncertain"],
            "full_data_kME": [0.70, 0.95, 0.80, 0.90, 0.60, 1.0],
        }
    ))

    arms, audit = build_membership_arms(assignments)

    assert set(arms) == {
        ARM_FULL,
        ARM_MATCHED,
        ARM_CORE,
        ARM_CONFIDENCE,
    }
    assert set(arms[ARM_FULL]["event_key"]) == {"A", "B", "C", "D", "E"}
    assert set(arms[ARM_CORE]["event_key"]) == {"A", "C", "E"}
    assert set(arms[ARM_MATCHED]["event_key"]) == {"B", "C", "D"}
    assert len(arms[ARM_CONFIDENCE]) == len(arms[ARM_CORE])
    assert audit[ARM_MATCHED].tolist() == audit[ARM_CORE].tolist()


def test_membership_arms_exclude_grey_even_when_core_labeled() -> None:
    assignments = _with_confidence_evidence(pd.DataFrame(
        {
            "event_key": ["A", "B"],
            "consensus_module": ["M01", "grey"],
            "membership_role": ["core", "core"],
            "full_data_kME": [0.8, 0.9],
        }
    ))

    arms, _audit = build_membership_arms(assignments)

    assert arms[ARM_CORE]["event_key"].tolist() == ["A"]
    assert arms[ARM_MATCHED]["event_key"].tolist() == ["A"]


def test_continuous_confidence_uses_read_evidence_when_available() -> None:
    assignments = pd.DataFrame(
        {
            "event_key": ["A", "B"],
            "consensus_module": ["M01", "M01"],
            "consensus_membership_confidence": [0.8, 0.8],
            "consensus_confidence_margin": [0.7, 0.7],
            "bootstrap_q10_abs_kME": [0.7, 0.7],
            "bootstrap_kME_sign_consistency": [0.9, 0.9],
            "split_half_membership_role_a": ["core", "core"],
            "split_half_membership_role_b": ["core", "core"],
        }
    )
    reliability = pd.Series({"A": 0.95, "B": 0.20})

    scored = add_continuous_membership_confidence(
        assignments, reliability
    ).set_index("event_key")

    assert scored.loc["A", "continuous_membership_confidence"] > scored.loc[
        "B", "continuous_membership_confidence"
    ]
    assert scored["continuous_membership_confidence"].between(0, 1).all()
    assert scored["confidence_evidence_domains"].eq(4).all()
