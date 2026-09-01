from __future__ import annotations

import pandas as pd

from katana_cosplice.formal_preservation import (
    prepare_formal_preservation_inputs,
)


def test_formal_preservation_keeps_common_background() -> None:
    reference = pd.DataFrame(
        {
            "s1": [0.1, 0.2, 0.3, 0.4],
            "s2": [0.2, 0.3, 0.5, 0.7],
            "s3": [0.3, 0.5, 0.7, 0.9],
        },
        index=["A", "B", "C", "D"],
    )
    target = pd.DataFrame(
        {
            "t1": [0.2, 0.1, 0.4, 0.8],
            "t2": [0.3, 0.2, 0.6, 0.6],
            "t3": [0.5, 0.4, 0.8, 0.3],
        },
        index=reference.index,
    )
    assignments = pd.DataFrame(
        {
            "event_key": ["A", "B"],
            "consensus_module": ["M01", "M01"],
        }
    )

    aligned_reference, aligned_target, labels = (
        prepare_formal_preservation_inputs(
            reference, target, assignments, min_module_size=2
        )
    )

    assert aligned_reference.index.tolist() == ["A", "B", "C", "D"]
    assert aligned_target.index.tolist() == ["A", "B", "C", "D"]
    assert labels["wgcna_module"].tolist() == ["M01", "M01", "grey", "grey"]
    assert labels["netrep_module"].tolist() == ["M01", "M01", "0", "0"]
