import numpy as np
import pandas as pd
import pytest

from katana_next.science.single_cell_benchmark import validate_psi_observations


@pytest.mark.parametrize("value", [-1, 1.5, 50, np.inf, "not_measured"])
def test_invalid_psi_is_not_silently_used(value):
    frame = pd.DataFrame({"cell_id": ["c"], "event_id": ["e"], "psi": [value]})
    with pytest.raises((ValueError, TypeError)):
        validate_psi_observations(frame, "psi", ["cell_id", "event_id"])


def test_missing_allowed_duplicates_rejected():
    frame = pd.DataFrame({"cell_id": ["c", "d"], "event_id": ["e", "e"], "psi": [.3, np.nan]})
    result = validate_psi_observations(frame, "psi", ["cell_id", "event_id"])
    assert pd.isna(result.psi.iloc[1])
    with pytest.raises(ValueError, match="Duplicate"):
        validate_psi_observations(pd.concat([frame, frame]), "psi", ["cell_id", "event_id"])
