"""Rendering must preserve the producer's exact recorded float64 Hub order."""
import numpy as np
import pandas as pd
import pytest

from katana_next.reports.network_views import _read_observed_kme_hubs


def _write_family(tmp_path):
    # These are distinct adjacent float64 values. Pandas' default CSV parser
    # rounds both to the lower value, inventing a tie that reverses e1/e2.
    lower = np.float64(0.9238421084414568)
    higher = np.nextafter(lower, np.inf)
    frame = pd.DataFrame({
        "event_key": ["e1", "e2", "e3", "e4", "e5", "e6"],
        "module": ["M01", "M01", "M02", "M02", "M03", "M03"],
        "module_membership_kME": [lower, -higher, .8, -.8, np.nan, np.nan],
        "hub_event_score": [lower, higher, .8, .8, np.nan, np.nan],
        "hub_event_rank": [2, 1, 1, 2, 1, 2],
        "hub_metric": ["abs_kME"] * 6,
    })
    path = tmp_path / "hub-events.tsv"
    frame.to_csv(path, sep="\t", index=False)
    members = pd.concat([frame[["event_key", "module"]],
                         pd.DataFrame({"event_key": ["unassigned"], "module": ["grey"]})])
    return path, frame, members


def test_exact_float64_roundtrip_preserves_near_neighbor_order_and_real_ties(tmp_path):
    path, original, members = _write_family(tmp_path)
    default = pd.read_csv(path, sep="\t")
    assert default.hub_event_score.iloc[0] == default.hub_event_score.iloc[1]
    assert original.hub_event_score.iloc[0] < original.hub_event_score.iloc[1]
    before = path.read_bytes()
    loaded = _read_observed_kme_hubs(path, members)
    np.testing.assert_array_equal(loaded.module_membership_kME, original.module_membership_kME)
    np.testing.assert_array_equal(loaded.hub_event_score, original.hub_event_score)
    assert loaded.hub_event_rank.tolist() == [2, 1, 1, 2, 1, 2]
    assert path.read_bytes() == before


@pytest.mark.parametrize("indices", [(0, 1), (2, 3), (4, 5)])
def test_corrupt_rank_rejected_for_near_neighbors_exact_ties_and_missing_values(tmp_path, indices):
    path, frame, members = _write_family(tmp_path)
    a, b = indices
    frame.loc[[a, b], "hub_event_rank"] = frame.loc[[b, a], "hub_event_rank"].to_numpy()
    frame.to_csv(path, sep="\t", index=False)
    with pytest.raises(ValueError, match="ranks differ"):
        _read_observed_kme_hubs(path, members)


@pytest.mark.parametrize("problem", ["score", "membership", "module", "missing_event", "duplicate_event"])
def test_rank_parser_preserves_strict_family_contract(tmp_path, problem):
    path, frame, members = _write_family(tmp_path)
    if problem == "score":
        frame.loc[0, "hub_event_score"] = .5
    elif problem == "membership":
        frame.loc[0, "module_membership_kME"] = np.inf
    elif problem == "module":
        frame.loc[0, "module"] = "M02"
    elif problem == "missing_event":
        frame = frame.iloc[1:]
    else:
        frame = pd.concat([frame, frame.iloc[[0]]])
    frame.to_csv(path, sep="\t", index=False)
    with pytest.raises(ValueError):
        _read_observed_kme_hubs(path, members)
