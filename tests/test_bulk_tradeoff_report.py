import pandas as pd
import pytest

from scripts.report_bulk_splicing_tradeoffs import paired_differences


def test_pair_identity_is_required():
    rows = []
    for i in range(3):
        for arm in ("CoSplice", "SUPPA_OPTICS_tutorial"):
            rows.append(dict(iteration=i, arm=arm, draw_sha256=str(i), ari_all_events=.5,
                ari_jointly_assigned=.8, coassignment_recall=.6, coassignment_precision=.7, candidate_coverage=.6))
    frame = pd.DataFrame(rows)
    result = paired_differences(frame)
    assert result.mean_difference.eq(0).all()
    pd.testing.assert_frame_equal(result, paired_differences(frame.sample(frac=1, random_state=17)))
    missing = frame.copy()
    missing["draw_sha256"] = None
    with pytest.raises(ValueError, match="identities"):
        paired_differences(missing)
    invalid = frame.copy()
    invalid.loc[0, "ari_all_events"] = float("inf")
    with pytest.raises(ValueError, match="Infinite"):
        paired_differences(invalid)
    frame.loc[1, "draw_sha256"] = "not-paired"
    with pytest.raises(ValueError, match="exactly paired"):
        paired_differences(frame)
