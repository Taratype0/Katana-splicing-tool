import hashlib

import numpy as np
import pandas as pd
import pytest

from katana_cosplice.motif_scanner import Motif
from katana_next.science.se_rbp_density import (
    average_odds_scores, canonical_rbp, consensus_coverage, consensus_pwm,
    fixed_pattern, paired_sign_rows,
)


def catalog(matrices):
    return pd.DataFrame([dict(RBP="NOVA1", evidence_class="direct_experiment",
        motif_token=f"p{i}", pwm_sha256=hashlib.sha256(m.astype("<f8").tobytes()).hexdigest())
        for i, m in enumerate(matrices)])


@pytest.mark.parametrize("text", ["A.*C", "(U)(n)", "R(U)(1-3)AGU", "[AZ]", ""])
def test_reject_unbounded_or_invalid_patterns(text):
    with pytest.raises(ValueError):
        fixed_pattern(text)


def test_rna_and_native_alias():
    assert fixed_pattern("YCAY") == ("[CT]CA[CT]", 4)
    assert fixed_pattern("UGCU")[0] == "TGCT"
    assert canonical_rbp("SRp20") == "SRSF3"
    assert np.allclose(consensus_pwm("YCAY").sum(axis=1), 1)


def test_overlapping_consensus_union_and_junction_filter():
    patterns = pd.DataFrame([dict(pattern="AAA", width=3)])
    windows = pd.DataFrame(dict(sequence=["AAAAA", "AAACCC", "AA"],
                                junction_offset=[np.nan, 3., np.nan]))
    result = consensus_coverage(patterns, windows)[0]
    assert result[0] == 1
    assert result[1] == 0  # hit ending at boundary is not crossing
    assert np.isnan(result[2])


def test_vectorized_average_odds_matches_full_naive_scan():
    matrices = [consensus_pwm("YCAY"), np.full((4, 4), .25)]
    motifs = [Motif("NOVA1", str(i), "fixture", pwm=m) for i, m in enumerate(matrices)]
    windows = pd.DataFrame(dict(sequence=["TCATACCA", "TCAANTCAT", "TCATACCA", "AA"],
                                junction_offset=[np.nan, np.nan, 4., np.nan]))
    bg = np.full(4, .25)
    actual = average_odds_scores(motifs, catalog(matrices), windows, bg)
    for i, matrix in enumerate(matrices):
        smooth = (matrix*20+.1*bg)/20.1
        for j, row in enumerate(windows.itertuples(index=False)):
            odds = []
            for start in range(len(row.sequence)-3):
                word = row.sequence[start:start+4]
                if "N" in word or (np.isfinite(row.junction_offset) and not start < row.junction_offset < start+4):
                    continue
                odds.append(np.prod([smooth[k, "ACGT".index(x)]/.25 for k, x in enumerate(word)]))
            expected = np.log2(1+np.mean(odds)) if odds else np.nan
            assert actual[i, j] == pytest.approx(expected, nan_ok=True)
    assert actual[0, 0] > 0  # short motif contributes without impossible site P filter
    assert np.allclose(actual[1, :3], 1.)


def test_fixed_paired_sign_test_retains_all_hypotheses_and_ties():
    windows = pd.DataFrame(dict(event_id=["f1", "b1", "f2", "b2", "f3", "b3"],
                                region=["cassette_exon"]*6))
    pairs = pd.DataFrame(dict(foreground_event=["f1", "f2", "f3"],
        background_event=["b1", "b2", "b3"], direction=["inclusion_down"]*3, status=["matched"]*3))
    result = paired_sign_rows(np.array([[2, 1, 3, 1, 1, 1]]),
        catalog([np.full((4, 4), .25)]), windows, pairs, "fixed", "PWM_average_odds")
    row = result[result.region.eq("cassette_exon") & result.direction.eq("inclusion_down")].iloc[0]
    assert row.pvalue == .25
    assert row.ties == 1
    assert row.matched_pairs == 3
    assert len(result) == 12
