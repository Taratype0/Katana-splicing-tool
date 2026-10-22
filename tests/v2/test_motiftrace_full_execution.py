import numpy as np
import pandas as pd
import pytest

from katana_next.science.motiftrace_clusters import _profile, _signed_ranks
from katana_next.science.motiftrace_full_execution import profiles, signed_ranks, fast_match


@pytest.mark.parametrize("pattern", ["AAA", "[CT]CT[CT]", "AC", "A", "[ACGT]A"])
def test_vectorized_profiles_preserve_full_windows_short_NA_and_greedy_counts(pattern):
    rng = np.random.default_rng(13)
    for word in ("", "AA", "AAAAAAA", "ACNNAC", "".join(rng.choice(list("ACGTN"), 95))):
        model = {"kind": "consensus", "pattern": pattern}
        actual = profiles(word, pattern)
        for width in (0, 10, 20, 50):
            for method, arm in (("finite_nonoverlap_counts", 0), ("overlapping_union_coverage", 1)):
                np.testing.assert_allclose(actual[width][arm], _profile(word, model, method, width), rtol=0, atol=1e-15, equal_nan=True)


def test_vectorized_signed_ranks_preserve_zero_ties_NA_and_denominators():
    rng = np.random.default_rng(22)
    matrix = rng.choice([0., 1., -1., np.nan, .2, -.2], size=(28, 91))
    for first, second in zip(signed_ranks(matrix), _signed_ranks(matrix, 20)):
        np.testing.assert_allclose(first, second, atol=1e-15, rtol=0, equal_nan=True)


def test_vectorized_null_search_preserves_NA_gaps_and_cluster_maxima():
    from katana_next.science.motiftrace_full_execution import null_cluster_max
    from katana_next.science.motiftrace_clusters import CLUSTER_Z, cluster_runs
    rng = np.random.default_rng(8)
    z = rng.uniform(0, 5, (99, 201))
    z[:, 50:54] = np.nan
    actual = null_cluster_max(z, CLUSTER_Z)
    expected = [max([m for _, _, m in cluster_runs(row)] or [0]) for row in z]
    np.testing.assert_allclose(actual, expected, atol=1e-12, rtol=0)


def test_vectorized_matching_preserves_original_greedy_gene_choice(tmp_path):
    from test_motiftrace_consequences import event, reference
    from katana_next.science.motiftrace_events import extract_event_rna, match_event_controls
    fasta, _ = reference(tmp_path)
    rows = [event(event_id=f"f{i}", gene_id=f"fg{i//2}", IncLevel2=f"{.3+i*.003}") for i in range(8)]
    rows += [event(event_id=f"b{i}", gene_id=f"bg{i//2}", stable_significant=False, fdr=.9,
        dpsi=.001, IncLevel2=f"{.3+i*.002}") for i in range(30)]
    frame = pd.DataFrame(rows)
    selected = {r["event_id"] for r in rows[:8]}
    rna = extract_event_rna(frame, fasta, feature_only=True)
    pd.testing.assert_frame_equal(fast_match(frame, selected, rna), match_event_controls(frame, selected, rna))
