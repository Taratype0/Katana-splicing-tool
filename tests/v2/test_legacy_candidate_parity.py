import pandas as pd
import pytest

from katana_next.science.bulk import build_candidate_evidence, pairwise_candidate_membership, usage_gene_support, candidate_display_labels, candidate_pair_plot_counts


def test_downstream_usage_does_not_relabel_exon_significance_and_labels_are_unique():
    source = pd.DataFrame([
        dict(comparison_id="A", gene_id="G", source_kind="dtu_result", dtu_fdr=.9),
        dict(comparison_id="A", gene_id="G", source_kind="dexseq_result", dtu_fdr=.001),
        dict(comparison_id="A", gene_id="H", source_kind="dexseq_result", dtu_fdr=.01),
    ])
    support = usage_gene_support(source).set_index("gene_id")
    assert support.loc["G", "best_dtu_fdr"] == .9
    assert support.loc["G", "best_dexseq_fdr"] == .001
    assert pd.isna(support.loc["H", "best_dtu_fdr"])
    labels = pd.DataFrame([
        dict(comparison_id="A", gene_id="G", gene_symbol="same", row_id="A::G"),
        dict(comparison_id="A", gene_id="H", gene_symbol="same", row_id="A::H"),
        dict(comparison_id="B", gene_id="G", gene_symbol="same", row_id="B::G"),
    ])
    assert candidate_display_labels(labels, {"A": "same contrast", "B": "same contrast"}).is_unique


def test_dexseq_is_not_misreported_as_dtu_and_single_evidence_is_retained():
    usage = pd.DataFrame([
        {"comparison_id": "A", "gene_id": "exon-only", "transcript_id": "E1", "dtu_fdr": .01, "source_kind": "dexseq_result"},
        {"comparison_id": "A", "gene_id": "tx-only", "transcript_id": "T1", "dtu_fdr": .01, "source_kind": "dtu_result"},
        {"comparison_id": "A", "gene_id": "both", "transcript_id": "E1", "dtu_fdr": .01, "source_kind": "dexseq_result"},
        {"comparison_id": "A", "gene_id": "both", "transcript_id": "T1", "dtu_fdr": .01, "source_kind": "dtu_result"},
    ])
    table = build_candidate_evidence(pd.DataFrame(), pd.DataFrame(), usage, {}).set_index("gene_id")
    assert table.loc["exon-only", "dexseq_supported"]
    assert not table.loc["exon-only", "dtu_available"]
    assert table.loc["tx-only", "dtu_supported"]
    assert not table.loc["tx-only", "dexseq_available"]
    assert table.loc["both", "ranking_score"] == 2
    assert list(table.columns)[-1] == "ranking_score"


def test_pair_membership_and_gene_not_event_direction():
    events = pd.DataFrame([
        {"comparison_id": c, "gene_id": "G", "gene_symbol": "Gene", "event_id": e,
         "significant": True, "stable_significant": True, "dpsi": d, "abs_dpsi": abs(d), "fdr": .01}
        for c, e, d in [("A", "different-event-1", .2), ("B", "different-event-2", -.3)]
    ])
    deg = pd.DataFrame([{"comparison_id": "A", "gene_id": "A-only", "log2fc": 2., "deg_fdr": .01},
                        {"comparison_id": "B", "gene_id": "B-only", "log2fc": -2., "deg_fdr": .01}])
    table = build_candidate_evidence(events, deg, pd.DataFrame(), {})
    pairs = pairwise_candidate_membership(table, ["A", "B"]).set_index("gene_id")
    assert pairs.loc["G", "n_comparisons_candidate"] == 2
    assert pairs.loc["G", "direction_pattern"] == "gene_level_direction_changed"
    assert pairs.loc["G", "dominant_as_event_id_A"] != pairs.loc["G", "dominant_as_event_id_B"]
    assert pairs.loc["A-only", "membership"] == "A_specific"
    assert pairs.loc["B-only", "membership"] == "B_specific"
    assert pairs.loc["A-only", "direction_pattern"] == "direction_not_comparable"
    assert pairwise_candidate_membership(table, ["A"]).empty
    with pytest.raises(ValueError, match="one row"):
        pairwise_candidate_membership(pd.concat([table, table]), ["A", "B"])


def test_four_distinct_supported_analyses_and_all_pairs():
    events = pd.DataFrame([dict(comparison_id=c, gene_id="G", gene_symbol="Gene", event_id="e",
                               significant=True, stable_significant=True, abs_dpsi=.2, dpsi=.2, fdr=.01) for c in "ABC"])
    deg = pd.DataFrame([dict(comparison_id=c, gene_id="G", log2fc=2., deg_fdr=.01) for c in "ABC"])
    usage = pd.DataFrame([dict(comparison_id=c, gene_id="G", transcript_id=k, source_kind=k, dtu_fdr=.01)
                          for c in "ABC" for k in ["dtu_result", "dexseq_result"]])
    table = build_candidate_evidence(events, deg, usage, {})
    assert table.ranking_score.eq(4).all()
    pair = pairwise_candidate_membership(table, list("ABC"))
    assert len(pair) == 3 and pair.membership.eq("shared").all()
    assert pair.direction_pattern.eq("comparable_directions_agree").all()
    plotted = candidate_pair_plot_counts(pair, {c: "same label" for c in "ABC"})
    assert len(plotted) == 3 and plotted.pair.nunique() == 3
    assert plotted.genes.eq(1).all()
