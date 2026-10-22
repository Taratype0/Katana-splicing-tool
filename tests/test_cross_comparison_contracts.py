from app.pages.cross_comparison_controller import CrossComparisonController


def test_primary_cross_comparison_contracts_keep_gene_and_event_levels_separate():
    assert CrossComparisonController.contract("matrix").result_level == "gene"
    assert CrossComparisonController.contract("shared").result_level == "gene"
    for subbranch in ("interaction", "venn", "venn_multi", "venn_multi_similarity"):
        assert CrossComparisonController.contract(subbranch).result_level == "event"


def test_cross_comparison_pair_tab_requirements_follow_page_semantics():
    assert CrossComparisonController.needs_pair_tabs("matrix")
    assert CrossComparisonController.needs_pair_tabs("interaction")
    assert not CrossComparisonController.needs_pair_tabs("interaction", interaction_view="batch_summary")
    assert not CrossComparisonController.needs_pair_tabs("venn_multi")
    assert not CrossComparisonController.needs_pair_tabs("venn_multi_similarity")


def test_cross_comparison_actions_are_independent_and_stable():
    assert CrossComparisonController.action("matrix") == "candidate_matrix"
    assert CrossComparisonController.action("shared") == "candidate_matrix"
    assert CrossComparisonController.action("interaction") == "interaction"
    assert CrossComparisonController.action("venn") == "venn"
    assert CrossComparisonController.action("venn_multi") == "venn_multi"
    assert CrossComparisonController.action("venn_multi_similarity") == "venn_multi"
    assert CrossComparisonController.run_label("venn_multi_similarity") == "Run 5.4.6 Similarity Heatmap"
