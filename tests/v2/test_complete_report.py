import json
from katana_next.science.report_complete import complete_sections


def test_gene_cards_and_cross_comparison_are_in_body(tmp_path):
    out=tmp_path/'24_Narrative_report';out.mkdir()
    for folder,name in [('12_Gene_cards_heatmaps','gene_card_example.png'),('07_Cross_comparison','selected_contrast_intersections.png')]:
        root=tmp_path/folder/'run'/'plots';root.mkdir(parents=True)
        (root/name).write_bytes(b'fixture')
    body=complete_sections(tmp_path,out,'EN')
    assert 'src="../12_Gene_cards_heatmaps/run/plots/gene_card_example.png"' in body
    assert 'src="../07_Cross_comparison/run/plots/selected_contrast_intersections.png"' in body
    assert len(json.loads((out/'MODULE_COVERAGE.json').read_text()))==15


def test_missing_figures_not_reported_as_zero_results(tmp_path):
    out=tmp_path/'report';out.mkdir()
    body=complete_sections(tmp_path,out,'EN')
    assert 'This is not a zero result.' in body
    assert 'first 8 rows' not in body


def test_figure_link_quotes_are_escaped_without_losing_the_figure(tmp_path):
    out=tmp_path/'report';out.mkdir()
    root=tmp_path/'07_Cross_comparison'/'run'/'plots';root.mkdir(parents=True)
    (root/'selected&contrast.png').write_bytes(b'fixture')
    body=complete_sections(tmp_path,out,'EN')
    assert 'selected&amp;contrast.png' in body
    assert '<figure>' in body
