"""English-only rendering and preservation of cached scientific results."""
import hashlib
import json
import re
import zipfile
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest
from bs4 import BeautifulSoup

from katana_next.science.narrative_report import build_report, package_report
from katana_next.science.report_complete import complete_sections
from katana_next.science.report_language import copy_english_reader, is_retired_report
from katana_next.workflows.narrative_report import NarrativeReportWorkflow


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_table(root, relative, data):
    path = root / relative
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(data).to_csv(path, sep='\t', index=False)
    return path


def test_english_reader_does_not_require_chinese_and_preserves_original(tmp_path):
    source = tmp_path / 'prepared'
    (source / 'tables').mkdir(parents=True)
    (source / 'figures').mkdir()
    (source / 'figures/test.png').write_bytes(b'figure')
    (source / 'tables/all.tsv').write_text('P\tFDR\n0.04\t0.08\n', encoding='utf-8')
    english = source / 'REPORT_EN.html'
    english.write_text('<html lang="en"><section id="results"><img src="figures/test.png"><a href="tables/all.tsv">All results</a></section></html>', encoding='utf-8')
    manifest = tmp_path / 'manifest.json'
    manifest.write_text('{}')
    config_path = tmp_path / 'config.json'
    config_path.write_text(json.dumps({'study': 'Fixture', 'approved_reader_root': str(source), 'manifest_path': str(manifest)}))
    workflow = NarrativeReportWorkflow()
    assert workflow._config({'config_path': str(config_path)})['study'] == 'Fixture'
    assert manifest in workflow.input_paths(None, {'config_path': str(config_path)})
    chinese = source / 'REPORT_ZH.html'
    chinese.write_text('<p>Historical Chinese report</p>', encoding='utf-8')
    english.write_text(english.read_text().replace('</section>', '<a href="REPORT_ZH.html">Previous language</a></section>'), encoding='utf-8')
    before = {p: digest(p) for p in source.rglob('*') if p.is_file()}
    rendered = copy_english_reader(source, tmp_path / 'reader')
    assert not (rendered.parent / 'REPORT_ZH.html').exists()
    soup = BeautifulSoup(rendered.read_text(encoding='utf-8'), 'html.parser')
    assert all((rendered.parent / a['href']).is_file() for a in soup.find_all('a'))
    assert before == {p: digest(p) for p in before}
    with pytest.raises(FileExistsError):
        copy_english_reader(source, rendered.parent)


def test_language_guard_and_filename_policy(tmp_path):
    with pytest.raises(ValueError, match='Only English'):
        complete_sections(tmp_path, tmp_path, 'ZH')
    assert is_retired_report('REPORT_ZH.html')
    assert is_retired_report('study_REPORT_ZH.html')
    assert is_retired_report('study_Bilingual_complete_results.zip')
    assert is_retired_report('REAL_CASE_CONCLUSION_ZH.md')
    assert not is_retired_report('ZH_gene_counts.tsv')
    assert not is_retired_report('REPORT_EN.html')


def test_report_generation_and_zip_need_no_chinese_first_pass(tmp_path):
    base = tmp_path / 'result'
    write_table(base, '20_Presentation_report/comparison_totals.tsv', {
        'comparison_id': ['KO_vs_WT'], 'significant_events': [1], 'DEG': [1]})
    write_table(base, '20_Presentation_report/DEG_DAS_gene_intersection.tsv', {
        'comparison_id': ['KO_vs_WT'], 'gene_id': ['ENSMUSG00000000001'], 'P': [.031], 'FDR': [.062]})
    write_table(base, '22_Expression_interaction_pathways/design/run/interaction.tsv', {
        'gene_id': ['ENSMUSG00000000001'], 'gene_symbol': ['Example'], 'padj': [.062]})
    write_table(base, '22_Expression_interaction_pathways/design/run/ranked_pathways.tsv', {
        'contrast': ['interaction'], 'pathway': ['KEGG_EXAMPLE'], 'study_family_FDR': [.2], 'NES': [-1.2], 'leadingEdge': ['Example']})
    write_table(base, '16_Interaction/design/run/dexseq_interaction.tsv', {
        'groupID': ['ENSMUSG00000000001'], 'bin_id': ['E001'], 'interaction_FDR': [.08], 'pvalue': [.04]})
    write_table(base, '23_Module_explanation/run/module_cards.tsv', {'module': ['C1'], 'n_events': [2]})
    chart = base / '04_DAS/run/differential_as_counts.vl.json'
    chart.parent.mkdir(parents=True)
    chart.write_text(json.dumps({'data': {'values': [{
        'comparison_id': 'KO_vs_WT', 'event_type': 'SE',
        'tested_events': 7, 'significant_events': 1, 'stable_significant_events': 1}]}}))
    manifest = tmp_path / 'project.json'
    manifest.write_text(json.dumps({'samples': [{'sample_id': 's1', 'condition': 'WT', 'biological_sample_id': 'mouse1'}]}))
    before = {p: digest(p) for p in base.rglob('*') if p.is_file()}
    output = base / '24_Narrative_report'
    build_report('Fixture', base, manifest, output)
    english = output / 'REPORT_EN.html'
    text = english.read_text(encoding='utf-8')
    assert not re.search('[\u4e00-\u9fff]', text)
    assert 'REPORT_ZH.html' not in text
    assert not (output / 'REPORT_ZH.html').exists()
    assert (output / 'figures/AS_type_distribution.png').is_file()
    frame = pd.read_csv(output / 'figures/AS_type_distribution_source.tsv', sep='\t')
    assert frame.tested_events.tolist() == [7]
    assert frame.significant_events.tolist() == [1]
    assert before == {p: digest(p) for p in before}
    # A historical report may still exist in source results; do not erase it or
    # let the new English download package silently include it.
    stale = output / 'REPORT_ZH.html'
    stale.write_text('<p>Historical Chinese report</p>', encoding='utf-8')
    archive = tmp_path / 'English_results.zip'
    package_report('Fixture', base, output, archive)
    with zipfile.ZipFile(archive) as z:
        assert z.testzip() is None
        assert any(name.endswith('/REPORT_EN.html') for name in z.namelist())
        assert not any('REPORT_ZH' in name for name in z.namelist())
        assert any(name.endswith('/DEG_DAS_gene_intersection.tsv') for name in z.namelist())
    assert stale.exists()


def test_english_filter_footnote_keeps_counts_units_and_statistical_boundary(monkeypatch):
    from katana_next.science import gene_filter
    details = dict(active=True, excluded_annotation_genes=23, annotation_genes=100,
                   config=dict(exclude_mitochondrial=True, exclude_pseudogenes=True,
                               exclude_rrna=False, protein_coding_only=False,
                               excluded_genes=['example'], min_expression=1, min_samples=3))
    monkeypatch.setattr(gene_filter, 'gene_filter_details', lambda manifest: details)
    text = gene_filter.gene_filter_footnote(None)
    assert 'mitochondrial-encoded genes' in text and 'pseudogenes' in text
    assert '23 / 100 annotated genes' in text
    assert 'not tested-gene counts' in text and 'models were not refitted' in text
    assert 'supplied matrix units' in text
    assert not re.search('[\u4e00-\u9fff]', text)
