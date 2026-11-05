"""Report regressions: real absence must differ from a hidden completed result."""
import json
from pathlib import Path
import re
import shutil
import subprocess

import pandas as pd
import pytest

from katana_next.reports.native_enrichment import native_enrichment_scope
from katana_next.reports.cosplice_presentation import refresh_cosplice_copy, refresh_cosplice_gene_labels
from katana_next.reports.gene_names import annotate_report_genes
from katana_next.reports.mapp_activity import activity_axis, refresh_activity_panel


def test_maei_labels_follow_declared_conditions_without_averaging_or_reordering(tmp_path):
    from hashlib import sha256
    project = {'samples': [{'sample_id': 'R4', 'condition': 'MUT_DEL_KA'},
                           {'sample_id': 'R2', 'condition': 'MUT_DEL_PBS'}]}
    columns = ['activities_R4', 'activities_R2']
    axis = activity_axis(columns, project)
    assert [x['label'] for x in axis] == ['MUT_DEL_KA\nR4', 'MUT_DEL_PBS\nR2']
    assert [x['column'] for x in axis] == columns
    assert activity_axis(['activities_unknown'], project)[0]['label'] == 'unknown'
    path = tmp_path/'native.tsv'
    frame = pd.DataFrame({'site':['3ss'], 'native_model_id':['model1'],
                          'activities_R4':[1.25], 'activities_R2':[-.75], 'pvalue':[.2]})
    frame.to_csv(path,sep='\t',index=False)
    before = path.read_bytes()
    panel = {'contrast_id':'KA_vs_PBS','metadata':{'component':'mapp_native',
             'native_plot_data_sha256':sha256(before).hexdigest()},
             'data_links':[{'role':'plot_data','path':str(path)}]}
    rendered, audit = refresh_activity_panel(panel, project, tmp_path/'report')
    assert Path(rendered['path']).is_file()
    assert audit['activity_axis'] == axis
    assert audit['new_fits'] == 0 and audit['native_values_unchanged']
    assert path.read_bytes() == before


def test_completed_ora_selection_requires_measured_background_and_current_rs():
    a = {'analysis_id': 'prioritize.enrichment', 'contrast_ids': ['B_vs_A'],
         'metadata': {'foreground_source': 'das', 'background_source': 'das',
                      'native_background_provenance': {'artifacts': ['recorded']},
                      'report_completion_evidence': 'completed_workflow_manifest',
                      'das_selection': {'event_set': 'stable_significant',
                                        'read_support_policy': 'bulk-all-samples-counts-v1'}}}
    assert native_enrichment_scope(a)['policy'] == 'stable_significant'
    a['metadata']['das_selection']['read_support_policy'] = 'historical'
    assert native_enrichment_scope(a) is None
    a['metadata']['foreground_source'] = 'deg'
    assert native_enrichment_scope(a)['policy'] is None
    a['contrast_ids'] = ['B_vs_A', 'C_vs_A']
    assert native_enrichment_scope(a) is None


def test_report_gene_labels_preserve_statistics_and_unambiguous_reference():
    native = pd.DataFrame({'gene_id': ['ENSMUSG1', 'ENSMUSG2', 'ENSMUSG3'],
                           'gene_symbol': ['ENSMUSG1.16', 'ENSMUSG2.8', 'ENSMUSG3.1'],
                           'log2fc': [5.0996, 3.9781, 1.0], 'deg_fdr': [1.4888e-145, 1.7439e-143, .8]})
    annotation = pd.DataFrame({'gene_id': ['ENSMUSG1', 'ENSMUSG2', 'ENSMUSG2'],
                               'gene_symbol': ['Ythdc1', 'NameA', 'NameB']})
    resolved = annotate_report_genes(native, annotation)
    assert resolved.gene_symbol.tolist() == ['Ythdc1', 'ENSMUSG2', 'ENSMUSG3']
    assert resolved.source_gene_symbol.tolist() == native.gene_symbol.tolist()
    pd.testing.assert_frame_equal(resolved[['gene_id', 'log2fc', 'deg_fdr']], native[['gene_id', 'log2fc', 'deg_fdr']])
    assert resolved.gene_symbol_source.tolist() == ['project_gtf', 'unresolved', 'unresolved']


def test_cosplice_copy_preserves_native_payload_and_browses_nonsignificant_modules(tmp_path):
    root = Path(__file__).resolve().parents[2]
    template = (root/'katana_cosplice/_study_core/science/study_report39.html').read_text(encoding='utf-8')
    native = json.dumps({'modules': [{'module': 'M1'}, {'module': 'M2'}],
                         'associations': [{'module': 'M1', 'q_value': .1322}, {'module': 'M2', 'q_value': None}]})
    original = template.replace('__REPORT_DATA__', native).replace('__STUDY_TITLE__', 'Study')
    original = original.replace('<b>01</b>', '<b>6.1</b>').replace('CO-SPLICING ANALYSIS', 'Chapter 06')
    source = tmp_path/'source.html'; source.write_text(original, encoding='utf-8')
    target = tmp_path/'copy.html'; shutil.copyfile(source, target)
    receipt = refresh_cosplice_copy(target)
    html = target.read_text(encoding='utf-8')
    assert receipt['display_chapter'] == 7
    assert source.read_text(encoding='utf-8') == original
    assert re.search(r'id="reportData">(.*?)</script>', html, re.S)[1] == native
    assert '<b>7.1</b>' in html and 'Chapter 06' not in html
    assert 'const browse=scopeRows();' in html
    assert 'r=scopeRows().find(x=>x.module===selectedModule)' in html
    assert "sig=all.filter(r=>r.status==='estimated'&&num(r.q_value)&&r.q_value<=.05)" in html
    assert refresh_cosplice_copy(target) is None
    node = shutil.which('node')
    if node:
        script = re.findall(r'<script>(.*?)</script>', html, re.S)[-1]
        js = tmp_path/'report.js'; js.write_text(script, encoding='utf-8')
        subprocess.run([node, '--check', str(js)], check=True, capture_output=True)


def test_unknown_cosplice_template_is_not_silently_partially_rewritten(tmp_path):
    target = tmp_path/'report.html'
    text = '<script type="application/json" id="reportData">{"modules":[]}</script>'
    target.write_text(text, encoding='utf-8')
    with pytest.raises(ValueError, match='requires review'):
        refresh_cosplice_copy(target)
    assert target.read_text(encoding='utf-8') == text


def test_zero_eligible_motif_tests_are_not_labelled_nonsignificant():
    from katana_next.reports.motif_sequences import summary_html
    result = summary_html(dict(regional_counts={'affected_sequence':0,'regulatory_pre_mRNA':0},
        regional_available=False,window_test_count=0,window_estimable_models=0,
        window_significant_models=0,displayed_models=0,eligible_models=0,selection_rule='Recorded native tests only'))
    assert 'not estimable for this comparison' in result
    assert 'completed, no significant' not in result


def test_browser_can_inspect_modules_when_no_association_is_significant(tmp_path):
    if not __import__('os').environ.get('KATANA_REPORT_BROWSER_TEST'):
        pytest.skip('Run in the dedicated Linux Chromium job')
    from playwright.sync_api import sync_playwright
    root = Path(__file__).resolve().parents[2]
    template = (root/'katana_cosplice/_study_core/science/study_report39.html').read_text(encoding='utf-8')
    modules = [{'module': m, 'events': 1, 'pc1_variance_fraction': .8,
                'latent_factor': 'unassigned', 'loading_direction': 'positive'} for m in ('M1', 'M2')]
    associations = [dict(module=m, kind='omnibus', factor='condition', status='estimated' if i == 0 else 'NA',
                         p_value=.0005 if i == 0 else None, q_value=.1322 if i == 0 else None,
                         qualified_candidate=False, qualification_reason='association_FDR_above_0_05',
                         reason='Insufficient coverage' if i else '', statistic=1)
                    for i, m in enumerate(('M1', 'M2'))]
    samples = [dict(sample_id=s, biological_sample_id=s, condition=c) for s,c in [('S1','A'),('S2','B'),('S3','C')]]
    members = [dict(module=m, event_key='E'+m, gene_id='G'+m, gene_symbol='Symbol'+m, kME=.8, observed_samples=3,
                    within_module_strength=.5, event_type='SE') for m in ('M1', 'M2')]
    data = dict(modules=modules, associations=associations, samples=samples, members=members,
                scores=[dict(sample_id=s['sample_id'], M1=i, M2=-i) for i,s in enumerate(samples)],
                genes=[dict(module=m, gene_id='G'+m, events=1, max_abs_kME=.8) for m in ('M1','M2')],
                graphs={m:dict(nodes=['E'+m], edges=[], total_nodes=1,total_edges=0) for m in ('M1','M2')},
                event_usage={'EM1':[.1,.2,.3],'EM2':[.4,.5,.6]},
                module_files={m:dict(genes=m+'-genes.tsv',events=m+'-events.tsv') for m in ('M1','M2')},
                manifest=dict(study_id='Synthetic',samples=3,biological_units=3,network_id='fixture',assignment_sha256='fixture'),
                result=dict(input_events=2,eligible_events=2,assigned_events=2,
                            pathway_enrichment=dict(background_genes=2,network_events=2,tests=2,tested_modules=2,
                              libraries=[dict(module=i,library='GO',module_genes=1,all_module_genes=1,background_genes=2,
                                              network_genes=2,status='tested') for i in (1,2)])),
                analysis={}, assembly={}, qualification=dict(status='not_measured'),pathways=[],downloads=[])
    target = tmp_path/'report.html'
    target.write_text(template.replace('__REPORT_DATA__',json.dumps(data)).replace('__STUDY_TITLE__','Synthetic'),encoding='utf-8')
    refresh_cosplice_copy(target)
    refresh_cosplice_gene_labels(target)
    with sync_playwright() as p:
        browser = p.chromium.launch()
        page = browser.new_page()
        errors=[]; page.on('pageerror', lambda e: errors.append(str(e)))
        page.goto(target.as_uri())
        assert page.locator('#moduleSelect option').count() == 2
        page.locator('#moduleSelect').select_option('M2')
        assert page.locator('#moduleBody').is_visible()
        assert 'M2' in page.locator('#moduleTitle').inner_text()
        assert page.locator('#networkPlot svg').count() == 1
        assert 'GM2' in page.locator('#independentModuleLists').inner_text()
        assert 'SymbolM2' in page.locator('#independentModuleLists').inner_text()
        assert 'SymbolM2' in page.locator('#networkPlot').inner_text()
        page.locator('#geneSearch').fill('SymbolM2')
        assert 'GM2' in page.locator('#geneTable').inner_text()
        page.locator('#geneTable tbody tr').first.click()
        assert 'SymbolM2' in page.locator('#eventPlot').inner_text()
        assert 'No significant terms' in page.locator('#pathwayTable').inner_text()
        assert 'No module passes' in page.locator('#topModules').inner_text()
        assert not errors, errors
        browser.close()


def test_cosplice_labels_keep_native_data_and_hardlinked_downloads_unchanged(tmp_path):
    import os
    import zipfile
    template = (Path(__file__).resolve().parents[2]/'katana_cosplice/_study_core/science/study_report39.html').read_text('utf-8')
    data = dict(modules=[{'module':'M1'}], genes=[{'gene_id':'G1','module':'M1'}],
        members=[{'gene_id':'G1','gene_symbol':'Alpha'}, {'gene_id':'G2','gene_symbol':'Beta'},
                 {'gene_id':'G2','gene_symbol':'Other'}], module_files={'M1':{'genes':'module_M1_genes.tsv'}})
    source = tmp_path/'original'; source.mkdir()
    target = tmp_path/'copy'; (target/'report_tables').mkdir(parents=True)
    native = json.dumps(data)
    (source/'report.html').write_text(template.replace('__REPORT_DATA__',native),encoding='utf-8')
    original_table = b'gene_id\tmodule\tmax_abs_kME\nG1\tM1\t0.987654321\nG2\tM1\t0.2\n'
    (source/'module_M1_genes.tsv').write_bytes(original_table)
    os.link(source/'report.html', target/'report.html')
    os.link(source/'module_M1_genes.tsv', target/'report_tables/module_M1_genes.tsv')
    with zipfile.ZipFile(target/'report_tables/all_module_lists.zip','w') as z:
        z.writestr('module_M1_genes.tsv',original_table)
    receipt = refresh_cosplice_gene_labels(target/'report.html')
    assert receipt['mapped_gene_ids'] == 1 and receipt['ambiguous_gene_ids'] == 1
    html = (target/'report.html').read_text('utf-8')
    assert re.search(r'id="reportData">(.*?)</script>',html,re.S)[1] == native
    updated = pd.read_csv(target/'report_tables/module_M1_genes.tsv',sep='\t',dtype=str)
    assert updated.gene_symbol.tolist() == ['Alpha','G2']
    assert updated.max_abs_kME.tolist() == ['0.987654321','0.2']
    assert (source/'module_M1_genes.tsv').read_bytes() == original_table
    assert 'katana-cosplice-gene-labels-v1' not in (source/'report.html').read_text('utf-8')
    with zipfile.ZipFile(target/'report_tables/all_module_lists.zip') as z:
        assert z.read('module_M1_genes.tsv') == (target/'report_tables/module_M1_genes.tsv').read_bytes()
    assert refresh_cosplice_gene_labels(target/'report.html') is None
