import base64
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import shutil

import pytest

from katana_next.reports.compact import compact_report, reconcile_compact_metadata, recover_compact_report


def _report(tmp_path, *, rows=None):
    root = tmp_path / 'report'
    (root / 'outputs' / 'results').mkdir(parents=True)
    (root / 'outputs' / 'comparisons' / 'B_vs_A').mkdir(parents=True)
    (root / 'outputs' / 'views').mkdir(parents=True)
    source = root / 'outputs' / 'results' / 'native.tsv'
    source.write_bytes(b'gene\tp\tfdr\teffect\nG1\t1e-200\t0.01\t-Infinity\n')
    alias = root / 'outputs' / 'comparisons' / 'B_vs_A' / 'same.tsv'
    alias.write_bytes(source.read_bytes())
    sha = hashlib.sha256(source.read_bytes()).hexdigest()
    (root / 'REPORT.html').write_text('<html><a href="outputs/comparisons/B_vs_A/same.tsv?download=1#full">raw</a></html>', encoding='utf-8')
    (root / 'report-manifest.json').write_text(json.dumps({'complete_results': [{'source': 'outputs/comparisons/B_vs_A/same.tsv'}], 'panels': [{'data_links': [{'path':'outputs/comparisons/B_vs_A/same.tsv','source_sha256':sha,'portable_sha256':sha,'sha256':sha}]}]}), encoding='utf-8')
    (root / 'outputs' / 'views' / 'native.json').write_text('{"schema":"katana.portable_native_MAPP_fit_index/v1","fits":[{"path":"../comparisons/B_vs_A/same.tsv","p":1e-200,"source_sha256":"'+sha+'"}]}',encoding='utf-8')
    template = Path(__file__).resolve().parents[2] / 'katana_next' / 'reports' / 'assets' / 'table-template.html'
    data = json.dumps({'columns':['gene','p','effect','missing','source_path'], 'rows':rows or [['G1','1e-200','-Infinity',None,'../comparisons/B_vs_A/same.tsv'],['G2','0','Infinity','', 'x']]}, separators=(',', ':'))
    view = root / 'outputs' / 'views' / 'table.html'
    view.write_text(template.read_text(encoding='utf-8').replace('__TITLE__','Native result').replace('__RAW__','../comparisons/B_vs_A/same.tsv').replace('__DATA__',data),encoding='utf-8')
    return root, source, alias, view, data.encode('utf-8'), sha


def test_dedup_rewrites_portable_links_preserves_native_scalar_tokens(tmp_path):
    root, source, alias, view, raw, sha = _report(tmp_path,rows=[['G1','1e-200','-Infinity',None,'../comparisons/B_vs_A/same.tsv']]*100)
    original = source.read_bytes()
    audit = compact_report(root, html_min_bytes=0)
    assert source.read_bytes() == original
    assert not alias.exists()
    assert 'outputs/results/native.tsv?download=1#full' in (root/'REPORT.html').read_text(encoding='utf-8')
    native_index=(root/'outputs/views/native.json').read_text(encoding='utf-8')
    assert '"path":"../results/native.tsv"' in native_index
    assert '"p":1e-200' in native_index
    manifest=json.loads((root/'report-manifest.json').read_text(encoding='utf-8'))
    binding=manifest['panels'][0]['data_links'][0]
    assert binding['source_sha256'] == binding['sha256'] == binding['portable_sha256'] == sha
    assert binding['path'] == 'outputs/results/native.tsv'
    assert manifest['complete_results'][0]['source'] == 'outputs/results/native.tsv'
    content=view.read_text(encoding='utf-8')
    assert 'href="../results/native.tsv"' in content
    assert '(async()=>{try{' in content and 'DecompressionStream' in content
    compressed=re.search(r'<script id="table-data"[^>]*>(.*?)</script>',content,re.S).group(1)
    assert gzip.decompress(base64.b64decode(compressed)) == raw
    assert audit['native_statistics_recomputed'] is False
    assert audit['complete_table_payloads'][0]['rows'] == 100
    assert len(audit['deleted_duplicate_aliases']) == 1


def test_atomic_compaction_never_changes_hardlinked_source_html(tmp_path):
    root, source, alias, view, raw, sha = _report(tmp_path)
    original = tmp_path/'original-table.html'
    os.link(view, original)
    old=original.read_bytes()
    compact_report(root, html_min_bytes=0)
    assert original.read_bytes() == old
    assert view.read_bytes() != old


def test_rejects_previously_archived_report_before_mutation(tmp_path):
    root, source, alias, view, raw, sha = _report(tmp_path)
    (root/'REPORT.zip').write_bytes(b'existing delivered archive')
    old=view.read_bytes()
    with pytest.raises(ValueError, match='unarchived'):
        compact_report(root)
    assert alias.exists() and view.read_bytes()==old


def test_supplied_aliases_must_be_identical_and_inside_root(tmp_path):
    root, source, alias, view, raw, sha = _report(tmp_path)
    with pytest.raises(ValueError, match='escapes'):
        compact_report(root, alias_map={'outputs/results/native.tsv':'../outside.tsv'})
    alias.write_text('different native bytes',encoding='utf-8')
    with pytest.raises(ValueError, match='differ'):
        compact_report(root, alias_map={'outputs/comparisons/B_vs_A/same.tsv':'outputs/results/native.tsv'})


def test_position_dependent_html_not_deduplicated(tmp_path):
    root, source, alias, view, raw, sha = _report(tmp_path)
    duplicate = root/'outputs'/'views'/'same.html'
    duplicate.write_bytes(view.read_bytes())
    compact_report(root, html_min_bytes=0)
    assert view.exists() and duplicate.exists()


def test_preserves_every_table_row_and_missing_infinity_strings(tmp_path):
    rows=[['G'+str(i),'1e-200','-Infinity' if i%2 else 'Infinity',None,'../comparisons/B_vs_A/same.tsv'] for i in range(501)]
    root, source, alias, view, raw, sha = _report(tmp_path,rows=rows)
    compact_report(root, html_min_bytes=0)
    compressed=re.search(r'<script id="table-data"[^>]*>(.*?)</script>',view.read_text(encoding='utf-8'),re.S).group(1)
    actual=json.loads(gzip.decompress(base64.b64decode(compressed)))
    assert actual['rows']==rows


def test_provenance_paths_and_json_file_index_keys_are_native_identities(tmp_path):
    root, source, alias, view, raw, sha = _report(tmp_path)
    (root/'outputs/provenance/native').mkdir(parents=True)
    provenance=root/'outputs/provenance/native/provenance.json'
    original='{"outputs/comparisons/B_vs_A/same.tsv":"'+sha+'","source_path":"outputs/comparisons/B_vs_A/same.tsv"}'
    provenance.write_text(original,encoding='utf-8')
    index=root/'outputs/views/file-index.json'
    index.write_text('{"outputs/comparisons/B_vs_A/same.tsv":"native identity"}',encoding='utf-8')
    compact_report(root,html_min_bytes=0)
    assert provenance.read_text(encoding='utf-8')==original
    assert '"outputs/comparisons/B_vs_A/same.tsv":' in index.read_text(encoding='utf-8')


def test_generated_report_provenance_portable_bindings_are_rebound(tmp_path):
    root, source, alias, view, raw, sha = _report(tmp_path)
    provenance=root/'provenance.json'
    provenance.write_text(json.dumps({'schema':'katana.report_provenance/v1',
        'source_files':[{'source_name':'same.tsv','source_sha256':sha,'bytes':source.stat().st_size,'report_path':'outputs/comparisons/B_vs_A/same.tsv'}],
        'figure_data_links':{'bindings':[{'path':'outputs/comparisons/B_vs_A/same.tsv','source_sha256':sha,'portable_sha256':sha}]}}),encoding='utf-8')
    compact_report(root,html_min_bytes=0)
    actual=json.loads(provenance.read_text(encoding='utf-8'))
    binding=actual['figure_data_links']['bindings'][0]
    assert binding['path']=='outputs/results/native.tsv'
    assert binding['source_sha256']==binding['portable_sha256']==sha
    record=actual['source_files'][0]
    assert record['report_path']=='outputs/results/native.tsv'
    assert record['source_name']=='same.tsv' and record['source_sha256']==sha
    assert record['source_bytes']==record['bytes']==source.stat().st_size


def test_equal_bytes_with_different_decoder_extensions_are_not_deduplicated(tmp_path):
    root, source, alias, view, raw, sha = _report(tmp_path)
    csv=root/'outputs/results/native.csv'
    csv.write_bytes(source.read_bytes())
    compact_report(root)
    assert csv.exists() and source.exists()


def test_scientific_json_path_cells_remain_native_bytes(tmp_path):
    root, source, alias, view, raw, sha = _report(tmp_path)
    result=root/'outputs/results/native-measurements.json'
    original=b'{"rows":[{"gene":"G1","path":"outputs/comparisons/B_vs_A/same.tsv","P":1e-200,"fdr":0.01}]}'
    result.write_bytes(original)
    compact_report(root,html_min_bytes=0)
    assert result.read_bytes()==original


def test_metadata_reconciliation_after_aliases_removed_and_old_native_json_restore(tmp_path):
    root, source, alias, view, raw, sha = _report(tmp_path)
    original_root=tmp_path/'immutable-source'
    (original_root/'outputs/results').mkdir(parents=True)
    native_relative='outputs/results/native-measurements.json'
    native_bytes=b'{"path":"outputs/comparisons/B_vs_A/same.tsv","p":1e-200}'
    (original_root/native_relative).write_bytes(native_bytes)
    (root/native_relative).write_bytes(native_bytes)
    audit=compact_report(root,html_min_bytes=0)
    assert not alias.exists()
    # Simulate an older pass that left generated root provenance stale and
    # rewrote a native JSON cell; the recorded before-SHA enables exact repair.
    (root/'provenance.json').write_text(json.dumps({'source_files':[{'source_name':'same.tsv','source_sha256':sha,'bytes':source.stat().st_size,'report_path':'outputs/comparisons/B_vs_A/same.tsv'}],
        'figure_data_links':{'bindings':[{'path':'outputs/comparisons/B_vs_A/same.tsv','source_sha256':sha,'portable_sha256':sha}]}}),encoding='utf-8')
    (root/native_relative).write_bytes(b'{"path":"outputs/results/native.tsv","p":1e-200}')
    audit['rewritten_files'].append({'path':native_relative,'previous_sha256':hashlib.sha256(native_bytes).hexdigest()})
    receipt=reconcile_compact_metadata(root,audit,source_report_dir=original_root)
    assert (root/native_relative).read_bytes()==native_bytes
    binding=json.loads((root/'provenance.json').read_text(encoding='utf-8'))['figure_data_links']['bindings'][0]
    assert binding['path']=='outputs/results/native.tsv'
    assert binding['source_sha256']==binding['portable_sha256']==sha
    assert receipt['restored_native_json'][0]['restored_source_sha256']==hashlib.sha256(native_bytes).hexdigest()


def test_all_embedded_native_json_and_octet_stream_payloads_are_immutable(tmp_path):
    root, source, alias, view, raw, sha = _report(tmp_path)
    payload='{"members":[{"gene":"G1","path":"../comparisons/B_vs_A/same.tsv","P":1e-200,"fdr":0.01}],"identity":"gene\\u0027s original value"}'
    opaque='["../comparisons/B_vs_A/same.tsv","Infinity",null]'
    content=view.read_text(encoding='utf-8')
    content += '<script id="data" type="application/json">'+payload+'</script>'
    content += "<script id='module-display-data' type='APPLICATION/OCTET-STREAM'>"+opaque+'</script>'
    content += '<script>const navigation="../comparisons/B_vs_A/same.tsv";</script>'
    view.write_text(content,encoding='utf-8')
    compact_report(root,html_min_bytes=0)
    actual=view.read_text(encoding='utf-8')
    assert '<script id="data" type="application/json">'+payload+'</script>' in actual
    assert "<script id='module-display-data' type='APPLICATION/OCTET-STREAM'>"+opaque+'</script>' in actual
    assert 'const navigation="../results/native.tsv"' in actual


def test_recovery_proves_prior_encoded_tables_and_handles_already_removed_alias(tmp_path):
    source, raw_source, alias, view, raw, sha = _report(tmp_path/'immutable',rows=[['G1','1e-200','-Infinity',None,'../comparisons/B_vs_A/same.tsv']]*501)
    new=tmp_path/'partial-report'
    shutil.copytree(source,new)
    compact_report(new,html_min_bytes=0)
    old_html=(source/'REPORT.html').read_bytes()
    events=[]
    audit=recover_compact_report(new,source,{'outputs/comparisons/B_vs_A/same.tsv':'outputs/results/native.tsv'},progress=events.append)
    assert (source/'REPORT.html').read_bytes()==old_html
    assert (source/'outputs/comparisons/B_vs_A/same.tsv').exists()
    assert audit['recovery']['already_encoded_tables']==1
    table=audit['complete_table_payloads'][0]
    assert table['rows']==501 and table['immutable_source_payload_verified']
    assert table['json_sha256']==hashlib.sha256(raw).hexdigest()
    assert audit['recovery']['scientific_results_rerun']==0
    assert events[-1]['phase']=='completed'


def test_custom_mapp_and_rich_browser_portable_hashes_refresh(tmp_path):
    root, source, alias, view, raw, sha = _report(tmp_path,rows=[['G1','1e-200','-Infinity',None,'../comparisons/B_vs_A/same.tsv']]*501)
    native_sha='a'*64
    catalog=root/'outputs/figure-data/abc-complete_native_fits_portable.json'
    catalog.parent.mkdir()
    catalog.write_text('{"schema":"katana.portable_native_MAPP_fit_index/v1","fits":[{"path":"../comparisons/B_vs_A/same.tsv","sha256":"'+sha+'"}]}',encoding='utf-8')
    catalog_before=hashlib.sha256(catalog.read_bytes()).hexdigest()
    browser_before=hashlib.sha256(view.read_bytes()).hexdigest()
    manifest=json.loads((root/'report-manifest.json').read_text(encoding='utf-8'))
    manifest['panels'].append({'metadata':{'portable_fit_index_sha256':catalog_before,'original_native_fit_index_sha256':native_sha},'data_links':[{'path':'outputs/figure-data/abc-complete_native_fits_portable.json','source_sha256':catalog_before,'portable_sha256':catalog_before}]})
    manifest['complete_results'].append({'source':'outputs/views/table.html','metadata':{'portable_sha256':browser_before,'source_sha256':native_sha}})
    (root/'report-manifest.json').write_text(json.dumps(manifest),encoding='utf-8')
    compact_report(root,html_min_bytes=0)
    actual=json.loads((root/'report-manifest.json').read_text(encoding='utf-8'))
    metadata=actual['panels'][-1]['metadata']
    assert metadata['portable_fit_index_sha256']==hashlib.sha256(catalog.read_bytes()).hexdigest()
    assert metadata['pre_storage_portable_fit_index_sha256']==catalog_before
    assert metadata['original_native_fit_index_sha256']==native_sha
    rich=actual['complete_results'][-1]['metadata']
    assert rich['portable_sha256']==hashlib.sha256(view.read_bytes()).hexdigest()
    assert rich['pre_storage_portable_sha256']==browser_before
    assert rich['source_sha256']==native_sha


def test_metadata_only_reconciliation_reuses_completed_unchanged_table_proof(tmp_path,monkeypatch):
    root, source, alias, view, raw, sha = _report(tmp_path)
    compact_report(root)
    import katana_next.reports.compact as module
    original_sha=module._sha
    def checked_sha(path):
        assert Path(path).suffix!='.tsv', 'Metadata reconciliation must not rehash already verified unchanged native GB tables'
        return original_sha(path)
    monkeypatch.setattr(module,'_sha',checked_sha)
    receipt=reconcile_compact_metadata(root,rehash_unchanged=False)
    assert receipt['unchanged_native_hashes'].startswith('retained')
    binding=json.loads((root/'report-manifest.json').read_text(encoding='utf-8'))['panels'][0]['data_links'][0]
    assert binding['source_sha256']==binding['portable_sha256']==sha
