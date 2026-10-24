import base64
import gzip
import hashlib
import json
import os
from pathlib import Path
import re

import pytest

from katana_next.reports.compact import compact_report, reconcile_compact_metadata


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
    assert 'outputs/results/native.tsv?download=1#full' in (root/'REPORT.html').read_text()
    native_index=(root/'outputs/views/native.json').read_text()
    assert '"path":"../results/native.tsv"' in native_index
    assert '"p":1e-200' in native_index
    manifest=json.loads((root/'report-manifest.json').read_text())
    binding=manifest['panels'][0]['data_links'][0]
    assert binding['source_sha256'] == binding['sha256'] == binding['portable_sha256'] == sha
    assert binding['path'] == 'outputs/results/native.tsv'
    assert manifest['complete_results'][0]['source'] == 'outputs/results/native.tsv'
    content=view.read_text()
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
    compressed=re.search(r'<script id="table-data"[^>]*>(.*?)</script>',view.read_text(),re.S).group(1)
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
    assert provenance.read_text()==original
    assert '"outputs/comparisons/B_vs_A/same.tsv":' in index.read_text()


def test_generated_report_provenance_portable_bindings_are_rebound(tmp_path):
    root, source, alias, view, raw, sha = _report(tmp_path)
    provenance=root/'provenance.json'
    provenance.write_text(json.dumps({'schema':'katana.report_provenance/v1',
        'source_files':[{'source_name':'same.tsv','source_sha256':sha,'bytes':source.stat().st_size,'report_path':'outputs/comparisons/B_vs_A/same.tsv'}],
        'figure_data_links':{'bindings':[{'path':'outputs/comparisons/B_vs_A/same.tsv','source_sha256':sha,'portable_sha256':sha}]}}),encoding='utf-8')
    compact_report(root,html_min_bytes=0)
    actual=json.loads(provenance.read_text())
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
    binding=json.loads((root/'provenance.json').read_text())['figure_data_links']['bindings'][0]
    assert binding['path']=='outputs/results/native.tsv'
    assert binding['source_sha256']==binding['portable_sha256']==sha
    assert receipt['restored_native_json'][0]['restored_source_sha256']==hashlib.sha256(native_bytes).hexdigest()
