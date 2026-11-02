"""Check a synthetic report's files and presentation contracts, not its science."""
from pathlib import Path
from html.parser import HTMLParser
from urllib.parse import urlsplit, unquote
import argparse
import hashlib
import json
import zipfile


class Links(HTMLParser):
    def __init__(self):
        super().__init__(); self.refs=[]; self.ids=[]
    def handle_starttag(self, tag, attrs):
        a=dict(attrs)
        if a.get('id'): self.ids.append(a['id'])
        for key in ('src','href'):
            if a.get(key): self.refs.append(a[key])


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--directory',type=Path,required=True)
    ap.add_argument('--receipt',type=Path,required=True)
    a=ap.parse_args(); root=a.directory.resolve()
    manifest=json.loads((root/'report-manifest.json').read_text('utf8'))
    assert manifest['dummy'] is True
    assert manifest['figure_data_links']['main_figures_with_plot_data']==len(manifest['panels'])
    assert not manifest['figure_data_links']['unbound_main_figures']
    assert len(manifest['contrasts'])==3
    checked=0; missing=[]; escaping=[]
    for page in root.rglob('*.html'):
        if page.relative_to(root).as_posix()=='assets/table-template.html':
            continue  # renderer template, not a published navigation target
        parser=Links(); parser.feed(page.read_text('utf8'))
        if page.name=='REPORT.html' and page.parent==root:
            assert all(f'chapter-{i}' in parser.ids for i in range(1,8))
            assert 'SYNTHETIC DISPLAY REVIEW' in page.read_text('utf8')
        for ref in parser.refs:
            u=urlsplit(ref)
            if u.scheme or u.netloc or not u.path: continue
            target=(page.parent/unquote(u.path)).resolve(); checked+=1
            if not target.is_relative_to(root): escaping.append((str(page),ref))
            elif not target.is_file(): missing.append((str(page),ref))
    assert not missing,missing
    assert not escaping,escaping
    summaries={m['contrast_id']:m['sequence_report'] for m in manifest['motif_summaries']}
    assert [summaries[c]['displayed_models'] for c in manifest['contrasts']]==[10,5,0]
    assert summaries['C_vs_A']['window_significant_models']==0
    assert summaries['C_vs_A']['regional_significant_regions']==5
    assert all(m['regional_FWER_is_window_FDR'] is False for m in manifest['motif_summaries'])
    assert all(m['native_statistics_unchanged'] for m in summaries.values())
    for t in manifest['top10']:
        assert len(t['rows'])==10
        assert (root/t['full_path']).is_file()
    assert len(manifest['network_browsers'])==1 and manifest['network_browsers'][0]['whole_study']
    archive=root/'REPORT.zip'
    with zipfile.ZipFile(archive) as z:
        assert z.testzip() is None
        files=[p for p in root.rglob('*') if p.is_file() and p!=archive]
        assert set(z.namelist())=={p.relative_to(root).as_posix() for p in files}
        assert all(z.read(p.relative_to(root).as_posix())==p.read_bytes() for p in files)
        members=len(files)
    result=dict(status='passed',scope='synthetic presentation and portable files only',
        science_acceptance=False,native_pipeline_acceptance=False,
        chapters=7,comparisons=3,panels=len(manifest['panels']),
        panels_with_data=len(manifest['panels']),top10_tables=len(manifest['top10']),
        complete_output_records=len(manifest['complete_results']),local_links_checked=checked,
        missing_local_files=missing,links_outside_report=escaping,
        motif_models_displayed=[10,5,0],network_browsers=1,zip_members=members,
        zip_bytes=archive.stat().st_size,zip_sha256=hashlib.sha256(archive.read_bytes()).hexdigest())
    a.receipt.write_text(json.dumps(result,indent=2),encoding='utf8')
    print(json.dumps(result))


if __name__=='__main__': main()
