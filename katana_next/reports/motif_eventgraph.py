"""Presentation only: keep regional FWER and window BH results separate."""
import html
import json
from pathlib import Path


def render_eventgraph_report(summary, clusters, destination):
    destination = Path(destination)
    from .motif_sequences import standalone_fragment
    sequence_section = standalone_fragment(summary, clusters, destination)
    files = sorted(p.name for p in destination.parent.iterdir() if p.is_file() and p.suffix != '.html')
    files.append('COMPLETE.json')
    sig = clusters[clusters.cluster_P_FWER_new_family.le(.05)]
    sections = []
    columns = ['model_id', 'RBP', 'event_type', 'direction', 'region', 'window_nt',
               'covered_region_start0', 'covered_region_end0', 'minimum_events',
               'minimum_gene_blocks', 'cluster_P_FWER_new_family']
    for layer, label in [('affected_sequence', 'Motif consequences in affected RNA'),
                         ('regulatory_pre_mRNA', 'Regulatory pre-mRNA flanks')]:
        rows = sig[sig.layer.eq(layer)]
        body = ('<p>No region meets whole-family FWER P ≤ 0.05.</p>' if rows.empty else
                rows[columns].to_html(index=False, escape=True, float_format=lambda x: f'{x:.6g}', border=0))
        sections.append(f'<h3>{label}</h3><div class="scroll">{body}</div>')
    payload = json.dumps({'summary': summary, 'available_files': files}, ensure_ascii=True, allow_nan=False).replace('<', '\\u003c')
    page = '''<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>MotifTrace | Event-level consequences</title><style>
body{font:16px/1.6 Arial,sans-serif;color:#1e2630;background:#fff;max-width:1180px;margin:32px auto;padding:0 24px}
h1,h2,h3,a{color:#2863A5}.note{border-left:4px solid #C63F4A;padding:12px}.scroll{overflow:auto}
td,th{padding:9px;border-bottom:1px solid #d9dfe7;text-align:left}table{border-collapse:collapse;font-size:13px;white-space:nowrap}
section{margin:28px 0}pre{white-space:pre-wrap;overflow-wrap:anywhere;background:#f5f7fa;padding:16px}</style>
<h1>MotifTrace: event-level RNA consequences</h1><h2>''' + html.escape(str(summary['comparison_id'])) + '''</h2>
<p>All eligible foreground events remain in the matching ledger. Matched events from the same gene share one randomization sign. Unmatched events and reasons remain downloadable.</p>
<p class="note">Regional significance does not imply that individual windows or bases are significant. Sequence opportunity is not measured binding or evidence of causal regulation.</p>
''' + sequence_section + '''<section id="window-significance"><h2>Window / full-sequence BH FDR</h2><p>''' + str(summary['new_window_q_lt_0_05']) + ''' tests have BH FDR ≤ 0.05 out of ''' + str(summary['eligible_new_window_tests']) + ''' estimable tests. The correction family contains all declared models, event types, directions, regions and scalar/window tests in this comparison.</p>
<p><a download data-file="NEW_WINDOW_AND_FULL_SEQUENCE_TESTS.tsv.gz">Complete window and full-sequence P / BH FDR table</a>. Column: <code>q_new_window_family</code>. Regional P values are never substituted into this column.</p></section>
<section id="regional-significance"><h2>Regional significance: whole-family FWER P</h2><p>''' + str(summary['cluster_P_FWER_lt_0_05']) + ''' regions meet FWER P ≤ 0.05. The maximum-cluster correction spans both biological layers together; the sections below are display groups, not separate corrections.</p>
<p><a download data-file="NEW_CLUSTER_FWER.tsv">Complete regional FWER P table</a>. Column: <code>cluster_P_FWER_new_family</code>; this is a corrected P value, not FDR.</p>''' + ''.join(sections) + '''</section>
<section><h2>Sequences, event evidence and complete downloads</h2><p>Actual RNA forms, gained/lost/retained motif loci, PSI-weighted opportunities, matching losses and all statistical results are retained separately below. Motif sites and per-event scores do not have individual significance inferred from regional results.</p><div id="files"></div></section>
<section><h2>Analysis record</h2><pre id="summary"></pre></section><script>
const D=''' + payload + ''';
document.getElementById('summary').textContent=JSON.stringify(D.summary,null,2);
for(const a of document.querySelectorAll('a[data-file]')){const n=a.dataset.file;a.href=encodeURIComponent(n)}
for(const img of document.querySelectorAll('img[data-file]')){const n=img.dataset.file;img.src=encodeURIComponent(n)}
for(const n of D.available_files){const p=document.createElement('p'),a=document.createElement('a');a.href=encodeURIComponent(n);a.download=n;a.textContent=n;p.append(a);document.getElementById('files').append(p)}
</script></html>'''
    destination.write_text(page, encoding='utf-8')
