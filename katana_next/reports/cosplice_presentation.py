"""Repair presentation in an independent copy of a completed CoSplice report.

Association filtering remains unchanged. Module inspection does not require a
significant association. No network, scores, P values or FDR are recalculated.
"""
from __future__ import annotations
import json
import re
from pathlib import Path


def refresh_cosplice_copy(target):
    target = Path(target)
    text = target.read_text(encoding='utf-8')
    if 'katana-all-module-browser-v1' in text:
        return None
    payload = re.search(r'<script type="application/json" id="reportData">(.*?)</script>', text, re.S)
    if not payload:
        raise ValueError('Completed CoSplice browser lacks its recorded result payload')
    data = json.loads(payload[1])
    # Keep persisted filenames and fragment IDs compatible with old downloads.
    text = text.replace('Chapter 06', 'Chapter 07').replace('CHAPTER 06', 'CHAPTER 07')
    for n in range(1, 5):
        text = text.replace(f'<b>6.{n}</b>', f'<b>7.{n}</b>')
        text = text.replace(f'<span class="chapter-number">6.{n}</span>',
                            f'<span class="chapter-number">7.{n}</span>')
    text = text.replace('Module in selected evidence view', 'Module (all discovered modules)')
    text = text.replace('<details><summary>Optional: inspect connections between splicing events</summary>',
                        '<details open><summary>Network connections between splicing events</summary>')
    text = text.replace('How does a selected module change over time?',
                        'How does a selected module vary across conditions?')
    text = text.replace('Choose a module to see its overall time response.',
                        'Choose any discovered module to inspect its sample-level measurements.')
    # Show the actual saved score coverage. Missing scores are not null effects.
    old_sample = "renderTable('sampleTable',D.samples,[{key:'sample_id',label:'Sample'},{key:'biological_sample_id',label:'Biological unit'},{key:'condition',label:'Condition'}])"
    new_sample = "renderTable('sampleTable',D.samples.map(s=>({...s,score_coverage:D.modules.filter(m=>num(bySample[s.sample_id]?.[m.module])).length+' / '+D.modules.length})),[{key:'sample_id',label:'Sample'},{key:'biological_sample_id',label:'Biological unit'},{key:'condition',label:'Condition'},{key:'score_coverage',label:'Available module scores'}])"
    text = text.replace(old_sample, new_sample)
    text = text.replace('<div id="sampleTable"', '<p class="thin">Score coverage is counted from saved module scores. Missing scores can make associations not estimable; they are not evidence of no difference.</p><div id="sampleTable"')
    start = text.find("$('moduleSelect').replaceChildren();for(let r of significant)")
    end = text.find('\nfunction network()', start)
    if start < 0 or end < 0:
        raise ValueError('Unknown CoSplice module selector; presentation repair requires review')
    text = text[:start] + """const browse=scopeRows();
$('moduleSelect').replaceChildren();for(const r of browse){const o=document.createElement('option');o.value=r.module;o.textContent=r.module+' · association FDR '+fmt(r.q_value);$('moduleSelect').append(o)}
$('moduleSelect').disabled=!browse.length;$('noModules').hidden=!!browse.length;$('noModules').textContent='No discovered module is available in the recorded scope.';$('moduleBody').hidden=!browse.length;
if(browse.length)showModule(browse.some(r=>r.module===selectedModule)?selectedModule:browse[0].module);
else{renderTable('pathwayTable',[],[],{empty:'No recorded module available.'})}}
""" + text[end:]
    text = text.replace('r=significant.find(x=>x.module===selectedModule)',
                        'r=scopeRows().find(x=>x.module===selectedModule)')
    # Gene lists and enrichment follow the selected module, avoiding 122 long
    # sections. The association tables/heatmap still use the original gate.
    selected = 'scopeRows().filter(r=>r.module===selectedModule)'
    text = text.replace('for(const r of significant)', 'for(const r of '+selected+')')
    text = text.replace("if(!significant.length)host.innerHTML=", "if(!"+selected+".length)host.innerHTML=")
    text = text.replace("renderModule();if(scroll)",
                        "renderModule();renderIndependentLists();renderEnrichment39();if(scroll)")
    text = text.replace('Each module has its own Top 10 gene list below.',
                        'The selected module has its own Top 10 gene list below. All discovered modules remain accessible regardless of association significance.')
    text = text.replace('in the selected evidence view. Contributing genes',
                        'for the selected module. Contributing genes')
    text = text.replace('</body>', '<!-- katana-all-module-browser-v1 --></body>')
    if 'katana-all-module-browser-v1' not in text:
        text += '\n<!-- katana-all-module-browser-v1 -->'
    after = re.search(r'<script type="application/json" id="reportData">(.*?)</script>', text, re.S)
    if after[1] != payload[1]:
        raise ValueError('Presentation must preserve the complete native data payload byte-for-byte')
    from .compact import _atomic_write
    _atomic_write(target, text.encode('utf-8'))
    return {'version': 'katana-all-module-browser-v1', 'display_chapter': 7,
            'modules': len(data['modules']), 'native_data_unchanged': True,
            'association_filter_unchanged': True, 'new_scientific_computation': False}


def refresh_cosplice_gene_labels(target):
    """Use saved member annotations in display and copied module downloads.

    Native JSON, object IDs, graph edges, ranking and tests remain untouched.
    Atomic replacement also protects source reports copied using hard links.
    """
    from collections import defaultdict
    import io
    import zipfile
    import pandas as pd
    from .compact import _atomic_write
    target = Path(target)
    text = target.read_text('utf-8')
    marker = 'katana-cosplice-gene-labels-v1'
    if marker in text: return None
    payload = re.search(r'<script type="application/json" id="reportData">(.*?)</script>', text, re.S)
    if not payload: raise ValueError('CoSplice gene display requires its completed native payload')
    data = json.loads(payload[1])
    candidates = defaultdict(set)
    for row in [*data.get('members', []), *data.get('genes', [])]:
        gid, symbol = str(row.get('gene_id') or ''), str(row.get('gene_symbol') or '')
        if gid and symbol and symbol.lower() not in {'nan', 'none', 'na', 'unannotated'} and symbol != gid and not symbol.startswith('ENS'):
            candidates[gid].add(symbol)
    names = {gid: next(iter(values)) for gid, values in candidates.items() if len(values) == 1}
    # A display overlay is separate from reportData, which is kept byte exact.
    anchor = "const moduleMap="
    position = text.find(anchor, payload.end())
    if position < 0: raise ValueError('Unknown CoSplice gene browser template')
    overlay = ('const katanaGeneNames=' + json.dumps(names, ensure_ascii=True).replace('<', '\\u003c') + ';\n'
               "const geneLabel=r=>katanaGeneNames[r?.gene_id]||r?.gene_id||'Unannotated';\n"
               "D.genes=D.genes.map(r=>({...r,gene_symbol:geneLabel(r)}));\n"
               "D.members=D.members.map(r=>({...r,gene_symbol:geneLabel(r)}));\n")
    text = text[:position] + overlay + text[position:]
    # Replace display expressions only. Clicks, joins, ordering and activeGene
    # continue to use stable IDs; searching accepts either symbol or ID.
    changes = {
        "{key:'gene_id',label:'Gene'}": "{key:'gene_symbol',label:'Gene'},{key:'gene_id',label:'Gene ID'}",
        "g.gene_id.toLowerCase().includes(q)": "(geneLabel(g)+' '+g.gene_id).toLowerCase().includes(q)",
        "(e.gene_id+' '+e.event_key).toLowerCase()": "(geneLabel(e)+' '+e.gene_id+' '+e.event_key).toLowerCase()",
        "b.textContent=g.gene_id+' · '": "b.textContent=geneLabel(g)+' · '",
        "E(e?.gene_id||'Unannotated')": "E(geneLabel(e))",
        "(e?.gene_id||'Unannotated')+' · kME '": "geneLabel(e)+' · '+(e?.gene_id||'Unannotated')+' · kME '",
        "r=>r.gene_id+' · E'": "r=>geneLabel(r)+' · E'",
        "o.textContent=e.gene_id+' · '": "o.textContent=geneLabel(e)+' · '+e.gene_id+' · '",
        "e.gene_id+' observed PSI'": "geneLabel(e)+' observed PSI'",
        "textContent=e.gene_id+' · '": "textContent=geneLabel(e)+' · '+e.gene_id+' · '",
    }
    # Apply only to executable script, never the stored scientific payload.
    head, script = text[:payload.end()], text[payload.end():]
    for old, new in changes.items(): script = script.replace(old, new)
    script = script.replace("map(g=>g.gene_id).join(', ')", "map(g=>geneLabel(g)).join(', ')")
    text = head + script + '\n<!-- ' + marker + ' -->'
    updated = {}
    download_sets = [*data.get('module_files', {}).values(),
                     {'genes': 'module_gene_summary.tsv', 'events': 'gene_event_membership.tsv'}]
    for files in download_sets:
        for role in ('genes', 'events'):
            if not files.get(role): continue
            path = (target.parent / 'report_tables' / files[role]).resolve()
            if not path.is_relative_to(target.parent.resolve()):
                raise ValueError('CoSplice member download escapes report directory')
            if not path.is_file(): continue
            frame = pd.read_csv(path, sep='\t', dtype=str, keep_default_na=False)
            if 'gene_id' not in frame: continue
            frame['gene_symbol'] = frame.gene_id.map(names).fillna(frame.gene_id)
            content = frame.to_csv(sep='\t', index=False).encode('utf-8')
            _atomic_write(path, content)
            updated[path.name] = content
    # Existing ZIPs are report downloads; augment only matching tables and
    # preserve every other entry. The original source ZIP remains unchanged.
    for path in (target.parent / 'report_tables').glob('*module_lists.zip'):
        buffer = io.BytesIO()
        with zipfile.ZipFile(path) as original, zipfile.ZipFile(buffer, 'w', zipfile.ZIP_DEFLATED) as archive:
            for entry in original.infolist():
                archive.writestr(entry, updated.get(Path(entry.filename).name, original.read(entry)))
        _atomic_write(path, buffer.getvalue())
    if re.search(r'id="reportData">(.*?)</script>', text, re.S)[1] != payload[1]:
        raise ValueError('CoSplice native payload changed while labelling genes')
    _atomic_write(target, text.encode('utf-8'))
    return {'version': marker, 'mapped_gene_ids': len(names),
            'module_gene_rows': len(data.get('genes', [])),
            'mapped_module_gene_rows': sum(r.get('gene_id') in names for r in data.get('genes', [])),
            'ambiguous_gene_ids': sum(len(v) > 1 for v in candidates.values()),
            'download_tables': len(updated), 'native_data_unchanged': True,
            'new_scientific_computation': False}
