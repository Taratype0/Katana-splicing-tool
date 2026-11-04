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
    target.write_text(text, encoding='utf-8')
    return {'version': 'katana-all-module-browser-v1', 'display_chapter': 7,
            'modules': len(data['modules']), 'native_data_unchanged': True,
            'association_filter_unchanged': True, 'new_scientific_computation': False}
