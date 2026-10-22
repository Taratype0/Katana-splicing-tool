"""Display-only RBP evidence taxonomy; never combines or recalculates results."""
from __future__ import annotations
import html
import json
from pathlib import Path
CONTRACT_PATH = Path(__file__).with_suffix('.json')
CONTRACT = json.loads(CONTRACT_PATH.read_text(encoding='utf-8'))
LAYERS = tuple(CONTRACT['layers'])

def evidence_layer(workflow_id, parameters=None, metadata=None):
    parameters, metadata = (parameters or {}, metadata or {})
    if workflow_id == 'prioritize.rbp_mechanism':
        scope = parameters.get('sequence_scope', 'boundary_association')
        key = CONTRACT['sequence_scope_layers'].get(scope)
    else:
        key = CONTRACT['workflow_layers'].get(workflow_id)
    return next((layer for layer in LAYERS if layer['id'] == key), None)

def artifact_layer(artifact):
    return evidence_layer(artifact.analysis_id, artifact.parameters, artifact.metadata)

def report_groups(artifacts):
    """Do not merge modes/contrasts under the same durable RBP entry."""
    groups = {}
    for artifact in artifacts:
        layer = artifact_layer(artifact)
        key = (artifact.analysis_id, layer['id'] if layer else None, tuple(artifact.contrast_ids) if layer else (), artifact.parameters.get('sequence_scope') if layer else None)
        groups.setdefault(key, []).append(artifact)
    return [(key[0], next((l for l in LAYERS if l['id'] == key[1]), None), values) for key, values in groups.items()]

def layer_inventory(artifacts):
    artifacts = list(artifacts)
    return [dict(layer, artifact_ids=[a.artifact_id for a in artifacts if (artifact_layer(a) or {}).get('id') == layer['id']]) for layer in LAYERS]

def overview_html(artifacts=(), language='EN', inventory=True):
    if str(language).upper() != 'EN':
        raise ValueError('Only English reports are supported')
    suffix = 'en'
    parts = [f"""<section id="rbp-evidence-layers"><h2>{html.escape(CONTRACT['title_' + suffix])}</h2>"""]
    for layer in layer_inventory(artifacts):
        parts.append(f'''<article data-rbp-layer="{layer['id']}"><h3>{layer['number']}. {html.escape(layer['title_' + suffix])}</h3>''')
        for field in ('question_', 'value_', 'unit_', 'boundary_'):
            parts.append('<p>' + html.escape(layer[field + suffix]) + '</p>')
        if inventory:
            count = len(layer['artifact_ids'])
            status = f'{count} saved outputs included; output files are not gene/event counts.' if count else 'No saved output of this layer is included; not a negative result and no analysis was started.'
            parts.append('<p class="rbp-layer-status">' + status + '</p>')
        parts.append('</article>')
    parts.append('<p>' + html.escape(CONTRACT['shared_boundary_' + suffix]) + '</p></section>')
    return ''.join(parts)

def layer_note_html(layer, language='EN'):
    if str(language).upper() != 'EN':
        raise ValueError('Only English reports are supported')
    if layer is None:
        return ''
    suffix = 'en'
    return f'''<p class="rbp-layer-label" data-rbp-layer="{layer['id']}"><strong>{layer['number']}. {html.escape(layer['title_' + suffix])}</strong></p><p>{html.escape(layer['boundary_' + suffix])}</p>'''
