"""Render verified completed MotifTrace families; never scan RNA or run tests.

Numeric PWM logos and count/position summaries are display transformations.
Native raw P, BH q and cluster FWER families remain complete and unchanged.
"""
from __future__ import annotations

from datetime import datetime, timezone
import html
import json
from pathlib import Path

import numpy as np
import pandas as pd

from .comparisons import BLUE, RED, _safe
from .motif_display import motif_display_summary, significant_motif_ranking
from .sources import digest

VERSION = 'completed-motif-report-v3'
SEQUENCE_FILES = ('ANALYSIS_PLAN_BEFORE_SCORING.json', 'FROZEN_SELECTED_EVENTS.tsv',
                  'PAIRED_FORM_INSTANCE_CHANGES.tsv.gz', 'ACTUAL_RNA_FORMS.tsv.gz',
                  'ACTUAL_LOCUS_HITS.tsv.gz', 'PAIRED_FORM_NORMALIZED_SCORES.tsv.gz')
STAT_FILES = ('NEW_WINDOW_AND_FULL_SEQUENCE_TESTS.tsv.gz', 'NEW_CLUSTER_FWER.tsv', 'MAX_CLUSTER_NULL.npy')
TEST_COLUMNS = ['event_type', 'direction', 'model_id', 'RBP', 'region', 'method', 'window_nt',
                'window_start0', 'n_matched_gene_pairs', 'paired_rank_Z', 'mean_paired_score_difference',
                'raw_P_new_paired_swaps', 'q_new_window_family', 'status']
RANK_COLUMNS = ['q_new_window_family', 'raw_P_new_paired_swaps', 'event_type',
                'direction', 'region', 'method', 'window_nt', 'window_start0']


def _best_recorded_tests(chunk, existing):
    """Select one actual row per model, keeping its P/q/effect paired.

    This is a descriptive ranking across already corrected tests, not a new
    model-level hypothesis or a new correction family.
    """
    for column in ('q_new_window_family', 'raw_P_new_paired_swaps'):
        values = chunk.loc[chunk[column].notna(), column]
        if not np.isfinite(values).all() or not values.between(0., 1.).all():
            raise ValueError('Recorded native motif P/FDR must be finite in [0, 1] or explicitly missing: ' + column)
    valid = chunk.loc[chunk.q_new_window_family.notna() & chunk.raw_P_new_paired_swaps.notna()]
    for token, part in valid.groupby('model_id', sort=False):
        candidate = part.sort_values(RANK_COLUMNS, kind='stable').iloc[0]
        key = tuple(candidate[column] for column in RANK_COLUMNS)
        if token not in existing or key < existing[token][0]:
            existing[token] = key, candidate.to_dict()


def _motif_ranking(summary, best, models):
    records = []
    for item in summary.to_dict('records'):
        token = item['model_id']
        record = dict(item)
        record.update({column: np.nan for column in TEST_COLUMNS if column not in record})
        if token in best:
            record.update(best[token][1])
        record['model_id'] = token
        record['RBP'] = str(models[token].get('RBPs', token))
        matrix = _matrix(models[token])
        record['motif_sequence'] = (''.join('ACGU'[i] for i in matrix.argmax(axis=1))
                                    if matrix is not None else _recorded_pattern(models[token]) or 'not_available')
        record['sequence_representation'] = ('maximum-probability PWM consensus' if models[token].get('matrix') is not None
                                              else 'literal' if matrix is not None else 'recorded IUPAC/fixed-pattern consensus')
        record['test_scope'] = ('not_estimable' if token not in best else
                                'full_sequence_test' if record['window_nt'] == 0 else 'window_test')
        q = record['q_new_window_family']
        record['significance'] = ('not_estimable' if pd.isna(q) else
                                  'recorded_test_FDR_le_0.05' if q <= .05 else 'not_significant')
        record['correction_scope'] = 'complete comparison window/full-sequence test family; not model-level FDR'
        records.append(record)
    result = pd.DataFrame(records).sort_values(
        ['q_new_window_family', 'raw_P_new_paired_swaps', 'model_id'], kind='stable', na_position='last')
    result.insert(0, 'display_rank', range(1, len(result) + 1))
    return result


def verify_completed_motif_family(source, completion_proof, *, comparison_id, proof_kind='native', native_job_id=None):
    """Verify every manifest member and an explicit terminal execution receipt.

    Native receipts bind actual sacct output, Nextflow trace, terminal claim,
    and the downloaded collection ledger. A caller must also supply the
    expected native_job_id. The explicit verified_local route accepts the
    original local completion receipt (passed, source_unchanged, zero upstream
    refits, exact summary comparison, and the same complete-manifest hash).
    Queue absence, an output directory or ``passed`` alone is insufficient.
    """
    source, proof_path = Path(source).resolve(), Path(completion_proof).resolve()
    marker = source / 'COMPLETE.json'
    complete = json.loads(marker.read_text('utf-8-sig'))
    proof = json.loads(proof_path.read_text('utf-8-sig'))
    identity = digest(marker)
    binding = proof.get('complete_manifest', {})
    if not complete.get('complete') or ('statistical_family_split' in complete and complete['statistical_family_split'] is not False):
        raise ValueError('A complete unsplit scientific MotifTrace family is required.')
    if binding.get('sha256', binding.get('SHA256')) != identity:
        raise ValueError('Completion receipt does not bind this exact scientific manifest.')
    if proof_kind == 'native':
        job = proof.get('native_job', {})
        if (proof.get('schema') != 'katana.native_motiftrace_completion/v1' or proof.get('passed') is not True
                or native_job_id is None or str(proof.get('native_job_id')) != str(native_job_id)
                or str(job.get('job_id')) != str(native_job_id)
                or job.get('state') != 'COMPLETED' or job.get('exit_code') != '0:0'
                or proof.get('comparison_id') != comparison_id):
            raise ValueError('Exact native COMPLETED/0:0 job and comparison proof are required.')
        sacct, trace, claim = (proof.get(key, {}) for key in ('sacct', 'nextflow_trace', 'claim_completed'))
        trace_row = trace.get('actual_row', {})
        command = sacct.get('actual_command', {})
        rows = [line.strip().split('|') for line in str(command.get('stdout', '')).splitlines() if line.strip()]
        exact = [row for row in rows if row[0].strip() == str(native_job_id)]
        if (sacct.get('state') != 'COMPLETED' or sacct.get('exit_code') != '0:0'
                or command.get('returncode') != 0 or not command.get('argv') or command['argv'][0] != 'sacct'
                or len(exact) != 1 or len(exact[0]) < 3 or exact[0][1:3] != ['COMPLETED', '0:0']
                or str(sacct.get('main_row', {}).get('job_id')) != str(native_job_id)
                or sacct['main_row'].get('state') != 'COMPLETED' or sacct['main_row'].get('exit_code') != '0:0'
                or trace.get('status') != 'COMPLETED' or trace.get('exit') != 0
                or str(trace.get('native_id')) != str(native_job_id)
                or trace_row.get('status') != 'COMPLETED' or str(trace_row.get('exit')) != '0'
                or str(trace_row.get('native_id')) != str(native_job_id)
                or not str(trace_row.get('name', '')).rstrip().endswith('(' + comparison_id + ')')
                or claim.get('passed') is not True or str(claim.get('job_id')) != str(native_job_id)
                or claim.get('comparison_id') != comparison_id):
            raise ValueError('Native scientific completion requires matching actual sacct, trace and terminal claim evidence.')
        if (proof.get('scope') != 'collection_only' or proof.get('collector_RNA_rescans') != 0
                or proof.get('upstream_refits') != 0 or proof.get('family_split') is not False
                or binding.get('relative_path') != 'COMPLETE.json' or binding.get('bytes') != marker.stat().st_size
                or binding.get('members_verified') is not True):
            raise ValueError('Native collection scope and complete unsplit family must be explicitly verified.')
        collection_path = source / 'COLLECTION.receipt.json'
        collection = json.loads(collection_path.read_text('utf-8-sig'))
        collected_marker = collection.get('complete_manifest', {})
        collected_proof = collection.get('native_completion_receipt', {})
        if (collection.get('schema') != 'katana.motiftrace_collection/v1' or collection.get('passed') is not True
                or collection.get('comparison_id') != comparison_id or str(collection.get('native_job_id')) != str(native_job_id)
                or collection.get('all_scientific_members_unchanged') is not True
                or collected_marker.get('relative_path') != 'COMPLETE.json'
                or any(collected_marker.get(key) != identity for key in ('sha256', 'sha256_before', 'sha256_after'))
                or collected_marker.get('members_verified') is not True or collected_marker.get('bytes') != marker.stat().st_size
                or collected_proof.get('relative_path') != proof_path.name
                or collected_proof.get('sha256') != digest(proof_path) or collected_proof.get('bytes') != proof_path.stat().st_size
                or collection.get('members') != complete.get('members')):
            raise ValueError('Downloaded native collection ledger does not bind every original scientific member and completion receipt.')
    elif proof_kind == 'verified_local':
        if (proof.get('passed') is not True or proof.get('source_unchanged') is not True
                or proof.get('upstream_refits') != 0 or proof.get('summary', {}).get('comparison_id') != comparison_id):
            raise ValueError('Verified local scientific completion proof is incomplete or belongs to another comparison.')
    elif proof_kind == 'flow':
        from katana_next.flow.presentation import validate_motif_view_proof
        validate_motif_view_proof(proof_path, source, comparison_id)
    else:
        raise ValueError('Choose an explicit native, verified_local or flow completion-proof route.')
    members = complete.get('members', {})
    missing = set((*SEQUENCE_FILES, *STAT_FILES)) - members.keys()
    if missing:
        raise ValueError('Scientific manifest omits required outputs: ' + ', '.join(sorted(missing)))
    audit = []
    for name, record in members.items():
        path = (source / name).resolve()
        if not path.is_relative_to(source) or not path.is_file():
            raise ValueError('Missing or unbounded completed MotifTrace member: ' + name)
        stat = path.stat()
        if stat.st_size != record['bytes'] or digest(path) != record['sha256']:
            raise ValueError('Completed MotifTrace member SHA/size mismatch: ' + name)
        audit.append({'path': str(path), 'bytes': stat.st_size, 'mtime_ns': stat.st_mtime_ns, 'sha256': record['sha256']})
    plan = json.loads((source / 'ANALYSIS_PLAN_BEFORE_SCORING.json').read_text('utf-8-sig'))
    # The frozen native writer stores this declaration in the hashed plan;
    # older verified completion markers also repeat it at their top level.
    declared_split = complete.get('statistical_family_split', plan.get('statistical_family_split', plan.get('provenance', {}).get('statistical_family_split')))
    if declared_split is not False:
        raise ValueError('The hashed completed scientific plan must explicitly establish an unsplit family.')
    if plan.get('comparison_id') != comparison_id or not isinstance(plan.get('models'), dict) or not plan['models']:
        raise ValueError('Completed scientific plan has another comparison or no declared model catalogue.')
    if proof_kind == 'native' and (proof.get('model_count') != len(plan['models']) or proof.get('permutations') != plan.get('permutations')):
        raise ValueError('Native receipt model/permutation family does not match the complete scientific plan.')
    return {'manifest_sha256': identity, 'manifest': str(marker), 'proof': str(proof_path),
            'proof_sha256': digest(proof_path), 'proof_kind': proof_kind, 'native_job_id': native_job_id,
            'members': audit, 'plan': plan}


def _read(path, **kwargs):
    return pd.read_csv(path, sep='\t', float_precision='round_trip', dtype={'model_id': str, 'event_id': str}, **kwargs)


def verify_empty_completed_statistics(source, members):
    """A complete zero-pair family has sequences but no estimable tests.

    The native writer retains a minimal header for an empty test table. Accept
    it only with its authenticated full matching audit and native result
    summary; absence of a file or a positive matched-pair count is not N/A.
    """
    source = Path(source)
    needed = {'RESULTS.json', 'MATCHING_ALL_SELECTED_WITH_LOSSES.tsv'}
    if not needed <= set(members):
        raise ValueError('Empty native MotifTrace tests require their authenticated summary and complete matching audit')
    summary = json.loads((source / 'RESULTS.json').read_text('utf-8-sig'))
    matching = _read(source / 'MATCHING_ALL_SELECTED_WITH_LOSSES.tsv')
    if (summary.get('matched_independent_gene_pairs') != 0 or summary.get('eligible_new_window_tests') != 0 or
            summary.get('inference_status') != 'not_estimable_insufficient_independent_matched_pairs' or
            'status' not in matching or matching.status.eq('matched').any() or
            len(matching) != summary.get('selected_events')):
        raise ValueError('Empty native MotifTrace tests contradict their complete gene-pair matching evidence')


def _test_chunks(source, members):
    path = Path(source) / 'NEW_WINDOW_AND_FULL_SEQUENCE_TESTS.tsv.gz'
    header = _read(path, nrows=0)
    if not set(TEST_COLUMNS) <= set(header):
        frame = _read(path)
        if not frame.empty:
            raise ValueError('Nonempty native MotifTrace tests lack required statistical identity columns')
        verify_empty_completed_statistics(source, members)
        return iter([pd.DataFrame({name: pd.Series(dtype='float64' if name in {
            'window_nt', 'window_start0', 'n_matched_gene_pairs', 'paired_rank_Z',
            'mean_paired_score_difference', 'raw_P_new_paired_swaps', 'q_new_window_family'} else 'object')
                                  for name in TEST_COLUMNS})])
    return _read(path, chunksize=50000, usecols=TEST_COLUMNS, low_memory=False)


def _recorded_pattern(model):
    return str(model.get('sequence', model.get('motif', model.get('pattern', '')))).upper().replace('T', 'U')


def _matrix(model):
    value = model.get('matrix')
    if value is None:
        literal = _recorded_pattern(model)
        if not literal or any(c not in 'ACGU' for c in literal):
            return None
        return np.array([[float(base == c) for base in 'ACGU'] for c in literal])
    matrix = np.asarray(value, dtype=float)
    if matrix.ndim != 2 or matrix.shape[1] != 4 or not np.isfinite(matrix).all() or (matrix < 0).any() or (matrix.sum(axis=1) <= 0).any():
        raise ValueError('Recorded PWM must be a finite nonnegative position-by-four matrix.')
    return matrix / matrix.sum(axis=1, keepdims=True)


def _logo(ax, model):
    from matplotlib.textpath import TextPath
    from matplotlib.transforms import Affine2D
    from matplotlib.patches import PathPatch
    matrix = _matrix(model)
    if matrix is None:
        pattern = _recorded_pattern(model)
        label = (pattern + '\nRecorded IUPAC/fixed-pattern consensus; no base probabilities supplied'
                 if pattern else 'No numeric PWM/literal logo representation supplied')
        ax.text(.5, .5, label, ha='center', va='center', transform=ax.transAxes)
        ax.axis('off')
        return
    information = 2 + (matrix * np.log2(np.maximum(matrix, 1e-300))).sum(axis=1)
    colors = {'A': '#238B45', 'C': BLUE, 'G': '#D49A25', 'U': RED}
    for i, row in enumerate(matrix):
        bottom = 0.
        for j in np.argsort(row):
            height = float(row[j] * information[i])
            if height <= 1e-9:
                continue
            glyph = TextPath((0, 0), 'ACGU'[j], size=1, prop={'weight': 'bold'})
            bounds = glyph.get_extents()
            transform = Affine2D().translate(-bounds.x0, -bounds.y0).scale(.85 / bounds.width, height / bounds.height).translate(i + .075, bottom)
            ax.add_patch(PathPatch(glyph, transform=transform + ax.transData, color=colors['ACGU'[j]], linewidth=0))
            bottom += height
    ax.set(xlim=(0, len(matrix)), ylim=(0, 2), ylabel='bits')
    ax.set_xticks(np.arange(len(matrix)) + .5, np.arange(len(matrix)) + 1, fontsize=7)


def _catalogue_browser(models):
    """All models, paginated locally; no CDN, downloads, motif scan or Top N."""
    records = [{**model, 'model_id': key} for key, model in models.items()]
    data = json.dumps(records, allow_nan=False).replace('<', '\\u003c')
    return '''<!doctype html><html lang="en"><meta charset="utf-8"><title>All MotifTrace models</title>
<style>body{font:15px Arial;color:#1e2630;margin:28px}article{border:1px solid #ddd;padding:14px;margin:15px 0}pre{white-space:pre-wrap}button,input{padding:7px;margin:4px}svg{max-width:100%}</style>
<h1>Complete MotifTrace model catalogue</h1><p>Every declared model, native assignment and available recorded-test summary. The displayed order is descriptive: selected native test FDR, paired P, then model ID. Test-family FDR is not whole-model FDR. Logos describe recorded sequence preference, not measured binding.</p>
<input id="query" placeholder="Model or RBP"><button id="previous">Previous</button><button id="next">Next</button><span id="count"></span><main id="models"></main>
<script id="catalogue" type="application/json">''' + data + '''</script><script>
const all=JSON.parse(document.getElementById('catalogue').textContent);let page=0;const size=20;
function text(tag,value){let n=document.createElement(tag);n.textContent=value;return n}
function draw(){let query=document.getElementById('query').value.toLowerCase(),rows=all.filter(x=>JSON.stringify(x).toLowerCase().includes(query));page=Math.max(0,Math.min(page,Math.ceil(rows.length/size)-1));let root=document.getElementById('models');root.replaceChildren();
for(const model of rows.slice(page*size,(page+1)*size)){let box=document.createElement('article');box.append(text('h2',(model.RBPs||model.model_id)+' · '+model.model_id));
if(Array.isArray(model.matrix)){const ns='http://www.w3.org/2000/svg',svg=document.createElementNS(ns,'svg');svg.setAttribute('viewBox','0 0 '+(model.matrix.length*34)+' 130');svg.style.width=(model.matrix.length*34)+'px';
model.matrix.forEach((raw,i)=>{let sum=raw.reduce((a,b)=>a+b,0),row=raw.map(x=>x/sum),bits=2+row.reduce((a,b)=>a+(b?b*Math.log2(b):0),0),bottom=110;[0,1,2,3].sort((a,b)=>row[a]-row[b]).forEach(j=>{let height=row[j]*bits*45;if(height<.01)return;let t=document.createElementNS(ns,'text');t.textContent='ACGU'[j];t.setAttribute('x',i*34+3);t.setAttribute('y',bottom);t.setAttribute('font-family','Arial');t.setAttribute('font-weight','bold');t.setAttribute('font-size',height);t.setAttribute('fill',['#238B45','#2863A5','#D49A25','#C63F4A'][j]);svg.append(t);bottom-=height});});box.append(svg)}
box.append(text('pre',JSON.stringify(model,null,2)));root.append(box)}document.getElementById('count').textContent=rows.length+' / '+all.length+' models; page '+(page+1)+' / '+Math.max(1,Math.ceil(rows.length/size));}
document.getElementById('query').oninput=()=>{page=0;draw()};document.getElementById('previous').onclick=()=>{page--;draw()};document.getElementById('next').onclick=()=>{page++;draw()};draw();</script></html>'''


def build_completed_motif_views(source, destination, *, comparison_id, project_id,
                                source_project_manifest_sha256, completion_proof,
                                proof_kind='native', native_job_id=None, example_max_nt=450):
    """Return data-bound sequence/statistical panels from a verified full family.

    This renderer accepts completed outputs only. It never imports a scientific
    MotifTrace kernel, fits BH/FWER, scores a sequence, or submits a command.
    Profile/example selection is descriptive and explicitly recorded.
    """
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    evidence = verify_completed_motif_family(source, completion_proof, comparison_id=comparison_id,
                                               proof_kind=proof_kind, native_job_id=native_job_id)
    if not isinstance(example_max_nt, int) or example_max_nt <= 0:
        raise ValueError('RNA example display length must be a positive integer.')
    proof = json.loads(Path(completion_proof).read_text('utf-8-sig'))
    if proof.get('project_id', project_id) != project_id:
        raise ValueError('Native MotifTrace completion belongs to another project.')
    if (proof_kind == 'native'
            and proof.get('study_lineage', {}).get('source_project_manifest_SHA256') != source_project_manifest_sha256):
        raise ValueError('Native MotifTrace completion belongs to another source-project manifest lineage.')
    if proof_kind == 'flow' and proof.get('render_project_manifest_sha256') != source_project_manifest_sha256:
        raise ValueError('Flow MotifTrace completion belongs to another collected project manifest.')
    source, out = Path(source).resolve(), Path(destination).resolve()
    if out == source or out.is_relative_to(source) or source.is_relative_to(out):
        raise ValueError('Motif report output must be independent of the immutable scientific source.')
    if out.exists():
        raise FileExistsError('Existing MotifTrace report outputs are preserved; use a new destination.')
    out.mkdir(parents=True)
    plan, models = evidence['plan'], evidence['plan']['models']
    native_members = json.loads(Path(evidence['manifest']).read_text('utf-8-sig'))['members']
    panels, tables, unavailable = [], [], []
    stem = _safe(comparison_id)
    scope = {'mode': 'fixed_completed_result', 'recorded_event_set': 'stable_significant',
             'read_support_policy': plan.get('provenance', {}).get('read_support_policy')}
    methods = out / 'completed-motif-methods.json'
    methods.write_text(json.dumps({'schema': VERSION, 'comparison_id': comparison_id,
                                  'complete_manifest_sha256': evidence['manifest_sha256'],
                                  'completion_proof_sha256': evidence['proof_sha256'], 'proof_kind': proof_kind,
                                  'native_P_q_FWER_unchanged': True, 'new_scans': 0, 'new_statistical_tests': 0,
                                  'profile_ranking': 'minimum native window-family q, paired raw P from that same test, model ID; no model-level FDR',
                                  'main_motif_limit': 10,
                                  'main_motif_ranking': 'only recorded tests with native q <= 0.05; one test per distinct model, ordered by that q, paired P, model ID; all models remain in outputs',
                                  'significance_rule': 'native selected-test FDR <= 0.05; not a whole-model test',
                                  'RNA_example_rule': 'for each ranked model, earliest stable changing event with complete local forms within display limit; unavailable examples explicitly labelled',
                                  'RNA_example_max_nt': example_max_nt}, indent=2), encoding='utf-8')

    def link(role, path, label, unit):
        return {'role': role, 'path': str(Path(path).resolve()), 'label': label,
                'contrast_id': comparison_id, 'policy': 'stable_significant', 'unit': unit}

    def table(frame, name, unit):
        frame = frame.copy()
        for key, value in [('comparison_id', comparison_id), ('report_policy', 'stable_significant')]:
            if key in frame:
                if not frame[key].dropna().eq(value).all():
                    raise ValueError('Completed motif table has a different recorded scope.')
            else:
                frame.insert(0 if key == 'comparison_id' else 1, key, value)
        path = out / (name + '.tsv')
        frame.to_csv(path, sep='\t', index=False, float_format='%.17g')
        tables.append({'name': name + ' | ' + comparison_id, 'path': str(path), 'chapter': 7, 'unit': unit, 'contrast_id': comparison_id})
        return path

    def figure(name, title, caption, frame, unit, full=(), *, size=(10, 5), links=(), output_only=False):
        fig = Figure(figsize=size); FigureCanvasAgg(fig)
        data = table(frame, name + '-plot-data', unit)
        def save():
            path = out / (_safe(stem + '-' + name) + '.png')
            fig.savefig(path, dpi=145, bbox_inches='tight', facecolor='white')
            fig.savefig(path.with_suffix('.svg'), bbox_inches='tight', facecolor='white')
            panels.append({'id': _safe(stem + '-' + name), 'title': title + ' | ' + comparison_id,
                           'path': str(path), 'chapter': 7, 'contrast_id': comparison_id, 'policy': None,
                           'caption': caption, 'output_only': output_only,
                           'metadata': {'das_scope': scope, 'statistical_family_split': False,
                                        'component': 'motiftrace_sequence',
                                        'complete_manifest_sha256': evidence['manifest_sha256'],
                                        'completion_proof_sha256': evidence['proof_sha256'],
                                        'motif_image_sha256': digest(path)},
                           'data_links': [link('plot_data', data, 'Plot data', unit),
                                          *[link('full_results', path, label, source_unit) for path, label, source_unit in full],
                                          *links, link('methods', methods, 'Methods', 'method_record')]})
            fig.clear()
        return fig, save

    changes = _read(source / 'PAIRED_FORM_INSTANCE_CHANGES.tsv.gz', low_memory=False)
    selected = _read(source / 'FROZEN_SELECTED_EVENTS.tsv', low_memory=False)
    forms = _read(source / 'ACTUAL_RNA_FORMS.tsv.gz', low_memory=False)
    hits = _read(source / 'ACTUAL_LOCUS_HITS.tsv.gz', low_memory=False)
    required = {'event_id', 'event_type', 'model_id', 'instance_state', 'consequence_direction'}
    if not required <= set(changes):
        if not changes.empty:
            raise ValueError('Completed sequence-instance table lacks its declared identity/state columns.')
        # Native zero-hit results use a minimal empty header. Preserve that
        # source; add only empty display columns, never synthetic hit rows.
        changes = changes.reindex(columns=sorted(set(changes) | required))
    if not set(changes.model_id.astype(str)) <= set(models):
        raise ValueError('Sequence results contain a model absent from the completed catalogue.')
    types = [t for t in ['SE', 'A3SS', 'A5SS', 'MXE', 'RI'] if t in set(selected.event_type)]
    types += sorted(set(selected.event_type) - set(types))
    full_instances = [(source / 'PAIRED_FORM_INSTANCE_CHANGES.tsv.gz', 'Complete sequence instances', 'motif_instance')]
    colors = {'inclusion_only': RED, 'exclusion_only': BLUE, 'retained': '#CDD1D6',
              'increased': RED, 'decreased': BLUE, 'unchanged': '#CDD1D6'}
    for field, name, title in [('instance_state', 'instance-states', 'Motif presence in local RNA forms'),
                              ('consequence_direction', 'sequence-opportunity', 'PSI-weighted local sequence opportunity')]:
        counts = changes.groupby(['event_type', field]).size().unstack(fill_value=0).reindex(types, fill_value=0)
        fig, save = figure(name, title, 'Complete recorded instance counts; descriptive local sequence availability, not measured binding or occupancy.',
                           counts.rename_axis('event_type').reset_index(), 'motif_instance_count', full_instances,
                           output_only=field == 'consequence_direction')
        ax = fig.subplots(); bottom = np.zeros(len(counts))
        for state in counts:
            y = counts[state].to_numpy(); ax.bar(counts.index, y, bottom=bottom, color=colors.get(state, '#6993BF'), label=state); bottom += y
        ax.set_ylabel('Exact local RNA motif instances'); ax.legend(frameon=False); save()
    counts = changes.groupby(['model_id', 'instance_state']).size().unstack(fill_value=0).reindex(list(models), fill_value=0)
    catalogue_path = out / 'ALL_MOTIF_MODELS.html'
    catalogue_path.write_text(_catalogue_browser(models), encoding='utf-8')
    tables.append({'name': 'Complete model catalogue and native identities | ' + comparison_id, 'path': str(catalogue_path),
                   'chapter': 7, 'kind': 'html', 'contrast_id': comparison_id})
    fig, save = figure('all-models', 'All declared motif models: sequence instances',
                       'All models are retained in catalogue order, without model-level significance claims. For large catalogues, axis labels are sampled for readability; every model remains in the full browser and count table.',
                       counts.rename_axis('model_id').reset_index(), 'model_instance_count', full_instances,
                       size=(11, min(18, max(5, .28 * len(counts) + 2))),
                       links=[link('full_results', catalogue_path, 'All models browser', 'motif_model')], output_only=True)
    ax = fig.subplots(); left = np.zeros(len(counts))
    for state in counts:
        values = counts[state].to_numpy(); ax.barh(range(len(counts)), values, left=left, color=colors.get(state, '#CDD1D6'), label=state); left += values
    positions = np.arange(len(counts)) if len(counts) <= 60 else np.unique(np.linspace(0, len(counts) - 1, 60).astype(int))
    ax.set_yticks(positions, [str(models[counts.index[i]].get('RBPs', counts.index[i])) for i in positions], fontsize=7)
    ax.invert_yaxis(); ax.legend(frameon=False); ax.set_xlabel('Recorded local RNA motif instances'); save()

    test_source = source / 'NEW_WINDOW_AND_FULL_SEQUENCE_TESTS.tsv.gz'
    aggregate = {token: {'model_id': token, 'RBP': str(model.get('RBPs', token)), 'all_test_rows': 0,
                        'estimable_tests': 0, 'window_q_lt_0_05': 0,
                        'eligible_recorded_windows': 0,
                        'minimum_full_family_window_q': np.nan, 'minimum_raw_window_P': np.nan} for token, model in models.items()}
    selected_tests = {}
    for chunk in _test_chunks(source, native_members):
        _best_recorded_tests(chunk, selected_tests)
        if not set(chunk.model_id.astype(str)) <= set(models):
            raise ValueError('Window test family contains undeclared models.')
        for token, part in chunk.groupby('model_id', sort=False):
            item = aggregate[str(token)]; p, q = part.raw_P_new_paired_swaps, part.q_new_window_family
            item['all_test_rows'] += len(part); item['estimable_tests'] += int(p.notna().sum()); item['window_q_lt_0_05'] += int(q.lt(.05).sum())
            item['eligible_recorded_windows'] += int((part.window_nt.gt(0) & q.notna()).sum())
            for name, values in [('minimum_full_family_window_q', q), ('minimum_raw_window_P', p)]:
                if values.notna().any():
                    item[name] = float(values.min()) if pd.isna(item[name]) else min(item[name], float(values.min()))
    summary = pd.DataFrame(aggregate.values()).sort_values(['minimum_full_family_window_q', 'minimum_raw_window_P', 'model_id'], kind='stable', na_position='last')
    ranking = _motif_ranking(summary, selected_tests, models)
    ranking_path = table(ranking, 'all-motif-ranked-recorded-tests', 'motif_model_recorded_test')
    top10 = significant_motif_ranking(ranking)
    display_summary = motif_display_summary(ranking, comparison_id)
    display_summary.update(ranking_path=str(ranking_path), ranking_sha256=digest(ranking_path))
    if display_summary['status'] != 'significant':
        # A completed negative result gets one concise status in the body.
        # Descriptive instance-count figures remain complete attachments.
        for panel in panels:
            panel['output_only'] = True
    top10_path = table(top10, 'top10-motif-recorded-tests', 'motif_model_recorded_test')
    # The complete browser contains every model and its actual selected native
    # test, including non-estimable models. No min-q value is relabelled as a
    # new model-wide FDR and no P is taken from a different test row.
    catalogue_models = {}
    for row in ranking.astype(object).where(ranking.notna(), None).to_dict('records'):
        catalogue_models[row['model_id']] = {**models[row['model_id']], 'recorded_test_summary': row}
    catalogue_path.write_text(_catalogue_browser(catalogue_models), encoding='utf-8')
    mature = forms.loc[forms.role.eq('foreground') & forms.layer.eq('mature_local_form')]
    compact_forms = {}
    for event, choices in mature.groupby('event_id', sort=True):
        choices = choices.sort_values('region', kind='stable')
        complete_forms = choices.sequence_RNA.notna() & choices.length.gt(0) & choices.length.le(example_max_nt)
        complete_forms &= choices.sequence_RNA.astype(str).str.len().eq(choices.length)
        if len(choices) >= 2 and choices.region.nunique() >= 2 and complete_forms.all():
            compact_forms[event] = choices
    example_count = 0
    for record in top10.to_dict('records'):
        token = record['model_id']
        rank = record['display_rank']
        name = 'top10-motif-' + str(rank).zfill(2) + '-' + _safe(token)
        model_changes = changes.loc[changes.model_id.eq(token)]
        changing_events = sorted(set(model_changes.loc[model_changes.instance_state.ne('retained'), 'event_id']) & compact_forms.keys())
        all_events = sorted(set(model_changes.event_id) & compact_forms.keys())
        event = next(iter(changing_events or all_events), None)
        chosen_forms = compact_forms.get(event, mature.iloc[:0])
        chosen_hits = hits.loc[hits.event_id.eq(event) & hits.model_id.eq(token)].copy()
        hit_path = table(chosen_hits, name + '-exact-sites', 'motif_site')
        form_path = table(chosen_forms, name + '-RNA-forms', 'RNA_form')
        matrix = _matrix(models[token])
        pwm = pd.DataFrame(matrix, columns=list('ACGU')) if matrix is not None else pd.DataFrame(columns=list('ACGU'))
        pwm.insert(0, 'position_1based', range(1, len(pwm) + 1)); pwm.insert(0, 'model_id', token)
        pwm_path = table(pwm, name + '-PWM-matrix', 'PWM_position')
        pvalue, qvalue = record['raw_P_new_paired_swaps'], record['q_new_window_family']
        stats = ('P = not estimable; FDR = not estimable' if pd.isna(qvalue) else
                 f'P = {pvalue:.4g}; test-family FDR = {qvalue:.4g}; ' +
                 ('significant (FDR <= 0.05)' if qvalue <= .05 else 'not significant (FDR > 0.05)'))
        test_label = ('No estimable native test is available for this recorded motif.' if pd.isna(qvalue) else
                      f"{record['event_type']} / {record['direction']} / {record['region']} / {record['method']}; "
                      f"{record['test_scope']}, start {record['window_start0']}, width {record['window_nt']} nt")
        caption = (f'Top {rank} / {len(top10)} displayed models. {stats}. '
                   'Ranking uses the lowest recorded test-family FDR per model, its paired P, then model ID; '
                   'the displayed FDR belongs to that window/full-sequence test, not a whole-model hypothesis. '
                   'Logo/consensus describes sequence preference, not measured binding. The local RNA example is descriptive, '
                   'and does not itself carry the aggregate test P/FDR. All models and complete tests are in the outputs.')
        fig, save = figure(name, f'Top motif {rank}: {record["RBP"]} | {token}', caption,
                           pd.DataFrame([record]), 'motif_model_recorded_test',
                           [(ranking_path, 'All ranked models and native test statistics', 'motif_model_recorded_test'),
                            (test_source, 'Complete native test family', 'motif_window_test'),
                            (source / 'ACTUAL_RNA_FORMS.tsv.gz', 'All actual RNA forms', 'RNA_form'),
                            (source / 'ACTUAL_LOCUS_HITS.tsv.gz', 'All native motif sites', 'motif_site')],
                           size=(12, max(6, 3 + sum(np.ceil(chosen_forms.length / 80)) * .25)),
                           links=[link('plot_data', hit_path, 'Displayed motif positions', 'motif_site'),
                                  link('plot_data', form_path, 'Displayed RNA forms', 'RNA_form'),
                                  link('plot_data', pwm_path, 'Target motif sequence / PWM', 'PWM_position'),
                                  link('full_results', top10_path, 'Top 10 statistics', 'motif_model_recorded_test')])
        axes = fig.subplots(2, 1, gridspec_kw={'height_ratios': [1, 2]})
        _logo(axes[0], models[token])
        axes[0].set_title(f'{record["RBP"]} | {record["motif_sequence"]}\n{stats}', fontsize=10)
        axes[1].axis('off')
        lines = [test_label, 'Exact native test scope; no model-level FDR.', '']
        if event is None:
            lines += ['No eligible complete local RNA example within the display length limit.',
                      'The recorded target motif is shown above; all RNA forms and hits remain in outputs.']
            unavailable.append({'chapter': 7, 'title': f'Compact RNA example | {token} | {comparison_id}',
                                'reason': 'No eligible complete local RNA forms for this ranked model within the display length limit; its recorded motif and available statistical evidence are retained.'})
        else:
            example_count += 1
            lines.append('Illustrative local RNA event: ' + str(event))
            for row in chosen_forms.itertuples():
                lines.append(str(row.region) + ' | ' + str(row.length) + ' nt')
                sequence = str(row.sequence_RNA)
                lines.extend(f'{start + 1:4d}  {sequence[start:start + 80]}' for start in range(0, len(sequence), 80))
                exact = chosen_hits.loc[chosen_hits.form.eq(row.region)]
                lines.append('RNA half-open sites: ' + (', '.join(f'{int(r.RNA_start0)}-{int(r.RNA_end0)}' for r in exact.itertuples()) or 'none recorded'))
        # Native event IDs and dense hit lists can be much longer than an RNA
        # line. Wrap the display only; complete identities/sites stay in TSVs.
        import textwrap
        display_lines = [part for line in lines for part in
                         (textwrap.wrap(line, width=100, break_on_hyphens=False) or [''])]
        text_height = max(2, len(display_lines) * .13)
        fig.set_size_inches(12, max(6, text_height + 2.2))
        axes[0].get_gridspec().set_height_ratios([1.2, text_height])
        axes[1].text(0, 1, '\n'.join(display_lines), va='top', transform=axes[1].transAxes,
                     fontfamily='DejaVu Sans Mono', fontsize=7)
        fig.tight_layout(); save()
        panels[-1]['metadata'].update(motif_view_kind='top_motif', model_id=token, display_rank=int(rank),
                                     main_display_fdr_threshold=.05, main_display_eligible=True)
    fig, save = figure('window-statistics-all-models', 'All-model completed window statistics',
                       'Minimum native window-family q/P defines descriptive catalogue order, not model-level FDR. All declared models and complete original test families are retained; no P/q is recalculated.',
                       summary, 'motif_model_summary', [(test_source, 'Complete native test family', 'motif_window_test')],
                       size=(11, min(18, max(5, .28 * len(summary) + 2))), output_only=True)
    ax = fig.subplots(); shown = summary.iloc[::-1]
    ax.barh(range(len(shown)), shown.estimable_tests, color=BLUE); ax.barh(range(len(shown)), shown.window_q_lt_0_05, color=RED)
    positions = np.arange(len(shown)) if len(shown) <= 60 else np.unique(np.linspace(0, len(shown) - 1, 60).astype(int))
    ax.set_yticks(positions, shown.iloc[positions].model_id, fontsize=7); ax.set_xlabel('Native window/full-sequence tests'); save()

    ranked = summary.loc[summary.minimum_full_family_window_q.notna() & summary.eligible_recorded_windows.gt(0)].head(2).model_id.tolist()
    best = {}
    keys = ['q_new_window_family', 'raw_P_new_paired_swaps', 'event_type', 'direction', 'region', 'window_nt', 'window_start0']
    for chunk in _test_chunks(source, native_members):
        valid = chunk.loc[chunk.model_id.isin(ranked) & chunk.window_nt.gt(0) & chunk.q_new_window_family.notna()]
        for token, part in valid.groupby('model_id', sort=False):
            candidate = part.sort_values(keys, kind='stable').iloc[0]
            key = tuple(candidate[c] for c in keys)
            if token not in best or key < best[token][0]:
                best[token] = key, candidate.to_dict()
    parts = {token: [] for token in best}
    for chunk in _test_chunks(source, native_members):
        for token, (_, chosen) in best.items():
            keep = chunk.model_id.eq(token)
            for key in ('event_type', 'direction', 'region', 'method', 'window_nt'):
                keep &= chunk[key].eq(chosen[key])
            if keep.any():
                parts[token].append(chunk.loc[keep].copy())
    for token, chunks in parts.items():
        profile = pd.concat(chunks, ignore_index=True).sort_values('window_start0', kind='stable')
        profile['plotted_minus_log10_q'] = -np.log10(np.maximum(profile.q_new_window_family, np.nextafter(0., 1.)))
        fig, save = figure('window-profile-' + _safe(token), 'Complete recorded window profile | ' + token,
                           'One complete native position grid, selected descriptively by minimum full-family q/P and stable coordinates. Position gaps/NA remain; native P/q are unmodified and no model-level test is created.',
                           profile, 'motif_window_test', [(test_source, 'Complete native test family', 'motif_window_test')], size=(11, 8), output_only=True)
        axes = fig.subplots(3, 1, sharex=True)
        for ax, column, color in zip(axes, ['mean_paired_score_difference', 'paired_rank_Z', 'plotted_minus_log10_q'], [RED, BLUE, BLUE]):
            ax.plot(profile.window_start0, profile[column], color=color); ax.set_ylabel(column, fontsize=8)
        axes[-1].set_xlabel('Native zero-based window start'); save()
    if len(parts) < 2:
        unavailable.append({'chapter': 7, 'title': 'Additional completed window profiles',
                            'reason': 'Fewer than two declared models have eligible recorded finite-q window grids; missing profiles are not fabricated.'})
    null_path = source / 'MAX_CLUSTER_NULL.npy'
    null = np.load(null_path, allow_pickle=False)
    if null.ndim != 1 or not len(null) or not np.isfinite(null).all():
        raise ValueError('A finite complete one-dimensional maximum-cluster null array is required.')
    counts, edges = np.histogram(null, bins=min(35, len(null)))
    histogram = pd.DataFrame({'bin_left': edges[:-1], 'bin_right': edges[1:], 'count': counts})
    fig, save = figure('whole-family-cluster-null', 'Whole-family maximum-cluster null',
                       f'All {len(null)} native paired-swap maximum-cluster values from one complete family. Original cluster FWER and common null are retained, without per-model recalculation.',
                       histogram, 'null_histogram_bin', [(null_path, 'All native null values', 'permutation'),
                                                          (source / 'NEW_CLUSTER_FWER.tsv', 'Complete cluster FWER results', 'motif_cluster')], output_only=True)
    ax = fig.subplots(); ax.bar(edges[:-1], counts, width=np.diff(edges), align='edge', color=BLUE)
    ax.set(xlabel='Maximum cluster mass across the complete comparison family', ylabel='Paired swaps'); save()
    # Keep the established sequence/statistical attachments first, then retain
    # every additional member authenticated by the scientific completion
    # manifest. Supporting RNA sequences, matching losses and replicate
    # opportunity tables must travel with the completed results too. Never
    # collect an unlisted sidecar merely because it exists in the directory.
    native_paths = [source / name for name in (*SEQUENCE_FILES, *STAT_FILES)]
    native_paths.extend(Path(row['path']) for row in evidence['members'])
    native_paths.append(source / 'COMPLETE.json')
    for path in dict.fromkeys(path.resolve() for path in native_paths):
        name = path.relative_to(source).as_posix()
        tables.append({'name': 'Complete native MotifTrace output | ' + name + ' | ' + comparison_id,
                       'path': str(path), 'kind': 'file', 'chapter': 7, 'contrast_id': comparison_id})
    tables.extend([{'name': 'Scientific completion proof | ' + comparison_id, 'path': str(Path(completion_proof).resolve()),
                    'kind': 'file', 'chapter': 7, 'contrast_id': comparison_id},
                   {'name': 'MotifTrace plotting methods | ' + comparison_id, 'path': str(methods),
                    'kind': 'file', 'chapter': 7, 'contrast_id': comparison_id}])
    for row in evidence['members']:
        path = Path(row['path'])
        if (path.stat().st_size, path.stat().st_mtime_ns, digest(path)) != (row['bytes'], row['mtime_ns'], row['sha256']):
            raise RuntimeError('Original completed MotifTrace family changed during rendering.')
    audit = {'schema': VERSION, 'comparison_id': comparison_id, 'models': len(models),
             'selected_events': len(selected), 'statistical_test_rows': int(summary.all_test_rows.sum()),
             'null_values': len(null), 'figures': len(panels), 'all_declared_models_retained': True,
             'top_motif_cards': len(top10), 'top_motif_RNA_examples': example_count,
             'significant_models': display_summary['significant_models'],
             'main_motif_status': display_summary['status'],
             'main_motif_limit': 10, 'model_level_FDR_created': False,
             'native_P_q_FWER_unchanged': True, 'statistical_family_split': False,
             'new_scans': 0, 'new_statistical_tests': 0, 'originals_unchanged': True,
             'complete_manifest_sha256': evidence['manifest_sha256'], 'proof_sha256': evidence['proof_sha256'],
             'proof_kind': proof_kind, 'native_job_id': native_job_id,
             'observed_at': datetime.now(timezone.utc).isoformat()}
    (out / 'motif-view.receipt.json').write_text(json.dumps(audit, indent=2), encoding='utf-8')
    return {'project_id': project_id, 'source_project_manifest_sha256': source_project_manifest_sha256,
            'status': 'complete', 'panels': panels, 'tables': tables, 'motif_summaries': [display_summary],
            'unavailable': unavailable, 'audit': audit}
