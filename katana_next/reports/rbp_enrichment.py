"""Compact, comparison-specific RBP enrichment from recorded result tables."""
import hashlib
import math
import re
from pathlib import Path

from .sources import read_table


def build_rbp_summaries(artifacts):
    records, seen = [], set()
    for artifact in artifacts:
        path = Path(artifact.get('path', ''))
        if artifact.get('analysis_id') != 'prioritize.rbp_support' or path.name != 'rbp_enrichment.parquet' or not path.is_file():
            continue
        if artifact.get('status', 'complete') not in ('complete', 'completed', 'available'):
            continue
        identity = (str(path.resolve()), artifact.get('run_fingerprint', ''))
        if identity in seen:
            continue
        seen.add(identity)
        table = read_table(path)
        required = {'term', 'genes', 'pvalue', 'fdr'}
        if not required <= set(table.columns):
            raise ValueError('RBP enrichment requires native genes, P and FDR; do not infer contributors from motif names')
        if 'comparison_id' not in table:
            contrasts = artifact.get('contrast_ids', [])
            if len(contrasts) != 1:
                raise ValueError('RBP enrichment table requires an exact comparison')
            table['comparison_id'] = contrasts[0]
        groups = ['comparison_id'] + [c for c in ('cell_type', 'assay', 'background_scope_id') if c in table]
        for group, frame in table.groupby(groups, sort=False, dropna=False):
            values = group if isinstance(group, tuple) else (group,)
            scope = dict(zip(groups, values))
            rows = []
            for row in frame.to_dict('records'):
                p, q = float(row['pvalue']), float(row['fdr'])
                if math.isfinite(q) and not 0 <= q <= 1 or math.isfinite(p) and not 0 <= p <= 1:
                    raise ValueError('Invalid native RBP enrichment P/FDR')
                if not math.isfinite(q) or q > .05:
                    continue
                if not math.isfinite(p):
                    raise ValueError('Significant RBP enrichment has no native P value')
                genes = [x.strip() for x in re.split(r';|,\s*', str(row['genes'])) if x.strip()]
                if not genes or str(row['genes']).lower() in ('nan', 'none') or len(set(genes)) != len(genes):
                    raise ValueError('Significant RBP enrichment must retain its actual unique hit genes')
                if 'overlap_genes' in row and int(row['overlap_genes']) != len(genes):
                    raise ValueError('Native RBP hit count disagrees with its recorded gene list')
                rows.append({'term': str(row['term']), 'hit_genes': genes, 'hits': len(genes),
                             'pvalue': p, 'fdr': q,
                             'fold_enrichment': row.get('fold_enrichment'),
                             'selected_genes': row.get('selected_genes'), 'background_genes': row.get('background_genes'),
                             'fdr_family': str(row.get('fdr_family', 'Recorded native correction family'))})
            rows.sort(key=lambda r: (r['fdr'], r['pvalue'], r['term']))
            token = hashlib.sha256((str(identity) + str(scope)).encode()).hexdigest()[:16]
            records.append({'id': 'rbp-enrichment-' + token, 'chapter': 7,
                            'contrast_id': str(scope['comparison_id']), 'cell_type': str(scope.get('cell_type', 'all')),
                            'assay': str(scope.get('assay', '')), 'policy': None,
                            'source_run': artifact.get('run_fingerprint', ''), 'source_artifact_id': artifact.get('artifact_id'),
                            'full_path': str(path), 'rows': rows, 'tested_sets': len(frame),
                            'foreground_source': artifact.get('metadata', {}).get('foreground_source', 'recorded foreground'),
                            'background_source': artifact.get('metadata', {}).get('background_definition', 'recorded detected/testable background'),
                            'resource_semantics': artifact.get('metadata', {}).get('resource_semantics', 'RBP membership and declared target sets retain their separate resource definitions.'),
                            'no_statistical_recomputation': True})
    return records


def build_motif_hit_summaries(artifacts, output=None):
    """Join saved matched pairs to saved FIMO hits; never scan or retest."""
    registered = {str(Path(a['path']).resolve()): a for a in artifacts if a.get('path')}
    records = []
    for artifact in registered.values():
        path = Path(artifact['path'])
        if artifact.get('analysis_id') not in {'prioritize.event_motifs', 'prioritize.rbp_mechanism'} or path.name != 'motif_enrichment.parquet':
            continue
        inputs = {}
        for name in ('fimo_sites.parquet', 'motif_matching.parquet', 'motif_catalog.parquet'):
            sibling = path.with_name(name).resolve()
            declared = registered.get(str(sibling))
            if not declared or declared.get('run_fingerprint') != artifact.get('run_fingerprint'):
                raise ValueError('Motif hit genes require matched registered native hits, catalog and pairs')
            inputs[name] = read_table(sibling)
        tests = read_table(path)
        comparison = tests.comparison_id.dropna().astype(str).unique()
        if len(comparison) != 1:
            raise ValueError('Motif enrichment must belong to one comparison')
        pairs = inputs['motif_matching.parquet']
        pairs = pairs.loc[pairs.status.eq('matched')]
        if not pairs.empty and pairs.comparison_id.dropna().astype(str).ne(comparison[0]).any():
            raise ValueError('Motif matching belongs to a different comparison')
        sites, catalog = inputs['fimo_sites.parquet'], inputs['motif_catalog.parquet']
        rows, all_rows = [], []
        for row in tests.to_dict('records'):
            q, p = float(row['enrichment_fdr']), float(row['enrichment_pvalue'])
            if math.isfinite(q) and not 0 <= q <= 1 or math.isfinite(p) and not 0 <= p <= 1:
                raise ValueError('Invalid native motif enrichment P/FDR')
            group = catalog.loc[catalog.RBP.str.casefold().eq(str(row['RBP']).casefold()) & catalog.evidence_class.eq(row['evidence_class'])]
            hit_ids = set(sites.loc[sites.motif_token.isin(group.motif_token), 'event_id'])
            hit_pairs = pairs.loc[pairs.foreground_event.isin(hit_ids)]
            genes = sorted(set(hit_pairs.foreground_gene.astype(str)))
            if len(genes) != int(row['foreground_hit_genes']):
                raise ValueError('Motif hit-gene join does not reproduce its native hit count')
            events = sorted(set(hit_pairs.foreground_event.astype(str)))
            all_rows.append({**row, 'foreground_hit_gene_ids': ';'.join(genes), 'foreground_hit_event_ids': ';'.join(events)})
            if not math.isfinite(q) or q > .05:
                continue
            if not math.isfinite(p):
                raise ValueError('Significant motif enrichment has no native P value')
            rows.append({'term': str(row['RBP']) + ' / ' + str(row['evidence_class']), 'hit_genes': genes,
                         'hit_events': events, 'hits':len(genes),
                         'pvalue':p, 'fdr':q, 'fold_enrichment':float(row['hit_fraction_difference']),
                         'fdr_family':str(row.get('FDR_family', 'Native RBP-by-motif-evidence family'))})
        rows.sort(key=lambda row:(row['fdr'], row['pvalue'], row['term']))
        records.append({'id':'motif-hits-' + hashlib.sha256(str(path).encode()).hexdigest()[:16],
                        'chapter':7, 'title':'Boundary motif enrichment', 'effect_label':'Hit fraction difference',
                        'contrast_id':comparison[0], 'cell_type':'all', 'policy':None,
                        'source_run':artifact.get('run_fingerprint',''), 'source_artifact_id':artifact.get('artifact_id'),
                        'full_path':str(path), 'rows':rows, 'tested_sets':len(tests),
                        'foreground_source':'fixed recorded RS-DAS matched foreground genes',
                        'background_source':'recorded sequence-matched measurable control genes',
                        'resource_semantics':'Motif matches are sequence evidence, not experimental RBP binding.',
                        'no_statistical_recomputation':True})
        if output is not None:
            import pandas as pd
            target = Path(output) / 'outputs/enrichment-hits' / (records[-1]['id'] + '.tsv')
            target.parent.mkdir(parents=True, exist_ok=True)
            pd.DataFrame(all_rows).to_csv(target, sep='\t', index=False, float_format='%.17g')
            records[-1]['all_hits_path'] = str(target.resolve())
    return records
