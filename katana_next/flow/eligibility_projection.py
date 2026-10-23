"""Prepare current read-supported inputs from immutable historical AS results.

This is input preparation, not a new rMATS/DAS fit or a completed-workflow
version upgrade. Native effects, significance, P/FDR and legacy flags remain
available unchanged. New analyses receive an explicitly named projection.
"""
from pathlib import Path
from types import SimpleNamespace

from .model import file_sha256, write_json


def project_existing_as(source, project, outdir, *, contrast_ids):
    import pandas as pd
    from katana_next.infrastructure.manifest_store import ManifestStore
    from katana_next.science.bulk import apply_stable_event_qc, READ_SUPPORT_POLICY_VERSION
    source, project, outdir = Path(source).resolve(), Path(project).resolve(), Path(outdir).resolve()
    _, manifest = ManifestStore().load(project)
    requested = list(contrast_ids)
    contrasts = {c.contrast_id: c for c in manifest.contrasts}
    if not requested or len(set(requested)) != len(requested) or set(requested)-set(contrasts):
        raise ValueError('Choose exact unique declared comparison IDs for historical input preparation')
    if outdir.exists():
        raise FileExistsError('Preserve existing projection; choose a new output directory')
    original = pd.read_parquet(source)
    required = {'row_id','event_id','comparison_id','significant','gene_id','chr','strand',
                'IncLevel1','IncLevel2','IJC_SAMPLE_1','IJC_SAMPLE_2','SJC_SAMPLE_1','SJC_SAMPLE_2'}
    if not required <= set(original):
        raise ValueError('Historical result lacks complete native event/count vectors: '+str(sorted(required-set(original))))
    selected = original.loc[original.comparison_id.isin(requested)].copy()
    if set(selected.comparison_id) != set(requested) or selected.duplicated(['comparison_id','event_id']).any():
        raise ValueError('Historical measurable universe is missing or duplicates a requested comparison/event')
    if not pd.api.types.is_bool_dtype(selected.significant) or selected.significant.isna().any():
        raise ValueError('Native caller significance must be a complete boolean flag; no implicit recaller test')
    if 'reverse_direction' in selected:
        if not pd.api.types.is_bool_dtype(selected.reverse_direction) or selected.reverse_direction.isna().any():
            raise ValueError('Historical count-vector direction must be explicit boolean')
    else:
        raise ValueError('Historical count-vector direction is not declared')
    for cid in requested:
        mask = selected.comparison_id.eq(cid)
        contrast = contrasts[cid]
        for event_type, rows in selected.loc[mask].groupby('event_type'):
            declarations = [a for a in manifest.assets if a.kind == 'rmats_result'
                            and a.metadata.get('contrast_id') == cid
                            and a.metadata.get('event_type') == event_type]
            if not declarations:
                raise ValueError('Count-vector sample order has no declared rMATS source: '+cid+'/'+event_type)
            orders = set()
            for asset in declarations:
                metadata = asset.metadata
                if (tuple(metadata.get('numerator_sample_ids', ())) != contrast.numerator.sample_ids
                    or tuple(metadata.get('denominator_sample_ids', ())) != contrast.denominator.sample_ids
                    or not isinstance(metadata.get('reverse_direction'), bool)):
                    raise ValueError('Declared count-vector sample order differs from the comparison: '+cid)
                orders.add(metadata['reverse_direction'])
            if len(orders) != 1 or set(rows.reverse_direction) != orders:
                raise ValueError('Historical vector direction is ambiguous or differs from its declared source: '+cid)
        for column, label in (('numerator', contrast.numerator.label), ('denominator', contrast.denominator.label)):
            if column in selected and not selected.loc[mask, column].eq(label).all():
                raise ValueError('Historical effect direction differs from the declared comparison: '+cid)
        selected.loc[mask,'numerator_expected_samples'] = len(contrast.numerator.sample_ids)
        selected.loc[mask,'denominator_expected_samples'] = len(contrast.denominator.sample_ids)
    from katana_next.science.gene_filter import apply_gene_filter, gene_filter_details
    filter_details = gene_filter_details(manifest)
    measured_before_exclusions = len(selected)
    eligible_genes = apply_gene_filter(SimpleNamespace(manifest=manifest), selected)
    excluded_genes = selected.loc[~selected.index.isin(eligible_genes.index),
                                  ['row_id','comparison_id','event_id','gene_id']].copy()
    selected = eligible_genes
    revised = apply_stable_event_qc(selected, {'min_reads':10})
    projected = selected.copy()
    # Keep original policies/flags under source_* before publishing new input
    # eligibility. Caller statistics and measurements are never overwritten.
    changed = set(revised)^set(selected) | {c for c in revised if c in selected and not revised[c].equals(selected[c])}
    changed.add('stable_significant')
    for column in changed:
        if column in selected:
            projected['source_'+column] = selected[column]
        if column in revised:
            projected[column] = revised[column]
        elif column in projected:
            projected = projected.drop(columns=column)
    projected['stable_significant'] = selected.significant & revised.measurement_eligible
    immutable = [c for c in selected if c not in changed and c not in {'numerator_expected_samples','denominator_expected_samples'}]
    pd.testing.assert_frame_equal(projected[immutable],selected[immutable],check_exact=True)
    for column in ('pvalue','p_value','PValue','fdr','FDR','dpsi','IncLevelDifference','significant'):
        if column in selected:
            pd.testing.assert_series_equal(projected[column],selected[column],check_exact=True)
    outdir.mkdir(parents=True)
    destination = outdir/'measured_as_input_projection.parquet'
    projected.to_parquet(destination,index=False)
    summary = projected.groupby('comparison_id',sort=False).agg(
        measured_events=('event_id','size'), caller_significant=('significant','sum'),
        read_supported=('measurement_eligible','sum'), rs_das=('stable_significant','sum')).reset_index()
    summary.to_csv(outdir/'eligibility_counts.tsv',sep='\t',index=False)
    excluded_genes.to_csv(outdir/'gene_excluded_input_rows.tsv',sep='\t',index=False)
    receipt = {'schema':'katana.flow.as_input_projection/v1','status':'prepared',
               'analysis_id':'flow.as_input_projection','scientific_model_fitted':False,
               'source_result_modified':False,'native_statistics_reestimated':False,
               'source':{'path':str(source),'sha256':file_sha256(source)},
               'project_manifest':{'path':str(project),'sha256':file_sha256(project)},
               'contrast_ids':requested,'read_support_policy':READ_SUPPORT_POLICY_VERSION,
               'minimum_counts':10,'aggregation':'every declared biological sample in both arms',
               'retired_direction_and_SE_gates_applied':False,
               'count_vector_order_verified_against_declared_samples':True,
               'gene_filter':filter_details,'measured_before_gene_exclusions':measured_before_exclusions,
               'gene_excluded_rows':len(excluded_genes),
               'source_rows':len(original),'projected_rows':len(projected),
               'native_columns_preserved':immutable,'legacy_flags_preserved_as_source_columns':sorted(changed & set(selected)),
               'output':{'path':str(destination),'sha256':file_sha256(destination)}}
    write_json(outdir/'projection.receipt.json',receipt)
    return receipt
