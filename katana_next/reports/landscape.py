"""Descriptive views of saved PSI measurements, never a differential refit."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from .sources import read_table

BLUE, RED, GRAY = '#2863A5', '#C63F4A', '#CDD1D6'


def build_landscape_views(project, artifacts, destination):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap
    from scipy.cluster.hierarchy import linkage, leaves_list
    destination = Path(destination).resolve()
    destination.mkdir(parents=True, exist_ok=True)
    panels, tables, unavailable = [], [], []
    cmap = LinearSegmentedColormap.from_list('katana_b', [BLUE, 'white', RED])
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10, 'axes.spines.top': False,
                         'axes.spines.right': False, 'svg.fonttype': 'none'})

    def export(frame, stem, title, chapter, unit):
        target = destination / (stem + '.tsv')
        frame.to_csv(target, sep='\t', index=False)
        tables.append(dict(name=title, path=str(target), chapter=chapter, unit=unit))
        return target

    def data_link(role, path, label, unit):
        return dict(role=role, path=str(path), label=label, contrast_id='all', policy=None, unit=unit)

    def save(fig, key, title, chapter, caption, *, plot_data, unit, full_results=(), parameters=None, **extra):
        target = destination / (key + '.png')
        fig.savefig(target, dpi=150, bbox_inches='tight', facecolor='white')
        fig.savefig(target.with_suffix('.svg'), bbox_inches='tight', facecolor='white')
        plt.close(fig)
        method = destination / (key + '.methods.json')
        method.write_text(json.dumps({'schema': 'katana.report_plot_methods/v1', 'figure_id': key,
                                     'method': caption, 'parameters': parameters or {},
                                     'contrast_id': 'all', 'policy': None, 'unit': unit,
                                     'new_inferential_fit': False}, indent=2), encoding='utf-8')
        links = [data_link('plot_data', plot_data, 'Plot data', unit),
                 *[data_link('full_results', path, label, source_unit) for path, label, source_unit in full_results],
                 data_link('methods', method, 'Methods', 'method_record')]
        panels.append(dict(id=key, path=str(target), chapter=chapter, title=title,
                           caption=caption, contrast_id='all', data_links=links, **extra))

    def find(name):
        return next((r for r in artifacts if r.get('name') == name and Path(r['path']).is_file() and Path(r['path']).suffix in ('.parquet', '.tsv', '.csv')), None)

    samples = pd.DataFrame([{'sample_id': s['sample_id'],
                             'biological_sample_id': s.get('biological_sample_id', s['sample_id']),
                             'condition': s.get('condition', ''), **s.get('attributes', {})}
                            for s in project.get('samples', [])])
    if not samples.empty:
        counts = samples.groupby('condition', sort=False).biological_sample_id.nunique()
        fig, ax = plt.subplots(figsize=(11, 4))
        ax.bar(counts.index, counts.values, color=BLUE)
        ax.set(ylabel='Independent biological samples', title='Declared biological sample design')
        ax.tick_params(axis='x', rotation=30)
        for x, y in enumerate(counts.values): ax.text(x, y, str(y), ha='center', va='bottom')
        sample_path = export(samples, 'biological-samples', 'Biological samples', 1, 'biological_sample')
        design_path = export(counts.rename('independent_biological_samples').rename_axis('condition').reset_index(),
                             'biological-design-plot-data', 'Exact design plot counts', 1, 'condition_count')
        save(fig, 'biological-design', 'Biological samples and experimental conditions', 1,
             'Sample identities and independence come from the project manifest. Replicate numbers do not establish pairing.',
             plot_data=design_path, unit='condition_count', full_results=[(sample_path, 'Full sample metadata', 'biological_sample')])

    matrix_record = find('Sample-level PSI matrix')
    if not matrix_record:
        unavailable.append(dict(chapter=2, title='PSI matrix, heatmap and consistency',
                                reason='No completed sample-level PSI matrix was supplied.'))
        return dict(panels=panels, tables=tables, unavailable=unavailable)
    frame = read_table(matrix_record['path'])
    columns = [s['sample_id'] for s in project.get('samples', []) if s['sample_id'] in frame]
    if len(columns) < 2:
        return dict(panels=panels, tables=tables, unavailable=[dict(chapter=2, title='PSI views', reason='Fewer than two declared samples have saved PSI measurements.')])
    values = frame[columns].apply(pd.to_numeric, errors='coerce').to_numpy(float)
    if (np.isfinite(values) & ((values < 0) | (values > 1))).any():
        raise ValueError('Saved PSI measurements must be fractions in [0, 1]; no implicit percent conversion or clipping is applied.')
    labels = frame.get('event_id', pd.Series(frame.index.astype(str))).astype(str).tolist()
    finite = np.isfinite(values)
    with np.errstate(invalid='ignore'):
        variance = np.nanvar(values, axis=1)
        means = np.nanmean(values, axis=1)
    order = np.argsort(-np.nan_to_num(variance, nan=-1), kind='stable')
    selected = order[(finite.sum(axis=1)[order] >= 2)][:200]
    centered = values[selected] - means[selected, None]
    filled = np.nan_to_num(centered)
    if len(selected) > 1:
        row_order = leaves_list(linkage(filled, method='average', metric='euclidean'))
        selected = selected[row_order]
    fig, ax = plt.subplots(figsize=(12, 7))
    image = ax.imshow(values[selected], aspect='auto', cmap=cmap, vmin=0, vmax=1, interpolation='nearest')
    ax.set_xticks(range(len(columns)), columns, rotation=90, fontsize=7)
    ax.set(ylabel='Variable measured events (overview)', title='Saved sample-level PSI patterns')
    ax.set_yticks([]); fig.colorbar(image, ax=ax, label='PSI')
    plotted = frame.iloc[selected][[c for c in ('event_id', 'gene_id', 'gene_symbol') if c in frame] + columns].copy()
    plotted.insert(0, 'plot_row', range(len(plotted)))
    plotted_path = export(plotted, 'psi-heatmap-plot-data', 'Exact PSI heatmap rows in plotted order', 2, 'event_by_sample_PSI_fraction')
    save(fig, 'psi-heatmap', 'PSI heatmap and splicing pattern clustering', 2,
         'Overview: up to 200 highest-variance measured events, average-linkage ordering of mean-centered PSI; colors show native PSI. '
         'Every saved event is available in the full matrix browser. This measurement view is independent of differential-significance selection.',
         plot_data=plotted_path, unit='event_by_sample_PSI_fraction',
         full_results=[(matrix_record['path'], 'Full native PSI matrix', 'event_by_sample_PSI_fraction'),
                       (destination / 'PSI_MATRIX.html', 'Interactive full PSI matrix', 'event_by_sample_PSI_fraction')],
         parameters={'maximum_displayed_events': 200, 'displayed_events': len(selected), 'samples_in_order': columns,
                     'ordering': 'average-linkage Euclidean on mean-centered PSI; missing centered values filled only for ordering'},
         full_matrix=True, selection_scope='fixed_measured_matrix')

    # A portable all-row matrix browser, outside the report body. Values remain
    # native JSON numbers and missing measurements remain null, not zero.
    matrix_data = {'samples': columns, 'events': labels,
                   'values': [[float(x) if np.isfinite(x) else None for x in row] for row in values],
                   'conditions': [next((s.get('condition', '') for s in project.get('samples', []) if s['sample_id'] == c), '') for c in columns]}
    template = Path(__file__).with_name('assets').joinpath('matrix-browser.html').read_text('utf-8')
    browser = destination / 'PSI_MATRIX.html'
    browser.write_text(template.replace('__DATA__', json.dumps(matrix_data, separators=(',', ':'), allow_nan=False).replace('<', '\\u003c')), encoding='utf-8')
    tables.append(dict(name='Full PSI matrix browser: every saved event', path=str(browser), kind='html', chapter=2))

    pca_record = find('Sample PSI PCA')
    if pca_record and Path(pca_record['path']).suffix == '.parquet':
        pca = read_table(pca_record['path'])
    else:
        pca = None
    if pca is not None and {'PC1', 'PC2', 'sample_id'}.issubset(pca):
        fig, ax = plt.subplots(figsize=(11, 6))
        conditions = pca.get('condition', pd.Series('sample', index=pca.index)).fillna('Unknown')
        colors = [BLUE, RED, '#6993BF', '#D9777F', '#173C66', '#862B33', '#99B6D4', '#E7ADB2']
        for i, condition in enumerate(conditions.unique()):
            selected_pca = pca.loc[conditions == condition]
            ax.scatter(selected_pca.PC1, selected_pca.PC2, label=condition, color=colors[i % len(colors)], s=45)
            for _, r in selected_pca.iterrows(): ax.annotate(str(r.sample_id), (r.PC1, r.PC2), xytext=(3, 3), textcoords='offset points', fontsize=6)
        exp1 = float(pca.get('explained_PC1', pd.Series([np.nan])).iloc[0])
        exp2 = float(pca.get('explained_PC2', pd.Series([np.nan])).iloc[0])
        ax.set(xlabel=f'PC1 ({exp1:.1%})', ylabel=f'PC2 ({exp2:.1%})', title='Saved PSI principal components')
        ax.legend(loc='upper left', bbox_to_anchor=(1.01, 1), fontsize=8, frameon=False)
        pca_path = export(pca, 'psi-pca-plot-data', 'Exact saved PCA coordinates and identities', 2, 'sample_PCA_coordinates')
        save(fig, 'psi-pca', 'PSI PCA: saved measurements', 2,
             'Completed source PCA coordinates, explained variance and sample identities are retained. No differential model or PCA is refitted. '
             'Its original measured-event population is recorded separately from ordinary DAS/RS-DAS selection.',
             plot_data=pca_path, unit='sample_PCA_coordinates', full_results=[(pca_record['path'], 'Full saved PCA results', 'sample_PCA_coordinates')],
             selection_scope='fixed_completed_pca')

    similarity_record = find('Sample PSI similarity')
    if similarity_record and Path(similarity_record['path']).suffix == '.parquet':
        sim = read_table(similarity_record['path'])
        condition_by_sample = {s['sample_id']: s.get('condition', '') for s in project.get('samples', [])}
        rows = sim.loc[sim.sample_a != sim.sample_b].copy()
        rows['pair_type'] = ['Within condition' if condition_by_sample.get(a) == condition_by_sample.get(b) else 'Between conditions' for a, b in zip(rows.sample_a, rows.sample_b)]
        fig, ax = plt.subplots(figsize=(8, 4))
        groups = [rows.loc[rows.pair_type == label, 'spearman_r'].dropna().values for label in ['Within condition', 'Between conditions']]
        parts = ax.boxplot(groups, patch_artist=True, showfliers=True)
        ax.set_xticks([1, 2], ['Within condition', 'Between conditions'])
        for patch, color in zip(parts['boxes'], [BLUE, GRAY]): patch.set_facecolor(color)
        ax.set(ylabel='Saved pairwise Spearman correlation', title='Replicate consistency (descriptive)')
        pairs_path = export(rows, 'replicate-consistency-plot-data', 'Every plotted within/between sample pair', 2, 'sample_pair_spearman')
        save(fig, 'replicate-consistency', 'Within-condition consistency', 2,
             'All saved off-diagonal sample pairs; both orientations are retained consistently. Pairwise values are descriptive and dependent; no independent-pair significance test is claimed.',
             plot_data=pairs_path, unit='sample_pair_spearman', full_results=[(similarity_record['path'], 'Full saved pairwise results', 'sample_pair_spearman')])
        diag = []
        for sample in columns:
            s = rows.loc[rows.sample_a == sample]
            within = s.loc[s.pair_type == 'Within condition', 'spearman_r']
            diag.append(dict(sample_id=sample, condition=condition_by_sample.get(sample, ''), median_within_condition_similarity=float(within.median()) if len(within) else np.nan,
                             diagnostic_only=True, removed=False))
        diagnostics = pd.DataFrame(diag)
        p = destination / 'sample-consistency-diagnostics.tsv'; diagnostics.to_csv(p, sep='\t', index=False)
        tables.append(dict(name='Sample consistency diagnostics (no sample exclusion)', path=str(p), chapter=2))
        fig, ax = plt.subplots(figsize=(12, 4))
        ax.bar(range(len(diag)), diagnostics.median_within_condition_similarity, color=BLUE)
        ax.set_xticks(range(len(diag)), diagnostics.sample_id, rotation=90, fontsize=7)
        ax.set(ylabel='Median within-condition similarity', title='Sample-level consistency diagnostics')
        save(fig, 'sample-consistency-diagnostic', 'Sample consistency and potential outliers', 2,
             'Continuous within-condition similarity is shown for every sample. No arbitrary outlier cutoff or automatic sample removal is applied.',
             plot_data=p, unit='sample_median_spearman', full_results=[(similarity_record['path'], 'Full saved pairwise results', 'sample_pair_spearman')])
        # The complete similarity image belongs in output files, not the report.
        pivot = sim.pivot(index='sample_a', columns='sample_b', values='spearman_r').reindex(index=columns, columns=columns)
        fig, ax = plt.subplots(figsize=(9, 8)); image = ax.imshow(pivot, cmap=cmap, vmin=0, vmax=1)
        ax.set_xticks(range(len(columns)), columns, rotation=90, fontsize=6); ax.set_yticks(range(len(columns)), columns, fontsize=6)
        fig.colorbar(image, ax=ax, label='Spearman correlation')
        matrix_path = export(pivot.rename_axis('sample_a').reset_index(), 'sample-similarity-matrix-plot-data',
                             'Complete similarity matrix in plotted sample order', 2, 'sample_by_sample_spearman')
        save(fig, 'sample-similarity-complete', 'Sample similarity: complete matrix', 2,
             'Full completed pairwise similarity; output-only as requested.', plot_data=matrix_path, unit='sample_by_sample_spearman',
             full_results=[(similarity_record['path'], 'Full saved pairwise results', 'sample_pair_spearman')], output_only=True)
    return dict(panels=panels, tables=tables, unavailable=unavailable)
