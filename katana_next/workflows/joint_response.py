"""Cached, contrast-resolved expression and splicing interpretation."""
from pathlib import Path
import html
import pandas as pd
from katana_next.workflows.base import WorkflowBase, WorkflowDefinition, ArtifactWriter
from katana_next.workflows.bulk import _latest_table

def classify_response(frame):
    result = frame.copy()
    if result.duplicated(['comparison_id', 'gene_id']).any():
        raise ValueError('Expected one gene per comparison')
    both = result.as_available.fillna(False) & result.deg_available.fillna(False)
    result['joint_class'] = 'not jointly assessable'
    for a, d, label in ((True, True, 'DEG and stable DAS'), (True, False, 'stable DAS; not DEG'), (False, True, 'DEG; not stable DAS'), (False, False, 'neither threshold passed')):
        mask = both & result.as_supported.eq(a) & result.deg_supported.eq(d)
        result.loc[mask, 'joint_class'] = label
    summary = result.groupby(['comparison_id', 'joint_class']).size().rename('genes').reset_index()
    totals = result.groupby('comparison_id').size()
    eligible = result.loc[both].groupby('comparison_id').size()
    summary['all_evidence_genes'] = summary.comparison_id.map(totals)
    summary['jointly_assessable_genes'] = summary.comparison_id.map(eligible).fillna(0).astype(int)
    summary['row_id'] = summary.comparison_id + '::' + summary.joint_class
    return (result, summary)

class JointResponseWorkflow(WorkflowBase):
    version = '1.0.1+english.20261021'
    definition = WorkflowDefinition(workflow_id='interpretation.joint_response', workspace_id='prioritize', label='Integrated expression and splicing response', question='Does each contrast affect gene abundance, splicing, or both?', method='Reuse cached gene evidence and its thresholds; no upstream refitting or combined P-value.', output_summary='All-gene table, assessability denominators, PNG/SVG/PDF and English interpretation.', analysis_unit='gene within an explicitly named comparison', supports_contrast_selection=True, limitations='Not significant is not unchanged. Correlated assays are not independent validation. No causal or NMD inference; no factorial refitting in this module.')

    def parent_artifacts(self, context, parameters):
        p = _latest_table(context, 'prioritize.evidence', 'Candidate evidence matrix', parameters)
        return [p] if p else []

    def input_paths(self, context, parameters):
        return [p.path for p in self.parent_artifacts(context, parameters)]

    def run(self, context, job, run_fingerprint, parameters):
        parents = self.parent_artifacts(context, parameters)
        if not parents:
            raise ValueError('Run candidate evidence explicitly first; no implicit upstream execution')
        frame = pd.read_parquet(parents[0].path)
        selected = parameters.get('contrast_ids')
        if selected:
            if set(selected) - set(frame.comparison_id):
                raise ValueError('Selected comparisons missing from cached evidence')
            frame = frame[frame.comparison_id.isin(selected)]
        if frame.empty:
            raise ValueError('No cached gene evidence')
        frame, summary = classify_response(frame)
        # The unchanged classifier reads policy-specific as_supported values;
        # display labels must name that actual selection rather than always RS.
        from katana_next.science.das_selection import EVENT_SET_LABELS, normalize_event_set
        das_label = EVENT_SET_LABELS[normalize_event_set(parameters.get('event_set'))]
        frame['joint_class'] = frame.joint_class.str.replace('stable DAS', das_label, regex=False)
        summary['joint_class'] = summary.joint_class.str.replace('stable DAS', das_label, regex=False)
        summary['row_id'] = summary.comparison_id + '::' + summary.joint_class
        writer = ArtifactWriter(context, self.definition.workflow_id, run_fingerprint, parameters, self.version, common_metadata={'limitations': self.definition.limitations, 'threshold_policy': 'Inherited from parent candidate and DAS artifacts; not recalculated', 'source_artifact': parents[0].artifact_id})
        ids = (parents[0].artifact_id,)
        writer.table('Integrated all-gene response', 'joint_response.parquet', frame, entity_type='gene', parents=ids)
        writer.table('Response classes and denominators', 'joint_summary.parquet', summary, entity_type='gene_summary', parents=ids)
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        labels = {c.contrast_id: c.label for c in context.manifest.contrasts}
        categories = [f'DEG and {das_label}', f'{das_label}; not DEG', f'DEG; not {das_label}']
        matrix = summary.pivot(index='comparison_id', columns='joint_class', values='genes').reindex(columns=categories).fillna(0)
        fig, ax = plt.subplots(figsize=(12, max(4, len(matrix) * 1.3)))
        matrix.rename(index=labels).plot.barh(ax=ax, width=0.8)
        for container in ax.containers:
            ax.bar_label(container, padding=3, fontsize=8)
        ax.set_xlabel('Genes (not events)')
        ax.set_ylabel('')
        ax.set_title('Threshold-passing evidence among jointly assessable genes\nUnavailable and neither-passing genes are listed in the complete table')
        ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), ncol=2, frameon=False)
        ax.spines[['top', 'right']].set_visible(False)
        fig.tight_layout()
        for ext in ('png', 'svg', 'pdf'):
            path = writer.directory / ('joint_response.' + ext)
            fig.savefig(path, dpi=300, bbox_inches='tight')
            writer.existing_file('Joint response ' + ext.upper(), path, entity_type='joint_response_plot', parents=ids, metadata={'mime_type': {'png': 'image/png', 'svg': 'image/svg+xml', 'pdf': 'application/pdf'}[ext]})
        plt.close(fig)
        path = writer.directory / ('REPORT_' + 'EN' + '.html')
        path.write_text('<!doctype html><meta charset="utf-8"><title>' + 'Integrated expression and splicing' + '</title><style>body{font:17px sans-serif;max-width:1100px;margin:40px auto;line-height:1.7}img{width:100%}td,th{padding:8px;border-bottom:1px solid #ddd}</style><h1>' + 'Integrated expression and splicing' + '</h1><p>' + 'Expression measures abundance; splicing measures usage. Each row is a gene within one contrast, not a splice event. Classification requires assessability in both layers. Non-significance does not establish absence. DTU and ordinary DEXSeq columns are complementary, not factorial interaction results. The plot omits the neither-threshold-passed class; the full table retains every gene. These categories do not establish causality or NMD.' + '</p><img src="joint_response.png"><h2>All classes</h2>' + summary.to_html(index=False, escape=True) + '<p>Source: ' + html.escape(parents[0].artifact_id) + '</p>', encoding='utf-8')
        writer.existing_file('Joint interpretation ' + 'EN', path, entity_type='report', parents=ids, metadata={'mime_type': 'text/html'})
        job.report(0.95, 'Joint response exported; no upstream analyses rerun')
        return writer.finish()
