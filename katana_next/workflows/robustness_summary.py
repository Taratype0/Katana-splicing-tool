"""Readable cached replicate reliability without silently rerunning models."""
import pandas as pd
from katana_next.workflows.base import WorkflowBase, WorkflowDefinition, ArtifactWriter
from katana_next.workflows.sources import resolve_table

class RobustnessSummaryWorkflow(WorkflowBase):
    version = '1.0.0+english.20261021'
    definition = WorkflowDefinition('validate.replicate_robustness', 'validate', 'Replicate robustness evidence', 'Which event effects have consistent direction and leave-one-out support?', 'Expose existing Stable-AS per-event replicate, read and leave-one-out diagnostics without changing its thresholds or rerunning significance.', 'All event diagnostics, comparison summaries and English explanation; no automatic sample deletion.', 'exact event within comparison', supports_contrast_selection=True, limitations='Cached event-effect leave-one-out diagnostics are not DEG refits or module-discovery resampling. Small samples limit sensitivity. A stable event is not automatically significant; diagnostic flags do not justify deleting a sample.')

    def parent_artifacts(self, context, parameters):
        source = resolve_table(context, 'quality.stable_as', 'Stable-AS event QC', parameters)
        if source is None:
            raise ValueError('Run Stable-AS quality explicitly first')
        return [source]

    def input_paths(self, context, parameters):
        return [a.path for a in self.parent_artifacts(context, parameters)]

    def run(self, context, job, run_fingerprint, parameters):
        source = self.parent_artifacts(context, parameters)[0]
        frame = pd.read_parquet(source.path)
        if parameters.get('contrast_ids'):
            frame = frame[frame.comparison_id.isin(parameters['contrast_ids'])]
        columns = [c for c in frame if c in ('row_id', 'comparison_id', 'gene_id', 'gene_symbol', 'event_id', 'event_type', 'dpsi', 'fdr', 'stable_event', 'stable_significant', 'measurement_eligible', 'first_exclusion_stage') or any((t in c for t in ('loo', 'replicate', 'consistency', 'reads', 'measurement_se')))]
        evidence = frame[columns].copy()
        summary = frame.groupby('comparison_id').agg(events=('row_id', 'size'), measurement_eligible=('measurement_eligible', 'sum'), repeat_supported=('stable_event', 'sum'), stable_significant=('stable_significant', 'sum')).reset_index()
        summary['row_id'] = summary.comparison_id
        writer = ArtifactWriter(context, self.definition.workflow_id, run_fingerprint, parameters, self.version, common_metadata={'limitations': self.definition.limitations, 'inherited_parameters': source.parameters})
        writer.table('All replicate diagnostics', 'replicate_diagnostics.parquet', evidence, entity_type='splice_event', parents=(source.artifact_id,))
        writer.table('Robustness denominators', 'robustness_summary.parquet', summary, entity_type='qc_summary', parents=(source.artifact_id,))
        message = self.definition.limitations
        p = writer.directory / f"REPORT_{'EN'}.html"
        p.write_text('<meta charset="utf-8"><h1>Replicate robustness</h1><p>' + message + '</p>' + summary.to_html(index=False) + '<a href="replicate_diagnostics.tsv">All event evidence</a>', encoding='utf-8')
        writer.existing_file(p.name, p, entity_type='report', parents=(source.artifact_id,))
        return writer.finish()
