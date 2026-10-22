"""Published-result import entry points for the single-cell workspace."""
from pathlib import Path
from katana_next.science.single_cell_imports import read_config, input_files, run_import
from katana_next.workflows.base import WorkflowBase, WorkflowDefinition, ParameterDefinition, ArtifactWriter

class SingleCellResultImportWorkflow(WorkflowBase):
    version = '1.0.0+english.20261021'

    def __init__(self, mode):
        self.mode = mode
        self.definition = WorkflowDefinition(workflow_id=f'single_cell.import_{mode}', workspace_id='single-cell', label='MARVEL result analysis' if mode == 'marvel' else 'Long-read isoform result analysis', question='What do completed source-tool outputs show across declared cells, donors and conditions?', method='Validated source-output adapter with preserved statistics, design audit, overview and English report.', output_summary='Complete imported tables, overview PNG/SVG, English HTML and input provenance.', analysis_unit='source event / quantified transcript; donor-aware interpretation', parameters=(ParameterDefinition('config_path', 'Import configuration JSON', 'text', '', 'Source paths, method/version, target/reference, metadata and annotation.'),), limitations='No original tool or raw-read pipeline execution. MARVEL cell tests remain source-level evidence. Long-read DIU uses existing Katana donor modeling, not the paper implementation.')

    def _config(self, parameters):
        config = read_config(parameters.get('config_path'))
        if config.get('mode') != self.mode:
            raise ValueError('Import mode disagrees with selected workflow')
        return config

    def input_paths(self, context, parameters):
        config = self._config(parameters)
        required = ('metadata_path', 'results_path') if self.mode == 'marvel' else ('metadata_path', 'counts_path', 'annotation_path')
        for key in required:
            if not config.get(key):
                raise ValueError(f'Import configuration requires {key}')
        paths = input_files(config)
        for path in paths:
            if not path.is_file():
                raise ValueError(f'Import input must be an existing regular file: {path}')
        return [Path(parameters['config_path']), *paths]

    def run(self, context, job, run_fingerprint, parameters):
        config = self._config(parameters)
        if config.get('project_id') != context.manifest.project_id:
            raise ValueError('Declare matching project_id in import configuration')
        writer = ArtifactWriter(context, self.definition.workflow_id, run_fingerprint, parameters, self.version)
        job.report(0.1, 'Validating completed source results; no upstream tool execution')
        tables = run_import(config, writer.directory / 'reader')
        for name, frame in tables.items():
            writer.table(name, name + '.parquet', frame, entity_type='event' if self.mode == 'marvel' else 'transcript_usage')
        for name in ('REPORT_EN.html', 'overview.png', 'overview.svg', 'provenance.json'):
            writer.existing_file(name, writer.directory / 'reader' / name, entity_type='report')
        return writer.finish()
