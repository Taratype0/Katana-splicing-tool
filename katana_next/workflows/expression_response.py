"""Explicit factorial expression and ranked pathway workflow."""
from pathlib import Path
import json
import subprocess
import time
import pandas as pd
from katana_next.workflows.base import WorkflowBase, WorkflowDefinition, ParameterDefinition, ArtifactWriter


class ExpressionResponseWorkflow(WorkflowBase):
    version='1.0.1'
    definition = WorkflowDefinition(
        workflow_id='differential.expression_response', workspace_id='differential',
        label='DESeq2 interaction and ranked pathways',
        question='Does gene abundance respond differently to treatment across genotypes?',
        method='Official DESeq2 tximport-aware factorial Wald tests; fgseaMultilevel signed-rank pathway analysis.',
        output_summary='Four simple effects, one interaction, all genes, pathway FDR and leading-edge genes.',
        analysis_unit='gene; independent biological sample',
        parameters=(ParameterDefinition('config_path','Explicit expression configuration JSON','text','','Paths: design, tximport, gene_names, gmt, r_library; optional seed.'),),
        limitations='No causal inference. Independent samples only. Normalized counts cannot replace tximport counts. Cross-cohort confounding cannot be repaired by statistical significance.')

    def _config(self, context, parameters):
        path=Path(parameters.get('config_path',''))
        if not path.is_file(): raise ValueError('Provide expression config JSON')
        config=json.loads(path.read_text(encoding='utf-8'))
        for key in ('design','tximport','gene_names','gmt'):
            if not Path(config[key]).is_file(): raise ValueError('Missing '+key)
        design=pd.read_csv(config['design'],sep='\t')
        required={'sample_id','factor_a','factor_b','biological_sample_id'}
        if not required <= set(design): raise ValueError('Explicit sample and biological identity required')
        if design[list(required)].isna().any().any() or design.sample_id.duplicated().any() or design.biological_sample_id.duplicated().any():
            raise ValueError('Missing or repeated biological/sample identities')
        if not set(design.sample_id)<= {s.sample_id for s in context.manifest.samples}: raise ValueError('Unknown samples')
        for factor in ('factor_a','factor_b'):
            if set(design[factor])!={'reference','target'}: raise ValueError('Use explicit reference/target coding')
        if len(design.groupby(['factor_a','factor_b']))!=4 or design.groupby(['factor_a','factor_b']).size().min()<2: raise ValueError('Incomplete 2x2 design')
        labels=config.get('labels',{})
        if any(not isinstance(labels.get(k),str) or not labels[k].strip() for k in ('reference_a','target_a','reference_b','target_b')):
            raise ValueError('Explicit reference and target labels required for both factors')
        if labels['reference_a']==labels['target_a'] or labels['reference_b']==labels['target_b']:
            raise ValueError('Reference and target labels must differ')
        return config

    def input_paths(self, context, parameters):
        c=self._config(context,parameters)
        return [Path(parameters['config_path']),*[Path(c[k]) for k in ('design','tximport','gene_names','gmt')],Path(__file__).parents[1]/'science/expression_response.R']

    def run(self, context, job, run_fingerprint, parameters):
        c=self._config(context,parameters)
        from katana_next.flow.stages import _executable, _r_environment
        r = _executable(context.manifest.metadata.get('rscript_path'), 'Rscript')
        r_environment = _r_environment(r, c)
        writer=ArtifactWriter(context,self.definition.workflow_id,run_fingerprint,parameters,self.version,
            common_metadata={'method':self.definition.method,'limitations':self.definition.limitations,'design_labels':c.get('labels',{}),'interpretation_limit':c.get('interpretation_limit','')})
        writer.directory.mkdir(parents=True,exist_ok=True)
        cmd=[str(r),str(Path(__file__).parents[1]/'science/expression_response.R'),c['design'],c['tximport'],c['gene_names'],c['gmt'],str(writer.directory),c.get('r_library',''),str(c.get('seed',20260925))]
        if context.manifest.metadata.get('gene_filter'):
            from katana_next.science.gene_filter import expression_filter_rules
            rules, policy = expression_filter_rules(context, writer.directory/'downstream_gene_exclusion_rules.tsv')
            cmd.extend([str(rules), policy])
        with (writer.directory/'stdout.log').open('w') as stdout, (writer.directory/'stderr.log').open('w') as stderr:
            process=subprocess.Popen(cmd,env=r_environment,stdout=stdout,stderr=stderr)
            try:
                while process.poll() is None:
                    job.check_cancelled();time.sleep(1)
            finally:
                if process.poll() is None: process.terminate();process.wait()
        if process.returncode: raise RuntimeError('Expression model failed; see '+str(writer.directory/'stderr.log'))
        from katana_next.science.expression_figures import expression_figures
        for path in expression_figures(writer.directory,c['labels']):
            writer.existing_file(path.name,path,entity_type='plot')
        for path in sorted(writer.directory.glob('*.tsv')):
            if path.name.startswith('unfiltered_model_') or path.name=='downstream_gene_exclusion_rules.tsv':
                continue  # Preserve native model evidence locally; filtered reports use scoped tables only.
            frame=pd.read_csv(path,sep='\t')
            frame['row_id']=[str(i) for i in range(len(frame))]
            writer.table(path.stem,path.stem+'.parquet',frame,entity_type='gene' if 'gene_id' in frame else 'pathway')
        for name in ('METHODS.txt','expression_model.rds','R_sessionInfo.txt'):
            writer.existing_file(name,writer.directory/name,entity_type='provenance')
        return writer.finish()
