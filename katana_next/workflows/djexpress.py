from pathlib import Path
import pandas as pd
from katana_next.workflows.base import WorkflowBase,WorkflowDefinition,ParameterDefinition,ArtifactWriter
from katana_next.science.djexpress import run_djexpress,sources

class DJExpressWorkflow(WorkflowBase):
    version='1.0.0'
    definition=WorkflowDefinition('programs.djexpress','programs','DJExpress: author co-splicing engine',
        'Which measured junctions vary together, with optional gene-expression-aware second pass?',
        'Pinned original DJExpress JCNA1pass; optional JCNAgenePrepare and JCNA2pass. Katana handles preparation and result registration.',
        'Original modules, sample scores, trait associations, RDS, preparation audit and method provenance','splice junction',
        expected_artifacts=7,parameters=(
            ParameterDefinition('counts','Measured junction counts TSV','text','','Feature-by-sample integer counts, chr:start:end:strand IDs.'),
            ParameterDefinition('traits','Observed traits TSV','text','','Sample-by-trait numeric table, >=2 observed nonconstant traits. Blank uses the declared condition indicators.'),
            ParameterDefinition('gene_expression','Normalized gene expression TSV (second pass)','text','','Gene-symbol-by-sample measured normalized expression. Supply together with the matching GTF.'),
            ParameterDefinition('gtf','Matching GTF (second pass)','text','','Must match the junction assembly and gene expression annotation.'),
            ParameterDefinition('rscript','Rscript executable','text','','R runtime with WGCNA, limma, dplyr, magrittr, tidyr and jsonlite.'),
            ParameterDefinition('input_kind','Count evidence identity','select','junction_counts','rMATS IJC/SJC form evidence is a core-function adaptation, not a raw-junction reproduction.',options=('junction_counts','rmats_form_evidence'))),
        limitations='Uses original author core functions and defaults. The fixed author JCNAprepare raw entry has an undefined min.expressed; this workflow supplies prepared objects explicitly. No new Katana network score or superiority claim.')
    def input_paths(self,context,parameters):
        root,expected=sources()
        if not parameters.get('counts'):raise ValueError('Select a measured junction count table for DJExpress')
        paths=[Path(parameters[k]).resolve() for k in ['counts','traits','gene_expression','gtf'] if parameters.get(k)]
        for path in paths:
            if not path.is_file():raise FileNotFoundError(path)
        return paths+[root/name for name in expected['files']]+[Path(__file__).resolve().parents[1]/'science/djexpress_driver.R',Path(__file__).resolve().parents[1]/'science/djexpress.py']
    def run(self,context,job,run_fingerprint,parameters):
        counts=pd.read_csv(parameters['counts'],sep='\t',index_col=0)
        samples=context.manifest.samples
        if set(counts.columns)!={s.sample_id for s in samples}:raise ValueError('Count columns must match declared project sample IDs exactly')
        if len({s.biological_sample_id for s in samples})<15:raise ValueError('Require at least 15 declared biological samples; cells/repeated observations are not independent samples')
        if parameters.get('traits'):traits=pd.read_csv(parameters['traits'],sep='\t',index_col=0)
        else:traits=pd.get_dummies(pd.Series({s.sample_id:s.condition for s in samples}),dtype=int)
        genes=pd.read_csv(parameters['gene_expression'],sep='\t',index_col=0) if parameters.get('gene_expression') else None
        writer=ArtifactWriter(context,self.definition.workflow_id,run_fingerprint,parameters,self.version,common_metadata={'engine':'DJExpress','algorithm_source':'MauerLab/DJExpress','novel_algorithm':False,'primary_stage':'JCNA1pass'})
        job.report(.1,'Running original DJExpress author core functions')
        out=run_djexpress(counts,traits,writer.directory/'author_results',rscript=parameters.get('rscript') or None,gene_expression=genes,gtf=parameters.get('gtf') or None,input_kind=parameters.get('input_kind','junction_counts'))
        assignment=pd.read_csv(out/'JCNA1pass_assignments.tsv',sep='\t',dtype={'module':str})
        writer.table('consensus_modules__assignments','consensus_modules__assignments.parquet',assignment,entity_type='splice_event',metadata={'consensus_refits':0,'assignment_source':'unmodified original JCNA1pass; name retained for Katana module browser compatibility'})
        for name,path in [('network_observed__scores','JCNA1pass_scores.tsv'),('djexpress__preparation','EVERY_FEATURE_PREPARATION.tsv')]:
            frame=pd.read_csv(out/path,sep='\t');writer.table(name,name+'.parquet',frame,entity_type='module')
        for p in sorted(out.glob('JCNA2pass_*_assignments.tsv')):
            writer.table(p.stem,p.stem+'.parquet',pd.read_csv(p,sep='\t',dtype={'module':str}),entity_type='splice_event')
        for name in ['ENGINE_PROVENANCE.json','ENGINE_COMPLETE.json','METHOD_STATUS.json','R_SESSION.txt','R_STDOUT.log','R_STDERR.log']:
            writer.existing_file('djexpress__'+name, out/name,entity_type='module')
        for p in sorted(out.glob('*.rds')):writer.existing_file('djexpress__'+p.stem,p,entity_type='module')
        for p in sorted(out.glob('*_module_trait_*.tsv')):writer.table(p.stem,p.stem+'.parquet',pd.read_csv(p,sep='\t',index_col=0).reset_index(),entity_type='module')
        job.report(.95,'Registering original first-pass and separate second-pass results')
        return writer.finish()
