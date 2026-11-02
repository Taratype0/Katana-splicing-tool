"""Production CoSplice: one sample-level network per complete project study."""
from pathlib import Path
import json,os,sys,subprocess,time
import pandas as pd
from katana_next.workflows.base import WorkflowBase,WorkflowDefinition,ParameterDefinition,ArtifactWriter

PARAMETERS=(ParameterDefinition('paired','Paired observations','boolean',False,'Optional. Use biological_sample_id to link repeated observations; declared paired comparisons retain their design.'),
    ParameterDefinition('association_contrast_ids','Module association comparisons','json',None,'Optional comparison ID list; network discovery always retains all project samples.'),
    ParameterDefinition('min_module_size','Minimum module events','number',20,'WGCNA minimum module size.',minimum=2),
    ParameterDefinition('merge_cut_height','Module merge height','number',.25,'WGCNA eigengene merge cut height.',minimum=.001,maximum=.999),
    ParameterDefinition('event_limit_mode','Maximum variable events','select','5000','Choose 2000, 5000, 10000 or Custom. After QC, retain up to this many events by variance; every study sample remains. Never filter by DAS significance.',options=('2000','5000','10000','Custom')),
    ParameterDefinition('max_events','Event count for Custom','number',5000,'Used only when Maximum variable events is Custom. Enter an integer from 40 to 100000. Fewer eligible events means all eligible events are used.',minimum=40,maximum=100000,step=1),
    ParameterDefinition('bootstrap_iterations','Optional stability refits','number',0,'Off by default: one whole-study WGCNA + kME fit. Opt-in unit-deletion refits assess stability separately; kME is not stability validation.',minimum=0,maximum=500))

class CoSpliceWorkflow(WorkflowBase):
    version='4.3.3'
    definition=WorkflowDefinition('programs.cosplice','programs','CoSplice: whole-study network','Which splicing modules vary across all comparable study samples?',
        'One WGCNA/TOM + kME network using all declared project samples, followed by module association and stability analysis.',
        'Study modules, per-module genes/events, sample scores, continuous kME and study report','individual sample',parameters=PARAMETERS,optional_asset_kinds=('cosplice_config','psi_matrix','event_annotation','rmats_result'),supports_contrast_selection=False,
        limitations='Repeated-donor inference requires matching calibration; otherwise associations remain exploratory. A pooled response is not within-condition rewiring.')

    def normalize_parameters(self,supplied):
        values=dict(supplied or {})
        if isinstance(values.get('association_contrast_ids'),str):
            values['association_contrast_ids']=json.loads(values['association_contrast_ids']) if values['association_contrast_ids'].strip() else None
        # Retain compatibility with numeric API requests made before the presets.
        if 'event_limit_mode' not in values and 'max_events' in values:
            values['event_limit_mode']='Custom'
        mode=values.get('event_limit_mode','5000')
        if mode in {'2000','5000','10000'}:
            values['max_events']=int(mode)
        normalized=super().normalize_parameters(values)
        selected=normalized.get('association_contrast_ids')
        if selected is not None and (not isinstance(selected,list) or any(not isinstance(c,str) for c in selected)):
            raise ValueError('Module association comparisons must be a JSON list of comparison IDs')
        if normalized['max_events'] is None:
            raise ValueError('Custom event limit must be an integer')
        return normalized

    def validate_assets(self,context):
        super().validate_assets(context)
        kinds={a.kind for a in context.manifest.assets}
        if not ('cosplice_config' in kinds or {'psi_matrix','event_annotation'} <= kinds or 'rmats_result' in kinds):
            raise ValueError('CoSplice requires a PSI matrix plus event annotation, native rMATS results, or a CoSplice input configuration')

    def input_paths(self,context,parameters):
        from katana_cosplice.study_design import project_study_design
        project_study_design(context.manifest,paired=parameters['paired'],association_contrast_ids=parameters.get('association_contrast_ids'))
        self.validate_assets(context)
        kinds={a.kind for a in context.manifest.assets}
        if 'rmats_result' in kinds and not kinds.intersection({'psi_matrix','cosplice_config'}):
            from katana_next.flow.network_inputs import resolve_network_rmats_mode
            resolve_network_rmats_mode(context.manifest,parameters)
        paths=super().input_paths(context,parameters)
        from katana_cosplice import _study_core
        paths.extend(sorted((Path(_study_core.__file__).parent/'science').glob('*.py')))
        paths.extend(sorted((Path(_study_core.__file__).parent/'science').glob('*.R')))
        paths.extend(sorted((Path(_study_core.__file__).parent/'science').glob('*.html')))
        paths.extend(sorted((Path(_study_core.__file__).parent/'science').glob('*.css')))
        paths.extend(sorted((Path(_study_core.__file__).parent/'science').glob('*.js')))
        paths.extend([Path(__file__),Path(__file__).parents[2]/'katana_cosplice/study_design.py'])
        paths.extend(a.path for a in context.manifest.assets if a.kind in {'read_support','gene_sets_gmt'})
        asset=next((a for a in context.manifest.assets if a.kind=='cosplice_config'),None)
        if asset:
            from katana_cosplice.config import load_config
            from katana_next.workflows.programs import configuration_input_files
            paths.extend(configuration_input_files(load_config(asset.path)))
        return list(dict.fromkeys(paths))

    def run(self,context,job,run_fingerprint,parameters):
        writer=ArtifactWriter(context,self.definition.workflow_id,run_fingerprint,parameters,self.version,common_metadata={'analysis_scope':'whole_study','comparison_id':'all','displayed_comparisons':[],'cohort_samples':len(context.manifest.samples),'method':'WGCNA + kME','stability_refits_requested':int(parameters['bootstrap_iterations'])})
        target=writer.directory/'whole-study'
        if target.exists():
            import uuid
            target=writer.directory/('whole-study-'+uuid.uuid4().hex[:10])
        target.mkdir(exist_ok=False)
        inputs=target/'source_inputs';inputs.mkdir()
        assets={a.kind:a.path for a in context.manifest.assets}
        samples=pd.DataFrame([{'sample_id':s.sample_id,'biological_sample_id':s.biological_sample_id,'condition':s.condition,**{k:v for k,v in s.attributes.items() if k not in {'sample_id','biological_sample_id','condition'}}} for s in context.manifest.samples])
        samples.to_csv(inputs/'samples.tsv',sep='\t',index=False)
        from katana_cosplice.study_design import project_study_design
        sample_policy, contrasts=project_study_design(context.manifest,paired=parameters['paired'],association_contrast_ids=parameters.get('association_contrast_ids'))
        support=None
        if 'psi_matrix' in assets:
            psi=assets['psi_matrix'];annotation=assets['event_annotation'];support=assets.get('read_support')
        else:
            from katana_cosplice.config import load_config
            from katana_cosplice.pipeline import CoSplicePipeline
            from katana_next.services.design_editor import sync_cosplice_design
            if 'cosplice_config' in assets:
                cfg=load_config(assets['cosplice_config']);cfg.output_dir=target/'native_inputs';cfg.output_dir.mkdir()
                cfg.gene_exclusion=dict(context.manifest.metadata.get('gene_filter') or {})
                declared=context.manifest.model_copy(update={'metadata':{**context.manifest.metadata,'design_revision':'whole-study-43'}})
                sync_cosplice_design(cfg,declared,cfg.output_dir)
                cfg.selected_conditions=list(dict.fromkeys(samples.condition))
                built=CoSplicePipeline(cfg).build_matrix()
            else:
                from katana_next.flow.network_inputs import prepare_project_network,matrix_display_labels
                from katana_next.flow.model import file_sha256
                cfg,transport,_=prepare_project_network(context.manifest_path,target/'native_inputs',
                    {'network_scope':'whole_study','paired':parameters['paired'], **({'rmats_mode':parameters['rmats_mode']} if parameters.get('rmats_mode') else {})},file_sha256)
                with matrix_display_labels(transport['condition_aliases']):
                    built=CoSplicePipeline(cfg).build_matrix()
            psi=built['psi'];annotation=built['events']
            ijc=pd.read_csv(built['ijc'],sep='\t',index_col=0);sjc=pd.read_csv(built['sjc'],sep='\t',index_col=0)
            if not ijc.index.equals(sjc.index) or not ijc.columns.equals(sjc.columns):raise ValueError('Junction-support matrix identities differ')
            support=inputs/'read_support.tsv.gz';(ijc+sjc).to_csv(support,sep='\t',index_label='event_key')
        columns=pd.read_csv(psi,sep='\t',nrows=0).columns[1:]
        if len(columns)!=len(set(columns)) or set(columns)!=set(samples.sample_id):raise ValueError('Whole-study PSI must contain every project sample exactly once; a contrast subset cannot replace it')
        # Apply the project's already-declared gene exclusion policy to imported inputs too.
        if context.manifest.metadata.get('gene_filter'):
            from katana_next.science.gene_filter import apply_gene_filter
            events=pd.read_csv(annotation,sep='\t');kept=apply_gene_filter(context,events)
            if len(kept)!=len(events):
                annotation=inputs/'events.tsv';kept.to_csv(annotation,sep='\t',index=False)
                matrix=pd.read_csv(psi,sep='\t',index_col=0);matrix=matrix.loc[matrix.index.isin(kept.event_key)]
                psi=inputs/'psi.tsv.gz';matrix.to_csv(psi,sep='\t',index_label='event_key')
                if support:
                    table=pd.read_csv(support,sep='\t',index_col=0).reindex(matrix.index);support=inputs/'filtered_support.tsv.gz';table.to_csv(support,sep='\t',index_label='event_key')
        config={'method':'wgcna','feature_cap':int(parameters['max_events']),'bootstrap':int(parameters['bootstrap_iterations']),'qc_profile':'unit_robust','merge_modules':False,'value_kind':'usage','study_id':context.manifest.project_id,'cohort_scope':'All declared comparable project samples; one shared network across every condition','cohort_complete':True,'cohort_design':context.manifest.metadata.get('cosplice_cohort_design','comparable_samples'),'psi':str(psi),'samples':str(inputs/'samples.tsv'),'annotation':str(annotation),'read_support':str(support) if support else None,'gene_sets':str(assets['gene_sets_gmt']) if 'gene_sets_gmt' in assets else None,'expression':None,'condition_a':None,'condition_b':None,'genome_build':context.manifest.annotation.genome_build}
        config.update(paired=sample_policy['paired'],sample_size_policy=sample_policy,contrasts=contrasts,min_module_size=parameters['min_module_size'],merge_cut_height=parameters['merge_cut_height'])
        if 'study_id' in samples and samples.study_id.nunique()==1:config['study_id']=str(samples.study_id.iloc[0])
        (target/'config.json').write_text(json.dumps(config,indent=2),encoding='utf8')
        from katana_cosplice import _study_core
        core=Path(_study_core.__file__).parent
        threads=str(max(1,int(os.environ.get('KATANA_CPUS','2'))))
        env=os.environ.copy();env.update(PYTHONPATH=os.pathsep.join([str(core/'science'),str(core.parents[1])]),PYTHONUTF8='1',PYTHONDONTWRITEBYTECODE='1',OPENBLAS_NUM_THREADS=threads,MKL_NUM_THREADS=threads,OMP_NUM_THREADS=threads,KATANA_CPUS=threads,KATANA_MIN_MODULE_SIZE=str(parameters['min_module_size']),KATANA_MERGE_CUT_HEIGHT=str(parameters['merge_cut_height']))
        job.report(.1,f'Building ONE network across {len(samples)} samples and {samples.condition.nunique()} conditions')
        with (target/'worker.log').open('w',encoding='utf8') as log:
            proc=subprocess.Popen([sys.executable,'-B','-X','utf8',str(core/'science/production43.py'),str(target)],env=env,stdout=log,stderr=subprocess.STDOUT,creationflags=subprocess.CREATE_NO_WINDOW if os.name=='nt' else 0)
            try:
                while proc.poll() is None:
                    job.check_cancelled();time.sleep(.3)
            except BaseException:
                import psutil
                owner=psutil.Process(proc.pid)
                for child in owner.children(recursive=True):
                    try:child.terminate()
                    except psutil.NoSuchProcess:pass
                proc.terminate();proc.wait(timeout=15);raise
        if proc.returncode:raise RuntimeError('Whole-study CoSplice failed; '+(target/'worker.log').read_text(encoding='utf8')[-2500:])
        status=json.loads((target/'status.json').read_text())
        if status.get('status')!='completed':raise ValueError('Whole-study worker did not complete')
        for name,file in [('CoSplice whole-study modules, genes and functions','cosplice_chapter6.html'),('CoSplice complete study report','report.html'),('CoSplice complete outputs','results.zip')]:
            writer.existing_file(name,target/file,entity_type='cosplice_report',metadata={'completion_verified':True,'preserve_relative_report_links':file.endswith('.html'),'display_in_main':file=='cosplice_chapter6.html','report_role':'whole_study_chapter' if file=='cosplice_chapter6.html' else 'complete_outputs','portable_root':str(writer.directory),'contrast_id':'all'})
        job.report(.98,'Completed one whole-study network; Shared network retains every sample and condition')
        return writer.finish()

class PsiNetworkWorkflow(CoSpliceWorkflow):
    definition=WorkflowDefinition('programs.psi_network','programs','CoSplice: whole-study PSI network',CoSpliceWorkflow.definition.question,CoSpliceWorkflow.definition.method,CoSpliceWorkflow.definition.output_summary,'individual sample',parameters=PARAMETERS,required_asset_kinds=('psi_matrix','event_annotation'),supports_contrast_selection=False,limitations=CoSpliceWorkflow.definition.limitations)
