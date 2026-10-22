"""Explicit immutable hub-event analysis; the dashboard only reads its outputs."""
from pathlib import Path
from dataclasses import asdict
import json
import pandas as pd
from katana_next.workflows.base import ArtifactWriter, ParameterDefinition as P, WorkflowBase, WorkflowDefinition
from katana_next.workflows.sources import validate_source, source_complete
from katana_next.workflows.reliability import read_matrix
from katana_next.science.hub_events import HubSettings, run_hub_analysis


def hub_sources(context, parameters):
    source=context.catalog.get_artifact(parameters.get('source_artifact_id',''))
    validate_source(context,source)
    if source.analysis_id not in {'programs.cosplice','programs.psi_network'} or Path(source.path).stem!='consensus_modules__assignments':
        raise ValueError('Choose completed frozen consensus event-module assignments')
    artifacts=[a for a in context.catalog.list_artifacts(context.manifest.project_id)
               if a.analysis_id==source.analysis_id and a.run_fingerprint==source.run_fingerprint
               and a.metadata.get('project_snapshot_fingerprint')==context.snapshot.manifest_fingerprint and source_complete(context,a)]
    def find(*names):
        a=next((a for a in artifacts if Path(a.path).stem in names),None)
        validate_source(context,a)
        return a
    matrix=find('matrix__psi','stable_qc__matrix')
    annotation=find('matrix__events')
    config=find('resolved_psi_configuration','resolved_cosplice')
    samples=next((a for a in artifacts if Path(a.path).stem=='matrix__samples'),None)
    if samples: validate_source(context,samples)
    return source,matrix,annotation,config,samples


def source_network(config):
    if Path(config.path).suffix=='.json':
        value=json.loads(Path(config.path).read_text(encoding='utf-8'))
    else:
        import yaml
        value=yaml.safe_load(Path(config.path).read_text(encoding='utf-8'))
    network=dict(value.get('network',{}))
    # Resolved NetworkConfig uses soft_power/minimum_network_samples;
    # legacy configuration files use power/min_network_independent_samples.
    if 'power' not in network and 'soft_power' in network:
        network['power']=network['soft_power']
    if 'power' in network and 'soft_power' in network and float(network['power'])!=float(network['soft_power']):
        raise ValueError('Completed discovery configuration has conflicting power values')
    network['min_network_independent_samples']=max(int(network.get('min_network_independent_samples',15)),int(network.get('minimum_network_samples',15)))
    if not {'method','power','signed'} <= set(network):
        raise ValueError('Completed discovery configuration must declare correlation method, power and signedness')
    return network


class HubEventsWorkflow(WorkflowBase):
    version='1.0.1'
    definition=WorkflowDefinition('programs.hub_events','programs','cosplice: validate hub events',
        'Which events remain highly connected, module-representative and stable after biological-block resampling and removal of same-gene edges?',
        'Official R WGCNA correlation and PC1; complete within-module weighted connectivity without export truncation; fixed membership; condition-profile biological-block bootstrap. All-edge and cross-gene ranks are computed separately.',
        'All-event metrics, candidate/rejection reasons, full bootstrap metrics, same-gene sensitivity and exact input/config/session provenance.',
        'splice event within one fixed module and explicitly labelled sample scope',expected_artifacts=8,
        parameters=(
            P('source_artifact_id','Completed module assignments','artifact','','Select one completed immutable discovery run.',options=('programs.cosplice|programs.psi_network','consensus_modules__assignments')),
            P('condition_scope','Sample scopes','select','all_and_conditions','All is the marginal discovery-wide reference; each condition is analysed separately. These are not independent replications.',options=('all_and_conditions','all','conditions')),
            P('bootstrap_iterations','Biological-block bootstrap iterations','number',100,'Preserve complete donor blocks and condition profiles. Invalid replicates count against selection frequency.',minimum=20,maximum=1000,step=10),
            P('seed','Bootstrap seed','number',20261018,'Recorded before resampling; no validation-set tuning.',minimum=0),
            P('top_fraction','Top connectivity fraction','number',.10,'Average tied ranks; separately applied to complete and no-same-gene connectivity.',minimum=.01,maximum=1),
            P('minimum_kme','Minimum module membership (kME)','number',.70,'Signed networks require positive kME; unsigned networks use absolute kME. This is a descriptive eligibility rule, not a P/FDR.',minimum=0,maximum=1),
            P('minimum_selection_frequency','Minimum bootstrap selection frequency','number',.70,'Denominator is all requested bootstraps, including invalid replicates.',minimum=.01,maximum=1),
            P('minimum_valid_bootstrap_fraction','Minimum valid bootstrap fraction','number',.80,'Require estimable kME/connectivity and sufficient pair coverage.',minimum=.01,maximum=1),
            P('minimum_pair_coverage','Minimum estimable pair fraction','number',.80,'Unknown or insufficiently observed pairs are not imputed as known zero connections.',minimum=.01,maximum=1),
        ),limitations='Candidate hubs are descriptive, not causal regulators or statistically significant edges. Observed PSI is reassessed on fixed modules, with source correlation settings. Same-gene removal leaves membership and module PC1 unchanged. Missing gene annotation cannot establish a cross-gene hub. No discovery or trait-based selection is rerun.')

    def parent_artifacts(self,context,parameters):
        return [a for a in hub_sources(context,parameters) if a is not None]

    def input_paths(self,context,parameters):
        return [Path(a.path) for a in self.parent_artifacts(context,parameters)]

    def run(self,context,job,run_fingerprint,parameters):
        source,matrix,annotation,config,sample_source=hub_sources(context,parameters)
        network=source_network(config)
        network_min=int(network.get('min_network_independent_samples',15))
        settings=HubSettings(method=str(network['method']),power=float(network['power']),signed=bool(network['signed']),
            edge_threshold=float(network.get('edge_threshold',.10)),minimum_samples=max(15,network_min),
            **{key:parameters[key] for key in ('bootstrap_iterations','seed','top_fraction','minimum_kme','minimum_selection_frequency','minimum_valid_bootstrap_fraction','minimum_pair_coverage')})
        if sample_source:
            samples=pd.read_parquet(sample_source.path).drop(columns='row_id',errors='ignore')
        else:
            samples=pd.DataFrame([dict(sample_id=s.sample_id,biological_sample_id=s.biological_sample_id,condition=s.condition) for s in context.manifest.samples])
        psi=read_matrix(matrix)
        if set(samples.sample_id)!=set(psi.columns):
            # Cached samples may include unused samples, but no PSI sample may be guessed.
            samples=samples.loc[samples.sample_id.isin(psi.columns)].copy()
        rscript=context.manifest.metadata.get('rscript_path') or network.get('wgcna_rscript')
        if not rscript or not Path(rscript).is_file():
            release=Path(__file__).resolve().parents[3]
            lock=release/'RUNTIME_LOCK.json'
            rscript=json.loads(lock.read_text(encoding='utf-8')).get('rscript') if lock.is_file() else None
        if not rscript or not Path(rscript).is_file(): raise ValueError('Configure an existing Rscript with official WGCNA')
        writer=ArtifactWriter(context,self.definition.workflow_id,run_fingerprint,parameters,self.version,
            common_metadata={'source_discovery_run':source.run_fingerprint,'source_assignment_id':source.artifact_id,
                             'hub_settings':asdict(settings),'method':self.definition.method,'assessment_space':'observed PSI; fixed module membership'})
        metrics,boot,plan=run_hub_analysis(psi,pd.read_parquet(source.path),samples,pd.read_parquet(annotation.path),settings,
            writer.directory,Path(rscript),scope=parameters['condition_scope'],progress=lambda msg:job.report(.3,msg))
        parent_ids=tuple(a.artifact_id for a in (source,matrix,annotation,config,sample_source) if a is not None)
        summary=metrics.groupby(['module','condition'],as_index=False).agg(module_events=('event_key','size'),
            candidates_all_edges=('candidate_all_edges','sum'),candidates_no_same_gene=('candidate_no_same_gene','sum'),
            robust_hub_candidates=('hub_candidate','sum'),biological_units=('biological_units','first'))
        edges=pd.read_csv(writer.directory/'hub_display_edges.tsv',sep='\t')
        if len(edges):
            pairs=[tuple(sorted((str(a),str(b)))) for a,b in zip(edges.event_a,edges.event_b)]
            edges[['event_a','event_b']]=pd.DataFrame(pairs,index=edges.index)
            edges=edges.drop_duplicates(['module','condition','event_a','event_b'])
        for name,frame,entity in [('hub_events__metrics',metrics,'module_hub'),('hub_events__bootstrap',boot,'hub_bootstrap'),('hub_events__summary',summary,'module'),('hub_events__display_edges',edges,'module')]:
            frame=frame.copy(); frame.insert(0,'row_id',[f'{name}:{i}' for i in range(len(frame))])
            writer.table(name,name+'.parquet',frame,entity_type=entity,parents=parent_ids)
        for name in ('hub_plan.json','hub_input_psi.tsv','hub_input_events.tsv','hub_biological_units.tsv','technical_repeat_audit.json','R_session.txt','R_stdout.log','R_stderr.log','hub_raw_metrics.tsv','hub_display_edges.tsv'):
            writer.existing_file(name,writer.directory/name,entity_type='hub_provenance',parents=parent_ids)
        job.report(.95,'Registering full-module hub metrics and same-gene sensitivity; discovery membership unchanged')
        return writer.finish()
