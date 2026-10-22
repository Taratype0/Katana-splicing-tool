"""Explicit immutable hub-event analysis; the dashboard only reads its outputs."""
from pathlib import Path
from dataclasses import asdict
import json
import shutil
import pandas as pd
from katana_next.workflows.base import ArtifactWriter, ParameterDefinition as P, WorkflowBase, WorkflowDefinition
from katana_next.workflows.sources import validate_source, source_complete
from katana_next.workflows.reliability import read_matrix
from katana_next.science.hub_events import HubSettings, network_label_policy, run_hub_analysis
from katana_next.science.hub_events import summarize_hubs
from katana_next.infrastructure.fingerprints import fingerprint_file, fingerprint_value


HUB_SCIENTIFIC_FILES = (
    'hub_plan.json','hub_input_psi.tsv','hub_input_events.tsv','hub_biological_units.tsv',
    'technical_repeat_audit.json','R_session.txt','R_stdout.log','R_stderr.log',
    'hub_raw_metrics.tsv','hub_display_edges.tsv',
)


def _scientific_fingerprint(parents, settings, rscript, scope):
    # Rank/preset labels are a view of the same raw measurements. Their switches
    # and cutoffs must not resubmit R or the biological-block bootstrap.
    r_settings = {key:getattr(settings,key) for key in (
        'method','power','signed','edge_threshold','bootstrap_iterations','seed',
        'minimum_samples','block_size')}
    return fingerprint_value({
        'contract':'fixed-module-hub-measurements-v1',
        'inputs':[{ 'artifact_id':a.artifact_id, 'fingerprint':a.fingerprint,
                    'content':fingerprint_file(a.path,force_full=True).model_dump(mode='json')}
                  for a in parents if a is not None],
        'rscript':fingerprint_file(Path(rscript),force_full=True).model_dump(mode='json'),
        'script':fingerprint_file(Path(__file__).parents[1]/'science'/'hub_metrics.R',force_full=True).model_dump(mode='json'),
        'settings':r_settings,'scope':scope,
    })


def _completed_scientific_source(context, scientific_fingerprint):
    """Only complete, unchanged raw execution bundles can supply a new view."""
    groups = {}
    for artifact in context.catalog.list_artifacts(context.manifest.project_id,analysis_id='programs.hub_events'):
        if artifact.metadata.get('hub_scientific_fingerprint') != scientific_fingerprint:
            continue
        if artifact.path.name not in HUB_SCIENTIFIC_FILES:
            continue
        if str(artifact.status) != 'complete' or artifact.metadata.get('artifact_contract_version') != '3.0-correctness':
            continue
        try:
            validate_source(context,artifact)
            expected=fingerprint_value({'run':artifact.run_fingerprint,'name':artifact.name,
                'path':artifact.path.name,
                'content':fingerprint_file(artifact.path,force_full=True).model_dump(mode='json')})
            if expected != artifact.fingerprint:
                continue
        except (ValueError, FileNotFoundError):
            continue
        groups.setdefault(artifact.run_fingerprint,{})[artifact.path.name]=artifact
    for run, artifacts in groups.items():
        if set(HUB_SCIENTIFIC_FILES) <= set(artifacts):
            return run,artifacts
    return None


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
    version='2.0.0'
    definition=WorkflowDefinition('programs.hub_events','programs','cosplice: hub metrics and optional labels',
        'How do all fixed-module events rank in connectivity and module representation, and how sensitive are these descriptions to resampling and same-gene edges?',
        'Official R WGCNA correlation and PC1; complete within-module weighted connectivity without export truncation; fixed membership; condition-profile biological-block bootstrap. All events and all-edge/cross-gene ranks are retained. Heuristic candidate labels are opt-in and never filter rows.',
        'All-event metrics and complete ranks, preset-selection sensitivity frequencies, optional candidate labels or not-applied status, full bootstrap metrics, same-gene sensitivity and input/config/session provenance.',
        'splice event within one fixed module and explicitly labelled sample scope',expected_artifacts=8,
        parameters=(
            P('source_artifact_id','Completed module assignments','artifact','','Select one completed immutable discovery run.',options=('programs.cosplice|programs.psi_network','consensus_modules__assignments')),
            P('condition_scope','Sample scopes','select','all_and_conditions','All is the marginal discovery-wide reference; each condition is analysed separately. These are not independent replications.',options=('all_and_conditions','all','conditions')),
            P('apply_threshold_labels','Apply optional descriptive candidate labels','boolean',False,'Off retains all continuous metrics/ranks and leaves candidate labels not applied. On annotates the recorded preset without filtering events or changing modules.'),
            P('bootstrap_iterations','Biological-block bootstrap iterations','number',100,'Preserve complete donor blocks and condition profiles. Preset-selection frequency is a threshold-dependent sensitivity metric, not a probability of truth.',minimum=20,maximum=1000,step=10),
            P('seed','Bootstrap seed','number',20261018,'Recorded before resampling; no validation-set tuning.',minimum=0),
            P('top_fraction','Optional-label connectivity fraction','number',.10,'Preset sensitivity and opt-in labels only; average tied ranks, separately for all and cross-gene connectivity. Complete ranks remain available.',minimum=.01,maximum=1),
            P('minimum_kme','Optional-label minimum kME','number',.70,'Preset sensitivity and opt-in labels only. Signed networks use positive kME; unsigned networks use absolute kME. Not P/FDR.',minimum=0,maximum=1),
            P('minimum_selection_frequency','Optional-label preset selection frequency','number',.70,'Opt-in labels only; denominator is all requested bootstraps, including invalid replicates.',minimum=.01,maximum=1),
            P('minimum_valid_bootstrap_fraction','Optional-label valid bootstrap fraction','number',.80,'Opt-in labels only; requires estimable kME/connectivity and the recorded pair-coverage preset.',minimum=.01,maximum=1),
            P('minimum_pair_coverage','Optional-label estimable pair fraction','number',.80,'Preset sensitivity and opt-in labels only. Missing pairs remain unknown; continuous coverage is always reported.',minimum=.01,maximum=1),
        ),limitations='Optional candidate labels and preset selection frequencies are Katana descriptive heuristics, not P/FDR, calibrated probabilities or causal regulation. Labels are not applied by default. Mathematical non-estimability remains explicit. Observed PSI, fixed source membership and full ranks are preserved. Missing host-gene annotation cannot establish cross-gene evidence. No discovery or trait-based selection is rerun.')

    def normalize_parameters(self, supplied):
        supplied=dict(supplied or {})
        enabled=supplied.get('apply_threshold_labels',False)
        if enabled is None:
            enabled=False
        elif isinstance(enabled,str):
            if enabled.strip().lower() not in {'true','false'}:
                raise ValueError('apply_threshold_labels must be true or false')
            enabled=enabled.strip().lower()=='true'
        elif type(enabled) is not bool:
            raise ValueError('apply_threshold_labels must be a boolean')
        supplied['apply_threshold_labels']=enabled
        return super().normalize_parameters(supplied)

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
            apply_threshold_labels=parameters.get('apply_threshold_labels',False),
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
        scientific_fingerprint=_scientific_fingerprint((source,matrix,annotation,config,sample_source),settings,rscript,parameters['condition_scope'])
        reused=_completed_scientific_source(context,scientific_fingerprint)
        scientific_source_run=(reused[1]['hub_raw_metrics.tsv'].metadata.get('hub_scientific_source_run',reused[0])
                               if reused else run_fingerprint)
        writer=ArtifactWriter(context,self.definition.workflow_id,run_fingerprint,parameters,self.version,
            common_metadata={'source_discovery_run':source.run_fingerprint,'source_assignment_id':source.artifact_id,
                             'hub_settings':asdict(settings),'network_label_policy':network_label_policy(settings),
                             'hub_scientific_fingerprint':scientific_fingerprint,'hub_scientific_source_run':scientific_source_run,
                             'hub_raw_source_view_run':reused[0] if reused else run_fingerprint,
                             'scientific_measurements_reused':bool(reused),
                             'method':self.definition.method,'assessment_space':'observed PSI; fixed module membership'})
        if reused:
            job.report(.3,'Reusing a verified complete raw hub execution; only labels and descriptive summaries are refreshed')
            for filename,artifact in reused[1].items():
                shutil.copy2(artifact.path,writer.directory/filename)
            metrics,boot=summarize_hubs(pd.read_csv(writer.directory/'hub_raw_metrics.tsv',sep='\t'),settings)
            # The original measured event/host-gene table is copied unchanged;
            # do not infer new annotation while rebuilding a label-only view.
            measured_events=pd.read_csv(writer.directory/'hub_input_events.tsv',sep='\t')
            metrics=metrics.merge(measured_events[['event_key','gene','host_gene']],on='event_key',how='left')
        else:
            metrics,boot,plan=run_hub_analysis(psi,pd.read_parquet(source.path),samples,pd.read_parquet(annotation.path),settings,
                writer.directory,Path(rscript),scope=parameters['condition_scope'],progress=lambda msg:job.report(.3,msg))
        view_policy=writer.directory/'hub_view_policy.json'
        view_policy.write_text(json.dumps({'network_label_policy':network_label_policy(settings),
            'view_settings':asdict(settings),'scientific_source_run':scientific_source_run,
            'scientific_measurements_reused':bool(reused),
            'source_plan_role':'Original raw execution plan, retained unchanged; view_settings control this result view.'},indent=2),encoding='utf-8')
        parent_ids=tuple(a.artifact_id for a in (source,matrix,annotation,config,sample_source) if a is not None)
        if reused:
            parent_ids=tuple(dict.fromkeys((*parent_ids,*(a.artifact_id for a in reused[1].values()))))
        summary=metrics.groupby(['module','condition'],as_index=False).agg(module_events=('event_key','size'),
            candidates_all_edges=('candidate_all_edges',lambda values:values.sum(min_count=1)),
            candidates_no_same_gene=('candidate_no_same_gene',lambda values:values.sum(min_count=1)),
            robust_hub_candidates=('hub_candidate',lambda values:values.sum(min_count=1)),
            label_status=('label_status',lambda values:';'.join(sorted(set(values)))),
            estimable_events=('estimation_status',lambda values:int(values.eq('estimable').sum())),
            biological_units=('biological_units','first'))
        edges=pd.read_csv(writer.directory/'hub_display_edges.tsv',sep='\t')
        if len(edges):
            pairs=[tuple(sorted((str(a),str(b)))) for a,b in zip(edges.event_a,edges.event_b)]
            edges[['event_a','event_b']]=pd.DataFrame(pairs,index=edges.index)
            edges=edges.drop_duplicates(['module','condition','event_a','event_b'])
        for name,frame,entity in [('hub_events__metrics',metrics,'module_hub'),('hub_events__bootstrap',boot,'hub_bootstrap'),('hub_events__summary',summary,'module'),('hub_events__display_edges',edges,'module')]:
            frame=frame.copy(); frame.insert(0,'row_id',[f'{name}:{i}' for i in range(len(frame))])
            writer.table(name,name+'.parquet',frame,entity_type=entity,parents=parent_ids)
        for name in (*HUB_SCIENTIFIC_FILES,'hub_view_policy.json'):
            writer.existing_file(name,writer.directory/name,entity_type='hub_provenance',parents=parent_ids)
        job.report(.95,'Registering full-module hub metrics and same-gene sensitivity; discovery membership unchanged')
        return writer.finish()
