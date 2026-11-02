"""Separate complete gene/event exports for every module; no inference changes."""
from pathlib import Path
import json,re,zipfile
import pandas as pd

def module_exports(exports,genes,events,associations):
    exports=Path(exports);files={};all_names=[];qualified=[]
    ids=sorted(set(genes.get('module',[]))|set(events.get('module',[])),key=lambda m:int(m[1:]))
    selected=set(associations.loc[associations.qualified_candidate,'module']) if 'qualified_candidate' in associations else set()
    for module in ids:
        if re.fullmatch(r'M[1-9][0-9]*',module) is None:raise ValueError('Invalid module export ID')
        g=genes[genes.module.eq(module)].copy();e=events[events.module.eq(module)].copy()
        g=g.sort_values(['gene_id'],kind='stable');e=e.sort_values(['gene_id','event_key'],kind='stable')
        names=dict(genes=f'module_{module}_genes.tsv',events=f'module_{module}_events.tsv')
        g.to_csv(exports/names['genes'],sep='\t',index=False);e.to_csv(exports/names['events'],sep='\t',index=False)
        files[module]=dict(**names,gene_labels=len(g),events_count=len(e));all_names.extend(names.values())
        if module in selected:qualified.extend(names.values())
    for name,members in [('all_module_lists.zip',all_names),('qualified_module_lists.zip',qualified)]:
        with zipfile.ZipFile(exports/name,'w',zipfile.ZIP_DEFLATED) as z:
            for entry in members:z.write(exports/entry,entry)
            z.writestr('README.txt','One genes file and one events file per module. Shared latent programs do not merge memberships. Qualification depends on the current interpretation; see module_qualification.tsv. A gene may occur in multiple modules through different events.\n')
            z.write(exports/'module_qualification.tsv','module_qualification.tsv')
    return files,all_names+['all_module_lists.zip','qualified_module_lists.zip']

def template():
    root=Path(__file__).parent;s=(root/'study_report28.html').read_text(encoding='utf8')
    s=s.replace('CoSplice Next 2.8','CoSplice Next 2.9')
    s=s.replace("programs.size+' latent programs from '","programs.size+' latent program'+(programs.size===1?'':'s')+' from '")
    needle='<div class="controls"><label for="moduleSelect">'
    section='''<p class="chapter-note">Every module in the selected evidence view has its own complete gene list below. Two modules belonging to one latent program still have separate memberships. A shared gene is retained in each relevant module.</p><div class="link-row"><a download href="report_tables/qualified_module_lists.zip">Qualified module lists for this interpretation (separate files in ZIP)</a></div><div id="independentModuleLists"></div><h3 class="subhead">Explore one module: activity, network and individual events</h3>'''
    s=s.replace(needle,section+needle)
    s=s.replace("searchModules();$('moduleHeatmapFigure')","searchModules();renderIndependentLists();$('moduleHeatmapFigure')")
    js='''
function renderIndependentLists(){
 const host=$('independentModuleLists');host.replaceChildren();
 if(!significant.length){host.innerHTML='<p class="empty">No module passes the selected evidence view.</p>';return}
 for(const r of significant){
  const m=moduleMap[r.module],files=D.module_files?.[r.module],genes=D.genes.filter(g=>g.module===r.module).sort((a,b)=>a.gene_id.localeCompare(b.gene_id));
  const card=document.createElement('article');card.className='module-gene-card';card.dataset.module=r.module;
  const id='independentGenes_'+r.module;
  card.innerHTML='<div class="exhibit-head"><span class="exhibit-no">'+E(r.module)+'</span><h3>'+E(r.module)+' — '+genes.filter(g=>g.gene_id!=='Unannotated').length+' genes / '+m.events+' events</h3><span class="scope-tag">'+E(m.latent_factor||'')+' '+E(m.loading_direction||'')+' · FDR '+E(fmt(r.q_value))+'</span></div><div class="link-row">'+(files?'<a download href="report_tables/'+E(files.genes)+'">'+E(r.module)+' genes (TSV)</a><a download href="report_tables/'+E(files.events)+'">'+E(r.module)+' events (TSV)</a>':'')+'</div><div class="controls"><input aria-label="Search '+E(r.module)+' genes" placeholder="Search '+E(r.module)+' genes"><span class="thin" data-count></span></div><div id="'+E(id)+'" class="table-scroll"></div>';
  host.append(card);const input=card.querySelector('input'),count=card.querySelector('[data-count]');
  const draw=()=>{let query=input.value.trim().toLowerCase(),rows=genes.filter(g=>!query||g.gene_id.toLowerCase().includes(query));count.textContent=rows.length+' / '+genes.length+' gene labels';renderTable(id,rows,[{key:'gene_id',label:'Gene'},{key:'events',label:'Events in '+r.module},{key:'max_abs_kME',label:'Max |kME|'}],{pageSize:15,onClick:g=>{showModule(r.module);chooseGene(g.gene_id);$('moduleTitle').scrollIntoView({behavior:'smooth'})}})};input.oninput=draw;draw();
 }
}
'''
    s=s.replace('function updateScope(){',js+'\nfunction updateScope(){')
    # Generic study reports also use the independent module lists, without
    # inventing latent factors for graph-derived modules.
    s=s.replace("qualified.length+' stability-qualified signed modules / '+programs.size+' latent program'+(programs.size===1?'':'s')+' from '","qualified.length+' stability-qualified modules'+(programs.size?' / '+programs.size+' latent program'+(programs.size===1?'':'s'):'')+' from '")
    s=s.replace("E(selectedModule)+' · '+E(m.latent_factor)+' ('+E(m.loading_direction)+')'","E(selectedModule)+(m.latent_factor?' · '+E(m.latent_factor)+' ('+E(m.loading_direction)+')':'')")
    s=s.replace("V.program_note||'Modules are signed factor cores.","V.program_note||(D.manifest.discovery_identity.method==='factor_core'?'Modules are signed factor cores.")
    s=s.replace("all eligible events participate in fitting.';","all eligible events participate in fitting.':'Module memberships are discovered across the whole study. Module count does not establish the number of independent biological mechanisms.');")
    start=s.index('<h3 class="subhead">Matched deletion fits')
    end=s.index('<h3 class="subhead">Qualified modules:',start)
    historical=s[start:end]
    section29='''<div id="resolutionEvidence"><h3 class="subhead">Resolution and coverage controlled comparison</h3><p id="resolutionConclusion" class="signal"></p><p id="resolutionScope" class="chapter-note"></p><div id="resolutionTable" class="table-scroll"></div><p id="resolutionLimit" class="fig-notes"></p><figure class="exhibit"><img id="resolutionFigure" style="width:100%;height:auto" alt="Matched recovery, held-out coherence and fixed-coverage resolution curves"><p class="fig-notes">These are diagnostic transformations of method outputs, not unmodified native algorithms. Actual assigned-event deficits and failed fits remain in the complete tables.</p></figure><details><summary>Common-panel MEGENA coverage diagnosis</summary><div id="resolutionPanel" class="table-scroll"></div><p id="resolutionMegen" class="fig-notes"></p></details><details><summary>All development candidates, including unsuccessful refinements</summary><div id="resolutionCandidates" class="table-scroll"></div><p class="fig-notes">Development pilots: three outer refits and one held-out fold. No candidate is qualified from this pilot; unchanged qualification requires at least 20 full-QC refits. Coverage loss, score redundancy and biological coherence remain visible.</p></details></div>'''
    s=s[:start]+section29+'<details id="historicalComparison"><summary>Historical 2.8 comparisons: module resolution was not matched</summary>'+historical+'</details>'+s[end:]
    auditjs='''
const R=D.resolution_validation||{};
$('resolutionEvidence').hidden=!R.primary_rows;
if(R.primary_rows){
 $('resolutionConclusion').textContent=R.conclusion;
 $('resolutionScope').textContent=R.scope;
 $('resolutionLimit').textContent=R.limit;
 $('resolutionMegen').textContent=R.megen_note;
 $('resolutionFigure').src='report_tables/resolution_audit.png';
 const cols=[{key:'label',label:'Method + diagnostic adapter'},{key:'assigned',label:'Assigned events'},{key:'modules',label:'Modules'},{key:'paired_median_member_recovery',label:'Paired member recovery'},{key:'heldout_coherence',label:'Held-out coherence'},{key:'heldout_absolute_separation',label:'Held-out within–between |r|'},{key:'effective_score_rank',label:'Score effective rank'},{key:'enriched_modules',label:'Enriched modules'}];
 renderTable('resolutionTable',R.primary_rows,cols);
 renderTable('resolutionPanel',R.panel_rows,cols);
 renderTable('resolutionCandidates',R.candidate_rows,[{key:'method',label:'Candidate'},{key:'modules',label:'Modules'},{key:'assigned',label:'Assigned events'},{key:'median_recovery',label:'Pilot recovery'},{key:'pilot_stable_modules',label:'Pass recovery metrics only'},{key:'heldout_coherence',label:'Held-out coherence'},{key:'effective_score_rank',label:'Score effective rank'}],{pageSize:25});
}
'''
    s=s.replace("const V=D.scientific_validation||{};",auditjs+"\nconst V=D.scientific_validation||{};")
    s=s.replace('</style>','#independentModuleLists{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:24px}.module-gene-card{border-top:2px solid var(--navy);margin:22px 0 32px;padding-top:14px;min-width:0}.module-gene-card .exhibit-head{margin-top:0;flex-wrap:wrap}.module-gene-card .metadata td:first-child{font-weight:bold;color:var(--blue)}@media(max-width:1050px){#independentModuleLists{grid-template-columns:1fr}}\n</style>')
    (root/'study_report29.html').write_text(s,encoding='utf8')

if __name__=='__main__':template()
