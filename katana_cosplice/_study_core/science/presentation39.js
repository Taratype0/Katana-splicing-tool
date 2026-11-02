// Presentation only: frozen memberships, tests and enrichment are never altered.
const geneRank39=(a,b)=>(b.max_abs_kME??-1)-(a.max_abs_kME??-1)||a.gene_id.localeCompare(b.gene_id);
function renderIndependentLists(){
 const host=$('independentModuleLists');host.replaceChildren();
 for(const r of significant){
  const genes=D.genes.filter(g=>g.module===r.module).sort(geneRank39),files=D.module_files[r.module],id='independentGenes_'+r.module;
  const card=document.createElement('article');card.className='module-gene-card';card.dataset.module=r.module;
  card.innerHTML='<h3>'+E(r.module)+' — Top 10 genes</h3><p class="thin">'+genes.length+' genes · '+moduleMap[r.module].events+' events · FDR '+E(fmt(r.q_value))+'</p><p class="thin">Ranked by maximum absolute event–module correlation (|kME|); ties by gene name. Membership strength, not regulatory importance.</p><div id="'+id+'"></div><div class="link-row"><a download href="report_tables/'+E(files.genes)+'">Download all '+genes.length+' genes ('+E(r.module)+')</a><a download href="report_tables/'+E(files.events)+'">Download all events ('+E(r.module)+')</a></div>';
  host.append(card);renderTable(id,genes.slice(0,10),[{key:'gene_id',label:'Gene'},{key:'events',label:'Events'},{key:'max_abs_kME',label:'Max |kME|'}],{onClick:g=>{showModule(r.module);chooseGene(g.gene_id);$('moduleTitle').scrollIntoView({behavior:'smooth'})}});
 }
 renderEnrichment39();
}
function heatmap(id,rows,keys,labels,values,limit,click,rawPsi=false){
 const host=$(id);host.innerHTML='';if(!rows.length)return;
 const controls=document.createElement('div');controls.className='controls';
 controls.innerHTML='<label>Heatmap columns <select aria-label="'+E(id)+' column view"><option value="condition">Time / condition — equal donor means</option><option value="sample">Individual samples — scroll horizontally</option></select></label>';
 const panel=document.createElement('div');panel.className='heat39-scroll';host.append(controls,panel);
 function draw(){
  const individual=controls.querySelector('select').value==='sample';
  const groups=individual?D.samples.map((s,i)=>({label:'S'+(i+1),indices:[i],title:s.sample_id+' · '+s.condition+' · donor '+s.biological_sample_id})):conditions.map(c=>({label:c.replace(/^differentiation day /,'Day '),indices:D.samples.map((s,i)=>s.condition===c?i:-1).filter(i=>i>=0),title:c}));
  panel.innerHTML='<table class="heat heat39 '+(individual?'individual':'')+'"><thead><tr><th>'+(rawPsi?'Gene / event':'Module')+'</th>'+groups.map(g=>'<th title="'+E(g.title)+'">'+E(g.label)+'</th>').join('')+'</tr></thead><tbody>'+rows.map((r,ri)=>{
   const vs=values(r);return '<tr><td data-row="'+ri+'" title="'+E(keys(r))+'">'+E(labels(r))+'</td>'+groups.map(g=>{
    const units=new Map();for(const i of g.indices)if(num(vs[i])){const key=D.samples[i].biological_sample_id??D.samples[i].sample_id;if(!units.has(key))units.set(key,[]);units.get(key).push(vs[i])}
    const means=[...units.values()].map(a=>a.reduce((s,v)=>s+v,0)/a.length),v=means.length?means.reduce((s,v)=>s+v,0)/means.length:null;
    const count=g.indices.filter(i=>num(vs[i])).length;
    return '<td><span class="cell" style="background:'+colour(v,limit)+'" title="'+E(labels(r)+' · '+g.title+' · '+(rawPsi?'PSI ':'score ')+fmt(rawPsi&&num(v)?v+.5:v)+' · '+count+' observed samples / '+units.size+' donors')+'"></span></td>';
   }).join('')+'</tr>';
  }).join('')+'</tbody></table>';
  panel.querySelectorAll('[data-row]').forEach(e=>e.onclick=()=>click?.(rows[+e.dataset.row]));
 }
 controls.querySelector('select').onchange=draw;draw();
 const note=document.createElement('p');note.className='fig-notes';note.innerHTML='<span class="heat39-legend"></span> '+(rawPsi?'PSI: blue 0 · white 0.5 · red 1':'Module score (study SD): blue −2 · white 0 · red +2')+' · grey = missing. Default columns average within each donor, then across donors; no missing values are imputed. Display aggregation only: inference still uses individual samples with donor-aware statistics. Hover for values and sample counts.';host.append(note);
}
// Presentation uses the complete, jointly corrected enrichment result.
// Presentation uses the complete, jointly corrected enrichment result.
function renderEnrichment39(){
 const info=D.result.pathway_enrichment;if(!info?.libraries){$('pathwayContext').textContent='No species-matched gene-set resource supplied; enrichment not computed.';$('pathwayTable').innerHTML='';return;}
 $('pathwayContext').textContent='Each module uses ALL its unique annotated genes, not only the displayed Top 10. Background: '+info.background_genes.toLocaleString()+' unique genes from all '+info.network_events.toLocaleString()+' network events (including unassigned events), intersected with each database. One BH correction covers all '+info.tests.toLocaleString()+' tests across all '+info.tested_modules+' modules and all libraries. GO terms overlap; their counts do not represent independent mechanisms.';
 const host=$('pathwayTable');host.replaceChildren();
 for(const r of significant){
  const section=document.createElement('article');section.className='enrichment39';section.dataset.module=r.module;
  const categories=info.libraries.some(a=>a.library==='GO')?['GO BP','GO MF','GO CC','KEGG','Reactome']:[...new Set(info.libraries.map(a=>a.library))];
  section.innerHTML='<h3>'+E(r.module)+' — functional enrichment</h3><div class="controls"><label>Annotation library <select aria-label="'+r.module+' enrichment library">'+categories.map(c=>'<option>'+c+'</option>').join('')+'</select></label></div><div class="ora42-content"></div>';
  host.append(section);const select=section.querySelector('select'),body=section.querySelector('.ora42-content');let page=0;
  function draw(){
   const paths=D.pathways.filter(p=>p.module===r.module&&p.category===select.value&&p.q_value<=.05).sort((a,b)=>a.q_value-b.q_value||a.pathway.localeCompare(b.pathway));
   const library=select.value.startsWith('GO')?'GO':select.value,audit=info.libraries.find(a=>'M'+a.module===r.module&&a.library===library);
   body.innerHTML='<p class="thin">'+audit.module_genes+' / '+audit.all_module_genes+' module genes annotated; '+audit.background_genes+' / '+audit.network_genes+' background genes annotated. '+paths.length+' terms pass joint FDR ≤ 0.05. Ranked by FDR, then term name.</p>';
   if(!paths.length){body.innerHTML+='<p class="empty">No significant terms in this library. No nonsignificant terms substituted.</p>';return}
   const shown=paths.slice(page*12,(page+1)*12);
   body.innerHTML+='<div class="pager"><button data-dir="-1" '+(!page?'disabled':'')+'>Previous terms</button><span>Terms '+(page*12+1)+'–'+Math.min((page+1)*12,paths.length)+' of '+paths.length+'</span><button data-dir="1" '+((page+1)*12>=paths.length?'disabled':'')+'>Next terms</button></div>'+katanaEnrichment40(shown,r.module)+'<p class="fig-notes">X: percentage of annotated module genes in the term. Bubble size: contributing gene count. Colour: joint FDR (red = smaller). Legends use the displayed terms. Gene ratio and FDR are different quantities.</p><details><summary>All significant terms and contributing genes</summary><div id="ora42_'+r.module+'"></div></details>';
   renderTable('ora42_'+r.module,paths,[{key:'pathway',label:'Term'},{key:'q_value',label:'Joint FDR'},{key:'overlap',label:'Genes'},{key:'genes',label:'Contributing genes'}],{pageSize:12});
   body.querySelectorAll('[data-dir]').forEach(b=>b.onclick=()=>{page+=Number(b.dataset.dir);draw()});
  }
  select.onchange=()=>{page=0;draw()};draw();
 }
 if(!significant.length)host.innerHTML='<p class="empty">No qualifying modules in this evidence view.</p>';
}

// Katana approved report reference: katana_next/reports/enrichment_views.py
// Same encodings: gene ratio on X, contributor size, blue/red -log10(FDR).
function katanaEnrichment40(paths,module){
 const width=920,left=365,right=735,top=62,step=52,bottom=top+(paths.length-1)*step+35,height=bottom+155;
 const ratios=paths.map(p=>p.overlap/p.module_genes),scores=paths.map(p=>-Math.log10(Math.max(p.q_value,Number.MIN_VALUE))),maxRatio=Math.max(...ratios)*1.13,minScore=Math.min(...scores),maxScore=Math.max(...scores),maxHits=Math.max(...paths.map(p=>p.overlap));
 const x=v=>left+(right-left)*v/maxRatio,color=v=>{const t=maxScore===minScore?.5:(v-minScore)/(maxScore-minScore);return 'rgb('+[40,99,165].map((v,i)=>Math.round(v+([198,63,74][i]-v)*t)).join(',')+')'};
 const radius=hits=>Math.sqrt((35+350*hits/maxHits)/Math.PI);
 const gradient='enrichmentFdr40_'+module;
 let svg='<svg class="katana-ora40" viewBox="0 0 '+width+' '+height+'" role="img" aria-label="'+E(module)+' Katana pathway enrichment: gene ratio, gene count and FDR"><defs><linearGradient id="'+gradient+'" x1="0" y1="1" x2="0" y2="0"><stop offset="0%" stop-color="#2863A5"/><stop offset="100%" stop-color="#C63F4A"/></linearGradient></defs><text x="550" y="24" text-anchor="middle" font-size="17" fill="#1e2630">'+E(paths[0].category||paths[0].library||'Reactome')+' | '+E(module)+'</text>';
 paths.forEach((p,i)=>{const y=top+i*step,label=p.pathway.replace(/^(REACTOME_|GOBP_|GOMF_|GOCC_)/,'').replaceAll('_',' ').toLowerCase(),lines=[''];for(const word of label.split(' ')){if((lines.at(-1)+' '+word).trim().length>43)lines.push(word);else lines[lines.length-1]=(lines.at(-1)+' '+word).trim()}
  svg+='<line x1="'+left+'" x2="'+right+'" y1="'+y+'" y2="'+y+'" stroke="#e6ebf0"/><text x="'+(left-12)+'" y="'+(y+4-(lines.length-1)*7)+'" text-anchor="end" font-size="12" fill="#1e2630">'+lines.map((s,k)=>'<tspan x="'+(left-12)+'" dy="'+(k?15:0)+'">'+E(s)+'</tspan>').join('')+'</text><circle data-pathway="'+E(p.pathway)+'" cx="'+x(ratios[i])+'" cy="'+y+'" r="'+radius(p.overlap)+'" fill="'+color(scores[i])+'" stroke="white" stroke-width=".5"><title>'+E(p.pathway+' | Gene ratio '+(ratios[i]*100).toFixed(2)+'% ('+p.overlap+'/'+p.module_genes+') | FDR '+fmt(p.q_value)+' | Genes: '+p.genes)+'</title></circle>';
 });
 svg+='<path d="M'+left+' '+(top-25)+'V'+bottom+'H'+right+'" fill="none" stroke="#202020"/>';
 for(let i=0;i<=4;i++){const value=maxRatio*i/4,cx=x(value);svg+='<path d="M'+cx+' '+bottom+'v5" stroke="#202020"/><text x="'+cx+'" y="'+(bottom+22)+'" text-anchor="middle" font-size="12">'+(value*100).toFixed(1)+'%</text>'}
 svg+='<text x="550" y="'+(bottom+49)+'" text-anchor="middle" font-size="13">Module genes in this pathway (%)</text>';
 const barTop=top,barHeight=Math.max(65,(bottom-top)*.75),barX=788;
 svg+='<rect x="'+barX+'" y="'+barTop+'" width="16" height="'+barHeight+'" fill="url(#'+gradient+')" stroke="#202020"/>';
 for(let i=0;i<=4;i++){const y=barTop+barHeight*i/4,q=10**(-(maxScore-(maxScore-minScore)*i/4));svg+='<text x="'+(barX+24)+'" y="'+(y+4)+'" font-size="11">'+q.toPrecision(2)+'</text>'}
 svg+='<text x="788" y="'+(barTop-16)+'" font-size="13">FDR</text><text transform="translate(904 '+(barTop+barHeight/2)+') rotate(-90)" text-anchor="middle" font-size="12">Red = smaller FDR</text>';
 const legendY=bottom+99,counts=[...new Set([Math.min(...paths.map(p=>p.overlap)),Math.round((Math.min(...paths.map(p=>p.overlap))+maxHits)/2),maxHits])];
 svg+='<g class="bubble-size-legend"><text x="'+left+'" y="'+(legendY-22)+'" font-size="13">Bubble size: contributing gene count</text>';
 counts.forEach((count,i)=>{const cx=left+25+i*125;svg+='<circle cx="'+cx+'" cy="'+legendY+'" r="'+radius(count)+'" fill="#6b7280" stroke="white" stroke-width=".5"/><text x="'+(cx+20)+'" y="'+(legendY+4)+'" font-size="12">'+count+' genes</text>'});
 svg+='</g><text x="'+left+'" y="'+(legendY+35)+'" font-size="12">Larger bubble = more genes; redder bubble = stronger enrichment evidence.</text></svg>';
 return svg;
}
