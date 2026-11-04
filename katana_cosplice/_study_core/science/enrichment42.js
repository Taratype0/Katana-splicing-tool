// Presentation uses the complete, jointly corrected enrichment result.
function renderEnrichment39(){
 const info=D.result.pathway_enrichment;
 $('pathwayContext').textContent='Each module uses ALL its unique annotated genes, not only the displayed Top 10. Background: '+info.background_genes.toLocaleString()+' unique genes from all '+info.network_events.toLocaleString()+' network events (including unassigned events), intersected with each database. One BH correction covers all '+info.tests.toLocaleString()+' tests across all '+info.tested_modules+' modules and all libraries. GO terms overlap; their counts do not represent independent mechanisms.';
 const host=$('pathwayTable');host.replaceChildren();
 if(info.status==='not_estimable'){
  $('pathwayContext').textContent='Module enrichment was not estimable: '+(info.reasons||[]).join(', ')+'. No P values or FDR were computed. This is not a nonsignificant result. Full gene mapping and background records are provided with the output tables.';
  host.innerHTML='<p class="empty">No eligible enrichment tests. Network modules remain available.</p>';return;
 }
 for(const r of significant){
  const section=document.createElement('article');section.className='enrichment39';section.dataset.module=r.module;
  const collections=info.supplied_libraries||[...new Set(info.libraries.map(a=>a.library))].map(library=>({library,categories:[library]}));
  const categories=[...new Set(collections.flatMap(a=>a.categories))];
  section.innerHTML='<h3>'+E(r.module)+' — functional enrichment</h3><div class="controls"><label>Annotation library <select aria-label="'+r.module+' enrichment library">'+categories.map(c=>'<option>'+c+'</option>').join('')+'</select></label></div><div class="ora42-content"></div>';
  host.append(section);const select=section.querySelector('select'),body=section.querySelector('.ora42-content');let page=0;
  function draw(){
   const paths=D.pathways.filter(p=>p.module===r.module&&p.category===select.value&&p.q_value<=.05).sort((a,b)=>a.q_value-b.q_value||a.pathway.localeCompare(b.pathway));
   const library=collections.find(a=>a.categories.includes(select.value))?.library||select.value,audit=info.libraries.find(a=>'M'+a.module===r.module&&a.library===library);
   if(!audit){body.textContent='No enrichment assessment recorded for this module and library.';return}
   if(audit.status!=='tested'){body.textContent='Enrichment not estimable: '+audit.status+'. No P value or FDR was computed for this module and library.';return}
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
