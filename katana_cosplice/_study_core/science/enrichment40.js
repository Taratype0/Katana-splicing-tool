// Katana approved report reference: katana_next/reports/enrichment_views.py
// Same encodings: gene ratio on X, contributor size, blue/red -log10(FDR).
function katanaEnrichment40(paths,module){
 const width=920,left=365,right=735,top=62,step=52,bottom=top+(paths.length-1)*step+35,height=bottom+155;
 const ratios=paths.map(p=>p.overlap/p.module_genes),scores=paths.map(p=>-Math.log10(Math.max(p.q_value,Number.MIN_VALUE))),maxRatio=Math.max(...ratios)*1.13,minScore=Math.min(...scores),maxScore=Math.max(...scores),maxHits=Math.max(...paths.map(p=>p.overlap));
 const x=v=>left+(right-left)*v/maxRatio,color=v=>{const t=maxScore===minScore?.5:(v-minScore)/(maxScore-minScore);return 'rgb('+[40,99,165].map((v,i)=>Math.round(v+([198,63,74][i]-v)*t)).join(',')+')'};
 const radius=hits=>Math.sqrt((35+350*hits/maxHits)/Math.PI);
 const gradient='enrichmentFdr40_'+module;
 let svg='<svg class="katana-ora40" viewBox="0 0 '+width+' '+height+'" role="img" aria-label="'+E(module)+' Katana pathway enrichment: gene ratio, gene count and FDR"><defs><linearGradient id="'+gradient+'" x1="0" y1="1" x2="0" y2="0"><stop offset="0%" stop-color="#2863A5"/><stop offset="100%" stop-color="#C63F4A"/></linearGradient></defs><text x="550" y="24" text-anchor="middle" font-size="17" fill="#1e2630">Reactome | '+E(module)+'</text>';
 paths.forEach((p,i)=>{const y=top+i*step,label=p.pathway.replace(/^REACTOME_/,'').replaceAll('_',' ').toLowerCase(),lines=[''];for(const word of label.split(' ')){if((lines.at(-1)+' '+word).trim().length>43)lines.push(word);else lines[lines.length-1]=(lines.at(-1)+' '+word).trim()}
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
