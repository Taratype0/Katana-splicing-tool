"""Portable HTML representation of the exact event-module dashboard payload."""
import json


def dashboard_section(payloads, *, identifier="katana-module-explorer"):
    data = json.dumps(payloads, ensure_ascii=True, allow_nan=False).replace("<", "\\u003c")
    return TEMPLATE.replace("__ID__", identifier).replace("__DATA__", data)


def render_dashboard_png(payload, target):
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap
    d = payload
    fig, axes = plt.subplots(2, 2, figsize=(13, 9), layout="constrained")
    nodes = {n["event_key"]: n for n in d["nodes"]}
    for e in d["edges"]:
        a, b = nodes[e["event_a"]], nodes[e["event_b"]]
        axes[0,0].plot([a["x"],b["x"]], [a["y"],b["y"]], color="#bf9095" if e["same_gene"] else "#9bb6bd", alpha=.45, lw=.5+e["weight"])
    maximum = max([n["connectivity"] for n in nodes.values()] or [1]) or 1
    for i,n in enumerate(nodes.values()):
        axes[0,0].scatter(n["x"], n["y"], s=20+100*n["connectivity"]/maximum, c="#287a98", edgecolor="white")
        if i<6:axes[0,0].annotate(n["gene"] or f"Event {i+1}",(n["x"],n["y"]),xytext=(5,5),textcoords="offset points",fontsize=8)
    axes[0,0].set_title("Event association network")
    axes[0,0].axis("off")
    heat=np.array([[np.nan if x is None else x for x in row] for row in d["heatmap"]["z"]])
    cmap=LinearSegmentedColormap.from_list("soft_blue_red",["#7aabcd","#fafafa","#dd979e"])
    cmap.set_bad("#d8dcde")
    image=axes[0,1].imshow(heat,aspect="auto",cmap=cmap,vmin=-2,vmax=2,interpolation="nearest")
    axes[0,1].set(title="Sample PSI heatmap (row z-score)",xlabel="Samples ordered by condition",ylabel="Splice events",xticks=[],yticks=[])
    fig.colorbar(image,ax=axes[0,1],label="PSI row z-score",shrink=.6)
    for ax,rows,field,title in [(axes[1,0],d["scores"],"value","Module activity by condition"),(axes[1,1],d["event_profile"],"psi","Selected event inclusion (PSI)")]:
        for i,condition in enumerate(d["conditions"]):
            values=[r[field] for r in rows if r["condition"]==condition and r[field] is not None]
            ax.scatter(i+np.linspace(-.15,.15,len(values)),values,s=14,alpha=.55,c=["#287a98","#b85362","#3a8c78","#a88337"][i%4])
        ax.set_xticks(range(len(d["conditions"])),d["conditions"],rotation=15)
        ax.set_title(title)
        ax.spines[["top","right"]].set_visible(False)
    axes[1,1].set_ylim(-.03,1.03)
    axes[1,1].set_ylabel("Measured PSI")
    fig.suptitle(f"Katana-Splice-Network | {d['module']} | {d['condition']} | {d['total_events']} event nodes",fontsize=16)
    fig.savefig(target,dpi=170)
    plt.close(fig)


TEMPLATE = r'''
<style>
.event-explorer{font:15px 'IBM Plex Sans',sans-serif;color:#243f4b;margin:24px 0;background:#fffefa;border:1px solid #d6e2df;padding:24px}
.event-explorer h2{font-family:Georgia,serif;font-size:28px}.event-explorer .controls{display:flex;flex-wrap:wrap;gap:16px;margin:16px 0}
.event-explorer select{padding:8px;max-width:100%;border:1px solid #adc7c9}.event-explorer .views{display:grid;grid-template-columns:repeat(auto-fit,minmax(min(100%,380px),1fr));gap:20px}
.event-explorer article{background:#f3f7f6;padding:16px;min-width:0}.event-explorer svg{width:100%;max-height:460px}.event-explorer table{border-collapse:collapse;width:100%;font-size:12px}.event-explorer td,.event-explorer th{padding:8px;border-bottom:1px solid #d6e2df;text-align:left;overflow-wrap:anywhere}
.event-explorer .heat{overflow:auto;max-height:440px}.event-explorer .heat svg{min-width:360px;max-height:none}.event-explorer .counts{font-weight:600;margin:15px 0}.event-explorer p{line-height:1.55}.event-explorer .warning{background:#faf0d9;padding:12px}.event-explorer button{cursor:pointer;border:0;background:none;color:#08756d;text-align:left;overflow-wrap:anywhere}
</style>
<section id="__ID__" class="event-explorer"><h2>Katana-Splice-Network: Event Network Explorer</h2><p>Each node is one splice event, not a gene or protein. Gene names label events. Select a condition to examine associations within that group; click a hub or node to inspect the actual measurements.</p>
<div class="controls"><label>Module <select class="module"></select></label><label>Network condition <select class="condition"></select></label><label><input type="checkbox" class="hide"/> Hide same-gene edges in the drawing</label></div>
<p class="counts"></p><p class="warning" hidden></p><div class="views"><article><h3>Splice-Event Network</h3><div class="network"></div><p class="method"></p></article><article><h3>Sample Heatmap: Scaled PSI</h3><div class="heat"></div><p>Row z-score: pale blue = lower inclusion, pale red = higher inclusion; grey = missing. Values are scaled per event, not gene expression.</p><p class="display"></p></article><article><h3>Module Activity by Condition</h3><div class="scores"></div><p>Cached module score. Score orientation is arbitrary and does not mean more exon inclusion.</p></article><article><h3>Selected Splice Event</h3><p class="event-label"></p><div class="profile"></div><p>Each dot is one PSI measurement (0-1). Hover for sample identity. A measurement is not automatically an independent biological replicate.</p></article></div>
<h3>Hub Events: Full-Module Validation</h3><p>Complete connectivity, module membership and biological-block bootstrap; same-gene sensitivity is reported separately. Candidate status does not imply causality.</p><div class="hubs"></div><h3>Functions and Contributing Genes</h3><p class="enrichment-scope"></p><div class="enrichment"></div><details><summary>Algorithm, source and display rules</summary><p class="rules"></p><p class="interpretation"></p><pre class="source" style="white-space:pre-wrap;overflow-wrap:anywhere"></pre></details></section>
<script>(()=>{
const root=document.getElementById('__ID__'), payloads=__DATA__, $=s=>root.querySelector(s), ns='http://www.w3.org/2000/svg';
const key=d=>(d.discovery_run||'')+'|'+d.module, modules=[...new Set(payloads.map(key))], palette=['#287a98','#b85362','#3a8c78','#a88337','#646c9f'];
const element=(tag,attrs={},text)=>{const node=document.createElementNS(ns,tag);for(const[k,v]of Object.entries(attrs))node.setAttribute(k,String(v));if(text!=null)node.textContent=String(text);return node};
const option=(value,label)=>{const o=document.createElement('option');o.value=value;o.textContent=label;return o};
modules.forEach(m=>{const d=payloads.find(d=>key(d)===m);$('.module').append(option(m,d.module+(new Set(payloads.map(d=>d.discovery_run)).size>1?' / '+d.discovery_run.slice(0,8):'')))});
function chart(rows,score=false,conditionOrder=[]){const svg=element('svg',{viewBox:'0 0 410 270'}), values=rows.filter(r=>r.value!=null).map(r=>r.value),levels=conditionOrder.filter(c=>rows.some(r=>r.condition===c)), low=score?Math.min(0,...values):0,high=score?Math.max(.001,...values):1,y=v=>205-(v-low)/(high-low)*155;svg.append(element('line',{x1:45,y1:50,x2:45,y2:205,stroke:'#a3b6bb'}));[low,(low+high)/2,high].forEach(v=>svg.append(element('text',{x:40,y:y(v)+4,'text-anchor':'end','font-size':11},v.toFixed(2))));levels.forEach((c,j)=>{const group=rows.filter(r=>r.condition===c),x=60+(j+.5)*320/levels.length;group.forEach((r,i)=>{if(r.value==null)return;const point=element('circle',{cx:x+((i*17)%31-15),cy:y(r.value),r:3,fill:palette[conditionOrder.indexOf(c)%5],opacity:.6});point.append(element('title',{},r.sample_id+': '+r.value.toFixed(3)));svg.append(point)});svg.append(element('text',{x,y:230,'text-anchor':'middle','font-size':12},c));svg.append(element('text',{x,y:250,'text-anchor':'middle','font-size':10},group.length+' measurements'))});return svg}
function draw(){const m=$('.module').value,c=$('.condition').value||'all',d=payloads.find(d=>key(d)===m&&d.condition===c)||payloads.find(d=>key(d)===m);if(!d)return;
$('.counts').textContent=d.total_events+' events | '+d.biological_samples+' biological samples | '+d.total_edges+' edges | '+d.same_gene_edges+' same-gene edges';$('.method').textContent=d.method;$('.warning').hidden=d.status==='available';$('.warning').textContent=d.method;$('.rules').textContent=d.display_rule;$('.interpretation').textContent=d.interpretation;$('.source').textContent=JSON.stringify({run:d.discovery_run,sources:d.sources},null,2);
const svg=element('svg',{viewBox:'0 0 430 370'}), nodes=new Map(d.nodes.map(n=>[n.event_key,n])), pos=n=>[215+n.x*165,175+n.y*145], max=Math.max(1,...d.nodes.map(n=>n.connectivity));
d.edges.filter(e=>!$('.hide').checked||!e.same_gene).forEach(e=>{const a=pos(nodes.get(e.event_a)),b=pos(nodes.get(e.event_b));svg.append(element('line',{x1:a[0],y1:a[1],x2:b[0],y2:b[1],stroke:e.same_gene?'#bf9095':'#9bb6bd','stroke-width':.5+e.weight*2,opacity:.5}))});
d.nodes.forEach((n,i)=>{const [x,y]=pos(n),point=element('circle',{cx:x,cy:y,r:3+6*Math.sqrt(n.connectivity/max),fill:'#287a98',stroke:'white'});point.append(element('title',{},n.gene+' | '+n.event_key+' | connectivity '+n.connectivity.toFixed(2)));point.style.cursor='pointer';point.onclick=()=>profile(n.event_key,d);svg.append(point);if(i<6)svg.append(element('text',{x:x+10,y:y-5,'font-size':11},n.gene||'Event '+(i+1)))});$('.network').replaceChildren(svg);
const h=d.heatmap,hm=element('svg',{viewBox:`0 0 ${Math.max(380,h.samples.length*7)} ${h.events.length*7+18}`});h.conditions.forEach((c,j)=>hm.append(element('rect',{x:j*7,y:0,width:7,height:9,fill:palette[d.conditions.indexOf(c)%5]})));h.z.forEach((row,i)=>row.forEach((v,j)=>{let color='#d8dcde';if(v!=null){const t=Math.min(Math.abs(v)/2,1),end=v<0?[122,171,205]:[221,151,158];color='rgb('+end.map(c=>Math.round(250+(c-250)*t)).join(',')+')'}const r=element('rect',{x:j*7,y:i*7+15,width:7,height:7,fill:color});r.append(element('title',{},h.events[i]+' / '+h.samples[j]+': '+(v==null?'missing':'z='+v.toFixed(2))));hm.append(r)}));$('.heat').replaceChildren(hm);$('.display').textContent='Shown '+d.nodes.length+'/'+d.total_events+' events; '+h.samples.length+' samples. Display thinning does not change module membership.';
$('.scores').replaceChildren(d.scores.length?chart(d.scores,true,d.conditions):document.createTextNode('No cached score available.'));
const table=document.createElement('table'), header=document.createElement('tr');['Gene label','Event','Full kWithin','kME','Bootstrap all / cross-gene','Cross-gene rank'].forEach(t=>{const th=document.createElement('th');th.textContent=t;header.append(th)});table.append(header);d.hub_events.forEach(h=>{const tr=document.createElement('tr');[h.gene||'Not annotated',h.event_key,(h.kWithin??h.connectivity).toFixed(2),h.kME?.toFixed(3)??'Not computed',(h.bootstrap_frequency_all_edges?.toFixed(2)??'NA')+' / '+(h.bootstrap_frequency_no_same_gene?.toFixed(2)??'NA'),h.rank_no_same_gene??'NA'].forEach((t,i)=>{const td=document.createElement('td');if(i===1){const b=document.createElement('button');b.textContent=t;b.onclick=()=>profile(h.event_key,d);td.append(b)}else td.textContent=String(t);tr.append(td)});table.append(tr)});const notice=document.createElement('p');notice.textContent=d.hub_method||'Run the explicit hub-event analysis first.';$('.hubs').replaceChildren(notice,table);if(!d.hub_events.length){const empty=document.createElement('p');empty.textContent='No validated hub candidates: '+(d.hub_status||'not_computed');$('.hubs').append(empty)}if(d.hub_comparison?.length){const detail=document.createElement('details'),summary=document.createElement('summary'),pre=document.createElement('pre');summary.textContent='All-event ranks, same-gene comparison and rejection reasons';pre.style.whiteSpace='pre-wrap';pre.textContent=JSON.stringify(d.hub_comparison,null,2);detail.append(summary,pre);$('.hubs').append(detail)}profile(d.selected_event,d);
$('.enrichment-scope').textContent=d.enrichment_scope||'Frozen module genes; enrichment is annotation association, not causal evidence.';const terms=d.enrichment||[], list=document.createElement('div');if(!terms.length)list.textContent='No cached functional enrichment for this module.';else{const bubble=element('svg',{viewBox:`0 0 950 ${terms.length*28+24}`});terms.forEach((e,i)=>{bubble.append(element('text',{x:4,y:24+i*28,'font-size':12},e.term.slice(0,85)));const q=e.FDR??1,p=element('circle',{cx:700+Math.min(10,-Math.log10(Math.max(q,1e-300)))*18,cy:20+i*28,r:3+Math.sqrt(e.overlap_count),fill:q<.05?'#b85362':'#83abc3'});p.append(element('title',{},e.term+' | FDR='+q+' | '+e.overlap_genes));bubble.append(p)});list.append(bubble);const p=document.createElement('p');p.textContent='Bubble size: contributing genes; position: -log10(FDR). Red FDR<0.05; blue other terms. Hover to inspect contributors.';list.append(p);terms.forEach(e=>{const p=document.createElement('p');p.textContent=e.term+': '+e.overlap_genes+' | P='+e.pvalue+' | FDR='+e.FDR;list.append(p)})}$('.enrichment').replaceChildren(list);
}
function profile(event,d){$('.event-label').textContent=(d.event_options.find(e=>e.event_key===event)?.gene||'')+' | '+event;const measurements=d.event_measurements,index=measurements?.events.indexOf(event);const rows=measurements&&index>=0?measurements.samples.map((s,i)=>({...s,psi:measurements.psi[index][i]})):(event===d.selected_event?d.event_profile:[]);$('.profile').replaceChildren(rows.length?chart(rows.map(r=>({...r,value:r.psi})),false,d.conditions):document.createTextNode('No cached measurements for this event.'))}
function conditions(){const m=$('.module').value;$('.condition').replaceChildren(...payloads.filter(d=>key(d)===m).map(d=>option(d.condition,d.condition==='all'?'All samples (discovery)':d.condition)));draw()}
$('.module').onchange=conditions;$('.condition').onchange=draw;$('.hide').onchange=draw;conditions();
})();</script>
'''
