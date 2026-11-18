import {useEffect, useRef, useState} from 'react';
import type {Result} from 'vega-embed';

export type Population = {mode:'points'|'density'; total:number; plotted:number; unplottable:number;
  x_min:number;x_max:number;y_max:number;points?:[number,number,string,string,boolean?,number?][];bins?:[number,number,number][];grid?:[number,number]};
export type PlotRange = {effect_low:number;effect_high:number;fdr_min:number;fdr_max:number};

// Use the application's Vega-Lite grammar and export controls for imported data.
export function volcanoSpec(data:Population, effect:string, fdr:string, fdrCutoff:number, effectCutoff:number, labels:boolean, zoom=false, readSupport=false) {
  const extent=Math.max(Math.abs(data.x_min),Math.abs(data.x_max),effectCutoff*1.2,1);
  const yLine=-Math.log10(Math.max(fdrCutoff,1e-300));
  const ymax=Math.max(data.y_max*1.04,yLine*1.15,1);
  const rows=(data.points||[]).map(([x,y,row_id,gene,supported,rawP])=>({x,y,row_id,gene,p:rawP??10**-y,
    status:(!readSupport||supported!==false)&&y>=yLine&&Math.abs(x)>=effectCutoff?(x>=0?'Higher in numerator':'Lower in numerator'):'Other measured rows'}));
  const x={field:'x',type:'quantitative',scale:{domain:[-extent,extent],nice:false},axis:{title:effect==='log2fc'||effect.toLowerCase().includes('fold')?'log2 fold change (numerator / denominator)':effect,values:[-extent,-extent/2,0,extent/2,extent],format:'.2~f',grid:true}};
  const y={field:'y',type:'quantitative',scale:{domain:[0,ymax],nice:false},axis:{title:`−log10(${fdr})`,tickMinStep:1,grid:true}};
  const points=data.mode==='points'?{data:{values:rows},params:[{name:'region',select:{type:'interval',encodings:['x','y'],clear:'dblclick'},...(zoom?{bind:'scales'}:{})},{name:'target',select:{type:'point',fields:['row_id'],on:'click',clear:false}}],mark:{type:'point',filled:true,opacity:.68,size:34,clip:true},encoding:{x,y,color:{field:'status',type:'nominal',scale:{domain:['Higher in numerator','Lower in numerator','Other measured rows'],range:['#2251FF','#B83A4B','#CBD5E1']},legend:{title:null,orient:'top'}},tooltip:[{field:'gene',title:'Gene'},{field:'p',title:fdr,format:'.4g'},{field:'row_id',title:'Row'},{field:'x',title:effect,format:'.4g'},{field:'y',title:`−log10(${fdr})`,format:'.4g'}]}}:
    {data:{values:(data.bins||[]).map(([bx,by,n])=>({x:data.x_min+bx/240*(data.x_max-data.x_min),x2:data.x_min+(bx+1)/240*(data.x_max-data.x_min),y:by/120*data.y_max,y2:(by+1)/120*data.y_max,n}))},mark:'rect',encoding:{x,y,x2:{field:'x2'},y2:{field:'y2'},color:{field:'n',type:'quantitative',scale:{type:'log',scheme:'blues'},legend:{title:'All rows / bin'}},tooltip:[{field:'n',title:'Rows'}]}};
  const seen=new Set<string>();const labelPositions:{x:number;y:number;gene:string}[]=[];const selectedLabels=rows.filter(r=>r.status!=='Other measured rows').sort((a,b)=>b.y-a.y||Math.abs(b.x)-Math.abs(a.x)).filter(r=>{if(seen.has(r.gene))return false;seen.add(r.gene);if(labelPositions.some(p=>Math.abs(p.y-r.y)/ymax*390<16&&Math.abs(p.x-r.x)/(2*extent)*900<(p.gene.length+r.gene.length)*3.2+8))return false;labelPositions.push(r);return true;}).slice(0,12);
  return {$schema:'https://vega.github.io/schema/vega-lite/v5.json',width:'container',height:390,
    config:{font:'IBM Plex Sans, Arial, sans-serif',view:{stroke:null},axis:{labelFontSize:11.5,titleFontSize:12.5,labelColor:'#1E293B',titleColor:'#051C2C',gridColor:'#F1F5F9'},legend:{labelFontSize:11.5,labelLimit:280}},
    layer:[points,
      {data:{values:[{x:0}]},mark:{type:'rule',color:'#94A3B8',strokeWidth:1},encoding:{x}},
      {data:{values:[{x:-effectCutoff},{x:effectCutoff}]},mark:{type:'rule',color:'#64748B',strokeDash:[5,4]},encoding:{x}},
      {data:{values:[{y:yLine}]},mark:{type:'rule',color:'#64748B',strokeDash:[5,4]},encoding:{y}},
      {data:{values:[{x:extent*.98,y:yLine,threshold:`FDR = ${fdrCutoff}`}]},mark:{type:'text',align:'right',dy:-8,fontSize:11,color:'#051C2C'},encoding:{x,y,text:{field:'threshold'}}},
      ...(labels?[{data:{values:selectedLabels},mark:{type:'text',dy:-9,fontSize:11,fontWeight:600,color:'#051C2C',clip:true},encoding:{x,y,text:{field:'gene'}}}]:[])]};
}

export function PopulationPlot({data,effect,fdr,onRange,onSelect,fdrCutoff=.05,effectCutoff=1,readSupport=false}:{data:Population;effect:string;fdr:string;
  onRange:(range:PlotRange)=>void;onSelect:(rowId:string)=>void;fdrCutoff?:number;effectCutoff?:number;readSupport?:boolean}) {
  const host=useRef<HTMLDivElement>(null),embedded=useRef<Result|null>(null);
  const callbacks=useRef({onRange,onSelect});callbacks.current={onRange,onSelect};
  const [zoom,setZoom]=useState(false),[revision,setRevision]=useState(0);
  const [labels,setLabels]=useState(true),[error,setError]=useState('');const [region,setRegion]=useState<PlotRange>();
  useEffect(()=>{
    let disposed=false;const element=host.current;if(!element)return;
    setError('');setRegion(undefined);element.replaceChildren();
    void import('vega-embed').then(({default:embed})=>embed(element,volcanoSpec(data,effect,fdr,fdrCutoff,effectCutoff,labels,zoom,readSupport) as never,
      {actions:{export:true,source:false,compiled:false,editor:false},renderer:data.plotted>30000?'canvas':'svg'})).then(result=>{
      if(disposed){result.finalize();return;}embedded.current=result;
      if(data.mode==='points'){
        result.view.addSignalListener('region',(_name,value)=>{const v=value as {x?:number[];y?:number[]};if(v.x&&v.y)setRegion({effect_low:v.x[0],effect_high:v.x[1],fdr_min:v.y[1]>=data.y_max?0:10**-v.y[1],fdr_max:Math.min(1,10**-v.y[0])});else setRegion(undefined);});
        result.view.addEventListener('click',(_event,item)=>{const id=item?.datum?.row_id;if(id)callbacks.current.onSelect(String(id));});
      }
    }).catch(e=>{if(!disposed)setError(String(e));});
    const observer=new ResizeObserver(()=>{void embedded.current?.view.resize().runAsync();});observer.observe(element);
    return()=>{disposed=true;observer.disconnect();embedded.current?.finalize();embedded.current=null;element.replaceChildren();};
  },[data,effect,fdr,fdrCutoff,effectCutoff,labels,zoom,revision,readSupport]);
  return <div className="plot-stage">
    <div className="preview-toolbar"><strong>Volcano plot</strong><button aria-pressed={zoom} onClick={()=>setZoom(!zoom)}>{zoom?"Zoom / pan mode":"Brush selection mode"}</button><button onClick={()=>setRevision(v=>v+1)}>Reset view</button><label><input type="checkbox" checked={labels} onChange={e=>setLabels(e.target.checked)}/> Label leading genes</label><button disabled={!region} onClick={()=>region&&onRange(region)}>Filter to brushed region</button></div>
    <p className="plot-summary">{data.plotted.toLocaleString()} plotted / {data.total.toLocaleString()} rows · {data.unplottable.toLocaleString()} missing or invalid coordinates · FDR ≤ {fdrCutoff} · |effect| ≥ {effectCutoff}</p>
    {error&&<p role="alert">{error}</p>}<div ref={host} className="plot-host" aria-label="Full-population volcano plot"/>
    <p className="plot-summary">All valid results in the current comparison and biological scope remain in the plot. Dashed lines show significance cutoffs; the horizontal axis is symmetric around 0. Hover for gene names; click to select a row. {readSupport?'PS-DAS highlights require recorded read support; grey points also include unsupported events.':''} Drag to brush, then apply the region. {data.mode==='density'?'Full-population density above 150,000 rows; no sampling. ':''}Zero FDR is displayed at −log10 = 300. Up to 12 leading genes are labelled; the data are not truncated.</p>
  </div>;
}
