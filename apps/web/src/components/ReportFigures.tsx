import {useEffect,useRef,useState} from 'react';
import {useQuery} from '@tanstack/react-query';
import type {Result} from 'vega-embed';
export type ReportPanel={id:string;chapter:number;title:string;image?:string;svg?:string;contrast_id?:string;caption?:string;data_links?:{path:string;label?:string;role:string}[];metadata?:Record<string,unknown>};
export const reportUrl=(project:string,path:string)=>`/api/projects/${encodeURIComponent(project)}/report-source/${path.split('/').map(encodeURIComponent).join('/')}`;
export async function reportGet<T>(project:string,route:string):Promise<T>{const response=await fetch(`/api/projects/${encodeURIComponent(project)}/${route}`);if(!response.ok)throw new Error(await response.text());return response.json();}
export function VegaFigure({spec,title,onPoint}:{spec:Record<string,unknown>;title:string;onPoint?:(row:Record<string,unknown>)=>void}){
 const host=useRef<HTMLDivElement>(null);const callback=useRef(onPoint);callback.current=onPoint;const [error,setError]=useState('');
 useEffect(()=>{let disposed=false;let result:Result|undefined;const node=host.current;if(!node)return;node.replaceChildren();setError('');
  const updated=structuredClone(spec);
  function recolor(value:unknown){if(!value||typeof value!=='object')return;for(const [key,item] of Object.entries(value)){if(['data','datasets','values'].includes(key))continue;if(key==='range'&&Array.isArray(item)&&item.length&&item.every(c=>typeof c==='string'&&/^#[a-f0-9]{3,8}$/i.test(c))){(value as Record<string,unknown>)[key]=item.length===3?['#162e51','#f1f5fa','#0077b6']:item.map((_,i)=>['#162e51','#0077b6','#76b6d9','#526b86','#b8d0e3'][i%5]);}else if(Array.isArray(item))item.forEach(recolor);else recolor(item);}}
  recolor(updated);updated.config={...(updated.config as object||{}),background:'white',font:'Arial',range:{category:['#2863a5','#2d7e9b','#8bb9d7','#1e2630','#5b6470']}};
  void import('vega-embed').then(({default:embed})=>embed(node,updated as never,{renderer:'canvas',actions:{export:true,source:false,compiled:false,editor:false}})).then(r=>{if(disposed){r.finalize();return;}result=r;r.view.addEventListener('click',(_event,item)=>{if(item?.datum)callback.current?.(item.datum);});}).catch(e=>{if(!disposed)setError(String(e));});
  const observer=new ResizeObserver(()=>{void result?.view.resize().runAsync();});observer.observe(node);
  return()=>{disposed=true;observer.disconnect();result?.finalize();node.replaceChildren();};
 },[spec]);
 return <>{error&&<p role="alert">{error}</p>}<div className="plot-host native-figure" ref={host} aria-label={title}/></>;
}
export function ReportFigures({projectId,panels}:{projectId:string;panels:ReportPanel[]}){
 const [chosen,setChosen]=useState('');const unique=panels.filter((p,i)=>!panels.slice(0,i).some(a=>(a.image&&a.image===p.image)||a.title.replace(/ — .*$/,'')===p.title.replace(/ — .*$/,'')));const panel=unique.find(p=>p.id===chosen)||unique[0];
 const query=useQuery({queryKey:['delivered-figure',projectId,panel?.id],queryFn:()=>reportGet<{spec:Record<string,unknown>|null;reason?:string}>(projectId,`report-figures/${encodeURIComponent(panel!.id)}`),enabled:!!panel});
 if(!panel)return null;
 return <section className="result-card delivered-figures"><h3>Visualizations</h3><div className="preview-toolbar"><label>Figure<select aria-label="Delivered visualization" value={panel.id} onChange={e=>setChosen(e.target.value)}>{unique.map(p=><option key={p.id} value={p.id}>{p.title}</option>)}</select></label><span>{unique.length} delivered figures</span></div>
  <h4>{panel.title}</h4>{query.isLoading?<p>Loading figure…</p>:query.data?.spec?<VegaFigure spec={query.data.spec} title={panel.title}/>:panel.image?<img className="delivered-image" src={reportUrl(projectId,panel.image)} alt={panel.title}/>:<p>No figure delivered.</p>}
  {query.error&&<p role="alert">{query.error.message}</p>}<p>{panel.caption}</p><div className="target-actions">{panel.svg&&<a href={reportUrl(projectId,panel.svg)}>Download SVG</a>}{panel.image&&<a href={reportUrl(projectId,panel.image)}>Download figure</a>}{panel.data_links?.map((link,i)=><a key={i} href={reportUrl(projectId,link.path)}>{link.label||link.role}</a>)}</div>
 </section>;
}
