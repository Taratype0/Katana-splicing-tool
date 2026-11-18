import {render,screen,fireEvent,waitFor} from '@testing-library/react';
import {QueryClient,QueryClientProvider} from '@tanstack/react-query';
import {vi,test,expect} from 'vitest';
import {ReportExplorer} from './ReportExplorer';
import {api} from '../api';
import type {Artifact} from '../types';

vi.mock('./PopulationPlot',()=>({PopulationPlot:()=> <div>Application volcano</div>}));

test('threshold filtering preserves the full volcano background and shows mapped names in the table',async()=>{
 const browse=vi.spyOn(api,'browseReport').mockResolvedValue({columns:['row_id','gene_name','gene_id','gene_symbol','deg_fdr','log2fc'],rows:[{row_id:'r1',gene_name:'Fos',gene_symbol:'ENSMUSG1.2',gene_id:'ENSMUSG1',deg_fdr:.01,log2fc:2}],total:1,offset:0,limit:100});
 const plot=vi.spyOn(api,'plotReport').mockResolvedValue({mode:'points',total:1,plotted:1,unplottable:0,x_min:-2,x_max:2,y_max:2,points:[[2,2,'r1','Fos']]});
 const artifact={artifact_id:'names',analysis_id:'imported.report',name:'Complete DEG',row_count:1,metadata:{},parameters:{},contrast_ids:[]} as unknown as Artifact;
 render(<QueryClientProvider client={new QueryClient({defaultOptions:{queries:{retry:false}}})}><ReportExplorer projectId="names" artifacts={[artifact]} onAdd={vi.fn()}/></QueryClientProvider>);
 await screen.findByText('Application volcano');
 expect(screen.getByLabelText('Maximum FDR')).toHaveValue(.05);
 expect(screen.getByLabelText('Minimum absolute effect')).toHaveValue(1);
 fireEvent.change(screen.getByLabelText('Maximum FDR'),{target:{value:'.05'}});
 await waitFor(()=>expect(browse.mock.calls.at(-1)?.[2]).toMatchObject({fdr_max:.05}));
 expect(plot.mock.calls.at(-1)?.[2]).not.toHaveProperty('fdr_max');
 fireEvent.click(screen.getByRole('tab',{name:'Complete result table'}));
 expect(await screen.findByText('Fos')).toBeInTheDocument();
 expect(screen.queryByText('ENSMUSG1.2')).not.toBeInTheDocument();
 fireEvent.click(screen.getByLabelText('Show all original columns and mapping provenance'));
 expect(screen.getByText('ENSMUSG1.2')).toBeInTheDocument();
 browse.mockRestore();plot.mockRestore();
});

test('a delayed full selection cannot restore rows from a previous filter',async()=>{
 const browse=vi.spyOn(api,'browseReport').mockResolvedValue({columns:['row_id','gene_id','deg_fdr'],rows:[{row_id:'r1',gene_id:'G1',deg_fdr:.01}],total:1,offset:0,limit:100});
 let finish!:(value:{rows:Record<string,unknown>[];total:number})=>void;
 const select=vi.spyOn(api,'selectReport').mockImplementation(()=>new Promise(resolve=>{finish=resolve;}));
 const artifact={artifact_id:'race',analysis_id:'imported.report',name:'Real DEG schema',row_count:1,metadata:{},parameters:{},contrast_ids:[]} as unknown as Artifact;
 render(<QueryClientProvider client={new QueryClient({defaultOptions:{queries:{retry:false}}})}><ReportExplorer projectId="race-study" artifacts={[artifact]} onAdd={vi.fn()}/></QueryClientProvider>);
 await screen.findByText('0.01');
 expect(screen.getByLabelText('FDR column')).toHaveValue('deg_fdr');
 fireEvent.click(screen.getByRole('button',{name:'Select all 1 matching rows'}));
 fireEvent.change(screen.getByLabelText('Search'),{target:{value:'different'}});
 finish({rows:[{row_id:'r1',gene_id:'G1'}],total:1});
 await waitFor(()=>expect(screen.queryByText('Selecting…')).not.toBeInTheDocument());
 expect(screen.getByRole('button',{name:'Add genes to targets'})).toBeDisabled();
 browse.mockRestore();select.mockRestore();
});
test('complete table selection retains exact non-significant rows and browsing context', async()=>{
 const browse=vi.spyOn(api,'browseReport').mockResolvedValue({columns:['row_id','gene_id','event_id','fdr'],rows:[{row_id:'r1',gene_id:'G1',event_id:'SE:1',fdr:.9}],total:1,offset:0,limit:100});
 const add=vi.fn();
 const artifact={artifact_id:'a',analysis_id:'imported.report',name:'Complete events',row_count:1,metadata:{report_hash_verified:true},parameters:{},contrast_ids:['B_vs_A']} as unknown as Artifact;
 render(<QueryClientProvider client={new QueryClient({defaultOptions:{queries:{retry:false}}})}><ReportExplorer projectId="study" artifacts={[artifact]} onAdd={add}/></QueryClientProvider>);
 await screen.findByText('0.9');
 fireEvent.click(screen.getByRole('button',{name:'Select this page'}));
 fireEvent.click(screen.getByRole('button',{name:'Add events to targets'}));
 expect(add.mock.calls[0][0]).toMatchObject({source_artifact_id:'a',target_kind:'event',row_ids:['r1']});
 fireEvent.change(screen.getByLabelText('Gene IDs, symbols or event IDs'),{target:{value:'G1 G2'}});
 await waitFor(()=>expect(browse.mock.calls.at(-1)?.[2]).toMatchObject({identifiers:['G1','G2']}));
 expect(screen.getByRole('button',{name:'Add events to targets'})).toBeDisabled();
 browse.mockRestore();
});

test('select all matching rows exports the full population rather than the page',async()=>{
 const all=Array.from({length:305},(_,i)=>({row_id:`r${i}`,gene_id:`G${i}`}));
 const browse=vi.spyOn(api,'browseReport').mockResolvedValue({columns:['row_id','gene_id'],rows:all.slice(0,100),total:305,offset:0,limit:100});
 const select=vi.spyOn(api,'selectReport').mockResolvedValue({rows:all,total:305});
 const add=vi.fn();const artifact={artifact_id:'full',analysis_id:'imported.report',name:'Complete genes',row_count:305,metadata:{},parameters:{},contrast_ids:[]} as unknown as Artifact;
 render(<QueryClientProvider client={new QueryClient({defaultOptions:{queries:{retry:false}}})}><ReportExplorer projectId="full-study" artifacts={[artifact]} onAdd={add}/></QueryClientProvider>);
 fireEvent.click(await screen.findByRole('button',{name:'Select all 305 matching rows'}));
 await waitFor(()=>expect(screen.getByRole('button',{name:'Add genes to targets'})).toBeEnabled());
 fireEvent.click(screen.getByRole('button',{name:'Add genes to targets'}));
 expect(add.mock.calls[0][0].row_ids).toHaveLength(305);
 expect(add.mock.calls[0][0].row_ids).toContain('r304');
 browse.mockRestore();select.mockRestore();
});
