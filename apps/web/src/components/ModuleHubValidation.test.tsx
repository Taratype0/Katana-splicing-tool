import {render,screen,fireEvent} from "@testing-library/react";
import {QueryClient,QueryClientProvider} from "@tanstack/react-query";
import {it,expect,vi} from "vitest";
import {api,type ModuleDashboard as Dashboard} from "../api";
import {ModuleDashboard} from "./ModuleDashboard";

function data(): Dashboard {
  return {module:"M1",condition:"all",conditions:["A"],status:"available",method:"Cached hub-analysis edges",
    biological_samples:30,minimum_biological_samples:15,total_events:2,total_edges:0,same_gene_edges:0,
    nodes:[{event_key:"e1",gene:"G1",connectivity:0,x:0,y:0}],edges:[],hub_events:[],
    hub_status:"no_supported_edges",hub_method:"Full connectivity and biological-block bootstrap",
    event_options:[{event_key:"e1",gene:"G1"}],selected_event:"e1",event_profile:[],scores:[],
    heatmap:{events:["e1"],samples:["s1"],conditions:["A"],z:[[0]],scale:"Row z-score"},
    score_definition:"Cached PC1",display_rule:"Drawing only",interpretation:"Descriptive",discovery_run:"fixed",
    connectivity_definition:"Full kWithin",enrichment:[],enrichment_total:0,enrichment_scope:"Whole module",sources:[]};
}

it("does not show zero-connectivity hubs when there are no edges",async()=>{
  vi.spyOn(api,"moduleDashboard").mockResolvedValue(data());
  const run=vi.spyOn(api,"runWorkflow");
  render(<QueryClientProvider client={new QueryClient()}><ModuleDashboard projectId="p" artifactId="a" module="M1"/></QueryClientProvider>);
  expect(await screen.findByText(/No optional labelled hub candidates/)).toHaveTextContent("no_supported_edges");
  expect(screen.queryByRole("columnheader",{name:"Full kWithin"})).not.toBeInTheDocument();
  expect(run).not.toHaveBeenCalled();vi.restoreAllMocks();
});

it("shows persisted full metrics and same-gene rank comparison without launching analysis",async()=>{
  const d=data();d.total_edges=10;d.hub_status="available";
  const hub={event_key:"e1",gene:"G1",connectivity:7.125,kWithin:7.125,kWithin_no_same_gene:5.25,kME:.93,
    rank_all_edges:1,rank_no_same_gene:2,rank_change_no_same_gene:1,bootstrap_frequency_all_edges:.9,
    bootstrap_frequency_no_same_gene:.8,selection_reason:"passed",hub_candidate:true};
  d.hub_events=[hub];d.hub_comparison=[hub];
  vi.spyOn(api,"moduleDashboard").mockResolvedValue(d);const run=vi.spyOn(api,"runWorkflow");
  render(<QueryClientProvider client={new QueryClient()}><ModuleDashboard projectId="p" artifactId="a" module="M1"/></QueryClientProvider>);
  expect((await screen.findAllByRole("columnheader",{name:"Full kWithin"})).length).toBeGreaterThanOrEqual(1);
  expect(screen.getByText("0.90 / 0.80 (all / cross-gene)")).toBeInTheDocument();
  expect(screen.getByRole("columnheader",{name:"Rank change"})).toBeInTheDocument();
  expect(screen.getByText(/not recorded; passed/)).toBeInTheDocument();expect(run).not.toHaveBeenCalled();vi.restoreAllMocks();
});

it("shows every continuous metric and sortable rank while optional hub labels remain NA",async()=>{
  const d=data();d.total_events=25;d.hub_status="labels_not_applied";d.hub_label_status="not_applied";
  d.hub_comparison=Array.from({length:25},(_,i)=>({event_key:`event-${i+1}`,gene:`Gene${i+1}`,connectivity:25-i,
    kWithin:25-i,kWithin_no_same_gene:20-i/2,kME:.9-i/100,rank_all_edges:i+1,rank_no_same_gene:i+1,
    hub_candidate:null,label_status:"not_applied",estimation_status:"estimable"}));
  vi.spyOn(api,"moduleDashboard").mockResolvedValue(d);const run=vi.spyOn(api,"runWorkflow");
  render(<QueryClientProvider client={new QueryClient()}><ModuleDashboard projectId="p" artifactId="a" module="M1"/></QueryClientProvider>);
  expect(await screen.findByText(/Optional threshold labels were not applied/)).toBeInTheDocument();
  const section=screen.getByRole("region",{name:"Complete hub metrics"});
  expect(section.querySelectorAll("tbody tr")).toHaveLength(25);
  expect(section.querySelector("tbody tr")?.textContent).toContain("event-1");
  fireEvent.click(screen.getByRole("button",{name:/^Full rank$/}));
  expect(section.querySelector("tbody tr")?.textContent).toContain("event-25");
  expect(section.querySelectorAll("tbody tr")).toHaveLength(25);
  expect(screen.getByRole("button",{name:"Optional labelled hub candidates TSV"})).toBeDisabled();
  expect(screen.queryByText(/No optional labelled hub candidates/)).not.toBeInTheDocument();
  expect(screen.getAllByText("not_applied; estimable")).toHaveLength(25);
  expect(screen.getByText(/not a probability, P value or FDR/)).toBeInTheDocument();
  expect(run).not.toHaveBeenCalled();vi.restoreAllMocks();
});

it("uses inclusive FDR 0.05 for the functional-label boundary without changing statistics",async()=>{
  const d=data();d.enrichment=[{term:"Boundary term",overlap_count:2,overlap_genes:"G1;G2",pvalue:.01,FDR:.05}];d.enrichment_total=1;
  vi.spyOn(api,"moduleDashboard").mockResolvedValue(d);
  render(<QueryClientProvider client={new QueryClient()}><ModuleDashboard projectId="p" artifactId="a" module="M1"/></QueryClientProvider>);
  const plot=await screen.findByRole("img",{name:"Module functional enrichment bubble plot"});
  expect(plot.querySelector("circle")).toHaveAttribute("fill","#b85362");
  expect(screen.getByText(/Red: FDR≤0.05/)).toBeInTheDocument();
  expect(d.enrichment[0].FDR).toBe(.05);vi.restoreAllMocks();
});
