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
  expect(await screen.findByText(/No validated hub candidates/)).toHaveTextContent("no_supported_edges");
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
  fireEvent.click(screen.getByText(/All events: same-gene sensitivity/));
  expect(screen.getByRole("columnheader",{name:"Rank change"})).toBeInTheDocument();
  expect(screen.getByText("passed")).toBeInTheDocument();expect(run).not.toHaveBeenCalled();vi.restoreAllMocks();
});
