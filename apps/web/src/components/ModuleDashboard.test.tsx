import { render, screen, fireEvent, waitFor } from "@testing-library/react";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { expect, it, vi } from "vitest";
import { api, type ModuleDashboard as Dashboard } from "../api";
import { ModuleDashboard } from "./ModuleDashboard";

it("links condition, event, heatmap, hubs and enrichment without launching analysis", async () => {
  const data = { module: "M1", condition: "all", conditions: ["normal", "tumor"], node_unit: "splice_event", status: "available",
    biological_samples: 40, total_events: 2, total_edges: 1, same_gene_edges: 0, method: "Cached discovery edges", discovery_run: "frozen",
    nodes: [{event_key:"e1",gene:"A",connectivity:1,x:0,y:0},{event_key:"e2",gene:"B",connectivity:1,x:.5,y:.5}],
    edges: [{event_a:"e1",event_b:"e2",weight:.9,correlation:.95,pair_n:40,same_gene:false}],
    hub_events: [{event_key:"e1",gene:"A",connectivity:1}], event_options: [{event_key:"e1",gene:"A"},{event_key:"e2",gene:"B"}],
    heatmap:{events:["e1"],samples:["s1","s2"],conditions:["normal","tumor"],z:[[-1,1]],scale:"Row z-score"},
    scores:[{sample_id:"s2",condition:"tumor",value:1},{sample_id:"s1",condition:"normal",value:-1}],
    event_profile:[{sample_id:"s2",biological_sample_id:"d2",condition:"tumor",psi:.8},{sample_id:"s1",biological_sample_id:"d1",condition:"normal",psi:.2}],
    selected_event:"e1", score_definition:"Cached score",display_rule:"Display only",interpretation:"Not PPI",sources:[],
    enrichment:[{term:"RNA processing",overlap_count:2,overlap_genes:"A;B",pvalue:.001,FDR:.01}],enrichment_total:20,enrichment_scope:"Measured network background"} as unknown as Dashboard;
  const read = vi.spyOn(api,"moduleDashboard").mockImplementation(async (_p,_a,_m,c,e)=>({...data,condition:c||"all",selected_event:e||"e1"}));
  const run = vi.spyOn(api,"runWorkflow");
  render(<QueryClientProvider client={new QueryClient({defaultOptions:{queries:{retry:false}}})}><ModuleDashboard projectId="p" artifactId="a" module="M1"/></QueryClientProvider>);
  expect(await screen.findByRole("img",{name:"Row-scaled PSI heatmap"})).toBeInTheDocument();
  expect(screen.getByRole("img",{name:"Module functional enrichment bubble plot"})).toBeInTheDocument();
  expect(screen.getByRole("img",{name:"Module score by condition"}).querySelectorAll("circle")[0]).toHaveAttribute("fill", "#287a98");
  expect(screen.getByRole("img",{name:"Event PSI by condition"}).querySelectorAll("circle")[1]).toHaveAttribute("fill", "#b85362");
  fireEvent.change(screen.getByLabelText("Network condition"),{target:{value:"tumor"}});
  await waitFor(()=>expect(read).toHaveBeenLastCalledWith("p","a","M1","tumor",undefined,false));
  fireEvent.click(screen.getByRole("button",{name:"Select event e2"}));
  await waitFor(()=>expect(read).toHaveBeenLastCalledWith("p","a","M1","tumor","e2",false));
  fireEvent.click(screen.getByLabelText("Hide same-gene edges in drawing"));
  await waitFor(()=>expect(read).toHaveBeenLastCalledWith("p","a","M1","tumor","e2",true));
  expect(run).not.toHaveBeenCalled();
  vi.restoreAllMocks();
});
