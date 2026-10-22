import { fireEvent, render, screen } from "@testing-library/react";
import { expect, it, vi } from "vitest";
import type { WorkflowDefinition } from "../types";
import { WorkflowControls } from "./WorkflowControls";

it("shows restored JSON resource mappings instead of object coercion", () => {
  const workflow: WorkflowDefinition = {
    workflow_id: "programs.module_functions", workspace_id: "programs",
    label: "Module functions", question: "Which functions?", method: "Hypergeometric",
    output_summary: "Contributing genes", analysis_unit: "module", expected_artifacts: 2,
    required_asset_kinds: [], optional_asset_kinds: [], supports_contrast_selection: false,
    enabled: true, limitations: "Annotation association",
    parameters: [{key: "libraries", label: "GMT resources", kind: "text", default: "{}", help: "Explicit paths", options: []}],
  };
  const change = vi.fn();
  render(<WorkflowControls workflow={workflow} contrasts={[]} selectedContrasts={new Set()}
    parameters={{libraries: {GO_BP: "E:/resources/GO.gmt"}}} selections={[]} conditions={[]}
    onToggleContrast={vi.fn()} onSelectAllContrasts={vi.fn()} onClearContrasts={vi.fn()}
    onParameterChange={change} onReset={vi.fn()} />);
  const field = screen.getByRole("textbox");
  expect(field).toHaveValue('{"GO_BP":"E:/resources/GO.gmt"}');
  fireEvent.change(field, {target: {value: '{"KEGG":"E:/resources/KEGG.gmt"}'}});
  expect(change).toHaveBeenCalledWith("libraries", '{"KEGG":"E:/resources/KEGG.gmt"}');
});

it("MotifTrace exposes custom-library controls, not unused legacy inference thresholds", () => {
  const workflow: WorkflowDefinition = {
    workflow_id: "prioritize.rbp_mechanism", workspace_id: "prioritize", label: "RBP motifs",
    question: "Where are sequence matches?", method: "Completed frozen results", output_summary: "RNA maps",
    analysis_unit: "event", expected_artifacts: 1, required_asset_kinds: [], optional_asset_kinds: [],
    supports_contrast_selection: false, enabled: true, limitations: "Not measured binding",
    parameters: [
      {key:"sequence_scope", label:"RNA sequence analysis", kind:"select", default:"boundary_association", help:"Method", options:["boundary_association","completed_sequence_maps"]},
      {key:"site_pvalue", label:"FIMO site P threshold", kind:"number", default:.0001, help:"Legacy only", options:[]},
      {key:"binding_mode", label:"Experimental binding evidence", kind:"select", default:"required", help:"Legacy only", options:["required"]},
      {key:"bundle_asset_id", label:"Completed result asset ID", kind:"text", default:"", help:"Registered result", options:[]},
      {key:"window_bundle_asset_id", label:"Matched count / coverage result asset ID (optional)", kind:"text", default:"", help:"Reuse separate original FDR families", options:[]},
      {key:"view_scope", label:"Display comparisons", kind:"select", default:"selected_comparison", help:"No FDR changes", options:["selected_comparison","all_frozen_comparisons"]},
      {key:"custom_motif_path", label:"User motif library path (optional)", kind:"text", default:"", help:"Explicit local file", options:[]},
      {key:"custom_motif_mode", label:"User motif analysis", kind:"select", default:"preview_only", help:"Separate family only when explicit", options:["preview_only","fixed_enrichment"]},
    ],
  };
  render(<WorkflowControls workflow={workflow} contrasts={[]} selectedContrasts={new Set()}
    parameters={{sequence_scope:"completed_sequence_maps",custom_motif_mode:"preview_only"}} selections={[]} conditions={[]}
    onToggleContrast={vi.fn()} onSelectAllContrasts={vi.fn()} onClearContrasts={vi.fn()}
    onParameterChange={vi.fn()} onReset={vi.fn()} />);
  expect(screen.queryByText("FIMO site P threshold")).not.toBeInTheDocument();
  expect(screen.queryByText("Experimental binding evidence")).not.toBeInTheDocument();
  expect(screen.getByText("Completed result asset ID")).toBeInTheDocument();
  expect(screen.getByText("Matched count / coverage result asset ID (optional)")).toBeInTheDocument();
  expect(screen.getByText("Display comparisons")).toBeInTheDocument();
  expect(screen.getByRole("option", {name:"Katana MotifTrace: completed SE RNA and motifs"})).toBeInTheDocument();
  expect(screen.getByText("User motif library path (optional)")).toBeInTheDocument();
  expect(screen.getByRole("option", {name:"User library: sequence positions only (no new P/FDR)"})).toBeInTheDocument();
});

it("MotifTrace consequence mode shows only the new fixed-family controls", () => {
  const workflow: WorkflowDefinition = {
    workflow_id:"prioritize.rbp_mechanism",workspace_id:"prioritize",label:"Katana MotifTrace",
    question:"RNA consequences",method:"Paired swaps",output_summary:"Five types",analysis_unit:"gene pair",
    expected_artifacts:1,required_asset_kinds:[],optional_asset_kinds:[],supports_contrast_selection:false,enabled:true,limitations:"Not binding",
    parameters:[
      {key:"sequence_scope",label:"RNA sequence analysis",kind:"select",default:"boundary_association",help:"Method",options:["boundary_association","event_consequences"]},
      ...["custom_motif_path","model_ids","bundle_asset_id","custom_motif_mode","site_pvalue","binding_mode"].map(key=>({key,label:key,kind:"text" as const,default:"",help:"Method",options:[]})),
      {key:"cluster_permutations",label:"Paired permutations",kind:"number" as const,default:999,help:"New family",options:[]},
    ],
  };
  render(<WorkflowControls workflow={workflow} contrasts={[]} selectedContrasts={new Set()} parameters={{sequence_scope:"event_consequences"}}
    selections={[]} conditions={[]} onToggleContrast={vi.fn()} onSelectAllContrasts={vi.fn()} onClearContrasts={vi.fn()} onParameterChange={vi.fn()} onReset={vi.fn()}/>);
  expect(screen.getByText("custom_motif_path")).toBeInTheDocument();
  expect(screen.getByText("model_ids")).toBeInTheDocument();
  expect(screen.getByText("Paired permutations")).toBeInTheDocument();
  for(const key of ["bundle_asset_id","custom_motif_mode","site_pvalue","binding_mode"])expect(screen.queryByText(key)).not.toBeInTheDocument();
});
