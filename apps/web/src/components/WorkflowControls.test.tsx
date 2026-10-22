import { fireEvent, render, screen } from "@testing-library/react";
import { expect, it, vi } from "vitest";
import type { Artifact, WorkflowDefinition, WorkflowParameter } from "../types";
import { WorkflowControls } from "./WorkflowControls";

const readSupportWorkflow: WorkflowDefinition = {
  workflow_id: "quality.stable_as", workspace_id: "quality", label: "AS read support",
  question: "Which events have read support in every selected sample?",
  method: "IJC + SJC of at least 10 in every selected sample across both groups.",
  output_summary: "Read-support audit and separately labelled differential results",
  analysis_unit: "event within the selected contrast", expected_artifacts: 3,
  required_asset_kinds: ["rmats_result"], optional_asset_kinds: [],
  supports_contrast_selection: true, enabled: true, limitations: "Read support is not independent validation.",
  parameters: [
    {key: "fdr", label: "rMATS FDR", kind: "number", default: 0.05, help: "Adjusted caller P-value.", minimum: 0, maximum: 1, step: 0.01, options: []},
    {key: "min_abs_dpsi", label: "Minimum |dPSI|", kind: "number", default: 0.1, help: "Absolute effect size.", minimum: 0, maximum: 1, step: 0.01, options: []},
  ],
};

it("explains support in every selected sample and keeps differential cutoffs editable", () => {
  const change = vi.fn();
  const contrast = {contrast_id: "treated-vs-control", label: "Treated vs control", numerator: {label: "Treated", sample_ids: ["t1", "t2"]}, denominator: {label: "Control", sample_ids: ["c1", "c2"]}};
  render(<WorkflowControls workflow={readSupportWorkflow} contrasts={[contrast]} selectedContrasts={new Set([contrast.contrast_id])}
    parameters={{fdr: 0.05, min_abs_dpsi: 0.1}} selections={[]} conditions={[]}
    onToggleContrast={vi.fn()} onSelectAllContrasts={vi.fn()} onClearContrasts={vi.fn()}
    onParameterChange={change} onReset={vi.fn()} />);
  const policy = screen.getByRole("note", {name: "AS read-support policy"});
  expect(policy).toHaveTextContent("IJC + SJC must be at least 10 in every selected sample across both groups.");
  expect(policy).toHaveTextContent("Missing counts fail this check.");
  expect(policy).toHaveTextContent("FDR and |dPSI| thresholds are applied separately.");
  expect(screen.getByText("1 of 1 included")).toBeInTheDocument();
  const fdr = screen.getByRole("spinbutton", {name: /rMATS FDR/});
  const effect = screen.getByRole("spinbutton", {name: /Minimum \|dPSI\|/});
  expect(fdr).toHaveValue(0.05);
  expect(effect).toHaveValue(0.1);
  fireEvent.change(fdr, {target: {value: "0.01"}});
  fireEvent.change(effect, {target: {value: "0.2"}});
  expect(change).toHaveBeenCalledWith("fdr", 0.01);
  expect(change).toHaveBeenCalledWith("min_abs_dpsi", 0.2);
});

it("does not restore retired QC preset, replicate, direction or SE controls", () => {
  const retired: WorkflowParameter[] = [
    {key: "stability_profile", label: "Stable-AS profile", kind: "select", default: "balanced", help: "Retired profile", options: ["exploratory", "balanced", "stringent", "custom"]},
    {key: "min_reads", label: "Custom reads per replicate", kind: "number", default: 10, help: "Retired read setting", options: []},
    {key: "min_supported_replicates", label: "Custom supported replicates", kind: "number", default: 2, help: "Retired replicate setting", options: []},
    {key: "min_direction_consistency", label: "Custom direction consistency", kind: "number", default: 0.75, help: "Retired consistency setting", options: []},
    {key: "apply_measurement_se_gate", label: "Experimental SE gate", kind: "boolean", default: false, help: "Retired gate", options: []},
    {key: "max_measurement_se", label: "Custom maximum PSI SE", kind: "number", default: 0.12, help: "Retired uncertainty setting", options: []},
  ];
  render(<WorkflowControls workflow={{...readSupportWorkflow, parameters: [...readSupportWorkflow.parameters, ...retired]}}
    contrasts={[]} selectedContrasts={new Set()} parameters={{stability_profile: "custom", min_reads: 20, min_direction_consistency: 0.8, apply_measurement_se_gate: true}}
    selections={[]} conditions={[]} onToggleContrast={vi.fn()} onSelectAllContrasts={vi.fn()} onClearContrasts={vi.fn()}
    onParameterChange={vi.fn()} onReset={vi.fn()} />);
  expect(screen.getAllByRole("spinbutton")).toHaveLength(2);
  expect(screen.queryByRole("combobox")).not.toBeInTheDocument();
  expect(screen.queryByRole("switch")).not.toBeInTheDocument();
  for (const parameter of retired) expect(screen.queryByText(parameter.label)).not.toBeInTheDocument();
  expect(screen.queryByText(/Choose Custom|75% direction|2 supporting replicates/)).not.toBeInTheDocument();
  expect(screen.getByRole("note", {name: "AS read-support policy"})).toHaveTextContent("at least 10 in every selected sample");
});

it("shows the saved read-support result identity rather than a retired preset name", () => {
  const workflow: WorkflowDefinition = {...readSupportWorkflow, workflow_id: "differential.as", workspace_id: "differential", parameters: [
    {key: "quality_artifact_id", label: "Read support from Quality", kind: "artifact", default: "", help: "Completed compatible result", options: ["quality.stable_as", "AS read-support event QC"]},
  ]};
  const artifact: Artifact = {
    artifact_id: "qc-result", project_id: "toy", analysis_id: "quality.stable_as", name: "AS read-support event QC", kind: "table", entity_type: "splice_event", status: "complete",
    fingerprint: "content", run_fingerprint: "1234567890", path: "results/qc.parquet", row_count: 1,
    parent_artifact_ids: [], contrast_ids: [], parameters: {stability_profile: "balanced"}, metadata: {}, warnings: [], created_at: "2026-10-22T12:00:00Z",
  };
  render(<WorkflowControls workflow={workflow} contrasts={[]} selectedContrasts={new Set()} parameters={{}} artifacts={[artifact]}
    selections={[]} conditions={[]} onToggleContrast={vi.fn()} onSelectAllContrasts={vi.fn()} onClearContrasts={vi.fn()}
    onParameterChange={vi.fn()} onReset={vi.fn()} />);
  expect(screen.getByRole("option", {name: /AS read-support event QC/})).toBeInTheDocument();
  expect(screen.queryByRole("option", {name: /balanced|75%|supporting replicates/})).not.toBeInTheDocument();
  expect(screen.queryByRole("note", {name: "AS read-support policy"})).not.toBeInTheDocument();
});

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
