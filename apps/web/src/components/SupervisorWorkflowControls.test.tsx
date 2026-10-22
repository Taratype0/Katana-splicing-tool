import { fireEvent, render, screen } from "@testing-library/react";
import { expect, it, vi } from "vitest";
import type { WorkflowDefinition, WorkflowParameter } from "../types";
import { WorkflowControls } from "./WorkflowControls";

function controls(parameters: WorkflowParameter[], workflowId: string) {
  const change = vi.fn();
  const workflow: WorkflowDefinition = {
    workflow_id: workflowId, workspace_id: "prioritize", label: "Workflow",
    question: "What is measured?", method: "Declared method", output_summary: "Auditable outputs",
    analysis_unit: "Declared unit", expected_artifacts: 1, parameters,
    required_asset_kinds: [], optional_asset_kinds: [], supports_contrast_selection: false,
    enabled: true, limitations: "No pooling or silent substitution",
  };
  render(<WorkflowControls workflow={workflow} contrasts={[]} selectedContrasts={new Set()}
    parameters={{}} selections={[]} conditions={[]} onToggleContrast={vi.fn()}
    onSelectAllContrasts={vi.fn()} onClearContrasts={vi.fn()} onParameterChange={change} onReset={vi.fn()}/>);
  return change;
}

it("allows GO namespace selection without changing the statistical settings", () => {
  const change = controls([{key:"ontology_category",label:"Show GO / pathway category",kind:"select",
    default:"All",help:"Display only; all P/FDR retained",options:["All","GO BP","GO MF","GO CC","KEGG","Unclassified"]}], "prioritize.enrichment");
  fireEvent.change(screen.getByRole("combobox"), {target:{value:"GO MF"}});
  expect(change).toHaveBeenCalledExactlyOnceWith("ontology_category", "GO MF");
  expect(screen.getByText("Display only; all P/FDR retained")).toBeInTheDocument();
});

it("shows editable settings and puts methods in a separate disclosure", () => {
  controls([{key:"fdr",label:"FDR threshold",kind:"number",default:.05,help:"Adjusted P-value cutoff",options:[],minimum:0,maximum:1,step:.01}], "prioritize.enrichment");
  expect(screen.getByText("Analysis settings")).toBeInTheDocument();
  expect(screen.getByRole("spinbutton")).toHaveValue(.05);
  expect(screen.getByText("Statistical method and interpretation")).toBeInTheDocument();
  expect(screen.queryByText("Question")).not.toBeInTheDocument();
});

it("exposes cell feature source, read support and donor evaluation independently", () => {
  const change = controls([
    {key:"mode",label:"Feature source",kind:"select",default:"splicing",help:"Counts, not expression substitution",options:["splicing","expression"]},
    {key:"min_reads",label:"Reads per cell-event",kind:"number",default:10,help:"Unmeasured PSI is missing, not zero",options:[],minimum:1,step:1},
    {key:"evaluate_cell_types",label:"Evaluate held-out donors",kind:"boolean",default:true,help:"Cells are not biological replicates",options:[]},
  ], "single_cell.cell_embedding");
  fireEvent.change(screen.getByRole("combobox"), {target:{value:"expression"}});
  fireEvent.change(screen.getByRole("spinbutton"), {target:{value:"15"}});
  fireEvent.click(screen.getByRole("switch"));
  expect(change.mock.calls).toEqual([["mode","expression"],["min_reads",15],["evaluate_cell_types",false]]);
  expect(screen.getByText("Cells are not biological replicates")).toBeInTheDocument();
});
