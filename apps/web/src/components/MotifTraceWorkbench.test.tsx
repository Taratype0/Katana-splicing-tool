import { fireEvent, render, screen } from "@testing-library/react";
import { expect, it, vi } from "vitest";
import { api, type MotifTraceCatalog } from "../api";
import type { Artifact, WorkflowDefinition } from "../types";
import { MotifTraceWorkbench, motifTraceModeParameters } from "./MotifTraceWorkbench";

const catalog: MotifTraceCatalog = {package_version:"1.0.0.post1",scientific_release:"1.0.0-frozen.20261017",modes:[
  {mode_id:"boundary",label:"Boundary associations",workflow_id:"toy.boundary",fixed_parameters:{sequence_scope:"boundary_association"},description:"Matched splice-boundary analysis",operation:"new_analysis",requires_selection:true},
  {mode_id:"se_impact",label:"SE impact",workflow_id:"toy.se",fixed_parameters:{},description:"SE RNA sequences",operation:"new_analysis",requires_selection:true},
  {mode_id:"completed_maps",label:"Completed sequence maps",workflow_id:"toy.maps",fixed_parameters:{},description:"Read completed RNA and motif positions",operation:"completed_result_import",requires_selection:false},
  {mode_id:"consequences",label:"Event consequences",workflow_id:"toy.consequences",fixed_parameters:{},description:"Event-level sequence consequences",operation:"new_analysis",requires_selection:true},
  {mode_id:"usage",label:"Transcript usage",workflow_id:"toy.usage",fixed_parameters:{},description:"Transcript-usage consequences",operation:"existing_result_integration",requires_selection:false},
  {mode_id:"parallel",label:"Completed MAPP / MotifTrace",workflow_id:"toy.parallel",fixed_parameters:{},description:"Original MAPP fits with separate MotifTrace evidence",operation:"completed_result_import",requires_selection:false},
]};
const workflows=catalog.modes.map(mode=>({workflow_id:mode.workflow_id,enabled:true,parameters:[]}) as unknown as WorkflowDefinition);
const browser={artifact_id:"saved-browser",project_id:"p",analysis_id:"toy.maps",name:"Sequence browser",status:"complete",path:"results/browser.html",parameters:{},run_fingerprint:"exact-source-run"} as Artifact;

it("exposes all server-catalog modes and selects without running or modifying a saved result",()=>{
  const select=vi.fn(), run=vi.spyOn(api,"runWorkflow"), plan=vi.spyOn(api,"workflowPlan");
  render(<MotifTraceWorkbench projectId="p" catalog={catalog} workflows={workflows} artifacts={[]}
    activeWorkflowId="toy.boundary" parameters={{sequence_scope:"boundary_association",selection_id:"saved-selection"}} onSelect={select}/>);
  const choice=screen.getByRole("combobox",{name:"MotifTrace mode"});
  expect(choice).toHaveValue("boundary");
  for(const mode of catalog.modes) expect(screen.getByRole("option",{name:mode.label})).toBeInTheDocument();
  fireEvent.change(choice,{target:{value:"parallel"}});
  expect(select).toHaveBeenCalledWith(catalog.modes[5]);
  expect(run).not.toHaveBeenCalled();expect(plan).not.toHaveBeenCalled();
  expect(screen.getByText(/1.0.0-frozen.20261017/)).toBeInTheDocument();vi.restoreAllMocks();
});

it("opens an existing completed browser and keeps completed-result mode distinct from new analysis",()=>{
  const run=vi.spyOn(api,"runWorkflow");
  render(<MotifTraceWorkbench projectId="p" catalog={catalog} workflows={workflows} artifacts={[browser,{...browser,artifact_id:"other",analysis_id:"toy.boundary"}]}
    activeWorkflowId="toy.maps" parameters={{}} onSelect={vi.fn()}/>);
  expect(screen.getByText("Completed-result mode.")).toBeInTheDocument();
  expect(screen.getByText(/no model fitting or upstream quantification is started/)).toBeInTheDocument();
  expect(screen.getByRole("link",{name:"Open Sequence browser"})).toHaveAttribute("href","/api/projects/p/artifacts/saved-browser/file");
  expect(screen.getByText(/Saved outputs for this mode: 1/)).toBeInTheDocument();
  expect(run).not.toHaveBeenCalled();vi.restoreAllMocks();
});

it("retains shared current settings and source selection while applying only catalog mode patch",()=>{
  const target={...workflows[0],parameters:[
    {key:"selection_id",default:""},{key:"custom_motif_path",default:""},{key:"sequence_scope",default:"legacy"},
  ]} as WorkflowDefinition;
  const current={selection_id:"current-selection",custom_motif_path:"user-models.meme",sequence_scope:"other",unrelated_threshold:999};
  const parameters=motifTraceModeParameters(catalog.modes[0],target,current,{selection_id:"old-selection"});
  expect(parameters).toEqual({selection_id:"current-selection",custom_motif_path:"user-models.meme",sequence_scope:"boundary_association"});
  expect(current.selection_id).toBe("current-selection");
});

it("reports catalog failure and unavailable modes without inventing alternate thresholds",()=>{
  render(<MotifTraceWorkbench projectId="p" error="Catalog unavailable" workflows={workflows} artifacts={[]}
    activeWorkflowId="toy.maps" parameters={{}} onSelect={vi.fn()}/>);
  expect(screen.getByRole("alert")).toHaveTextContent("Catalog unavailable");
  expect(screen.getByRole("combobox",{name:"MotifTrace mode"})).toBeDisabled();
});

it("labels usage as integration of existing results rather than upstream quantification",()=>{
  render(<MotifTraceWorkbench projectId="p" catalog={catalog} workflows={workflows} artifacts={[]}
    activeWorkflowId="toy.usage" parameters={{}} onSelect={vi.fn()}/>);
  expect(screen.getByText("Existing-result integration.")).toBeInTheDocument();
  expect(screen.getByText(/No upstream quantification is started/)).toBeInTheDocument();
  expect(screen.queryByText("New analysis mode.")).not.toBeInTheDocument();
});

it("does not keep the boundary mode label or reset settings after manual scope changes",()=>{
  const select=vi.fn();
  const props={projectId:"p",catalog,workflows,artifacts:[],activeWorkflowId:"toy.boundary",onSelect:select};
  const view=render(<MotifTraceWorkbench {...props} parameters={{sequence_scope:"boundary_association"}}/>);
  expect(screen.getByRole("combobox",{name:"MotifTrace mode"})).toHaveValue("boundary");
  view.rerender(<MotifTraceWorkbench {...props} parameters={{sequence_scope:"completed_sequence_maps"}}/>);
  expect(screen.getByRole("combobox",{name:"MotifTrace mode"})).toHaveValue("");
  expect(screen.getByRole("option",{name:"Choose a MotifTrace mode"})).toHaveProperty("selected",true);
  expect(select).not.toHaveBeenCalled();
});
