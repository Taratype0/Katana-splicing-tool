import { api, type MotifTraceCatalog, type MotifTraceMode } from "../api";
import { displayArtifactName } from "../branding";
import type { Artifact, WorkflowDefinition } from "../types";

type Props = {
  projectId: string; catalog?: MotifTraceCatalog; loading?: boolean; error?: string;
  workflows: WorkflowDefinition[]; artifacts: Artifact[];
  activeWorkflowId: string; parameters: Record<string, unknown>;
  onSelect: (mode: MotifTraceMode) => void;
};

export function motifTraceModeParameters(mode: MotifTraceMode, target: WorkflowDefinition,
  current: Record<string, unknown>, saved: Record<string, unknown> = {}) {
  const defaults=Object.fromEntries(target.parameters.map(parameter=>[parameter.key,parameter.default]));
  const shared=Object.fromEntries(target.parameters.filter(parameter=>parameter.key in current)
    .map(parameter=>[parameter.key,current[parameter.key]]));
  return {...defaults,...saved,...shared,...mode.fixed_parameters};
}

export function MotifTraceWorkbench({projectId,catalog,loading,error,workflows,artifacts,activeWorkflowId,parameters,onSelect}: Props) {
  const modes=catalog?.modes ?? [];
  const active=modes.find(mode=>mode.workflow_id===activeWorkflowId && Object.entries(mode.fixed_parameters).every(([key,value])=>JSON.stringify(parameters[key])===JSON.stringify(value)));
  const outputs=active ? artifacts.filter(artifact=>artifact.status==="complete" && artifact.analysis_id===active.workflow_id
    && Object.entries(active.fixed_parameters).every(([key,value])=>JSON.stringify(artifact.parameters[key])===JSON.stringify(value))) : [];
  const browsers=outputs.filter(artifact=>/\.html?$/i.test(artifact.path));
  return <section className="run-plan-panel" aria-label="Katana MotifTrace workbench">
    <h3>Katana MotifTrace</h3>
    <p>Choose an analysis or open completed sequence results. All modes use the same project, comparison controls and analysis engine.</p>
    {loading && <p role="status">Loading MotifTrace modes...</p>}
    {error && <p role="alert">Could not load MotifTrace modes: {error}</p>}
    <label className="parameter-field"><span>MotifTrace mode</span>
      <select aria-label="MotifTrace mode" value={active?.mode_id ?? ""} disabled={!modes.length}
        onChange={event=>{const mode=modes.find(item=>item.mode_id===event.target.value);if(mode)onSelect(mode);}}>
        <option value="" disabled>Choose a MotifTrace mode</option>
        {modes.map(mode=><option key={mode.mode_id} value={mode.mode_id} disabled={!workflows.some(workflow=>workflow.workflow_id===mode.workflow_id && workflow.enabled)}>{mode.label}</option>)}
      </select>
    </label>
    {active && <>
      <p>{active.description}</p>
      <p><strong>{active.operation==="completed_result_import"?"Completed-result mode":active.operation==="existing_result_integration"?"Existing-result integration":"New analysis mode"}.</strong> {active.operation==="completed_result_import"?"Reuse an explicitly supplied completed result; no model fitting or upstream quantification is started.":active.operation==="existing_result_integration"?"Integrate explicitly supplied existing results. No upstream quantification is started.":"Check the required inputs below before running this analysis."} {active.requires_selection?"A compatible saved selection is required.":"Review this mode's input and comparison requirements below."}</p>
      <p>Saved outputs for this mode: {outputs.length}. These are output files, not event or model counts.</p>
      {browsers.length>0 ? <ul aria-label="Completed MotifTrace browsers">{browsers.map(artifact=><li key={artifact.artifact_id}>
        <a href={api.artifactFileUrl(projectId,artifact.artifact_id)} target="_blank" rel="noreferrer">Open {displayArtifactName(artifact.name)}</a>
        <small> · source run {artifact.run_fingerprint.slice(0,12)} · {artifact.contrast_ids?.length?artifact.contrast_ids.join(", "):"scope recorded in browser"}</small>
      </li>)}</ul> : <p>No completed browser is listed for this mode. Missing output is not a negative biological result.</p>}
    </>}
    <small>Scientific release: {catalog?.scientific_release ?? "not loaded"}. Selecting a mode or opening a saved browser does not run an analysis. Use Check prerequisites and the existing Run analysis button.</small>
  </section>;
}
