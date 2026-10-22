import { displayParameterOption } from "../branding";
import { DEFAULT_DAS_EVENT_SET } from "../dasSelection";
import type { WorkflowDefinition } from "../types";

export function DasSelectionNote({workflow, parameters}: {
  workflow: WorkflowDefinition; parameters: Record<string, unknown>;
}) {
  const policy = workflow.das_selection;
  const fixedResult = ["reports.interpreted_bundle", "prioritize.rbp_sequence_results", "prioritize.motiftrace_parallel", "prioritize.motiftrace_usage"].includes(workflow.workflow_id);
  if (!policy || (policy.mode === "not_applicable" && !fixedResult)) return null;
  return <div className="run-plan-panel" role="note" aria-label="DAS event selection policy">
    <strong>{policy.label}</strong>
    <p>{policy.description}</p>
    {fixedResult && <p><strong>Completed result scope.</strong> This view retains its saved input population and original statistics. It is not recalculated by the ordinary DAS / RS-DAS switch.</p>}
    {policy.mode === "switchable" && <p>
      Current selection: <strong>{displayParameterOption("event_set", String(parameters.event_set ?? policy.default_event_set ?? DEFAULT_DAS_EVENT_SET))}</strong>.
      This selection is shared with other switchable downstream modules. Changed inputs require their own matching downstream result; existing results are not relabelled.
    </p>}
    {policy.mode === "rs_required" && <p><strong>RS-DAS required.</strong> Switching other modules to ordinary DAS does not change this module's eligible input or reuse a result under a different label.</p>}
    {policy.mode === "measurement_population" && <p><strong>Independent measurement scope.</strong> The significant DAS / RS-DAS foreground does not replace this module's declared measured-input population. Its own measurement and read-support requirements remain in effect.</p>}
  </div>;
}
