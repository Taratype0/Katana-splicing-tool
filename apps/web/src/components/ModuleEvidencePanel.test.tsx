import { render, screen } from "@testing-library/react";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { expect, it, vi } from "vitest";
import { api } from "../api";
import { ModuleEvidencePanel } from "./ModuleEvidencePanel";

it("reads bounded module evidence without launching analysis", async () => {
  const read = vi.spyOn(api,"moduleEvidence").mockResolvedValue({module:"M1",discovery_run:"frozen",scope:"Cached evidence only",tables:[
    {artifact_id:"events",name:"Module members",entity_type:"event_diagnostic",row_count:100,preview_limit:30,
      rows:[{module:"M1",event_key:"event1"}],run_fingerprint:"frozen"},
  ]});
  const run = vi.spyOn(api,"runWorkflow");
  render(<QueryClientProvider client={new QueryClient()}><ModuleEvidencePanel projectId="p" artifactId="a" module="M1" /></QueryClientProvider>);
  expect(await screen.findByText("Module members (100 rows)")).toBeInTheDocument();
  expect(screen.getByText(/first 30 matching records/)).toBeInTheDocument();
  expect(read).toHaveBeenCalledWith("p","a","M1");
  expect(run).not.toHaveBeenCalled();
  vi.restoreAllMocks();
});
