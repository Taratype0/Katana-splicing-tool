import { render, screen, fireEvent } from "@testing-library/react";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { expect, it, vi } from "vitest";
import { api } from "../api";
import { RunPlanPanel } from "./RunPlanPanel";

it("checks prerequisites only on request and never submits an analysis", async () => {
  const plan = vi.spyOn(api,"workflowPlan").mockResolvedValue({ready:false,errors:["Missing source"],warnings:[],scope:"read-only",input_bytes:0,memory_estimate:"unknown",parent_artifacts:[]});
  const run = vi.spyOn(api,"runWorkflow");
  render(<QueryClientProvider client={new QueryClient()}><RunPlanPanel projectId="p" workflowId="programs.diagnostics" parameters={{}} /></QueryClientProvider>);
  expect(plan).not.toHaveBeenCalled();
  fireEvent.click(screen.getByText("Check prerequisites"));
  expect(await screen.findByText("Missing source")).toBeInTheDocument();
  expect(run).not.toHaveBeenCalled();
  vi.restoreAllMocks();
});
