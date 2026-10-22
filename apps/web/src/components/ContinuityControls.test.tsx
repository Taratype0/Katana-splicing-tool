import { fireEvent, render, screen } from "@testing-library/react";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { afterEach, expect, it, vi } from "vitest";
import { api } from "../api";
import { JobDrawer } from "./JobDrawer";
import { LineagePanel } from "./LineagePanel";
import type { Job } from "../types";

afterEach(() => vi.restoreAllMocks());

it("does not retry persisted failures until explicitly clicked", async () => {
  const failed = {job_id:"j", analysis_id:"quality.stable_as", status:"failed", message:"Interrupted", progress:.4, cache_status:"rebuild"} as Job;
  const retry = vi.spyOn(api,"retryJob").mockResolvedValue({...failed,status:"queued"});
  render(<QueryClientProvider client={new QueryClient()}><JobDrawer jobs={[failed]}/></QueryClientProvider>);
  expect(retry).not.toHaveBeenCalled();
  fireEvent.click(screen.getByText("Activity"));
  fireEvent.click(screen.getByText("Retry / resume saved run"));
  await vi.waitFor(() => expect(retry).toHaveBeenCalledWith("j", expect.anything()));
});

it("reads ancestry without launching analysis or turning missing counts into zero", async () => {
  const read = vi.spyOn(api,"lineage").mockResolvedValue({scope:"Do not subtract different units",nodes:[{artifact_id:"a",name:"Source",workflow:"quality.stable_as",rows:null,entity_type:"splice_event",state:"missing_file"}]});
  const run = vi.spyOn(api,"runWorkflow");
  render(<QueryClientProvider client={new QueryClient()}><LineagePanel projectId="p" artifactId="a"/></QueryClientProvider>);
  expect(read).not.toHaveBeenCalled();
  fireEvent.click(screen.getByText("Read recorded analysis path"));
  expect(await screen.findByText(/unavailable rows/)).toBeInTheDocument();
  expect(run).not.toHaveBeenCalled();
});
