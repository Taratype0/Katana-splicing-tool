import { render, screen } from "@testing-library/react";
import { describe, expect, it } from "vitest";
import { ContextBar } from "./ContextBar";
import type { ProjectPayload } from "../types";

describe("ContextBar", () => {
  it("does not display all comparisons when none are selected", () => {
    const project = { snapshot: { name: "Long project name", contrast_count: 4,
      output_root: "E:/results", genome_build: "GRCm39", manifest_fingerprint: "abcdefgh123456" } } as ProjectPayload;
    render(<ContextBar project={project} selectedContrastCount={0} workflowLabel="Candidate evidence" />);
    expect(screen.getByText("0 selected")).toBeTruthy();
    expect(screen.getByTitle("Long project name")).toBeTruthy();
    expect(screen.getByTitle("Candidate evidence")).toBeTruthy();
  });
});
