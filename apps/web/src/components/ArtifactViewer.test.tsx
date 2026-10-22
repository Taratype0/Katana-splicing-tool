import { useState } from "react";
import { fireEvent, render, screen } from "@testing-library/react";
import { describe, expect, it, vi } from "vitest";
import { ArtifactViewer } from "./ArtifactViewer";
import type { Artifact } from "../types";

vi.mock("./ArtifactTable", () => ({
  ArtifactTable: () => {
    const [value, setValue] = useState("");
    return <input aria-label="table state" value={value} onChange={event => setValue(event.target.value)} />;
  },
}));
vi.mock("./PlotViewer", () => ({ PlotViewer: () => null }));
vi.mock("./ImagePreview", () => ({ ImagePreview: () => null }));

describe("artifact table identity", () => {
  it("retains state for the same result and resets it across tables", () => {
    const first = { artifact_id: "first", kind: "table" } as Artifact;
    const props = { projectId: "project", selectedIds: new Set<string>(), onToggle: vi.fn(), onFocus: vi.fn() };
    const { rerender } = render(<ArtifactViewer {...props} artifact={first} />);
    fireEvent.change(screen.getByLabelText("table state"), { target: { value: "old sort/search" } });
    rerender(<ArtifactViewer {...props} artifact={first} />);
    expect(screen.getByLabelText("table state")).toHaveValue("old sort/search");
    rerender(<ArtifactViewer {...props} artifact={{ ...first, artifact_id: "second" }} />);
    expect(screen.getByLabelText("table state")).toHaveValue("");
  });
});
