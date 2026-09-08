import { render, screen } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import { describe, expect, it, vi } from "vitest";
import { WorkspaceNav } from "./WorkspaceNav";

describe("WorkspaceNav", () => {
  it("allows enabled workspaces and blocks planned ones", async () => {
    const user = userEvent.setup();
    const onChange = vi.fn();
    render(
      <WorkspaceNav
        activeId="project"
        onChange={onChange}
        workspaces={[
          { workspace_id: "project", label: "Project", purpose: "Inputs", icon: "folder", enabled: true, phase: 1 },
          { workspace_id: "programs", label: "Programs", purpose: "Modules", icon: "network", enabled: false, phase: 3 },
        ]}
      />,
    );

    await user.click(screen.getByRole("button", { name: /Project/ }));
    expect(onChange).toHaveBeenCalledWith("project");
    expect(screen.getByRole("button", { name: /Programs/ })).toBeDisabled();
  });
});
