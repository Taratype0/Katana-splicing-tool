import { render, screen } from "@testing-library/react";
import { expect, it } from "vitest";
import { WorkspaceNav } from "./WorkspaceNav";

it("uses the original Katana artwork rather than a placeholder mark", () => {
  render(<WorkspaceNav workspaces={[]} activeId="project" onChange={() => {}} />);
  expect(screen.getByRole("img", { name: "Katana logo" })).toHaveAttribute("src", "/katana_icon.png");
});
