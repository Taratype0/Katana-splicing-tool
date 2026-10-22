import { render, screen } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import { describe, expect, it, vi } from "vitest";
import { ProjectGate } from "./ProjectGate";

describe("ProjectGate", () => {
  it("attaches a resource bundle only on explicit submit", async () => {
    const user = userEvent.setup();
    const onOpen = vi.fn();
    const onConfigureResources = vi.fn();
    render(<ProjectGate defaultPath="E:/Analysis/project.json" loading={false} onOpen={onOpen} onConfigureResources={onConfigureResources} />);
    await user.click(screen.getByRole("button", { name: "GO / MotifTrace resources" }));
    expect(onConfigureResources).not.toHaveBeenCalled();
    await user.type(screen.getByLabelText("Prepared resources.json path"), "E:/Annotations/resources.json");
    await user.type(screen.getByLabelText("New empty project directory"), "E:/Analysis/new");
    await user.click(screen.getByRole("button", { name: "Attach resources and open" }));
    expect(onConfigureResources).toHaveBeenCalledWith({ manifest_path: "E:/Analysis/project.json", resource_bundle_path: "E:/Annotations/resources.json", output_root: "E:/Analysis/new" });
    expect(onOpen).not.toHaveBeenCalled();
  });
  it("submits the explicit manifest path", async () => {
    const user = userEvent.setup();
    const onOpen = vi.fn();
    render(<ProjectGate defaultPath="" loading={false} onOpen={onOpen} />);

    await user.type(
      screen.getByLabelText("Project manifest path"),
      "E:\\Analysis\\katana.project.json",
    );
    await user.click(screen.getByRole("button", { name: "Open project" }));

    expect(onOpen).toHaveBeenCalledWith("E:\\Analysis\\katana.project.json");
  });

  it("shows API validation errors without hiding the input", () => {
    render(
      <ProjectGate
        defaultPath={"E:\\bad.json"}
        loading={false}
        error="Project manifest not found"
        onOpen={() => undefined}
      />,
    );

    expect(screen.getByRole("alert")).toHaveTextContent("Project manifest not found");
    expect(screen.getByLabelText("Project manifest path")).toHaveValue("E:\\bad.json");
  });
});
