import { render, screen } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import { describe, expect, it, vi } from "vitest";
import { ProjectGate } from "./ProjectGate";

describe("ProjectGate", () => {
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
