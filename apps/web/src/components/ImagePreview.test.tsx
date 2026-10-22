import { render, screen, fireEvent } from "@testing-library/react";
import { expect, it } from "vitest";
import { ImagePreview } from "./ImagePreview";
import type { Artifact } from "../types";

it("offers readable zoom modes and both cached figure exports", () => {
  const artifact = { artifact_id: "image-1", name: "Donor Sashimi", metadata: { exports: { png: "a.png", svg: "a.svg" } } } as unknown as Artifact;
  render(<ImagePreview projectId="p1" artifact={artifact} />);
  expect(screen.getByRole("button", {name: "Fit page"})).toHaveAttribute("aria-pressed", "true");
  fireEvent.click(screen.getByRole("button", {name: "Fit width"}));
  expect(screen.getByLabelText("Image viewport")).toHaveClass("image-stage-width");
  fireEvent.click(screen.getByRole("button", {name: "100%"}));
  expect(screen.getByLabelText("Image viewport")).toHaveClass("image-stage-actual");
  expect(screen.getByRole("link", {name: "SVG"})).toHaveAttribute("href", "/api/projects/p1/artifacts/image-1/export/svg");
});
