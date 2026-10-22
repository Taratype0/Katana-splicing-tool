import { describe, expect, it } from "vitest";
import { displayModuleName } from "./branding";

describe("Katana-Splice-Network public naming", () => {
  it("renames public legacy labels", () => {
    expect(displayModuleName("CoSplice / Katana-CoSplice / Katana-Splice")).toBe("Katana-Splice-Network / Katana-Splice-Network / Katana-Splice-Network");
  });
  it("preserves durable identifiers and comparator names", () => {
    const original = "CoSplicEx CoSpliceNet programs.cosplice 15_CoSplice CoSpliceConfig";
    expect(displayModuleName(original)).toBe(original);
  });
  it("is idempotent and preserves paths and CLI aliases", () => {
    const title = "Katana-Splice-Network / CoSplice完整结果";
    expect(displayModuleName(displayModuleName(title))).toBe("Katana-Splice-Network / Katana-Splice-Network完整结果");
    const paths = String.raw`E:/results/Katana-Splice/CoSplice.tsv E:\results\CoSplice\tables katana-splice`;
    expect(displayModuleName(paths)).toBe(paths);
  });
});
