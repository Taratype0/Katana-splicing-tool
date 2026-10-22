import {render, screen} from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import {describe, expect, it, vi} from "vitest";
import {GeneFilterPanel} from "./GeneFilterPanel";
import type {ProjectPayload} from "../types";

const project = {manifest:{metadata:{}}} as ProjectPayload;
describe("Gene exclusion",()=>{
  it("applies annotation categories and a deduplicated custom list only on submit",async()=>{
    const user=userEvent.setup(), onApply=vi.fn();
    render(<GeneFilterPanel project={project} busy={false} onApply={onApply}/>);
    await user.click(screen.getByText(/Gene exclusion/));
    await user.click(screen.getByLabelText(/Exclude pseudogenes/));
    await user.type(screen.getByLabelText(/Additional gene IDs/),"GeneA, GeneB\nGeneA");
    expect(onApply).not.toHaveBeenCalled();
    await user.click(screen.getByRole("button",{name:/Apply exclusions/}));
    expect(onApply).toHaveBeenCalledWith(expect.objectContaining({exclude_pseudogenes:true,exclude_mitochondrial:false,excluded_genes:["GeneA","GeneB"]}));
  });
  it("restores saved exclusions and shows actionable API errors",()=>{
    const filtered={manifest:{metadata:{gene_filter:{config:{exclude_mitochondrial:true}}}}} as unknown as ProjectPayload;
    render(<GeneFilterPanel project={filtered} busy={false} error="Chromosome annotation missing" onApply={()=>{}}/>);
    expect(screen.getByLabelText(/Exclude mitochondrial/)).toBeChecked();
    expect(screen.getByRole("alert")).toHaveTextContent("Chromosome annotation missing");
  });
});
