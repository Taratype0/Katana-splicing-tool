import { fireEvent, render, screen } from "@testing-library/react";
import { beforeEach, expect, it, vi } from "vitest";
import { PathPicker } from "./PathPicker";
import { api } from "../api";

vi.mock("../api",()=>({api:{browseLocalFiles:vi.fn()}}));
beforeEach(()=>{vi.clearAllMocks();localStorage.clear();HTMLDialogElement.prototype.showModal = function(){this.setAttribute("open", "");};});

it("navigates and chooses a folder without entering a path; cancel preserves selection",async()=>{
  vi.mocked(api.browseLocalFiles).mockImplementation(async path=>({path:path||"/home",parent:"/",locations:[{name:"Home",path:"/home"}],
    entries:path==="/home/report"?[]:[{name:"report",path:"/home/report",directory:true}]}));
  const change=vi.fn(); render(<PathPicker label="Report" value="" onChange={change}/>);
  fireEvent.click(screen.getByRole("button",{name:"Browse Report"}));
  fireEvent.click(await screen.findByRole("button",{name:"report"}));
  await screen.findByText("/home/report");
  expect(change).not.toHaveBeenCalled();
  fireEvent.click(screen.getByRole("button",{name:"Choose folder"}));
  expect(change).toHaveBeenCalledWith("/home/report");
  fireEvent.click(screen.getByRole("button",{name:"Browse Report"}));
  fireEvent.click(await screen.findByRole("button",{name:"Cancel"}));
  expect(change).toHaveBeenCalledTimes(1);
});

it("file mode requires a selected file and filters unsupported extensions",async()=>{
  vi.mocked(api.browseLocalFiles).mockResolvedValue({path:"/home",parent:"/",locations:[],entries:[
    {name:"report-manifest.json",path:"/home/report-manifest.json",directory:false},
    {name:"notes.txt",path:"/home/notes.txt",directory:false}]});
  const change=vi.fn();render(<PathPicker label="Manifest" value="" onChange={change} kind="file" extensions={[".json"]}/>);
  fireEvent.click(screen.getByRole("button",{name:"Browse Manifest"}));
  const file=await screen.findByRole("button",{name:"report-manifest.json"});
  expect(screen.queryByRole("button",{name:"notes.txt"})).not.toBeInTheDocument();
  expect(screen.getByRole("button",{name:"Choose file"})).toBeDisabled();
  fireEvent.click(file);fireEvent.click(screen.getByRole("button",{name:"Choose file"}));
  expect(change).toHaveBeenCalledWith("/home/report-manifest.json");
});
