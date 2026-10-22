import { expect, it } from "vitest";
import type { ModuleDashboard } from "../api";
import { networkSVG, tableTSV } from "./moduleDownloads";

it("keeps identifiers, numeric zero, missing values and multiline reasons in TSV exports",()=>{
  const tsv=tableTSV([{event_key:'chr1:12-34',kME:0,reason:'missing\tgene\n"annotation"',unknown:null}],['event_key','kME','reason','unknown']);
  expect(tsv).toBe('event_key\tkME\treason\tunknown\nchr1:12-34\t0\t"missing\tgene\n""annotation"""\t\n');
  expect(tableTSV([],['event_key','condition'])).toBe('event_key\tcondition\n');
});

it("exports a standalone labelled SVG with the current scope and unchanged visible graph",()=>{
  const svg=document.createElementNS('http://www.w3.org/2000/svg','svg');
  svg.setAttribute('viewBox','0 0 380 340');
  const circle=document.createElementNS(svg.namespaceURI,'circle');circle.setAttribute('cx','42');svg.appendChild(circle);
  const data={module:'M1',condition:'A & B',discovery_run:'frozen-run',nodes:[{}],edges:[],total_events:80,display_rule:'Shown subset',method:'Observed PSI',interpretation:'Association only'}as unknown as ModuleDashboard;
  const output=networkSVG(svg,data);
  const parsed=new DOMParser().parseFromString(output,'image/svg+xml');
  expect(parsed.querySelector('parsererror')).toBeNull();
  expect(parsed.querySelector('title')?.textContent).toContain('M1, A & B');
  expect(parsed.querySelector('desc')?.textContent).toContain('frozen-run');
  expect(parsed.querySelector('circle')?.getAttribute('cx')).toBe('42');
  expect(parsed.documentElement.textContent).toContain('1/80 events shown');
  expect(svg.hasAttribute('width')).toBe(false);
});
