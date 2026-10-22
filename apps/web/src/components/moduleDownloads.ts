import type { ModuleDashboard } from "../api";

export function tableTSV(rows: Array<Record<string, unknown>>, columns: string[]) {
  const cell = (value: unknown) => {
    const text = value == null ? "" : String(value);
    return /[\t\r\n"]/.test(text) ? `"${text.replaceAll('"', '""')}"` : text;
  };
  return [columns.map(cell).join("\t"), ...rows.map(row => columns.map(key => cell(row[key])).join("\t"))].join("\n") + "\n";
}

export function saveBlob(blob: Blob, name: string) {
  const url = URL.createObjectURL(blob);
  const link = document.createElement("a");
  link.href = url; link.download = name;
  document.body.appendChild(link); link.click(); link.remove();
  window.setTimeout(() => URL.revokeObjectURL(url), 1000);
}

export function moduleFileStem(data: ModuleDashboard) {
  const safe = (value: string) => value.replace(/[^\p{L}\p{N}_.-]+/gu, "_").slice(0, 80) || "scope";
  return `cosplice-${safe(data.module)}-${safe(data.condition)}`;
}

export function networkSVG(element: SVGSVGElement, data: ModuleDashboard, hideSameGene = false) {
  const ns = "http://www.w3.org/2000/svg";
  const outer = document.createElementNS(ns, "svg");
  outer.setAttribute("width", "1140"); outer.setAttribute("height", "1200");
  outer.setAttribute("viewBox", "0 0 380 400"); outer.setAttribute("font-family", "Arial, sans-serif");
  const title = document.createElementNS(ns, "title");
  title.textContent = `cosplice event network: ${data.module}, ${data.condition}`; outer.appendChild(title);
  const description = document.createElementNS(ns, "desc");
  description.textContent = JSON.stringify({module:data.module,condition:data.condition,discovery_run:data.discovery_run,
    hub_analysis_run:data.hub_analysis_run,selected_event:data.selected_event,hide_same_gene_in_drawing:hideSameGene,
    display_rule:data.display_rule,method:data.method,interpretation:data.interpretation});
  outer.appendChild(description);
  const background = document.createElementNS(ns, "rect");
  background.setAttribute("width", "380"); background.setAttribute("height", "400"); background.setAttribute("fill", "white"); outer.appendChild(background);
  for (const [y, text, size] of [[17, `Module: ${data.module} | Scope: ${data.condition}`, 11],
    [32, `${data.nodes.length}/${data.total_events} events shown; ${data.edges.length} edges shown`, 9],
    [389, "PSI associations; drawing is a subset. Hub metrics use the full module.", 8]] as const) {
    const label = document.createElementNS(ns, "text"); label.setAttribute("x", "10"); label.setAttribute("y", String(y));
    label.setAttribute("font-size", String(size)); label.setAttribute("fill", "#243f4b"); label.textContent = text; outer.appendChild(label);
  }
  const graph = element.cloneNode(true) as SVGSVGElement;
  graph.setAttribute("x", "0"); graph.setAttribute("y", "39"); graph.setAttribute("width", "380"); graph.setAttribute("height", "340");
  outer.appendChild(graph);
  return new XMLSerializer().serializeToString(outer);
}

export async function saveNetwork(element: SVGSVGElement, data: ModuleDashboard, format: "svg" | "png", hideSameGene = false) {
  const blob = new Blob([networkSVG(element, data, hideSameGene)], {type:"image/svg+xml;charset=utf-8"});
  const filename = `${moduleFileStem(data)}-${hideSameGene?"cross-gene-view-":""}network.${format}`;
  if (format === "svg") { saveBlob(blob, filename); return; }
  const url = URL.createObjectURL(blob);
  try {
    const image = new Image();
    await new Promise<void>((resolve, reject) => {image.onload=()=>resolve(); image.onerror=()=>reject(new Error("Network image could not be rendered")); image.src=url;});
    const canvas = document.createElement("canvas"); canvas.width=1140; canvas.height=1200;
    const context = canvas.getContext("2d"); if (!context) throw new Error("PNG export is unavailable in this browser; use SVG");
    context.drawImage(image, 0, 0, canvas.width, canvas.height);
    const png = await new Promise<Blob>((resolve, reject)=>canvas.toBlob(value=>value?resolve(value):reject(new Error("PNG export failed; use SVG")), "image/png"));
    saveBlob(png, filename);
  } finally { URL.revokeObjectURL(url); }
}
