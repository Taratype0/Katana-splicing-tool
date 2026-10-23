import { useEffect, useState } from "react";
import type { Artifact, TargetBasketEntry, TargetKind } from "./types";

export function browsingSnapshot(artifact: Artifact, settings: Record<string, unknown>, table: Record<string, unknown>) {
  return {
    source_parameters: structuredClone(artifact.parameters),
    source_fingerprint: artifact.fingerprint,
    source_run_fingerprint: artifact.run_fingerprint,
    source_comparison_ids: [...artifact.contrast_ids],
    browse_settings: { ...structuredClone(settings), table: structuredClone(table) },
  };
}

export function makeTargetEntry(artifact: Artifact, rows: Record<string, unknown>[], kind: TargetKind,
  filters: Record<string, unknown>): TargetBasketEntry {
  if (!rows.length || rows.some(row => !row.row_id || !row[kind === "gene" ? "gene_id" : "event_id"])) {
    throw new Error(`Selected rows do not contain exact ${kind} identities. Choose an annotated member/result table.`);
  }
  const unique = [...new Map(rows.map(row => [String(row.row_id), row])).values()];
  const rowIds = unique.map(row => String(row.row_id));
  return {
    basket_id: crypto.randomUUID(), source_artifact_id: artifact.artifact_id,
    source_name: artifact.name, row_ids: rowIds, target_kind: kind, filters: structuredClone(filters),
    labels: unique.map(row => kind === "gene" ? String(row.gene_symbol || row.gene_id) :
      `${row.gene_symbol || row.gene_id || ""} ${row.event_type || ""} ${row.event_id}`.trim()),
  };
}

export function mergeTargetEntries(entries: TargetBasketEntry[], entry: TargetBasketEntry) {
  const match = entries.find(item => item.source_artifact_id === entry.source_artifact_id &&
    item.target_kind === entry.target_kind && JSON.stringify(item.filters) === JSON.stringify(entry.filters));
  if (!match) return [...entries, entry];
  const labels = new Map(match.row_ids.map((id, index) => [id, match.labels[index]]));
  entry.row_ids.forEach((id, index) => labels.set(id, entry.labels[index]));
  return entries.map(item => item === match ? { ...item, row_ids: [...labels.keys()], labels: [...labels.values()] } : item);
}

function load(key: string): TargetBasketEntry[] {
  try {
    const saved = JSON.parse(localStorage.getItem(key) || "[]");
    return Array.isArray(saved) ? saved.filter(item => item && ["gene", "event"].includes(item.target_kind) &&
      typeof item.basket_id === "string" && typeof item.source_artifact_id === "string" &&
      Array.isArray(item.row_ids) && item.row_ids.every((id: unknown) => typeof id === "string") &&
      Array.isArray(item.labels) && item.filters && typeof item.filters === "object") : [];
  } catch { return []; }
}

export function useTargetBasket(projectId: string, fingerprint: string) {
  const key = `katana:target-basket:v1:${projectId}:${fingerprint}`;
  const [saved, setSaved] = useState<{key: string; entries: TargetBasketEntry[]}>(() => ({key, entries: load(key)}));
  const [storageError, setStorageError] = useState("");
  const entries = saved.key === key ? saved.entries : load(key);
  useEffect(() => {
    if (saved.key !== key || !projectId) return;
    try { localStorage.setItem(key, JSON.stringify(saved.entries)); setStorageError(""); }
    catch { setStorageError("Target basket persistence is unavailable. Export the request before closing this tab."); }
  }, [key, projectId, saved]);
  const change = (update: (current: TargetBasketEntry[]) => TargetBasketEntry[]) => {
    setSaved(current => ({key, entries: update(current.key === key ? current.entries : load(key))}));
  };
  return {entries, storageError, add: (entry: TargetBasketEntry) => change(current => mergeTargetEntries(current, entry)),
    remove: (id: string) => change(current => current.filter(entry => entry.basket_id !== id)), clear: () => change(() => [])};
}
