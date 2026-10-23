import { useDeferredValue, useEffect, useState } from "react";
import { useQuery } from "@tanstack/react-query";
import { ArrowDown, ArrowUp, ChevronLeft, ChevronRight, Search } from "lucide-react";
import { api } from "../api";
import type { Artifact } from "../types";

type Props = {
  projectId: string;
  artifact?: Artifact;
  selectedIds: Set<string>;
  onToggle: (rowId: string, row?: Record<string, unknown>) => void;
  onFocus: (row: Record<string, unknown>) => void;
  onBrowseSettingsChange?: (settings: Record<string, unknown>) => void;
};

const PAGE_SIZE = 50;

function cellText(value: unknown) {
  if (value === null || value === undefined || value === "") return "NA";
  if (typeof value === "boolean") return value ? "Yes" : "No";
  if (typeof value === "number" && !Number.isInteger(value)) {
    return Math.abs(value) > 0 && Math.abs(value) < 0.001 ? value.toExponential(2) : value.toPrecision(4);
  }
  return String(value);
}

export function ArtifactTable({ projectId, artifact, selectedIds, onToggle, onFocus, onBrowseSettingsChange }: Props) {
  const [search, setSearch] = useState("");
  const deferredSearch = useDeferredValue(search);
  const [offset, setOffset] = useState(0);
  const [sortBy, setSortBy] = useState<string>();
  const [descending, setDescending] = useState(false);

  useEffect(() => setOffset(0), [deferredSearch, artifact?.artifact_id]);
  useEffect(() => {onBrowseSettingsChange?.({search: deferredSearch, sort_by: sortBy || null,
    descending, page_offset: offset, page_size: PAGE_SIZE});}, [deferredSearch, sortBy, descending, offset]);
  const pageQuery = useQuery({
    queryKey: ["artifact-rows", projectId, artifact?.artifact_id, offset, deferredSearch, sortBy, descending],
    queryFn: () => api.artifactRows(projectId, artifact!.artifact_id, {
      offset,
      limit: PAGE_SIZE,
      search: deferredSearch,
      sortBy,
      descending,
    }),
    enabled: Boolean(artifact),
    placeholderData: (previous) => previous,
  });

  if (!artifact) return <div className="table-empty">No cached table is available.</div>;
  const columns = pageQuery.data?.columns ?? [];
  const rows = pageQuery.data?.rows ?? [];
  const total = pageQuery.data?.total ?? 0;

  function sort(column: string) {
    if (sortBy === column) setDescending((current) => !current);
    else {
      setSortBy(column);
      setDescending(false);
    }
    setOffset(0);
  }

  return (
    <div className="table-panel">
      <div className="table-toolbar">
        <label className="table-search">
          <Search size={16} />
          <input
            value={search}
            onChange={(event) => setSearch(event.target.value)}
            placeholder="Search event, gene, state, or value"
            aria-label="Search table"
          />
        </label>
        <span>{total.toLocaleString()} rows · 50 per page</span>
      </div>
      {pageQuery.error ? (
        <div className="table-empty error">{pageQuery.error.message}</div>
      ) : (
        <div className="data-table-scroll">
          <table className="data-table">
            <thead>
              <tr>
                <th className="select-cell" aria-label="Selection" />
                {columns.map((column) => (
                  <th key={column}>
                    <button type="button" onClick={() => sort(column)}>
                      {column.replaceAll("_", " ")}
                      {sortBy === column && (descending ? <ArrowDown size={13} /> : <ArrowUp size={13} />)}
                    </button>
                  </th>
                ))}
              </tr>
            </thead>
            <tbody>
              {rows.map((row, index) => {
                const rowId = String(row.row_id ?? `${offset + index}`);
                return (
                  <tr key={rowId} onClick={() => onFocus(row)}>
                    <td className="select-cell" onClick={(event) => event.stopPropagation()}>
                      <input
                        type="checkbox"
                        checked={selectedIds.has(rowId)}
                        onChange={() => onToggle(rowId, row)}
                        aria-label={`Select ${rowId}`}
                      />
                    </td>
                    {columns.map((column) => <td key={column} title={cellText(row[column])}>{cellText(row[column])}</td>)}
                  </tr>
                );
              })}
            </tbody>
          </table>
          {!rows.length && !pageQuery.isLoading && <div className="table-empty">No rows match the current search.</div>}
        </div>
      )}
      <div className="table-pagination">
        <span>{total ? `${offset + 1}-${Math.min(offset + PAGE_SIZE, total)} of ${total}` : "0 rows"}</span>
        <div>
          <button type="button" disabled={offset === 0} onClick={() => setOffset(Math.max(0, offset - PAGE_SIZE))}>
            <ChevronLeft size={16} />Previous
          </button>
          <button type="button" disabled={offset + PAGE_SIZE >= total} onClick={() => setOffset(offset + PAGE_SIZE)}>
            Next<ChevronRight size={16} />
          </button>
        </div>
      </div>
    </div>
  );
}
