import { useDeferredValue, useEffect, useRef, useState } from "react";
import { useQuery } from "@tanstack/react-query";
import {
  type ColumnDef,
  type SortingState,
  flexRender,
  getCoreRowModel,
  useReactTable,
} from "@tanstack/react-table";
import { useVirtualizer } from "@tanstack/react-virtual";
import { ArrowDown, ArrowUp, ChevronLeft, ChevronRight, Search } from "lucide-react";
import { api } from "../api";
import type { Artifact } from "../types";

type Props = {
  projectId: string;
  artifact?: Artifact;
  selectedIds: Set<string>;
  onToggle: (rowId: string) => void;
  onFocus: (row: Record<string, unknown>) => void;
};

const PAGE_SIZE = 50;

function cellText(value: unknown) {
  if (value === null || value === undefined || value === "") return "NA";
  if (typeof value === "number" && !Number.isInteger(value)) return value.toPrecision(4);
  return String(value);
}

export function ArtifactTable({ projectId, artifact, selectedIds, onToggle, onFocus }: Props) {
  const [search, setSearch] = useState("");
  const deferredSearch = useDeferredValue(search);
  const [offset, setOffset] = useState(0);
  const [sorting, setSorting] = useState<SortingState>([]);
  const scrollHost = useRef<HTMLDivElement>(null);

  useEffect(() => setOffset(0), [deferredSearch, artifact?.artifact_id]);

  const pageQuery = useQuery({
    queryKey: ["artifact-rows", projectId, artifact?.artifact_id, offset, deferredSearch, sorting],
    queryFn: () =>
      api.artifactRows(projectId, artifact!.artifact_id, {
        offset,
        limit: PAGE_SIZE,
        search: deferredSearch,
        sortBy: sorting[0]?.id,
        descending: sorting[0]?.desc,
      }),
    enabled: Boolean(artifact),
    placeholderData: (previous) => previous,
  });

  const columns: ColumnDef<Record<string, unknown>>[] = (pageQuery.data?.columns ?? []).map((column) => ({
    id: column,
    accessorKey: column,
    header: column.replaceAll("_", " "),
    cell: (context) => cellText(context.getValue()),
  }));
  const rows = pageQuery.data?.rows ?? [];
  const table = useReactTable({
    data: rows,
    columns,
    state: { sorting },
    manualSorting: true,
    onSortingChange: setSorting,
    getCoreRowModel: getCoreRowModel(),
    getRowId: (row, index) => String(row.row_id ?? `${offset + index}`),
  });
  const tableRows = table.getRowModel().rows;
  const virtualizer = useVirtualizer({
    count: tableRows.length,
    getScrollElement: () => scrollHost.current,
    estimateSize: () => 43,
    overscan: 8,
  });
  const gridTemplate = `46px repeat(${Math.max(1, columns.length)}, minmax(138px, 1fr))`;
  const gridWidth = `${46 + Math.max(1, columns.length) * 138}px`;
  const total = pageQuery.data?.total ?? 0;

  if (!artifact) {
    return <div className="table-empty">No cached table is available.</div>;
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
        <span>{total.toLocaleString()} rows</span>
      </div>
      {pageQuery.error ? (
        <div className="table-empty error">{pageQuery.error.message}</div>
      ) : (
        <div className="data-grid-shell">
          <div className="data-grid-header" style={{ gridTemplateColumns: gridTemplate, minWidth: gridWidth }}>
            <span className="selection-column">Select</span>
            {table.getHeaderGroups()[0]?.headers.map((header) => {
              const direction = header.column.getIsSorted();
              return (
                <button
                  key={header.id}
                  type="button"
                  className="column-header"
                  onClick={header.column.getToggleSortingHandler()}
                >
                  {flexRender(header.column.columnDef.header, header.getContext())}
                  {direction === "asc" && <ArrowUp size={13} />}
                  {direction === "desc" && <ArrowDown size={13} />}
                </button>
              );
            })}
          </div>
          <div ref={scrollHost} className="data-grid-scroll">
            <div
              className="virtual-grid"
              style={{ height: `${virtualizer.getTotalSize()}px`, minWidth: gridWidth }}
            >
              {virtualizer.getVirtualItems().map((virtualRow) => {
                const row = tableRows[virtualRow.index];
                return (
                  <div
                    key={row.id}
                    className="data-grid-row"
                    style={{
                      gridTemplateColumns: gridTemplate,
                      transform: `translateY(${virtualRow.start}px)`,
                    }}
                    onClick={() => onFocus(row.original)}
                  >
                    <label className="row-checkbox" onClick={(event) => event.stopPropagation()}>
                      <input
                        type="checkbox"
                        checked={selectedIds.has(row.id)}
                        onChange={() => onToggle(row.id)}
                        aria-label={`Select ${row.id}`}
                      />
                      <span />
                    </label>
                    {row.getVisibleCells().map((cell) => (
                      <span key={cell.id} className="data-cell" title={cellText(cell.getValue())}>
                        {flexRender(cell.column.columnDef.cell, cell.getContext())}
                      </span>
                    ))}
                  </div>
                );
              })}
            </div>
          </div>
        </div>
      )}
      <div className="table-pagination">
        <span>
          {total ? `${offset + 1}-${Math.min(offset + PAGE_SIZE, total)} of ${total}` : "0 rows"}
        </span>
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
