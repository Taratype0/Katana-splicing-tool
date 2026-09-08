from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Any, Iterable

import duckdb
import pyarrow as pa
import pyarrow.parquet as pq


def _quote_identifier(value: str) -> str:
    return '"' + value.replace('"', '""') + '"'


class ResultStore:
    def __init__(self, output_root: str | Path) -> None:
        self.output_root = Path(output_root).resolve()

    def artifact_directory(self, analysis_id: str, run_fingerprint: str) -> Path:
        safe_analysis = "".join(
            character if character.isalnum() or character in "-_." else "_"
            for character in analysis_id
        )
        path = self.output_root / "artifacts" / safe_analysis / run_fingerprint[:16]
        path.mkdir(parents=True, exist_ok=True)
        return path

    def write_table(self, rows: Iterable[dict[str, Any]], path: str | Path) -> int:
        target = Path(path).resolve()
        target.parent.mkdir(parents=True, exist_ok=True)
        records = list(rows)
        table = pa.Table.from_pylist(records)
        temporary = target.with_suffix(target.suffix + ".tmp")
        pq.write_table(table, temporary, compression="zstd")
        os.replace(temporary, target)
        return table.num_rows

    def write_plot_spec(self, specification: dict[str, Any], path: str | Path) -> Path:
        target = Path(path).resolve()
        target.parent.mkdir(parents=True, exist_ok=True)
        temporary = target.with_suffix(target.suffix + ".tmp")
        temporary.write_text(
            json.dumps(specification, ensure_ascii=False, indent=2), encoding="utf-8"
        )
        os.replace(temporary, target)
        return target

    def table_page(
        self,
        path: str | Path,
        *,
        offset: int = 0,
        limit: int = 100,
        search: str = "",
        sort_by: str | None = None,
        descending: bool = False,
    ) -> dict[str, Any]:
        source = Path(path).resolve()
        if not source.is_file():
            raise FileNotFoundError(source)
        if source.suffix.lower() != ".parquet":
            raise ValueError(f"Unsupported table format: {source.suffix}")

        with duckdb.connect() as connection:
            schema_rows = connection.execute(
                "DESCRIBE SELECT * FROM read_parquet(?)", [str(source)]
            ).fetchall()
            columns = [str(row[0]) for row in schema_rows]
            where = ""
            parameters: list[Any] = [str(source)]
            if search.strip() and columns:
                expression = " || ' ' || ".join(
                    f"COALESCE(CAST({_quote_identifier(column)} AS VARCHAR), '')"
                    for column in columns
                )
                where = f" WHERE ({expression}) ILIKE ?"
                parameters.append(f"%{search.strip()}%")

            count_row = connection.execute(
                f"SELECT COUNT(*) FROM read_parquet(?) {where}", parameters
            ).fetchone()
            total = int(count_row[0]) if count_row is not None else 0
            order = ""
            if sort_by:
                if sort_by not in columns:
                    raise ValueError(f"Unknown sort column: {sort_by}")
                order = f" ORDER BY {_quote_identifier(sort_by)}"
                if descending:
                    order += " DESC"
                order += " NULLS LAST"
            page_parameters = [
                *parameters,
                max(1, min(int(limit), 500)),
                max(0, int(offset)),
            ]
            cursor = connection.execute(
                f"SELECT * FROM read_parquet(?) {where}{order} LIMIT ? OFFSET ?",
                page_parameters,
            )
            records = [
                dict(zip(columns, row, strict=True)) for row in cursor.fetchall()
            ]

        normalized = [
            {
                key: value.item() if hasattr(value, "item") else value
                for key, value in record.items()
            }
            for record in records
        ]
        return {
            "columns": columns,
            "rows": normalized,
            "offset": max(0, int(offset)),
            "limit": max(1, min(int(limit), 500)),
            "total": int(total),
        }

    @staticmethod
    def read_plot_spec(path: str | Path) -> dict[str, Any]:
        return json.loads(Path(path).read_text(encoding="utf-8"))
