from __future__ import annotations

from pathlib import Path
import sqlite3
from tempfile import TemporaryDirectory

import pandas as pd


class MoleculeStore:
    """Bounded-memory junction deduplication on the selected output volume."""

    def __init__(self, output_directory: Path):
        output_directory.mkdir(parents=True, exist_ok=True)
        self.temporary = TemporaryDirectory(prefix=".junction-dedup-", dir=output_directory)
        self.connection = sqlite3.connect(str(Path(self.temporary.name) / "molecules.sqlite"))
        self.connection.execute("PRAGMA temp_store=MEMORY")
        self.connection.execute("PRAGMA cache_size=-16384")
        self.connection.execute("CREATE TABLE molecules (sample_id TEXT, cell_id TEXT, chrom TEXT, start INTEGER, end INTEGER, strand TEXT, molecule TEXT, PRIMARY KEY(sample_id, cell_id, chrom, start, end, strand, molecule)) WITHOUT ROWID")
        self.pending: list[tuple] = []

    def __enter__(self):
        return self

    def add(self, molecule: tuple) -> None:
        self.pending.append(molecule)
        if len(self.pending) >= 20000:
            self.flush()

    def flush(self) -> None:
        self.connection.executemany("INSERT OR IGNORE INTO molecules VALUES (?, ?, ?, ?, ?, ?, ?)", self.pending)
        self.connection.commit()
        self.pending.clear()

    def counts(self) -> pd.DataFrame:
        self.flush()
        return pd.read_sql_query("SELECT sample_id, cell_id, chrom, start, end, strand, COUNT(*) AS junction_count FROM molecules GROUP BY sample_id, cell_id, chrom, start, end, strand", self.connection)

    def __exit__(self, *_):
        self.connection.close()
        self.temporary.cleanup()
