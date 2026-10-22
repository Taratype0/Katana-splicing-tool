from __future__ import annotations

from dataclasses import dataclass
from contextlib import nullcontext
from pathlib import Path
import re
from zipfile import ZipFile

import numpy as np
import pandas as pd


IUPAC = {
    "A": "A", "C": "C", "G": "G", "T": "T", "U": "T",
    "R": "[AG]", "Y": "[CT]", "S": "[GC]", "W": "[AT]",
    "K": "[GT]", "M": "[AC]", "B": "[CGT]", "D": "[AGT]",
    "H": "[ACT]", "V": "[ACG]", "N": "[ACGT]",
}


@dataclass(frozen=True)
class Motif:
    rbp: str
    motif_id: str
    source: str
    consensus: str | None = None
    pwm: np.ndarray | None = None
    evidence_class: str = "unspecified"
    target_species: str = ""


class IndexedFasta:
    """Minimal faidx reader; avoids loading a multi-gigabyte genome into memory."""

    def __init__(self, fasta: Path) -> None:
        self.fasta = fasta
        fai = Path(str(fasta) + ".fai")
        if not fai.exists():
            raise FileNotFoundError(f"Missing FASTA index: {fai}")
        self.index: dict[str, tuple[int, int, int, int]] = {}
        for line in fai.read_text(encoding="utf-8").splitlines():
            fields = line.split("\t")
            if len(fields) >= 5:
                self.index[fields[0]] = tuple(map(int, fields[1:5]))

    def fetch(self, chrom: str, start: int, end: int, *, stream=None) -> str:
        if chrom not in self.index:
            return ""
        length, offset, line_bases, line_width = self.index[chrom]
        start = max(0, min(int(start), length))
        end = max(start, min(int(end), length))
        if end <= start:
            return ""
        byte_start = offset + (start // line_bases) * line_width + start % line_bases
        byte_end = offset + ((end - 1) // line_bases) * line_width + (end - 1) % line_bases + 1
        with (self.fasta.open("rb") if stream is None else nullcontext(stream)) as handle:
            handle.seek(byte_start)
            raw = handle.read(byte_end - byte_start)
        return raw.replace(b"\n", b"").replace(b"\r", b"").decode("ascii").upper()[: end - start]


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def _parse_attract_pwm(text: str) -> dict[str, np.ndarray]:
    matrices: dict[str, np.ndarray] = {}
    motif_id = ""
    expected = 0
    rows: list[list[float]] = []
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if motif_id and len(rows) == expected:
                matrices[motif_id] = np.asarray(rows, dtype=float)
            fields = line[1:].split()
            motif_id = fields[0]
            expected = int(fields[1]) if len(fields) > 1 else 0
            rows = []
            continue
        values = [float(value) for value in line.split()]
        if len(values) == 4:
            rows.append(values)
    if motif_id and len(rows) == expected:
        matrices[motif_id] = np.asarray(rows, dtype=float)
    return matrices


def _attract_motifs(metadata: pd.DataFrame, matrices: dict[str, np.ndarray], organism: str) -> list[Motif]:
    if organism and "Organism" in metadata:
        metadata = metadata.loc[metadata["Organism"].astype(str).str.casefold().eq(organism.casefold())]
    if "Mutated" in metadata:
        metadata = metadata.loc[metadata.Mutated.astype(str).str.casefold().eq("no")]
    if "Gene_id" in metadata and organism in {"Mus_musculus", "Homo_sapiens"}:
        prefix = "ENSMUSG" if organism == "Mus_musculus" else "ENSG"
        metadata = metadata.loc[metadata.Gene_id.astype(str).str.startswith(prefix)]
    motifs: list[Motif] = []
    for _, row in metadata.drop_duplicates(["Gene_name", "Matrix_id"]).iterrows():
        motif_id = str(row["Matrix_id"])
        pwm = matrices.get(motif_id)
        if pwm is None:
            continue
        motifs.append(Motif(
            rbp=str(row["Gene_name"]), motif_id=motif_id,
            source=f"ATtRACT:{row.get('Database', '')}", pwm=pwm,
            evidence_class="direct_experiment", target_species=organism,
        ))
    return motifs


def load_motif_database(path: Path, *, organism: str = "") -> list[Motif]:
    if path.suffix.lower() == ".zip":
        with ZipFile(path) as archive:
            names = set(archive.namelist())
            if {"ATtRACT_db.txt", "pwm.txt"}.issubset(names):
                metadata = pd.read_csv(archive.open("ATtRACT_db.txt"), sep="\t", low_memory=False)
                matrices = _parse_attract_pwm(archive.read("pwm.txt").decode("utf-8"))
                return _attract_motifs(metadata, matrices, organism)
            information = next((n for n in sorted(names) if n.endswith("RBP_Information_all_motifs.txt")), None)
            if information:
                return _cisbp_motifs(archive, information, organism)
        raise ValueError("Unsupported motif ZIP. Expected native ATtRACT or CISBP-RNA metadata and matrices.")
    if path.name == "ATtRACT_db.txt" and (path.parent / "pwm.txt").exists():
        metadata = pd.read_csv(path, sep="\t", low_memory=False)
        return _attract_motifs(metadata, _parse_attract_pwm((path.parent / "pwm.txt").read_text(encoding="utf-8")), organism)
    frame = pd.read_csv(path, sep="," if path.suffix.lower() == ".csv" else "\t", low_memory=False)
    if organism and "Organism" in frame:
        frame = frame.loc[frame["Organism"].astype(str).str.casefold().eq(organism.casefold())]
    rbp_col = next((c for c in ("RBP", "rbp", "gene_symbol", "protein", "Gene_name", "RBP_Name") if c in frame), None)
    id_col = next((c for c in ("motif_id", "Motif_ID", "Matrix_id", "PWM_id", "id") if c in frame), None)
    source_col = next((c for c in ("source", "database", "Database") if c in frame), None)
    if rbp_col is None:
        raise ValueError("Motif database requires an RBP/rbp/gene_symbol/protein column.")
    motifs: list[Motif] = []
    consensus_col = next((c for c in ("motif", "consensus", "sequence", "Motif") if c in frame), None)
    if consensus_col is not None:
        for index, row in frame.dropna(subset=[consensus_col]).iterrows():
            consensus = re.sub(r"[^ACGTURYSWKMBDHVN]", "", str(row[consensus_col]).upper())
            if len(consensus) >= 3:
                motifs.append(Motif(
                    rbp=str(row[rbp_col]),
                    motif_id=str(row[id_col]) if id_col else f"motif_{index + 1}",
                    source=str(row[source_col]) if source_col else path.stem,
                    consensus=consensus,
                ))
        return motifs
    required = {"position", "A", "C", "G", "T"}
    if not required.issubset(frame.columns) or id_col is None:
        raise ValueError(
            "Motif database must contain either a motif/consensus sequence or long-format PWM columns "
            "motif_id, position, A, C, G, T."
        )
    for motif_id, part in frame.groupby(id_col, sort=False, observed=True):
        part = part.sort_values("position")
        pwm = part[["A", "C", "G", "T"]].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
        row_sums = pwm.sum(axis=1, keepdims=True)
        if len(pwm) < 3 or not np.isfinite(pwm).all() or (row_sums <= 0).any():
            continue
        pwm = (pwm + 1e-6) / (row_sums + 4e-6)
        motifs.append(Motif(
            rbp=str(part.iloc[0][rbp_col]), motif_id=str(motif_id),
            source=str(part.iloc[0][source_col]) if source_col else path.stem, pwm=pwm,
            evidence_class=str(part.iloc[0].get("evidence_class", "unspecified")),
            target_species=organism,
        ))
    return motifs


def _cisbp_motifs(archive: ZipFile, information: str, organism: str) -> list[Motif]:
    metadata = pd.read_csv(archive.open(information), sep="\t", dtype=str, keep_default_na=False)
    required = {"RBP_ID", "RBP_Name", "RBP_Species", "RBP_Status", "Motif_ID", "Motif_Type"}
    if not required <= set(metadata):
        raise ValueError("Incomplete CISBP-RNA motif metadata")
    if not organism:
        raise ValueError("CISBP-RNA requires an explicit target species")
    metadata = metadata.loc[metadata.RBP_Species.str.replace(" ", "_").eq(organism.replace(" ", "_"))]
    matrices = {Path(n).stem: n for n in archive.namelist() if "/" in n and n.endswith(".txt") and Path(n).stem.startswith("M")}
    result = []
    for row in metadata.drop_duplicates(["RBP_ID", "Motif_ID"]).to_dict("records"):
        if row["Motif_ID"] in {"", "."}:
            continue
        if row["Motif_ID"] not in matrices:
            continue
        member = matrices[row["Motif_ID"]]
        # Native CISBP archives include unavailable-matrix placeholders.
        # Audit these separately instead of treating them as negative hits.
        if archive.getinfo(member).file_size == 0:
            continue
        frame = pd.read_csv(archive.open(member), sep="\t")
        if not {"Pos", "A", "C", "G", "U"} <= set(frame):
            raise ValueError("CISBP-RNA matrix needs Pos/A/C/G/U columns")
        frame = frame.sort_values("Pos")
        pwm = frame[["A", "C", "G", "U"]].to_numpy(dtype=float)
        if len(pwm) < 3 or not np.isfinite(pwm).all() or (pwm < 0).any() or (pwm.sum(axis=1) <= 0).any():
            raise ValueError(f"Invalid CISBP-RNA PWM: {row['Motif_ID']}")
        if row["Motif_Type"].casefold() == "jple":
            evidence = "model_predicted"
        elif row["RBP_Status"] == "D":
            evidence = "direct_experiment"
        elif row["RBP_Status"] == "I":
            evidence = "homology_inferred"
        else:
            raise ValueError(f"Unknown CISBP-RNA assignment: {row['RBP_Status']}")
        result.append(Motif(row["RBP_Name"], f"{row['RBP_ID']}::{row['Motif_ID']}",
            f"CISBP-RNA:2.00:{row['Motif_Type']}:PMID={row.get('PMID', '')}",
            pwm=pwm / pwm.sum(axis=1, keepdims=True), evidence_class=evidence, target_species=organism))
    return result


def motif_database_audit(path: Path, *, organism: str) -> pd.DataFrame:
    """Expose native records that cannot participate in scanning."""
    if path.suffix.lower() != ".zip":
        return pd.DataFrame(columns=["RBP", "motif_id", "status"])
    with ZipFile(path) as archive:
        info = next((n for n in archive.namelist() if n.endswith("RBP_Information_all_motifs.txt")), None)
        if info is None:
            return pd.DataFrame(columns=["RBP", "motif_id", "status"])
        records = pd.read_csv(archive.open(info), sep="\t", dtype=str, keep_default_na=False)
        records = records.loc[records.RBP_Species.str.replace(" ", "_").eq(organism.replace(" ", "_"))]
        matrices = {Path(n).stem: n for n in archive.namelist() if "/" in n and Path(n).stem.startswith("M")}
        rows = []
        for row in records.drop_duplicates(["RBP_ID", "Motif_ID"]).to_dict("records"):
            motif_id = row["Motif_ID"]
            member = matrices.get(motif_id)
            status = "no_assigned_motif" if motif_id in {"", "."} else "missing_matrix" if member is None else "empty_matrix" if archive.getinfo(member).file_size == 0 else "available"
            rows.append({"RBP": row["RBP_Name"], "RBP_ID": row["RBP_ID"], "motif_id": motif_id,
                         "native_assignment": row["RBP_Status"], "native_method": row["Motif_Type"], "status": status})
        return pd.DataFrame(rows)


def _event_windows(row: pd.Series, window_bp: int) -> list[tuple[int, int]]:
    coordinates = []
    for value in str(row.get("event_coordinates", "")).split("|"):
        try:
            coordinates.append(int(float(value)))
        except ValueError:
            continue
    windows = sorted({(max(0, coordinate - window_bp), coordinate + window_bp) for coordinate in coordinates})
    merged: list[tuple[int, int]] = []
    for start, end in windows:
        if not merged or start > merged[-1][1]:
            merged.append((start, end))
        else:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
    return merged


def _consensus_count(sequence: str, consensus: str) -> tuple[int, float]:
    pattern = "".join(IUPAC.get(base, re.escape(base)) for base in consensus)
    count = len(list(re.finditer(f"(?=({pattern}))", sequence)))
    return count, 1.0 if count else 0.0


def _encode_sequence(sequence: str) -> np.ndarray:
    lookup = np.full(256, -1, dtype=np.int8)
    for base, index in (("A", 0), ("C", 1), ("G", 2), ("T", 3)):
        lookup[ord(base)] = index
    return lookup[np.frombuffer(sequence.encode("ascii"), dtype=np.uint8)]


def _encoded_windows(encoded: np.ndarray, length: int) -> np.ndarray:
    if len(encoded) < length:
        return np.empty((0, length), dtype=np.int8)
    windows = np.lib.stride_tricks.sliding_window_view(encoded, length)
    return windows[np.all(windows >= 0, axis=1)]


def _pwm_count(windows: np.ndarray, pwm: np.ndarray, threshold: float) -> tuple[int, float]:
    length = len(pwm)
    if windows.size == 0:
        return 0, np.nan
    log_odds = np.log2(np.clip(pwm, 1e-9, 1.0) / 0.25)
    minimum = float(log_odds.min(axis=1).sum())
    maximum = float(log_odds.max(axis=1).sum())
    raw = log_odds[np.arange(length)[None, :], windows].sum(axis=1)
    relative = (raw - minimum) / (maximum - minimum) if maximum > minimum else np.zeros_like(raw)
    return int(np.sum(relative >= threshold)), float(np.max(relative))


def scan_event_motifs(
    events: pd.DataFrame,
    genome_fasta: Path,
    motif_database: Path,
    *,
    window_bp: int,
    score_threshold: float,
    organism: str = "",
) -> pd.DataFrame:
    fasta = IndexedFasta(genome_fasta)
    motifs = load_motif_database(motif_database, organism=organism)
    rows: list[dict[str, object]] = []
    for _, event in events.drop_duplicates("event_key").iterrows():
        chrom = str(event.get("chr", ""))
        strand = str(event.get("strand", "+"))
        windows = _event_windows(event, window_bp)
        if not windows:
            continue
        sequences = []
        for start, end in windows:
            sequence = fasta.fetch(chrom, start, end)
            sequences.append(reverse_complement(sequence) if strand == "-" else sequence)
        sequence = "N".join(sequences)
        encoded = _encode_sequence(sequence)
        encoded_window_cache: dict[int, np.ndarray] = {}
        for motif in motifs:
            if motif.consensus is not None:
                count, maximum = _consensus_count(sequence, motif.consensus)
                method = "IUPAC consensus scan"
            else:
                motif_length = len(motif.pwm)
                if motif_length not in encoded_window_cache:
                    encoded_window_cache[motif_length] = _encoded_windows(encoded, motif_length)
                count, maximum = _pwm_count(encoded_window_cache[motif_length], motif.pwm, score_threshold)
                method = "PWM relative log-likelihood scan"
            if count:
                rows.append({
                    "event_key": event["event_key"], "RBP": motif.rbp,
                    "motif_id": motif.motif_id, "motif_source": motif.source,
                    "site_count": count, "max_relative_score": maximum,
                    "window_bp": window_bp, "score_threshold": score_threshold,
                    "strand_oriented": True, "scan_method": method,
                })
    columns = [
        "event_key", "RBP", "motif_id", "motif_source", "site_count",
        "max_relative_score", "window_bp", "score_threshold", "strand_oriented", "scan_method",
    ]
    return pd.DataFrame(rows, columns=columns)
