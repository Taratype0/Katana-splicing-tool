from __future__ import annotations

from pathlib import Path

import pandas as pd

from katana_cosplice.audit import resolve_rmats_root
from katana_cosplice.config import CoSpliceConfig


def _metadata_candidate(project_root: Path) -> Path:
    candidates = list(project_root.glob("**/pipeline_info/samplesheet.valid.csv"))
    ranked: list[tuple[int, Path]] = []
    for path in candidates:
        try:
            columns = set(pd.read_csv(path, nrows=2).columns)
        except Exception:
            continue
        if {"sample", "condition"}.issubset(columns) or {"sample_id", "condition"}.issubset(columns):
            score = 10 if "Splicing_genomebam" in str(path) else 0
            score += len(pd.read_csv(path, usecols=[next(iter(columns.intersection({"sample", "sample_id"})))]))
            ranked.append((score, path))
    if not ranked:
        raise FileNotFoundError(f"No sample metadata with sample/condition columns found under {project_root}")
    return max(ranked, key=lambda item: item[0])[1]


def _find_reference(search_root: Path, names: tuple[str, ...]) -> Path | None:
    for name in names:
        direct = search_root / name
        if direct.exists():
            return direct
    for name in names:
        matches = list(search_root.glob(f"**/{name}"))
        if matches:
            return matches[0]
    return None


def discover_config(project_root: str | Path, output_dir: str | Path | None = None) -> CoSpliceConfig:
    root = Path(project_root).resolve()
    rmats_root = resolve_rmats_root(root)
    metadata = _metadata_candidate(root)
    quant_candidates = [path.parent.parent for path in root.glob("**/quant.sf")]
    quant_root = None
    if quant_candidates:
        counts: dict[Path, int] = {}
        for candidate in quant_candidates:
            counts[candidate] = counts.get(candidate, 0) + 1
        quant_root = max(counts, key=counts.get)

    reference_root = root
    for parent in root.parents:
        if parent.name.lower() == "rnaseq":
            reference_root = parent
            break
    output = Path(output_dir).resolve() if output_dir else root / "katana_output" / "09_cosplice_network"
    support_root = output.parent / "02_candidate_gene_screening"
    return CoSpliceConfig(
        project_name=root.name,
        rmats_root=rmats_root,
        sample_metadata=metadata,
        output_dir=output,
        quant_root=quant_root,
        gtf=_find_reference(reference_root, ("gencode.vM34.annotation.gtf", "gencode.annotation.gtf")),
        genome_fasta=_find_reference(reference_root, ("GRCm39.primary_assembly.genome.fa",)),
        transcript_fasta=_find_reference(reference_root, ("gencode.vM34.transcripts.fa",)),
        support_root=support_root if support_root.exists() else None,
    )
