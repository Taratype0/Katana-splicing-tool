from __future__ import annotations

from functools import lru_cache
import gzip
import hashlib
import json
from pathlib import Path
import re
from typing import TYPE_CHECKING

import pandas as pd

from src.analysis.common import normalize_gene_id
from src.domain.models.gene_filter import GeneFilterConfig

if TYPE_CHECKING:
    from src.domain.models.project import ProjectConfig

ID_COLUMNS = ("gene_id", "GeneID", "geneID", "groupID", "normalized_gene_id", "gene_id_A", "gene_id_B")
SYMBOL_COLUMNS = ("gene_symbol", "geneSymbol", "gene_name", "symbol", "gene", "gene_symbol_A", "gene_symbol_B")
TYPE_COLUMNS = ("gene_type", "gene_biotype", "biotype")
CHROM_COLUMNS = ("chr", "chromosome", "chrom", "seqname", "seqnames")


def _text(value: object) -> str:
    return "" if value is None or pd.isna(value) else str(value).strip().strip('\"')


def _key(value: object) -> str:
    return (normalize_gene_id(_text(value)) or "").casefold()


def _first(row: dict, columns: tuple[str, ...]) -> str:
    return next((_text(row.get(column)) for column in columns if _text(row.get(column))), "")


def _file_stamp(value: str | Path | None) -> tuple[str, int, int]:
    if not value:
        return ("", 0, 0)
    path = Path(value).expanduser().resolve()
    if not path.is_file():
        return (str(path), 0, 0)
    stat = path.stat()
    return (str(path), stat.st_size, stat.st_mtime_ns)


@lru_cache(maxsize=12)
def _sha256(stamp: tuple[str, int, int]) -> str:
    if not stamp[0] or not Path(stamp[0]).is_file():
        return "missing"
    digest = hashlib.sha256()
    with Path(stamp[0]).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _source_stamps(project: ProjectConfig) -> tuple:
    return (_file_stamp(project.tool_paths.get("annotation_tsv")),
            _file_stamp(project.tool_paths.get("gtf")),
            _file_stamp(project.counts_path if project.gene_filter.min_expression > 0 else None))


def gene_filter_signature(project: ProjectConfig) -> dict:
    config = project.gene_filter.to_dict()
    if not project.gene_filter.active:
        return {"schema": 1, "config": config}
    return {"schema": 1, "config": config, "sources": [
        {"path": stamp[0], "sha256": _sha256(stamp)}
        for stamp in _source_stamps(project) if stamp[0]
    ]}


def gene_filter_scope(project: ProjectConfig) -> str:
    encoded = json.dumps(gene_filter_signature(project), sort_keys=True).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()[:16]


@lru_cache(maxsize=4)
def _read_annotation(stamp: tuple[str, int, int]) -> pd.DataFrame:
    if not stamp[0] or not Path(stamp[0]).is_file():
        return pd.DataFrame()
    path = Path(stamp[0])
    if path.name.lower().endswith((".gtf", ".gtf.gz")):
        rows = []
        opener = gzip.open if path.suffix.lower() == ".gz" else open
        with opener(path, "rt", encoding="utf-8") as handle:
            for line in handle:
                if line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) != 9 or parts[2] != "gene":
                    continue
                attrs = dict(re.findall(r'(\S+)\s+"([^"]*)"', parts[8]))
                rows.append({"gene_id": attrs.get("gene_id", ""),
                             "gene_symbol": attrs.get("gene_name", ""),
                             "gene_type": attrs.get("gene_type") or attrs.get("gene_biotype", ""),
                             "chr": parts[0]})
        return pd.DataFrame(rows)
    return pd.read_csv(path, sep="," if path.suffix.lower() == ".csv" else "\t", low_memory=False)


class GeneFilter:
    def __init__(self, config: GeneFilterConfig, annotation: pd.DataFrame,
                 expression: pd.DataFrame | None = None) -> None:
        self.config = config
        self.by_id: dict[str, dict] = {}
        self.by_symbol: dict[str, list[dict]] = {}
        for row in annotation.to_dict("records"):
            gene_id = _first(row, ID_COLUMNS)
            if not gene_id:
                continue
            entry = {"gene_id": gene_id, "gene_symbol": _first(row, SYMBOL_COLUMNS),
                     "gene_type": _first(row, TYPE_COLUMNS), "chr": _first(row, CHROM_COLUMNS)}
            key = _key(gene_id)
            current = self.by_id.setdefault(key, entry)
            for name, value in entry.items():
                if not current.get(name) and value:
                    current[name] = value
        for entry in self.by_id.values():
            if entry["gene_symbol"]:
                self.by_symbol.setdefault(entry["gene_symbol"].casefold(), []).append(entry)
        self.custom = {_key(gene) for gene in config.excluded_genes}
        self.expression_keep: set[str] | None = None
        self.expression_columns: list[str] = []
        if config.min_expression > 0:
            if expression is None or expression.empty or "gene_id" not in expression.columns:
                raise ValueError("Expression filtering requires a gene_id expression matrix (Counts input).")
            metadata = set(ID_COLUMNS + SYMBOL_COLUMNS + TYPE_COLUMNS + CHROM_COLUMNS) | {
                "length", "gene_length", "start", "end", "strand", "width", "GC", "gc_content"}
            numeric = expression.drop(columns=[c for c in metadata if c in expression], errors="ignore")
            numeric = numeric.apply(pd.to_numeric, errors="coerce").dropna(axis=1, how="all")
            if numeric.shape[1] < config.min_samples:
                raise ValueError("Minimum samples exceeds the available expression sample columns.")
            self.expression_columns = list(numeric.columns)
            keep = numeric.ge(config.min_expression).sum(axis=1).ge(config.min_samples)
            self.expression_keep = {_key(gene) for gene in expression.loc[keep, "gene_id"]}

    def validate_annotation(self) -> None:
        entries = self.by_id.values()
        if (self.config.exclude_pseudogenes or self.config.exclude_rrna or self.config.protein_coding_only):
            if not any(entry["gene_type"] for entry in entries):
                raise ValueError("Gene-type filtering requires gene_type/gene_biotype annotation in an annotation TSV or GTF.")
        if self.config.exclude_mitochondrial and not any(entry["chr"] for entry in entries):
            raise ValueError("Mitochondrial filtering requires chromosome annotation in an annotation TSV or GTF.")

    def _reasons(self, row: dict) -> list[str]:
        ids = [_text(row.get(column)) for column in ID_COLUMNS if _text(row.get(column))]
        symbols = [_text(row.get(column)) for column in SYMBOL_COLUMNS if _text(row.get(column))]
        tokens = [part for value in ids for part in re.split(r"_and_|\+", value)]
        keys = {_key(value) for value in tokens + symbols}
        entries = [self.by_id[_key(value)] for value in tokens + symbols if _key(value) in self.by_id]
        # Prefer stable IDs. Symbol aliases are used only when IDs did not map.
        if not entries:
            entries = [entry for symbol in symbols for entry in self.by_symbol.get(symbol.casefold(), [])]
        keys.update(_key(entry["gene_id"]) for entry in entries)
        keys.update(_key(entry["gene_symbol"]) for entry in entries if entry["gene_symbol"])
        types = {entry["gene_type"].casefold() for entry in entries if entry["gene_type"]}
        chroms = {entry["chr"].casefold() for entry in entries if entry["chr"]}
        native_type, native_chr = _first(row, TYPE_COLUMNS), _first(row, CHROM_COLUMNS)
        if native_type:
            types.add(native_type.casefold())
        if native_chr:
            chroms.add(native_chr.casefold())
        reasons = []
        if keys & self.custom:
            reasons.append("custom_list")
        if self.config.exclude_mitochondrial and chroms & {"m", "mt", "chrm", "chrmt", "mitochondrion_genome"}:
            reasons.append("mitochondrial_encoded")
        if self.config.exclude_pseudogenes and any("pseudogene" in kind for kind in types):
            reasons.append("pseudogene")
        if self.config.exclude_rrna and any(kind in {"rrna", "mt_rrna"} or kind.startswith("rrna_") for kind in types):
            reasons.append("rrna")
        if self.config.protein_coding_only and types != {"protein_coding"}:
            reasons.append("not_protein_coding" if types else "unannotated_not_protein_coding")
        if self.expression_keep is not None and not keys & self.expression_keep:
            reasons.append("expression_threshold_or_missing")
        return reasons

    def decisions(self, frame: pd.DataFrame) -> pd.DataFrame:
        if frame.empty or not self.config.active or not any(c in frame for c in ID_COLUMNS + SYMBOL_COLUMNS):
            return pd.DataFrame({"excluded": False, "reasons": ""}, index=frame.index)
        columns = [c for c in ID_COLUMNS + SYMBOL_COLUMNS + TYPE_COLUMNS + CHROM_COLUMNS if c in frame]
        reasons = [self._reasons(row) for row in frame[columns].to_dict("records")]
        return pd.DataFrame({"excluded": [bool(value) for value in reasons],
                             "reasons": [";".join(value) for value in reasons]}, index=frame.index)

    def apply(self, frame: pd.DataFrame) -> pd.DataFrame:
        if not self.config.active or frame.empty:
            return frame.copy()
        return frame.loc[~self.decisions(frame)["excluded"]].copy()

    def annotation_audit(self) -> pd.DataFrame:
        frame = pd.DataFrame(self.by_id.values(), columns=["gene_id", "gene_symbol", "gene_type", "chr"])
        decisions = self.decisions(frame)
        return frame.join(decisions)


@lru_cache(maxsize=4)
def _engine(config_json: str, stamps: tuple) -> GeneFilter:
    config = GeneFilterConfig.from_mapping(json.loads(config_json))
    annotations = [_read_annotation(stamp) for stamp in stamps[:2]] if config.active else []
    annotations = [frame for frame in annotations if not frame.empty]
    annotation = pd.concat(annotations, ignore_index=True) if annotations else pd.DataFrame()
    expression = None
    if config.min_expression > 0 and stamps[2][0] and Path(stamps[2][0]).is_file():
        expression = pd.read_csv(stamps[2][0], sep="\t", low_memory=False)
    engine = GeneFilter(config, annotation, expression)
    if config.active:
        engine.validate_annotation()
    return engine


def project_gene_filter(project: ProjectConfig) -> GeneFilter:
    return _engine(json.dumps(project.gene_filter.to_dict(), sort_keys=True), _source_stamps(project))


def filter_project_frame(project: ProjectConfig, frame: pd.DataFrame) -> pd.DataFrame:
    if not project.gene_filter.active:
        return frame
    return project_gene_filter(project).apply(frame)
