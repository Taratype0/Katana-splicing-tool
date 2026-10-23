#!/usr/bin/env python3
"""Export genuine upstream results as an independently verifiable Katana bundle.

Only Python's standard library is required. No Katana analysis code is imported.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import re
import shutil
from datetime import datetime, timezone
from pathlib import Path, PurePosixPath

SCHEMA = "katana-pipeline-bundle/v1"
EVENTS = ("SE", "RI", "A3SS", "A5SS", "MXE")
MODES = ("JC", "JCEC")
RMATS_COLUMNS = {
    "ID", "GeneID", "geneSymbol", "chr", "strand", "PValue", "FDR",
    "IncLevelDifference", "IncLevel1", "IncLevel2", "IJC_SAMPLE_1",
    "IJC_SAMPLE_2", "SJC_SAMPLE_1", "SJC_SAMPLE_2", "IncFormLen", "SkipFormLen",
}
COORDINATES = {
    "SE": {"exonStart_0base", "exonEnd", "upstreamES", "upstreamEE", "downstreamES", "downstreamEE"},
    "RI": {"riExonStart_0base", "riExonEnd", "upstreamES", "upstreamEE", "downstreamES", "downstreamEE"},
    "A3SS": {"longExonStart_0base", "longExonEnd", "shortES", "shortEE", "flankingES", "flankingEE"},
    "A5SS": {"longExonStart_0base", "longExonEnd", "shortES", "shortEE", "flankingES", "flankingEE"},
    "MXE": {"1stExonStart_0base", "1stExonEnd", "2ndExonStart_0base", "2ndExonEnd", "upstreamES", "upstreamEE", "downstreamES", "downstreamEE"},
}
TABLE_COLUMNS = {
    "deg": {"gene_id", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"},
    "normalized_counts": {"gene_id"},
    "raw_counts": {"gene_id"},
    "quant": {"Name", "Length", "EffectiveLength", "TPM", "NumReads"},
    "dexseq": {"groupID", "padj"},
    "dtu": {"groupID", "padj"},
    "stager": {"geneID", "txID", "gene", "transcript"},
    "annotation": {"gene_id", "gene_name"},
}
KINDS = set(TABLE_COLUMNS) | {
    "rmats", "suppa", "gtf", "genome_fasta", "transcript_fasta", "bam", "bai",
    "rmats_bamlist", "matrix", "provenance", "qc", "auxiliary",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def safe_relative(value: str) -> str:
    value = str(value).replace("\\", "/")
    path = PurePosixPath(value)
    if not value or path.is_absolute() or ".." in path.parts or ":" in value or path.parts[0] == ".":
        raise ValueError(f"Unsafe bundle relative path: {value!r}")
    return path.as_posix()


def inside(root: Path, relative: str) -> Path:
    destination = root / safe_relative(relative)
    if not destination.resolve().is_relative_to(root.resolve()):
        raise ValueError(f"Bundle path escapes its root: {relative}")
    return destination


def read_metadata(path: Path, aliases: dict[str, str], required: set[str]) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t" if path.suffix == ".tsv" else ",")
        rows = [{aliases.get(k, k): (v or "").strip() for k, v in row.items()} for row in reader]
        columns = {aliases.get(k, k) for k in (reader.fieldnames or [])}
    if required - columns or not rows:
        raise ValueError(f"Metadata missing {sorted(required - columns)} or contains no rows: {path}")
    if any(not row.get(key) for row in rows for key in required):
        raise ValueError(f"Metadata has blank required values: {path}")
    return rows


def identifier(value: str) -> str:
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", value):
        raise ValueError(f"Use a portable identifier of letters, digits, '.', '_' or '-': {value!r}")
    return value


def write_csv(path: Path, rows: list[dict], first: list[str]) -> None:
    columns = first + sorted(set().union(*(row.keys() for row in rows)) - set(first))
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


def suppa_info(path: Path) -> dict:
    """Read SUPPA 2.3's native unnamed event index without rewriting the file.

    Native merge_temp_output_files emits only statistical column names in the
    header, then pandas.to_csv(header=False, index=True) adds event IDs to rows.
    """
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, [])
        implicit_index = bool(header) and len(header) % 2 == 0
        value_columns = header if implicit_index else header[1:]
        if len(value_columns) < 2 or len(value_columns) % 2 or len(header) != len(set(header)):
            raise ValueError(f"Invalid SUPPA result header: {path}")
        for position, column in enumerate(value_columns):
            lower = column.lower()
            valid = lower.endswith("_dpsi") or lower == "dpsi" if position % 2 == 0 else lower.endswith("_p-val") or lower in {"pvalue", "pval", "p-val"}
            if not valid:
                raise ValueError(f"Invalid SUPPA statistical column {column!r}: {path}")
        count = 0
        for fields in reader:
            count += 1
            if len(fields) != len(value_columns) + 1 or not fields[0]:
                raise ValueError(f"Malformed SUPPA row {count}, including its event index: {path}")
            for position, value in enumerate(fields[1:]):
                if value.lower() in {"", "na", "nan"}:
                    continue
                number = float(value)
                lower_bound = -1 if position % 2 == 0 else 0
                if not math.isfinite(number) or not lower_bound <= number <= 1:
                    raise ValueError(f"Invalid SUPPA statistic {value!r} in row {count}: {path}")
    return {"columns": header, "rows": count, "status": "available" if count else "available_empty",
            "event_id_representation": "implicit_index" if implicit_index else "explicit_column", "data_fields": len(value_columns) + 1}


def table_info(path: Path, kind: str, samples: list[dict], event: str | None = None) -> dict:
    if kind == "suppa":
        return suppa_info(path)
    delimiter = "," if path.suffix.lower() == ".csv" else "\t"
    required = RMATS_COLUMNS | COORDINATES[event] if kind == "rmats" else TABLE_COLUMNS.get(kind, set())
    if kind in {"dexseq", "dtu"} and "DEXSeqResults." in path.name:
        required = required | {"featureID"}
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle, delimiter=delimiter)
        columns = next(reader, [])
        duplicates = sorted({column for column in columns if columns.count(column) > 1})
        # Native rMATS joins event annotations and count/statistic tables, both
        # containing ID. Validate the two original values before a DictReader
        # could silently discard one; all other duplicate headers are invalid.
        native_rmats_ids = kind == "rmats" and duplicates == ["ID"] and columns.count("ID") == 2
        if not columns or (duplicates and not native_rmats_ids) or required - set(columns):
            raise ValueError(f"Invalid {kind} header; missing {sorted(required - set(columns))}; duplicate columns {duplicates}: {path}")
        id_positions = [i for i, column in enumerate(columns) if column == "ID"] if native_rmats_ids else []
        if kind in {"normalized_counts", "raw_counts"}:
            missing = {sample["sample_id"] for sample in samples} - set(columns)
            if missing:
                raise ValueError(f"Count matrix missing sample columns {sorted(missing)}: {path}")
        count = 0
        for fields in reader:
            if not fields:  # Match DictReader's handling of blank physical lines.
                continue
            count += 1
            if len(fields) != len(columns):
                raise ValueError(f"Malformed {kind} row {count}: {path}")
            if id_positions and fields[id_positions[0]] != fields[id_positions[1]]:
                raise ValueError(f"Mismatched native rMATS ID values in row {count}: {path}")
    return {"columns": columns, "rows": count, "status": "available" if count else "available_empty"}


def vector(value: str, expected: int, *, probability: bool, path: Path) -> list[str]:
    values = value.split(",")
    if len(values) != expected:
        raise ValueError(f"rMATS vector has {len(values)} values, expected {expected} samples: {path}")
    for item in values:
        if item.lower() in {"na", "nan", ""}:
            continue
        number = float(item)
        if not math.isfinite(number) or number < 0 or (probability and number > 1):
            raise ValueError(f"Invalid rMATS {'PSI' if probability else 'read count'} value {item!r}: {path}")
    return values


def make_matrices(root: Path, artifacts: list[dict], comparisons: list[dict], order: dict[str, list[str]]) -> list[str]:
    """Expand original sample vectors without averaging or statistical testing."""
    generated = []
    by_comparison = {row["comparison_id"]: row for row in comparisons}
    groups: dict[tuple[str, str], list[dict]] = {}
    for artifact in artifacts:
        if artifact["kind"] == "rmats" and artifact.get("relative_path") and artifact["status"].startswith("available"):
            groups.setdefault((artifact["comparison_id"], artifact["mode"]), []).append(artifact)
    for (comparison_id, mode), sources in sorted(groups.items()):
        comparison = by_comparison[comparison_id]
        group1, group2 = comparison["rmats_group1"], comparison["rmats_group2"]
        ids1, ids2 = order[group1], order[group2]
        paths = {key: root / "matrices" / comparison_id / mode / f"{key}.tsv" for key in ("psi", "ijc", "sjc")}
        handles = {}
        try:
            writers = {}
            for key, path in paths.items():
                path.parent.mkdir(parents=True, exist_ok=True)
                handles[key] = path.open("w", encoding="utf-8", newline="")
                writers[key] = csv.writer(handles[key], delimiter="\t")
                writers[key].writerow(["event_key", *ids1, *ids2])
                generated.append(path.relative_to(root).as_posix())
            seen = set()
            for source in sorted(sources, key=lambda item: item["event_type"]):
                path = inside(root, source["relative_path"])
                with path.open(encoding="utf-8-sig", newline="") as handle:
                    for row in csv.DictReader(handle, delimiter="\t"):
                        event = source["event_type"]
                        keys = [event, row["GeneID"], row["chr"], row["strand"], *(row[key] for key in sorted(COORDINATES[event]))]
                        event_key = "|".join(keys)
                        if event_key in seen:
                            raise ValueError(f"Duplicate rMATS event coordinate key: {path}: {event_key}")
                        seen.add(event_key)
                        for key, prefix in (("psi", "IncLevel"), ("ijc", "IJC_SAMPLE_"), ("sjc", "SJC_SAMPLE_")):
                            left = vector(row[f"{prefix}1"], len(ids1), probability=key == "psi", path=path)
                            right = vector(row[f"{prefix}2"], len(ids2), probability=key == "psi", path=path)
                            writers[key].writerow([event_key, *left, *right])
        finally:
            for handle in handles.values():
                handle.close()
    return generated


def record(root: Path, relative: str, kind: str, **extra) -> dict:
    path = inside(root, relative)
    return {"kind": kind, "relative_path": relative, "bytes": path.stat().st_size, "sha256": sha256(path), "status": "available", **extra}


def export_bundle(source_map_path: Path, output_dir: Path, samples_path: Path, contrasts_path: Path, profile: str = "core", allow_partial: bool = False) -> dict:
    source_map_path = source_map_path.resolve()
    data = json.loads(source_map_path.read_text(encoding="utf-8-sig"))
    samples = read_metadata(samples_path, {"sample": "sample_id", "group": "condition"}, {"sample_id", "condition"})
    contrasts = read_metadata(contrasts_path, {"comparison_id": "contrast", "target": "treatment", "reference": "control"}, {"contrast", "treatment", "control"})
    sample_ids = [identifier(row["sample_id"]) for row in samples]
    if len(set(sample_ids)) != len(sample_ids):
        raise ValueError("Sample metadata contains duplicate sample IDs")
    condition_ids = {identifier(row["condition"]) for row in samples}
    if len({row["contrast"] for row in contrasts}) != len(contrasts):
        raise ValueError("Duplicate contrast IDs")
    overrides = {row["comparison_id"]: row for row in data.get("comparisons", [])}
    comparisons = []
    for contrast in contrasts:
        comparison_id = identifier(contrast["contrast"])
        treatment, control = contrast["treatment"], contrast["control"]
        if treatment == control or {treatment, control} - condition_ids:
            raise ValueError(f"Contrast refers to invalid conditions: {contrast}")
        override = overrides.get(comparison_id, {})
        comparison = {"comparison_id": comparison_id, "treatment": treatment, "control": control,
                      "rmats_group1": override.get("rmats_group1", treatment), "rmats_group2": override.get("rmats_group2", control),
                      "deg_target": override.get("deg_target", treatment), "deg_reference": override.get("deg_reference", control)}
        if {comparison["rmats_group1"], comparison["rmats_group2"]} != {treatment, control} or {comparison["deg_target"], comparison["deg_reference"]} != {treatment, control}:
            raise ValueError(f"Comparison source direction does not match the declared biological contrast: {comparison_id}")
        comparisons.append(comparison)
    order = data.get("rmats_sample_order", {})
    for condition in condition_ids:
        expected = {row["sample_id"] for row in samples if row["condition"] == condition}
        actual = order.get(condition, [])
        if len(actual) != len(set(actual)) or set(actual) != expected:
            raise ValueError(f"Explicit upstream rMATS sample order is missing or inconsistent for {condition}")
    sources = data.get("artifacts", [])
    if not sources:
        raise ValueError("Source map must contain explicit artifacts")
    output_dir = output_dir.resolve()
    if output_dir.exists():
        raise ValueError(f"Refusing to overwrite an existing bundle directory: {output_dir}")
    comparison_ids = {row["comparison_id"] for row in comparisons}
    reserved = {"manifest.json", "SHA256SUMS", ".katana_project.yaml", "pipeline_info/samplesheet.valid.csv", "pipeline_info/contrastsheet.valid.csv", "provenance/source_map.json"}
    destinations = set(reserved)
    checked = []
    for item in sources:
        item = dict(item)
        kind = item.get("kind")
        if kind not in KINDS:
            raise ValueError(f"Unknown artifact kind: {kind}")
        status = item.get("status", "available")
        if status not in {"available", "available_empty", "not_computed", "failed"}:
            raise ValueError(f"Unknown artifact status: {status}")
        if not status.startswith("available"):
            if not item.get("reason"):
                raise ValueError(f"{status} requires a reason: {kind}")
            if not allow_partial:
                raise ValueError(f"Incomplete source artifact {kind} ({status}): {item['reason']}; partial export requires explicit --allow-partial")
            checked.append(item)
            continue
        relative = safe_relative(item["relative_path"])
        if relative in destinations or relative.startswith("matrices/") or relative.startswith("rmats/bamlist/"):
            raise ValueError(f"Duplicate or reserved bundle destination: {relative}")
        destinations.add(relative)
        path = Path(item["source"])
        if not path.is_absolute():
            path = source_map_path.parent / path
        path = path.resolve()
        if not path.is_file() or path.stat().st_size == 0:
            raise ValueError(f"Available artifact is missing or zero bytes: {path}")
        item.update(source=str(path), relative_path=relative)
        if kind in {"deg", "dexseq", "dtu", "stager", "suppa"} and item.get("comparison_id") not in comparison_ids:
            raise ValueError(f"{kind} artifact requires a declared comparison_id")
        if kind in {"quant", "bam", "bai"} and item.get("sample_id") not in sample_ids:
            raise ValueError(f"{kind} artifact requires a declared sample_id")
        if item.get("sample_id") and item["sample_id"] not in sample_ids:
            raise ValueError(f"Artifact refers to unknown sample: {item['sample_id']}")
        if item.get("comparison_id") and item["comparison_id"] not in comparison_ids:
            raise ValueError(f"Artifact refers to unknown comparison: {item['comparison_id']}")
        if kind == "deg" and PurePosixPath(relative).name != f"{item['comparison_id']}.deseq2.results.tsv":
            raise ValueError(f"DEG destination filename must preserve its declared comparison ID: {relative}")
        if kind == "normalized_counts" and PurePosixPath(relative).name != "all.normalised_counts.tsv":
            raise ValueError(f"Normalized counts need the current scanner's exact filename: {relative}")
        if kind == "quant" and not relative.endswith(f"/{item['sample_id']}/quant.sf"):
            raise ValueError(f"quant.sf destination parent must equal sample_id: {relative}")
        if kind in {"dexseq", "dtu", "stager"}:
            suffix = ".csv" if kind == "dexseq" else ".tsv"
            prefixes = {"dexseq": ("perGeneQValue", "DEXSeqResults"), "dtu": ("perGeneQValue", "DEXSeqResults"), "stager": ("getAdjustedPValues",)}[kind]
            if PurePosixPath(relative).name not in {f"{prefix}.{item['comparison_id']}{suffix}" for prefix in prefixes}:
                raise ValueError(f"{kind} filename is not recognized by current Katana: {relative}")
        if kind == "suppa" and (item["comparison_id"] not in PurePosixPath(relative).parts[:-1] or not relative.endswith("_local_diffsplice.dpsi")):
            raise ValueError(f"SUPPA requires a comparison-named directory and native .dpsi filename: {relative}")
        if kind == "gtf" and path.suffix == ".gz":
            raise ValueError("Current Katana requires an uncompressed GTF; archive compressed sources as auxiliary and export the genuine decompressed GTF")
        if kind == "rmats":
            match = re.fullmatch(r"(SE|RI|A3SS|A5SS|MXE)\.MATS\.(JC|JCEC)\.txt", path.name)
            if not match or item.get("comparison_id") not in comparison_ids:
                raise ValueError(f"rMATS artifact requires comparison_id and a native event filename: {path}")
            item.update(event_type=match[1], mode=match[2])
            if relative != f"rmats/{item['comparison_id']}/rmats_post/{path.name}":
                raise ValueError(f"rMATS destination must match the current Katana layout: {relative}")
        if kind in TABLE_COLUMNS or kind in {"rmats", "suppa"}:
            item.update(table_info(path, kind, samples, item.get("event_type")))
        elif kind == "auxiliary" and PurePosixPath(relative).name.startswith(("perGeneQValue.", "DEXSeqResults.")):
            # The workflow may classify native supporting tables as auxiliary;
            # validate their scientific table contract without changing bytes.
            inferred_kind = "dexseq" if relative.endswith(".csv") else "dtu"
            item.update(table_info(path, inferred_kind, samples))
        elif kind in {"gtf", "genome_fasta", "transcript_fasta"}:
            opener = gzip.open if path.suffix == ".gz" else open
            with opener(path, "rt", encoding="utf-8") as handle:
                first = next((line for line in handle if line.strip() and not line.startswith("#")), "")
            if kind == "gtf" and len(first.rstrip().split("\t")) != 9:
                raise ValueError(f"Reference is not a 9-column GTF: {path}")
            if kind != "gtf" and not first.startswith(">"):
                raise ValueError(f"Reference does not have a FASTA header: {path}")
        checked.append(item)
    available = [item for item in checked if item.get("status", "available").startswith("available")]

    def require_available(kind: str, *, comparison_id: str | None = None, sample_id: str | None = None, filename: str | None = None) -> None:
        def in_scope(item: dict) -> bool:
            return (comparison_id is None or item.get("comparison_id") == comparison_id) and (sample_id is None or item.get("sample_id") == sample_id)
        exists = any(in_scope(item) and (item["kind"] == kind if filename is None else PurePosixPath(item["relative_path"]).name == filename) for item in available)
        if exists:
            return
        declared = any(item["kind"] == kind and item.get("status") in {"failed", "not_computed"}
                       and (item.get("comparison_id") in {None, comparison_id}) and (item.get("sample_id") in {None, sample_id}) for item in checked)
        if allow_partial and declared:
            return
        scope = f" for {comparison_id or sample_id}" if comparison_id or sample_id else ""
        raise ValueError(f"{profile.capitalize()} bundle missing {filename or kind}{scope}; intentional omissions require --allow-partial and a scoped unavailable source record")

    for comparison_id in comparison_ids:
        observed = {(item.get("event_type"), item.get("mode")) for item in available if item["kind"] == "rmats" and item.get("comparison_id") == comparison_id}
        if observed != {(event, mode) for event in EVENTS for mode in MODES}:
            raise ValueError(f"Incomplete rMATS JC/JCEC event set for {comparison_id}")
        require_available("deg", comparison_id=comparison_id)
    require_available("normalized_counts")
    if profile == "full":
        for kind in ("annotation", "gtf", "genome_fasta", "transcript_fasta"):
            require_available(kind)
        for sample_id in sample_ids:
            for kind in ("quant", "bam", "bai"):
                require_available(kind, sample_id=sample_id)
        for comparison_id in comparison_ids:
            for kind in ("dexseq", "dtu", "stager", "suppa"):
                require_available(kind, comparison_id=comparison_id)
            for kind, filename in (("dexseq", f"DEXSeqResults.{comparison_id}.csv"), ("dexseq", f"perGeneQValue.{comparison_id}.csv"),
                                   ("dtu", f"DEXSeqResults.{comparison_id}.tsv"), ("dtu", f"perGeneQValue.{comparison_id}.tsv")):
                require_available(kind, comparison_id=comparison_id, filename=filename)
    output_dir.mkdir(parents=True)
    exported = []
    for item in checked:
        if not item.get("status", "available").startswith("available"):
            exported.append({key: value for key, value in item.items() if key != "source"})
            continue
        destination = inside(output_dir, item["relative_path"])
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(item["source"], destination)
        source_hash = sha256(Path(item["source"]))
        result = record(output_dir, item["relative_path"], item["kind"], **{key: value for key, value in item.items() if key not in {"kind", "relative_path", "source"}})
        if result["sha256"] != source_hash:
            raise ValueError(f"Source/copy checksum mismatch: {destination}")
        result["source_sha256"] = source_hash
        exported.append(result)
    # These normalized metadata files are adapters; original source sheets are also archived unchanged.
    sample_rows = [{"sample": row["sample_id"], **{key: value for key, value in row.items() if key != "sample_id"}} for row in samples]
    write_csv(output_dir / "pipeline_info/samplesheet.valid.csv", sample_rows, ["sample", "condition"])
    write_csv(output_dir / "pipeline_info/contrastsheet.valid.csv", contrasts, ["contrast", "treatment", "control"])
    for relative, source in (("provenance/source_samplesheet" + samples_path.suffix, samples_path), ("provenance/source_contrastsheet" + contrasts_path.suffix, contrasts_path)):
        destination = inside(output_dir, relative)
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, destination)
        exported.append(record(output_dir, relative, "provenance"))
    for relative in ("pipeline_info/samplesheet.valid.csv", "pipeline_info/contrastsheet.valid.csv"):
        exported.append(record(output_dir, relative, "metadata"))
    source_map_destination = output_dir / "provenance/source_map.json"
    shutil.copyfile(source_map_path, source_map_destination)
    exported.append(record(output_dir, "provenance/source_map.json", "provenance"))
    bam_paths = {item["sample_id"]: item["relative_path"] for item in exported if item["kind"] == "bam" and item.get("sample_id")}
    for condition, sample_order in sorted(order.items()):
        relative = f"rmats/bamlist/{condition}_bamlist.txt"
        path = inside(output_dir, relative)
        path.parent.mkdir(parents=True, exist_ok=True)
        # Without downloaded BAMs, names retain sample identity for PSI/rMATS import.
        path.write_text(",".join(bam_paths.get(sample_id, f"bam/{sample_id}.bam") for sample_id in sample_order) + "\n", encoding="utf-8")
        exported.append(record(output_dir, relative, "rmats_bamlist", condition=condition, sample_ids=sample_order, generated=True))
    for relative in make_matrices(output_dir, exported, comparisons, order):
        exported.append(record(output_dir, relative, "matrix", derived_from="unmodified rMATS per-sample vectors"))
    manifest = {"schema": SCHEMA, "created_at": datetime.now(timezone.utc).isoformat(), "project_name": data.get("project_name", "Katana project"),
                "profile": profile, "allow_partial": allow_partial, "status": "partial" if any(item["status"] in {"failed", "not_computed"} for item in exported) else "complete",
                "analysis_purpose": data.get("analysis_purpose"), "design": data.get("design"), "quantification_method": data.get("quantification_method"),
                "reference": data.get("reference", {}), "software_versions": data.get("software_versions", {}), "samples": samples,
                "comparisons": comparisons, "rmats_sample_order": order, "artifacts": exported,
                "import_instructions": "Run katana_prepare_import.py --bundle <downloaded-directory> on the destination computer before selecting this directory in Katana."}
    (output_dir / "manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    sums = [f"{item['sha256']}  {item['relative_path']}\n" for item in exported if item.get("sha256")]
    sums.append(f"{sha256(output_dir / 'manifest.json')}  manifest.json\n")
    (output_dir / "SHA256SUMS").write_text("".join(sums), encoding="utf-8")
    return manifest


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-map", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--samples", type=Path, required=True)
    parser.add_argument("--contrasts", type=Path, required=True)
    parser.add_argument("--profile", choices=("core", "full"), default="core")
    parser.add_argument("--allow-partial", action="store_true", help="Deliberately export scoped failed/not-computed optional branches; never fabricate their outputs")
    args = parser.parse_args(argv)
    try:
        result = export_bundle(args.source_map, args.output_dir, args.samples, args.contrasts, args.profile, args.allow_partial)
    except (ValueError, OSError, KeyError, json.JSONDecodeError) as error:
        parser.exit(2, f"Bundle export failed: {error}\n")
    print(json.dumps({"status": result["status"], "artifacts": len(result["artifacts"]), "bundle": str(args.output_dir.resolve())}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
