"""Read-only experiment metadata adapter around existing upstream samplesheets.

nf-core owns read merging and processing. This module validates declared identity;
it neither reconstructs biological replication nor generates comparison jobs.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from copy import deepcopy
import csv
import ntpath
from pathlib import Path, PureWindowsPath
import re
from urllib.parse import urlsplit

from .validation import load_document, validate_document

SAMPLE_IDENTITY_VERSION = "1.0.0"
_ID = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.:-]*$")
_FORMATS = {"fastq", "bam", "sample_ids"}
_ASSAYS = {"bulk_rna", "single_cell_rna", "single_nucleus_rna"}
_STRANDS = {"auto", "forward", "reverse", "unstranded"}
_META_REQUIRED = {"sample", "condition", "biological_unit_id", "donor_id"}
_META_OPTIONAL = {"display_name", "batch", "cell_type", "technical_replicate_of"}


def _issue(report, code, path, message, severity="error"):
    report["issues"].append({"code": code, "path": path, "message": message, "severity": severity})


def _report(input_format, assay_mode):
    return {"sample_identity_version": SAMPLE_IDENTITY_VERSION,
            "scope": "samplesheet_and_experiment_metadata_only", "valid": False,
            "input_format": input_format, "assay_mode": assay_mode,
            "execution_authorized": False, "scientific_tools_executed": 0,
            "merges_executed": 0, "contrasts_generated": 0,
            "upstream_processing_owner": "nf-core", "method_eligibility_checked": False,
            "independent_replicates_inferred": False, "samples": [],
            "upstream_groups": [], "counts": {}, "related_samples": [], "related_donors": [], "issues": []}


def _finish(report):
    report["valid"] = not any(x["severity"] == "error" for x in report["issues"])
    return report


def _rows_shape(rows, name, report):
    if not isinstance(rows, list) or any(not isinstance(row, dict) or any(
            not isinstance(k, str) or not isinstance(v, str) for k, v in row.items()) for row in rows):
        _issue(report, "INVALID_ROWS", "/" + name, "Expected a list of CSV row dictionaries with string keys and values.")
        return False
    if not rows:
        _issue(report, "EMPTY_SAMPLESHEET", "/" + name, "At least one sample row is required.")
        return False
    return True


def _identity(value, report, path, *, nullable=False):
    if nullable and value == "":
        return None
    if not _ID.fullmatch(value):
        _issue(report, "INVALID_SAMPLE_ID", path, "Use an explicit stable ID without whitespace; display names are separate.")
    return value


def _locator_key(value):
    """Compare declared locators without opening files or guessing sample names."""
    try:
        uri = urlsplit(value)
    except ValueError:
        return ("malformed", value)
    if "://" in value and uri.scheme:
        return ("uri", uri.scheme.lower(), uri.netloc.lower(), uri.path, uri.query, uri.fragment)
    win = PureWindowsPath(value)
    if win.drive or "\\" in value:
        return ("windows", ntpath.normpath(value).casefold())
    # Lexically normalise relative paths, retaining POSIX case and avoiding stat.
    import posixpath
    return ("posix", posixpath.normpath(value))


def _technical_links(samples, report):
    for sid, sample in samples.items():
        parent_id = sample["technical_replicate_of"]
        if parent_id is None:
            continue
        path = "/metadata/" + sid + "/technical_replicate_of"
        if parent_id == sid:
            _issue(report, "TECHNICAL_REPLICATE_CYCLE", path, "A technical replicate cannot refer to itself.")
        elif parent_id not in samples:
            _issue(report, "UNKNOWN_TECHNICAL_PARENT", path, "Technical replicate parent must be an existing sample ID.")
        else:
            parent = samples[parent_id]
            for key in ["group", "biological_unit_id", "donor_id", "cell_type"]:
                if sample[key] != parent[key]:
                    _issue(report, "TECHNICAL_REPLICATE_CONFLICT", path,
                           f"Technical replicas must retain the same condition, biological unit, donor and cell type ({key} disagrees).")
            if sample["biological_unit_id"] is None:
                _issue(report, "TECHNICAL_UNIT_UNKNOWN", path,
                       "Technical replication requires an explicitly identified shared biological unit.")
    # Iterative graph walk avoids recursion limits on large sidecars.
    finished = set()
    for start in samples:
        seen, current = set(), start
        while current in samples and current not in finished:
            if current in seen:
                _issue(report, "TECHNICAL_REPLICATE_CYCLE", "/metadata/" + start,
                       "Technical replicate relationships must be acyclic.")
                break
            seen.add(current)
            current = samples[current]["technical_replicate_of"]
        finished.update(seen)


def validate_sample_rows(upstream_rows, metadata_rows, *, input_format="fastq", assay_mode="bulk_rna"):
    """Validate copied string rows; never mutate records, merge reads or run science.

    Header checks for empty CSVs and duplicate header names are additionally made
    by validate_samplesheets. Row identity checks apply equally to API callers.
    """
    report = _report(input_format, assay_mode)
    if not isinstance(input_format, str) or input_format not in _FORMATS:
        _issue(report, "UNSUPPORTED_INPUT_FORMAT", "/input_format", "Select fastq, bam or processed sample_ids explicitly.")
    if not isinstance(assay_mode, str) or assay_mode not in _ASSAYS:
        _issue(report, "UNSUPPORTED_ASSAY", "/assay_mode", "Select bulk_rna, single_cell_rna or single_nucleus_rna explicitly.")
    elif assay_mode != "bulk_rna" and input_format == "fastq":
        _issue(report, "RAW_SINGLE_CELL_UNSUPPORTED", "/input_format", "Raw scRNA/snRNA FASTQ processing is outside the bulk upstream scope.")
    upstream_ok = _rows_shape(upstream_rows, "upstream", report)
    metadata_ok = _rows_shape(metadata_rows, "metadata", report)
    if not upstream_ok or not metadata_ok or report["issues"]:
        return _finish(report)
    samples = {}
    for i, row in enumerate(metadata_rows):
        path = f"/metadata/{i}"
        if not _META_REQUIRED.issubset(row):
            _issue(report, "MISSING_METADATA_COLUMN", path, "Metadata requires sample, condition, biological_unit_id and donor_id columns (unknown identity stays blank).")
            continue
        sid = _identity(row["sample"], report, path + "/sample")
        if sid in samples:
            _issue(report, "DUPLICATE_METADATA_SAMPLE", path, "A canonical sample ID occurs more than once in metadata; lane rows belong only in the upstream sheet.")
            continue
        if not row["condition"].strip() or row["condition"] != row["condition"].strip():
            _issue(report, "MISSING_OR_PADDED_CONDITION", path + "/condition", "Record an explicit nonempty condition without leading/trailing whitespace.")
        unit = _identity(row["biological_unit_id"], report, path + "/biological_unit_id", nullable=True)
        donor = _identity(row["donor_id"], report, path + "/donor_id", nullable=True)
        parent = _identity(row.get("technical_replicate_of", ""), report, path + "/technical_replicate_of", nullable=True)
        if unit is None:
            _issue(report, "UNKNOWN_BIOLOGICAL_UNIT", path + "/biological_unit_id", "Biological unit is unknown; no independent replicate count or analysis eligibility is inferred.", "warning")
        if donor is None:
            _issue(report, "UNKNOWN_DONOR", path + "/donor_id", "Donor is unspecified; donor-level inference requires explicit identity and later method/design qualification.", "warning")
        if assay_mode != "bulk_rna" and not row.get("cell_type", "").strip():
            _issue(report, "CELL_TYPE_UNSPECIFIED", path + "/cell_type", "Cell type is unspecified; cell-type/donor aggregation must be qualified before pseudobulk inference.", "warning")
        extra = set(row) - _META_REQUIRED - _META_OPTIONAL
        covariates = {k.removeprefix("covariate_"): v for k, v in row.items() if k.startswith("covariate_")}
        if extra - {k for k in extra if k.startswith("covariate_")}:
            _issue(report, "UNMAPPED_METADATA_COLUMNS", path, "Additional metadata columns are preserved as metadata; they do not silently change a statistical model.", "warning")
        sample = {"sample_id": sid, "display_name": row.get("display_name", "") or sid,
                  "group": row["condition"], "condition": row["condition"],
                  "biological_unit_id": unit, "donor_id": donor,
                  "cell_type": row.get("cell_type", "") or None, "batch": row.get("batch", "") or None,
                  "technical_replicate_of": parent, "source_assay": assay_mode,
                  "covariates": covariates, "source_metadata": dict(row)}
        samples[sid] = sample
    labels = Counter(s["display_name"] for s in samples.values())
    for label, count in labels.items():
        if count > 1:
            _issue(report, "DUPLICATE_DISPLAY_NAME", "/metadata", f"Display name is shared by {count} sample IDs ({label}); these records are kept separate.", "warning")
    unit_donors = defaultdict(set)
    for sample in samples.values():
        if sample["biological_unit_id"] and sample["donor_id"]:
            unit_donors[sample["biological_unit_id"]].add(sample["donor_id"])
    for unit, donors in unit_donors.items():
        if len(donors) > 1:
            _issue(report, "BIOLOGICAL_UNIT_DONOR_CONFLICT", "/metadata", f"One biological unit ({unit}) has conflicting donor IDs; pooled donors need a separate explicit adapter.")
    _technical_links(samples, report)

    grouped, layouts, strands, paths, exact_rows = {}, defaultdict(set), defaultdict(set), {}, set()
    for i, row in enumerate(upstream_rows):
        path = f"/upstream/{i}"
        required = {"sample", "fastq_1", "fastq_2", "strandedness"} if input_format == "fastq" else {"sample"}
        if not required.issubset(row):
            _issue(report, "MISSING_UPSTREAM_COLUMN", path, f"Missing required {input_format} samplesheet columns: {', '.join(sorted(required - set(row)))}.")
            continue
        sid = _identity(row["sample"], report, path + "/sample")
        if sid not in samples:
            _issue(report, "MISSING_SAMPLE_METADATA", path + "/sample", f"Upstream sample {sid} has no canonical experiment metadata row.")
        row_key = tuple(sorted(row.items()))
        if row_key in exact_rows:
            _issue(report, "DUPLICATE_UPSTREAM_ROW", path, "The same upstream row is repeated; duplicate lane/file rows must not be counted or merged twice.")
        exact_rows.add(row_key)
        grouped.setdefault(sid, []).append(dict(row))
        sample = samples.get(sid)
        if sample:
            for column, field in [("condition", "group"), ("group", "group"), ("biological_unit_id", "biological_unit_id"), ("donor_id", "donor_id"), ("cell_type", "cell_type")]:
                if row.get(column) and row[column] != sample[field]:
                    _issue(report, "UPSTREAM_METADATA_CONFLICT", path + "/" + column, "Upstream metadata and experiment sidecar disagree; neither is silently overwritten.")
        files = []
        if input_format == "fastq":
            if not row["fastq_1"].strip():
                _issue(report, "READ1_REQUIRED", path + "/fastq_1", "R1 FASTQ path is required; R2 may be blank only for single-end reads.")
            files = [("fastq_1", row["fastq_1"]), ("fastq_2", row["fastq_2"])]
            layouts[sid].add("paired" if row["fastq_2"] else "single")
        elif input_format == "bam":
            genome, alias = row.get("genome_bam", ""), row.get("bam", "")
            if genome and alias and _locator_key(genome) != _locator_key(alias):
                _issue(report, "AMBIGUOUS_BAM_COLUMN", path, "genome_bam and bam identify different files; select one explicit genome-aligned input.")
            if not (genome or alias):
                _issue(report, "GENOME_BAM_REQUIRED", path, "A genome_bam or bam path is required; transcriptome_bam is not a genome BAM substitute.")
            files = [("genome_bam", genome or alias), ("bai", row.get("bai", ""))]
            layouts[sid].add("bam")
        else:
            layouts[sid].add("sample_ids")
            if len(grouped[sid]) > 1:
                _issue(report, "DUPLICATE_PROCESSED_SAMPLE", path, "Processed sample-ID mapping is canonical; repeated rows are not sequencing lanes.")
        if "strandedness" in row:
            value = row["strandedness"]
            if value not in _STRANDS:
                _issue(report, "INVALID_STRANDEDNESS", path + "/strandedness", "Use auto, forward, reverse or unstranded explicitly.")
            strands[sid].add(value)
        for column, value in files:
            if not value:
                continue
            if value != value.strip():
                _issue(report, "PADDED_INPUT_PATH", path + "/" + column, "Input paths must not contain leading/trailing whitespace.")
            key = _locator_key(value)
            if key[0] == "malformed":
                _issue(report, "INVALID_INPUT_LOCATOR", path + "/" + column, "Malformed input URI; the adapter never fetches or repairs it.")
            if key in paths:
                _issue(report, "DUPLICATE_INPUT_FILE", path + "/" + column, f"An input locator is reused (first at {paths[key]}); mates and samples must not reuse a file.")
            else:
                paths[key] = path + "/" + column
        # Optional run/lane labels are provenance. Lane L001 in different runs is legal.
    for sid, rows in grouped.items():
        if input_format == "bam" and len(rows) > 1:
            _issue(report, "MULTIPLE_BAM_ROWS", "/upstream/" + sid, "Multiple BAM records for one sample require a qualified upstream adapter; FASTQ multi-lane concatenation rules do not authorise BAM merging.")
        if len(layouts[sid]) > 1:
            _issue(report, "MIXED_READ_LAYOUT", "/upstream/" + sid, "One sample cannot mix single-end and paired-end lane rows in this adapter.")
        if len(strands[sid]) > 1:
            _issue(report, "MIXED_STRANDEDNESS", "/upstream/" + sid, "Strandedness differs across rows for one upstream sample.")
        if len({r.get("library_id") for r in rows if r.get("library_id")}) > 1:
            _issue(report, "MULTIPLE_LIBRARY_PREPARATIONS", "/upstream/" + sid, "Multiple declared library preparations share one upstream sample ID; nf-core owns grouping, and they are not independent biological replicates.", "warning")
        report["upstream_groups"].append({"sample_id": sid, "rows": rows, "row_count": len(rows),
                                           "layout": sorted(layouts[sid])[0] if len(layouts[sid]) == 1 else "mixed",
                                           "merge_owner": "nf-core", "merge_performed": False})
    for sid in sorted(samples.keys() - grouped.keys()):
        _issue(report, "UNMAPPED_METADATA_SAMPLE", "/metadata/" + sid, "Experiment sample is absent from the declared upstream sheet; choose the matching sheet instead of silently dropping it.")
    related = defaultdict(list)
    for sample in samples.values():
        if sample["biological_unit_id"]:
            related[sample["biological_unit_id"]].append(sample["sample_id"])
    for unit, ids in related.items():
        if len(ids) > 1:
            report["related_samples"].append({"biological_unit_id": unit, "sample_ids": ids})
            _issue(report, "REPEATED_BIOLOGICAL_UNIT", "/metadata", f"Biological unit {unit} occurs in multiple sample records; these are related measurements, not additional independent replicates. T11 must qualify their design.", "warning")
    donors = defaultdict(list)
    for sample in samples.values():
        if sample["donor_id"]:
            donors[sample["donor_id"]].append(sample["sample_id"])
    for donor, ids in donors.items():
        units = sorted({samples[sid]["biological_unit_id"] for sid in ids if samples[sid]["biological_unit_id"]})
        if len(ids) > 1:
            report["related_donors"].append({"donor_id": donor, "sample_ids": ids, "biological_unit_ids": units})
        if len(units) > 1:
            _issue(report, "SHARED_DONOR_ACROSS_UNITS", "/metadata", f"Donor {donor} appears under multiple declared units; distinct unit IDs do not establish independent donors.", "warning")
    report["samples"] = list(samples.values())
    report["counts"] = {"sample_count": len(samples), "upstream_row_count": len(upstream_rows),
                         "biological_unit_count": len(related),
                         "donor_count": len({s["donor_id"] for s in samples.values() if s["donor_id"]}),
                         "technical_replicate_count": sum(s["technical_replicate_of"] is not None for s in samples.values()),
                         "unknown_biological_unit_count": sum(s["biological_unit_id"] is None for s in samples.values())}
    return _finish(report)


def _read_csv(path, required, report, name):
    try:
        with Path(path).open("r", encoding="utf-8-sig", newline="") as stream:
            reader = csv.reader(stream, strict=True)
            header = next(reader, None)
            if header is None:
                _issue(report, "EMPTY_SAMPLESHEET", "/" + name, "CSV has no header or sample rows.")
                return []
            if len(header) != len(set(header)) or any(not value or value != value.strip() for value in header):
                _issue(report, "INVALID_CSV_HEADER", "/" + name, "CSV headers must be nonempty, unpadded and unique.")
                return []
            if not required.issubset(header):
                _issue(report, "MISSING_CSV_COLUMN", "/" + name, "Missing required CSV columns: " + ", ".join(sorted(required - set(header))))
                return []
            rows = []
            for line, cells in enumerate(reader, start=2):
                if not cells or all(not value for value in cells):
                    continue
                if len(cells) != len(header):
                    _issue(report, "CSV_ROW_WIDTH", f"/{name}/{line}", "CSV row has a different number of fields from the header; use CSV quoting for commas.")
                    continue
                rows.append(dict(zip(header, cells)))
            return rows
    except (OSError, UnicodeError, csv.Error) as error:
        _issue(report, "CSV_LOAD", "/" + name, str(error))
        return []


def check_project_sample_mapping(project_document, sample_report):
    """Check an existing T09 project's declared identities without rewriting it."""
    report = deepcopy(sample_report)
    report["project_mapping_checked"] = True
    validation = validate_document(project_document)
    if not validation["valid"] or project_document.get("kind") != "project":
        _issue(report, "PROJECT_CONTRACT_INVALID", "/project", "Project must pass the existing T09 project contract before identity comparison.")
        return _finish(report)
    if not sample_report.get("valid"):
        _issue(report, "SAMPLE_MAPPING_INVALID", "/project", "Upstream/metadata mapping must be valid before it is matched to a project.")
        return _finish(report)
    expected = {s["sample_id"]: s for s in report["samples"]}
    actual = {s["sample_id"]: s for s in project_document["samples"]}
    if expected.keys() != actual.keys():
        _issue(report, "PROJECT_SAMPLE_SCOPE", "/project/samples", "Project and experiment metadata must contain the same explicit sample IDs.")
    if project_document["assay_mode"] != report["assay_mode"]:
        _issue(report, "PROJECT_ASSAY_MISMATCH", "/project/assay_mode", "Project and metadata source assays disagree.")
    for sid in sorted(expected.keys() & actual.keys()):
        for field in ["group", "biological_unit_id", "donor_id", "cell_type"]:
            if expected[sid][field] != actual[sid][field]:
                _issue(report, "PROJECT_SAMPLE_IDENTITY", "/project/samples/" + sid + "/" + field,
                       "Project sample identity conflicts with experiment metadata; no field is overwritten.")
    return _finish(report)


def validate_samplesheets(upstream_path, metadata_path, *, input_format="fastq", assay_mode="bulk_rna", project_path=None):
    """Read two CSVs and optional T09 project. Never read FASTQ/BAM data."""
    precheck = _report(input_format, assay_mode)
    required = {"sample", "fastq_1", "fastq_2", "strandedness"} if input_format == "fastq" else {"sample"}
    upstream = _read_csv(upstream_path, required, precheck, "upstream")
    metadata = _read_csv(metadata_path, _META_REQUIRED, precheck, "metadata")
    if precheck["issues"]:
        return _finish(precheck)
    report = validate_sample_rows(upstream, metadata, input_format=input_format, assay_mode=assay_mode)
    if project_path is not None:
        try:
            project = load_document(project_path)
        except (ValueError, OSError) as error:
            _issue(report, "PROJECT_LOAD", "/project", str(error))
        else:
            report = check_project_sample_mapping(project, report)
    return _finish(report)
