"""Audit original SpliceSeq event definitions, independently of measured PSI."""
import csv
import hashlib
import json
from pathlib import Path


TABLE_COLUMNS = {
    "graph": ("graph_id", "species", "symbol", "chr", "strand", "start", "stop"),
    "as_ref": ("as_id", "graph_id", "splice_type", "start_ex", "end_ex",
               "novel_splice", "exons", "exons_hash", "avg_exp_pct"),
    "as_ref_exon": ("as_id", "exon_id"),
    "exon": ("graph_id", "exon_id", "chr_start", "chr_stop", "exon_name",
             "graph_start", "graph_stop", "model_occur", "nucleic_seq"),
    "version": ("current",),
}


def canonical_sha(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, ensure_ascii=True,
                                    separators=(",", ":"), allow_nan=False).encode()).hexdigest()


def read_table(folder, table, *, native_dump):
    columns = TABLE_COLUMNS[table]
    path = Path(folder) / (table + (".txt" if native_dump else ".tsv"))
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t", quoting=csv.QUOTE_NONE)
        if not native_dump and tuple(next(reader, ())) != columns:
            raise ValueError("Unexpected original native definition columns: " + table)
        for line, fields in enumerate(reader, 1 if native_dump else 2):
            if len(fields) != len(columns):
                raise ValueError(f"Malformed native definition {table}:{line}")
            yield {key: None if value in {"\\N", "NULL"} else value
                   for key, value in zip(columns, fields)}


def read_identity_map(path):
    with Path(path).open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    fields = {"event_key", "native_as_id", "original_as_id",
              "composite_event", "exact_full_path_identity", "reasons"}
    if not rows or not fields <= set(rows[0]):
        raise ValueError("Complete original identity mapping required")
    if len({row["event_key"] for row in rows}) != len(rows):
        raise ValueError("Duplicate original event identity")
    ids = [row["native_as_id"] for row in rows]
    if len(set(ids)) != len(ids) or any(not value.isdigit() for value in ids):
        raise ValueError("Unique integer original native event IDs required")
    for row in rows:
        if row["exact_full_path_identity"] not in {"True", "False"}:
            raise ValueError("Identity flags must remain explicit booleans")
    return rows


def build_native_definition_snapshot(folder, identity_rows, *, native_dump=False):
    """Retain complete original paths, including previously recorded mismatches.

    Exon sequence, genomic/graph coordinates, flags, ordered native path text,
    inclusive-exon links, gene identity and reference version are checked.
    Sample-dependent avg_exp_pct is not an event definition. No PSI is read.
    """
    wanted = {row["native_as_id"] for row in identity_rows}
    events = {}
    for row in read_table(folder, "as_ref", native_dump=native_dump):
        if row["as_id"] not in wanted:
            continue
        if row["as_id"] in events:
            raise ValueError("Duplicate original native as_ref key")
        events[row["as_id"]] = {k: v for k, v in row.items() if k != "avg_exp_pct"}
    missing = wanted - set(events)
    links = {event: [] for event in wanted}
    for row in read_table(folder, "as_ref_exon", native_dump=native_dump):
        if row["as_id"] in wanted:
            if row["exon_id"] in links[row["as_id"]]:
                raise ValueError("Duplicate original inclusive-exon link")
            links[row["as_id"]].append(row["exon_id"])
    exon_ids = {x for row in events.values() for x in (row["start_ex"], row["end_ex"])
                if x is not None}
    exon_ids.update(x for xs in links.values() for x in xs if x is not None)
    exons = {}
    for row in read_table(folder, "exon", native_dump=native_dump):
        if row["exon_id"] not in exon_ids:
            continue
        if row["exon_id"] in exons:
            raise ValueError("Duplicate original native exon key")
        sequence = row["nucleic_seq"]
        exons[row["exon_id"]] = {k: v for k, v in row.items()
                                if k not in {"model_occur", "nucleic_seq"}}
        exons[row["exon_id"]]["nucleic_seq_sha256"] = (
            hashlib.sha256(sequence.encode()).hexdigest() if sequence is not None else None)
    graph_ids = {row["graph_id"] for row in events.values()}
    graphs = {}
    for row in read_table(folder, "graph", native_dump=native_dump):
        if row["graph_id"] not in graph_ids:
            continue
        if row["graph_id"] in graphs:
            raise ValueError("Duplicate original native graph key")
        graphs[row["graph_id"]] = row
    versions = [row["current"] for row in read_table(folder, "version", native_dump=native_dump)]
    if len(versions) != 1 or not versions[0]:
        raise ValueError("One explicit native reference version required")
    fingerprints = {}
    for identity in identity_rows:
        event_id = identity["native_as_id"]
        event = events.get(event_id)
        ids = sorted(set(links[event_id] + ([event["start_ex"], event["end_ex"]] if event else [])),
                     key=lambda value: "" if value is None else value)
        definition = dict(native_as_ref=event,
            native_graph=graphs.get(event["graph_id"]) if event else None,
            inclusive_exon_ids=sorted(links[event_id]),
            all_original_exons={str(x): exons.get(x) for x in ids},
            native_reference_version=versions[0])
        fingerprints[identity["event_key"]] = canonical_sha(definition)
    return dict(native_reference_version=versions[0], original_event_count=len(identity_rows),
        missing_native_as_ids=sorted(missing),
        missing_native_exon_ids=sorted(exon_ids - set(exons)),
        missing_native_graph_ids=sorted(graph_ids - set(graphs)),
        event_definition_sha256=fingerprints, measurements_opened=False)


def compare_native_definition_snapshots(baseline, current, identity_rows):
    keys = {row["event_key"] for row in identity_rows}
    for snapshot in (baseline, current):
        if set(snapshot["event_definition_sha256"]) != keys:
            raise ValueError("Definition audit must retain every original event")
    changed = [row["event_key"] for row in identity_rows
               if baseline["event_definition_sha256"][row["event_key"]]
               != current["event_definition_sha256"][row["event_key"]]]
    return dict(original_event_count=len(keys), changed_original_event_keys=changed,
        changed_original_event_count=len(changed), all_original_definitions_unchanged=not changed,
        native_reference_version_unchanged=(baseline["native_reference_version"]
                                           == current["native_reference_version"]),
        inherited_exact_identity_events=sum(row["exact_full_path_identity"] == "True"
                                            for row in identity_rows),
        inherited_identity_mismatches=sum(row["exact_full_path_identity"] == "False"
                                         for row in identity_rows),
        missing_values_or_mismatches_repaired=False, full_cohort_quantified=False,
        external_preservation_tested=False)
