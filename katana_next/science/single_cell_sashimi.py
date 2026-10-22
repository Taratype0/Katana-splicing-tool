"""Regional BAM evidence with explicit cell-to-donor mapping and count semantics."""
from collections import defaultdict
import re

import numpy as np
import pandas as pd

from katana_next.science.single_cell import _normalize_cigar_item


def junction_coordinates(event_id):
    match = re.fullmatch(r"(.+):(\d+)-(\d+):([+.-])", str(event_id))
    if not match or int(match[3]) <= int(match[2]):
        raise ValueError(f"Sashimi needs exact 0-based half-open junction IDs (chrom:start-end:strand), not gene names: {event_id}")
    return match[1], int(match[2]), int(match[3]), match[4]


def cigar_evidence(read, left, right):
    reference = int(read.reference_start)
    blocks, junctions = [], []
    for item in getattr(read, "cigartuples", None) or getattr(read, "cigar", []) or []:
        operation, length = _normalize_cigar_item(item)
        if operation in {0, 7, 8}:
            a, b = max(left, reference), min(right, reference+length)
            if b > a:
                blocks.append((a, b))
        if operation == 3:
            junctions.append((reference, reference+length))
        if operation in {0, 2, 3, 7, 8}:
            reference += length
    return blocks, junctions


def indexed_reads(handle, chrom, left, right):
    """Respect reference span, including long introns, on the Python BAM backend."""
    if type(handle).__module__.startswith("bamnostic"):
        if not getattr(handle, "_index", None):
            raise ValueError("Indexed BAM is required for regional splicing evidence")
        offset = handle._index.query(handle.ref2tid[chrom], left, right)
        if offset is None:
            return
        intervals = getattr(handle._index.current_ref, "intervals", ())
        linear = intervals[min(left >> 14, len(intervals)-1)] if intervals else 0
        handle.seek(max(offset, linear))
        for read in handle:
            if read.reference_name != chrom or read.reference_start >= right:
                break
            span = sum(n for op, n in read.cigartuples if op in {0, 2, 3, 7, 8})
            if read.reference_start + span > left:
                yield read
    else:
        yield from handle.fetch(contig=chrom, start=left, stop=right)


def regional_evidence(samples, metadata, event_ids, *, flank=150, min_mapq=20,
                      barcode_tag="CB", umi_tag="UB", check_cancelled=lambda: None,
                      opener=None):
    coordinates = {event: junction_coordinates(event) for event in event_ids}
    chroms = {v[0] for v in coordinates.values()}
    strands = {v[3] for v in coordinates.values()}
    if len(chroms) != 1 or len(strands) != 1:
        raise ValueError("Each junction cluster must have one chromosome and one strand convention.")
    chrom, strand = next(iter(chroms)), next(iter(strands))
    left = max(0, min(v[1] for v in coordinates.values())-flank)
    right = max(v[2] for v in coordinates.values())+flank
    if right-left > 2000000:
        raise ValueError("Selected cluster spans >2 Mb; select a local junction cluster rather than unrelated events.")
    if opener is None:
        try:
            import pysam
            opener = pysam.AlignmentFile
        except ImportError:
            from bamnostic import AlignmentFile
            opener = AlignmentFile
    keys = ["donor_id", "condition", "cell_type"]
    track_cells = {tuple(k): set(g.cell_id) for k, g in metadata.groupby(keys, sort=True)}
    coverage = {key: np.zeros(right-left+1, dtype=np.int64) for key in track_cells}
    molecules, regional_cells = defaultdict(set), defaultdict(set)
    coordinate_events = defaultdict(list)
    for event, (_, a, b, _) in coordinates.items():
        coordinate_events[(a, b)].append(event)
    audits = []
    known_samples = {s.sample_id for s in samples}
    if not set(metadata.sample_id).issubset(known_samples):
        raise ValueError("Single-cell metadata references sample IDs absent from the project.")
    for sample in samples:
        cells = metadata[metadata.sample_id.eq(sample.sample_id)]
        if cells.empty:
            continue
        bam = sample.assets.get("single_cell_bam") or sample.assets.get("bam")
        if bam is None or not bam.is_file():
            raise ValueError(f"Missing single-cell BAM for {sample.sample_id}. Junction counts cannot reconstruct coverage.")
        full = cells.assay.str.lower().isin({"full_length", "smart_seq2", "smart_seq3", "smart-seq2", "smart-seq3"}).all()
        modes = set(cells.get("count_mode", pd.Series(dtype=str)).dropna().astype(str))
        if len(modes) > 1 or modes - {"umi", "read_name"}:
            raise ValueError(f"Declare one count_mode (umi/read_name) for sample {sample.sample_id}.")
        mode = next(iter(modes), "read_name" if full else "umi")
        lookup = {row.cell_id: (row.donor_id, row.condition, row.cell_type) for row in cells.itertuples()}
        raw_lookup = {}
        if "barcode" in cells:
            if cells.barcode.isna().any() or cells.barcode.duplicated().any():
                raise ValueError(f"BAM barcode mapping must be unique and nonmissing within {sample.sample_id}.")
            raw_lookup = dict(zip(cells.barcode.astype(str), cells.cell_id))
        counters = dict(sample_id=sample.sample_id, bam=str(bam), count_mode=mode,
                        fetched_reads=0, accepted_reads=0, unmatched_barcode=0, missing_molecule_id=0, strand_filtered=0)
        handle = opener(str(bam), "rb")
        try:
            try:
                reads = indexed_reads(handle, chrom, left, right)
                for i, read in enumerate(reads):
                    if i % 10000 == 0:
                        check_cancelled()
                    counters["fetched_reads"] += 1
                    if any(getattr(read, flag, False) for flag in ("is_unmapped", "is_secondary", "is_supplementary", "is_qcfail")) or getattr(read, "mapping_quality", 0) < min_mapq:
                        continue
                    try:
                        barcode = str(read.get_tag(barcode_tag))
                        cell = raw_lookup.get(barcode, barcode if barcode in lookup else f"{sample.sample_id}:{barcode}")
                    except (KeyError, AttributeError):
                        cell = next(iter(lookup)) if full and len(lookup) == 1 else ""
                    if cell not in lookup:
                        counters["unmatched_barcode"] += 1
                        continue
                    try:
                        read_strand = str(read.get_tag("XS"))
                    except (KeyError, AttributeError):
                        read_strand = "."
                    if strand != "." and read_strand != strand:
                        counters["strand_filtered"] += 1
                        continue
                    try:
                        molecule = str(read.get_tag(umi_tag)) if mode == "umi" else str(getattr(read, "query_name", "") or "")
                    except (KeyError, AttributeError):
                        molecule = ""
                    if not molecule:
                        counters["missing_molecule_id"] += 1
                        continue
                    key = lookup[cell]
                    blocks, junctions = cigar_evidence(read, left, right)
                    for a, b in blocks:
                        coverage[key][a-left] += 1
                        coverage[key][b-left] -= 1
                    for junction in junctions:
                        for event in coordinate_events[junction]:
                            molecules[(key, event)].add((sample.sample_id, cell, molecule))
                    regional_cells[key].add(cell)
                    counters["accepted_reads"] += 1
            except (ValueError, AssertionError, OSError) as exc:
                raise ValueError(f"Cannot fetch {chrom}:{left}-{right} from {bam}. Supply a matching coordinate-sorted BAM index and chromosome naming. {exc}") from exc
        finally:
            from katana_cosplice.bam_io import close_alignment
            close_alignment(handle)
        audits.append(counters)
    rows = []
    for key in track_cells:
        total = sum(len(molecules[(key, event)]) for event in event_ids)
        for event in event_ids:
            count = len(molecules[(key, event)])
            rows.append(dict(zip(keys, key), event_id=event, junction_molecules=count,
                        cluster_junction_total=total, junction_usage=count/total if total else np.nan,
                        cells_in_metadata=len(track_cells[key]), cells_with_regional_reads=len(regional_cells[key])))
        coverage[key] = np.cumsum(coverage[key][:-1])
    return (chrom, left, right, strand), coverage, pd.DataFrame(rows), pd.DataFrame(audits)
