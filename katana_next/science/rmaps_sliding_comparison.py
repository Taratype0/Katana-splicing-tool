"""Event-level RNA-map audits; author source and unmodified outputs stay separate."""
from __future__ import annotations

import json
import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

from katana_cosplice.motif_scanner import IndexedFasta
from katana_next.science.se_rbp import _blocks, sequence_from_blocks

REGION_SPECS = (
    ("UpstreamExon_3prime", "UpstreamExon", "negative", 50),
    ("UpstreamExonIntron", "UpstreamExonIntron", "positive", 200),
    ("UpstreamIntron", "UpstreamIntron", "negative", 200),
    ("TargetExon_5prime", "TargetExon", "positive", 50),
    ("TargetExon-3prime", "TargetExon", "negative", 50),
    ("DownstreamIntron", "DownstreamIntron", "positive", 200),
    ("DownstreamExonIntron", "DownstreamExonIntron", "negative", 200),
    ("DownstreamExon_5prime", "DownstreamExon", "positive", 50),
)


def window_design():
    rows = []
    for index, (name, fasta_region, side, span) in enumerate(REGION_SPECS):
        for position in range(span):
            a = position if side == "positive" else position-span-49
            b = a+50
            intron = "Intron" in name
            mask = not intron or (a >= 6 if side == "positive" else b <= -20)
            rows.append(dict(region_index=index, region=name, fasta_region=fasta_region,
                side=side, native_position=position, window_start_RNA_relative=a,
                window_end_RNA_relative=b, window_center=(a+b)/2,
                canonical_mask_eligible=mask))
    return pd.DataFrame(rows)


def covered_windows(sequence, pattern, side, span):
    """Literal 50nt windows, all overlapping words, covered-base fraction."""
    covered = np.zeros(len(sequence), dtype=np.uint8)
    for hit in re.finditer("(?=("+pattern+"))", sequence):
        covered[hit.start(1):hit.end(1)] = 1
    starts = np.arange(span) if side == "positive" else len(sequence)-span-49+np.arange(span)
    stops = starts+50
    result = np.full(span, np.nan)
    good = (starts >= 0) & (stops <= len(sequence))
    valid_bases = np.array([base in "ACGT" for base in sequence], dtype=np.uint8)
    valid_sum = np.r_[0, np.cumsum(valid_bases)]
    locations = np.flatnonzero(good)
    good[locations] &= valid_sum[stops[locations]]-valid_sum[starts[locations]] == 50
    cumulative = np.r_[0, np.cumsum(covered)]
    result[good] = (cumulative[stops[good]]-cumulative[starts[good]])/50
    return result


def mw_windows(first, second):
    """Fixed greater alternative, asymptotic tied-rank MW; NA stays unavailable."""
    nfirst, nsecond = np.isfinite(first).sum(axis=0), np.isfinite(second).sum(axis=0)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        p = np.asarray(mannwhitneyu(first, second, axis=0, alternative="greater",
            method="asymptotic", use_continuity=True, nan_policy="omit").pvalue, float)
    p[(nfirst < 20) | (nsecond < 20)] = np.nan
    mean1 = np.divide(np.nansum(first, axis=0), nfirst, out=np.full(first.shape[1], np.nan), where=nfirst > 0)
    mean2 = np.divide(np.nansum(second, axis=0), nsecond, out=np.full(second.shape[1], np.nan), where=nsecond > 0)
    return p, mean1, mean2, nfirst, nsecond


def read_native_counts(path, n_events):
    table = pd.read_csv(path, sep="\t")
    arrays = np.asarray([json.loads(value) for value in table["values"]], dtype=float).T
    if arrays.shape != (n_events, 1000):
        raise ValueError("Native count distributions do not bind the frozen event rows")
    return arrays


def author_intervals(row):
    """Exactly the author's unmodified extraction intervals, in RNA orientation."""
    up, exon, down = _blocks(row)
    before = lambda block: (block[0]-250, block[0]) if row.strand == "+" else (block[1], block[1]+250)
    after = lambda block: (block[1], block[1]+250) if row.strand == "+" else (block[0]-250, block[0])
    return dict(zip(("UpstreamExon", "UpstreamExonIntron", "UpstreamIntron", "TargetExon",
                     "DownstreamIntron", "DownstreamExonIntron", "DownstreamExon"),
                    (up, after(up), before(exon), exon, after(exon), before(down), down)))


def load_and_verify_native_sequences(case, fasta_path):
    """All author-produced strings independently reloaded from the real genome."""
    case = Path(case)
    fasta = IndexedFasta(Path(fasta_path))
    groups, sequences, checks = {}, {}, 0
    geometric_audit = []
    with fasta.fasta.open("rb") as genome_stream:
        for group in ("up", "dn", "bg"):
            frame = pd.read_parquet(case/"input"/(group+".events.parquet"))
            groups[group] = frame
            sequences[group] = {}
            for region in dict.fromkeys(spec[1] for spec in REGION_SPECS):
                path = case/"author_native/fasta"/(group+"."+region+".fasta")
                strings = []
                with path.open() as stream:
                    for row in frame.itertuples(index=False):
                        header, sequence = stream.readline().strip(), stream.readline().strip()
                        if not header.startswith(">") or not sequence:
                            raise ValueError("Native FASTA event row missing")
                        coords = [str(row.chr), str(row.strand), *[str(int(getattr(row, key))) for key in
                            ("exonStart_0base", "exonEnd", "upstreamES", "upstreamEE", "downstreamES", "downstreamEE")]]
                        if header != ">"+":".join(coords):
                            raise ValueError("Native FASTA header no longer matches the ordered event")
                        interval = author_intervals(row)[region]
                        actual = sequence_from_blocks(fasta, genome_stream, row.chr, row.strand, [interval])
                        if actual != sequence:
                            raise ValueError("Native FASTA silently padded or used a different genome/strand")
                        strings.append(sequence)
                        checks += 1
                    if stream.readline():
                        raise ValueError("Unexpected extra native FASTA events")
                sequences[group][region] = strings
            for row in frame.itertuples(index=False):
                up, exon, down = _blocks(row)
                upstream_len = exon[0]-up[1] if row.strand == "+" else up[0]-exon[1]
                downstream_len = down[0]-exon[1] if row.strand == "+" else exon[0]-down[1]
                geometric_audit.append(dict(event_id=row.event_id, group=group, gene_id=row.gene_id,
                    cassette_length=exon[1]-exon[0], upstream_intron_length=upstream_len,
                    downstream_intron_length=downstream_len,
                    author_250nt_intron_extract_can_cross_neighbouring_exon=bool(min(upstream_len, downstream_len) < 250)))
    return groups, sequences, pd.DataFrame(geometric_audit), checks
