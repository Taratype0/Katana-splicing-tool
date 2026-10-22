"""Whole-gene alignment evidence; never an automatic knockout classifier."""
import re
from collections import Counter

import numpy as np
import pandas as pd

from katana_cosplice.bam_io import closing_alignment


def gene_exons(gtf, query):
    rows = []
    with open(gtf, encoding='utf-8') as stream:
        for line in stream:
            if line.startswith('#') or query not in line:
                continue
            fields = line.rstrip().split('\t')
            if len(fields) != 9 or fields[2] != 'exon':
                continue
            attrs = dict(re.findall(r'(\w+) "([^"]*)"', fields[8]))
            if query not in {attrs.get('gene_id', '').split('.')[0], attrs.get('gene_id'), attrs.get('gene_name')}:
                continue
            number = re.search(r'exon_number\s+"?(\d+)', fields[8])
            rows.append(dict(gene_id=attrs['gene_id'], gene_symbol=attrs.get('gene_name', ''),
                chrom=fields[0], start=int(fields[3])-1, end=int(fields[4]), strand=fields[6],
                transcript=attrs.get('transcript_id', ''), exon_number=int(number.group(1)) if number else None))
    result = pd.DataFrame(rows)
    if result.empty or result.gene_id.nunique() != 1 or result.chrom.nunique() != 1:
        raise ValueError('Gene must resolve unambiguously to one annotated locus')
    return result


def alignment_evidence(reads, start, end, min_mapq=20, min_anchor=8):
    delta = np.zeros(end-start+1, dtype=np.int64)
    junctions = Counter()
    qc = Counter()
    for read in reads:
        qc['fetched_alignments'] += 1
        if read.flag & (4 | 256 | 512 | 2048):
            qc['excluded_flags'] += 1
            continue
        tags = dict(read.get_tags())
        if read.mapq < min_mapq or tags.get('NH', 1) != 1:
            qc['excluded_mapping'] += 1
            continue
        qc['retained_alignments'] += 1
        qc['NH_missing'] += int('NH' not in tags)
        qc['marked_duplicates'] += int(bool(read.flag & 1024))
        pos = read.pos
        cigar = read.cigar
        for i, (op, length) in enumerate(cigar):
            if op in (0, 7, 8):
                a, b = max(start, pos)-start, min(end, pos+length)-start
                if b > a:
                    delta[a] += 1
                    delta[b] -= 1
                pos += length
            elif op == 3:
                left = cigar[i-1][1] if i and cigar[i-1][0] in (0, 7, 8) else 0
                right = cigar[i+1][1] if i+1 < len(cigar) and cigar[i+1][0] in (0, 7, 8) else 0
                junctions[(pos, pos+length, 'primary')] += 1
                if not read.flag & 1024 and min(left, right) >= min_anchor:
                    junctions[(pos, pos+length, 'nonduplicate_anchored')] += 1
                pos += length
            elif op == 2:
                qc['CIGAR_deleted_bases'] += length
                pos += length
            elif op not in (1, 4, 5, 6):
                raise ValueError(f'Unsupported CIGAR operation: {op}')
    return np.cumsum(delta[:-1]), junctions, dict(qc)


def audit_gene(exons, samples, min_mapq=20, min_anchor=8, progress=None):
    import bamnostic
    start, end = int(exons.start.min()), int(exons.end.max())
    if end-start > 5_000_000:
        raise ValueError('Gene span exceeds the declared 5 Mb memory guard; use a scoped event view')
    chrom = str(exons.chrom.iloc[0])
    mask = np.zeros(end-start, dtype=bool)
    for row in exons.itertuples():
        mask[row.start-start:row.end-start] = True
    coverage, junctions, qc = [], [], []
    for i, (sample, bam) in enumerate(samples.items()):
        if progress:
            progress(i, len(samples), sample)
        with closing_alignment(bamnostic.AlignmentFile(str(bam), 'rb')) as handle:
            depth, counts, audit = alignment_evidence(handle.fetch(chrom, start, end), start, end, min_mapq, min_anchor)
        coverage.append(depth)
        qc.append(dict(sample_id=sample, **audit, mean_exon_depth=float(depth[mask].mean())))
        for (a, b, mode), count in counts.items():
            junctions.append(dict(sample_id=sample, chrom=chrom, intron_start_1based=a+1,
                intron_end_1based=b, filter=mode, alignment_count=count))
    matrix = np.array(coverage)
    denominator = matrix[:, mask].mean(axis=1)
    rows = []
    edges = sorted(set(exons.start) | set(exons.end))
    for a, b in zip(edges[:-1], edges[1:]):
        if not mask[a-start:b-start].any():
            continue
        raw = matrix[:, a-start:b-start].mean(axis=1)
        for i, sample in enumerate(samples):
            rows.append(dict(sample_id=sample, gene_id=exons.gene_id.iloc[0], chrom=chrom,
                start_1based=a+1, end_1based=b, mean_depth=raw[i],
                relative_exon_depth=raw[i]/denominator[i] if denominator[i] > 0 else np.nan))
    return pd.DataFrame(rows), pd.DataFrame(junctions, columns=['sample_id','chrom','intron_start_1based','intron_end_1based','filter','alignment_count']), pd.DataFrame(qc), matrix
