"""Exact event identity and measured-count transport for frozen external validation."""
from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
import re

import numpy as np
import pandas as pd


def full_frozen_member_gate(matrix, members, samples, exact_identity,
                            minimum_patients=15, minimum_members=20,
                            variance_floor=1e-8):
    """Qualify a full measured module without selecting a smaller event core.

    Call on one biological scope (for example tumors), with one observation
    per patient. This is an input gate, not a preservation or significance test.
    """
    keys = pd.Index(members)
    if keys.empty or keys.has_duplicates:
        raise ValueError('Frozen membership must be nonempty and unique')
    if matrix.index.has_duplicates or matrix.columns.has_duplicates:
        raise ValueError('Duplicate event or sample matrix identity')
    if not {'sample_id', 'patient_id'} <= set(samples):
        raise ValueError('Explicit sample and biological patient identities required')
    if samples[['sample_id', 'patient_id']].isna().any().any():
        raise ValueError('Missing biological sample identity')
    if samples.sample_id.duplicated().any() or samples.patient_id.duplicated().any():
        raise ValueError('One biological observation per patient per structural scope required')
    if minimum_patients < 2 or minimum_members < 2 or variance_floor < 0:
        raise ValueError('Invalid predeclared full-module input gates')
    absent_samples = set(samples.sample_id) - set(matrix.columns)
    if absent_samples:
        raise ValueError('Declared samples missing from matrix: ' + str(sorted(absent_samples)))
    absent = keys.difference(matrix.index).tolist()
    exact = [k for k in keys if isinstance(exact_identity.get(k), (bool, np.bool_))
             and bool(exact_identity[k])]
    values = matrix.reindex(keys).loc[:, samples.sample_id].to_numpy(dtype=float)
    finite = np.isfinite(values)
    if np.any(finite & ((values < 0) | (values > 1))):
        raise ValueError('Measured PSI outside [0, 1]')
    complete = finite.all(axis=0)
    patients = samples.loc[complete, 'patient_id'].astype(str).tolist()
    variable = 0
    if complete.sum() >= 2:
        variable = int((np.var(values[:, complete], axis=1, ddof=1) > variance_floor).sum())
    reasons = []
    if len(keys) < minimum_members:
        reasons.append('original_module_below_minimum_members')
    if absent:
        reasons.append('missing_original_event_rows')
    if len(exact) != len(keys):
        reasons.append('incomplete_exact_event_identity')
    if len(patients) < minimum_patients:
        reasons.append('insufficient_full_module_complete_patients')
    if variable != len(keys):
        reasons.append('not_all_original_members_variable_in_full_cases')
    return dict(full_module_input_eligible=not reasons,
                original_members=len(keys), exact_identity_members=len(exact),
                measured_members=int(finite.any(axis=1).sum()),
                complete_case_patients=len(patients),
                variable_members_in_complete_cases=variable,
                complete_case_patient_ids=patients,
                absent_original_events=absent, reasons=reasons,
                dropped_original_members=0, imputed_cells=0,
                minimum_patients=minimum_patients, minimum_members=minimum_members,
                variance_floor=variance_floor, preservation_pvalue=None,
                preservation_fdr=None)


@dataclass(frozen=True)
class ExonInterval:
    chromosome: str
    start: int
    end: int
    strand: str

    def __post_init__(self):
        if self.start < 0 or self.end <= self.start or self.strand not in {'+', '-'}:
            raise ValueError('Expected a zero-based half-open exon interval')


def chain_blocks(lines):
    """UCSC chain target is the old build; query is the new build."""
    result = defaultdict(list)
    active = None
    for line in lines:
        fields = line.split()
        if not fields:
            continue
        if fields[0] == 'chain':
            if fields[4] != '+':
                raise ValueError('Unsupported reverse source orientation in chain')
            active = (fields[2], fields[7], int(fields[8]), fields[9])
            t = int(fields[5])
            q = int(fields[10])
        elif active is not None:
            size = int(fields[0])
            chrom, dest, qsize, orientation = active
            result[chrom].append((t, t+size, dest, q, qsize, orientation))
            if len(fields) == 3:
                t += size + int(fields[1])
                q += size + int(fields[2])
            elif len(fields) == 1:
                active = None
            else:
                raise ValueError('Malformed UCSC chain block')
    return dict(result)


def translate_interval(exon, blocks):
    hits = set()
    for a, b, dest, q, qsize, orientation in blocks.get(exon.chromosome, ()):
        if a <= exon.start < exon.end <= b:
            if orientation == '+':
                hits.add(ExonInterval(dest, q+exon.start-a, q+exon.end-a, exon.strand))
            else:
                hits.add(ExonInterval(dest, qsize-(q+exon.end-a),
                    qsize-(q+exon.start-a), '+' if exon.strand == '-' else '-'))
    return next(iter(hits)) if len(hits) == 1 else None


def normalize_exon_id(value):
    return str(value).strip().removesuffix('.0')


def splice_graph(table):
    required = {'Symbol', 'Chromosome', 'Strand', 'Exon', 'Chr_Start', 'Chr_Stop'}
    if not required <= set(table):
        raise ValueError('Official SpliceSeq graph fields missing')
    result = {}
    for row in table.itertuples(index=False):
        key = (str(row.Symbol), normalize_exon_id(row.Exon))
        a, b = sorted((int(row.Chr_Start), int(row.Chr_Stop)))
        chrom = str(row.Chromosome).removesuffix('.0')
        chrom = 'chrM' if chrom in {'M', 'MT'} else 'chr'+chrom.removeprefix('chr')
        exon = ExonInterval(chrom, a-1, b, row.Strand)
        if key in result and result[key] != exon:
            raise ValueError('Ambiguous official gene/exon graph key: ' + str(key))
        result[key] = exon
    return result


def graph_interval(gene, token, graph):
    parts = [graph.get((gene, normalize_exon_id(t))) for t in str(token).split(':')]
    if not parts or any(p is None for p in parts):
        return None
    parts = sorted(set(parts), key=lambda p: p.start)
    if len({(p.chromosome, p.strand) for p in parts}) != 1:
        return None
    if any(a.end != b.start for a, b in zip(parts, parts[1:])):
        return None
    return ExonInterval(parts[0].chromosome, parts[0].start, parts[-1].end, parts[0].strand)


def junction_triplet(exons):
    """Return exact one-based closed introns in genomic order, on either strand."""
    if len(exons) != 3 or len(set(exons)) != 3:
        return None
    if len({(p.chromosome, p.strand) for p in exons}) != 1:
        return None
    left, middle, right = sorted(exons, key=lambda p: p.start)
    if not left.end < middle.start < middle.end < right.start:
        return None
    chrom, strand = middle.chromosome, middle.strand
    return ((chrom, left.end+1, middle.start, strand),
            (chrom, middle.end+1, right.start, strand),
            (chrom, left.end+1, right.start, strand))


def geo_samples(text):
    records = []
    for section in text.split('^SAMPLE = ')[1:]:
        lines = section.splitlines()
        row = {'geo_accession': lines[0].strip()}
        for line in lines[1:]:
            if line.startswith('!Sample_title = '):
                row['sample_id'] = line.split(' = ', 1)[1].strip()
            elif line.startswith('!Sample_characteristics_ch1 = '):
                key, _, value = line.split(' = ', 1)[1].partition(': ')
                row[key] = value.strip()
            elif line.startswith('!Sample_relation = SRA: '):
                match = re.search(r'(SRX\d+)', line)
                if match:
                    row['experiment_acc'] = match.group(1)
        title = row.get('sample_id', '')
        row['provider_title'] = title
        canonical = title.removeprefix('matched sample_')
        match = re.fullmatch(r'(L\d+)([TN])', canonical)
        if not match:
            raise ValueError('Undeclared biological sample alias: ' + title)
        row['sample_id'] = canonical
        row['patient_id'] = match.group(1)
        row['condition'] = 'tumor' if match.group(2) == 'T' else 'normal'
        records.append(row)
    data = pd.DataFrame(records)
    if data.sample_id.duplicated().any() or data.experiment_acc.duplicated().any():
        raise ValueError('Duplicate biological sample/experiment identity')
    cases = set(data.loc[data.condition.eq('tumor') & data.histology.eq('2'), 'patient_id'])
    data['included'] = data.patient_id.isin(cases)
    data['exclusion_reason'] = np.where(data.included, '', 'patient_without_confirmed_adenocarcinoma')
    return data


def aggregate_junction_runs(counts, columns, runs, samples):
    """Add technical partitions at count level, never count them as patients."""
    if counts.shape[1] != len(columns) or runs.rail_id.duplicated().any():
        raise ValueError('Junction column/run identity mismatch')
    lookup = samples.set_index('experiment_acc').sample_id
    run_lookup = runs.set_index('rail_id').experiment_acc
    assignments = run_lookup.reindex(columns).map(lookup)
    if assignments.isna().any():
        raise ValueError('Unmapped public junction column')
    result = np.zeros((counts.shape[0], len(samples)), dtype=np.int64)
    positions = {s: i for i, s in enumerate(samples.sample_id)}
    for j, sample in enumerate(assignments):
        result[:, positions[sample]] += counts[:, j]
    return result


def inclusion_proportion(left, right, skip, minimum_support=10):
    left, right, skip = (np.asarray(x, dtype=float) for x in (left, right, skip))
    if left.shape != right.shape or left.shape != skip.shape:
        raise ValueError('Junction count shapes differ')
    if any((~np.isfinite(x) | (x < 0)).any() for x in (left, right, skip)):
        raise ValueError('Counts must be finite and nonnegative')
    inclusion = (left+right)/2
    total = inclusion+skip
    psi = np.divide(inclusion, total, out=np.full_like(total, np.nan), where=total >= minimum_support)
    return psi, total


def frozen_scores(matrix, loadings, minimum_fraction=.8):
    """Project only frozen weights; retain unavailable modules and coverage."""
    matrix = matrix.reindex(loadings.event_key)
    present = matrix.notna().any(axis=1).to_numpy()
    weights = loadings.loading.to_numpy(float)
    coverage = dict(member_fraction=float(present.mean()),
                    loading_energy_fraction=float(np.sum(weights[present]**2)/np.sum(weights**2)))
    if min(coverage.values()) < minimum_fraction:
        return None, coverage
    values = matrix.to_numpy(float)
    center, scale = loadings.center.to_numpy(float), loadings.scale.to_numpy(float)
    observed = np.isfinite(values).sum(axis=0)/len(loadings)
    z = (np.where(np.isfinite(values), values, center[:, None])-center[:, None])/scale[:, None]
    scores = weights @ z
    scores[observed < minimum_fraction] = np.nan
    return pd.Series(scores, index=matrix.columns), coverage
