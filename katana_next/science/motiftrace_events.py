"""Local RNA consequences for all five rMATS event types, without quantification.

Coordinates are 0-based half-open; blocks always follow RNA 5'-3' order.
Regulatory pre-mRNA maps are separate from the two post-splicing local forms.
"""
from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

from katana_cosplice.motif_scanner import IndexedFasta
from katana_next.science.motiftrace_custom import fixed_pattern
from katana_next.science.se_rbp import _chromosome, mapped_blocks, sequence_from_blocks, tx_slice

EVENT_TYPES = ("SE", "A3SS", "A5SS", "MXE", "RI")


def _interval(row, left, right):
    values = np.asarray([row[left], row[right]], dtype=float)
    if not np.isfinite(values).all() or (values != np.floor(values)).any():
        raise ValueError(f"Noninteger/missing coordinates: {left}/{right}")
    a, b = map(int, values)
    if a < 0 or b <= a:
        raise ValueError(f"Empty/reversed interval: {left}/{right}")
    return a, b


def coalesce(blocks, strand):
    result = []
    for a, b in sorted(blocks):
        if result and a < result[-1][1]:
            raise ValueError("Overlapping blocks within a local RNA form")
        if result and a == result[-1][1]:
            result[-1] = (result[-1][0], b)
        else:
            result.append((a, b))
    return result if strand == "+" else result[::-1]


def subtract_blocks(first, second, strand):
    output = []
    for a, b in first:
        pieces = [(a, b)]
        for c, d in second:
            pieces = [(x, min(y, c)) for x, y in pieces if x < min(y, c)] + [
                (max(x, d), y) for x, y in pieces if max(x, d) < y]
        output.extend(pieces)
    return coalesce(output, strand)


def event_forms(row, anchor_nt=100):
    """Respect rMATS PSI semantics, including the negative-strand MXE reversal."""
    kind, strand = str(row["event_type"]), str(row["strand"])
    if kind not in EVENT_TYPES or strand not in {"+", "-"}:
        raise ValueError("Explicit supported event type and +/- strand required")
    if not isinstance(anchor_nt, int) or anchor_nt < 1:
        raise ValueError("Positive integer external exon anchor required")
    get = lambda a, b: _interval(row, a, b)
    if kind in {"SE", "MXE", "RI"}:
        up, down = get("upstreamES", "upstreamEE"), get("downstreamES", "downstreamEE")
        if up[1] >= down[0]:
            raise ValueError("Flanking exons must have a nonempty intervening interval")
        u, d = (max(up[0], up[1]-anchor_nt), up[1]), (down[0], min(down[1], down[0]+anchor_nt))
        if kind == "SE":
            target = get("exonStart_0base", "exonEnd")
            if not up[1] < target[0] < target[1] < down[0]:
                raise ValueError("SE target must lie strictly between its flanks")
            first, second = [u, target, d], [u, d]
        elif kind == "MXE":
            a, b = get("1stExonStart_0base", "1stExonEnd"), get("2ndExonStart_0base", "2ndExonEnd")
            if not up[1] < a[0] < a[1] <= b[0] < b[1] < down[0]:
                raise ValueError("MXE exon order/overlap is invalid")
            inc, exc = (a, b) if strand == "+" else (b, a)
            first, second = [u, inc, d], [u, exc, d]
        else:
            retained = get("riExonStart_0base", "riExonEnd")
            if retained != (up[0], down[1]):
                raise ValueError("RI retained-exon extent disagrees with the flanking exons")
            first, second = [u, (up[1], down[0]), d], [u, d]
    else:
        long = get("longExonStart_0base", "longExonEnd")
        short = get("shortES", "shortEE")
        flank = get("flankingES", "flankingEE")
        suffix = (kind == "A5SS") == (strand == "+")
        if suffix:
            if not long[0] == short[0] < short[1] < long[1] < flank[0]:
                raise ValueError("Alternative site geometry conflicts with type/strand")
            flank = flank[0], min(flank[1], flank[0]+anchor_nt)
        else:
            if not flank[1] < long[0] < short[0] < short[1] == long[1]:
                raise ValueError("Alternative site geometry conflicts with type/strand")
            flank = max(flank[0], flank[1]-anchor_nt), flank[1]
        first, second = [long, flank], [short, flank]
    first, second = coalesce(first, strand), coalesce(second, strand)
    return {"inclusion": first, "exclusion": second,
            "affected_inclusion": subtract_blocks(first, second, strand),
            "affected_exclusion": subtract_blocks(second, first, strand)}


def extract_event_rna(events, fasta_path, *, anchor_nt=100, intron_nt=200,
                      max_total_bases=None, feature_only=False, check_cancelled=lambda: None):
    """Reconstruct every requested local RNA; the legacy base-budget argument is ignored."""
    if events.empty or events.event_id.duplicated().any() or events.comparison_id.nunique() != 1:
        raise ValueError("One nonempty comparison with unique event IDs is required")
    fasta, rows = IndexedFasta(Path(fasta_path)), []
    with fasta.fasta.open("rb") as stream:
        for row in events.sort_values("event_id").to_dict("records"):
            check_cancelled()
            forms = event_forms(row, anchor_nt)
            chrom, strand = _chromosome(fasta, str(row["chr"])), str(row["strand"])
            specs = [(name, blocks, "mature_local_form" if name in {"inclusion", "exclusion"}
                      else "affected_sequence") for name, blocks in forms.items()]
            if feature_only:
                specs = [s for s in specs if s[2] == "affected_sequence"]
            for form in (() if feature_only else ("inclusion", "exclusion")):
                for i, (a, b) in enumerate(zip(forms[form], forms[form][1:]), 1):
                    alen, blen = a[1]-a[0], b[1]-b[0]
                    specs.append((f"{form}_junction_{i}", [tx_slice(a, strand, max(0, alen-anchor_nt), alen),
                        tx_slice(b, strand, 0, min(anchor_nt, blen))], "splice_junction"))
                    gap = (a[1], b[0]) if strand == "+" else (b[1], a[0])
                    size = gap[1]-gap[0]
                    # Keep the same declared canonical-site masks as the existing SE module.
                    if size > 6:
                        specs.append((f"{form}_junction_{i}_donor_intron", [tx_slice(gap, strand, 6, min(size, intron_nt))], "regulatory_pre_mRNA"))
                    if size > 20:
                        specs.append((f"{form}_junction_{i}_acceptor_intron", [tx_slice(gap, strand, max(0, size-intron_nt), size-20)], "regulatory_pre_mRNA"))
            for region, blocks, layer in specs:
                if not blocks:
                    continue
                sequence = sequence_from_blocks(fasta, stream, chrom, strand, blocks)
                joins = np.cumsum([b-a for a, b in blocks])[:-1].tolist()
                rows.append({"event_id": str(row["event_id"]), "gene_id": str(row["gene_id"]),
                    "gene_symbol": str(row.get("gene_symbol", "")), "comparison_id": str(row["comparison_id"]),
                    "event_type": row["event_type"], "chr": chrom, "strand": strand,
                    "role": row.get("role", "foreground"), "region": region, "layer": layer,
                    "blocks_json": json.dumps(blocks), "junctions_json": json.dumps(joins),
                    "sequence": sequence, "sequence_RNA": sequence.replace("T", "U"),
                    "sequence_sha256": hashlib.sha256(sequence.encode()).hexdigest(),
                    "length": len(sequence), "dpsi": row.get("dpsi", np.nan), "DAS_FDR": row.get("fdr", np.nan)})
    return pd.DataFrame(rows)


def validate_models(models):
    if not models:
        raise ValueError("Require at least one explicitly selected model; no silent library truncation")
    for model in models.values():
        if model.get("kind") == "consensus":
            _, width = fixed_pattern(model["pattern"])
        elif model.get("kind") == "PWM":
            matrix = np.asarray(model["matrix"], dtype=float)
            if matrix.ndim != 2 or matrix.shape[1] != 4 or not 1 <= len(matrix) <= 25 or not np.isfinite(matrix).all() or (matrix < 0).any() or not np.allclose(matrix.sum(axis=1), 1):
                raise ValueError("Expected native finite probability PWM, width 1-25")
            width = len(matrix)
        else:
            raise ValueError("Unknown motif representation")
        if model.get("width", width) != width:
            raise ValueError("Motif width disagrees with its representation")


def score_positions(sequence, model, background=(.25, .25, .25, .25)):
    """All valid starts; short/ambiguous opportunities are NA, never zero-filled."""
    sequence = sequence.upper().replace("U", "T")
    if model["kind"] == "consensus":
        pattern, width = fixed_pattern(model["pattern"])
    else:
        width = len(model["matrix"])
    n = max(0, len(sequence)-width+1)
    valid = np.array([set(sequence[a:a+width]) <= set("ACGT") for a in range(n)], bool)
    scores = np.full(n, np.nan)
    if model["kind"] == "consensus":
        scores[valid] = 0
        for match in re.finditer("(?=(" + pattern + "))", sequence):
            scores[match.start(1)] = 1
        return scores, scores == 1, width
    bg = np.asarray(background, float)
    if bg.shape != (4,) or (bg <= 0).any() or not np.isclose(bg.sum(), 1):
        raise ValueError("A positive fixed four-base background is required")
    matrix = (20*np.asarray(model["matrix"], float)+.1*bg)/20.1
    odds = matrix/bg
    relative = matrix/matrix.max(axis=1, keepdims=True)
    candidate = np.zeros(n, bool)
    starts = np.flatnonzero(valid)
    values, rel = np.ones(len(starts)), np.ones(len(starts))
    lookup = np.full(256, 4, dtype=np.uint8)
    lookup[np.frombuffer(b"ACGT", dtype=np.uint8)] = np.arange(4)
    encoded = lookup[np.frombuffer(sequence.encode("ascii"), dtype=np.uint8)]
    for i in range(width):
        values *= odds[i, encoded[starts+i]]
        rel *= relative[i, encoded[starts+i]]
    scores[starts] = values
    candidate[starts] = rel >= .8
    return scores, candidate, width


class ScoringCache:
    """Per-analysis bounded memoization; eviction never removes scientific tests."""
    def __init__(self, max_bytes=32 * 1024 * 1024):
        from collections import OrderedDict
        self.values = OrderedDict()
        self.bytes = 0
        self.max_bytes = max_bytes

    def __call__(self, sequence, model):
        key = (sequence, json.dumps(model, sort_keys=True, separators=(',', ':')))
        if key in self.values:
            self.values.move_to_end(key)
            return self.values[key][0]
        value = score_positions(sequence, model)
        size = value[0].nbytes + value[1].nbytes + len(sequence) + len(key[1])
        if size <= self.max_bytes:
            while self.values and self.bytes + size > self.max_bytes:
                _, (_, removed) = self.values.popitem(last=False)
                self.bytes -= removed
            self.values[key] = (value, size)
            self.bytes += size
        return value


def trace_instances(rna, events, models, *, check_cancelled=lambda: None, scorer=None):
    """Identify loci, not net hit-count changes; PWM candidates are not P-values."""
    validate_models(models)
    records, changes, weighted = [], [], []
    lookup = events.set_index("event_id")
    for event_id, group in rna[rna.layer.eq("mature_local_form")].groupby("event_id", sort=False):
        check_cancelled()
        source = lookup.loc[event_id]
        for token, model in models.items():
            loci = {}
            for form in group.to_dict("records"):
                scores, candidate, width = (scorer or score_positions)(form["sequence"], model)
                for start in np.flatnonzero(candidate):
                    start = int(start)
                    blocks = coalesce(mapped_blocks(json.loads(form["blocks_json"]), form["strand"], start, start+width), form["strand"])
                    key = json.dumps(blocks)
                    hit = {"event_id": event_id, "event_type": form["event_type"], "gene_id": form["gene_id"],
                        "comparison_id": form["comparison_id"], "model_id": token, "RBP": model.get("RBPs", ""),
                        "form": form["region"], "RNA_start0": start, "RNA_end0": start+width,
                        "hit_RNA": form["sequence_RNA"][start:start+width], "chr": form["chr"], "strand": form["strand"],
                        "genomic_blocks_json": key, "score": float(scores[start]),
                        "crosses_splice_join": any(start < j < start+width for j in json.loads(form["junctions_json"])),
                        "evidence": "finite_pattern_match" if model["kind"] == "consensus" else "relative_likelihood_candidate_ge_0.8_not_site_P"}
                    records.append(hit)
                    loci.setdefault(key, {})[form["region"]] = hit
            psi = psi_values(source)
            for key, forms in loci.items():
                inc, exc = "inclusion" in forms, "exclusion" in forms
                state = "retained" if inc and exc else "inclusion_only" if inc else "exclusion_only"
                first = next(iter(forms.values()))
                delta = float(source.get("dpsi", np.nan)) * (int(inc)-int(exc))
                changes.append({k: first[k] for k in ("event_id", "event_type", "gene_id", "comparison_id", "model_id", "RBP", "chr", "strand", "genomic_blocks_json", "hit_RNA", "evidence")} |
                    {"instance_state": state, "delta_sequence_opportunity": delta,
                     "consequence_direction": "unchanged" if inc == exc else "increased" if delta > 0 else "decreased" if delta < 0 else "unknown_or_zero"})
                for value in psi:
                    weighted.append({"event_id": event_id, "model_id": token, "genomic_blocks_json": key,
                        "instance_state": state, **value,
                        "sequence_opportunity": int(inc)*value["PSI"]+int(exc)*(1-value["PSI"]),
                        "meaning": "local-form sequence availability, not expression-adjusted occupancy"})
    return pd.DataFrame(records), pd.DataFrame(changes), pd.DataFrame(weighted)


def psi_values(row):
    output = []
    reverse = row.get("reverse_direction", False)
    if pd.isna(reverse):
        raise ValueError("Unknown PSI comparison direction")
    if isinstance(reverse, str):
        if reverse.lower() not in {"true", "false"}:
            raise ValueError("Invalid PSI direction flag")
        reverse = reverse.lower() == "true"
    for group, field in (("numerator", "IncLevel2" if reverse else "IncLevel1"),
                          ("denominator", "IncLevel1" if reverse else "IncLevel2")):
        raw = row.get(field, "")
        for i, token in enumerate(str(raw).split(","), 1):
            if token.strip().lower() in {"", "na", "nan", "none"}:
                value = np.nan
            else:
                value = float(token)
                if not np.isfinite(value) or not 0 <= value <= 1:
                    raise ValueError("PSI outside [0,1]")
            output.append({"group": group, "replicate_index_not_sample_identity": i, "PSI": value})
    return output


def paired_form_scores(rna, events, models, *, check_cancelled=lambda: None, scorer=None):
    """Descriptive opportunity-normalized form scores, including zero-hit events."""
    rows = []
    event_lookup = events.set_index("event_id")
    for event_id, group in rna[rna.layer.eq("mature_local_form")].groupby("event_id", sort=False):
        check_cancelled()
        source = event_lookup.loc[event_id]
        for token, model in models.items():
            record = {"event_id": event_id, "event_type": source.event_type, "gene_id": source.gene_id,
                "model_id": token, "RBP": model.get("RBPs", ""), "dpsi": source.dpsi,
                "score_definition": "mean odds over valid starts" if model["kind"] == "PWM" else "overlapping hits per valid motif start",
                "interpretation": "descriptive local-form sequence score, not binding or a per-event statistical test"}
            for form in group.to_dict("records"):
                positions, _, _ = (scorer or score_positions)(form["sequence"], model)
                valid = np.isfinite(positions)
                record[form["region"]+"_valid_starts"] = int(valid.sum())
                record[form["region"]+"_normalized_score"] = float(np.mean(positions[valid])) if valid.any() else np.nan
            first, second = record["inclusion_normalized_score"], record["exclusion_normalized_score"]
            record["form_score_difference"] = first-second
            record["delta_PSI_weighted_normalized_score"] = source.dpsi*(first-second)
            rows.append(record)
    return pd.DataFrame(rows)


def match_event_controls(events, selected_ids, rna):
    """Motif-blind, one-to-one matching of unique genes within an AS type."""
    selected_ids = set(map(str, selected_ids))
    if not selected_ids or not selected_ids <= set(events.event_id):
        raise ValueError("Exact nonempty selected event identity required")
    features = {}
    for token, group in rna[rna.layer.eq("affected_sequence")].groupby("event_id"):
        seq = "".join(group.sequence)
        valid = sum(seq.count(base) for base in "ACGT")
        features[token] = (len(seq), (seq.count("G")+seq.count("C"))/valid if valid else np.nan,
                           1-valid/len(seq) if seq else 1)
    significant = events.stable_significant.fillna(False) | events.fdr.le(.05) | events.event_id.isin(selected_ids)
    forbidden = set(events.loc[significant, "gene_id"])
    pool = events[events.measurement_eligible.fillna(False) & events.fdr.gt(.1) & events.dpsi.abs().lt(.05)
                  & ~events.gene_id.isin(forbidden)].copy()
    def info(row):
        length, gc, ambiguous = features.get(row["event_id"], (0, np.nan, 1))
        vals = [p["PSI"] for p in psi_values(row) if p["group"] == "denominator" and np.isfinite(p["PSI"])]
        return length, gc, ambiguous, float(np.mean(vals)) if vals else np.nan, min(float(row["numerator_mean_reads"]), float(row["denominator_mean_reads"]))
    backgrounds = {row["event_id"]: (row, info(row)) for row in pool.to_dict("records")}
    fg_genes, bg_genes, result = set(), set(), []
    for row in events[events.event_id.isin(selected_ids)].sort_values("event_id").to_dict("records"):
        gene = row["gene_id"]
        record = {"foreground_event": row["event_id"], "foreground_gene": gene, "event_type": row["event_type"],
            "direction": "inclusion_up" if row["dpsi"] > 0 else "inclusion_down", "background_event": "", "background_gene": "", "status": "unmatched"}
        if pd.isna(gene) or not str(gene).strip():
            record["status"] = "missing_gene_id"
        elif gene in fg_genes:
            record["status"] = "same_gene_extra_event_trace_only"
        else:
            fg_genes.add(gene)
            length, gc, ambiguous, psi, depth = info(row)
            options = []
            if ambiguous <= .05 and length > 0 and np.isfinite([gc, psi, depth]).all() and depth > 0:
                for token, (other, feature) in backgrounds.items():
                    n, g, amb, p, d = feature
                    if other["event_type"] != row["event_type"] or other["gene_id"] in bg_genes or pd.isna(other["gene_id"]) or not str(other["gene_id"]).strip() or d <= 0:
                        continue
                    distance = (abs(gc-g), abs(length-n)/length, abs(psi-p), abs(np.log2(depth/d)))
                    if amb <= .05 and np.isfinite(distance).all() and all(a <= b for a, b in zip(distance, (.05, .2, .15, 1))):
                        options.append((distance, token, other))
            if options:
                _, token, other = min(options, key=lambda value: (value[0], value[1]))
                bg_genes.add(other["gene_id"])
                record.update(background_event=token, background_gene=other["gene_id"], status="matched")
        result.append(record)
    return pd.DataFrame(result)
