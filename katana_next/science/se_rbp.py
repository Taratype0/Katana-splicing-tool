"""SE sequence consequences and matched, direction-specific RBP hypotheses.

Only postprocess existing DAS results. Genomic blocks are half-open and listed
in RNA 5'-3' order; junction sites must not be exported as continuous DNA hits.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import binomtest

from katana_cosplice.motif_scanner import IndexedFasta, reverse_complement
from katana_next.science.motif_significance import calibrated_scan
from katana_next.science.statistics import bh_adjust

REGIONS = ("cassette_exon", "inclusion_junction_5prime", "inclusion_junction_3prime",
           "exclusion_junction", "upstream_intron", "downstream_intron")
POLICY = {
    "foreground": "existing stable_significant SE; no differential model refit",
    "control": "measurement_eligible SE, original FDR > 0.1 and abs(dPSI) < 0.05",
    "independence": "one lexicographic event per foreground gene across both directions; unique control genes",
    "matching": "cassette GC <=0.05, length <=20%, minimum group mean read depth within 2-fold, denominator mean PSI <=0.15",
    "site_pvalue": 1e-4, "intron_bp": 200,
    "splice_consensus_mask": "upstream intron last 20 nt and downstream intron first 6 nt excluded",
    "hypotheses": "RBP x evidence class x 6 regions x inclusion direction; all tested rows in one BH family",
    "test": "one-sided exact paired binomial/McNemar; motifs of same RBP/evidence unioned",
    "claim_boundary": "predicted motif availability/association, not binding, affinity change, or causal function",
}


def _blocks(row):
    names = (("upstreamES", "upstreamEE"), ("exonStart_0base", "exonEnd"),
             ("downstreamES", "downstreamEE"))
    blocks = []
    for left, right in names:
        a, b = float(getattr(row, left)), float(getattr(row, right))
        if not np.isfinite([a, b]).all() or not a.is_integer() or not b.is_integer() or a < 0 or b <= a:
            raise ValueError(f"Invalid SE exon coordinates: {row.event_id}")
        blocks.append((int(a), int(b)))
    if not (blocks[0][1] <= blocks[1][0] and blocks[1][1] <= blocks[2][0]):
        raise ValueError(f"Overlapping or misordered SE exons: {row.event_id}")
    if row.strand not in {"+", "-"}:
        raise ValueError("Explicit SE strand required")
    return blocks if row.strand == "+" else blocks[::-1]


def tx_slice(block, strand, start, stop):
    """Slice a genomic block using transcript-oriented relative offsets."""
    a, b = block
    if not 0 <= start <= stop <= b-a:
        raise ValueError("RNA slice outside genomic block")
    return (a+start, a+stop) if strand == "+" else (b-stop, b-start)


def mapped_blocks(blocks, strand, start, stop):
    """Project an RNA substring onto one or more genomic blocks, losslessly."""
    result, offset = [], 0
    for a, b in blocks:
        left, right = max(start, offset), min(stop, offset+b-a)
        if right > left:
            result.append(tx_slice((a, b), strand, left-offset, right-offset))
        offset += b-a
    if start < 0 or stop > offset or sum(b-a for a, b in result) != stop-start:
        raise ValueError("Motif coordinate projection is incomplete")
    return result


def _chromosome(fasta, chrom):
    if chrom in fasta.index:
        return chrom
    other = chrom.removeprefix("chr") if chrom.startswith("chr") else "chr"+chrom
    if other not in fasta.index:
        raise ValueError(f"Chromosome absent from FASTA: {chrom}")
    return other


def sequence_from_blocks(fasta, stream, chrom, strand, blocks):
    parts = []
    for a, b in blocks:
        seq = fasta.fetch(chrom, a, b, stream=stream)
        if len(seq) != b-a:
            raise ValueError("FASTA interval was truncated")
        parts.append(seq if strand == "+" else reverse_complement(seq))
    return "".join(parts)


def _mean_psi(value):
    values = pd.to_numeric(pd.Series(str(value).split(",")), errors="coerce").dropna()
    return float(values.mean()) if len(values) else np.nan


def prepare_se_background(frame, selected_ids, fasta_path):
    """Freeze motif-blind matching, auditing every selected event and loss."""
    if frame.comparison_id.nunique() != 1:
        raise ValueError("SE RBP analysis requires exactly one comparison")
    required = {"measurement_eligible", "stable_significant", "fdr", "dpsi", "IncLevel1", "IncLevel2",
                "numerator_mean_reads", "denominator_mean_reads", "gene_id", "event_id", "event_type"}
    if not required <= set(frame):
        raise ValueError(f"Missing frozen DAS/background columns: {sorted(required-set(frame))}")
    if frame.event_id.duplicated().any():
        raise ValueError("Ambiguous duplicate event identity in comparison")
    selected = frame[frame.event_id.isin(selected_ids)].copy()
    if set(selected.event_id) != set(selected_ids) or not selected.event_type.eq("SE").all():
        raise ValueError("Select exact SE events from this complete comparison")
    if not selected.stable_significant.fillna(False).all() or not selected.measurement_eligible.fillna(False).all():
        raise ValueError("Foreground must retain existing read-supported DAS qualification")
    if not selected.dpsi.ne(0).all():
        raise ValueError("Selected SE has no inclusion direction")
    # Exclude ALL foreground/significant genes, including opposite-direction DAS.
    excluded_genes = set(frame.loc[frame.stable_significant.fillna(False) | frame.fdr.le(.05), "gene_id"])
    candidates = frame[frame.event_type.eq("SE") & frame.measurement_eligible.fillna(False)
                       & frame.fdr.gt(.1) & frame.dpsi.abs().lt(.05) & ~frame.gene_id.isin(excluded_genes)]
    pool = pd.concat([selected, candidates]).drop_duplicates("event_id").sort_values("event_id")
    fasta, features, audit = IndexedFasta(Path(fasta_path)), [], []
    with fasta.fasta.open("rb") as stream:
        for row in pool.itertuples(index=False):
            try:
                blocks = _blocks(row)
                chrom = _chromosome(fasta, str(row.chr))
                seq = sequence_from_blocks(fasta, stream, chrom, row.strand, [blocks[1]])
                valid = sum(seq.count(x) for x in "ACGT")
                reverse = bool(getattr(row, "reverse_direction", False))
                features.append({"event_id": row.event_id, "gene_id": row.gene_id, "bases": len(seq),
                    "gc": (seq.count("G")+seq.count("C"))/max(valid, 1), "ambiguous_fraction": 1-valid/len(seq),
                    "baseline_psi": _mean_psi(row.IncLevel1 if reverse else row.IncLevel2),
                    "depth": min(row.numerator_mean_reads, row.denominator_mean_reads),
                    "direction": "inclusion_up" if row.dpsi > 0 else "inclusion_down"})
                audit.append({"event_id": row.event_id, "status": "sequence_valid"})
            except (ValueError, KeyError) as error:
                audit.append({"event_id": row.event_id, "status": "sequence_rejected", "reason": str(error)})
    features = pd.DataFrame(features)
    fg = selected.sort_values("event_id")
    lookup = features.set_index("event_id")
    control = features[features.event_id.isin(candidates.event_id) & features.ambiguous_fraction.le(.05)
                       & features.baseline_psi.notna() & features.depth.gt(0)]
    used_fg, used_bg, pairs = set(), set(), []
    for row in fg.itertuples(index=False):
        pair = {"foreground_event": row.event_id, "foreground_gene": row.gene_id,
                "direction": "inclusion_up" if row.dpsi > 0 else "inclusion_down",
                "background_event": "", "background_gene": "", "status": "unmatched"}
        if pd.isna(row.gene_id) or not str(row.gene_id):
            pair["status"] = "missing_gene_id"
        elif row.gene_id in used_fg:
            pair["status"] = "same_gene_additional_event_sites_only"
        else:
            used_fg.add(row.gene_id)
            if row.event_id not in lookup.index:
                pair["status"] = "invalid_sequence"
            else:
                f = lookup.loc[row.event_id]
                c = control[~control.gene_id.isin(used_bg)].copy()
                c["gc_distance"] = (c.gc-f.gc).abs()
                c["length_distance"] = (c.bases-f.bases).abs()/f.bases
                c["psi_distance"] = (c.baseline_psi-f.baseline_psi).abs()
                c["depth_log2_distance"] = np.abs(np.log2(c.depth/f.depth)) if f.depth > 0 else np.inf
                c = c[c.gc_distance.le(.05) & c.length_distance.le(.20) & c.psi_distance.le(.15)
                      & c.depth_log2_distance.le(1)]
                if f.ambiguous_fraction > .05:
                    pair["status"] = "ambiguous_foreground"
                elif len(c):
                    b = c.sort_values(["gc_distance", "length_distance", "psi_distance", "depth_log2_distance", "event_id"]).iloc[0]
                    used_bg.add(b.gene_id)
                    pair.update(background_event=b.event_id, background_gene=b.gene_id, status="matched")
                    for key in ("gc", "bases", "baseline_psi", "depth"):
                        pair["foreground_"+key], pair["background_"+key] = f[key], b[key]
        pairs.append(pair)
    return pd.DataFrame(pairs), features, pd.DataFrame(audit)


def extract_se_sequences(events, fasta_path, max_motif_width, intron_bp=200):
    fasta, rows = IndexedFasta(Path(fasta_path)), []
    flank = max_motif_width-1
    if flank < 1:
        raise ValueError("Motif widths must be at least two")
    with fasta.fasta.open("rb") as stream:
        for row in events.sort_values("event_id").itertuples(index=False):
            up, exon, down = _blocks(row)
            chrom, strand = _chromosome(fasta, str(row.chr)), row.strand
            tail = lambda b: tx_slice(b, strand, max(0, b[1]-b[0]-flank), b[1]-b[0])
            head = lambda b: tx_slice(b, strand, 0, min(flank, b[1]-b[0]))
            specifications = [("cassette_exon", [exon], None, "inclusion")]
            for region, a, b, form in (("inclusion_junction_5prime", tail(up), head(exon), "inclusion"),
                                       ("inclusion_junction_3prime", tail(exon), head(down), "inclusion"),
                                       ("exclusion_junction", tail(up), head(down), "exclusion")):
                specifications.append((region, [a, b], a[1]-a[0], form))
            upstream = (up[1], exon[0]) if strand == "+" else (exon[1], up[0])
            downstream = (exon[1], down[0]) if strand == "+" else (down[1], exon[0])
            ulen, dlen = upstream[1]-upstream[0], downstream[1]-downstream[0]
            if ulen > 20:
                specifications.append(("upstream_intron", [tx_slice(upstream, strand, max(0, ulen-intron_bp), ulen-20)], None, "pre_mRNA"))
            if dlen > 6:
                specifications.append(("downstream_intron", [tx_slice(downstream, strand, 6, min(dlen, intron_bp))], None, "pre_mRNA"))
            for region, blocks, junction, form in specifications:
                seq = sequence_from_blocks(fasta, stream, chrom, strand, blocks)
                if not seq:
                    continue
                rows.append({"window_id": f"s{len(rows):07d}", "event_id": row.event_id,
                    "gene_id": row.gene_id, "gene_symbol": getattr(row, "gene_symbol", ""),
                    "chr": chrom, "genomic_strand": strand, "region": region, "form": form,
                    "blocks_json": json.dumps(blocks), "junction_offset": junction,
                    "anchor_start": -min(ulen, intron_bp) if region == "upstream_intron" else 6 if region == "downstream_intron" else 0,
                    "dpsi": row.dpsi, "DAS_fdr": row.fdr, "sequence": seq,
                    "sequence_sha256": hashlib.sha256(seq.encode()).hexdigest(),
                    "sequence_rna": seq.replace("T", "U"), "role": getattr(row, "role", "foreground"),
                    "direction": getattr(row, "analysis_direction", "inclusion_up" if row.dpsi > 0 else "inclusion_down"),
                    # FIMO scans RNA-oriented strings. Coordinates are projected below.
                    "strand": "+", "start_0based": 0, "end_exclusive": len(seq)})
    return pd.DataFrame(rows)


def scan_se_sequences(motifs, windows, output, runtime, threshold=1e-4):
    """Identical PWMs are computed once, retaining every RBP/source assignment."""
    unique, keys, aliases = [], {}, []
    for motif in motifs:
        pwm = np.asarray(motif.pwm, float)
        pwm = pwm/pwm.sum(axis=1, keepdims=True)
        digest = hashlib.sha256(pwm.astype("<f8").tobytes()).hexdigest()
        if digest not in keys:
            keys[digest] = len(unique)
            unique.append(motif)
        aliases.append({"scan_token": f"m{keys[digest]:06d}", "motif_token": f"a{len(aliases):06d}",
                        "motif_id": motif.motif_id, "RBP": motif.rbp, "source": motif.source,
                        "evidence_class": motif.evidence_class, "target_species": motif.target_species,
                        "pwm_sha256": digest, "width": len(pwm)})
    catalog = pd.DataFrame(aliases)
    raw, _, provenance = calibrated_scan(unique, windows.drop(columns="sequence_rna"), Path(output), runtime, threshold)
    raw = raw.drop(columns=["motif_id", "RBP", "source", "evidence_class", "target_species"], errors="ignore")
    raw = raw.rename(columns={"motif_token": "scan_token"})
    # Only genuinely crossing hits are junction hypotheses; flank-only hits do
    # not become splicing-created sites just because they are near a junction.
    crossed = (raw.start_1based-1 < raw.junction_offset) & (raw.stop_1based > raw.junction_offset)
    raw = raw[raw.junction_offset.isna() | crossed].copy()
    by_window = windows.set_index("window_id")
    sequences, contexts, mappings, positions = [], [], [], []
    for hit in raw.itertuples(index=False):
        w = by_window.loc[hit.window_id]
        a, b = int(hit.start_1based)-1, int(hit.stop_1based)
        sequences.append(w.sequence[a:b].replace("T", "U"))
        contexts.append(w.sequence[max(0, a-15):a].lower().replace("t", "u") + "[" + w.sequence[a:b].replace("T", "U") + "]" + w.sequence[b:b+15].lower().replace("t", "u"))
        mappings.append(json.dumps(mapped_blocks(json.loads(w.blocks_json), w.genomic_strand, a, b)))
        positions.append((a+b)/2/len(w.sequence)*100 if w.region == "cassette_exon" else w.anchor_start+(a+b)/2 if w.form == "pre_mRNA" else (a+b)/2-w.junction_offset)
    raw["matched_rna_sequence"], raw["rna_context"], raw["hit_genomic_blocks_json"], raw["relative_position"] = sequences, contexts, mappings, positions
    raw = raw.drop(columns=["genomic_start_0based", "genomic_end_exclusive", "strand"], errors="ignore")
    raw["strand"] = raw.genomic_strand
    sites = raw.merge(catalog, on="scan_token", validate="many_to_many")
    sites["availability_interpretation"] = np.where(sites.form.eq("pre_mRNA"), "regulatory_association_only",
        np.where(sites.form.eq("inclusion") == sites.direction.eq("inclusion_up"),
                 "increased_representation_in_numerator", "decreased_representation_in_numerator"))
    sites.loc[sites.role.eq("background"), "availability_interpretation"] = "matched_background_not_DAS"
    candidates = catalog.groupby("pwm_sha256").RBP.agg(lambda x: ";".join(sorted(set(x))))
    sites["indistinguishable_PWM_RBPs"] = sites.pwm_sha256.map(candidates)
    provenance.update(motif_assignments=len(motifs), unique_pwms_scanned=len(unique),
                      exact_pwm_reuse="byte-identical normalized float64 matrices only; all original assignments retained",
                      junction_rule="start < boundary < stop; genomic blocks in RNA order",
                      site_q_scope="all unique PWM positions, before junction filter; not enrichment FDR")
    return sites, catalog, provenance


def region_enrichment(sites, catalog, pairs, windows):
    rows = []
    groups = catalog.assign(rbp_key=catalog.RBP.str.casefold()).groupby(["rbp_key", "evidence_class"], sort=True)
    for (rbp, evidence), motifs in groups:
        rbp_name = sorted(set(motifs.RBP), key=lambda x: (x.isupper(), x))[0]
        own = sites[sites.motif_token.isin(motifs.motif_token)]
        for region in REGIONS:
            valid = set(windows.loc[windows.region.eq(region), "event_id"])
            hit_ids = set(own.loc[own.region.eq(region), "event_id"])
            for direction in ("inclusion_up", "inclusion_down"):
                pairs0 = pairs[pairs.status.eq("matched") & pairs.direction.eq(direction)]
                paired = pairs0[pairs0.foreground_event.isin(valid) & pairs0.background_event.isin(valid)]
                a, b = paired.foreground_event.isin(hit_ids), paired.background_event.isin(hit_ids)
                only_a, only_b = int((a & ~b).sum()), int((~a & b).sum())
                n = len(paired)
                p = float(binomtest(only_a, only_a+only_b, alternative="greater").pvalue) if only_a+only_b else 1.
                difference = a.to_numpy(dtype=float)-b.to_numpy(dtype=float)
                effect = float(difference.mean()) if n else np.nan
                se = float(difference.std(ddof=1)/np.sqrt(n)) if n > 1 else np.nan
                rows.append({"RBP": rbp_name, "evidence_class": evidence, "region": region, "direction": direction,
                    "motif_ids": ";".join(motifs.motif_id), "motif_tokens": ";".join(motifs.motif_token),
                    "matched_gene_pairs": n, "region_ineligible_pairs": len(pairs0)-n,
                    "foreground_hit_genes": int(a.sum()), "background_hit_genes": int(b.sum()),
                    "foreground_only_pairs": only_a, "background_only_pairs": only_b,
                    "hit_fraction_difference": effect, "difference_CI95_low_approx": max(-1, effect-1.96*se),
                    "difference_CI95_high_approx": min(1, effect+1.96*se),
                    "enrichment_pvalue": p if n >= 2 else np.nan,
                    "status": "tested" if n >= 2 else "insufficient_matched_genes",
                    "foreground_hit_events": ";".join(paired.loc[a, "foreground_event"]),
                    "background_hit_events": ";".join(paired.loc[b, "background_event"])})
    result = pd.DataFrame(rows)
    result["enrichment_fdr"] = bh_adjust(result.enrichment_pvalue)
    result["fdr_family"] = "all RBP/evidence/region/direction hypotheses in this comparison"
    return result


def run_se_rbp(frame, selected_ids, motifs, fasta, runtime, directory, *, threshold=1e-4, progress=None):
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    report = progress or (lambda message: None)
    report("Matching SE backgrounds without inspecting motif hits")
    pairs, features, audit = prepare_se_background(frame, selected_ids, fasta)
    matched = pairs[pairs.status.eq("matched")]
    if matched.empty:
        raise ValueError("No eligible independent matched SE gene pairs")
    # Keep ALL valid selected SE sequences, even non-representative/unmatched.
    valid = set(audit.loc[audit.status.eq("sequence_valid"), "event_id"])
    ids = (set(selected_ids) & valid) | set(matched.background_event)
    events = frame[frame.event_id.isin(ids)].copy()
    events["role"] = np.where(events.event_id.isin(selected_ids), "foreground", "background")
    bg_direction = dict(zip(matched.background_event, matched.direction))
    events["analysis_direction"] = events.event_id.map(bg_direction).fillna(events.dpsi.map(lambda x: "inclusion_up" if x > 0 else "inclusion_down"))
    width = max(len(m.pwm) for m in motifs)
    windows = extract_se_sequences(events, fasta, width)
    for name, table in (("matching", pairs), ("sequence_features", features), ("sequence_audit", audit), ("sequences", windows)):
        table.to_parquet(directory/(name+".parquet"), index=False)
        table.to_csv(directory/(name+".tsv"), sep="\t", index=False)
    report(f"FIMO scanning {len(windows)} RNA-oriented sequences; one local process")
    sites, catalog, provenance = scan_se_sequences(motifs, windows, directory/"fimo", runtime, threshold)
    enrichment = region_enrichment(sites, catalog, pairs, windows)
    for name, table in (("sites", sites), ("motif_catalog", catalog), ("enrichment", enrichment)):
        table.to_parquet(directory/(name+".parquet"), index=False)
        table.to_csv(directory/(name+".tsv"), sep="\t", index=False)
    with (directory/"all_sequences.rna.fa").open("w", encoding="ascii") as out:
        for row in windows.itertuples(index=False):
            out.write(f">{row.window_id} {row.region} {row.gene_symbol} {row.role} genomic_strand={row.genomic_strand}\n{row.sequence_rna}\n")
    foreground = sites[sites.role.eq("foreground")]
    foreground.to_csv(directory/"DAS_actual_motif_hits.tsv", sep="\t", index=False)
    with (directory/"DAS_actual_motif_hits.rna.fa").open("w", encoding="ascii") as out:
        for i, row in enumerate(foreground.itertuples(index=False)):
            out.write(f">hit{i:07d} {row.RBP} {row.motif_token} {row.window_id}:{row.start_1based}-{row.stop_1based} {row.region}\n{row.matched_rna_sequence}\n")
    summary = {"comparison": str(frame.comparison_id.iloc[0]), "policy": {**POLICY, "site_pvalue": threshold},
        "selected_SE": len(selected_ids), "selected_genes": frame[frame.event_id.isin(selected_ids)].gene_id.nunique(),
        "matching_status": pairs.status.value_counts().to_dict(), "matched_by_direction": matched.direction.value_counts().to_dict(),
        "motif_assignments": len(catalog), "RBPs": catalog.RBP.str.casefold().nunique(),
        "hypotheses": len(enrichment), "significant_BH_005": int(enrichment.enrichment_fdr.lt(.05).sum()),
        "significant_by_evidence": enrichment[enrichment.enrichment_fdr.lt(.05)].evidence_class.value_counts().to_dict(),
        "foreground_hit_assignments": len(foreground), "foreground_distinct_sites": len(foreground.drop_duplicates(["scan_token", "window_id", "start_1based", "stop_1based"])),
        "fimo": provenance, "quantification_rerun": False, "CLIP_analysis": "not computed in this sequence-consequence run"}
    (directory/"RESULTS.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return summary
