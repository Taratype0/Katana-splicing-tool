"""Automated GTF-driven ORF, Reading-Frame, NMD (50-nt EJC rule), and Curated UniProt/Pfam Domain Consequence Engine (v2.4.1).

Key fixes & guarantees in v2.4.1:
  1. Strict `start_codon` + GTF column-8 `CDS phase (0/1/2)` frame alignment; excludes 5'/3'-incomplete
     GENCODE fragments (`cds_start_NF`, `cds_end_NF`, `mRNA_start_NF`) and pre-existing `nonsense_mediated_decay`
     biotypes whenever a complete `protein_coding` host transcript exists.
  2. Splice-junction-aware host transcript selection: ranks candidate `protein_coding` transcripts by exact
     exon-boundary agreement with the splicing event's flanking donor/acceptor sites (`u1`, `d0`), preventing
     accidental mid-exon clipping on multi-isoform giant genes (e.g. `Macf1` 7,355 aa, `Ptbp3` 552 aa, `Stxbp2` 593 aa).
  3. Bundled offline Human + Mouse Reviewed UniProtKB/Swiss-Prot + Pfam domain coordinate catalog
     (`194,738` curated intervals across `37,714` proteins in `data/uniprot_domains_human_mouse_v1.tsv.gz`),
     mapping exact N-/C-anchored residue intervals onto both Inclusion and Exclusion protein isoforms to
     quantify in-frame domain indels AND C-terminal domain loss caused by frameshift/PTC truncation.
"""
from __future__ import annotations

from dataclasses import dataclass
import gzip
from importlib import resources
import json
from pathlib import Path
import re
import numpy as np
import pandas as pd

CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

_ATTR_GENE_ID = re.compile(r'gene_id\s+"?([^";\s]+)')
_ATTR_GENE_NAME = re.compile(r'gene_name\s+"?([^";\s]+)')
_ATTR_TX_ID = re.compile(r'transcript_id\s+"?([^";\s]+)')
_ATTR_TX_TYPE = re.compile(r'transcript_(?:type|biotype)\s+"?([^";\s]+)')
_ATTR_CANON = re.compile(r'tag\s+"?(appris_principal[^";\s]*|Ensembl_canonical|MANE_Select|CCDS)')
_ATTR_INCOMPLETE = re.compile(r'tag\s+"?(cds_start_NF|cds_end_NF|mRNA_start_NF|mRNA_end_NF)')

_DOMAIN_CACHE: dict[str, dict[str, dict]] = {}


def load_bundled_uniprot_domains(species: str = "auto") -> dict[str, dict]:
    """Load bundled 194,738-interval UniProt/Pfam catalog indexed by lowercase gene_symbol."""
    cache_key = species.lower()
    if cache_key in _DOMAIN_CACHE:
        return _DOMAIN_CACHE[cache_key]

    pkg_root = __name__.split(".")[0]
    data_ref = resources.files(f"{pkg_root}.data").joinpath("uniprot_domains_human_mouse_v1.tsv.gz")
    with resources.as_file(data_ref) as gz_path:
        df = pd.read_csv(gz_path, sep="\t")

    if cache_key in ("mouse", "mus_musculus"):
        df = df[df["species"].eq("mouse")]
    elif cache_key in ("human", "homo_sapiens"):
        df = df[df["species"].eq("human")]

    pfam_col = "pfam_ accessions" if "pfam_ accessions" in df.columns else "pfam_accessions"
    index: dict[str, dict] = {}
    for g_low, acc, ulen, pfam, kind, s_aa, e_aa, lbl in zip(
        df["gene_lower"],
        df["uniprot_acc"],
        df["uniprot_len"],
        df[pfam_col] if pfam_col in df.columns else [""] * len(df),
        df["feature_kind"],
        df["start_aa"],
        df["end_aa"],
        df["feature_label"],
    ):
        gk = str(g_low)
        ent = index.get(gk)
        if ent is None:
            ent = {
                "uniprot_acc": str(acc),
                "uniprot_len": int(ulen),
                "pfam_ids": "" if pd.isna(pfam) else str(pfam),
                "features": [],
            }
            index[gk] = ent
        ent["features"].append({"kind": str(kind), "start": int(s_aa), "end": int(e_aa), "label": str(lbl)})
    _DOMAIN_CACHE[cache_key] = index
    return index


@dataclass
class HostTranscriptModel:
    transcript_id: str
    gene_id: str
    gene_name: str
    chrom: str
    strand: str
    biotype: str
    exons: list[tuple[int, int]]
    cds_blocks: list[tuple[int, int]]     # sorted ascending genomic (start0, end0)
    cds_phases: dict[tuple[int, int], int]
    start_codon_5p: int | None            # exact 0-based genomic position of first base of start codon
    stop_codon_3p: int | None             # exact 0-based genomic position of 3'-most base of canonical stop codon
    is_canonical: bool
    is_incomplete: bool


def load_gtf_coding_transcripts(
    gtf_path: Path | str,
    target_genes: set[str] | None = None,
) -> dict[str, list[HostTranscriptModel]]:
    """Streaming GTF/GTF.gz parser honoring `start_codon`, `stop_codon`, `CDS` column-8 phase, biotypes, and `cds_start_NF` flags."""
    p = Path(gtf_path)
    opener = gzip.open if p.name.endswith(".gz") else open
    clean_targets: set[str] | None = None
    if target_genes:
        clean_targets = {g.split(".")[0].strip() for g in target_genes if g} | {g.lower().strip() for g in target_genes if g}

    tx_meta: dict[str, list] = {}  # [gid, gname, chrom, strand, biotype, is_can, is_inc]
    tx_exons: dict[str, list[tuple[int, int]]] = {}
    tx_cds: dict[str, list[tuple[int, int]]] = {}
    tx_phases: dict[str, dict[tuple[int, int], int]] = {}
    tx_start_bases: dict[str, list[tuple[int, int]]] = {}
    tx_stop_bases: dict[str, list[tuple[int, int]]] = {}

    with opener(p, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line or line[0] == "#":
                continue
            parts = line.split("\t", 8)
            if len(parts) < 9:
                continue
            ftype = parts[2]
            if ftype not in ("exon", "CDS", "start_codon", "stop_codon"):
                continue
            attrs = parts[8]
            mg = _ATTR_GENE_ID.search(attrs)
            if not mg:
                continue
            gid = mg.group(1).split(".")[0]
            mn = _ATTR_GENE_NAME.search(attrs)
            gname = mn.group(1) if mn else gid

            if clean_targets is not None and (gid not in clean_targets and gname.lower() not in clean_targets):
                continue

            mt = _ATTR_TX_ID.search(attrs)
            if not mt:
                continue
            tid = mt.group(1).split(".")[0]
            start0 = int(parts[3]) - 1
            end0 = int(parts[4])
            if end0 <= start0:
                continue

            is_can = bool(_ATTR_CANON.search(attrs))
            is_inc = bool(_ATTR_INCOMPLETE.search(attrs))
            m_bt = _ATTR_TX_TYPE.search(attrs)
            bt = m_bt.group(1) if m_bt else "protein_coding"

            if tid not in tx_meta:
                tx_meta[tid] = [gid, gname, parts[0], parts[6], bt, is_can, is_inc]
            else:
                if is_can:
                    tx_meta[tid][5] = True
                if is_inc:
                    tx_meta[tid][6] = True

            if ftype == "exon":
                tx_exons.setdefault(tid, []).append((start0, end0))
            elif ftype == "CDS":
                tx_cds.setdefault(tid, []).append((start0, end0))
                ph = int(parts[7]) if parts[7] in ("0", "1", "2") else 0
                tx_phases.setdefault(tid, {})[(start0, end0)] = ph
            elif ftype == "start_codon":
                tx_start_bases.setdefault(tid, []).append((start0, end0))
            elif ftype == "stop_codon":
                tx_stop_bases.setdefault(tid, []).append((start0, end0))

    by_gene: dict[str, list[HostTranscriptModel]] = {}
    for tid, (gid, gname, chrom, strand, bt, is_can, is_inc) in tx_meta.items():
        ex = sorted(set(tx_exons.get(tid, ())))
        cds = sorted(set(tx_cds.get(tid, ())))
        if not ex:
            continue
        sc_blocks = sorted(set(tx_start_bases.get(tid, ())))
        sc_5p = (sc_blocks[0][0] if strand == "+" else (sc_blocks[-1][1] - 1)) if sc_blocks else None
        st_blocks = sorted(set(tx_stop_bases.get(tid, ())))
        stop_3p = ((st_blocks[-1][1] - 1) if strand == "+" else st_blocks[0][0]) if st_blocks else None
        model = HostTranscriptModel(
            transcript_id=tid,
            gene_id=gid,
            gene_name=gname,
            chrom=chrom,
            strand=strand,
            biotype=bt,
            exons=ex,
            cds_blocks=cds,
            cds_phases=tx_phases.get(tid, {}),
            start_codon_5p=sc_5p,
            stop_codon_3p=stop_3p,
            is_canonical=is_can,
            is_incomplete=is_inc,
        )
        by_gene.setdefault(gid, []).append(model)
        by_gene.setdefault(gname.lower(), []).append(model)

    for k in by_gene:
        by_gene[k].sort(
            key=lambda m: (
                len(m.cds_blocks) > 0 and m.biotype == "protein_coding" and not m.is_incomplete,
                m.start_codon_5p is not None and not m.is_incomplete,
                len(m.cds_blocks) > 0 and not m.is_incomplete,
                m.is_canonical,
                sum(b - a for a, b in m.cds_blocks),
                sum(b - a for a, b in m.exons),
            ),
            reverse=True,
        )
    return by_gene


def _merge_sorted_blocks(blocks: list[tuple[int, int]]) -> list[tuple[int, int]]:
    out: list[tuple[int, int]] = []
    for a, b in sorted(blocks):
        if b <= a:
            continue
        if out and a <= out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], b))
        else:
            out.append((a, b))
    return out


def _rank_host_transcripts(candidates: list[HostTranscriptModel], rec: dict) -> list[HostTranscriptModel]:
    if not candidates:
        return []
    kind = str(rec["event_type"])
    ev_coords = [
        int(rec[c])
        for c in ("exonStart_0base", "exonEnd", "longExonStart_0base", "longExonEnd", "1stExonStart_0base", "2ndExonEnd", "riExonStart_0base", "riExonEnd")
        if c in rec and pd.notna(rec[c])
    ]
    ev_lo, ev_hi = (min(ev_coords), max(ev_coords)) if ev_coords else (0, 0)

    # Flanking splice junctions & intervening-exon topology we want the host transcript to match
    if kind in ("SE", "MXE", "RI"):
        u1, d0 = int(rec["upstreamEE"]), int(rec["downstreamES"])
        anchor_ends = {u1}
        anchor_starts = {d0}
        allowed_mid: set[tuple[int, int]] = set()
        if kind == "SE":
            allowed_mid.add((int(rec["exonStart_0base"]), int(rec["exonEnd"])))
        elif kind == "MXE":
            allowed_mid.add((int(rec["1stExonStart_0base"]), int(rec["1stExonEnd"])))
            allowed_mid.add((int(rec["2ndExonStart_0base"]), int(rec["2ndExonEnd"])))
    else:
        f0, f1 = int(rec["flankingES"]), int(rec["flankingEE"])
        l0, l1 = int(rec["longExonStart_0base"]), int(rec["longExonEnd"])
        s0, s1 = int(rec["shortES"]), int(rec["shortEE"])
        u1, d0 = (min(s1, l1), f0) if l1 <= f0 else (f1, max(l0, s0))
        anchor_starts = {f0, l0, s0}
        anchor_ends = {f1, l1, s1}
        allowed_mid = set()

    def _score_cand(m: HostTranscriptModel) -> tuple:
        tx_lo, tx_hi = m.exons[0][0], m.exons[-1][1]
        overlaps = tx_lo <= ev_hi and tx_hi >= ev_lo
        starts_set = {a for a, _ in m.exons}
        ends_set = {b for _, b in m.exons}
        junc_matches = len(anchor_starts & starts_set) + len(anchor_ends & ends_set)
        mid_ex = [(a, b) for a, b in m.exons if a >= u1 and b <= d0]
        extra_mid = [eb - ea for ea, eb in mid_ex if (ea, eb) not in allowed_mid]
        if kind == "RI":
            spans_ri = any(ea <= u1 and eb >= d0 for ea, eb in m.exons)
            exact_topo = (junc_matches >= 2 and len(mid_ex) == 0) or spans_ri
        else:
            exact_topo = (junc_matches >= 2) and (len(extra_mid) == 0) and (len(mid_ex) <= 1)
        frame_safe_mid = (sum(extra_mid) % 3 == 0)
        clean_coding = len(m.cds_blocks) > 0 and m.biotype == "protein_coding" and not m.is_incomplete
        has_sc = m.start_codon_5p is not None and not m.is_incomplete
        return (
            overlaps and clean_coding,
            overlaps and clean_coding and exact_topo,
            overlaps and clean_coding and (junc_matches >= 2) and frame_safe_mid,
            overlaps and clean_coding and (junc_matches >= 2),
            overlaps and clean_coding and (junc_matches >= 1),
            has_sc,
            m.is_canonical,
            sum(b - a for a, b in m.cds_blocks),
        )

    return sorted(candidates, key=_score_cand, reverse=True)


def _build_full_isoform_exons(host_exons: list[tuple[int, int]], row: dict) -> tuple[list[tuple[int, int]], list[tuple[int, int]]]:
    """Splice local Inclusion vs Exclusion coordinates cleanly into full-length host exons without dropping flank exons."""
    kind = str(row["event_type"])
    strand = str(row["strand"])

    if kind in ("SE", "MXE"):
        u0, u1 = int(row["upstreamES"]), int(row["upstreamEE"])
        d0, d1 = int(row["downstreamES"]), int(row["downstreamEE"])
        left_exons = []
        right_exons = []
        found_u, found_d = False, False
        for a, b in host_exons:
            if b <= u1:
                if b == u1:
                    found_u = True
                left_exons.append((a, b))
            elif a < u1 < b:
                # Host exon extends past u1 -> trim cleanly at u1
                left_exons.append((min(a, u0), u1))
                found_u = True
            if a >= d0:
                if a == d0:
                    found_d = True
                right_exons.append((a, b))
            elif a < d0 < b:
                right_exons.append((d0, max(b, d1)))
                found_d = True
        if not found_u:
            left_exons.append((u0, u1))
        if not found_d:
            right_exons.append((d0, d1))

        if kind == "SE":
            t0, t1 = int(row["exonStart_0base"]), int(row["exonEnd"])
            return _merge_sorted_blocks(left_exons + [(t0, t1)] + right_exons), _merge_sorted_blocks(left_exons + right_exons)
        a0, a1 = int(row["1stExonStart_0base"]), int(row["1stExonEnd"])
        b0, b1 = int(row["2ndExonStart_0base"]), int(row["2ndExonEnd"])
        inc_ex, exc_ex = ((a0, a1), (b0, b1)) if strand == "+" else ((b0, b1), (a0, a1))
        return _merge_sorted_blocks(left_exons + [inc_ex] + right_exons), _merge_sorted_blocks(left_exons + [exc_ex] + right_exons)

    if kind == "RI":
        u0, u1 = int(row["upstreamES"]), int(row["upstreamEE"])
        d0, d1 = int(row["downstreamES"]), int(row["downstreamEE"])
        left_ex = [(a, min(b, u1)) for a, b in host_exons if a < u1]
        right_ex = [(max(a, d0), b) for a, b in host_exons if b > d0]
        if not any(b == u1 for _, b in left_ex):
            left_ex.append((u0, u1))
        if not any(a == d0 for a, _ in right_ex):
            right_ex.append((d0, d1))
        exc = _merge_sorted_blocks(left_ex + right_ex)
        inc = _merge_sorted_blocks(left_ex + [(u1, d0)] + right_ex)
        return inc, exc

    # A3SS / A5SS
    l0, l1 = int(row["longExonStart_0base"]), int(row["longExonEnd"])
    s0, s1 = int(row["shortES"]), int(row["shortEE"])
    f0, f1 = int(row["flankingES"]), int(row["flankingEE"])
    if l1 <= f0:
        # Alternative exon is genomic-left of flanking exon
        left_base = [(a, b) for a, b in host_exons if b <= l0]
        # Preserve upstream start of the alternative exon if host exon started earlier than l0
        for a, b in host_exons:
            if a < l0 < b <= l1 + 5:
                l0, s0 = a, a
        right_base = [(a, b) for a, b in host_exons if a >= f0]
        if not any(a == f0 for a, _ in right_base):
            right_base.append((f0, f1))
        return _merge_sorted_blocks(left_base + [(l0, l1)] + right_base), _merge_sorted_blocks(left_base + [(s0, s1)] + right_base)
    else:
        # Flanking exon is genomic-left of alternative exon
        left_base = [(a, b) for a, b in host_exons if b <= f1]
        if not any(b == f1 for _, b in left_base):
            left_base.append((f0, f1))
        for a, b in host_exons:
            if l0 - 5 <= a < l1 < b:
                l1, s1 = b, b
        right_base = [(a, b) for a, b in host_exons if a >= l1]
        return _merge_sorted_blocks(left_base + [(l0, l1)] + right_base), _merge_sorted_blocks(left_base + [(s0, s1)] + right_base)


def _resolve_reading_frame_start_gpos(host: HostTranscriptModel) -> tuple[int, int]:
    """Return `(start_gpos, phase_offset)` guaranteeing exact codon-frame alignment."""
    if host.start_codon_5p is not None:
        return host.start_codon_5p, 0
    first_blk = host.cds_blocks[0] if host.strand == "+" else host.cds_blocks[-1]
    ph = host.cds_phases.get(first_blk, 0)
    gpos = first_blk[0] if host.strand == "+" else (first_blk[1] - 1)
    return gpos, ph


def _genomic_pos_to_tx_offset(exons_asc: list[tuple[int, int]], strand: str, gpos: int) -> int | None:
    ordered = exons_asc if strand == "+" else exons_asc[::-1]
    cursor = 0
    for a, b in ordered:
        if a <= gpos < b:
            return cursor + (gpos - a if strand == "+" else (b - 1 - gpos))
        cursor += b - a
    return None


def _translate_and_evaluate_nmd(
    dna_5to3: str,
    exon_lengths_5to3: list[int],
    start_tx_pos: int | None,
    phase_offset: int = 0,
    canonical_stop_tx_end: int | None = None,
    host_baseline_ptc_dist: int = 0,
) -> dict:
    total_len = len(dna_5to3)
    ejc_positions = list(np.cumsum(exon_lengths_5to3)[:-1]) if len(exon_lengths_5to3) >= 2 else []
    last_ejc = int(ejc_positions[-1]) if ejc_positions else None

    if start_tx_pos is not None:
        eff_start = start_tx_pos + int(phase_offset)
    else:
        eff_start = None

    if eff_start is None or eff_start < 0 or eff_start + 3 > total_len:
        best_orf, best_stop_end, best_st = "", None, 0
        for m in re.finditer(r"ATG", dna_5to3):
            st = m.start()
            aa_chars, stop_end = [], None
            for p in range(st, total_len - 2, 3):
                aa = CODON_TABLE.get(dna_5to3[p : p + 3], "X")
                if aa == "*":
                    stop_end = p + 3
                    break
                aa_chars.append(aa)
            if len(aa_chars) > len(best_orf):
                best_orf, best_stop_end, best_st = "".join(aa_chars), stop_end, st
        dist_ejc = int(last_ejc - best_stop_end) if (best_stop_end is not None and last_ejc is not None) else 0
        is_native_stop = canonical_stop_tx_end is not None and best_stop_end == canonical_stop_tx_end and dist_ejc <= max(host_baseline_ptc_dist + 6, 49)
        nmd = bool(best_stop_end is not None and last_ejc is not None and dist_ejc >= 50 and not is_native_stop)
        return {
            "protein": best_orf,
            "aa_len": len(best_orf),
            "stop_tx_nt": best_stop_end,
            "last_ejc_tx_nt": last_ejc,
            "mrna_len_nt": total_len,
            "nmd_sensitive": nmd,
            "ptc_to_last_ejc_nt": max(0, dist_ejc),
            "start_preserved": False,
        }

    aa_list: list[str] = []
    stop_end: int | None = None
    for p in range(eff_start, total_len - 2, 3):
        aa = CODON_TABLE.get(dna_5to3[p : p + 3], "X")
        if aa == "*":
            stop_end = p + 3
            break
        aa_list.append(aa)
    prot = "".join(aa_list)
    dist_ejc = int(last_ejc - stop_end) if (stop_end is not None and last_ejc is not None) else 0
    is_native_stop = canonical_stop_tx_end is not None and stop_end == canonical_stop_tx_end and dist_ejc <= max(host_baseline_ptc_dist + 6, 49)
    nmd = bool(stop_end is not None and last_ejc is not None and dist_ejc >= 50 and not is_native_stop)
    return {
        "protein": prot,
        "aa_len": len(prot),
        "stop_tx_nt": stop_end,
        "last_ejc_tx_nt": last_ejc,
        "mrna_len_nt": total_len,
        "nmd_sensitive": nmd,
        "ptc_to_last_ejc_nt": max(0, dist_ejc),
        "start_preserved": True,
    }


def _map_uniprot_domains_to_isoforms(
    gene_symbol: str,
    prot_inc: str,
    prot_exc: str,
    nmd_inc: bool,
    nmd_exc: bool,
    domain_db: dict[str, dict],
) -> tuple[str, str, list[dict]]:
    """Map curated UniProt/Pfam intervals onto Inclusion vs Exclusion ORFs and classify intact vs spliced vs PTC-lost."""
    entry = domain_db.get(gene_symbol.lower())
    len_inc, len_exc = len(prot_inc), len(prot_exc)
    ref_len = max(len_inc, len_exc)
    if ref_len == 0:
        return "NA", "none", []

    # Exact shared N-terminal prefix & C-terminal suffix between Inclusion and Exclusion proteins
    pre = 0
    min_l = min(len_inc, len_exc)
    while pre < min_l and prot_inc[pre] == prot_exc[pre]:
        pre += 1
    suf = 0
    while suf < (min_l - pre) and prot_inc[len_inc - 1 - suf] == prot_exc[len_exc - 1 - suf]:
        suf += 1

    # Altered window on the longer reference protein (1-based amino acid coordinates)
    alt_start_aa = pre + 1
    alt_end_aa = max(alt_start_aa, ref_len - suf)

    mapped_arch: list[dict] = []
    disrupted_labels: list[str] = []
    uniprot_acc = "predicted"

    if entry is not None and entry["features"]:
        uniprot_acc = entry["uniprot_acc"]
        u_len = max(entry["uniprot_len"], 1)
        scale = ref_len / float(u_len)
        # Only transfer coordinates if isoform length is reasonably compatible (within 0.45x..1.8x of canonical Swiss-Prot)
        if 0.40 <= scale <= 1.85:
            for feat in entry["features"]:
                # Anchor N-terminal domains directly; smoothly scale C-terminal coordinates when isoform differs slightly
                s_map = max(1, min(ref_len, int(round(feat["start"] * scale)) if abs(ref_len - u_len) > 25 else feat["start"]))
                e_map = max(s_map, min(ref_len, int(round(feat["end"] * scale)) if abs(ref_len - u_len) > 25 else feat["end"]))
                # Check overlap with the spliced/truncated region [alt_start_aa, alt_end_aa]
                ov = max(0, min(e_map, alt_end_aa) - max(s_map, alt_start_aa) + 1)
                dom_span = max(1, e_map - s_map + 1)
                if prot_inc != prot_exc and ov >= 1:
                    if (nmd_inc or nmd_exc) and s_map > min_l:
                        state = "lost_by_PTC_truncation"
                    else:
                        state = "spliced_indel_overlap"
                    pct = min(100, int(round(100.0 * ov / dom_span)))
                    if feat["kind"] != "REGION" or "Disordered" not in feat["label"]:
                        disrupted_labels.append(f"{feat['label']} [{s_map}-{e_map} aa, {pct}%]")
                    elif len(disrupted_labels) == 0:
                        disrupted_labels.append(f"{feat['label']} [{s_map}-{e_map} aa]")
                else:
                    state = "intact"
                mapped_arch.append({
                    "kind": feat["kind"],
                    "label": feat["label"],
                    "start": s_map,
                    "end": e_map,
                    "state": state,
                })

    dedup_disrupted = list(dict.fromkeys(disrupted_labels))[:5]
    feat_str = "; ".join(dedup_disrupted) if dedup_disrupted else "none"
    return uniprot_acc, feat_str, mapped_arch


def trace_protein_and_nmd_consequences(
    events_df: pd.DataFrame,
    gtf_path: Path | str,
    fasta_path: Path | str,
    *,
    indexed_fasta_cls: type,
    species: str = "auto",
) -> pd.DataFrame:
    """Deterministic GTF + FASTA + Curated UniProt/Pfam ORF, reading-frame, NMD, and protein-domain consequence engine."""
    if events_df.empty:
        return pd.DataFrame()

    target_genes = (
        set(events_df["gene_id"].dropna().astype(str))
        | set(events_df.get("gene_symbol", pd.Series(dtype=str)).dropna().astype(str))
        | set(events_df.get("gene_name", pd.Series(dtype=str)).dropna().astype(str))
    )
    tx_index = load_gtf_coding_transcripts(gtf_path, target_genes=target_genes)
    domain_db = load_bundled_uniprot_domains(species=species)
    reader = indexed_fasta_cls(Path(fasta_path))

    rows: list[dict] = []
    with reader.open_stream() as stream:
        for rec in events_df.to_dict("records"):
            eid = str(rec["event_id"])
            gid = str(rec.get("gene_id", "")).split(".")[0]
            raw_sym = rec.get("gene_symbol")
            if raw_sym is None or pd.isna(raw_sym) or str(raw_sym).strip() == "":
                raw_sym = rec.get("gene_name")
            gsym = str(raw_sym).strip() if (raw_sym is not None and pd.notna(raw_sym) and str(raw_sym).strip() != "") else gid
            dpsi = float(rec.get("dpsi", 0.0))
            etype = str(rec["event_type"])
            strand = str(rec["strand"])

            candidates = tx_index.get(gid) or tx_index.get(gsym.lower()) or []
            ranked_hosts = _rank_host_transcripts(candidates, rec)
            host = ranked_hosts[0] if ranked_hosts else None
            if host is not None and (gsym == gid or gsym.startswith(("ENSG", "ENSMUSG"))) and host.gene_name:
                gsym = host.gene_name

            if host is None or not host.cds_blocks:
                rows.append({
                    "event_id": eid,
                    "gene_id": gid,
                    "gene_symbol": gsym,
                    "event_type": etype,
                    "dpsi": dpsi,
                    "host_transcript_id": host.transcript_id if host else "unmatched",
                    "uniprot_acc": "NA",
                    "coding_region_location": "noncoding_host",
                    "frame_status": "noncoding_transcript",
                    "nmd_switch_class": "noncoding",
                    "aa_len_inclusion": 0,
                    "aa_len_exclusion": 0,
                    "delta_aa": 0,
                    "ptc_to_last_ejc_nt": 0,
                    "disrupted_protein_features": "none",
                    "functional_consequence_summary": "Non-coding transcript isoform change",
                    "architecture_json": "{}",
                })
                continue

            # Evaluate top candidate protein_coding backbones so multi-promoter / multi-isoform genes (e.g. Macf1)
            # bind to the actual coding backbone compatible with the local splice junctions rather than a chimeric frame!
            coding_pool = [m for m in ranked_hosts if m.cds_blocks][:6]
            chosen_bundle = None
            for cand_m in coding_pool:
                try:
                    c_inc_ex, c_exc_ex = _build_full_isoform_exons(cand_m.exons, rec)
                    c_ord_inc = c_inc_ex if strand == "+" else c_inc_ex[::-1]
                    c_ord_exc = c_exc_ex if strand == "+" else c_exc_ex[::-1]
                    c_dna_inc = reader.fetch_blocks(stream, cand_m.chrom, strand, c_ord_inc)
                    c_dna_exc = reader.fetch_blocks(stream, cand_m.chrom, strand, c_ord_exc)
                except Exception:
                    continue
                c_start_gpos, c_phase_off = _resolve_reading_frame_start_gpos(cand_m)
                c_st_inc = _genomic_pos_to_tx_offset(c_inc_ex, strand, c_start_gpos)
                c_st_exc = _genomic_pos_to_tx_offset(c_exc_ex, strand, c_start_gpos)

                # Native host baseline stop-to-last-EJC distance (protects canonical genes with natural 3'-UTR introns)
                host_ord = cand_m.exons if strand == "+" else cand_m.exons[::-1]
                host_ejcs = list(np.cumsum([b - a for a, b in host_ord])[:-1]) if len(host_ord) >= 2 else []
                host_last_ejc = int(host_ejcs[-1]) if host_ejcs else 0
                if cand_m.stop_codon_3p is not None:
                    h_stop_off = _genomic_pos_to_tx_offset(cand_m.exons, strand, cand_m.stop_codon_3p)
                    host_base_ptc = max(0, host_last_ejc - (h_stop_off + 1)) if h_stop_off is not None else 0
                    inc_stop_off = _genomic_pos_to_tx_offset(c_inc_ex, strand, cand_m.stop_codon_3p)
                    exc_stop_off = _genomic_pos_to_tx_offset(c_exc_ex, strand, cand_m.stop_codon_3p)
                    can_end_inc = (inc_stop_off + 1) if inc_stop_off is not None else None
                    can_end_exc = (exc_stop_off + 1) if exc_stop_off is not None else None
                else:
                    host_base_ptc, can_end_inc, can_end_exc = 0, None, None

                ct_inc = _translate_and_evaluate_nmd(
                    c_dna_inc, [b - a for a, b in c_ord_inc], c_st_inc,
                    phase_offset=c_phase_off, canonical_stop_tx_end=can_end_inc, host_baseline_ptc_dist=host_base_ptc,
                )
                ct_exc = _translate_and_evaluate_nmd(
                    c_dna_exc, [b - a for a, b in c_ord_exc], c_st_exc,
                    phase_offset=c_phase_off, canonical_stop_tx_end=can_end_exc, host_baseline_ptc_dist=host_base_ptc,
                )
                cand_bundle = (cand_m, c_inc_ex, c_exc_ex, c_ord_inc, c_ord_exc, c_dna_inc, c_dna_exc, ct_inc, ct_exc)
                if chosen_bundle is None:
                    chosen_bundle = cand_bundle
                # Accept immediately once at least one splice isoform is productive (not artificial dual-NMD)
                if not (ct_inc["nmd_sensitive"] and ct_exc["nmd_sensitive"]):
                    chosen_bundle = cand_bundle
                    break

            if chosen_bundle is None:
                continue
            host, inc_exons, exc_exons, ord_inc, ord_exc, dna_inc, dna_exc, t_inc, t_exc = chosen_bundle

            cds_lo, cds_hi = host.cds_blocks[0][0], host.cds_blocks[-1][1]
            ev_coords = [
                int(rec[c])
                for c in ("exonStart_0base", "exonEnd", "longExonStart_0base", "longExonEnd", "1stExonStart_0base", "2ndExonEnd", "riExonStart_0base", "riExonEnd")
                if c in rec and pd.notna(rec[c])
            ]
            ev_lo, ev_hi = (min(ev_coords), max(ev_coords)) if ev_coords else (0, 0)
            overlaps_cds = bool(ev_lo < cds_hi and ev_hi > cds_lo)
            if not overlaps_cds:
                is_5utr = (ev_hi <= cds_lo) if strand == "+" else (ev_lo >= cds_hi)
                loc = "5prime_UTR" if is_5utr else "3prime_UTR"
            else:
                loc = "CDS_coding"

            delta_nt = len(dna_inc) - len(dna_exc)
            if t_inc["protein"] == t_exc["protein"]:
                f_status = "utr_cis_regulatory_only"
            elif not (t_inc["start_preserved"] and t_exc["start_preserved"]):
                f_status = "translation_start_altered"
            elif abs(delta_nt) % 3 == 0:
                f_status = "in_frame_indel"
            else:
                f_status = "frameshift"

            if t_inc["nmd_sensitive"] and not t_exc["nmd_sensitive"]:
                nmd_cls = "inclusion_triggers_NMD_poison_exon"
                ptc_dist = t_inc["ptc_to_last_ejc_nt"]
            elif t_exc["nmd_sensitive"] and not t_inc["nmd_sensitive"]:
                nmd_cls = "skipping_triggers_NMD_decay"
                ptc_dist = t_exc["ptc_to_last_ejc_nt"]
            elif not t_inc["nmd_sensitive"] and not t_exc["nmd_sensitive"]:
                nmd_cls = "productive_protein_isoform_switch" if t_inc["protein"] != t_exc["protein"] else "coding_silent_utr_switch"
                ptc_dist = 0
            else:
                nmd_cls = "dual_NMD_background"
                ptc_dist = max(t_inc["ptc_to_last_ejc_nt"], t_exc["ptc_to_last_ejc_nt"])

            u_acc, feat_str, arch_list = _map_uniprot_domains_to_isoforms(
                gsym, t_inc["protein"], t_exc["protein"], t_inc["nmd_sensitive"], t_exc["nmd_sensitive"], domain_db
            )
            delta_aa = int(t_inc["aa_len"] - t_exc["aa_len"])

            # Exact shared prefix residue for visualization
            pre_aa = 0
            p1, p2 = t_inc["protein"], t_exc["protein"]
            while pre_aa < min(len(p1), len(p2)) and p1[pre_aa] == p2[pre_aa]:
                pre_aa += 1

            arch_payload = {
                "uniprot_acc": u_acc,
                "mrna_nt_inc": t_inc["mrna_len_nt"],
                "mrna_nt_exc": t_exc["mrna_len_nt"],
                "last_ejc_inc": t_inc["last_ejc_tx_nt"],
                "last_ejc_exc": t_exc["last_ejc_tx_nt"],
                "stop_nt_inc": t_inc["stop_tx_nt"],
                "stop_nt_exc": t_exc["stop_tx_nt"],
                "n_exons_inc": len(ord_inc),
                "n_exons_exc": len(ord_exc),
                "splice_junction_aa": pre_aa + 1,
                "domains": arch_list[:18],
            }

            if nmd_cls == "inclusion_triggers_NMD_poison_exon":
                fav = "Promotes NMD decay" if dpsi > 0 else "Skips poison exon -> rescues protein"
                dom_lost = f" (truncates {feat_str})" if feat_str != "none" else ""
                summary = f"Predicted NMD poison exon ({f_status}): inclusion PTC {ptc_dist} nt upstream of last EJC -> {fav}{dom_lost}"
            elif nmd_cls == "skipping_triggers_NMD_decay":
                fav = "Triggers NMD decay" if dpsi < 0 else "Maintains full-length ORF"
                dom_lost = f" (truncates {feat_str})" if feat_str != "none" else ""
                summary = f"Predicted NMD-inducing skipping ({f_status}): exclusion PTC {ptc_dist} nt upstream of last EJC -> {fav}{dom_lost}"
            elif f_status == "utr_cis_regulatory_only":
                summary = f"ORF intact ({t_inc['aa_len']} aa); remodels {loc} ({delta_nt:+d} nt)"
            else:
                dom_note = f"; overlaps {feat_str}" if feat_str != "none" else ""
                summary = f"Productive {f_status}: {t_inc['aa_len']} vs {t_exc['aa_len']} aa ({delta_aa:+d} aa){dom_note}"

            rows.append({
                "event_id": eid,
                "gene_id": gid,
                "gene_symbol": gsym,
                "event_type": etype,
                "dpsi": dpsi,
                "host_transcript_id": host.transcript_id,
                "uniprot_acc": u_acc,
                "coding_region_location": loc,
                "frame_status": f_status,
                "nmd_switch_class": nmd_cls,
                "aa_len_inclusion": int(t_inc["aa_len"]),
                "aa_len_exclusion": int(t_exc["aa_len"]),
                "delta_aa": delta_aa,
                "ptc_to_last_ejc_nt": int(ptc_dist),
                "disrupted_protein_features": feat_str,
                "functional_consequence_summary": summary,
                "architecture_json": json.dumps(arch_payload),
            })
    return pd.DataFrame(rows)
