"""Threshold-free PWM preference and native-consensus density sensitivity.

Average motif odds follows the sequence score documented by MEME AME; the
paired sign test here is Katana's declared matched-gene test, NOT native AME
or rMAPS2 execution. Nothing in this module refits DAS or quantifies reads.
"""
from __future__ import annotations

import hashlib
import io
import json
import re
import tarfile
from pathlib import Path
from zipfile import ZipFile

import numpy as np
import pandas as pd
from scipy.stats import binomtest

from katana_cosplice.motif_scanner import IUPAC, Motif, load_motif_database
from katana_next.science.se_rbp import REGIONS
from katana_next.science.statistics import bh_adjust

ALIASES = {"srp20": "SRSF3", "hur": "ELAVL1", "ptb": "PTBP1", "sf2-asf": "SRSF1",
           "srp40": "SRSF5", "srp55": "SRSF6", "9g8": "SRSF7", "tra2-beta": "TRA2B",
           "brunol4": "CELF4", "brunol5": "CELF5", "brunol6": "CELF6", "a2bp1": "RBFOX1"}


def canonical_rbp(name):
    return ALIASES.get(str(name).casefold(), str(name))


def fixed_pattern(text):
    """Accept finite native IUPAC or bracket patterns, never arbitrary regex."""
    text = str(text).upper().replace("U", "T")
    tokens, i = [], 0
    while i < len(text):
        if text[i] == "[":
            j = text.find("]", i)
            if j < 0 or not text[i+1:j] or set(text[i+1:j])-set("ACGT"):
                raise ValueError("Invalid finite base class")
            tokens.append("["+"".join(sorted(set(text[i+1:j])))+"]")
            i = j+1
        elif text[i] in IUPAC:
            tokens.append(IUPAC[text[i]])
            i += 1
        else:
            raise ValueError("Variable-length/probe notation is not a fixed native motif")
    if not tokens:
        raise ValueError("Empty motif")
    return "".join(tokens), len(tokens)


def consensus_pwm(pattern):
    """Schematic only: equal weights within allowed bases, NOT measured PWM."""
    regex, width = fixed_pattern(pattern)
    tokens = re.findall(r"\[[ACGT]+\]|[ACGT]", regex)
    matrix = np.zeros((width, 4))
    for i, token in enumerate(tokens):
        bases = token.strip("[]")
        for base in bases:
            matrix[i, "ACGT".index(base)] = 1/len(bases)
    return matrix


def expanded_catalog(manifest, resources):
    """Use original numeric PWMs and separately identified public references."""
    resources = Path(resources)
    pwm, patterns, audit = [], [], []
    for asset in manifest["assets"]:
        if asset["kind"] == "rbp_motif_database":
            pwm.extend(load_motif_database(Path(asset["path"]), organism="Mus_musculus"))
    baseline_count = len(pwm)
    def add(rbp, native_id, pattern, source, evidence):
        record = {"RBP": canonical_rbp(rbp), "native_RBP": rbp, "motif_id": native_id,
                  "source": source, "evidence_class": evidence, "native_pattern": str(pattern)}
        try:
            regex, width = fixed_pattern(pattern)
            if not 3 <= width <= 25:
                raise ValueError("Outside frozen 3-25 nt motif representation; preserve as bound probe")
            record.update(pattern=regex, width=width)
            patterns.append(record)
            audit.append({**record, "status": "included_native_consensus_or_binding_sequence"})
        except ValueError as error:
            audit.append({**record, "status": "not_compatible_fixed_motif", "reason": str(error)})
    # Recover genuine mouse consensus entries, not just ATtRACT entries with PWM.
    for asset in manifest["assets"]:
        if asset["kind"] != "rbp_motif_database":
            continue
        with ZipFile(asset["path"]) as archive:
            if "ATtRACT_db.txt" not in archive.namelist():
                continue
            table = pd.read_csv(archive.open("ATtRACT_db.txt"), sep="\t", dtype=str).fillna("")
            table = table[table.Organism.eq("Mus_musculus") & table.Mutated.str.casefold().eq("no")
                          & table.Gene_id.str.startswith("ENSMUSG")]
            for index, row in table.iterrows():
                add(row.Gene_name, f"ATtRACT_native_row_{index}", row.Motif, "ATtRACT native mouse consensus; "+row.Database, "direct_mouse_native_sequence")
    rmaps = pd.read_csv(resources/"rmaps_known_human_mouse.txt", sep="\t")
    for i, row in rmaps.iterrows():
        add(row.Protein_name, f"rMAPS_author_{i}", row.regularExpression,
            "rMAPS author-maintained human/mouse reference; revision b9a9dce9", "curated_human_mouse_consensus")
    # Native RBPmap labels encode published IUPAC consensuses, not numeric PSSM.
    # Preserve that representation and never pretend that a logo filename is PWM.
    page = (resources/"RBPmap_native_motif_list.html").read_text(encoding="utf-8-sig")
    for name in re.findall(r"images/motifs_updated/([^\"<>]+)\.png", page):
        match = re.fullmatch(r"(.+)_([acgturyswkmbdhvn]+)_human", name, flags=re.I)
        if match:
            add(match[1], name, match[2], "RBPmap native human/mouse motif label; not numeric PSSM", "curated_RBPmap_consensus_reference")
    with ZipFile(resources/"RBPDB_mouse_tables.zip") as archive:
        read = lambda part: pd.read_csv(archive.open(next(n for n in archive.namelist() if part in n)), sep="\t", header=None, dtype=str, keep_default_na=False)
        proteins, experiments, links = read("proteins"), read("experiments"), read("protExp")
        prot = proteins.set_index(0)
        experiments = experiments.drop_duplicates(0).set_index(0)
        for _, link in links.drop_duplicates([0, 1]).iterrows():
            if link[0] not in prot.index or link[1] not in experiments.index:
                raise ValueError("RBPDB native mouse relationship missing")
            row = experiments.loc[link[1]]
            if row[4] not in {"", "\\N"}:
                add(prot.loc[link[0], 4], f"RBPDB_exp_{link[1]}", row[4], f"RBPDB mouse {row[2]}; PMID={row[1]}", "direct_mouse_native_sequence")
        with ZipFile(resources/"RBPDB_mouse_matrices.zip") as matrix_zip:
            for member in matrix_zip.namelist():
                if not member.endswith(".pfm"):
                    continue
                experiment_id = Path(member).stem.split("_")[0]
                counts = np.loadtxt(io.BytesIO(matrix_zip.read(member))).T
                if counts.shape[1] != 4 or (counts.sum(axis=1) <= 0).any():
                    raise ValueError("RBPDB native PFM alphabet/orientation invalid")
                matrix = counts/counts.sum(axis=1, keepdims=True)
                for pid in links.loc[links[1].eq(experiment_id), 0].unique():
                    pwm.append(Motif(str(prot.loc[pid, 4]), "RBPDB::"+Path(member).stem,
                                     "RBPDB native mouse PFM", pwm=matrix, evidence_class="direct_experiment", target_species="Mus_musculus"))
    # Reference-only: oRNAment mixes species, so do not label this direct mouse.
    mammal_names = {canonical_rbp(m.rbp).casefold() for m in pwm} | {p["RBP"].casefold() for p in patterns}
    mapping = pd.read_csv(resources/"oRNAment_RBP_id_encoding.csv.gz", header=None, skipinitialspace=True,
                          names=["id", "native_name"])
    names = {int(r.id): re.sub(r"\s*\(\d+\)$", "", str(r.native_name)) for r in mapping.itertuples(index=False)}
    with tarfile.open(resources/"oRNAment_PWMs.tgz") as archive:
        for member in archive.getmembers():
            if not member.isfile() or not member.name.endswith(".PWM"):
                continue
            identifier = int(Path(member.name).stem)
            if identifier not in names:
                audit.append({"source": "oRNAment", "motif_id": str(identifier), "status": "no_native_RBP_mapping"})
                continue
            rbp = canonical_rbp(names[identifier])
            if rbp.casefold() not in mammal_names:
                audit.append({"source": "oRNAment", "motif_id": str(identifier), "RBP": rbp, "status": "non_mammalian_or_unresolved_label_not_transferred"})
                continue
            table = pd.read_csv(archive.extractfile(member), sep="\t")
            matrix = table[["A", "C", "G", "U"]].to_numpy(float)
            matrix /= matrix.sum(axis=1, keepdims=True)
            pwm.append(Motif(rbp, f"oRNAment::{identifier:03d}", "oRNAment native RNAcompete/RBNS PWM; species transfer unverified",
                             pwm=matrix, evidence_class="cross_species_reference_not_mouse_binding", target_species="reference_only"))
    assignments = pd.DataFrame([{"motif_token": f"p{i:05d}", "RBP": canonical_rbp(m.rbp), "native_RBP": m.rbp,
        "motif_id": m.motif_id, "source": m.source, "evidence_class": m.evidence_class,
        "target_species": m.target_species, "width": len(m.pwm), "baseline": i < baseline_count,
        "pwm_sha256": hashlib.sha256((np.asarray(m.pwm)/np.asarray(m.pwm).sum(axis=1, keepdims=True)).astype("<f8").tobytes()).hexdigest()}
        for i, m in enumerate(pwm)])
    consensus = pd.DataFrame(patterns).drop_duplicates(["RBP", "evidence_class", "pattern", "source"]).reset_index(drop=True)
    consensus["motif_token"] = [f"c{i:05d}" for i in range(len(consensus))]
    return pwm, assignments, consensus, pd.DataFrame(audit)


def average_odds_scores(motifs, catalog, windows, background, progress=None):
    """RNA-forward average odds at every valid position, no site threshold."""
    sequences = list(windows.sequence)
    width_max = max(len(m.pwm) for m in motifs)
    lookup = np.full(256, 4, dtype=np.uint8)
    for i, base in enumerate(b"ACGT"):
        lookup[base] = i
    chunks, offsets, offset = [], [], 0
    for sequence in sequences:
        offsets.append(offset)
        chunk = np.r_[lookup[np.frombuffer(sequence.encode("ascii"), dtype=np.uint8)], np.full(width_max, 4, dtype=np.uint8)]
        chunks.append(chunk)
        offset += len(chunk)
    encoded = np.concatenate(chunks)
    n = len(windows)
    junctions = windows.junction_offset.to_numpy(float)
    cache = {}
    def positions(width):
        if width not in cache:
            starts, owners = [], []
            for i, sequence in enumerate(sequences):
                locations = np.arange(max(0, len(sequence)-width+1))
                if np.isfinite(junctions[i]):
                    locations = locations[(locations < junctions[i]) & (locations+width > junctions[i])]
                starts.extend(locations+offsets[i]); owners.extend([i]*len(locations))
            starts, owners = np.asarray(starts, int), np.asarray(owners, int)
            valid = np.ones(len(starts), dtype=bool)
            for j in range(width):
                valid &= encoded[starts+j] < 4
            cache[width] = (starts[valid], owners[valid])
        return cache[width]
    scores = np.full((len(motifs), n), np.nan)
    done = {}
    for i, motif in enumerate(motifs):
        key = catalog.pwm_sha256.iloc[i]
        if key in done:
            scores[i] = scores[done[key]]
            continue
        matrix = np.asarray(motif.pwm, float)
        matrix = matrix/matrix.sum(axis=1, keepdims=True)
        matrix = (matrix*20 + .1*np.asarray(background)[None, :])/20.1
        logodds = np.log2(matrix/np.asarray(background)[None, :])
        starts, owners = positions(len(matrix))
        values = np.zeros(len(starts))
        for j in range(len(matrix)):
            values += logodds[j, encoded[starts+j]]
        counts = np.bincount(owners, minlength=n)
        total = np.bincount(owners, weights=np.exp2(values), minlength=n)
        scores[i, counts > 0] = np.log2(1+total[counts > 0]/counts[counts > 0])
        done[key] = i
        if progress and i % 100 == 0:
            progress(f"Preference-score matrices {i}/{len(motifs)}; no site P filter")
    return scores


def consensus_coverage(patterns, windows):
    scores = np.zeros((len(patterns), len(windows)))
    cache = {}
    for i, row in enumerate(patterns.itertuples(index=False)):
        if row.pattern in cache:
            scores[i] = scores[cache[row.pattern]]
            continue
        regex = re.compile("(?=("+row.pattern+"))")
        for j, window in enumerate(windows.itertuples(index=False)):
            if len(window.sequence) < row.width:
                scores[i, j] = np.nan
                continue
            covered = np.zeros(len(window.sequence), dtype=bool)
            for hit in regex.finditer(window.sequence):
                a, b = hit.start(1), hit.end(1)
                if np.isfinite(window.junction_offset) and not a < window.junction_offset < b:
                    continue
                covered[a:b] = True
            scores[i, j] = covered.mean()
        cache[row.pattern] = i
    return scores


def paired_sign_rows(scores, catalog, windows, pairs, scope, method):
    """One fixed gene-matched test, using signs rather than selected score cutoffs."""
    rows = []
    data = catalog.assign(rbp_key=catalog.RBP.str.casefold())
    for (key, evidence), group in data.groupby(["rbp_key", "evidence_class"], sort=True):
        if method == "PWM_average_odds":
            group = group.drop_duplicates("pwm_sha256")
            combined = np.nanmean(scores[group.index], axis=0)
        else:
            # Union of nucleotide coverage for multiple native patterns would
            # need position masks. Mean of UNIQUE pattern coverage is declared
            # here to avoid database duplication inflating evidence.
            group = group.drop_duplicates("pattern")
            valid_counts = np.isfinite(scores[group.index]).sum(axis=0)
            combined = np.divide(np.nansum(scores[group.index], axis=0), valid_counts,
                                 out=np.full(len(windows), np.nan), where=valid_counts > 0)
        rbp = sorted(group.RBP.unique(), key=lambda x: (x.isupper(), x))[0]
        for region in REGIONS:
            select = windows.region.eq(region)
            score = pd.Series(combined[select], index=windows.loc[select, "event_id"])
            for direction in ("inclusion_up", "inclusion_down"):
                match = pairs[pairs.status.eq("matched") & pairs.direction.eq(direction)]
                a, b = match.foreground_event.map(score), match.background_event.map(score)
                valid = np.isfinite(a) & np.isfinite(b)
                a, b = a[valid].to_numpy(float), b[valid].to_numpy(float)
                delta = a-b
                tolerance = 1e-12*np.maximum(1, np.maximum(np.abs(a), np.abs(b)))
                up, down = int((delta > tolerance).sum()), int((delta < -tolerance).sum())
                nonzero = up+down
                p = float(binomtest(up, nonzero, .5, alternative="greater").pvalue) if nonzero else 1.
                rows.append({"catalog_scope": scope, "method": method, "RBP": rbp,
                    "evidence_class": evidence, "region": region, "direction": direction,
                    "motif_assignments": len(group), "motif_tokens": ";".join(group.motif_token),
                    "matched_pairs": len(a), "foreground_higher_pairs": up, "background_higher_pairs": down,
                    "ties": len(a)-nonzero, "mean_score_difference": float(delta.mean()) if len(delta) else np.nan,
                    "median_score_difference": float(np.median(delta)) if len(delta) else np.nan,
                    "pvalue": p if len(a) >= 2 else np.nan,
                    "status": "tested" if len(a) >= 2 else "insufficient_valid_pairs"})
    result = pd.DataFrame(rows)
    result["FDR_within_declared_scope"] = bh_adjust(result.pvalue)
    return result
