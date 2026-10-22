"""Calibrated FIMO sites and gene-independent matched-event motif enrichment."""
import json
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import binomtest

from katana_cosplice.motif_scanner import IndexedFasta, _event_windows, reverse_complement
from katana_next.science.statistics import bh_adjust


def measurable_event_background(frame):
    from katana_next.science.enrichment import enrichment_background
    measured, audit = enrichment_background(frame)
    if "measurement_eligible" in measured:
        eligible = measured.measurement_eligible.fillna(False).astype(bool)
        excluded = set(measured.loc[~eligible, "row_id"])
        audit.loc[audit.row_id.isin(excluded), "background_eligible"] = False
        audit.loc[audit.row_id.isin(excluded), "background_reason"] = "event_measurement_support_failed"
        measured = measured.loc[eligible].copy()
    return measured, audit


def extract_windows(events, fasta_path, window_bp, progress=None):
    fasta = IndexedFasta(fasta_path)
    with fasta.fasta.open("rb") as stream:
        return _extract_windows(events, fasta, window_bp, stream, progress)


def _extract_windows(events, fasta, window_bp, stream, progress):
    from katana_next.science.bulk import EVENT_COORDINATES
    windows, summaries = [], []
    coordinates = [c for c in EVENT_COORDINATES if c not in {"chr", "strand"} and c in events]
    for _, original in events.drop_duplicates("event_id").sort_values("event_id").iterrows():
        row = original.copy()
        row["event_coordinates"] = "|".join(str(int(row[c])) for c in coordinates if pd.notna(row[c]))
        chromosome = str(row.chr)
        if chromosome not in fasta.index:
            chromosome = chromosome.removeprefix("chr") if chromosome.startswith("chr") else "chr" + chromosome
        if chromosome not in fasta.index or row.strand not in {"+", "-"}:
            raise ValueError(f"Invalid FASTA chromosome/strand for {row.event_id}")
        sequence_parts = []
        for index, (start, end) in enumerate(_event_windows(row, window_bp)):
            sequence = fasta.fetch(chromosome, start, end, stream=stream)
            if row.strand == "-":
                sequence = reverse_complement(sequence)
            sequence_parts.append(sequence)
            windows.append({"window_id": f"w{len(windows):09d}", "event_id": row.event_id, "chr": chromosome,
                            "strand": row.strand, "start_0based": start,
                            "end_exclusive": start+len(sequence), "sequence": sequence})
        combined = "".join(sequence_parts)
        valid = sum(combined.count(b) for b in "ACGT")
        summaries.append({"event_id": row.event_id, "gene_id": row.get("gene_id", ""),
                          "event_type": row.event_type, "bases": len(combined),
                          "gc": (combined.count("G")+combined.count("C"))/valid if valid else np.nan,
                          "ambiguous_fraction": 1-valid/max(1, len(combined))})
        if progress and len(summaries) % 1000 == 0:
            progress(len(summaries), len(events))
    return pd.DataFrame(windows), pd.DataFrame(summaries)


def match_background(summary, foreground_ids, *, gc_tolerance=.05, length_tolerance=.20):
    """One lexicographically chosen event per foreground gene, unique control genes.

    Matching is deterministic and blind to PWM hits. All foreground genes are
    excluded from controls, including genes whose representative is unmatched.
    """
    summary = summary.sort_values("event_id").copy()
    fg = summary[summary.event_id.isin(foreground_ids)]
    fg_genes = set(fg.gene_id.dropna())
    controls = summary[~summary.gene_id.isin(fg_genes) & summary.gc.notna() & summary.ambiguous_fraction.le(.05)]
    used_fg, used_bg, rows = set(), set(), []
    for _, event in fg.iterrows():
        row = {"foreground_event": event.event_id, "foreground_gene": event.gene_id,
               "background_event": "", "background_gene": "", "status": "unmatched",
               "foreground_gc": event.gc, "foreground_bases": event.bases}
        if pd.isna(event.gene_id) or not str(event.gene_id):
            row["status"] = "missing_gene_id"
        elif event.gene_id in used_fg:
            row["status"] = "same_gene_additional_event_sites_only"
        else:
            used_fg.add(event.gene_id)
            pool = controls[~controls.gene_id.isin(used_bg) & controls.event_type.eq(event.event_type)].copy()
            pool["gc_distance"] = (pool.gc-event.gc).abs()
            pool["length_distance"] = (pool.bases-event.bases).abs()/max(1, event.bases)
            pool = pool[pool.gc_distance.le(gc_tolerance) & pool.length_distance.le(length_tolerance)]
            if event.ambiguous_fraction > .05 or not np.isfinite(event.gc):
                row["status"] = "ambiguous_foreground_sequence"
            elif not pool.empty:
                selected = pool.sort_values(["gc_distance", "length_distance", "event_id"]).iloc[0]
                used_bg.add(selected.gene_id)
                row.update(background_event=selected.event_id, background_gene=selected.gene_id, status="matched",
                           background_gc=selected.gc, background_bases=selected.bases,
                           gc_distance=selected.gc_distance, length_distance=selected.length_distance)
        rows.append(row)
    return pd.DataFrame(rows)


def write_meme(motifs, path, background):
    lines = ["MEME version 4", "", "ALPHABET= ACGT", "", "strands: +", "",
             "Background letter frequencies (from matched control windows):", " ".join(f"{base} {freq:.9g}" for base, freq in zip("ACGT", background)), ""]
    catalog = []
    for index, motif in enumerate(motifs):
        if motif.pwm is None:
            raise ValueError("Calibrated mode requires a numeric PWM, not a consensus string")
        pwm = np.asarray(motif.pwm, float)
        if pwm.ndim != 2 or pwm.shape[1] != 4 or not np.isfinite(pwm).all() or (pwm < 0).any() or (pwm.sum(axis=1) <= 0).any():
            raise ValueError(f"Invalid PWM: {motif.motif_id}")
        pwm = pwm/pwm.sum(axis=1, keepdims=True)
        pwm = (pwm*20 + .1*np.asarray(background)[None, :])/20.1
        token = f"m{index:06d}"
        catalog.append({"motif_token": token, "motif_id": motif.motif_id, "RBP": motif.rbp, "source": motif.source,
                        "evidence_class": motif.evidence_class, "target_species": motif.target_species})
        lines.extend([f"MOTIF {token}", f"letter-probability matrix: alength= 4 w= {len(pwm)} nsites= 20 E= 0", *[" ".join(f"{x:.10g}" for x in row) for row in pwm], ""])
    path.write_text("\n".join(lines), encoding="ascii")
    return pd.DataFrame(catalog)


def _wsl_path(path):
    path = Path(path).resolve().as_posix()
    if len(path) < 3 or path[1] != ":":
        raise ValueError("WSL FIMO requires files on an explicitly mapped Windows drive")
    return f"/mnt/{path[0].lower()}/{path[3:]}"


def calibrated_scan(motifs, windows, directory, runtime_config, threshold=1e-4):
    directory.mkdir(parents=True, exist_ok=True)
    controls = windows[windows.role.eq("background")]
    sequence = "".join(controls.sequence)
    frequencies = np.array([sequence.count(base)+1 for base in "ACGT"], float)
    frequencies /= frequencies.sum()
    meme, data, output = directory/"motifs.meme", directory/"windows.json", directory/"fimo_result.json"
    catalog = write_meme(motifs, meme, frequencies)
    data.write_text(windows.to_json(orient="records"), encoding="utf-8")
    worker = Path(__file__).resolve().parents[2]/"scripts"/"fimo_worker.py"
    config = json.loads(Path(runtime_config).read_text(encoding="utf-8"))
    if config.get("wsl_distribution"):
        command = ["wsl", "-d", config["wsl_distribution"], "--", config["python"], "-B", _wsl_path(worker),
                   _wsl_path(meme), _wsl_path(data), _wsl_path(output)]
    else:
        command = [config["python"], "-B", str(worker), str(meme), str(data), str(output)]
    command += ["--threshold", str(threshold)]
    completed = subprocess.run(command, capture_output=True, text=True, timeout=1800)
    (directory/"fimo_command.json").write_text(json.dumps({"command": command, "runtime": config, "background_ACGT": frequencies.tolist()}, indent=2), encoding="utf-8")
    (directory/"fimo.log").write_text(completed.stdout+completed.stderr, encoding="utf-8")
    if completed.returncode:
        raise RuntimeError("FIMO worker failed: " + completed.stderr[-2000:])
    payload = json.loads(output.read_text(encoding="utf-8"))
    sites = pd.DataFrame(payload.pop("sites"), columns=["motif_token", "window_id", "start_1based", "stop_1based", "score", "site_pvalue"])
    sites = sites.merge(catalog, on="motif_token", how="left", validate="many_to_one").merge(windows.drop(columns="sequence"), on="window_id", how="left", validate="many_to_one")
    sites["site_qvalue_global_conservative"] = np.nan
    if len(sites):
        order = np.argsort(sites.site_pvalue.to_numpy())
        q = sites.site_pvalue.to_numpy()[order] * payload["tested_positions_all_motifs"]/np.arange(1, len(sites)+1)
        sites.loc[sites.index[order], "site_qvalue_global_conservative"] = np.minimum(1, np.minimum.accumulate(q[::-1])[::-1])
        plus = sites.strand.eq("+")
        sites["genomic_start_0based"] = np.where(plus, sites.start_0based+sites.start_1based-1, sites.end_exclusive-sites.stop_1based)
        sites["genomic_end_exclusive"] = np.where(plus, sites.start_0based+sites.stop_1based, sites.end_exclusive-sites.start_1based+1)
    payload["site_q_rule"] = "conservative global BH across every motif-position test, omitted p-values replaced with 1; ambiguous sites counted conservatively"
    return sites, catalog, payload


def matched_enrichment(sites, catalog, pairs):
    matched = pairs[pairs.status.eq("matched")]
    rows = []
    catalog = catalog.copy()
    if "evidence_class" not in catalog:
        catalog["evidence_class"] = "unspecified"
    catalog["RBP_casefold"] = catalog.RBP.astype(str).str.casefold()
    groups = catalog.groupby(["RBP_casefold", "evidence_class"], sort=True)
    for (rbp, evidence_class), group in groups:
        native_names = sorted(set(group.RBP))
        rbp = sorted(native_names, key=lambda value: (value.isupper(), value))[0]
        hit_ids = set(sites.loc[sites.motif_token.isin(group.motif_token), "event_id"]) if "motif_token" in sites else set(sites.loc[sites.RBP.eq(rbp), "event_id"])
        a = matched.foreground_event.isin(hit_ids).to_numpy()
        b = matched.background_event.isin(hit_ids).to_numpy()
        only_a, only_b = int(np.sum(a & ~b)), int(np.sum(b & ~a))
        discordant = only_a+only_b
        p = float(binomtest(only_a, discordant, .5, alternative="greater").pvalue) if discordant else 1.
        if len(matched) < 2:
            p = np.nan
        rows.append({"RBP": rbp, "native_RBP_names": ";".join(native_names), "evidence_class": evidence_class, "motifs_tested": len(group), "matched_gene_pairs": len(matched), "foreground_hit_genes": int(a.sum()),
                     "background_hit_genes": int(b.sum()), "foreground_only_pairs": only_a, "background_only_pairs": only_b,
                     "hit_fraction_difference": float(a.mean()-b.mean()) if len(a) else np.nan,
                     "enrichment_pvalue": p, "test": "one-sided exact paired McNemar/binomial; one event per gene",
                     "status": "tested" if len(matched) >= 2 else "insufficient_matched_genes"})
    result = pd.DataFrame(rows)
    result["enrichment_fdr"] = bh_adjust(result.enrichment_pvalue)
    return result
