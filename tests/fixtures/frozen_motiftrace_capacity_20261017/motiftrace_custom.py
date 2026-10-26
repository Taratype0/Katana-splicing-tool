"""Bounded user motif libraries and separate, explicit cached-RNA analyses."""
from __future__ import annotations

import hashlib
import json
import re
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control

from katana_next.science.rbp_short_extension import (
    PositionIndex, consensus_scores, consensus_windows, pwm_windows,
    rank_test_matrix, window_validity,
)

NAME = "Katana MotifTrace"
MAX_MODELS = 256
FIELDS = ("motif_id", "RBP", "species", "evidence", "source")
EVIDENCE = {"direct_experiment", "homology_inferred", "model_prediction", "user_defined", "unknown"}
IUPAC = dict(zip("ACGTURYSWKMBDHVN", ("A", "C", "G", "T", "T", "AG", "CT", "CG", "AT", "GT", "AC", "CGT", "AGT", "ACT", "ACG", "ACGT")))


def sha(path):
    with Path(path).open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def _file(path):
    raw = Path(path).expanduser()
    if not raw.is_absolute() or raw.is_symlink() or not raw.is_file():
        raise ValueError("Custom motif input must be an existing absolute regular local file")
    if raw.stat().st_size > 16 * 1024 ** 2:
        raise ValueError("Custom library exceeds the explicit 16MiB local input budget")
    return raw.resolve(strict=True)


def fixed_pattern(value):
    """Allow finite words/character classes only, never executable general regex."""
    word = str(value).strip().upper()
    tokens = re.findall(r"\[[ACGTU]+\]|[ACGTURYSWKMBDHVN]", word)
    if "".join(tokens) != word or not 1 <= len(tokens) <= 25:
        raise ValueError("Motifs must be 1-25nt IUPAC words or fixed A/C/G/T/U classes; no gaps or variable regex")
    output = []
    for token in tokens:
        letters = "".join(b for b in "ACGT" if b in (token[1:-1].replace("U", "T") if token.startswith("[") else IUPAC[token]))
        output.append(letters if len(letters) == 1 else "[" + letters + "]")
    return "".join(output), len(tokens)


def _meme(text):
    if not re.match(r"\s*MEME version \d+", text):
        raise ValueError("Expected a Minimal MEME version header")
    alphabet = re.search(r"^ALPHABET\s*=\s*(\S+)\s*$", text, re.M)
    if not alphabet or alphabet[1] not in {"ACGT", "ACGU"}:
        raise ValueError("Declare MEME ALPHABET= ACGT or ACGU explicitly; protein/PSSM formats are not PWM input")
    matrices = {}
    for block in re.split(r"^MOTIF\s+", text, flags=re.M)[1:]:
        lines = block.splitlines()
        token = lines[0].split()[0]
        if token in matrices or not re.fullmatch(r"[^\s=]+", token):
            raise ValueError("MEME motif IDs must be unique")
        headers = [i for i, line in enumerate(lines) if line.startswith("letter-probability matrix:")]
        if len(headers) != 1:
            raise ValueError("Each MEME motif requires exactly one letter-probability matrix")
        start = headers[0]
        width = re.search(r"\bw\s*=\s*(\d+)", lines[start])
        alength = re.search(r"\balength\s*=\s*(\d+)", lines[start])
        if not width or not 1 <= int(width[1]) <= 25 or not alength or int(alength[1]) != 4:
            raise ValueError("MEME requires explicit alength=4 and width 1-25 for this cached-RNA contract")
        matrix = np.asarray([list(map(float, row.split())) for row in lines[start+1:start+1+int(width[1])]])
        if matrix.shape != (int(width[1]), 4) or not np.isfinite(matrix).all() or (matrix < 0).any() or not np.allclose(matrix.sum(axis=1), 1, rtol=0, atol=1e-5):
            raise ValueError("MEME matrix must contain finite nonnegative probabilities with rows summing to one")
        matrices[token] = matrix / matrix.sum(axis=1, keepdims=True)
    if not matrices:
        raise ValueError("No MEME numeric motifs found")
    return matrices


def load_custom_library(path, *, species):
    path = _file(path)
    inputs = [path]
    if path.suffix.lower() in {".meme", ".txt"}:
        matrices = _meme(path.read_text(encoding="utf-8-sig"))
        metadata = _file(str(path) + ".metadata.tsv")
        inputs.append(metadata)
        table = pd.read_csv(metadata, sep="\t", dtype=str, keep_default_na=False)
        if "motif_id" not in table or set(table.motif_id) != set(matrices) or table.motif_id.duplicated().any():
            raise ValueError("MEME sidecar must map every motif ID exactly once")
    elif path.suffix.lower() == ".tsv":
        matrices = {}
        table = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
        if "motif" not in table:
            raise ValueError("Custom consensus TSV needs a motif column")
    else:
        raise ValueError("Use consensus .tsv or Minimal MEME .meme plus .metadata.tsv")
    if not set(FIELDS).issubset(table) or table.empty or len(table) > 1024:
        raise ValueError("Require 1-1024 source rows with motif_id/RBP/species/evidence/source")
    if table[list(FIELDS)].apply(lambda column: column.str.strip().eq("")).any().any() or table.motif_id.duplicated().any():
        raise ValueError("Custom source metadata cannot be missing or have duplicate motif IDs")
    if not table.species.eq(species).all():
        raise ValueError("Custom source species must match the declared project species; do not relabel cross-species evidence")
    if not table.evidence.isin(EVIDENCE).all():
        raise ValueError("Unknown custom evidence class")
    models, sources = {}, {}
    for row in table.to_dict("records"):
        if matrices:
            matrix = matrices[row["motif_id"]]
            fingerprint = hashlib.sha256(matrix.astype("<f8").tobytes()).hexdigest()
            key = "u_p_" + fingerprint[:20]
            model = {"kind": "PWM", "matrix": matrix.tolist(), "width": len(matrix)}
        else:
            pattern, width = fixed_pattern(row["motif"])
            fingerprint = hashlib.sha256(pattern.encode()).hexdigest()
            key = "u_c_" + fingerprint[:20]
            model = {"kind": "consensus", "pattern": pattern, "width": width}
        models.setdefault(key, {**model, "origin": "user_custom", "model_SHA256": fingerprint})
        sources.setdefault(key, []).append({"model_id": key, "native_id": row["motif_id"], "RBP": row["RBP"],
            "species": row["species"], "evidence": "user_declared_" + row["evidence"], "source": row["source"]})
    if len(models) > MAX_MODELS:
        raise ValueError(f"Custom job exceeds the explicit {MAX_MODELS}-model budget; never silently truncate a library")
    for token, model in models.items():
        model.update(RBPs=";".join(sorted({row["RBP"] for row in sources[token]})),
                     native_ids=";".join(sorted({row["native_id"] for row in sources[token]})),
                     sources=";".join(sorted({row["source"] for row in sources[token]})),
                     evidence=";".join(sorted({row["evidence"] for row in sources[token]})))
    return {"models": models, "sources": sources, "inputs": [{"path": str(p), "sha256": sha(p), "bytes": p.stat().st_size} for p in inputs],
            "source_rows": len(table), "unique_models": len(models), "metadata_verified_as": "user declaration, not independent experimental validation"}


def add_custom_results(payload, library, *, mode, output, check_cancelled=lambda: None):
    """Keep the historical public family immutable and custom results separate."""
    if mode not in {"preview_only", "fixed_enrichment"}:
        raise ValueError("Unknown custom motif mode")
    output = Path(output)
    output.mkdir(parents=True, exist_ok=False)
    plan = {"tool": NAME, "mode": mode, "UTC_before_custom_scoring": datetime.now(timezone.utc).isoformat(),
            "cases": payload["cases"], "inputs": library["inputs"], "species": payload["species"], "genome_build": payload["genome_build"],
            "cached_RNA_identity_sha256": hashlib.sha256(json.dumps(payload["events"], sort_keys=True).encode()).hexdigest(),
            "model_ids": sorted(library["models"]), "minimum_per_group": 20, "window_sizes": [10, 20, 50], "step": 1,
            "score_and_test": "original PositionIndex PWM/overlapping consensus and one-sided asymptotic MW with continuity; original NA rules",
            "correction": "ONE separate BH family across all custom models, supplied cases, two directions, scalar scopes and all eligible windows",
            "original_public_family_unchanged": payload["family_size"],
            "interpretation": "User-selected library follow-up, not independent confirmation or evidence of improved sensitivity",
            "selection_posthoc_possible": True, "quantification_or_DAS_refit": False}
    (output / "CUSTOM_ANALYSIS_PLAN.json").write_text(json.dumps(plan, indent=2), encoding="utf-8")
    payload["models"].update(library["models"])
    payload["sources"].update(library["sources"])
    payload["model_count"] += library["unique_models"]
    payload["source_count"] += library["source_rows"]
    payload["custom"] = {"mode": mode, "family_size": None, "summary": [], "library": library, "plan": plan}
    (output / "CUSTOM_LIBRARY.json").write_text(json.dumps(library, indent=2), encoding="utf-8")
    if mode == "preview_only":
        return payload
    tests = []
    for case in payload["cases"]:
        check_cancelled()
        rows = payload["events"][case]
        roles = np.array([row["role"] for row in rows])
        length = np.array([row["cassette_length"] for row in rows])
        forms = {region: [next(form for form in row["forms"] if form["region"] == region) for row in rows]
                 for region in ("cassette_exon", "inclusion_junctions", "skipping_junction")}
        sequences = {region: [form["RNA"].replace("U", "T") for form in values] for region, values in forms.items()}
        indices = {region: PositionIndex(words, 25, None if region == "cassette_exon" else [form["junctions"] for form in forms[region]])
                   for region, words in sequences.items()}
        windows = {}
        for width in (10, 20, 50):
            eligible = [a for a in range(max(0, int(length.max()) - width + 1))
                        if ((length >= a+width) & (roles == "bg")).sum() >= 20
                        and any(((length >= a+width) & (roles == group)).sum() >= 20 for group in ("up", "dn"))]
            span = max(eligible)+1 if eligible else 0
            windows[width] = (span, window_validity(sequences["cassette_exon"], width, span))
        if len(library["models"]) * (14 + 2 * sum(span for span, _ in windows.values())) > 1_000_000:
            raise ValueError("Custom analysis exceeds one million planned tests per contrast; preserve inputs and use an explicitly budgeted job")
        for token, model in library["models"].items():
            check_cancelled()
            current = {}
            for region, words in sequences.items():
                if model["kind"] == "PWM":
                    values, _, sites = indices[region].pwm(model["matrix"], payload["backgrounds"][case])
                else:
                    joins = None if region == "cassette_exon" else [form["junctions"] for form in forms[region]]
                    values, _, sites = consensus_scores(words, model["pattern"], joins)
                current[region] = (values, sites)
            scopes = {region: (values[0], np.ones(len(rows), bool)) for region, values in current.items()}
            for scope, subset in (("short_exon_lt50", length < 50), ("short_exon_1_9", length < 10),
                                  ("short_exon_10_19", (length >= 10) & (length < 20)), ("short_exon_20_49", (length >= 20) & (length < 50))):
                scopes[scope] = (current["cassette_exon"][0], subset)
            def record(scope, score, subset, width=0):
                score = np.asarray(score)
                for role, direction in (("up", "inclusion_up"), ("dn", "inclusion_down")):
                    p, a, b, n1, n2 = rank_test_matrix(score[subset & (roles == role)], score[subset & (roles == "bg")], minimum=20)
                    tests.append(pd.DataFrame(dict(contrast=case, model_id=token, scope=scope, direction=direction,
                        window_nt=width, window_start0=np.arange(len(p)) if width else -1,
                        raw_P=p, mean_foreground=a, mean_background=b, mean_difference=a-b, n_foreground=n1, n_background=n2)))
            for scope, (score, subset) in scopes.items():
                record(scope, score, subset)
            for width, (span, valid) in windows.items():
                if not span:
                    continue
                sites = current["cassette_exon"][1]
                scores = (pwm_windows(sites, len(rows), model["width"], width, span, valid) if model["kind"] == "PWM"
                          else consensus_windows(sites, width, span, valid, model["width"]))
                record("full_exon_sliding", scores, np.ones(len(rows), bool), width)
    table = pd.concat(tests, ignore_index=True)
    table["q_custom_library_family"] = np.nan
    valid = table.raw_P.notna()
    if valid.any():
        table.loc[valid, "q_custom_library_family"] = false_discovery_control(table.loc[valid, "raw_P"].to_numpy())
    table.to_csv(output / "CUSTOM_ALL_TESTS_WITH_NA.tsv.gz", sep="\t", index=False)
    table.to_parquet(output / "CUSTOM_ALL_TESTS_WITH_NA.parquet", index=False)
    payload["custom"]["family_size"] = int(valid.sum())
    for case in payload["cases"]:
        scalar = table[table.contrast.eq(case) & table.window_nt.eq(0)]
        for token, frame in scalar.groupby("model_id", sort=False):
            fields = ["scope", "direction", "n_foreground", "n_background", "mean_difference", "raw_P", "q_custom_library_family"]
            payload["scalar"][case][token] = frame[fields].astype(object).where(pd.notna(frame[fields]), None).values.tolist()
    for (case, scope), frame in table.groupby(["contrast", "scope"], sort=False):
        q = frame.q_custom_library_family
        payload["custom"]["summary"].append(dict(contrast=case, scope=scope, estimable=int(frame.raw_P.notna().sum()),
            significant_tests=int(q.lt(.05).sum()), minimum_q=float(q.min()) if q.notna().any() else None))
    (output / "CUSTOM_RESULTS.json").write_text(json.dumps({"state": "complete", "family_size": int(valid.sum()),
        "summary": payload["custom"]["summary"], "original_family_unchanged": True, "not_independent_validation": True}, indent=2), encoding="utf-8")
    return payload
