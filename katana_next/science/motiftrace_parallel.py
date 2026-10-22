"""Read-only adapter for completed, independently calibrated MAPP and MotifTrace.

This module deliberately has no fitting, scanning, quantification or subprocess API.
Population statistics and descriptive transcript loci have different units of evidence.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path, PurePosixPath

import pandas as pd

ASSET_KIND = "motiftrace_parallel_results"
SCHEMA = "motiftrace.completed_parallel.v1"
KEY = ["model_id", "site", "start", "end"]
MT_P = "MotifTrace_conditional_family_P"
MAPP_P = "MAPP_empirical_P_if_estimable"


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def records(frame):
    # Do not use pandas' fixed-decimal JSON: very small original P values become 0.
    return frame.astype(object).where(pd.notna(frame), None).to_dict("records")


class ParallelModelBundle:
    def __init__(self, asset, *, species, genome_build):
        self.asset = asset
        self.manifest_path = Path(asset.path).resolve()
        if digest(self.manifest_path) != asset.metadata.get("sha256"):
            raise ValueError("Parallel result manifest SHA256 mismatch")
        self.root = self.manifest_path.parent
        self.plan = json.loads(self.manifest_path.read_text(encoding="utf-8"))
        p = self.plan
        if p.get("schema") != SCHEMA or p.get("state") != "complete":
            raise ValueError("A completed parallel-method result bundle is required")
        if (p.get("species"), p.get("genome_build")) != (species, genome_build):
            raise ValueError("Parallel result species/genome mismatch")
        if p.get("new_inference") is not False or set(p.get("methods", {})) != {"MotifTrace", "MAPP"}:
            raise ValueError("Both original methods must be retained without new inference")
        if not p.get("cases") or not p.get("family", {}).get("units_per_case"):
            raise ValueError("Missing complete comparison/family identities")
        self.files = {}
        for role, info in p["files"].items():
            name = PurePosixPath(info["path"])
            if name.is_absolute() or ".." in name.parts or "\\" in str(name) or ":" in str(name):
                raise ValueError("Unsafe result member path")
            target = self.root.joinpath(*name.parts)
            if target.is_symlink() or not target.resolve().is_relative_to(self.root):
                raise ValueError("Result member escapes the completed bundle")
            if not target.is_file() or target.stat().st_size != info["bytes"] or digest(target) != info["sha256"]:
                raise ValueError("Changed completed result member: " + role)
            self.files[role] = target

    def bind(self, contrast):
        case = self.plan["cases"].get(contrast.contrast_id)
        if case is None:
            raise ValueError("Comparison is outside the frozen result bundle")
        for arm in ("numerator", "denominator"):
            if set(case[arm]) != set(getattr(contrast, arm).sample_ids):
                raise ValueError("Original sample/contrast identity mismatch: " + arm)
        return contrast.contrast_id

    def table(self, role):
        path = self.files[role]
        return pd.read_parquet(path) if path.suffix == ".parquet" else pd.read_csv(path, sep="\t", low_memory=False, float_precision="round_trip")

    def payload(self, cases=None, check_cancelled=lambda: None):
        chosen = list(cases or self.plan["cases"])
        models = self.table("models")
        if models.model_id.duplicated().any() or len(models) != self.plan["family"]["models"]:
            raise ValueError("Motif model identities changed")
        motif_ids = set(models.model_id)
        model_patterns = models.set_index("model_id").literal_pattern.to_dict()
        result = dict(schema=SCHEMA, title="Katana MotifTrace / parallel MAPP", new_inference=False,
                      methods=self.plan["methods"], family=self.plan["family"], cases={},
                      limitations=self.plan["limitations"], sources=self.plan.get("sources", {}))
        all_sequences = self.table("window_sequences")
        all_hits = self.table("window_hits")
        all_availability = self.table("window_availability")
        for case in chosen:
            check_cancelled()
            identity = self.plan["cases"][case]
            windows = self.table(case + ":windows")
            branches = self.table(case + ":branches")
            usage = self.table(case + ":usage")
            if (len(windows) != self.plan["family"]["units_per_case"] or windows.duplicated(KEY).any()
                    or set(windows.model_id) != motif_ids or set(windows["case"]) != {case}):
                raise ValueError("Complete motif/window family changed")
            if not windows.pattern.eq(windows.model_id.map(model_patterns)).all():
                raise ValueError("Motif sequence/model identity mismatch")
            for col in (MT_P, MAPP_P):
                values = windows[col].dropna()
                if not values.between(0, 1).all():
                    raise ValueError("Invalid original P value")
            fitted = windows.MAPP_fitted
            if fitted.dtype != bool or windows.loc[~fitted, MAPP_P].notna().any():
                raise ValueError("MAPP non-applicability cannot be replaced by a P value")
            if usage.event_id.duplicated().any() or len(usage) != identity["event_count"]:
                raise ValueError("Original event universe changed")
            if int(usage.role.eq("background").sum()) != identity["background_count"]:
                raise ValueError("Original background changed")
            branch_key = KEY + ["channel", "score"]
            if (len(branches) != len(windows) * self.plan["family"]["MotifTrace_branches"]
                    or branches.duplicated(branch_key).any() or set(branches["case"]) != {case}):
                raise ValueError("Original MotifTrace branches changed")
            keys = set(map(tuple, windows[KEY].itertuples(index=False, name=None)))
            if set(map(tuple, branches[KEY].itertuples(index=False, name=None))) != keys:
                raise ValueError("MotifTrace and MAPP window identities differ")
            sequences = all_sequences[all_sequences["case"].eq(case)].copy()
            hits = all_hits[all_hits["case"].eq(case)].copy()
            available = all_availability[all_availability["case"].eq(case)].copy()
            events = set(usage.event_id)
            if not set(sequences.event_id) <= events or not set(available.event_id) <= events:
                raise ValueError("Sequence/usage event is outside the original family")
            if sequences.duplicated(["event_id", "pattern"]).any() or available.duplicated(["event_id", "pattern"]).any():
                raise ValueError("Ambiguous event/motif sequence join")
            seq_lookup = {(r.event_id, r.pattern): r for r in sequences.itertuples()}
            for hit in hits.itertuples():
                seq = seq_lookup.get((hit.event_id, hit.pattern_DNA))
                if seq is None or seq.sequence_RNA[hit.RNA_start0:hit.RNA_end0] != hit.matched_RNA:
                    raise ValueError("Saved motif hit does not match its actual RNA")
                expected = ((seq.window_start0 + hit.RNA_start0, seq.window_start0 + hit.RNA_end0)
                            if seq.strand == "+" else
                            (seq.window_end0 - hit.RNA_end0, seq.window_end0 - hit.RNA_start0))
                if expected != (hit.genomic_start0, hit.genomic_end0) or (seq.chr, seq.strand) != (hit.chr, hit.strand):
                    raise ValueError("Saved motif RNA/genomic coordinates differ")
            for seq in sequences.itertuples():
                if (seq.sequence_DNA.replace("T", "U") != seq.sequence_RNA
                        or len(seq.sequence_RNA) != seq.window_end0 - seq.window_start0):
                    raise ValueError("Saved window RNA/reference identity mismatch")
            data = dict(identity=identity, windows=records(windows), branches=records(branches),
                        events=records(usage), sequences=records(sequences), hits=records(hits),
                        availability=records(available),
                        calibration=json.loads(self.files[case + ":calibration"].read_text()),
                        design=records(self.table(case + ":design")))
            if case + ":rna_case" in self.files:
                data["rna_case"] = self._rna_case(case, events, motif_ids, model_patterns)
            result["cases"][case] = data
        return result

    def _rna_case(self, case, cohort_events, motif_ids, model_patterns):
        selected = json.loads(self.files[case + ":rna_case"].read_text())
        forms = self.table(case + ":rna_forms")
        transcript_samples = self.table(case + ":rna_transcript_samples")
        das = self.table(case + ":rna_das")
        hits = self.table(case + ":rna_hits")
        hits = hits[hits.model_id.isin(motif_ids)].copy()
        loci = self.table(case + ":rna_loci")
        loci = loci[loci.model_id.isin(motif_ids)].copy()
        per_sample = self.table(case + ":rna_locus_samples")
        per_sample = per_sample[per_sample.site_key.isin(set(loci.site_key))].copy()
        if selected["RBP_KD"] != case or forms.transcript_id.duplicated().any():
            raise ValueError("RNA case identity mismatch")
        if (len(das) != 1 or das.event_id.iloc[0] != selected["event_id"]
                or not forms.gene_id.eq(selected["gene_id"]).all()
                or not forms["chr"].eq(selected["chr"]).all()
                or not forms.strand.eq(selected["strand"]).all()):
            raise ValueError("Exact DAS event/transcript coordinate identity mismatch")
        expected_samples = set(self.plan["cases"][case]["numerator"] + self.plan["cases"][case]["denominator"])
        if (set(transcript_samples.sample_id) != expected_samples
                or transcript_samples.duplicated(["transcript_id", "sample_id"]).any()
                or not set(transcript_samples.transcript_id) <= set(forms.transcript_id)
                or not set(per_sample.sample_id) <= expected_samples
                or per_sample.duplicated(["site_key", "sample_id"]).any()):
            raise ValueError("RNA locus sample/transcript identity mismatch")
        if set(loci.site_key) != set(hits.site_key):
            raise ValueError("Saved RNA locus is missing its actual transcript hits")
        refs = forms.set_index("transcript_id").to_dict("index")
        for hit in hits.itertuples():
            ref = refs.get(hit.transcript_id)
            if ref is None or ref["gene_id"] != selected["gene_id"]:
                raise ValueError("RNA hit is not on the selected annotated transcript")
            if ref["sequence_RNA"][hit.RNA_start0:hit.RNA_end0] != hit.hit_RNA:
                raise ValueError("Transcript motif slice differs from the saved RNA")
            if hit.hit_RNA.replace("U", "T") != model_patterns[hit.model_id]:
                raise ValueError("Full-library model is not the same literal motif")
            blocks = json.loads(ref["blocks_json"])
            offset, mapped = 0, []
            for start, end in blocks:
                a, b = max(offset, hit.RNA_start0), min(offset + end - start, hit.RNA_end0)
                if a < b:
                    mapped.append([start + a - offset, start + b - offset] if ref["strand"] == "+"
                                  else [end - b + offset, end - a + offset])
                offset += end - start
            if (mapped != json.loads(hit.genomic_blocks_json)
                    or (ref["chr"], ref["strand"]) != (hit.chr, hit.strand)):
                raise ValueError("Saved transcript hit split coordinates differ from annotated RNA")
        for form in forms.itertuples():
            if (hashlib.sha256(form.sequence.encode("ascii")).hexdigest() != form.sequence_SHA256
                    or form.sequence.replace("T", "U") != form.sequence_RNA
                    or sum(b-a for a,b in json.loads(form.blocks_json)) != len(form.sequence)):
                raise ValueError("Saved transcript sequence SHA256 mismatch")
        present = selected["event_id"] in cohort_events
        return dict(selected=selected, in_calibrated_event_universe=present,
                    statistical_link="same_event_in_frozen_cohort" if present else "descriptive_case_outside_calibrated_event_universe_no_event_P",
                    forms=records(forms), hits=records(hits), loci=records(loci), locus_samples=records(per_sample),
                    transcript_samples=records(transcript_samples), das=records(das),
                    dtu=records(self.table(case + ":rna_dtu")), expression=records(self.table(case + ":rna_expression")))


def render_parallel(payload, output, attachments=()):
    template = Path(__file__).with_name("motiftrace_parallel.html").read_text(encoding="utf-8")
    def compact(value):
        if isinstance(value, list):
            if value and all(isinstance(r, dict) and list(r) == list(value[0]) for r in value):
                columns = list(value[0])
                return {"$columns": columns, "$rows": [[compact(r[k]) for k in columns] for r in value]}
            return [compact(v) for v in value]
        return {k: compact(v) for k, v in value.items()} if isinstance(value, dict) else value
    text = json.dumps(compact(dict(payload, attachments=list(attachments))), ensure_ascii=True, allow_nan=False, separators=(",", ":"))
    # JSON is data, not executable source. Prevent closing the enclosing script tag.
    text = text.replace("<", "\\u003c").replace(">", "\\u003e").replace("&", "\\u0026")
    Path(output).write_text(template.replace("__PARALLEL_PAYLOAD__", text), encoding="utf-8")
    return Path(output)
