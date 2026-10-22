"""Read-only, identity-bound display of completed count and coverage windows."""
import base64
import json
from pathlib import Path, PurePosixPath

import numpy as np
import pandas as pd

from katana_next.science.rbp_sequence_results import digest

ASSET_KIND = "motiftrace_matched_windows"
METHODS = ("author_rMAPS3_counts", "Katana_overlapping_coverage")
Q = "q_primary_per_arm_both_contrasts_all_positions"
FIELDS = ("raw_P", Q, "n_foreground", "n_background", "mean_foreground", "mean_background", "mean_difference")
KEYS = ("region_index", "native_position", "direction", "motif_index")
COORDS = ("exonStart_0base", "exonEnd", "upstreamES", "upstreamEE", "downstreamES", "downstreamEE")


def support_status(evidence, method_q):
    """Experimental support never replaces the sequence-test threshold."""
    if evidence is None:
        return "not_assessed"
    if evidence["dataset"] is None:
        return "no_qualified_reference"
    if evidence["level"] == "ambiguous_family":
        return "family_reference_only"
    if not evidence["motif_overlap"]:
        return "no_reference_motif_overlap"
    if method_q is None:
        return "reference_overlap_method_NA"
    if method_q < .05:
        return "reference_overlap_method_passed"
    if evidence["level"] == "primary_author_significant" and evidence["external_q"] is not None and evidence["external_q"] < .05:
        return "reference_enriched_method_not_passed"
    return "reference_overlap_method_not_passed"


class MatchedWindowBundle:
    def __init__(self, asset, *, species, genome_build, source_manifest_sha256):
        meta = asset.metadata
        self.manifest = Path(asset.path).resolve()
        self.root = Path(meta["bundle_root"]).resolve()
        if not self.manifest.is_relative_to(self.root) or digest(self.manifest) != meta.get("sha256"):
            raise ValueError("Matched-window manifest identity mismatch")
        self.registration = json.loads(self.manifest.read_text(encoding="utf-8"))
        reg = self.registration
        if reg.get("schema") != "motiftrace.matched_windows.v1" or reg.get("state") != "complete":
            raise ValueError("Matched-window registration is not complete")
        if reg.get("species") != species or reg.get("genome_build") != genome_build:
            raise ValueError("Matched-window species/build mismatch")
        if reg.get("source_bundle_manifest_sha256") != source_manifest_sha256:
            raise ValueError("Matched windows belong to another RNA bundle")
        if reg.get("threshold") != .05 or set(reg.get("families", {})) != set(METHODS):
            raise ValueError("Original method-specific threshold/families required")
        self.verified = {}

    def path(self, name):
        member = self.registration["members"][name]
        rel = PurePosixPath(member["path"])
        if rel.is_absolute() or ".." in rel.parts or "\\" in member["path"] or ":" in member["path"]:
            raise ValueError("Unsafe matched-window member")
        path = (self.root / Path(*rel.parts)).resolve()
        if not path.is_relative_to(self.root):
            raise ValueError("Matched-window member escapes registered root")
        if name not in self.verified:
            if path.stat().st_size != member["bytes"] or digest(path) != member["sha256"]:
                raise ValueError("Matched-window member hash mismatch: " + name)
            self.verified[name] = member["sha256"]
        return path

    def input_paths(self):
        return [self.manifest, *(self.path(name) for name in self.registration["members"])]

    def bind(self, payload):
        reg = self.registration
        if set(payload["cases"]) - set(reg["cases"]):
            raise ValueError("Unregistered matched-window comparison")
        for case in payload["cases"]:
            actual = {e["event_id"]: e for e in payload["events"][case]}
            expected = {}
            for role, name in reg["cases"][case].items():
                rows = pd.read_parquet(self.path(name))
                if rows.event_id.duplicated().any() or not rows.comparison_id.eq(case).all():
                    raise ValueError("Invalid original matched event identity")
                for row in rows.to_dict("records"):
                    token = str(row["event_id"])
                    if token in expected or token not in actual:
                        raise ValueError("Matched foreground/background event set differs")
                    event = actual[token]
                    if event["role"] != role or any(str(row[k]) != str(event[k]) for k in ("chr", "strand", "gene_id")):
                        raise ValueError("Matched event role/identity differs")
                    if not np.allclose([row["dpsi"], row["fdr"]], [event["dpsi"], event["DAS_FDR"]], atol=1e-12, rtol=0):
                        raise ValueError("Matched event effect/FDR differs")
                    # The exact event ID contains every splice coordinate; verify against the snapshot too.
                    for key in COORDS:
                        if f"{key}={int(row[key])}" not in token:
                            raise ValueError("Matched event coordinates differ")
                    expected[token] = row
            if set(actual) != set(expected):
                raise ValueError("Matched-window background or foreground is incomplete")

    def payload(self, base, check_cancelled=lambda: None):
        self.bind(base)
        reg = self.registration
        plan = json.loads(self.path("plan").read_text(encoding="utf-8"))
        if plan["threshold"] != .05 or plan["author_source_modified"] is not False:
            raise ValueError("Original author results or threshold changed")
        params = plan["parameters"]
        if any(params[k] != v for k, v in {"intron": 200, "exon": 50, "window": 50, "step": 1, "stat_method": "mannwhitney_greater"}.items()):
            raise ValueError("Unsupported original sliding-window design")
        frame = pd.read_parquet(self.path("tests"))
        frame = frame[frame.arm.isin(METHODS)].copy()
        if set(frame.contrast.astype(str)) != set(reg["cases"]):
            raise ValueError("Incomplete original multi-contrast test family")
        for method in METHODS:
            arm = frame[frame.arm.eq(method)]
            if arm.duplicated(["contrast", *KEYS]).any():
                raise ValueError("Duplicate original window test")
            if int(arm.primary_family_P.notna().sum()) != reg["families"][method]:
                raise ValueError("Original per-method FDR family size differs")
        for field in ("raw_P", Q):
            values = frame[field].dropna()
            if not values.between(0, 1).all():
                raise ValueError("Invalid original probability")
        design_cols = ["region_index", "region", "native_position", "window_start_RNA_relative", "window_end_RNA_relative"]
        design = frame[design_cols].drop_duplicates().sort_values(["region_index", "native_position"])
        if design.duplicated(["region_index", "native_position"]).any() or not (design.window_end_RNA_relative - design.window_start_RNA_relative).eq(50).all():
            raise ValueError("Inconsistent physical window coordinates")
        design_keys = pd.MultiIndex.from_frame(design[["region_index", "native_position"]])
        catalog = frame[["motif_index", "RBP", "native_RBP", "pattern"]].drop_duplicates().sort_values("motif_index")
        if catalog.motif_index.duplicated().any():
            raise ValueError("Conflicting original motif assignments")
        groups = frame.groupby(["contrast", "motif_index", "direction", "arm"], observed=True)
        if groups.ngroups != len(reg["cases"]) * len(catalog) * 2 * len(METHODS) or set(frame.direction.astype(str)) != {"inclusion_up", "inclusion_down"}:
            raise ValueError("Incomplete original method/direction/model grid")
        data, summaries = {}, []
        for (case, motif, direction, method), rows in groups:
            check_cancelled()
            rows = rows.set_index(["region_index", "native_position"])
            if len(rows) != len(design) or not rows.index.isin(design_keys).all():
                raise ValueError("Missing original negative/NA windows")
            if case in base["cases"]:
                array = rows.reindex(design_keys)[list(FIELDS)].to_numpy(dtype="<f8")
                data[f"{case}|{int(motif)}|{direction}|{method}"] = base64.b64encode(array.tobytes()).decode("ascii")
        for case in reg["cases"]:
            c = frame[frame.contrast.eq(case)]
            passing = [set(map(tuple, c[c.arm.eq(m) & c[Q].lt(.05)][list(KEYS)].to_numpy())) for m in METHODS]
            summaries.append(dict(contrast=case, counts=len(passing[0]), coverage=len(passing[1]),
                                  shared=len(passing[0] & passing[1]), counts_only=len(passing[0]-passing[1]), coverage_only=len(passing[1]-passing[0])))
        evidence = self._evidence(frame)
        return {"cases": base["cases"], "frozen_cases": list(reg["cases"]), "methods": list(METHODS),
                "families": reg["families"], "threshold": .05, "author_revision": plan["revision"],
                "plan": {k: plan[k] for k in ("parameters", "analysis_arms", "FDR_primary", "canonical_masks", "validity", "limitations")},
                "catalog": catalog.to_dict("records"), "design": design.values.tolist(), "fields": list(FIELDS),
                "data": data, "evidence": evidence, "sources": reg["external_sources"], "summary": summaries,
                "encoding": "base64 little-endian float64 row-major; NaN means NA, never zero", "new_inference": False}

    def _evidence(self, frame):
        evidence = {}
        master = pd.read_csv(self.path("external_windows"), sep="\t", float_precision="round_trip")
        if master.window_id.duplicated().any():
            raise ValueError("Duplicate external window evidence")
        for col in ("external_CLIP_raw_P", "external_CLIP_q_all_available_windows"):
            if not master[col].dropna().between(0, 1).all():
                raise ValueError("Invalid external descriptive probability")
        case = self.registration["external_contrast"]
        indexed = frame[frame.contrast.eq(case)].set_index([*KEYS, "arm"])
        for row in master.to_dict("records"):
            region, position, direction, motif = row["window_id"].split(":")
            key = (int(region), int(position), direction, int(motif))
            originals = []
            for method, prefix in zip(METHODS, ("author", "MotifTrace")):
                source = indexed.loc[(*key, method)]
                if source["pattern"] != row["pattern"] or source["RBP"] != row["RBP"]:
                    raise ValueError("External motif/window identity differs")
                if not np.allclose([source.raw_P, source[Q]], [row[f"original_{prefix}_P"], row[f"original_{prefix}_q"]], rtol=1e-10, atol=1e-15, equal_nan=True):
                    raise ValueError("External evidence changed original method P/FDR")
                if not np.allclose([source.window_start_RNA_relative, source.window_end_RNA_relative], [row["relative_start"], row["relative_stop"]], rtol=0, atol=0):
                    raise ValueError("External evidence window coordinates differ")
                originals.append(None if pd.isna(source[Q]) else float(source[Q]))
            dataset = None if pd.isna(row["dataset"]) else str(row["dataset"])
            source = self.registration["external_sources"].get(dataset) if dataset else None
            if dataset and (source is None or source["quality"] != row["quality"]):
                raise ValueError("External source quality identity differs")
            e = {"dataset": dataset, "level": source["level"] if source else "not_available",
                 "motif_overlap": bool(pd.notna(row["foreground_motif_CLIP"]) and row["foreground_motif_CLIP"] > 0),
                 "external_P": None if pd.isna(row["external_CLIP_raw_P"]) else float(row["external_CLIP_raw_P"]),
                 "external_q": None if pd.isna(row["external_CLIP_q_all_available_windows"]) else float(row["external_CLIP_q_all_available_windows"]),
                 "foreground_motif_CLIP": None if pd.isna(row["foreground_motif_CLIP"]) else int(row["foreground_motif_CLIP"]),
                 "foreground_CLIP": None if pd.isna(row["foreground_CLIP"]) else int(row["foreground_CLIP"]),
                 "background_CLIP": None if pd.isna(row["background_CLIP"]) else int(row["background_CLIP"]),
                 "foreground_N": None if pd.isna(row["foreground_N"]) else int(row["foreground_N"]),
                 "background_N": None if pd.isna(row["background_N"]) else int(row["background_N"])}
            e["status"] = [support_status(e, q) for q in originals]
            evidence[case + "|" + row["window_id"]] = e
        return evidence


def render_windows(data, destination, attachments=()):
    serialized = json.dumps({**data, "attachments": list(attachments)}, ensure_ascii=True, separators=(",", ":"), allow_nan=False).replace("<", "\\u003c")
    if len(serialized) > 96 * 1024 * 1024:
        raise ValueError("Portable matched-window display exceeds 96 MiB; no truncation applied")
    template = Path(__file__).with_name("motiftrace_windows.html").read_text(encoding="utf-8")
    destination.write_text(template.replace("__WINDOW_DATA__", serialized), encoding="utf-8")
    return destination
