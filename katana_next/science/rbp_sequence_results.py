"""Verified reuse of completed RNA-motif results; never fits or scores a cohort."""
from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path, PurePosixPath

import numpy as np
import pandas as pd


ASSET_KIND = "rbp_sequence_results"
SCALAR_Q = "q_primary_ALL_models_contrasts_scopes_windows"
BASE_MEMBERS = (
    "STATUS.json", "RESULTS.json", "ANALYSIS_PLAN_AND_FREEZE.json",
    "UNIQUE_TESTED_MODELS.tsv", "ALL_NATIVE_ASSIGNMENTS_AND_DEDUPLICATION.tsv",
    "NATIVE_NUMERIC_MATRICES.npz", "NATIVE_CONSENSUS_PATTERNS.json",
    "DELIVERY/ALL_COMPLETE_SEQUENCE_AND_SHORT_STRATUM_TESTS.tsv.gz",
)


def digest(path):
    result = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            result.update(block)
    return result.hexdigest()


def _relative(value):
    if not isinstance(value, str) or "\\" in value or ":" in value:
        raise ValueError("Bundle members must use safe relative POSIX paths")
    path = PurePosixPath(value)
    if path.is_absolute() or not path.parts or any(part in {"..", "."} for part in path.parts):
        raise ValueError("Bundle member escapes its registered root")
    return path


class CompletedRbpBundle:
    def __init__(self, asset, *, species, genome_build):
        self.asset = asset
        metadata = asset.metadata
        if not species or metadata.get("species") != species:
            raise ValueError("RBP result species must match the project explicitly")
        if not genome_build or metadata.get("genome_build") != genome_build:
            raise ValueError("RBP result genome build differs from the project")
        if not metadata.get("sha256") or digest(asset.path) != metadata["sha256"]:
            raise ValueError("Completed RBP manifest SHA256 mismatch")
        self.root = Path(metadata["bundle_root"]).resolve(strict=True)
        if not asset.path.resolve(strict=True).is_relative_to(self.root):
            raise ValueError("RBP manifest must be inside its declared bundle root")
        members = json.loads(asset.path.read_text(encoding="utf-8"))
        self.members = {}
        for item in members:
            key = str(_relative(item["path"]))
            if key in self.members or not re.fullmatch(r"[0-9a-f]{64}", item.get("sha256", "")):
                raise ValueError("Duplicate or invalid RBP member integrity record")
            self.members[key] = item
        self.verified = {}
        self.results = self.read_json("RESULTS.json")
        if self.read_json("STATUS.json").get("state") != "complete" or self.results.get("state") != "complete":
            raise ValueError("RBP analysis is not actually complete; partial results stay separate")
        self.plan = self.read_json("ANALYSIS_PLAN_AND_FREEZE.json")
        self.cases = self.plan["foreground_background"]
        if not self.cases or any(re.search(r"[/\\:]", case) for case in self.cases):
            raise ValueError("Invalid completed contrast identity")

    def path(self, relative, *, verify=True):
        key = str(_relative(relative))
        item = self.members.get(key)
        if item is None:
            raise ValueError("RBP member is not in the frozen manifest: " + key)
        path = self.root.joinpath(*PurePosixPath(key).parts).resolve(strict=True)
        if not path.is_relative_to(self.root) or not path.is_file():
            raise ValueError("RBP member is not a regular file inside the registered root")
        if verify and key not in self.verified:
            if path.stat().st_size != item["bytes"] or digest(path) != item["sha256"]:
                raise ValueError("RBP member SHA256/size mismatch: " + key)
            self.verified[key] = item["sha256"]
        return path

    def read_json(self, relative):
        return json.loads(self.path(relative).read_text(encoding="utf-8"))

    def input_paths(self):
        names = list(BASE_MEMBERS)
        for case in self.cases:
            names.extend(f"{case}/{name}" for name in ("ALL_FROZEN_EVENTS.parquet", "ACTUAL_AFFECTED_SEQUENCES.parquet"))
        names.extend(name for name in self.members if name.startswith("DELIVERY/figures/") and name.endswith(".png"))
        names.extend(["ACTUAL_AFFECTED_RNA_MOTIF_HITS.tsv.gz", "ALL_TESTS_WITH_NA_DENOMINATORS.tsv.gz"])
        return [self.asset.path, *[self.path(name, verify=False) for name in names]]

    def bind(self, complete, selected, *, view_scope):
        if view_scope not in {"selected_comparison", "all_frozen_comparisons"}:
            raise ValueError("Unknown display-only comparison scope")
        if selected.empty or selected.comparison_id.nunique() != 1 or not selected.event_type.eq("SE").all():
            raise ValueError("Completed SE maps require one exact SE-only DAS selection")
        case = str(selected.comparison_id.iloc[0])
        if case not in self.cases:
            raise ValueError("The saved contrast is not in this completed RBP bundle")
        expected = set(self.cases[case]["up"]) | set(self.cases[case]["dn"])
        if set(selected.event_id.astype(str)) != expected:
            raise ValueError("Use the full frozen DAS selection; event filters in the viewer do not refit enrichment")
        cases = list(self.cases) if view_scope == "all_frozen_comparisons" else [case]
        for current in cases:
            frozen = pd.read_parquet(self.path(f"{current}/ALL_FROZEN_EVENTS.parquet"))
            actual = complete.loc[complete.comparison_id.eq(current)].set_index("event_id", verify_integrity=True)
            if not set(frozen.event_id).issubset(actual.index):
                raise ValueError("Frozen RBP foreground/background rows are missing from the project")
            compare = actual.loc[frozen.event_id]
            if not compare.event_type.eq("SE").all():
                raise ValueError("SE result transport changed the event type")
            for name in ("chr", "strand", "gene_id"):
                if not np.array_equal(compare[name].astype(str).to_numpy(), frozen[name].astype(str).to_numpy()):
                    raise ValueError("Frozen event identity differs in " + name)
            for name in ("exonStart_0base", "exonEnd", "upstreamES", "upstreamEE", "downstreamES", "downstreamEE", "dpsi", "fdr"):
                if not np.allclose(compare[name].to_numpy(float), frozen[name].to_numpy(float), atol=1e-12, rtol=0, equal_nan=True):
                    raise ValueError("Frozen coordinates or DAS effects differ in " + name)
            if not compare.measurement_eligible.fillna(False).all():
                raise ValueError("Frozen RBP events lost original measurement eligibility")
            bg = compare.loc[frozen.role.eq("bg").to_numpy()]
            if not (bg.fdr.gt(.1) & bg.dpsi.abs().lt(.05)).all():
                raise ValueError("The frozen non-DAS background no longer meets its contract")
            fg = compare.loc[frozen.role.ne("bg").to_numpy()]
            if not fg.stable_significant.fillna(False).all():
                raise ValueError("Frozen foreground lost original supported DAS status")
        return cases

    def payload(self, cases):
        models = pd.read_csv(self.path("UNIQUE_TESTED_MODELS.tsv"), sep="\t").fillna("")
        sources = pd.read_csv(self.path("ALL_NATIVE_ASSIGNMENTS_AND_DEDUPLICATION.tsv"), sep="\t").fillna("")
        patterns = self.read_json("NATIVE_CONSENSUS_PATTERNS.json")
        required = {"model_id", "kind", "width", "RBPs"}
        if not required.issubset(models) or models.model_id.duplicated().any():
            raise ValueError("Invalid or duplicated native model catalog")
        if len(sources) != self.results["native_assignments"] or set(sources.model_id) != set(models.model_id):
            raise ValueError("Native source/model mappings are incomplete")
        output = {}
        with np.load(self.path("NATIVE_NUMERIC_MATRICES.npz"), allow_pickle=False) as matrices:
            for row in models.to_dict("records"):
                token = row["model_id"]
                if row["kind"] == "PWM":
                    matrix = matrices[token]
                    if matrix.shape != (int(row["width"]), 4) or not np.isfinite(matrix).all() or (matrix < 0).any() or not np.allclose(matrix.sum(axis=1), 1):
                        raise ValueError("Invalid original numeric matrix")
                    if token != "p_" + hashlib.sha256(matrix.astype("<f8").tobytes()).hexdigest()[:20]:
                        raise ValueError("Native numeric model fingerprint changed")
                    row["matrix"] = matrix.tolist()
                elif row["kind"] == "consensus":
                    pattern = patterns[token]
                    if not re.fullmatch(r"(?:[ACGT]|\[[ACGT]+\])+", pattern) or len(re.findall(r"\[[ACGT]+\]|[ACGT]", pattern)) != int(row["width"]):
                        raise ValueError("Native pattern is not a fixed-length DNA-encoded motif")
                    if token != "c_" + hashlib.sha256(pattern.encode()).hexdigest()[:20]:
                        raise ValueError("Native consensus model fingerprint changed")
                    row["pattern"] = pattern
                else:
                    raise ValueError("Unsupported native motif type")
                output[token] = row
        if len(patterns) != self.results["unique_consensus_patterns"]:
            raise ValueError("Consensus model denominator changed")
        source_groups = {token: group.to_dict("records") for token, group in sources.groupby("model_id", sort=False)}
        events, backgrounds, inventories = {}, {}, []
        for case in cases:
            frame = pd.read_parquet(self.path(f"{case}/ACTUAL_AFFECTED_SEQUENCES.parquet"))
            expected_ids = {value for values in self.cases[case].values() for value in values}
            if set(frame.event_id) != expected_ids or frame.duplicated(["event_id", "region"]).any():
                raise ValueError("RNA forms are missing or duplicated")
            rows = []
            for token, group in frame.groupby("event_id", sort=False):
                if set(group.region) != {"cassette_exon", "inclusion_junctions", "skipping_junction"}:
                    raise ValueError("Missing cassette or genuine splice-join sequence")
                first = group.iloc[0]
                forms = []
                for item in group.to_dict("records"):
                    if item["sequence_RNA"] != item["sequence"].replace("T", "U") or hashlib.sha256(item["sequence"].encode()).hexdigest() != item["sequence_sha256"]:
                        raise ValueError("Original reference RNA sequence hash changed")
                    blocks, joins = json.loads(item["blocks_json"]), json.loads(item["junctions_json"])
                    if sum(b - a for a, b in blocks) != len(item["sequence"]) or any(a < 0 or b <= a for a, b in blocks):
                        raise ValueError("RNA sequence and genomic blocks disagree")
                    if any(not 0 < join < len(item["sequence"]) for join in joins):
                        raise ValueError("Splice junction is outside its actual local RNA")
                    forms.append({"region": item["region"], "RNA": item["sequence_RNA"], "blocks": blocks, "junctions": joins})
                rows.append({"event_id": token, "gene_id": first.gene_id, "gene_symbol": first.gene_symbol,
                             "role": first.role, "dpsi": first.dpsi, "DAS_FDR": first.DAS_FDR,
                             "chr": first.chr, "strand": first.strand, "cassette_length": int(first.cassette_length), "forms": forms})
            events[case] = rows
            control = frame.loc[frame.region.eq("cassette_exon") & frame.role.eq("bg"), "sequence"]
            counts = np.array([sum(word.count(base) for word in control) + 1 for base in "ACGT"], float)
            backgrounds[case] = (counts / counts.sum()).tolist()
            inventories.append({"contrast": case, "events": len(rows), "foreground": sum(row["role"] != "bg" for row in rows),
                                "short_foreground_lt50": sum(row["role"] != "bg" and row["cassette_length"] < 50 for row in rows)})
        scalar = pd.read_csv(self.path("DELIVERY/ALL_COMPLETE_SEQUENCE_AND_SHORT_STRATUM_TESTS.tsv.gz"), sep="\t")
        scalar = scalar[scalar.contrast.isin(cases)]
        if set(scalar.model_id) - set(output) or scalar.duplicated(["contrast", "model_id", "scope", "direction"]).any():
            raise ValueError("Scalar results do not bind unique native hypotheses")
        for name in ("raw_P", SCALAR_Q):
            values = scalar[name].dropna()
            if not values.between(0, 1).all():
                raise ValueError("Invalid completed probability/FDR")
        if not scalar.raw_P.isna().equals(scalar[SCALAR_Q].isna()):
            raise ValueError("Unavailable tests must retain unavailable FDR")
        minimum = self.plan.get("minimum_valid_events_per_group", 20)
        insufficient = scalar.n_foreground.lt(minimum) | scalar.n_background.lt(minimum)
        if scalar.loc[insufficient, "raw_P"].notna().any():
            raise ValueError("Insufficient event denominators cannot have an inferential P value")
        fields = ["scope", "direction", "n_foreground", "n_background", "mean_difference", "raw_P", SCALAR_Q]
        scalar_data = {case: {model: group[fields].astype(object).where(pd.notna(group[fields]), None).values.tolist()
                             for model, group in scalar[scalar.contrast.eq(case)].groupby("model_id", sort=False)} for case in cases}
        return {"cases": cases, "events": events, "backgrounds": backgrounds, "models": output, "sources": source_groups,
                "scalar": scalar_data, "inventories": inventories, "summary": self.results["summary_by_contrast_scope"],
                "model_count": len(output), "source_count": len(sources), "frozen_cases": list(self.cases),
                "family_size": self.results["estimable_tests_in_one_BH_family"], "plan": {
                    name: self.plan.get(name) for name in ("PWM_score", "consensus_score", "statistic", "correction", "window_sizes", "minimum_valid_events_per_group", "short_exon_policy")},
                "genome_build": self.asset.metadata["genome_build"], "species": self.asset.metadata["species"]}


def render_browser(payload, destination, *, attachments=()):
    template = Path(__file__).with_name("rbp_sequence_browser.html").read_text(encoding="utf-8")
    data = {**payload, "attachments": list(attachments)}
    serialized = json.dumps(data, ensure_ascii=True, allow_nan=False, separators=(",", ":")).replace("<", "\\u003c").replace("\u2028", "\\u2028").replace("\u2029", "\\u2029")
    if len(serialized.encode()) > 64 * 1024 * 1024:
        raise ValueError("Sequence browser exceeds its 64MiB bounded presentation contract")
    Path(destination).write_text(template.replace("__RBP_DATA__", serialized), encoding="utf-8")
    return Path(destination)
