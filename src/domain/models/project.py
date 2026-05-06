from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
import re


@dataclass
class ComparisonDefinition:
    comparison_id: str
    display_name: str | None = None
    biological_question: str | None = None
    rmats_name: str | None = None
    deg_name: str | None = None
    suppa_name: str | None = None
    dexseq_name: str | None = None
    dtu_name: str | None = None
    quant_name: str | None = None
    sashimi_name: str | None = None
    bam_group_name: str | None = None
    rmats_path: Path | None = None
    deg_path: Path | None = None
    quant_dir: Path | None = None
    source_experiment_group: str | None = None
    source_control_group: str | None = None
    experiment_group: str | None = None
    control_group: str | None = None
    source_direction: str = "original"
    analysis_direction: str = "original"
    reverse_direction: bool = False
    reverse_splicing: bool = False
    reverse_deg: bool = False
    expected_direction_notes: str | None = None
    known_positive_control_genes: list[str] = field(default_factory=list)
    known_negative_control_genes: list[str] = field(default_factory=list)
    output_prefix: str | None = None
    detected_rmats_modes: list[str] = field(default_factory=list)
    has_deg: bool = False
    has_quant: bool = False

    @property
    def detected_rmats_modes_label(self) -> str:
        return ", ".join(self.detected_rmats_modes) if self.detected_rmats_modes else "missing"

    @property
    def resolved_name(self) -> str:
        if self.display_name:
            return self.display_name
        if self.experiment_group and self.control_group:
            return f"{self.experiment_group} vs {self.control_group}"
        return self.comparison_id

    @property
    def display_resolved_name(self) -> str:
        if self.display_name:
            return self.display_name
        if self.experiment_group and self.control_group:
            if self.reverse_splicing or self.reverse_deg:
                return f"{self.control_group} vs {self.experiment_group}"
            return f"{self.experiment_group} vs {self.control_group}"
        return self._swap_comparison_name(self.comparison_id, self.reverse_splicing or self.reverse_deg)

    @property
    def rmats_file_label(self) -> str:
        return str(self.rmats_path) if self.rmats_path else "missing"

    @property
    def deg_file_label(self) -> str:
        return str(self.deg_path) if self.deg_path else "missing"

    def source_name(self, source: str) -> str:
        mapping = {
            "rmats": self.rmats_name,
            "deg": self.deg_name,
            "suppa": self.suppa_name,
            "dexseq": self.dexseq_name,
            "dtu": self.dtu_name,
            "quant": self.quant_name,
            "sashimi": self.sashimi_name,
            "bam": self.bam_group_name,
        }
        return mapping.get(source) or self.comparison_id

    def display_source_name(self, source: str) -> str:
        raw_name = self.source_name(source)
        if source == "rmats":
            return self._swap_comparison_name(raw_name, self.reverse_splicing)
        if source == "deg":
            return self._swap_comparison_name(raw_name, self.reverse_deg)
        return raw_name

    def source_groups(self, source: str) -> tuple[str | None, str | None]:
        return self._split_comparison_name(self.source_name(source))

    def display_source_groups(self, source: str) -> tuple[str | None, str | None]:
        return self._split_comparison_name(self.display_source_name(source))

    def inferred_source_groups(self, source: str) -> tuple[str | None, str | None]:
        left, right = self.source_groups(source)
        if left and right:
            return left, right
        if self.source_experiment_group and self.source_control_group:
            return self.source_experiment_group, self.source_control_group
        return self.experiment_group, self.control_group

    def final_groups(self, source: str) -> tuple[str | None, str | None]:
        left, right = self.inferred_source_groups(source)
        if not left or not right:
            return None, None
        should_swap = self.reverse_splicing if source == "rmats" else self.reverse_deg if source == "deg" else False
        if should_swap:
            return right, left
        return left, right

    @property
    def final_as_groups(self) -> tuple[str | None, str | None]:
        return self.final_groups("rmats")

    @property
    def final_deg_groups(self) -> tuple[str | None, str | None]:
        return self.final_groups("deg")

    @property
    def final_direction_consistent(self) -> bool:
        as_groups = self.final_as_groups
        deg_groups = self.final_deg_groups
        as_ready = all(as_groups)
        deg_ready = all(deg_groups)
        if as_ready and deg_ready:
            return as_groups == deg_groups
        return as_ready or deg_ready

    @property
    def final_direction_groups(self) -> tuple[str | None, str | None]:
        as_groups = self.final_as_groups
        deg_groups = self.final_deg_groups
        if all(as_groups) and all(deg_groups):
            return as_groups if as_groups == deg_groups else (None, None)
        if all(as_groups):
            return as_groups
        if all(deg_groups):
            return deg_groups
        return None, None

    @property
    def final_direction_label(self) -> str:
        left, right = self.final_direction_groups
        if left and right:
            return f"{left} vs {right}"
        return "AS / DEG mismatch"

    @staticmethod
    def _swap_comparison_name(name: str, should_swap: bool) -> str:
        if not should_swap or not name:
            return name
        for pattern in (r"^(?P<left>.+?)(_vs_)(?P<right>.+)$", r"^(?P<left>.+?)(-vs-)(?P<right>.+)$", r"^(?P<left>.+?)(\s+vs\s+)(?P<right>.+)$", r"^(?P<left>.+?)(-)(?P<right>.+)$"):
            match = re.match(pattern, name, flags=re.IGNORECASE)
            if match:
                left = match.group("left")
                right = match.group("right")
                delimiter = match.group(2)
                return f"{right}{delimiter}{left}"
        return name

    @staticmethod
    def _split_comparison_name(name: str) -> tuple[str | None, str | None]:
        if not name:
            return None, None
        for pattern in (
            r"^(?P<left>.+?)(_vs_)(?P<right>.+)$",
            r"^(?P<left>.+?)(-vs-)(?P<right>.+)$",
            r"^(?P<left>.+?)(\s+vs\s+)(?P<right>.+)$",
            r"^(?P<left>.+?)(-)(?P<right>.+)$",
        ):
            match = re.match(pattern, name, flags=re.IGNORECASE)
            if match:
                return match.group("left"), match.group("right")
        return None, None


@dataclass
class IsoformSampleDefinition:
    sample_id: str
    quant_sf: Path
    sample_group: str | None = None
    comparison_id: str | None = None
    condition: str | None = None
    batch: str | None = None
    include: bool = True


@dataclass
class ComparisonPairDefinition:
    pair_id: str
    comparison_a: str | None = None
    comparison_b: str | None = None
    display_name: str | None = None
    experiment_group: str | None = None
    control_group: str | None = None
    enabled: bool = True

    @property
    def resolved_name(self) -> str:
        if self.display_name:
            return self.display_name
        if self.experiment_group and self.control_group:
            return f"{self.experiment_group} vs {self.control_group}"
        if self.comparison_a and self.comparison_b:
            return f"{self.comparison_a} vs {self.comparison_b}"
        return self.pair_id


@dataclass
class VisualizationGroupDefinition:
    group_id: str
    display_name: str | None = None
    comparison_ids: list[str] = field(default_factory=list)
    enabled: bool = True

    @property
    def resolved_name(self) -> str:
        if self.display_name:
            return self.display_name
        if self.comparison_ids:
            return " + ".join(self.comparison_ids)
        return self.group_id


@dataclass
class ProjectConfig:
    project_name: str
    project_root: Path
    config_path: Path | None = None
    output_root: Path | None = None
    rmats_root: Path | None = None
    deg_root: Path | None = None
    suppa_root: Path | None = None
    dexseq_root: Path | None = None
    dtu_root: Path | None = None
    quant_root: Path | None = None
    counts_path: Path | None = None
    contrastsheet_path: Path | None = None
    selected_rmats_mode: str = "JC"
    available_comparisons: list[ComparisonDefinition] = field(default_factory=list)
    scan_messages: list[str] = field(default_factory=list)
    tool_paths: dict[str, str] = field(default_factory=dict)
    input_paths: dict[str, str] = field(default_factory=dict)
    selected_comparison_ids: list[str] = field(default_factory=list)
    comparison_order_ids: list[str] = field(default_factory=list)
    program_comparison_ids: list[str] = field(default_factory=list)
    selected_analysis_modules: list[str] = field(default_factory=list)
    comparison_pairs: list[ComparisonPairDefinition] = field(default_factory=list)
    visualization_groups: list[VisualizationGroupDefinition] = field(default_factory=list)
    isoform_samples: list[IsoformSampleDefinition] = field(default_factory=list)
    shortlist_genes: list[str] = field(default_factory=list)
    blacklist_genes: list[str] = field(default_factory=list)
    candidate_selection_genes: dict[str, list[str]] = field(default_factory=dict)
    candidate_blacklist_genes: dict[str, list[str]] = field(default_factory=dict)
    removed_comparison_ids: list[str] = field(default_factory=list)
    confirmed: bool = False
    pairing_confirmed: bool = False
    comparison_sets_confirmed: bool = False
    visualization_groups_confirmed: bool = False
