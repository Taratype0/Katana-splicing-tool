from __future__ import annotations

import os
from pathlib import Path

import csv
import yaml
from dataclasses import dataclass

from src.domain.models.project import (
    ComparisonDefinition,
    ComparisonPairDefinition,
    IsoformSampleDefinition,
    ProjectConfig,
    VisualizationGroupDefinition,
)
from src.io.validators import ProjectValidator

EVENT_TYPES = ["SE", "RI", "A3SS", "A5SS", "MXE"]
SKIP_DIR_NAMES = {
    "__pycache__",
    ".git",
    ".next",
    ".snakemake",
    ".pytest_cache",
    "cache",
    "caches",
    "tmp",
    "temp",
    "logs",
    "log",
    "work",
    "workspace",
    "katana_output",
    "node_modules",
    "bamlist",
    "bamlists",
    "bam_list",
    "bam_lists",
    "figures",
    "figure",
    "plots",
    "plot",
    "reports",
    "report",
    "multiqc",
    "cards_fixed",
    "combined_python_analysis",
    "final_delivery",
    "delivery",
    "event_files",
}
SKIP_DIR_KEYWORDS = (
    "fastq",
    "fq",
    "bigwig",
    "bw",
    "bedgraph",
    "aligned",
    "alignment",
    "star_index",
    "salmon_index",
    "sashimi",
    "jutils",
    "card",
    "figure",
    "plot",
    "report",
)


@dataclass
class ScanIndex:
    directories: list[Path]
    files: list[Path]


class ProjectScanner:
    def __init__(self) -> None:
        self.validator = ProjectValidator()

    def scan(
        self,
        root_dir: str | Path,
        input_paths: dict[str, str] | None = None,
        force_auto_scan: bool = False,
        progress_callback=None,
    ) -> ProjectConfig:
        root = Path(root_dir).expanduser().resolve()
        saved_input_paths = {} if force_auto_scan else self._load_saved_input_paths(root)
        merged_input_paths = {**saved_input_paths, **(input_paths or {})}
        self._notify(progress_callback, "Preparing scan configuration...")
        index = self._build_index(root, merged_input_paths)
        self._notify(progress_callback, "Finding contrastsheets and input roots...")
        contrastsheets = (
            [Path(merged_input_paths["contrastsheet_path"]).expanduser().resolve()]
            if merged_input_paths.get("contrastsheet_path")
            else self._find_files(index, "contrastsheet.valid.csv")
        )
        project = ProjectConfig(
            project_name=root.name,
            project_root=root,
            config_path=root / ".katana_project.yaml",
            output_root=root / "katana_output",
            input_paths=merged_input_paths,
        )

        project.rmats_root = self._input_path(merged_input_paths, "rmats_root") or self._find_first(index, ["rmats"])
        project.deg_root = self._input_path(merged_input_paths, "deg_root") or self._find_first(index, ["differential", "deg", "deseq2"])
        project.suppa_root = self._input_path(merged_input_paths, "suppa_root") or self._find_first(index, ["suppa"])
        project.dexseq_root = self._input_path(merged_input_paths, "dexseq_root") or self._find_first(index, ["dexseq_exon", "dexseq"])
        project.dtu_root = self._input_path(merged_input_paths, "dtu_root") or self._find_first(index, ["dtu", "stageR", "dexseq_dtu"])
        project.quant_root = self._input_path(merged_input_paths, "quant_root")
        project.counts_path = self._input_path(merged_input_paths, "counts_path") or self._find_file(index, "all.normalised_counts.tsv")
        project.contrastsheet_path = (
            self._input_path(merged_input_paths, "contrastsheet_path")
            or self._select_primary_contrastsheet(contrastsheets)
        )
        project.tool_paths = self._discover_tools(index)
        annotation_tsv = self._input_path(merged_input_paths, "annotation_tsv") or self._find_annotation_tsv(index)
        if annotation_tsv is not None:
            project.tool_paths["annotation_tsv"] = str(annotation_tsv)
        self._notify(progress_callback, "Discovering comparisons...")
        project.available_comparisons = self._discover_comparisons(project)
        self._notify(progress_callback, "Applying contrastsheet metadata...")
        self._apply_contrastsheets(project, contrastsheets)
        self._notify(progress_callback, "Detecting quant.sf samples...")
        project.isoform_samples = self._discover_isoform_samples(project, index)
        self._attach_isoform_samples_to_comparisons(project)
        self._notify(progress_callback, "Merging aliases and applying saved config...")
        self._merge_alias_comparisons(project)
        self._apply_saved_config(project)
        self._reconcile_source_mappings(project)
        self._merge_alias_comparisons(project)
        self._resolve_missing_source_paths(project)
        self._notify(progress_callback, "Validating detected inputs...")
        project.scan_messages = self._build_scan_messages(project)
        project.scan_messages = self.validator.validate_project(project)
        self._notify(progress_callback, "Project scan completed.")
        return project

    def _build_index(self, root: Path, input_paths: dict[str, str] | None = None) -> ScanIndex:
        if input_paths:
            return self._build_targeted_index(root, input_paths)
        return self._walk_index(root)

    def _build_targeted_index(self, root: Path, input_paths: dict[str, str]) -> ScanIndex:
        directories: list[Path] = []
        files: list[Path] = []
        seen_dirs: set[Path] = set()
        seen_files: set[Path] = set()

        def add_path(path: Path) -> None:
            resolved = path.expanduser().resolve()
            if resolved.is_file():
                if resolved not in seen_files:
                    files.append(resolved)
                    seen_files.add(resolved)
                parent = resolved.parent
                if parent.exists() and parent not in seen_dirs:
                    directories.append(parent)
                    seen_dirs.add(parent)
                return
            if resolved.is_dir():
                subset = self._walk_index(resolved)
                for directory in subset.directories:
                    if directory not in seen_dirs:
                        directories.append(directory)
                        seen_dirs.add(directory)
                for file_path in subset.files:
                    if file_path not in seen_files:
                        files.append(file_path)
                        seen_files.add(file_path)

        for key in (
            "rmats_root",
            "deg_root",
            "suppa_root",
            "dexseq_root",
            "dtu_root",
            "quant_root",
            "counts_path",
            "contrastsheet_path",
            "annotation_tsv",
        ):
            value = input_paths.get(key)
            if value:
                add_path(Path(value))

        software_root = root / "software"
        if software_root.exists():
            add_path(software_root)

        return ScanIndex(directories=directories, files=files)

    def _walk_index(self, root: Path) -> ScanIndex:
        directories: list[Path] = []
        files: list[Path] = []

        for current_root, dirnames, filenames in os.walk(root):
            current_path = Path(current_root)
            filtered_dirnames = []
            for dirname in dirnames:
                lowered = dirname.lower()
                if lowered in SKIP_DIR_NAMES:
                    continue
                if any(keyword in lowered for keyword in SKIP_DIR_KEYWORDS):
                    continue
                filtered_dirnames.append(dirname)
            dirnames[:] = filtered_dirnames

            directories.append(current_path)
            for dirname in dirnames:
                directories.append(current_path / dirname)
            for filename in filenames:
                files.append(current_path / filename)
        return ScanIndex(directories=directories, files=files)

    def _notify(self, callback, message: str) -> None:
        if callback is not None:
            callback(message)

    def _load_saved_input_paths(self, root: Path) -> dict[str, str]:
        config_path = root / ".katana_project.yaml"
        if not config_path.exists():
            return {}
        with config_path.open("r", encoding="utf-8") as handle:
            data = yaml.safe_load(handle) or {}
        raw_paths = data.get("input_paths") or {}
        cleaned: dict[str, str] = {}
        for key, value in raw_paths.items():
            if value:
                cleaned[key] = str(Path(value).expanduser())
        return cleaned

    def _input_path(self, input_paths: dict[str, str], key: str) -> Path | None:
        value = input_paths.get(key)
        return Path(value).expanduser().resolve() if value else None

    def _find_first(self, index: ScanIndex, keywords: list[str]) -> Path | None:
        for keyword in keywords:
            lowered = keyword.lower()
            for path in index.directories:
                if lowered in path.name.lower():
                    return path
        return None

    def _find_file(self, index: ScanIndex, filename: str) -> Path | None:
        matches = [path for path in index.files if path.name == filename]
        return matches[0] if matches else None

    def _find_files(self, index: ScanIndex, filename: str) -> list[Path]:
        return [path for path in index.files if path.name == filename]

    def _find_annotation_tsv(self, index: ScanIndex) -> Path | None:
        exact = [path for path in index.files if path.name == "gencode.anno.tsv"]
        if exact:
            return exact[0]
        annotation_dir_matches = sorted(
            [
                path
                for path in index.files
                if path.suffix.lower() in {".tsv", ".csv"}
                and "anno" in path.name.lower()
                and "annotation" in str(path.parent).lower()
            ],
            key=lambda path: (len(path.parts), path.name.lower()),
        )
        return annotation_dir_matches[0] if annotation_dir_matches else None

    def _select_primary_contrastsheet(self, matches: list[Path]) -> Path | None:
        if not matches:
            return None
        preferred = sorted(
            matches,
            key=lambda path: (
                0 if "splicing_genomebam" in str(path).lower() else 1,
                len(path.parts),
            ),
        )
        return preferred[0]

    def _discover_tools(self, index: ScanIndex) -> dict[str, str]:
        tools: dict[str, str] = {}
        for name in ("Jutils", "rmats2sashimiplot", "SUPPA"):
            matches = [path for path in index.directories if path.name == name]
            if matches:
                tools[name.lower()] = str(matches[0])
        return tools

    def _discover_comparisons(self, project: ProjectConfig) -> list[ComparisonDefinition]:
        comparisons: dict[str, ComparisonDefinition] = {}

        if project.rmats_root and project.rmats_root.exists():
            for comp_dir in sorted(path for path in project.rmats_root.iterdir() if path.is_dir()):
                modes = self._detect_rmats_modes(comp_dir)
                if not modes:
                    continue
                comparisons[comp_dir.name] = ComparisonDefinition(
                    comparison_id=comp_dir.name,
                    rmats_path=comp_dir,
                    detected_rmats_modes=modes,
                )

        if project.deg_root and project.deg_root.exists():
            for result in project.deg_root.rglob("*.deseq2.results.tsv"):
                comp_id = result.stem.replace(".deseq2.results", "")
                comparison = comparisons.setdefault(comp_id, ComparisonDefinition(comparison_id=comp_id))
                comparison.deg_path = result
                comparison.has_deg = True

        return sorted(comparisons.values(), key=lambda item: item.comparison_id)

    def _apply_contrastsheets(self, project: ProjectConfig, paths: list[Path]) -> None:
        if not paths:
            return

        comparisons_by_id = {item.comparison_id: item for item in project.available_comparisons}
        for path in sorted(
            paths,
            key=lambda value: (
                0 if "splicing_genomebam" in str(value).lower() else 1,
                len(value.parts),
            ),
        ):
            with path.open("r", encoding="utf-8-sig", newline="") as handle:
                reader = csv.DictReader(handle)
                for row in reader:
                    contrast = (row.get("contrast") or "").strip()
                    if not contrast:
                        continue
                    treatment = (row.get("treatment") or "").strip() or None
                    control = (row.get("control") or "").strip() or None
                    comparison = comparisons_by_id.get(contrast)
                    if comparison is None:
                        comparison = ComparisonDefinition(comparison_id=contrast)
                        project.available_comparisons.append(comparison)
                        comparisons_by_id[contrast] = comparison
                    comparison.source_experiment_group = treatment or comparison.source_experiment_group
                    comparison.source_control_group = control or comparison.source_control_group
                    comparison.experiment_group = treatment or comparison.experiment_group
                    comparison.control_group = control or comparison.control_group
                    comparison.source_direction = "experiment_vs_control" if treatment and control else comparison.source_direction
                    comparison.analysis_direction = (
                        "control_vs_experiment"
                        if comparison.reverse_direction
                        else "experiment_vs_control"
                    ) if treatment and control else comparison.analysis_direction
                    comparison.reverse_splicing = comparison.reverse_direction
                    comparison.reverse_deg = comparison.reverse_direction
                    if comparison.display_name is None and treatment and control:
                        comparison.display_name = f"{treatment} vs {control}"
                    self._derive_source_names_from_groups(comparison, treatment, control, contrast)

        project.available_comparisons = sorted(project.available_comparisons, key=lambda item: item.comparison_id)

    def _derive_source_names_from_groups(
        self,
        comparison: ComparisonDefinition,
        treatment: str | None,
        control: str | None,
        contrast_name: str,
    ) -> None:
        if comparison.sashimi_name is None:
            comparison.sashimi_name = contrast_name
        if comparison.bam_group_name is None:
            comparison.bam_group_name = contrast_name
        if not treatment or not control:
            return

        forward_vs = f"{treatment}_vs_{control}"
        reverse_vs = f"{control}_vs_{treatment}"
        dashed = f"{treatment}-{control}"

        if comparison.rmats_name is None:
            comparison.rmats_name = dashed
        if comparison.suppa_name is None:
            comparison.suppa_name = forward_vs
        if comparison.dexseq_name is None:
            comparison.dexseq_name = forward_vs
        if comparison.dtu_name is None:
            comparison.dtu_name = forward_vs
        if comparison.quant_name is None:
            comparison.quant_name = forward_vs
        if comparison.deg_name is None:
            comparison.deg_name = reverse_vs

    def _discover_isoform_samples(self, project: ProjectConfig, index: ScanIndex) -> list[IsoformSampleDefinition]:
        samples: list[IsoformSampleDefinition] = []
        if project.quant_root and project.quant_root.exists():
            quant_files = sorted(project.quant_root.rglob("quant.sf"))
        else:
            quant_files = sorted(path for path in index.files if path.name == "quant.sf")
        for quant_file in quant_files:
            sample_id = quant_file.parent.name
            sample_group = self._normalize_sample_group(sample_id)
            comparison_id = self._infer_comparison_from_sample_name(sample_id, project)
            samples.append(
                IsoformSampleDefinition(
                    sample_id=sample_id,
                    quant_sf=quant_file,
                    sample_group=sample_group,
                    comparison_id=comparison_id,
                )
            )
        return samples

    def _attach_isoform_samples_to_comparisons(self, project: ProjectConfig) -> None:
        group_to_dirs: dict[str, list[Path]] = {}
        for sample in project.isoform_samples:
            if not sample.sample_group:
                continue
            group_to_dirs.setdefault(sample.sample_group, []).append(sample.quant_sf.parent)

        for comparison in project.available_comparisons:
            groups = {
                comparison.source_experiment_group or comparison.experiment_group,
                comparison.source_control_group or comparison.control_group,
            }
            groups = {group for group in groups if group}
            if not groups:
                continue
            if groups.issubset(group_to_dirs.keys()):
                comparison.has_quant = True
                if comparison.quant_dir is None:
                    first_group = sorted(groups)[0]
                    comparison.quant_dir = group_to_dirs[first_group][0]

    def _infer_comparison_from_sample_name(
        self,
        sample_id: str,
        project: ProjectConfig,
    ) -> str | None:
        normalized = self._normalize_sample_group(sample_id)

        comparisons = project.available_comparisons
        for comparison in comparisons:
            experiment = comparison.source_experiment_group or comparison.experiment_group
            control = comparison.source_control_group or comparison.control_group
            if normalized in {experiment, control}:
                return comparison.comparison_id
        return None

    def _normalize_sample_group(self, sample_id: str) -> str:
        normalized = sample_id
        for suffix in ("_rep_1", "_rep_2", "_rep_3", "_rep_4", "_rep_5", "_rep_6"):
            if normalized.endswith(suffix):
                normalized = normalized[: -len(suffix)]
                break
        return normalized

    def _apply_saved_config(self, project: ProjectConfig) -> None:
        if project.config_path is None or not project.config_path.exists():
            return
        with project.config_path.open("r", encoding="utf-8") as handle:
            data = yaml.safe_load(handle) or {}

        comparisons_by_id = {item.comparison_id: item for item in project.available_comparisons}
        for item in data.get("comparisons", []):
            comparison_id = item.get("comparison_id")
            if not comparison_id:
                continue
            comparison = comparisons_by_id.get(comparison_id)
            if comparison is None:
                has_meaningful_saved_state = any(
                    item.get(field)
                    for field in (
                        "display_name",
                        "rmats_name",
                        "deg_name",
                        "suppa_name",
                        "dexseq_name",
                        "dtu_name",
                        "quant_name",
                        "sashimi_name",
                        "bam_group_name",
                        "source_experiment_group",
                        "source_control_group",
                        "experiment_group",
                        "control_group",
                    )
                ) or bool(item.get("reverse_direction")) or bool(item.get("reverse_splicing")) or bool(item.get("reverse_deg"))
                if not has_meaningful_saved_state:
                    continue
                comparison = ComparisonDefinition(comparison_id=comparison_id)
                project.available_comparisons.append(comparison)
                comparisons_by_id[comparison_id] = comparison
            comparison.display_name = item.get("display_name")
            comparison.biological_question = item.get("biological_question") or comparison.biological_question
            comparison.rmats_name = item.get("rmats_name") or comparison.rmats_name
            comparison.deg_name = item.get("deg_name") or comparison.deg_name
            comparison.suppa_name = item.get("suppa_name") or comparison.suppa_name
            comparison.dexseq_name = item.get("dexseq_name") or comparison.dexseq_name
            comparison.dtu_name = item.get("dtu_name") or comparison.dtu_name
            comparison.quant_name = item.get("quant_name") or comparison.quant_name
            comparison.sashimi_name = item.get("sashimi_name") or comparison.sashimi_name
            comparison.bam_group_name = item.get("bam_group_name") or comparison.bam_group_name
            comparison.source_experiment_group = item.get("source_experiment_group") or comparison.source_experiment_group
            comparison.source_control_group = item.get("source_control_group") or comparison.source_control_group
            comparison.experiment_group = item.get("experiment_group") or comparison.experiment_group
            comparison.control_group = item.get("control_group") or comparison.control_group
            comparison.source_direction = item.get("source_direction") or comparison.source_direction
            comparison.analysis_direction = item.get("analysis_direction") or comparison.analysis_direction
            comparison.reverse_direction = bool(item.get("reverse_direction", comparison.reverse_direction))
            comparison.reverse_splicing = bool(item.get("reverse_splicing", comparison.reverse_direction))
            comparison.reverse_deg = bool(item.get("reverse_deg", comparison.reverse_direction))
            comparison.expected_direction_notes = item.get("expected_direction_notes") or comparison.expected_direction_notes
            comparison.known_positive_control_genes = list(item.get("known_positive_control_genes") or comparison.known_positive_control_genes)
            comparison.known_negative_control_genes = list(item.get("known_negative_control_genes") or comparison.known_negative_control_genes)
            comparison.output_prefix = item.get("output_prefix") or comparison.output_prefix

        sample_overrides = data.get("isoform_samples", []) or []
        sample_lookup = {
            (sample.sample_id, str(sample.quant_sf)): sample
            for sample in project.isoform_samples
        }
        for item in sample_overrides:
            sample_id = item.get("sample_id")
            quant_sf = item.get("quant_sf")
            if not sample_id or not quant_sf:
                continue
            sample = sample_lookup.get((sample_id, quant_sf))
            if sample is None:
                path = Path(quant_sf)
                sample = IsoformSampleDefinition(
                    sample_id=sample_id,
                    quant_sf=path,
                    sample_group=item.get("sample_group"),
                    comparison_id=item.get("comparison_id"),
                )
                project.isoform_samples.append(sample)
                sample_lookup[(sample_id, str(path))] = sample
            sample.sample_group = item.get("sample_group") or sample.sample_group
            sample.comparison_id = item.get("comparison_id") or sample.comparison_id
            sample.condition = item.get("condition") or sample.condition
            sample.batch = item.get("batch") or sample.batch
            sample.include = bool(item.get("include", sample.include))

        project.selected_rmats_mode = data.get("selected_rmats_mode", project.selected_rmats_mode)
        output_root = data.get("output_root")
        if output_root:
            project.output_root = Path(output_root)
        project.selected_comparison_ids = list(data.get("selected_comparison_ids", project.selected_comparison_ids))
        project.comparison_order_ids = list(data.get("comparison_order_ids", project.comparison_order_ids))
        project.program_comparison_ids = list(data.get("program_comparison_ids", project.program_comparison_ids))
        project.selected_analysis_modules = list(data.get("selected_analysis_modules", project.selected_analysis_modules))
        project.confirmed = bool(data.get("confirmed", False))
        project.pairing_confirmed = bool(data.get("pairing_confirmed", False))
        project.comparison_sets_confirmed = bool(data.get("comparison_sets_confirmed", False))
        project.visualization_groups_confirmed = bool(data.get("visualization_groups_confirmed", False))
        project.comparison_pairs = [
            ComparisonPairDefinition(
                pair_id=item.get("pair_id"),
                comparison_a=item.get("comparison_a"),
                comparison_b=item.get("comparison_b"),
                display_name=item.get("display_name"),
                experiment_group=item.get("experiment_group"),
                control_group=item.get("control_group"),
                enabled=bool(item.get("enabled", True)),
            )
            for item in (data.get("comparison_pairs") or [])
            if item.get("pair_id")
        ]
        project.visualization_groups = [
            VisualizationGroupDefinition(
                group_id=item.get("group_id"),
                display_name=item.get("display_name"),
                comparison_ids=list(item.get("comparison_ids") or []),
                enabled=bool(item.get("enabled", True)),
            )
            for item in (data.get("visualization_groups") or [])
            if item.get("group_id")
        ]
        project.shortlist_genes = list(data.get("shortlist_genes") or [])
        project.blacklist_genes = list(data.get("blacklist_genes") or [])
        project.candidate_selection_genes = {
            str(key): [str(item) for item in (value or []) if str(item).strip()]
            for key, value in (data.get("candidate_selection_genes") or {}).items()
        }
        project.candidate_blacklist_genes = {
            str(key): [str(item) for item in (value or []) if str(item).strip()]
            for key, value in (data.get("candidate_blacklist_genes") or {}).items()
        }
        project.removed_comparison_ids = list(data.get("removed_comparison_ids") or [])
        if project.removed_comparison_ids:
            removed = set(project.removed_comparison_ids)
            project.available_comparisons = [
                comparison
                for comparison in project.available_comparisons
                if comparison.comparison_id not in removed
            ]
            project.selected_comparison_ids = [
                comparison_id
                for comparison_id in project.selected_comparison_ids
                if comparison_id not in removed
            ]
            project.comparison_order_ids = [
                comparison_id
                for comparison_id in project.comparison_order_ids
                if comparison_id not in removed
            ]
            project.program_comparison_ids = [
                comparison_id
                for comparison_id in project.program_comparison_ids
                if comparison_id not in removed
            ]
            project.comparison_pairs = [
                pair
                for pair in project.comparison_pairs
                if pair.comparison_a not in removed and pair.comparison_b not in removed
            ]
            for group in project.visualization_groups:
                group.comparison_ids = [
                    comparison_id
                    for comparison_id in group.comparison_ids
                    if comparison_id not in removed
                ]
            project.visualization_groups = [
                group for group in project.visualization_groups if group.comparison_ids
            ]
            project.candidate_selection_genes = {
                comparison_id: values
                for comparison_id, values in project.candidate_selection_genes.items()
                if comparison_id not in removed
            }
            project.candidate_blacklist_genes = {
                comparison_id: values
                for comparison_id, values in project.candidate_blacklist_genes.items()
                if comparison_id not in removed
            }
        thresholds = data.get("thresholds") or {}
        if thresholds:
            project.tool_paths["saved_thresholds"] = thresholds
        saved_tool_paths = data.get("tool_paths") or {}
        if saved_tool_paths:
            project.tool_paths.update(saved_tool_paths)
        saved_input_paths = data.get("input_paths") or {}
        if saved_input_paths:
            project.input_paths.update({key: str(value) for key, value in saved_input_paths.items() if value})
        project.available_comparisons = sorted(project.available_comparisons, key=lambda entry: entry.comparison_id)

    def _reconcile_source_mappings(self, project: ProjectConfig) -> None:
        comparisons_by_id = {item.comparison_id: item for item in project.available_comparisons}

        for comparison in list(project.available_comparisons):
            self._copy_from_source(comparison, comparisons_by_id.get(comparison.source_name("rmats")), "rmats")
            self._copy_from_source(comparison, comparisons_by_id.get(comparison.source_name("deg")), "deg")
            self._copy_from_source(comparison, comparisons_by_id.get(comparison.source_name("quant")), "quant")

        referenced = {
            name
            for comparison in project.available_comparisons
            for name in [
                comparison.rmats_name,
                comparison.deg_name,
                comparison.suppa_name,
                comparison.dexseq_name,
                comparison.dtu_name,
                comparison.quant_name,
                comparison.sashimi_name,
                comparison.bam_group_name,
            ]
            if name and name != comparison.comparison_id
        }
        project.available_comparisons = [
            comparison
            for comparison in project.available_comparisons
            if not (
                comparison.comparison_id in referenced
                and comparison.display_name is None
                and comparison.rmats_name is None
                and comparison.deg_name is None
                and comparison.suppa_name is None
                and comparison.dexseq_name is None
                and comparison.dtu_name is None
                and comparison.quant_name is None
                and comparison.sashimi_name is None
                and comparison.bam_group_name is None
            )
        ]

    def _merge_alias_comparisons(self, project: ProjectConfig) -> None:
        canonical_candidates = [
            comparison
            for comparison in project.available_comparisons
            if (
                comparison.rmats_path is not None
                or comparison.source_experiment_group
                or comparison.source_control_group
                or comparison.display_name
            )
        ]
        if not canonical_candidates:
            project.available_comparisons = sorted(project.available_comparisons, key=lambda item: item.comparison_id)
            return

        alias_to_canonical: dict[str, ComparisonDefinition] = {}
        for comparison in sorted(canonical_candidates, key=self._comparison_priority, reverse=True):
            for alias in self._comparison_aliases(comparison):
                current = alias_to_canonical.get(alias)
                if current is None or self._comparison_priority(comparison) > self._comparison_priority(current):
                    alias_to_canonical[alias] = comparison

        merged: list[ComparisonDefinition] = []
        seen_ids: set[str] = set()

        for comparison in project.available_comparisons:
            target = None
            for alias in self._comparison_aliases(comparison):
                target = alias_to_canonical.get(alias)
                if target is not None:
                    break
            if target is None or target is comparison:
                if comparison.comparison_id not in seen_ids:
                    merged.append(comparison)
                    seen_ids.add(comparison.comparison_id)
                continue
            self._absorb_comparison(target, comparison)

        project.available_comparisons = sorted(merged, key=lambda item: item.comparison_id)

    def _comparison_priority(self, comparison: ComparisonDefinition) -> tuple[int, int, int, int]:
        return (
            1 if comparison.rmats_path is not None else 0,
            1 if comparison.source_experiment_group or comparison.source_control_group else 0,
            1 if comparison.has_deg else 0,
            1 if comparison.has_quant else 0,
        )

    def _comparison_aliases(self, comparison: ComparisonDefinition) -> set[str]:
        aliases = {
            value
            for value in [
                comparison.comparison_id,
                comparison.display_name,
                comparison.rmats_name,
                comparison.deg_name,
                comparison.suppa_name,
                comparison.dexseq_name,
                comparison.dtu_name,
                comparison.quant_name,
                comparison.sashimi_name,
                comparison.bam_group_name,
            ]
            if value
        }
        if comparison.experiment_group and comparison.control_group:
            treatment = comparison.experiment_group
            control = comparison.control_group
            aliases.update(
                {
                    f"{treatment}_vs_{control}",
                    f"{control}_vs_{treatment}",
                    f"{treatment}-{control}",
                    f"{control}-{treatment}",
                    f"{treatment} vs {control}",
                }
            )
        return aliases

    def _absorb_comparison(self, target: ComparisonDefinition, source: ComparisonDefinition) -> None:
        if target.rmats_path is None:
            target.rmats_path = source.rmats_path
        if target.deg_path is None:
            target.deg_path = source.deg_path
        if target.quant_dir is None:
            target.quant_dir = source.quant_dir
        if not target.detected_rmats_modes and source.detected_rmats_modes:
            target.detected_rmats_modes = list(source.detected_rmats_modes)
        target.has_deg = target.has_deg or source.has_deg
        target.has_quant = target.has_quant or source.has_quant
        target.rmats_name = target.rmats_name or source.rmats_name or source.comparison_id
        target.deg_name = target.deg_name or source.deg_name or source.comparison_id
        target.suppa_name = target.suppa_name or source.suppa_name or source.comparison_id
        target.dexseq_name = target.dexseq_name or source.dexseq_name or source.comparison_id
        target.dtu_name = target.dtu_name or source.dtu_name or source.comparison_id
        target.quant_name = target.quant_name or source.quant_name or source.comparison_id
        target.sashimi_name = target.sashimi_name or source.sashimi_name or source.comparison_id
        target.bam_group_name = target.bam_group_name or source.bam_group_name or source.comparison_id
        target.source_experiment_group = target.source_experiment_group or source.source_experiment_group
        target.source_control_group = target.source_control_group or source.source_control_group
        target.experiment_group = target.experiment_group or source.experiment_group
        target.control_group = target.control_group or source.control_group

    def _copy_from_source(
        self,
        target: ComparisonDefinition,
        source: ComparisonDefinition | None,
        source_type: str,
    ) -> None:
        if source is None or source is target:
            return
        if source_type == "rmats" and target.rmats_path is None:
            target.rmats_path = source.rmats_path
            target.detected_rmats_modes = list(source.detected_rmats_modes)
        if source_type == "deg" and target.deg_path is None:
            target.deg_path = source.deg_path
            target.has_deg = source.has_deg
        if source_type == "quant" and target.quant_dir is None:
            target.quant_dir = source.quant_dir
            target.has_quant = source.has_quant

    def _resolve_missing_source_paths(self, project: ProjectConfig) -> None:
        self._resolve_missing_deg_paths(project)

    def _resolve_missing_deg_paths(self, project: ProjectConfig) -> None:
        if project.deg_root is None or not project.deg_root.exists():
            return

        deg_files = {
            path.stem.replace(".deseq2.results", ""): path
            for path in project.deg_root.rglob("*.deseq2.results.tsv")
        }

        for comparison in project.available_comparisons:
            if comparison.deg_path is not None:
                continue
            candidate_names = []
            for name in (
                comparison.deg_name,
                comparison.comparison_id,
                f"{comparison.experiment_group}_vs_{comparison.control_group}"
                if comparison.experiment_group and comparison.control_group
                else None,
                f"{comparison.control_group}_vs_{comparison.experiment_group}"
                if comparison.experiment_group and comparison.control_group
                else None,
            ):
                if name and name not in candidate_names:
                    candidate_names.append(name)
                swapped = self._swap_vs_name(name) if name else None
                if swapped and swapped not in candidate_names:
                    candidate_names.append(swapped)

            for candidate in candidate_names:
                matched = deg_files.get(candidate)
                if matched is not None:
                    comparison.deg_path = matched
                    comparison.has_deg = True
                    if comparison.deg_name is None:
                        comparison.deg_name = candidate
                    break

    @staticmethod
    def _swap_vs_name(name: str | None) -> str | None:
        if not name or "_vs_" not in name:
            return None
        left, right = name.split("_vs_", 1)
        if not left or not right:
            return None
        return f"{right}_vs_{left}"

    def _detect_rmats_modes(self, comparison_dir: Path) -> list[str]:
        rmats_post = comparison_dir / "rmats_post"
        if not rmats_post.exists():
            return []
        modes: list[str] = []
        if all((rmats_post / f"{event}.MATS.JC.txt").exists() for event in EVENT_TYPES):
            modes.append("JC")
        if all((rmats_post / f"{event}.MATS.JCEC.txt").exists() for event in EVENT_TYPES):
            modes.append("JCEC")
        return modes

    def _build_scan_messages(self, project: ProjectConfig) -> list[str]:
        messages = []
        if project.rmats_root is None:
            messages.append("rMATS root not detected.")
        if project.deg_root is None:
            messages.append("DEG root not detected.")
        if project.counts_path is None:
            messages.append("Normalized counts file not detected.")
        if not project.isoform_samples:
            messages.append("No quant.sf files detected under the selected project root.")
        if project.input_paths:
            messages.append("Manual input paths are active; targeted scanning was used.")
        if project.contrastsheet_path is not None:
            messages.append(f"Contrast sheet detected: {project.contrastsheet_path}")
        if not project.available_comparisons:
            messages.append("No comparisons detected.")
        return messages
