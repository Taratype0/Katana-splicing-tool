from __future__ import annotations

import json
from pathlib import Path
from datetime import datetime

import yaml
import pandas as pd
import matplotlib.pyplot as plt
import importlib.util

from src.adapters.jutils_adapter import JutilsAdapter
from src.adapters.isoformswitch_runner import IsoformSwitchRunner
from src.adapters.cards_adapter import CardsAdapter
from src.adapters.sashimi_adapter import SashimiAdapter
from src.analysis.common import EVENT_TYPES, event_uid_from_row, load_annotation_gene_symbol_lookup, normalize_gene_id
from src.analysis.mechanism_support import MechanismSupportAnalyzer
from src.analysis.program_comparison import SplicingProgramComparator
from src.analysis.splicing_landscape import SplicingLandscapeAnalyzer
from src.analysis.tx_splicing_integration import TxSplicingIntegrator
from src.runtime_paths import bundled_configs_root, bundled_software_root
from src.services.shortlist_service import ShortlistService
from src.adapters.isoformswitch_adapter import IsoformSwitchAdapter
from src.domain.models.analysis import AnalysisRequest, AnalysisResults, AnalysisRunState, ThresholdConfig
from src.domain.models.project import ComparisonPairDefinition, IsoformSampleDefinition, ProjectConfig, VisualizationGroupDefinition
from src.io.scanner import ProjectScanner


class ProjectService:
    ANALYSIS_MODULES = [
        "splicing_landscape",
        "program_comparison",
        "mechanism_support",
        "tx_splicing_integration",
    ]

    def __init__(self) -> None:
        self.current_project: ProjectConfig | None = None
        self.scanner = ProjectScanner()
        self.splicing_landscape = SplicingLandscapeAnalyzer()
        self.program_comparison = SplicingProgramComparator()
        self.mechanism_support = MechanismSupportAnalyzer()
        self.tx_splicing_integration = TxSplicingIntegrator()
        self.isoform_adapter = IsoformSwitchAdapter()
        self.cards_adapter = CardsAdapter()
        self.jutils_adapter: JutilsAdapter | None = None
        self.sashimi_adapter: SashimiAdapter | None = None
        self.isoform_runner: IsoformSwitchRunner | None = None
        self.shortlist_service = ShortlistService()
        self.default_thresholds = self._load_default_thresholds()
        self.run_state = AnalysisRunState()
        self.current_thresholds = ThresholdConfig(**self.default_thresholds["presets"]["default"])
        self._last_jutils_report: dict[str, object] = {
            "status": "idle",
            "command_lines": [],
            "output_directory": "",
            "generated_files": [],
            "stdout": "",
            "stderr": "",
            "error_log": "",
        }
        self._last_sashimi_failures = pd.DataFrame()
        self._sashimi_event_catalog_cache: dict[str, pd.DataFrame] = {}
        self._cross_as_pattern_cache: dict[tuple[str, float, float], dict[str, pd.DataFrame]] = {}
        self._deg_expression_cache: dict[str, pd.DataFrame] = {}
        self._expression_support_cache = pd.DataFrame()
        self._expression_support_loaded = False
        self._module_states = self._new_module_states()

    def _emit_progress(self, callback, message: str) -> None:
        if callback is not None:
            callback(message)

    def _log_heavy(self, message: str) -> None:
        print(f"[katana-heavy] {message}")

    def _new_module_states(self) -> dict[str, dict[str, object]]:
        return {
            "splicing_landscape": {"status": "not_run", "message": "Not run."},
            "program_comparison": {"status": "not_run", "message": "Not run."},
            "mechanism_support": {"status": "not_run", "message": "Not run."},
            "tx_splicing_integration": {"status": "not_run", "message": "Not run."},
            "candidate_ranking": {"status": "not_run", "message": "Not run."},
            "candidate_selection": {"status": "not_run", "message": "Not run."},
            "candidate_heatmap": {"status": "not_run", "message": "Not run."},
            "candidate_followup": {"status": "not_run", "message": "Not run."},
            "cross_comparison_candidates": {"status": "not_run", "message": "Not run."},
            "cross_comparison_as_patterns": {"status": "not_run", "message": "Not run."},
            "sashimi_manifest": {"status": "not_run", "message": "Not run."},
            "sashimi_plot": {"status": "not_run", "message": "Not run."},
            "jutils_heatmap": {"status": "not_run", "message": "Not run."},
            "jutils_pca": {"status": "not_run", "message": "Not run."},
            "isoform": {"status": "not_run", "message": "Not run."},
        }

    def _reset_module_states(self) -> None:
        self._module_states = self._new_module_states()

    def _mark_module_state(self, module: str, status: str, message: str, **extra) -> None:
        current = dict(self._module_states.get(module, {}))
        current.update({"status": status, "message": message, **extra})
        self._module_states[module] = current

    def module_state(self, module: str) -> dict[str, object]:
        return dict(self._module_states.get(module, {"status": "not_run", "message": "Not run."}))

    def _copy_analysis_results(self, results: AnalysisResults | None = None) -> AnalysisResults:
        source = results or self.run_state.results
        return AnalysisResults(
            splicing_landscape=source.splicing_landscape.copy(),
            program_events=source.program_events.copy(),
            program_summary=source.program_summary.copy(),
            mechanism_support=source.mechanism_support.copy(),
            tx_splicing_gene_table=source.tx_splicing_gene_table.copy(),
            tx_splicing_summary=source.tx_splicing_summary.copy(),
            candidate_gene_table=source.candidate_gene_table.copy(),
            cross_comparison_candidate_matrix=source.cross_comparison_candidate_matrix.copy(),
            ranked_candidates=source.ranked_candidates.copy(),
            shortlist_dp=source.shortlist_dp.copy(),
            shortlist_ko=source.shortlist_ko.copy(),
            cards_shortlist=source.cards_shortlist.copy(),
            cards_expression_support=source.cards_expression_support.copy(),
            sashimi_manifest=source.sashimi_manifest.copy(),
            isoform_manifest=source.isoform_manifest.copy(),
            jutils_manifest=source.jutils_manifest.copy(),
        )

    def load_project(
        self,
        root_dir: str | Path,
        input_paths: dict[str, str] | None = None,
        force_auto_scan: bool = False,
        progress_callback=None,
    ) -> ProjectConfig:
        self._emit_progress(progress_callback, "load config")
        self.run_state = AnalysisRunState()
        self._reset_module_states()
        self._sashimi_event_catalog_cache = {}
        self._cross_as_pattern_cache = {}
        self._deg_expression_cache = {}
        self._expression_support_cache = pd.DataFrame()
        self._expression_support_loaded = False
        self.current_project = self.scanner.scan(
            root_dir,
            input_paths=input_paths,
            force_auto_scan=force_auto_scan,
            progress_callback=progress_callback,
        )
        self._emit_progress(progress_callback, "load comparison config")
        self.current_project.confirmed = False
        self.current_project.pairing_confirmed = False
        self.current_project.comparison_sets_confirmed = False
        self.current_project.visualization_groups_confirmed = False
        self._emit_progress(progress_callback, "load pairing direction")
        self._refresh_tool_adapters()
        if not self.current_project.selected_comparison_ids:
            self.current_project.selected_comparison_ids = [
                item.comparison_id for item in self.current_project.available_comparisons[:4]
            ]
        if not self.current_project.selected_analysis_modules:
            self.current_project.selected_analysis_modules = list(self.ANALYSIS_MODULES)
        self._ensure_default_comparison_pairs()
        self._ensure_default_visualization_groups()
        saved_thresholds = self.current_project.tool_paths.get("saved_thresholds")
        if isinstance(saved_thresholds, dict) and saved_thresholds:
            self.current_thresholds = ThresholdConfig(**saved_thresholds)
        self._emit_progress(progress_callback, "load cached candidate selection state")
        self._autodetect_embedded_tools()
        self._emit_progress(progress_callback, "finish load project")
        return self.current_project

    def confirm_project(self) -> None:
        if self.current_project is None:
            raise RuntimeError("No project loaded.")
        self.current_project.confirmed = True
        self.current_project.pairing_confirmed = False
        self.current_project.comparison_sets_confirmed = False
        self.current_project.visualization_groups_confirmed = False
        self.save_project_config()

    def invalidate_project_confirmation(self) -> None:
        if self.current_project is None:
            return
        if self.current_project.confirmed or self.current_project.pairing_confirmed:
            self.current_project.confirmed = False
            self.current_project.pairing_confirmed = False
            self.current_project.visualization_groups_confirmed = False
            self.save_project_config()

    def confirm_pairing(self) -> None:
        if self.current_project is None:
            raise RuntimeError("No project loaded.")
        if not self.current_project.confirmed:
            raise RuntimeError("Confirm the project before confirming pairing.")
        for comparison in self.current_project.available_comparisons:
            if comparison.experiment_group and comparison.control_group:
                comparison.analysis_direction = "finalized"
            comparison.reverse_direction = comparison.reverse_splicing or comparison.reverse_deg
        self.current_project.pairing_confirmed = True
        self.current_project.comparison_sets_confirmed = False
        self.current_project.visualization_groups_confirmed = False
        self.save_project_config()

    def invalidate_pairing_confirmation(self) -> None:
        if self.current_project is None:
            return
        if self.current_project.pairing_confirmed or self.current_project.comparison_sets_confirmed:
            self.current_project.pairing_confirmed = False
            self.current_project.comparison_sets_confirmed = False
            self.current_project.visualization_groups_confirmed = False
            self.save_project_config()

    def invalidate_analysis_results(self, reason: str = "Analysis results were cleared because pairing or source settings changed.") -> None:
        self.run_state = AnalysisRunState(
            status="idle",
            message=reason,
        )
        self._cross_as_pattern_cache = {}
        self._reset_module_states()

    def confirm_comparison_sets(self) -> None:
        if self.current_project is None:
            raise RuntimeError("No project loaded.")
        if not self.current_project.pairing_confirmed:
            raise RuntimeError("Confirm Pairing before confirming Comparison Sets.")
        self.current_project.comparison_sets_confirmed = True
        self.current_project.visualization_groups_confirmed = False
        self.save_project_config()

    def invalidate_comparison_sets_confirmation(self) -> None:
        if self.current_project is None:
            return
        if self.current_project.comparison_sets_confirmed:
            self.current_project.comparison_sets_confirmed = False
            self.current_project.visualization_groups_confirmed = False
            self.save_project_config()

    def confirm_visualization_groups(self) -> None:
        if self.current_project is None:
            raise RuntimeError("No project loaded.")
        if not self.current_project.comparison_sets_confirmed:
            raise RuntimeError("Confirm Comparison Sets before confirming Visual Groups.")
        self.current_project.visualization_groups_confirmed = True
        self.save_project_config()

    def invalidate_visualization_groups_confirmation(self) -> None:
        if self.current_project is None:
            return
        if self.current_project.visualization_groups_confirmed:
            self.current_project.visualization_groups_confirmed = False
            self.save_project_config()

    def pairing_direction_mismatches(self) -> list[str]:
        return []

    def rename_comparison(self, comparison_id: str, display_name: str) -> None:
        if self.current_project is None:
            raise RuntimeError("No project loaded.")
        for comparison in self.current_project.available_comparisons:
            if comparison.comparison_id == comparison_id:
                comparison.display_name = display_name.strip() or None
                self.save_project_config()
                return
        raise KeyError(f"Comparison not found: {comparison_id}")

    def update_comparison_mapping(
        self,
        comparison_id: str,
        *,
        rmats_name: str | None = None,
        deg_name: str | None = None,
        suppa_name: str | None = None,
        dexseq_name: str | None = None,
        dtu_name: str | None = None,
        quant_name: str | None = None,
        sashimi_name: str | None = None,
        bam_group_name: str | None = None,
    ) -> None:
        if self.current_project is None:
            raise RuntimeError("No project loaded.")
        comparison = next(
            (item for item in self.current_project.available_comparisons if item.comparison_id == comparison_id),
            None,
        )
        if comparison is None:
            raise KeyError(f"Comparison not found: {comparison_id}")
        if rmats_name is not None:
            comparison.rmats_name = rmats_name.strip() or None
        if deg_name is not None:
            comparison.deg_name = deg_name.strip() or None
        if suppa_name is not None:
            comparison.suppa_name = suppa_name.strip() or None
        if dexseq_name is not None:
            comparison.dexseq_name = dexseq_name.strip() or None
        if dtu_name is not None:
            comparison.dtu_name = dtu_name.strip() or None
        if quant_name is not None:
            comparison.quant_name = quant_name.strip() or None
        if sashimi_name is not None:
            comparison.sashimi_name = sashimi_name.strip() or None
        if bam_group_name is not None:
            comparison.bam_group_name = bam_group_name.strip() or None
        self.invalidate_analysis_results(
            "Analysis results were cleared because comparison file mappings changed. Re-run analysis to refresh downstream views."
        )
        self.save_project_config()

    def update_comparison_settings(
        self,
        comparison_id: str,
        *,
        experiment_group: str | None = None,
        control_group: str | None = None,
        source_direction: str | None = None,
        analysis_direction: str | None = None,
        reverse_direction: bool | None = None,
        reverse_splicing: bool | None = None,
        reverse_deg: bool | None = None,
    ) -> None:
        if self.current_project is None:
            raise RuntimeError("No project loaded.")
        comparison = next(
            (item for item in self.current_project.available_comparisons if item.comparison_id == comparison_id),
            None,
        )
        if comparison is None:
            raise KeyError(f"Comparison not found: {comparison_id}")
        if experiment_group is not None:
            comparison.experiment_group = experiment_group.strip() or None
        if control_group is not None:
            comparison.control_group = control_group.strip() or None
        if source_direction is not None:
            comparison.source_direction = source_direction.strip() or "original"
        if analysis_direction is not None:
            comparison.analysis_direction = analysis_direction.strip() or "original"
        if reverse_direction is not None:
            comparison.reverse_direction = reverse_direction
        if reverse_splicing is not None:
            comparison.reverse_splicing = reverse_splicing
        if reverse_deg is not None:
            comparison.reverse_deg = reverse_deg
        self._sashimi_event_catalog_cache = {}
        self._cross_as_pattern_cache = {}
        self.invalidate_analysis_results(
            "Analysis results were cleared because pairing direction settings changed. Re-run analysis so standardized log2FC and dPSI use the updated Pairing flips."
        )
        self.save_project_config()

    def update_isoform_sample(
        self,
        sample_id: str,
        quant_sf: str,
        *,
        sample_group: str | None = None,
        comparison_id: str | None = None,
        condition: str | None = None,
        batch: str | None = None,
        include: bool | None = None,
    ) -> None:
        if self.current_project is None:
            raise RuntimeError("No project loaded.")
        sample = next(
            (
                item
                for item in self.current_project.isoform_samples
                if item.sample_id == sample_id and str(item.quant_sf) == quant_sf
            ),
            None,
        )
        if sample is None:
            sample = IsoformSampleDefinition(
                sample_id=sample_id,
                quant_sf=Path(quant_sf),
                sample_group=sample_group,
                comparison_id=comparison_id,
            )
            self.current_project.isoform_samples.append(sample)
        if sample_group is not None:
            sample.sample_group = sample_group or None
        if comparison_id is not None:
            sample.comparison_id = comparison_id or None
        if condition is not None:
            sample.condition = condition or None
        if batch is not None:
            sample.batch = batch or None
        if include is not None:
            sample.include = include
        self.save_project_config()

    def import_isoform_metadata(self, metadata_path: str | Path) -> tuple[int, int]:
        if self.current_project is None:
            raise RuntimeError("No project loaded.")
        path = Path(metadata_path)
        if not path.exists():
            raise RuntimeError(f"Isoform metadata file not found: {path}")
        sep = "\t" if path.suffix.lower() in {".tsv", ".txt"} else ","
        frame = pd.read_csv(path, sep=sep, low_memory=False)
        if "sample_id" not in frame.columns:
            raise RuntimeError("Isoform metadata must contain a 'sample_id' column.")

        updated = 0
        created = 0
        sample_lookup = {(sample.sample_id, str(sample.quant_sf)): sample for sample in self.current_project.isoform_samples}
        by_sample_id: dict[str, list[IsoformSampleDefinition]] = {}
        for sample in self.current_project.isoform_samples:
            by_sample_id.setdefault(sample.sample_id, []).append(sample)

        for _, row in frame.iterrows():
            sample_id = str(row.get("sample_id", "")).strip()
            if not sample_id:
                continue
            quant_sf_value = str(row.get("quant_sf", "")).strip()
            target = None
            if quant_sf_value:
                target = sample_lookup.get((sample_id, quant_sf_value))
            if target is None:
                candidates = by_sample_id.get(sample_id, [])
                if len(candidates) == 1:
                    target = candidates[0]
            if target is None and quant_sf_value:
                target = IsoformSampleDefinition(
                    sample_id=sample_id,
                    quant_sf=Path(quant_sf_value),
                )
                self.current_project.isoform_samples.append(target)
                sample_lookup[(sample_id, quant_sf_value)] = target
                by_sample_id.setdefault(sample_id, []).append(target)
                created += 1
            if target is None:
                continue

            if "comparison_id" in frame.columns:
                value = str(row.get("comparison_id", "")).strip()
                target.comparison_id = value or target.comparison_id
            if "sample_group" in frame.columns:
                value = str(row.get("sample_group", "")).strip()
                target.sample_group = value or target.sample_group
            if "condition" in frame.columns:
                value = str(row.get("condition", "")).strip()
                target.condition = value or None
            if "batch" in frame.columns:
                value = str(row.get("batch", "")).strip()
                target.batch = value or None
            if "include" in frame.columns:
                raw_value = str(row.get("include", "")).strip().lower()
                target.include = raw_value not in {"0", "false", "no", "n", ""}
            updated += 1

        self.save_project_config()
        return updated, created

    def export_isoform_design_template(self, output_path: str | Path) -> Path:
        if self.current_project is None:
            raise RuntimeError("No project loaded.")
        samples = self.current_project.isoform_samples
        if not samples:
            raise RuntimeError("No isoform samples detected.")
        frame = self.isoform_adapter.build_import_manifest_frame(samples)
        if "condition" not in frame.columns:
            frame["condition"] = ""
        if "batch" not in frame.columns:
            frame["batch"] = ""
        if "include" not in frame.columns:
            include_lookup = {(sample.sample_id, str(sample.quant_sf)): sample.include for sample in samples}
            frame["include"] = [
                include_lookup.get((str(row["sample_id"]), str(row["quant_sf"])), True)
                for _, row in frame.iterrows()
            ]
        if "sample_group" not in frame.columns:
            sample_group_lookup = {(sample.sample_id, str(sample.quant_sf)): sample.sample_group for sample in samples}
            frame["sample_group"] = [
                sample_group_lookup.get((str(row["sample_id"]), str(row["quant_sf"])), "")
                for _, row in frame.iterrows()
            ]
        columns = [column for column in ["sample_id", "quant_sf", "sample_group", "comparison_id", "condition", "batch", "include"] if column in frame.columns]
        out = Path(output_path)
        out.parent.mkdir(parents=True, exist_ok=True)
        frame[columns].to_csv(out, sep="\t", index=False)
        return out

    def build_default_request(self) -> AnalysisRequest:
        comparison_ids = []
        program_comparison_ids = []
        if self.current_project is not None:
            comparison_ids = [item.comparison_id for item in self._selected_comparisons()]
            program_comparison_ids = list(self.current_project.program_comparison_ids[:2])
        return AnalysisRequest(
            comparison_ids=comparison_ids,
            program_comparison_ids=program_comparison_ids,
            analysis_modules=list(self.current_project.selected_analysis_modules) if self.current_project is not None else list(self.ANALYSIS_MODULES),
            thresholds=self.current_thresholds,
        )

    def set_selected_comparisons(self, comparison_ids: list[str]) -> None:
        if self.current_project is None:
            raise RuntimeError("No project loaded.")
        self.current_project.selected_comparison_ids = list(dict.fromkeys(comparison_ids))
        all_available_ids = [item.comparison_id for item in self.current_project.available_comparisons]
        ordered = [item for item in self.current_project.comparison_order_ids if item in all_available_ids]
        for comparison_id in all_available_ids:
            if comparison_id not in ordered:
                ordered.append(comparison_id)
        self.current_project.comparison_order_ids = ordered
        if self.current_project.program_comparison_ids:
            selected = set(self.current_project.selected_comparison_ids)
            self.current_project.program_comparison_ids = [
                item for item in self.current_project.program_comparison_ids if item in selected
            ]
        self._ensure_default_comparison_pairs()
        self._ensure_default_visualization_groups()
        self.invalidate_comparison_sets_confirmation()
        self.invalidate_analysis_results(
            "Analysis results were cleared because the selected comparisons changed. Re-run analysis to refresh per-comparison results."
        )
        self.save_project_config()

    def move_comparison_up(self, comparison_id: str) -> None:
        if self.current_project is None:
            raise RuntimeError("No project loaded.")
        ordered = self._ordered_available_comparison_ids()
        try:
            index = ordered.index(comparison_id)
        except ValueError:
            return
        if index <= 0:
            return
        ordered[index - 1], ordered[index] = ordered[index], ordered[index - 1]
        self.current_project.comparison_order_ids = ordered
        self._apply_comparison_order_to_available()
        self.save_project_config()

    def move_comparison_down(self, comparison_id: str) -> None:
        if self.current_project is None:
            raise RuntimeError("No project loaded.")
        ordered = self._ordered_available_comparison_ids()
        try:
            index = ordered.index(comparison_id)
        except ValueError:
            return
        if index >= len(ordered) - 1:
            return
        ordered[index + 1], ordered[index] = ordered[index], ordered[index + 1]
        self.current_project.comparison_order_ids = ordered
        self._apply_comparison_order_to_available()
        self.save_project_config()

    def set_program_comparison_ids(self, comparison_ids: list[str]) -> None:
        if self.current_project is None:
            raise RuntimeError("No project loaded.")
        unique = list(dict.fromkeys(comparison_ids))
        self.current_project.program_comparison_ids = unique[:2]
        self.invalidate_comparison_sets_confirmation()
        self.save_project_config()

    def set_selected_analysis_modules(self, modules: list[str]) -> None:
        if self.current_project is None:
            raise RuntimeError("No project loaded.")
        normalized = [module for module in modules if module in self.ANALYSIS_MODULES]
        self.current_project.selected_analysis_modules = normalized or list(self.ANALYSIS_MODULES)
        self.save_project_config()

    def remove_comparison_from_selection(self, comparison_id: str) -> None:
        if self.current_project is None:
            raise RuntimeError("No project loaded.")
        self.current_project.removed_comparison_ids = list(
            dict.fromkeys([*self.current_project.removed_comparison_ids, comparison_id])
        )
        self.current_project.available_comparisons = [
            item for item in self.current_project.available_comparisons if item.comparison_id != comparison_id
        ]
        self.current_project.selected_comparison_ids = [
            item for item in self.current_project.selected_comparison_ids if item != comparison_id
        ]
        self.current_project.comparison_order_ids = [
            item for item in self.current_project.comparison_order_ids if item != comparison_id
        ]
        self.current_project.program_comparison_ids = [
            item for item in self.current_project.program_comparison_ids if item != comparison_id
        ]
        self.current_project.candidate_selection_genes.pop(comparison_id, None)
        self.current_project.candidate_blacklist_genes.pop(comparison_id, None)
        self.current_project.comparison_pairs = [
            pair
            for pair in self.current_project.comparison_pairs
            if pair.comparison_a != comparison_id and pair.comparison_b != comparison_id
        ]
        for group in self.current_project.visualization_groups:
            group.comparison_ids = [item for item in group.comparison_ids if item != comparison_id]
        self.current_project.visualization_groups = [
            group for group in self.current_project.visualization_groups if group.comparison_ids
        ]
        for sample in self.current_project.isoform_samples:
            if sample.comparison_id == comparison_id:
                sample.comparison_id = None
        self.invalidate_comparison_sets_confirmation()
        self.invalidate_analysis_results(
            "Analysis results were cleared because a comparison was removed from selection. Re-run analysis to refresh downstream views."
        )
        self.save_project_config()

    def add_visualization_group(self) -> None:
        if self.current_project is None:
            raise RuntimeError("No project loaded.")
        selected = self._selected_comparisons()
        default_ids = [item.comparison_id for item in selected[: min(4, len(selected))]]
        group_id = f"viz_{len(self.current_project.visualization_groups) + 1}"
        self.current_project.visualization_groups.append(
            VisualizationGroupDefinition(
                group_id=group_id,
                comparison_ids=default_ids,
                display_name=self._default_visualization_group_name(default_ids),
            )
        )
        self.invalidate_visualization_groups_confirmation()
        self.save_project_config()

    def remove_visualization_group(self, group_id: str) -> None:
        if self.current_project is None:
            raise RuntimeError("No project loaded.")
        self.current_project.visualization_groups = [
            group for group in self.current_project.visualization_groups if group.group_id != group_id
        ]
        self.invalidate_visualization_groups_confirmation()
        self.save_project_config()

    def update_visualization_group(
        self,
        group_id: str,
        *,
        display_name: str | None = None,
        comparison_ids: list[str] | None = None,
        enabled: bool | None = None,
    ) -> None:
        if self.current_project is None:
            raise RuntimeError("No project loaded.")
        group = next((item for item in self.current_project.visualization_groups if item.group_id == group_id), None)
        if group is None:
            raise KeyError(f"Visualization group not found: {group_id}")
        if display_name is not None:
            group.display_name = display_name or None
        if comparison_ids is not None:
            resolved = [self._resolve_comparison_reference(item) for item in comparison_ids]
            group.comparison_ids = [item for item in resolved if item]
            if not group.display_name:
                group.display_name = self._default_visualization_group_name(group.comparison_ids)
        if enabled is not None:
            group.enabled = enabled
        self.invalidate_visualization_groups_confirmation()
        self.save_project_config()

    def blacklist_genes_for_display(self) -> list[str]:
        if self.current_project is None:
            return []
        return list(self.current_project.blacklist_genes)

    def available_comparisons_for_display(self) -> list:
        if self.current_project is None:
            return []
        lookup = {item.comparison_id: item for item in self.current_project.available_comparisons}
        return [lookup[item] for item in self._ordered_available_comparison_ids() if item in lookup]

    def shortlist_genes_for_display(self) -> list[str]:
        if self.current_project is None:
            return []
        return list(self.current_project.shortlist_genes)

    def add_shortlist_gene(self, gene: str) -> None:
        if self.current_project is None:
            raise RuntimeError("No project loaded.")
        value = gene.strip()
        if not value:
            return
        current = list(self.current_project.shortlist_genes)
        if value not in current:
            current.append(value)
            self.current_project.shortlist_genes = current
            self.save_project_config()

    def remove_shortlist_gene(self, gene: str) -> None:
        if self.current_project is None:
            raise RuntimeError("No project loaded.")
        self.current_project.shortlist_genes = [
            item for item in self.current_project.shortlist_genes if item != gene
        ]
        self.save_project_config()

    def add_blacklist_gene(self, gene: str) -> None:
        if self.current_project is None:
            raise RuntimeError("No project loaded.")
        value = gene.strip()
        if not value:
            return
        current = list(self.current_project.blacklist_genes)
        if value not in current:
            current.append(value)
            self.current_project.blacklist_genes = current
            self.save_project_config()

    def remove_blacklist_gene(self, gene: str) -> None:
        if self.current_project is None:
            raise RuntimeError("No project loaded.")
        self.current_project.blacklist_genes = [
            item for item in self.current_project.blacklist_genes if item != gene
        ]
        self.save_project_config()

    def candidate_selection_genes_for_display(self, comparison_id: str) -> list[str]:
        if self.current_project is None or not comparison_id:
            return []
        return list(self.current_project.candidate_selection_genes.get(comparison_id, []))

    def candidate_blacklist_genes_for_display(self, comparison_id: str) -> list[str]:
        if self.current_project is None or not comparison_id:
            return []
        return list(self.current_project.candidate_blacklist_genes.get(comparison_id, []))

    def add_candidate_selection_gene(self, comparison_id: str, gene: str) -> None:
        if self.current_project is None:
            raise RuntimeError("No project loaded.")
        value = gene.strip()
        if not comparison_id or not value:
            return
        blocked = {
            item.strip()
            for item in self.current_project.candidate_blacklist_genes.get(comparison_id, [])
            if item and item.strip()
        }
        if value in blocked:
            return
        current = list(self.current_project.candidate_selection_genes.get(comparison_id, []))
        if value not in current:
            current.append(value)
            self.current_project.candidate_selection_genes[comparison_id] = current
            self.save_project_config()

    def remove_candidate_selection_gene(self, comparison_id: str, gene: str) -> None:
        if self.current_project is None or not comparison_id:
            raise RuntimeError("No project loaded.")
        self.current_project.candidate_selection_genes[comparison_id] = [
            item for item in self.current_project.candidate_selection_genes.get(comparison_id, []) if item != gene
        ]
        self.save_project_config()

    def add_candidate_blacklist_gene(self, comparison_id: str, gene: str) -> None:
        if self.current_project is None:
            raise RuntimeError("No project loaded.")
        value = gene.strip()
        if not comparison_id or not value:
            return
        current = list(self.current_project.candidate_blacklist_genes.get(comparison_id, []))
        if value not in current:
            current.append(value)
            self.current_project.candidate_blacklist_genes[comparison_id] = current
            selected = [
                item for item in self.current_project.candidate_selection_genes.get(comparison_id, [])
                if item != value
            ]
            self.current_project.candidate_selection_genes[comparison_id] = selected
            self.save_project_config()

    def remove_candidate_blacklist_gene(self, comparison_id: str, gene: str) -> None:
        if self.current_project is None or not comparison_id:
            raise RuntimeError("No project loaded.")
        self.current_project.candidate_blacklist_genes[comparison_id] = [
            item for item in self.current_project.candidate_blacklist_genes.get(comparison_id, []) if item != gene
        ]
        self.save_project_config()

    def reset_candidate_selection(self, comparison_id: str, *, top_n: int = 20) -> None:
        if self.current_project is None or not comparison_id:
            raise RuntimeError("No project loaded.")
        self.current_project.candidate_selection_genes[comparison_id] = self._default_candidate_selection_genes(
            comparison_id,
            top_n=top_n,
        )
        self.save_project_config()

    def apply_gene_blacklist(self, frame: pd.DataFrame) -> pd.DataFrame:
        if self.current_project is None or frame.empty or not self.current_project.blacklist_genes:
            return frame
        blocked = {item.strip() for item in self.current_project.blacklist_genes if item and item.strip()}
        if not blocked:
            return frame
        filtered = frame.copy()
        keep_mask = pd.Series(True, index=filtered.index)
        for column in ("geneSymbol", "gene_symbol", "GeneID", "gene_id"):
            if column in filtered.columns:
                keep_mask &= ~filtered[column].astype(str).isin(blocked)
        return filtered.loc[keep_mask].copy()

    def add_comparison_pair(self) -> None:
        if self.current_project is None:
            raise RuntimeError("No project loaded.")
        selected = self._selected_comparisons()
        first = selected[0] if len(selected) >= 1 else None
        second = selected[1] if len(selected) >= 2 else None
        pair_id = f"pair_{len(self.current_project.comparison_pairs) + 1}"
        self.current_project.comparison_pairs.append(
            ComparisonPairDefinition(
                pair_id=pair_id,
                comparison_a=first.comparison_id if first else None,
                comparison_b=second.comparison_id if second else None,
                experiment_group=first.experiment_group if first else None,
                control_group=second.experiment_group if second else None,
            )
        )
        self._autofill_comparison_pair(self.current_project.comparison_pairs[-1])
        self.invalidate_comparison_sets_confirmation()
        self.save_project_config()

    def remove_comparison_pair(self, pair_id: str) -> None:
        if self.current_project is None:
            raise RuntimeError("No project loaded.")
        self.current_project.comparison_pairs = [
            pair for pair in self.current_project.comparison_pairs if pair.pair_id != pair_id
        ]
        self.invalidate_comparison_sets_confirmation()
        self.save_project_config()

    def update_comparison_pair(
        self,
        pair_id: str,
        *,
        comparison_a: str | None = None,
        comparison_b: str | None = None,
        display_name: str | None = None,
        experiment_group: str | None = None,
        control_group: str | None = None,
        enabled: bool | None = None,
    ) -> None:
        if self.current_project is None:
            raise RuntimeError("No project loaded.")
        pair = next((item for item in self.current_project.comparison_pairs if item.pair_id == pair_id), None)
        if pair is None:
            raise KeyError(f"Comparison pair not found: {pair_id}")
        previous_display = self._default_pair_display_name(pair)
        previous_experiment = self._default_pair_experiment_group(pair)
        previous_control = self._default_pair_control_group(pair)
        autofill_display = pair.display_name in {None, "", previous_display}
        autofill_experiment = pair.experiment_group in {None, "", previous_experiment}
        autofill_control = pair.control_group in {None, "", previous_control}
        if comparison_a is not None:
            pair.comparison_a = self._resolve_comparison_reference(comparison_a)
        if comparison_b is not None:
            pair.comparison_b = self._resolve_comparison_reference(comparison_b)
        if display_name is not None:
            pair.display_name = display_name or None
            autofill_display = False
        if experiment_group is not None:
            pair.experiment_group = experiment_group or None
            autofill_experiment = False
        if control_group is not None:
            pair.control_group = control_group or None
            autofill_control = False
        if enabled is not None:
            pair.enabled = enabled
        self._autofill_comparison_pair(
            pair,
            update_display=autofill_display,
            update_experiment=autofill_experiment,
            update_control=autofill_control,
        )
        self.invalidate_comparison_sets_confirmation()
        self.save_project_config()

    def update_thresholds(self, **values: float) -> None:
        for key, value in values.items():
            if hasattr(self.current_thresholds, key):
                setattr(self.current_thresholds, key, value)
        self.save_project_config()

    def apply_threshold_preset(self, preset_name: str) -> None:
        preset = self.default_thresholds["presets"].get(preset_name)
        if preset is None:
            raise KeyError(f"Unknown preset: {preset_name}")
        self.current_thresholds = ThresholdConfig(**preset)
        self.save_project_config()

    def preview_splicing_landscape(self):
        if self.current_project is None:
            raise RuntimeError("No project loaded.")
        if self.run_state.status == "completed" and not self.run_state.results.splicing_landscape.empty:
            return self.run_state.results.splicing_landscape
        request = self.build_default_request()
        request.rmats_mode = self.current_project.selected_rmats_mode
        return self.splicing_landscape.run(self.current_project, request)

    def preview_program_comparison(self, allow_generate: bool = True):
        if self.current_project is None:
            raise RuntimeError("No project loaded.")
        if self.run_state.status == "completed":
            return self.run_state.results.program_events, self.run_state.results.program_summary
        if not allow_generate:
            return pd.DataFrame(), pd.DataFrame()
        request = self.build_default_request()
        request.rmats_mode = self.current_project.selected_rmats_mode
        return self.preview_program_comparison_uncached(request)

    def preview_mechanism_support(self, allow_generate: bool = True):
        if self.current_project is None:
            raise RuntimeError("No project loaded.")
        self._log_heavy("preview_mechanism_support triggered")
        if self.run_state.status == "completed" and not self.run_state.results.mechanism_support.empty:
            return self.run_state.results.mechanism_support
        if not allow_generate:
            return pd.DataFrame()
        request = self.build_default_request()
        request.rmats_mode = self.current_project.selected_rmats_mode
        mechanism = self.mechanism_support.run(self.current_project, request)
        self.run_state.results.mechanism_support = mechanism.copy()
        self._mark_module_state(
            "mechanism_support",
            "finished",
            "Mechanism support loaded for candidate screening support layers.",
            output_folder=str(self.current_project.output_root / "03_mechanism_support") if self.current_project.output_root else "",
            from_cache=False,
        )
        return mechanism

    def preview_tx_splicing_integration(self, allow_generate: bool = True):
        if self.current_project is None:
            raise RuntimeError("No project loaded.")
        self._log_heavy("preview_tx_splicing_integration triggered")
        if self.run_state.status == "completed":
            return (
                self._enrich_tx_splicing_table(self.run_state.results.tx_splicing_gene_table),
                self.run_state.results.tx_splicing_summary,
            )
        if not allow_generate:
            return pd.DataFrame(), pd.DataFrame()
        request = self.build_default_request()
        request.rmats_mode = self.current_project.selected_rmats_mode
        gene_table, summary = self.tx_splicing_integration.run(self.current_project, request)
        return self._enrich_tx_splicing_table(gene_table), summary

    def preview_candidate_gene_screening(self, allow_generate: bool = True, *, force_rebuild: bool = False) -> pd.DataFrame:
        if self.current_project is None:
            raise RuntimeError("No project loaded.")
        self._log_heavy("preview_candidate_gene_screening triggered")
        if not force_rebuild and self.run_state.status == "completed" and not self.run_state.results.candidate_gene_table.empty:
            if allow_generate and (self.run_state.results.mechanism_support is None or self.run_state.results.mechanism_support.empty):
                try:
                    self.preview_mechanism_support(allow_generate=True)
                except Exception:
                    pass
            cached = self.run_state.results.candidate_gene_table.copy()
            if not self._candidate_support_table_is_stale(cached):
                return cached
        if not allow_generate:
            return pd.DataFrame()
        integration = self.run_state.results.tx_splicing_gene_table
        if integration is None or integration.empty:
            self._mark_module_state(
                "candidate_ranking",
                "not_run",
                "Missing input: transcript + splicing integration not run for this project. Please run 5.2.1 first.",
            )
            return pd.DataFrame()
        mechanism = self.run_state.results.mechanism_support
        if mechanism is None or mechanism.empty:
            mechanism = self.preview_mechanism_support(allow_generate=True)
        rebuilt = self._build_candidate_gene_table(integration, mechanism).copy()
        self.run_state.results.candidate_gene_table = rebuilt.copy()
        self._mark_module_state(
            "candidate_ranking",
            "finished",
            "Candidate ranking built from cached transcript + splicing integration and current support tables.",
            from_cache=False,
        )
        return rebuilt

    def preview_cross_comparison_candidate_matrix(self, allow_generate: bool = True) -> pd.DataFrame:
        if self.current_project is None:
            raise RuntimeError("No project loaded.")
        cached = self.run_state.results.cross_comparison_candidate_matrix
        if cached is not None and not cached.empty:
            self._mark_module_state(
                "cross_comparison_candidates",
                "cache_hit",
                "Cross-comparison candidate matrix loaded from cache.",
                from_cache=True,
                output_folder=str(self.cross_comparison_output_dir() or ""),
            )
            return cached.copy()
        candidate_table = self.run_state.results.candidate_gene_table
        if candidate_table is None or candidate_table.empty:
            self._mark_module_state(
                "cross_comparison_candidates",
                "not_run",
                "Missing input: candidate ranking not run for this project. Please run 5.3.1 first.",
            )
            return pd.DataFrame()
        if not allow_generate:
            self._mark_module_state(
                "cross_comparison_candidates",
                "not_run",
                "Not run. Click 'Run Current Branch' to build the cross-comparison candidate matrix from cached 5.3.1 results.",
            )
            return pd.DataFrame()
        matrix = self._build_cross_comparison_candidate_matrix(candidate_table).copy()
        self.run_state.results.cross_comparison_candidate_matrix = matrix.copy()
        self._mark_module_state(
            "cross_comparison_candidates",
            "finished",
            "Cross-comparison candidate matrix built from cached per-comparison candidate results.",
            from_cache=False,
            output_folder=str(self.cross_comparison_output_dir() or ""),
        )
        self._write_cross_comparison_outputs()
        return matrix

    def preview_cross_comparison_significant_as_patterns(
        self,
        pair_id: str | None,
        *,
        abs_dpsi_cutoff: float = 0.10,
        large_delta_dpsi_cutoff: float = 0.10,
        allow_generate: bool = False,
        force_rebuild: bool = False,
    ) -> dict[str, pd.DataFrame]:
        if self.current_project is None:
            raise RuntimeError("No project loaded.")
        pair = next((item for item in self._active_comparison_pairs() if item.pair_id == pair_id), None)
        if pair is None:
            self._mark_module_state(
                "cross_comparison_as_patterns",
                "not_run",
                "Missing input: no enabled comparison pair is selected for 5.4.6.",
            )
            return self._empty_cross_as_pattern_tables()

        abs_cutoff = max(float(abs_dpsi_cutoff), 0.0)
        delta_cutoff = max(float(large_delta_dpsi_cutoff), 0.0)
        cache_key = (pair.pair_id, round(abs_cutoff, 6), round(delta_cutoff, 6))
        cached = None if force_rebuild else self._cross_as_pattern_cache.get(cache_key)
        if cached is not None:
            self._mark_module_state(
                "cross_comparison_as_patterns",
                "cache_hit",
                "Loaded cached 5.4.6 cross-comparison significant AS event dPSI pattern analysis.",
                from_cache=True,
                output_folder=str(self.cross_comparison_pattern_output_dir(pair.pair_id) or ""),
            )
            return {key: value.copy() for key, value in cached.items()}

        if not allow_generate:
            self._mark_module_state(
                "cross_comparison_as_patterns",
                "not_run",
                "Not run. Click 'Apply / Recalculate' to build the 5.4.6 cross-comparison significant AS event dPSI pattern tables for the selected comparison pair.",
            )
            return self._empty_cross_as_pattern_tables()

        tables = self._build_cross_comparison_significant_as_patterns(
            pair,
            abs_dpsi_cutoff=abs_cutoff,
            large_delta_dpsi_cutoff=delta_cutoff,
        )
        self._cross_as_pattern_cache[cache_key] = {key: value.copy() for key, value in tables.items()}
        self._mark_module_state(
            "cross_comparison_as_patterns",
            "finished",
            "Built 5.4.6 cross-comparison significant AS event dPSI pattern tables.",
            from_cache=False,
            output_folder=str(self.cross_comparison_pattern_output_dir(pair.pair_id) or ""),
        )
        self._write_cross_comparison_pattern_outputs(pair, tables, abs_cutoff, delta_cutoff)
        return {key: value.copy() for key, value in tables.items()}

    def run_main_analysis(self) -> AnalysisRunState:
        return self.run_main_analysis_with_progress()

    def run_main_analysis_with_progress(self, progress_callback=None) -> AnalysisRunState:
        if self.current_project is None:
            raise RuntimeError("No project loaded.")
        request = self.build_default_request()
        request.rmats_mode = self.current_project.selected_rmats_mode
        selected_modules = set(request.analysis_modules or self.ANALYSIS_MODULES)
        previous_results = self._copy_analysis_results()
        self.run_state = AnalysisRunState(
            request=request,
            status="running",
            message="Running main analysis modules...",
            results=previous_results,
        )
        try:
            if "splicing_landscape" in selected_modules:
                self._mark_module_state("splicing_landscape", "running", "Running splicing landscape...")
                if progress_callback is not None:
                    progress_callback("Running splicing landscape...")
                splicing_landscape = self.splicing_landscape.run(self.current_project, request)
                self.run_state.results.splicing_landscape = splicing_landscape
                self._mark_module_state(
                    "splicing_landscape",
                    "finished",
                    "Finished 5.1 Overview branch.",
                    output_folder=str(self.current_project.output_root / "01_splicing_landscape") if self.current_project.output_root else "",
                    from_cache=False,
                )
            if "mechanism_support" in selected_modules:
                self._mark_module_state("mechanism_support", "running", "Running mechanism support...")
                if progress_callback is not None:
                    progress_callback("Running mechanism support...")
                mechanism_support = self.mechanism_support.run(self.current_project, request)
                self.run_state.results.mechanism_support = mechanism_support
                self._mark_module_state(
                    "mechanism_support",
                    "finished",
                    "Finished mechanism support.",
                    output_folder=str(self.current_project.output_root / "03_mechanism_support") if self.current_project.output_root else "",
                    from_cache=False,
                )
            if "tx_splicing_integration" in selected_modules:
                self._mark_module_state("tx_splicing_integration", "running", "Running transcript + splicing integration...")
                if progress_callback is not None:
                    progress_callback("Running transcript-splicing integration...")
                tx_gene_table, tx_summary = self.tx_splicing_integration.run(self.current_project, request)
                self.run_state.results.tx_splicing_gene_table = self._enrich_tx_splicing_table(tx_gene_table)
                self.run_state.results.tx_splicing_summary = tx_summary
                self._mark_module_state(
                    "tx_splicing_integration",
                    "finished",
                    "Finished transcript + splicing integration.",
                    output_folder=str(self.current_project.output_root / "04_tx_splicing_integration") if self.current_project.output_root else "",
                    from_cache=False,
                )
            if "program_comparison" in selected_modules:
                self._mark_module_state("program_comparison", "running", "Running program comparison...")
                if progress_callback is not None:
                    progress_callback("Running program comparison...")
                program_events, program_summary = self.preview_program_comparison_uncached(request)
                self.run_state.results.program_events = program_events if program_events is not None else pd.DataFrame()
                self.run_state.results.program_summary = program_summary if program_summary is not None else pd.DataFrame()
                self._mark_module_state(
                    "program_comparison",
                    "finished",
                    "Finished cross-comparison program comparison.",
                    output_folder=str(self.current_project.output_root / "02_program_comparison") if self.current_project.output_root else "",
                    from_cache=False,
                )
            if progress_callback is not None:
                progress_callback("Writing selected analysis outputs...")
            self._write_analysis_outputs_for_modules(selected_modules)
            self.run_state.status = "completed"
            self.run_state.message = "Completed modules: " + ", ".join(request.analysis_modules or self.ANALYSIS_MODULES)
        except Exception as exc:
            self.run_state.status = "failed"
            self.run_state.message = str(exc)
            for module in selected_modules:
                state = self._module_states.get(module)
                if state and state.get("status") == "running":
                    self._mark_module_state(module, "failed", str(exc))
            raise
        return self.run_state

    def preview_program_comparison_uncached(self, request: AnalysisRequest):
        if self.current_project is None:
            raise RuntimeError("No project loaded.")
        pair_definitions = self._active_comparison_pairs()
        if not pair_definitions:
            comparisons = request.program_comparison_ids[:2] or request.comparison_ids[:2]
            if len(comparisons) < 2:
                return None, None
            return self.program_comparison.run(self.current_project, request, comparisons[0], comparisons[1])

        event_frames = []
        summary_frames = []
        for pair in pair_definitions:
            if not pair.comparison_a or not pair.comparison_b:
                continue
            events, summary = self.program_comparison.run(
                self.current_project,
                request,
                pair.comparison_a,
                pair.comparison_b,
            )
            if events is not None and not events.empty:
                events = events.copy()
                events["pair_id"] = pair.pair_id
                events["pair_name"] = pair.resolved_name
                events["pair_experiment_group"] = pair.experiment_group
                events["pair_control_group"] = pair.control_group
                event_frames.append(events)
            if summary is not None and not summary.empty:
                summary = summary.copy()
                summary["pair_id"] = pair.pair_id
                summary["pair_name"] = pair.resolved_name
                summary["pair_experiment_group"] = pair.experiment_group
                summary["pair_control_group"] = pair.control_group
                summary_frames.append(summary)
        if not event_frames and not summary_frames:
            return None, None
        return (
            pd.concat(event_frames, ignore_index=True) if event_frames else pd.DataFrame(),
            pd.concat(summary_frames, ignore_index=True) if summary_frames else pd.DataFrame(),
        )

    def export_frame(self, frame: pd.DataFrame, output_path: str | Path) -> Path:
        out = Path(output_path)
        out.parent.mkdir(parents=True, exist_ok=True)
        frame.to_csv(out, sep="\t", index=False)
        return out

    def update_rmats_mode(self, mode: str) -> None:
        if self.current_project is None:
            raise RuntimeError("No project loaded.")
        self.current_project.selected_rmats_mode = mode
        self.save_project_config()

    def update_output_root(self, output_root: str | Path) -> None:
        if self.current_project is None:
            raise RuntimeError("No project loaded.")
        self.current_project.output_root = Path(output_root).expanduser()
        self.save_project_config()

    def update_input_path(self, key: str, value: str | Path | None) -> None:
        if self.current_project is None:
            raise RuntimeError("No project loaded.")
        normalized = str(Path(value).expanduser()) if value else ""
        if normalized:
            self.current_project.input_paths[key] = normalized
        else:
            self.current_project.input_paths.pop(key, None)

        path_value = Path(normalized).expanduser().resolve() if normalized else None
        if key == "rmats_root":
            self.current_project.rmats_root = path_value
        elif key == "deg_root":
            self.current_project.deg_root = path_value
        elif key == "suppa_root":
            self.current_project.suppa_root = path_value
        elif key == "dexseq_root":
            self.current_project.dexseq_root = path_value
        elif key == "dtu_root":
            self.current_project.dtu_root = path_value
        elif key == "quant_root":
            self.current_project.quant_root = path_value
        elif key == "counts_path":
            self.current_project.counts_path = path_value
        elif key == "contrastsheet_path":
            self.current_project.contrastsheet_path = path_value
        self.save_project_config()

    def run_jutils_pipeline(self) -> str:
        if self.current_project is None:
            raise RuntimeError("No project loaded.")
        if self.jutils_adapter is None:
            raise RuntimeError("Jutils script not configured.")
        if self.current_project.output_root is None:
            raise RuntimeError("Output directory is not configured.")

        base = self.current_project.output_root / "08_jutils"
        converted_root = base / "converted"
        meta_root = base / "meta"
        plots_root = base / "plots"
        venn_root = base / "venn"
        logs = []
        command_lines: list[str] = []

        selected = self._selected_comparisons()
        comparison_names = []
        self._last_jutils_report = {
            "status": "running",
            "command_lines": [],
            "output_directory": str(base),
            "generated_files": [],
            "stdout": "",
            "stderr": "",
            "error_log": "",
        }
        try:
            for comparison in selected:
                if comparison.rmats_path is None:
                    continue
                comparison_names.append(comparison.resolved_name)
                comp_dir_name = self._safe_output_name(comparison.resolved_name)
                converted_dir = converted_root / comp_dir_name
                meta_file = meta_root / f"{comp_dir_name}.meta.tsv"
                plots_dir = plots_root / comp_dir_name
                rmats_post_dir = comparison.rmats_path / "rmats_post"
                control_label = comparison.control_group or "control"
                experiment_label = comparison.experiment_group or "experiment"

                convert_result = self.jutils_adapter.run_convert_results(rmats_post_dir, converted_dir)
                command_lines.append(self._command_to_text(convert_result.args))
                self.jutils_adapter.build_meta_template(rmats_post_dir, meta_file, control_label, experiment_label)
                tsv_file = converted_dir / "rmats_JC_results.tsv"
                heatmap_result = self.jutils_adapter.run_heatmap(tsv_file, meta_file, plots_dir, f"{comp_dir_name}_JC")
                command_lines.append(self._command_to_text(heatmap_result.args))
                pca_result = self.jutils_adapter.run_pca(tsv_file, meta_file, plots_dir, f"{comp_dir_name}_JC")
                command_lines.append(self._command_to_text(pca_result.args))
                logs.extend(
                    [
                        convert_result.stdout.strip(),
                        heatmap_result.stdout.strip(),
                        pca_result.stdout.strip(),
                    ]
                )

            tsv_list = self.jutils_adapter.build_tsv_list(converted_root, comparison_names, venn_root)
            venn_result = self.jutils_adapter.run_venn(tsv_list, venn_root, "katana_JC")
            command_lines.append(self._command_to_text(venn_result.args))
            logs.append(venn_result.stdout.strip())
            self.run_state.results.jutils_manifest = self.build_jutils_manifest_preview()
            self._write_step_output("08_jutils", {"jutils_manifest.tsv": self.run_state.results.jutils_manifest})
            self._save_jutils_summary_plot()
            browser = self.build_jutils_output_browser()
            self._last_jutils_report = {
                "status": "finished",
                "command_lines": command_lines,
                "output_directory": str(base),
                "generated_files": browser.to_dict("records"),
                "stdout": "\n".join(line for line in logs if line),
                "stderr": "",
                "error_log": "",
            }
            return "\n".join(line for line in logs if line)
        except Exception as exc:
            stdout = ""
            stderr = ""
            if hasattr(exc, "stdout") and getattr(exc, "stdout"):
                stdout = str(getattr(exc, "stdout"))
            if hasattr(exc, "stderr") and getattr(exc, "stderr"):
                stderr = str(getattr(exc, "stderr"))
            self._last_jutils_report = {
                "status": "failed",
                "command_lines": command_lines,
                "output_directory": str(base),
                "generated_files": self.build_jutils_output_browser().to_dict("records") if base.exists() else [],
                "stdout": stdout,
                "stderr": stderr,
                "error_log": str(exc),
            }
            raise

    def run_sashimi_pipeline(self, manifest: pd.DataFrame | None = None) -> str:
        if self.current_project is None:
            raise RuntimeError("No project loaded.")
        if self.sashimi_adapter is None:
            raise RuntimeError("rmats2sashimiplot script is not configured.")
        if self.current_project.output_root is None:
            raise RuntimeError("Output directory is not configured.")
        active_manifest = manifest.copy() if manifest is not None else self.run_state.results.sashimi_manifest.copy()
        if active_manifest.empty:
            raise RuntimeError("No selected sashimi manifest is available. Please select events and generate the manifest first.")
        logs, failures = self.sashimi_adapter.run_manifest(active_manifest, self.current_project.output_root)
        self._last_sashimi_failures = failures.copy()
        if failures.empty:
            self._mark_module_state(
                "sashimi_plot",
                "finished",
                "Finished rmats2sashimi for the selected events.",
                output_folder=str(self.current_project.output_root / "06_sashimi"),
                from_cache=False,
            )
        else:
            self._mark_module_state(
                "sashimi_plot",
                "failed",
                f"rmats2sashimi finished with {len(failures)} failed event(s).",
                output_folder=str(self.current_project.output_root / "06_sashimi"),
                from_cache=False,
            )
        return logs

    def run_isoform_switch_pipeline(self) -> str:
        if self.current_project is None:
            raise RuntimeError("No project loaded.")
        if self.isoform_runner is None:
            raise RuntimeError("IsoformSwitch runner is not configured.")
        if self.current_project.output_root is None:
            raise RuntimeError("Output directory is not configured.")
        gtf_path = self.current_project.tool_paths.get("gtf")
        if not gtf_path:
            raise RuntimeError("GTF is not configured. Set it in Settings.")
        manifest = self.run_state.results.isoform_manifest
        if manifest.empty:
            manifest = self.isoform_adapter.build_import_manifest_frame(self.selected_isoform_samples())
        if manifest.empty:
            raise RuntimeError("No quant.sf files detected.")
        if "condition" not in manifest.columns:
            raise RuntimeError("Isoform conditions are missing. Configure them in the Isoform page.")
        included_conditions = [
            str(value).strip()
            for value in manifest["condition"].fillna("").tolist()
        ]
        if any(not value for value in included_conditions):
            raise RuntimeError("All included isoform samples must have a condition before running.")
        if len(set(included_conditions)) < 2:
            raise RuntimeError("IsoformSwitchAnalyzeR requires at least two distinct conditions.")
        step_dir = self.current_project.output_root / "07_isoform_switch"
        data_dir = step_dir / "data"
        data_dir.mkdir(parents=True, exist_ok=True)
        manifest_file = data_dir / "isoform_quant_manifest.tsv"
        manifest.to_csv(manifest_file, sep="\t", index=False)
        design_file = self._write_isoform_design_template(data_dir, manifest)
        gtf_file = Path(gtf_path)
        fasta_file = Path(self.current_project.tool_paths["fasta"]) if self.current_project.tool_paths.get("fasta") else None
        script_path = self.isoform_runner.write_script(step_dir, manifest_file, design_file, gtf_file, fasta_file)
        result = self.isoform_runner.run(script_path, manifest_file, design_file, gtf_file, fasta_file, step_dir / "results")
        return result.stdout.strip()

    def run_cards_pipeline(self) -> str:
        if self.current_project is None:
            raise RuntimeError("No project loaded.")
        if self.current_project.output_root is None:
            raise RuntimeError("Output directory is not configured.")
        shortlist = self.run_state.results.cards_shortlist
        if shortlist.empty:
            raise RuntimeError("No shortlist available. Run analysis and generate shortlist first.")
        expr = self.run_state.results.cards_expression_support
        if expr.empty:
            expr = self.load_first_matching_table("expression_support_group_means.tsv")
        output_dir = self.current_project.output_root / "05_cards" / "cards_fixed"
        generated = self.cards_adapter.render_cards(shortlist, expr, output_dir)
        return "\n".join(str(path) for path in generated)

    def load_first_matching_table(self, pattern: str) -> pd.DataFrame:
        if self.current_project is None:
            raise RuntimeError("No project loaded.")
        matches = list(self.current_project.project_root.rglob(pattern))
        if not matches:
            return pd.DataFrame()
        return pd.read_csv(matches[0], sep="\t", low_memory=False)

    def load_all_matching_tables(self, pattern: str) -> pd.DataFrame:
        if self.current_project is None:
            raise RuntimeError("No project loaded.")
        frames = []
        for path in self.current_project.project_root.rglob(pattern):
            try:
                frames.append(pd.read_csv(path, sep="\t", low_memory=False))
            except Exception:
                continue
        return pd.concat(frames, axis=0, ignore_index=True) if frames else pd.DataFrame()

    def runtime_status_lines(self) -> list[str]:
        lines = []
        pysam_available = importlib.util.find_spec("pysam") is not None
        lines.append(f"pysam: {'available' if pysam_available else 'missing'}")
        lines.append(f"Rscript: {self.current_project.tool_paths.get('rscript', 'not set') if self.current_project else 'not set'}")
        lines.append(f"GTF: {self.current_project.tool_paths.get('gtf', 'not set') if self.current_project else 'not set'}")
        lines.append(f"Annotation TSV: {self.current_project.tool_paths.get('annotation_tsv', 'not set') if self.current_project else 'not set'}")
        lines.append(f"FASTA: {self.current_project.tool_paths.get('fasta', 'not set') if self.current_project else 'not set'}")
        lines.append(f"BAM root: {self.current_project.tool_paths.get('bam_root', 'not set') if self.current_project else 'not set'}")
        lines.append(f"Jutils adapter: {'ready' if self.jutils_adapter is not None else 'not ready'}")
        lines.append(f"rmats2sashimiplot adapter: {'ready' if self.sashimi_adapter is not None else 'not ready'}")
        lines.append(f"IsoformSwitch runner: {'ready' if self.isoform_runner is not None and self.isoform_runner.rscript_path is not None else 'not ready'}")
        return lines

    def build_input_confirmation_report(self) -> str:
        if self.current_project is None:
            return "No project loaded."

        project = self.current_project
        lines = [
            f"Project root: {project.project_root}",
            f"Input mode: {'manual input mode' if project.input_paths else 'auto scan mode'}",
            "",
        ]

        def add_dir_summary(label: str, path: Path | None, patterns: list[str] | None = None, dirs_only: bool = False) -> None:
            lines.append(f"{label}: {path if path else 'not set'}")
            if path is None:
                lines.append("  Status: missing")
                lines.append("")
                return
            if not path.exists():
                lines.append("  Status: path does not exist")
                lines.append("")
                return
            if path.is_file():
                lines.append("  Status: file exists")
                lines.append("")
                return
            lines.append("  Status: directory exists")
            if dirs_only:
                subdirs = sorted([item.name for item in path.iterdir() if item.is_dir()])
                lines.append(f"  Subdirectories: {len(subdirs)}")
                if subdirs:
                    lines.append(f"  Examples: {', '.join(subdirs[:5])}")
            elif patterns:
                for pattern in patterns:
                    matches = list(path.rglob(pattern))
                    lines.append(f"  {pattern}: {len(matches)}")
                    if matches:
                        preview = ", ".join(match.name for match in matches[:3])
                        lines.append(f"    Examples: {preview}")
            lines.append("")

        add_dir_summary("rMATS root", project.rmats_root, dirs_only=True)
        add_dir_summary("DEG root", project.deg_root, patterns=["*.deseq2.results.tsv"])
        add_dir_summary("SUPPA root", project.suppa_root, patterns=["*_local_diffsplice.dpsi"])
        add_dir_summary("DEXSeq root", project.dexseq_root, patterns=["perGeneQValue.*.csv"])
        add_dir_summary("DTU root", project.dtu_root, patterns=["DEXSeqResults.*.tsv", "getAdjustedPValues.*.tsv"])
        add_dir_summary("quant root", project.quant_root, patterns=["quant.sf"])
        add_dir_summary("counts file", project.counts_path)
        add_dir_summary("contrastsheet", project.contrastsheet_path)
        add_dir_summary("annotation TSV", Path(project.tool_paths["annotation_tsv"]) if project.tool_paths.get("annotation_tsv") else None)

        lines.extend(
            [
                f"Detected comparisons: {len(project.available_comparisons)}",
                f"Selected comparisons: {', '.join(project.selected_comparison_ids) if project.selected_comparison_ids else 'none'}",
                f"Program comparison pair: {', '.join(project.program_comparison_ids) if project.program_comparison_ids else 'auto/first two selected'}",
                f"Defined comparison pairs: {len(project.comparison_pairs)}",
                f"Detected isoform samples: {len(project.isoform_samples)}",
                "",
                "Resolved comparison file mapping:",
                "",
            ]
        )
        for comparison in project.available_comparisons:
            lines.extend(
                [
                    f"[{comparison.comparison_id}]",
                    f"Display: {comparison.resolved_name}",
                    f"rMATS source: {comparison.rmats_name or comparison.comparison_id}",
                    f"DEG source: {comparison.deg_name or comparison.comparison_id}",
                    f"rMATS file/dir: {comparison.rmats_file_label}",
                    f"DEG file: {comparison.deg_file_label}",
                    f"quant: {'yes' if comparison.has_quant else 'no'} ({comparison.quant_dir if comparison.quant_dir else 'missing'})",
                    "",
                ]
            )
        return "\n".join(lines)

    def build_selected_input_confirmation_report(self, comparison_id: str | None = None) -> str:
        if self.current_project is None:
            return "No project loaded."

        project = self.current_project
        selected = self.selected_comparisons_for_display()
        if comparison_id:
            selected = [item for item in selected if item.comparison_id == comparison_id]

        lines = [
            f"Project root: {project.project_root}",
            f"Input mode: {'manual input mode' if project.input_paths else 'auto scan mode'}",
            f"Selected comparisons in scope: {len(selected)}",
            "",
            f"rMATS root      : {project.rmats_root or 'auto'}",
            f"DEG root        : {project.deg_root or 'auto'}",
            f"SUPPA root      : {project.suppa_root or 'auto'}",
            f"DEXSeq root     : {project.dexseq_root or 'auto'}",
            f"DTU root        : {project.dtu_root or 'auto'}",
            f"quant root      : {project.quant_root or 'auto'}",
            f"counts file     : {project.counts_path or 'auto'}",
            f"contrastsheet   : {project.contrastsheet_path or 'auto'}",
            f"annotation TSV  : {project.tool_paths.get('annotation_tsv', 'auto')}",
            "",
        ]

        if not selected:
            lines.append("No selected comparisons available.")
            return "\n".join(lines)

        for index, comparison in enumerate(selected, start=1):
            quant_files = self._comparison_quant_files(comparison)
            suppa_file = self._resolve_comparison_support_file(comparison, "suppa")
            dexseq_file = self._resolve_comparison_support_file(comparison, "dexseq")
            dtu_file = self._resolve_comparison_support_file(comparison, "dtu")
            stage_r_file = self._resolve_comparison_support_file(comparison, "stager")
            lines.extend(
                [
                    "=" * 88,
                    f"[{index}] {comparison.comparison_id}",
                    f"Display name          : {comparison.resolved_name}",
                    f"Analysis experiment   : {comparison.experiment_group or 'missing'}",
                    f"Analysis control      : {comparison.control_group or 'missing'}",
                    f"Read AS experiment    : {comparison.display_source_groups('rmats')[0] or 'missing'}",
                    f"Read AS control       : {comparison.display_source_groups('rmats')[1] or 'missing'}",
                    f"Read DEG experiment   : {comparison.display_source_groups('deg')[0] or 'missing'}",
                    f"Read DEG control      : {comparison.display_source_groups('deg')[1] or 'missing'}",
                    "",
                    f"rMATS source          : {comparison.display_source_name('rmats')}",
                    f"rMATS directory       : {comparison.rmats_file_label}",
                    f"DEG source            : {comparison.display_source_name('deg')}",
                    f"DEG file              : {comparison.deg_file_label}",
                    f"SUPPA source          : {comparison.source_name('suppa')}",
                    f"SUPPA file            : {suppa_file or 'missing'}",
                    f"DEXSeq source         : {comparison.source_name('dexseq')}",
                    f"DEXSeq file           : {dexseq_file or 'missing'}",
                    f"DTU source            : {comparison.source_name('dtu')}",
                    f"DTU file              : {dtu_file or 'missing'}",
                    f"stageR file           : {stage_r_file or 'missing'}",
                    f"Quant source          : {comparison.source_name('quant')}",
                    f"Quant group dir       : {comparison.quant_dir or 'missing'}",
                    f"Quant files detected  : {len(quant_files)}",
                ]
            )
            if quant_files:
                lines.extend(f"  - {path}" for path in quant_files)
            else:
                lines.append("  - missing")
            lines.extend(
                [
                    f"Counts file           : {project.counts_path or 'missing'}",
                    f"Contrastsheet         : {project.contrastsheet_path or 'missing'}",
                    f"Sashimi source        : {comparison.source_name('sashimi')}",
                    f"BAM group             : {comparison.source_name('bam')}",
                    f"Isoform ready         : {'yes' if comparison.has_quant else 'no'}",
                    "",
                ]
            )

        return "\n".join(lines)

    def build_jutils_manifest_preview(self) -> pd.DataFrame:
        if self.current_project is None:
            raise RuntimeError("No project loaded.")
        output_root = self.current_project.output_root or (self.current_project.project_root / "katana_output")
        converted_root = output_root / "08_jutils" / "converted"
        comparison_names = [comparison.resolved_name for comparison in self._selected_comparisons()]
        return self.jutils_adapter.build_tsv_list_frame(comparison_names, converted_root)

    def jutils_report(self) -> dict[str, object]:
        report = dict(self._last_jutils_report)
        report.setdefault("status", "idle")
        report.setdefault("command_lines", [])
        report.setdefault("output_directory", "")
        report.setdefault("generated_files", [])
        report.setdefault("stdout", "")
        report.setdefault("stderr", "")
        report.setdefault("error_log", "")
        return report

    def build_jutils_output_browser(self) -> pd.DataFrame:
        if self.current_project is None:
            raise RuntimeError("No project loaded.")
        output_root = self.current_project.output_root or (self.current_project.project_root / "katana_output")
        base = output_root / "08_jutils"
        if not base.exists():
            return pd.DataFrame(
                columns=[
                    "comparison_id",
                    "relative_path",
                    "absolute_path",
                    "file_type",
                    "file_size",
                    "last_modified",
                    "preview_kind",
                    "previewable",
                ]
            )
        safe_name_to_comparison = {
            self._safe_output_name(comparison.resolved_name): comparison.comparison_id
            for comparison in self.current_project.available_comparisons
        }
        rows: list[dict[str, object]] = []
        for path in sorted(base.rglob("*")):
            if not path.is_file():
                continue
            rel = path.relative_to(base)
            comparison_id = ""
            parts = rel.parts
            if len(parts) >= 2:
                token = Path(parts[1]).stem if parts[0] == "meta" else parts[1]
                comparison_id = safe_name_to_comparison.get(token, "")
            stat = path.stat()
            suffix = path.suffix.lower()
            preview_kind = self._preview_kind_for_path(path)
            rows.append(
                {
                    "comparison_id": comparison_id,
                    "relative_path": str(rel),
                    "absolute_path": str(path),
                    "file_type": suffix.lstrip(".") or "file",
                    "file_size": stat.st_size,
                    "last_modified": datetime.fromtimestamp(stat.st_mtime).strftime("%Y-%m-%d %H:%M:%S"),
                    "preview_kind": preview_kind,
                    "previewable": preview_kind != "binary",
                }
            )
        return pd.DataFrame(rows)

    def candidate_screening_output_dir(self, comparison_id: str | None) -> Path | None:
        if self.current_project is None or self.current_project.output_root is None or not comparison_id:
            return None
        comparison = next(
            (item for item in self.current_project.available_comparisons if item.comparison_id == comparison_id),
            None,
        )
        if comparison is None:
            return None
        return self.current_project.output_root / "02_candidate_gene_screening" / (comparison.output_prefix or self._safe_output_name(comparison.comparison_id))

    def candidate_screening_file_map(self, comparison_id: str | None) -> dict[str, Path]:
        base = self.candidate_screening_output_dir(comparison_id)
        if base is None:
            return {}
        return {
            "gene_level_integrated_candidates": base / "gene_level_integrated_candidates.tsv",
            "tier1": base / "tier_1_candidates.tsv",
            "tier2": base / "tier_2_candidates.tsv",
            "tier3": base / "tier_3_candidates.tsv",
            "tier4": base / "tier_4_candidates.tsv",
            "tier5": base / "tier_5_candidates.tsv",
            "deg": base / "DEG_significant_genes.tsv",
            "rmats_gene": base / "rMATS_gene_summary.tsv",
            "dexseq_gene": base / "DEXSeq_gene_summary.tsv",
            "dtu_gene": base / "DTU_gene_summary.tsv",
            "shortlist": base / "shortlist.tsv",
            "blacklist": base / "blacklist.tsv",
            "summary": base / "summary.tsv",
        }

    def cross_comparison_output_dir(self) -> Path | None:
        if self.current_project is None or self.current_project.output_root is None:
            return None
        return self.current_project.output_root / "03_cross_comparison_candidate_comparison"

    def cross_comparison_file_map(self) -> dict[str, Path]:
        base = self.cross_comparison_output_dir()
        if base is None:
            return {}
        return {
            "matrix": base / "cross_comparison_candidate_matrix.tsv",
            "summary": base / "cross_comparison_summary.tsv",
        }

    def cross_comparison_pattern_output_dir(self, pair_id: str | None) -> Path | None:
        base = self.cross_comparison_output_dir()
        if base is None or not pair_id:
            return None
        return base / self._safe_output_name(pair_id)

    def cross_comparison_pattern_file_map(self, pair_id: str | None) -> dict[str, Path]:
        base = self.cross_comparison_pattern_output_dir(pair_id)
        if base is None:
            return {}
        return {
            "same_direction_large_delta": base / "same_direction_large_delta_significant_as_events.tsv",
            "opposite_direction": base / "opposite_direction_significant_as_events.tsv",
            "all_shared": base / "all_shared_significant_as_pairwise_dpsi_delta.tsv",
            "condition_specific": base / "optional_condition_specific_as_events.tsv",
            "summary": base / "summary_counts.tsv",
        }

    def _write_cross_comparison_pattern_outputs(
        self,
        pair: ComparisonPairDefinition,
        tables: dict[str, pd.DataFrame],
        abs_cutoff: float,
        delta_cutoff: float,
    ) -> None:
        base = self.cross_comparison_pattern_output_dir(pair.pair_id)
        if base is None:
            return
        base.mkdir(parents=True, exist_ok=True)
        file_map = self.cross_comparison_pattern_file_map(pair.pair_id)
        for key in ("same_direction_large_delta", "opposite_direction", "all_shared", "condition_specific"):
            frame = tables.get(key, pd.DataFrame())
            path = file_map.get(key)
            if path is None:
                continue
            frame.to_csv(path, sep="\t", index=False)
        summary = pd.DataFrame(
            [
                {
                    "pair_id": pair.pair_id,
                    "pair_name": pair.resolved_name,
                    "comparison_A": pair.comparison_a,
                    "comparison_B": pair.comparison_b,
                    "n_same_direction_large_delta": len(tables.get("same_direction_large_delta", pd.DataFrame())),
                    "n_opposite_direction": len(tables.get("opposite_direction", pd.DataFrame())),
                    "n_all_shared": len(tables.get("all_shared", pd.DataFrame())),
                    "n_condition_specific": len(tables.get("condition_specific", pd.DataFrame())),
                    "abs_dPSI_significance_cutoff_used": abs_cutoff,
                    "large_delta_dPSI_cutoff_used": delta_cutoff,
                }
            ]
        )
        summary_path = file_map.get("summary")
        if summary_path is not None:
            summary.to_csv(summary_path, sep="\t", index=False)

    def selected_isoform_samples(self) -> list[IsoformSampleDefinition]:
        if self.current_project is None:
            return []
        selected = self._selected_comparisons()
        selected_ids = {item.comparison_id for item in selected}
        selected_groups = {
            group
            for comparison in selected
            for group in (
                comparison.source_experiment_group or comparison.experiment_group,
                comparison.source_control_group or comparison.control_group,
            )
            if group
        }
        samples = [
            sample
            for sample in self.current_project.isoform_samples
            if sample.include
            and sample.quant_sf.exists()
            and (
                (sample.sample_group in selected_groups if sample.sample_group else False)
                or (sample.comparison_id in selected_ids if sample.comparison_id else False)
                or (not sample.sample_group and not sample.comparison_id)
            )
        ]
        return samples

    def selected_comparisons_for_display(self) -> list:
        return list(self._selected_comparisons())

    def comparison_context(self, comparison_id: str | None) -> dict[str, str]:
        if self.current_project is None or not comparison_id:
            return {
                "comparison_id": "",
                "display_name": "",
                "group1": "group1",
                "group2": "group2",
                "log2fc_direction_label": "",
                "dpsi_direction_label": "",
                "rmats_canonical_direction": "",
                "dexseq_inherited_direction": "",
                "dtu_inherited_direction": "",
                "deg_final_direction": "",
                "output_folder": "",
            }
        comparison = next(
            (item for item in self.current_project.available_comparisons if item.comparison_id == comparison_id),
            None,
        )
        if comparison is None:
            return {
                "comparison_id": comparison_id,
                "display_name": comparison_id,
                "group1": "group1",
                "group2": "group2",
                "log2fc_direction_label": "",
                "dpsi_direction_label": "",
                "rmats_canonical_direction": "",
                "dexseq_inherited_direction": "",
                "dtu_inherited_direction": "",
                "deg_final_direction": "",
                "output_folder": "",
            }
        group1, group2 = self._comparison_group_labels(comparison)
        rmats_group1, rmats_group2 = self._comparison_splicing_direction_labels(comparison)
        deg_group1, deg_group2 = self._comparison_deg_direction_labels(comparison)
        output_folder = ""
        if self.current_project.output_root is not None:
            output_folder = str(
                self.current_project.output_root / "02_candidate_gene_screening" / (comparison.output_prefix or self._safe_output_name(comparison.comparison_id))
            )
        return {
            "comparison_id": comparison.comparison_id,
            "display_name": comparison.display_resolved_name,
            "group1": group1,
            "group2": group2,
            "log2fc_direction_label": f"log2FC > 0 means transcript expression is higher in {deg_group1} than {deg_group2}",
            "dpsi_direction_label": f"dPSI > 0 means event usage is higher in {rmats_group1} than {rmats_group2}",
            "rmats_canonical_direction": f"rMATS canonical direction: {rmats_group1} vs {rmats_group2}",
            "dexseq_inherited_direction": f"DEXSeq inherited direction: {rmats_group1} vs {rmats_group2}",
            "dtu_inherited_direction": f"DTU inherited direction: {rmats_group1} vs {rmats_group2}",
            "deg_final_direction": f"DEG selected/final direction: {deg_group1} vs {deg_group2}",
            "output_folder": output_folder,
        }

    def comparison_context_text(self, comparison_id: str | None) -> str:
        context = self.comparison_context(comparison_id)
        if not context["comparison_id"]:
            return "No comparison selected."
        return (
            f"{context['display_name']} [{context['comparison_id']}] | "
            f"{context['group1']} vs {context['group2']} | "
            f"{context['log2fc_direction_label']} | "
            f"{context['dpsi_direction_label']} | "
            f"{context['rmats_canonical_direction']} | "
            f"{context['dexseq_inherited_direction']} | "
            f"{context['dtu_inherited_direction']} | "
            f"{context['deg_final_direction']}"
        )

    def comparison_header_text(self, comparison_id: str | None, *, cache_status: str | None = None) -> str:
        if self.current_project is None or not comparison_id:
            return "No comparison selected."
        comparison = next(
            (item for item in self.current_project.available_comparisons if item.comparison_id == comparison_id),
            None,
        )
        if comparison is None:
            return self.comparison_context_text(comparison_id)
        context = self.comparison_context(comparison_id)
        thresholds = (
            f"DE_FDR<={self.current_thresholds.deg_padj:g}, "
            f"|log2FC|>={self.current_thresholds.deg_log2fc:g}, "
            f"rMATS_FDR<={self.current_thresholds.splicing_fdr:g}, "
            f"|dPSI|>={self.current_thresholds.splicing_dpsi:g}, "
            f"DEXSeq_q<={self.current_thresholds.dexseq_qvalue:g}, "
            f"DTU_q<={self.current_thresholds.dtu_qvalue:g}"
        )
        resolved_cache = cache_status or ("cache hit" if self.run_state.status == "completed" else "preview")
        lines = [
            f"{context['display_name']} [{context['comparison_id']}]",
            f"Biological question: {comparison.biological_question or 'not provided'}",
            f"Groups: {context['group1']} vs {context['group2']}",
            context["log2fc_direction_label"],
            context["dpsi_direction_label"],
            context["rmats_canonical_direction"],
            context["dexseq_inherited_direction"],
            context["dtu_inherited_direction"],
            context["deg_final_direction"],
            f"Thresholds: {thresholds}",
            f"Cache status: {resolved_cache}",
            f"Output folder: {context['output_folder'] or 'not available yet'}",
        ]
        return "\n".join(lines)

    def _enrich_tx_splicing_table(self, frame: pd.DataFrame) -> pd.DataFrame:
        if self.current_project is None or frame is None or frame.empty:
            return frame if frame is not None else pd.DataFrame()
        working = frame.copy()
        comparison_lookup = {
            item.comparison_id: item for item in self.current_project.available_comparisons
        }
        annotation_lookup = self._annotation_gene_symbol_lookup()
        normalized_gene_ids = self._normalized_gene_id_series(working)
        if "gene_symbol" not in working.columns:
            base_symbol = working.get("geneSymbol")
            if base_symbol is None:
                base_symbol = pd.Series(pd.NA, index=working.index, dtype="object")
            working["gene_symbol"] = base_symbol
        if annotation_lookup:
            mapped_symbols = normalized_gene_ids.map(annotation_lookup)
            working["gene_symbol"] = mapped_symbols.combine_first(working["gene_symbol"])
        working["gene_symbol"] = working["gene_symbol"].combine_first(working.get("gene_id"))
        if "standardized_log2FC" not in working.columns:
            working["standardized_log2FC"] = pd.to_numeric(working.get("log2FC"), errors="coerce")
        if "standardized_dPSI" not in working.columns:
            working["standardized_dPSI"] = pd.to_numeric(working.get("representative_dPSI"), errors="coerce")
        if "DE_FDR" not in working.columns:
            working["DE_FDR"] = pd.to_numeric(working.get("deg_padj"), errors="coerce")
        if "rMATS_FDR" not in working.columns:
            working["rMATS_FDR"] = pd.to_numeric(working.get("representative_event_FDR"), errors="coerce")
        working["group1"] = working["comparison_id"].map(
            lambda value: self._comparison_group_labels(comparison_lookup.get(str(value)))[0]
        )
        working["group2"] = working["comparison_id"].map(
            lambda value: self._comparison_group_labels(comparison_lookup.get(str(value)))[1]
        )
        working["rmats_group1"] = working["comparison_id"].map(
            lambda value: self._comparison_splicing_direction_labels(comparison_lookup.get(str(value)))[0]
        )
        working["rmats_group2"] = working["comparison_id"].map(
            lambda value: self._comparison_splicing_direction_labels(comparison_lookup.get(str(value)))[1]
        )
        working["deg_group1"] = working["comparison_id"].map(
            lambda value: self._comparison_deg_direction_labels(comparison_lookup.get(str(value)))[0]
        )
        working["deg_group2"] = working["comparison_id"].map(
            lambda value: self._comparison_deg_direction_labels(comparison_lookup.get(str(value)))[1]
        )
        working["log2fc_direction_label"] = working.apply(
            lambda row: (
                f"log2FC > 0 means transcript expression is higher in {row['deg_group1']} than {row['deg_group2']}"
                if pd.notna(row.get("deg_group1")) and pd.notna(row.get("deg_group2"))
                else "log2FC > 0 follows the configured transcript numerator group"
            ),
            axis=1,
        )
        working["dpsi_direction_label"] = working.apply(
            lambda row: (
                f"dPSI > 0 means event usage is higher in {row['rmats_group1']} than {row['rmats_group2']}"
                if pd.notna(row.get("rmats_group1")) and pd.notna(row.get("rmats_group2"))
                else "dPSI > 0 follows the configured splicing group1"
            ),
            axis=1,
        )
        working["transcript_direction_flipped"] = working["comparison_id"].map(
            lambda value: bool(comparison_lookup.get(str(value)).reverse_deg) if str(value) in comparison_lookup else False
        )
        working["splicing_direction_flipped"] = working["comparison_id"].map(
            lambda value: bool(comparison_lookup.get(str(value)).reverse_splicing) if str(value) in comparison_lookup else False
        )
        return working

    def _annotation_gene_symbol_lookup(self) -> dict[str, str]:
        if self.current_project is None:
            return {}
        return load_annotation_gene_symbol_lookup(self.current_project.tool_paths.get("annotation_tsv"))

    def _normalized_gene_id_series(self, frame: pd.DataFrame) -> pd.Series:
        if "gene_id" in frame.columns:
            return frame["gene_id"].map(normalize_gene_id)
        if "GeneID" in frame.columns:
            return frame["GeneID"].map(normalize_gene_id)
        if "groupID" in frame.columns:
            return frame["groupID"].map(normalize_gene_id)
        return pd.Series(pd.NA, index=frame.index, dtype="object")

    @staticmethod
    def _gene_lookup_key(gene_id: object, gene_symbol: object) -> str:
        normalized = normalize_gene_id(gene_id)
        if normalized:
            return f"id:{normalized}"
        symbol = str(gene_symbol or "").strip()
        if symbol:
            return f"sym:{symbol.casefold()}"
        return ""

    @staticmethod
    def _parse_inc_level_mean(value: object) -> float:
        if pd.isna(value):
            return float("nan")
        parts: list[float] = []
        for token in str(value).split(","):
            text = token.strip()
            if not text or text.upper() == "NA":
                continue
            try:
                parts.append(float(text))
            except ValueError:
                continue
        if not parts:
            return float("nan")
        return float(sum(parts) / len(parts))

    def _expression_support_frame(self) -> pd.DataFrame:
        if self._expression_support_loaded:
            return self._expression_support_cache.copy()
        self._expression_support_loaded = True
        if self.current_project is None:
            self._expression_support_cache = pd.DataFrame()
            return pd.DataFrame()
        frame = self.load_first_matching_table("expression_support_group_means.tsv")
        if frame is None or frame.empty:
            self._expression_support_cache = pd.DataFrame()
            return pd.DataFrame()
        working = frame.copy()
        working["normalized_gene_id"] = self._normalized_gene_id_series(working)
        if "gene_symbol" not in working.columns:
            base_symbol = working.get("geneSymbol")
            if base_symbol is None:
                base_symbol = pd.Series(pd.NA, index=working.index, dtype="object")
            working["gene_symbol"] = base_symbol
        annotation_lookup = self._annotation_gene_symbol_lookup()
        if annotation_lookup:
            mapped = working["normalized_gene_id"].map(annotation_lookup)
            working["gene_symbol"] = mapped.combine_first(working["gene_symbol"])
        if "expr" in working.columns:
            working["expr"] = pd.to_numeric(working["expr"], errors="coerce")
        working["gene_card_key"] = working.apply(
            lambda row: self._gene_lookup_key(row.get("normalized_gene_id") or row.get("gene_id"), row.get("gene_symbol")),
            axis=1,
        )
        self._expression_support_cache = working
        return working.copy()

    def _deg_expression_frame(self, comparison_id: str | None) -> pd.DataFrame:
        if self.current_project is None or not comparison_id:
            return pd.DataFrame()
        cached = self._deg_expression_cache.get(comparison_id)
        if cached is not None:
            return cached.copy()
        comparison = next(
            (item for item in self.current_project.available_comparisons if item.comparison_id == comparison_id),
            None,
        )
        if comparison is None or comparison.deg_path is None or not comparison.deg_path.exists():
            self._deg_expression_cache[comparison_id] = pd.DataFrame()
            return pd.DataFrame()
        try:
            frame = pd.read_csv(comparison.deg_path, sep="\t", low_memory=False)
        except Exception:
            self._deg_expression_cache[comparison_id] = pd.DataFrame()
            return pd.DataFrame()
        if frame.empty:
            self._deg_expression_cache[comparison_id] = pd.DataFrame()
            return pd.DataFrame()
        working = frame.copy()
        working["gene_id"] = self._normalized_gene_id_series(working)
        if "gene_symbol" not in working.columns:
            base_symbol = working.get("geneSymbol")
            if base_symbol is None:
                base_symbol = pd.Series(pd.NA, index=working.index, dtype="object")
            working["gene_symbol"] = base_symbol
        annotation_lookup = self._annotation_gene_symbol_lookup()
        if annotation_lookup:
            working["gene_symbol"] = working["gene_id"].map(annotation_lookup).combine_first(working["gene_symbol"])
        working["gene_symbol"] = working["gene_symbol"].combine_first(working["gene_id"])
        working["baseMean"] = pd.to_numeric(working.get("baseMean"), errors="coerce")
        working["standardized_log2FC"] = pd.to_numeric(working.get("log2FoldChange"), errors="coerce")
        if comparison.reverse_deg:
            working["standardized_log2FC"] = -working["standardized_log2FC"]
        working["DE_FDR"] = pd.to_numeric(working.get("padj"), errors="coerce")
        working["DEG_significant"] = (
            (working["DE_FDR"] <= self.current_thresholds.deg_padj)
            & (working["standardized_log2FC"].abs() >= self.current_thresholds.deg_log2fc)
        )
        working["gene_card_key"] = working.apply(
            lambda row: self._gene_lookup_key(row.get("gene_id"), row.get("gene_symbol")),
            axis=1,
        )
        result = working[
            [
                column
                for column in [
                    "gene_card_key",
                    "gene_id",
                    "gene_symbol",
                    "baseMean",
                    "standardized_log2FC",
                    "DE_FDR",
                    "DEG_significant",
                ]
                if column in working.columns
            ]
        ].dropna(subset=["gene_card_key"]).copy()
        self._deg_expression_cache[comparison_id] = result
        return result.copy()

    def _candidate_card_event_frame(self, comparison_id: str | None, gene_key: str) -> tuple[pd.DataFrame, str]:
        if self.current_project is None or not comparison_id or not gene_key:
            return pd.DataFrame(), ""
        comparison = next(
            (item for item in self.current_project.available_comparisons if item.comparison_id == comparison_id),
            None,
        )
        if comparison is None:
            return pd.DataFrame(), ""
        mode = self.current_project.selected_rmats_mode or "JC"
        annotation_lookup = self._annotation_gene_symbol_lookup()
        left_group, right_group = self._comparison_splicing_direction_labels(comparison)
        rows: list[pd.DataFrame] = []
        for event_type in EVENT_TYPES:
            source_path = self._resolve_rmats_event_table_path(comparison, event_type, mode)
            if source_path is None or not source_path.exists():
                continue
            try:
                frame = pd.read_csv(source_path, sep="\t", low_memory=False)
            except Exception:
                continue
            if frame.empty:
                continue
            working = frame.copy()
            working["gene_id"] = self._normalized_gene_id_series(working)
            if "gene_symbol" not in working.columns:
                base_symbol = working.get("geneSymbol")
                if base_symbol is None:
                    base_symbol = pd.Series(pd.NA, index=working.index, dtype="object")
                working["gene_symbol"] = base_symbol
            if annotation_lookup:
                working["gene_symbol"] = working["gene_id"].map(annotation_lookup).combine_first(working["gene_symbol"])
            working["gene_symbol"] = working["gene_symbol"].combine_first(working["gene_id"])
            working["gene_card_key"] = working.apply(
                lambda row: self._gene_lookup_key(row.get("gene_id"), row.get("gene_symbol")),
                axis=1,
            )
            working = working.loc[working["gene_card_key"].astype(str) == gene_key].copy()
            if working.empty:
                continue
            working["comparison_id"] = comparison.comparison_id
            working["comparison_display_name"] = comparison.display_resolved_name
            working["event_type"] = event_type
            working["event_id"] = working.apply(lambda row: event_uid_from_row(row, event_type), axis=1)
            working["event_matching_key"] = working["event_id"].astype(str)
            dpsi = pd.to_numeric(working.get("IncLevelDifference"), errors="coerce")
            if comparison.reverse_splicing:
                dpsi = -dpsi
            working["dPSI"] = dpsi
            working["FDR"] = pd.to_numeric(working.get("FDR"), errors="coerce")
            psi_left = working.get("IncLevel1", pd.Series(pd.NA, index=working.index)).map(self._parse_inc_level_mean)
            psi_right = working.get("IncLevel2", pd.Series(pd.NA, index=working.index)).map(self._parse_inc_level_mean)
            if comparison.reverse_splicing:
                working["psi_experiment"] = psi_right
                working["psi_control"] = psi_left
            else:
                working["psi_experiment"] = psi_left
                working["psi_control"] = psi_right
            working["psi_experiment_group"] = left_group
            working["psi_control_group"] = right_group
            working["significant_event"] = (
                (working["FDR"] <= self.current_thresholds.splicing_fdr)
                & (working["dPSI"].abs() >= self.current_thresholds.splicing_dpsi)
            )
            working["direction"] = working["dPSI"].map(
                lambda value: (
                    "positive in canonical AS direction"
                    if pd.notna(value) and float(value) > 0
                    else "negative in canonical AS direction"
                    if pd.notna(value) and float(value) < 0
                    else "zero / unresolved"
                )
            )
            working["inclusion_direction"] = working["dPSI"].map(
                lambda value: (
                    f"Inclusion higher in {left_group} than {right_group}"
                    if pd.notna(value) and float(value) > 0
                    else f"Inclusion lower in {left_group} than {right_group}"
                    if pd.notna(value) and float(value) < 0
                    else f"No inclusion difference between {left_group} and {right_group}"
                )
            )
            working["coordinates"] = working.apply(lambda row, et=event_type: self._rmats_event_coordinates(row, et), axis=1)
            rows.append(
                working[
                    [
                        "comparison_id",
                        "comparison_display_name",
                        "gene_symbol",
                        "gene_id",
                        "event_type",
                        "event_id",
                        "event_matching_key",
                        "dPSI",
                        "FDR",
                        "psi_experiment",
                        "psi_control",
                        "psi_experiment_group",
                        "psi_control_group",
                        "direction",
                        "inclusion_direction",
                        "coordinates",
                        "significant_event",
                    ]
                ].copy()
            )
        if not rows:
            return pd.DataFrame(), ""
        combined = pd.concat(rows, axis=0, ignore_index=True)
        dominant_source = combined.loc[combined["significant_event"].fillna(False)].copy()
        if dominant_source.empty:
            dominant_source = combined.copy()
            rule = "Dominant event selected by max abs(dPSI) among all events because no event passed the current significance cutoff."
        else:
            rule = "Dominant event selected by max abs(dPSI) among significant events."
        dominant = dominant_source.assign(abs_dPSI=dominant_source["dPSI"].abs()).sort_values(
            by=["abs_dPSI", "FDR", "event_type", "event_id"],
            ascending=[False, True, True, True],
            na_position="last",
        ).iloc[0]
        combined["is_dominant"] = combined["event_id"].astype(str) == str(dominant["event_id"])
        combined["dominant_event_rule"] = rule
        combined = combined.assign(abs_dPSI=combined["dPSI"].abs()).sort_values(
            by=["is_dominant", "significant_event", "abs_dPSI", "FDR", "event_type", "event_id"],
            ascending=[False, False, False, True, True, True],
            na_position="last",
        ).drop(columns=["abs_dPSI"]).reset_index(drop=True)
        return combined, rule

    def candidate_gene_card_payload(
        self,
        comparison_id: str | None,
        *,
        source: str = "selection",
        top_n: int = 20,
        custom_genes: list[str] | None = None,
        allow_generate: bool = False,
        selected_gene_key: str | None = None,
    ) -> tuple[dict[str, object], dict[str, object]]:
        if self.current_project is None:
            return {}, {}
        comparison_meta = self.comparison_context(comparison_id)
        selected_comparisons = self.selected_comparisons_for_display()
        option_rows: list[pd.DataFrame] = []
        current_selection_meta: dict[str, object] = {}
        for comparison in selected_comparisons:
            selected, meta = self.candidate_selection_frame(
                comparison.comparison_id,
                source=source,
                top_n=top_n,
                custom_genes=custom_genes,
                allow_generate=allow_generate,
            )
            if comparison.comparison_id == comparison_id:
                current_selection_meta = meta
            if selected.empty:
                continue
            working = selected.copy()
            working["gene_card_key"] = working.apply(
                lambda row: self._gene_lookup_key(row.get("gene_id"), row.get("gene_symbol")),
                axis=1,
            )
            working["comparison_display_name"] = comparison.display_resolved_name
            option_rows.append(working)
        if not option_rows:
            meta = {
                "comparison_id": comparison_id or "",
                "selection_rule": "No candidate genes selected for card rendering.",
            }
            return {
                "gene_options": [],
                "selected_gene_key": "",
                "comparison_context": comparison_meta,
            }, meta
        option_frame = (
            pd.concat(option_rows, axis=0, ignore_index=True)
            .dropna(subset=["gene_card_key"])
            .drop_duplicates(subset=["gene_card_key"], keep="first")
            .sort_values(by=["rank", "gene_symbol"], ascending=[True, True], na_position="last")
            .reset_index(drop=True)
        )
        gene_options = [
            {
                "key": str(row["gene_card_key"]),
                "gene_symbol": str(row.get("gene_symbol") or row.get("gene_id") or ""),
                "gene_id": str(row.get("gene_id") or ""),
                "label": (
                    f"{row.get('gene_symbol') or row.get('gene_id')} [{row.get('gene_id')}]"
                    if pd.notna(row.get("gene_id")) and str(row.get("gene_id")).strip()
                    else str(row.get("gene_symbol") or row.get("gene_id") or "")
                ),
            }
            for _, row in option_frame.iterrows()
        ]
        valid_keys = {item["key"] for item in gene_options}
        current_key = selected_gene_key if selected_gene_key in valid_keys else (gene_options[0]["key"] if gene_options else "")
        selected_option = next((item for item in gene_options if item["key"] == current_key), None)

        ranking = self._filter_frame_by_comparison(self.preview_candidate_gene_screening(allow_generate=False), comparison_id)
        candidate_row = pd.Series(dtype="object")
        if not ranking.empty:
            ranking = ranking.copy()
            ranking["gene_card_key"] = ranking.apply(
                lambda row: self._gene_lookup_key(row.get("gene_id"), row.get("gene_symbol")),
                axis=1,
            )
            subset = ranking.loc[ranking["gene_card_key"].astype(str) == current_key].copy()
            if not subset.empty:
                subset = subset.sort_values(by=["rank", "gene_symbol"], ascending=[True, True], na_position="last")
                candidate_row = subset.iloc[0]

        fallback_gene_id = selected_option["gene_id"] if selected_option else ""
        fallback_gene_symbol = selected_option["gene_symbol"] if selected_option else ""
        gene_id = str(candidate_row.get("gene_id") or fallback_gene_id or "")
        gene_symbol = str(candidate_row.get("gene_symbol") or fallback_gene_symbol or gene_id or "")

        deg_frame = self._deg_expression_frame(comparison_id)
        deg_row = pd.Series(dtype="object")
        if not deg_frame.empty:
            subset = deg_frame.loc[deg_frame["gene_card_key"].astype(str) == current_key].copy()
            if not subset.empty:
                deg_row = subset.iloc[0]

        expr_support = self._expression_support_frame()
        expr_subset = pd.DataFrame()
        if not expr_support.empty:
            expr_subset = expr_support.loc[expr_support["gene_card_key"].astype(str) == current_key].copy()
        event_frame, dominant_rule = self._candidate_card_event_frame(comparison_id, current_key)
        dominant_event = event_frame.loc[event_frame["is_dominant"].fillna(False)].head(1).copy() if not event_frame.empty else pd.DataFrame()

        meta = dict(current_selection_meta)
        meta.setdefault("comparison_id", comparison_id or "")
        meta["selected_gene_key"] = current_key
        meta["selected_gene_symbol"] = gene_symbol
        meta["selected_gene_id"] = gene_id
        meta["selection_rule"] = meta.get("selection_rule") or "Candidate selection driven gene card."

        return {
            "gene_options": gene_options,
            "selected_gene_key": current_key,
            "gene_symbol": gene_symbol,
            "gene_id": gene_id,
            "candidate_row": candidate_row.to_dict() if not candidate_row.empty else {},
            "deg_row": deg_row.to_dict() if not deg_row.empty else {},
            "expression_support": expr_subset.reset_index(drop=True),
            "event_frame": event_frame,
            "dominant_event": dominant_event.to_dict("records")[0] if not dominant_event.empty else {},
            "dominant_rule": dominant_rule,
            "comparison_context": comparison_meta,
        }, meta

    def candidate_selection_frame(
        self,
        comparison_id: str | None,
        *,
        source: str = "selection",
        top_n: int = 20,
        custom_genes: list[str] | None = None,
        allow_generate: bool = True,
    ) -> tuple[pd.DataFrame, dict[str, object]]:
        frame = self.preview_candidate_gene_screening(allow_generate=allow_generate)
        filtered = self._filter_frame_by_comparison(frame, comparison_id)
        if comparison_id:
            filtered = self._apply_candidate_blacklist(filtered, comparison_id)
        elif not filtered.empty and "blacklist_gene" in filtered.columns:
            filtered = filtered.loc[~filtered["blacklist_gene"].fillna(False)].copy()
        custom_genes = [item.strip() for item in (custom_genes or []) if item and item.strip()]
        source_key = (source or "selection").strip().lower()

        if filtered.empty:
            return pd.DataFrame(), {
                "comparison_id": comparison_id or "",
                "source": source_key,
                "selection_rule": "No candidate rows available.",
                "gene_count": 0,
                "custom_genes": custom_genes,
            }

        if source_key == "selection":
            selected = self._candidate_selected_gene_rows(filtered, comparison_id, top_n=top_n, allow_generate=allow_generate)
            selection_rule = f"Candidate gene selection for this comparison. Default top {top_n} by numeric rank unless manually edited."
        elif source_key == "tier1":
            selected = filtered.loc[filtered["candidate_tier"].astype(str) == "Tier 1"].copy()
            selection_rule = "Tier 1 candidates only."
        elif source_key in {"tier1_tier2", "tier1+tier2", "tier12"}:
            selected = filtered.loc[filtered["candidate_tier"].astype(str).isin(["Tier 1", "Tier 2"])].copy()
            selection_rule = "Tier 1 + Tier 2 candidates."
        elif source_key == "shortlist":
            selected = filtered.loc[filtered["shortlist_gene"].fillna(False)].copy()
            selection_rule = "Manual shortlist genes and auto-marked shortlist genes."
        elif source_key == "custom":
            selected = self._candidate_custom_gene_rows(filtered, comparison_id, custom_genes)
            selection_rule = "User-provided custom gene list."
        else:
            selected = filtered.copy()
            selection_rule = "Current filtered candidate rows for this comparison."

        if source_key not in {"selection", "shortlist", "custom"} and top_n > 0:
            selected = selected.sort_values(["rank", "gene_symbol"], ascending=[True, True], na_position="last").head(top_n).copy()
            selection_rule += f" Top N = {top_n}."

        return selected.reset_index(drop=True), {
            "comparison_id": comparison_id or "",
            "source": source_key,
            "selection_rule": selection_rule,
            "gene_count": int(len(selected)),
            "top_n": int(top_n),
            "custom_genes": custom_genes,
        }

    def candidate_selection_catalog(self, comparison_id: str | None, *, allow_generate: bool = True) -> pd.DataFrame:
        frame = self.preview_candidate_gene_screening(allow_generate=allow_generate)
        filtered = self._filter_frame_by_comparison(frame, comparison_id)
        if filtered.empty:
            return pd.DataFrame(columns=["comparison_id", "rank", "gene_symbol", "gene_id", "candidate_tier", "evidence_class", "direction_class"])
        filtered = filtered.sort_values(["rank", "gene_symbol"], ascending=[True, True], na_position="last").copy()
        columns = [
            column for column in [
                "comparison_id",
                "rank",
                "gene_symbol",
                "gene_id",
                "candidate_tier",
                "evidence_class",
                "direction_class",
                "DEG_significant",
                "rMATS_significant",
                "DEXSeq_significant",
                "DTU_significant",
                "DE_FDR",
                "best_rMATS_FDR",
                "best_DEXSeq_qvalue",
                "best_DTU_qvalue",
                "log2fc_direction_label",
                "dpsi_direction_label",
            ] if column in filtered.columns
        ]
        return filtered[columns].copy().reset_index(drop=True)

    def _apply_candidate_blacklist(self, frame: pd.DataFrame, comparison_id: str) -> pd.DataFrame:
        if self.current_project is None or frame.empty:
            return frame
        blocked = {
            item.strip()
            for item in self.current_project.candidate_blacklist_genes.get(comparison_id, [])
            if item and item.strip()
        }
        if not blocked:
            return frame
        filtered = frame.copy()
        keep_mask = pd.Series(True, index=filtered.index)
        for column in ("gene_symbol", "geneSymbol", "gene_id", "GeneID"):
            if column in filtered.columns:
                keep_mask &= ~filtered[column].astype(str).isin(blocked)
        return filtered.loc[keep_mask].copy()

    def _default_candidate_selection_genes(self, comparison_id: str, *, top_n: int = 20, allow_generate: bool = True) -> list[str]:
        frame = self._filter_frame_by_comparison(self.preview_candidate_gene_screening(allow_generate=allow_generate), comparison_id)
        if frame.empty:
            return []
        working = frame.sort_values(["rank", "gene_symbol"], ascending=[True, True], na_position="last").head(top_n).copy()
        genes: list[str] = []
        seen: set[str] = set()
        for _, row in working.iterrows():
            gene = str(row.get("gene_symbol") or row.get("gene_id") or "").strip()
            if gene and gene not in seen:
                genes.append(gene)
                seen.add(gene)
        return genes

    def _ensure_default_candidate_selection(self, comparison_id: str, *, top_n: int = 20, allow_generate: bool = True) -> None:
        if self.current_project is None or not comparison_id:
            return
        existing = self.current_project.candidate_selection_genes.get(comparison_id, [])
        if existing:
            return
        if not allow_generate:
            return
        self.current_project.candidate_selection_genes[comparison_id] = self._default_candidate_selection_genes(
            comparison_id,
            top_n=top_n,
            allow_generate=allow_generate,
        )
        self.save_project_config()

    def _candidate_selected_gene_rows(self, filtered: pd.DataFrame, comparison_id: str | None, *, top_n: int = 20, allow_generate: bool = True) -> pd.DataFrame:
        if self.current_project is None or not comparison_id:
            return filtered.sort_values(["rank", "gene_symbol"], ascending=[True, True], na_position="last").head(top_n).copy()
        self._ensure_default_candidate_selection(comparison_id, top_n=top_n, allow_generate=allow_generate)
        selected_genes = [
            item.strip()
            for item in self.current_project.candidate_selection_genes.get(comparison_id, [])
            if item and item.strip()
        ]
        if not selected_genes:
            return pd.DataFrame(columns=filtered.columns)
        selected = self._candidate_custom_gene_rows(filtered, comparison_id, selected_genes)
        return selected.sort_values(["rank", "gene_symbol"], ascending=[True, True], na_position="last").reset_index(drop=True)

    def _candidate_support_table_is_stale(self, frame: pd.DataFrame) -> bool:
        if frame is None or frame.empty:
            return False
        if "DEXSeq_significant" not in frame.columns or "DTU_significant" not in frame.columns:
            return True
        has_any_candidate_support = bool(
            frame["DEXSeq_significant"].fillna(False).astype(bool).any()
            or frame["DTU_significant"].fillna(False).astype(bool).any()
        )
        if has_any_candidate_support:
            return False
        mechanism = self.preview_mechanism_support(allow_generate=False)
        if mechanism is None or mechanism.empty:
            return False
        support_cols = [column for column in ("support_DEXSeq", "support_DTU") if column in mechanism.columns]
        if not support_cols:
            return False
        support_present = False
        for column in support_cols:
            values = pd.to_numeric(mechanism[column], errors="coerce").fillna(0)
            if bool((values > 0).any()):
                support_present = True
                break
        return support_present

    def candidate_isoform_followup_frame(
        self,
        comparison_id: str | None,
        *,
        source: str = "shortlist",
        top_n: int = 20,
        custom_genes: list[str] | None = None,
        allow_generate: bool = True,
    ) -> tuple[pd.DataFrame, dict[str, object]]:
        selected, meta = self.candidate_selection_frame(
            comparison_id,
            source=source,
            top_n=top_n,
            custom_genes=custom_genes,
            allow_generate=allow_generate,
        )
        if selected.empty:
            empty = pd.DataFrame(
                columns=[
                    "comparison_id",
                    "gene_symbol",
                    "gene_id",
                    "candidate_tier",
                    "DTU_significant",
                    "best_DTU_qvalue",
                    "n_DTU_significant_isoforms",
                    "transcript_id",
                    "isoform_id",
                    "DTU_qvalue",
                    "isoform_followup_status",
                ]
            )
            return empty, meta

        followup = selected[
            [
                "comparison_id",
                "gene_symbol",
                "gene_id",
                "candidate_tier",
                "DTU_significant",
                "best_DTU_qvalue",
                "n_DTU_significant_isoforms",
                "shortlist_gene",
                "blacklist_gene",
            ]
        ].copy()
        followup["transcript_id"] = followup["DTU_significant"].map(lambda value: "DTU-supported transcript set" if bool(value) else "no DTU evidence")
        followup["isoform_id"] = followup["DTU_significant"].map(lambda value: "DTU-supported isoform set" if bool(value) else "no DTU evidence")
        followup["DTU_qvalue"] = followup["best_DTU_qvalue"]
        followup["isoform_followup_status"] = followup["DTU_significant"].map(
            lambda value: "DTU evidence available" if bool(value) else "no DTU evidence"
        )
        return followup, meta

    def candidate_event_followup_frame(
        self,
        comparison_id: str | None,
        *,
        source: str = "selection",
        top_n: int = 20,
        custom_genes: list[str] | None = None,
        allow_generate: bool = True,
    ) -> tuple[pd.DataFrame, dict[str, object]]:
        selected, meta = self.candidate_selection_frame(
            comparison_id,
            source=source,
            top_n=top_n,
            custom_genes=custom_genes,
            allow_generate=allow_generate,
        )
        if selected.empty:
            empty = pd.DataFrame(
                columns=[
                    "comparison_id",
                    "gene_symbol",
                    "gene_id",
                    "candidate_tier",
                    "direction_class",
                    "best_rMATS_event_id",
                    "dominant_rMATS_event_type",
                    "dominant_rMATS_standardized_dPSI",
                    "best_rMATS_FDR",
                    "DE_FDR",
                    "standardized_log2FC",
                    "shortlist_gene",
                    "blacklist_gene",
                ]
            )
            return empty, meta
        columns = [
            column
            for column in [
                "comparison_id",
                "gene_symbol",
                "gene_id",
                "candidate_tier",
                "evidence_class",
                "direction_class",
                "best_rMATS_event_id",
                "dominant_rMATS_event_type",
                "dominant_rMATS_standardized_dPSI",
                "best_rMATS_FDR",
                "DE_FDR",
                "standardized_log2FC",
                "shortlist_gene",
                "blacklist_gene",
                "candidate_reason",
            ]
            if column in selected.columns
        ]
        return selected[columns].copy(), meta

    def available_sashimi_events(
        self,
        comparison_id: str | None,
        *,
        allow_generate: bool = False,
    ) -> pd.DataFrame:
        _ = allow_generate  # Explicitly unused: opening the page must not generate downstream inputs.
        return self._sashimi_event_catalog_frame(comparison_id)

    def _empty_sashimi_events_frame(self) -> pd.DataFrame:
        return pd.DataFrame(
            columns=[
                "comparison_id",
                "comparison_display_name",
                "gene_symbol",
                "gene_id",
                "event_type",
                "event_id",
                "event_matching_key",
                "dPSI",
                "FDR",
                "direction",
                "coordinates",
            ]
        )

    def _sashimi_event_catalog_frame(self, comparison_id: str | None) -> pd.DataFrame:
        if self.current_project is None or not comparison_id:
            return self._empty_sashimi_events_frame()
        comparison = next(
            (item for item in self.current_project.available_comparisons if item.comparison_id == comparison_id),
            None,
        )
        if comparison is None or comparison.rmats_path is None:
            return self._empty_sashimi_events_frame()
        mode = self.current_project.selected_rmats_mode or "JC"
        cache_key = f"{comparison.comparison_id}|{mode}|{int(bool(comparison.reverse_splicing))}"
        cached = self._sashimi_event_catalog_cache.get(cache_key)
        if cached is not None:
            return cached.copy()

        annotation_lookup = self._annotation_gene_symbol_lookup()
        rows: list[pd.DataFrame] = []
        display_name = comparison.display_resolved_name
        for event_type in EVENT_TYPES:
            source_path = self._resolve_rmats_event_table_path(comparison, event_type, mode)
            if source_path is None or not source_path.exists():
                continue
            try:
                frame = pd.read_csv(source_path, sep="\t", low_memory=False)
            except Exception:
                continue
            if frame.empty:
                continue
            working = frame.copy()
            working["comparison_id"] = comparison.comparison_id
            working["comparison_display_name"] = display_name
            working["event_type"] = event_type
            working["event_id"] = working.apply(lambda row: event_uid_from_row(row, event_type), axis=1)
            working["event_matching_key"] = working["event_id"].astype(str)

            if "GeneID" in working.columns:
                working["gene_id"] = working["GeneID"]
            elif "gene_id" not in working.columns:
                working["gene_id"] = pd.NA
            normalized_gene_ids = working["gene_id"].map(normalize_gene_id)

            if "geneSymbol" in working.columns:
                base_symbol = working["geneSymbol"].astype("object")
            elif "gene_symbol" in working.columns:
                base_symbol = working["gene_symbol"].astype("object")
            else:
                base_symbol = pd.Series(pd.NA, index=working.index, dtype="object")
            working["gene_symbol"] = base_symbol
            missing_symbol = working["gene_symbol"].isna() | (working["gene_symbol"].astype(str).str.strip() == "")
            if missing_symbol.any():
                working.loc[missing_symbol, "gene_symbol"] = normalized_gene_ids.loc[missing_symbol].map(annotation_lookup)
            working.loc[working["gene_symbol"].isna(), "gene_symbol"] = normalized_gene_ids.loc[working["gene_symbol"].isna()]

            dpsi_series = pd.to_numeric(
                working["IncLevelDifference"] if "IncLevelDifference" in working.columns else working.get("dPSI"),
                errors="coerce",
            )
            if comparison.reverse_splicing:
                dpsi_series = -dpsi_series
            working["dPSI"] = dpsi_series
            working["FDR"] = pd.to_numeric(working["FDR"] if "FDR" in working.columns else working.get("fdr"), errors="coerce")
            working["direction"] = working["dPSI"].map(
                lambda value: (
                    "positive in canonical AS direction"
                    if pd.notna(value) and float(value) > 0
                    else "negative in canonical AS direction"
                    if pd.notna(value) and float(value) < 0
                    else "zero / unresolved"
                )
            )
            working["coordinates"] = working.apply(
                lambda row, et=event_type: self._rmats_event_coordinates(row, et),
                axis=1,
            )
            rows.append(
                working[
                    [
                        "comparison_id",
                        "comparison_display_name",
                        "gene_symbol",
                        "gene_id",
                        "event_type",
                        "event_id",
                        "event_matching_key",
                        "dPSI",
                        "FDR",
                        "direction",
                        "coordinates",
                    ]
                ].copy()
            )

        if not rows:
            empty = self._empty_sashimi_events_frame()
            self._sashimi_event_catalog_cache[cache_key] = empty
            return empty.copy()

        combined = pd.concat(rows, ignore_index=True)
        combined = combined.drop_duplicates(subset=["comparison_id", "event_matching_key"]).reset_index(drop=True)
        combined = combined.sort_values(
            by=["FDR", "event_type", "gene_symbol", "event_id"],
            ascending=[True, True, True, True],
            na_position="last",
        ).reset_index(drop=True)
        self._sashimi_event_catalog_cache[cache_key] = combined.copy()
        return combined

    @staticmethod
    def _rmats_event_coordinates(row: pd.Series, event_type: str) -> str:
        ordered_columns = [
            "chr",
            "strand",
            "exonStart_0base",
            "exonEnd",
            "upstreamES",
            "upstreamEE",
            "downstreamES",
            "downstreamEE",
            "riExonStart_0base",
            "riExonEnd",
            "longExonStart_0base",
            "longExonEnd",
            "shortES",
            "shortEE",
            "flankingES",
            "flankingEE",
            "1stExonStart_0base",
            "1stExonEnd",
            "2ndExonStart_0base",
            "2ndExonEnd",
        ]
        parts = [event_type]
        for column in ordered_columns:
            if column in row.index:
                value = row[column]
                if pd.notna(value):
                    parts.append(f"{column}={value}")
        return " | ".join(parts)

    def build_sashimi_manifest_for_events(
        self,
        comparison_id: str | None,
        event_ids: list[str],
        *,
        allow_generate: bool = False,
    ) -> pd.DataFrame:
        if self.current_project is None:
            raise RuntimeError("No project loaded.")
        clean_ids = [str(item).strip() for item in event_ids if str(item).strip()]
        if not clean_ids:
            raise RuntimeError("No events selected. Please select one or more splicing events first.")
        events = self.available_sashimi_events(comparison_id, allow_generate=allow_generate)
        if events.empty:
            self._mark_module_state(
                "sashimi_manifest",
                "not_run",
                "Missing input: event follow-up source not available. Please run 5.3.2 / 5.3.4 first.",
            )
            return pd.DataFrame()
        source = events.loc[events["event_id"].astype(str).isin(clean_ids)].copy()
        if source.empty:
            raise RuntimeError("The selected splicing events are not available in the current candidate event list.")
        shortlist = source.copy()
        shortlist["event_uid"] = shortlist["event_id"].astype(str)
        shortlist["GeneID"] = shortlist["gene_id"].astype(str)
        shortlist["geneSymbol"] = shortlist["gene_symbol"].astype(str)
        shortlist["display_group"] = shortlist["comparison_display_name"].astype(str)
        shortlist["class_interpretation"] = shortlist.get("direction_class", pd.Series(dtype=str)).astype(str)
        raw_manifest = self.load_first_matching_table("rmats2sashimi_manifest.tsv")
        manifest = self._filter_sashimi_manifest(raw_manifest, shortlist)
        if manifest.empty:
            manifest = self._build_generated_sashimi_manifest(shortlist)
        manifest = manifest.reset_index(drop=True) if manifest is not None and not manifest.empty else pd.DataFrame()
        self.run_state.results.sashimi_manifest = manifest.copy()
        self._mark_module_state(
            "sashimi_manifest",
            "finished" if not manifest.empty else "failed",
            "Generated sashimi manifest for the selected events." if not manifest.empty else "No sashimi manifest rows could be generated for the selected events.",
            output_folder=str(self.current_project.output_root / "06_sashimi") if self.current_project.output_root else "",
            from_cache=False,
        )
        return manifest

    def last_sashimi_failures(self) -> pd.DataFrame:
        return self._last_sashimi_failures.copy()

    def write_candidate_selection_outputs(
        self,
        comparison_id: str | None,
        *,
        kind: str,
        selected: pd.DataFrame,
        parameters: dict[str, object],
    ) -> Path | None:
        if self.current_project is None or self.current_project.output_root is None or not comparison_id:
            return None
        comparison = next(
            (item for item in self.current_project.available_comparisons if item.comparison_id == comparison_id),
            None,
        )
        if comparison is None:
            return None
        comp_dir = self.current_project.output_root / "02_candidate_gene_screening" / (
            comparison.output_prefix or self._safe_output_name(comparison.comparison_id)
        )
        kind_dir = comp_dir / kind
        kind_dir.mkdir(parents=True, exist_ok=True)
        gene_columns = [column for column in ["comparison_id", "gene_symbol", "gene_id", "candidate_tier", "evidence_class", "direction_class"] if column in selected.columns]
        selected[gene_columns].to_csv(kind_dir / "heatmap_gene_list.tsv", sep="\t", index=False)
        self._write_json(kind_dir / "heatmap_parameters.json", parameters)
        summary = pd.DataFrame(
            [
                {
                    "comparison_id": comparison_id,
                    "kind": kind,
                    "source": parameters.get("source", ""),
                    "gene_count": len(selected),
                }
            ]
        )
        summary.to_csv(kind_dir / "summary.tsv", sep="\t", index=False)
        return kind_dir

    def active_visualization_groups_for_display(self) -> list[VisualizationGroupDefinition]:
        return list(self._active_visualization_groups())

    def selected_analysis_modules_for_display(self) -> list[str]:
        if self.current_project is None:
            return list(self.ANALYSIS_MODULES)
        return list(self.current_project.selected_analysis_modules or self.ANALYSIS_MODULES)

    def save_project_config(self) -> None:
        if self.current_project is None or self.current_project.config_path is None:
            return
        payload = {
            "selected_rmats_mode": self.current_project.selected_rmats_mode,
            "output_root": str(self.current_project.output_root) if self.current_project.output_root else None,
            "confirmed": self.current_project.confirmed,
            "pairing_confirmed": self.current_project.pairing_confirmed,
            "comparison_sets_confirmed": self.current_project.comparison_sets_confirmed,
            "visualization_groups_confirmed": self.current_project.visualization_groups_confirmed,
            "selected_comparison_ids": list(self.current_project.selected_comparison_ids),
            "comparison_order_ids": list(self.current_project.comparison_order_ids),
            "program_comparison_ids": list(self.current_project.program_comparison_ids),
            "selected_analysis_modules": list(self.current_project.selected_analysis_modules),
            "removed_comparison_ids": list(self.current_project.removed_comparison_ids),
            "comparison_pairs": [
                {
                    "pair_id": pair.pair_id,
                    "comparison_a": pair.comparison_a,
                    "comparison_b": pair.comparison_b,
                    "display_name": pair.display_name,
                    "experiment_group": pair.experiment_group,
                    "control_group": pair.control_group,
                    "enabled": pair.enabled,
                }
                for pair in self.current_project.comparison_pairs
            ],
            "visualization_groups": [
                {
                    "group_id": group.group_id,
                    "display_name": group.display_name,
                    "comparison_ids": list(group.comparison_ids),
                    "enabled": group.enabled,
                }
                for group in self.current_project.visualization_groups
            ],
            "shortlist_genes": list(self.current_project.shortlist_genes),
            "blacklist_genes": list(self.current_project.blacklist_genes),
            "candidate_selection_genes": {
                comparison_id: list(values)
                for comparison_id, values in self.current_project.candidate_selection_genes.items()
            },
            "candidate_blacklist_genes": {
                comparison_id: list(values)
                for comparison_id, values in self.current_project.candidate_blacklist_genes.items()
            },
            "thresholds": {
                "splicing_fdr": self.current_thresholds.splicing_fdr,
                "splicing_dpsi": self.current_thresholds.splicing_dpsi,
                "program_delta_dpsi": self.current_thresholds.program_delta_dpsi,
                "deg_padj": self.current_thresholds.deg_padj,
                "deg_log2fc": self.current_thresholds.deg_log2fc,
                "dexseq_qvalue": self.current_thresholds.dexseq_qvalue,
                "dtu_qvalue": self.current_thresholds.dtu_qvalue,
            },
            "tool_paths": {
                key: value
                for key, value in self.current_project.tool_paths.items()
                if key in {"rscript", "gtf", "annotation_tsv", "fasta", "jutils", "rmats2sashimiplot", "suppa", "bam_root"}
            },
            "input_paths": dict(self.current_project.input_paths),
            "comparisons": [
                {
                    "comparison_id": item.comparison_id,
                    "display_name": item.display_name,
                    "biological_question": item.biological_question,
                    "rmats_name": item.rmats_name,
                    "deg_name": item.deg_name,
                    "suppa_name": item.suppa_name,
                    "dexseq_name": item.dexseq_name,
                    "dtu_name": item.dtu_name,
                    "quant_name": item.quant_name,
                    "sashimi_name": item.sashimi_name,
                    "bam_group_name": item.bam_group_name,
                    "source_experiment_group": item.source_experiment_group,
                    "source_control_group": item.source_control_group,
                    "experiment_group": item.experiment_group,
                    "control_group": item.control_group,
                    "source_direction": item.source_direction,
                    "analysis_direction": item.analysis_direction,
                    "reverse_direction": item.reverse_direction,
                    "reverse_splicing": item.reverse_splicing,
                    "reverse_deg": item.reverse_deg,
                    "expected_direction_notes": item.expected_direction_notes,
                    "known_positive_control_genes": list(item.known_positive_control_genes),
                    "known_negative_control_genes": list(item.known_negative_control_genes),
                    "output_prefix": item.output_prefix,
                }
                for item in self.current_project.available_comparisons
            ],
            "isoform_samples": [
                {
                    "sample_id": item.sample_id,
                    "quant_sf": str(item.quant_sf),
                    "sample_group": item.sample_group,
                    "comparison_id": item.comparison_id,
                    "condition": item.condition,
                    "batch": item.batch,
                    "include": item.include,
                }
                for item in self.current_project.isoform_samples
            ],
        }
        with self.current_project.config_path.open("w", encoding="utf-8") as handle:
            yaml.safe_dump(payload, handle, sort_keys=False, allow_unicode=True)

    def _refresh_tool_adapters(self) -> None:
        if self.current_project is None:
            self.jutils_adapter = None
            self.sashimi_adapter = None
            self.isoform_runner = None
            return
        jutils_root = self.current_project.tool_paths.get("jutils")
        if jutils_root:
            script_path = Path(jutils_root) / "jutils.py"
            self.jutils_adapter = JutilsAdapter(script_path) if script_path.exists() else None
        else:
            self.jutils_adapter = None
        sashimi_root = self.current_project.tool_paths.get("rmats2sashimiplot")
        if sashimi_root:
            script_path = Path(sashimi_root) / "src" / "rmats2sashimiplot" / "rmats2sashimiplot.py"
            self.sashimi_adapter = SashimiAdapter(script_path) if script_path.exists() else None
        else:
            self.sashimi_adapter = None
        rscript_value = self.current_project.tool_paths.get("rscript")
        self.isoform_runner = IsoformSwitchRunner(Path(rscript_value)) if rscript_value else IsoformSwitchRunner(None)

    def _comparison_quant_files(self, comparison) -> list[Path]:
        if self.current_project is None:
            return []
        groups = {
            comparison.source_experiment_group or comparison.experiment_group,
            comparison.source_control_group or comparison.control_group,
        }
        groups = {group for group in groups if group}
        matched = [
            sample.quant_sf
            for sample in self.current_project.isoform_samples
            if sample.quant_sf.exists()
            and (
                (sample.sample_group in groups if sample.sample_group else False)
                or (sample.comparison_id == comparison.comparison_id if sample.comparison_id else False)
            )
        ]
        return sorted(set(matched))

    def _resolve_comparison_support_file(self, comparison, source_type: str) -> str | None:
        if self.current_project is None:
            return None

        root_map = {
            "suppa": self.current_project.suppa_root,
            "dexseq": self.current_project.dexseq_root,
            "dtu": self.current_project.dtu_root,
            "stager": self.current_project.dtu_root,
        }
        root = root_map.get(source_type)
        if root is None or not root.exists():
            return None

        for candidate in self._comparison_name_candidates(comparison, source_type):
            patterns = []
            if source_type == "suppa":
                patterns = [
                    f"{candidate}_local_diffsplice.dpsi",
                    f"*{candidate}*_local_diffsplice.dpsi",
                ]
            elif source_type == "dexseq":
                patterns = [
                    f"perGeneQValue.{candidate}.csv",
                    f"perGeneQValue.*{candidate}*.csv",
                ]
            elif source_type == "dtu":
                patterns = [
                    f"perGeneQValue.{candidate}.tsv",
                    f"perGeneQValue.*{candidate}*.tsv",
                    f"DEXSeqResults.{candidate}.tsv",
                    f"DEXSeqResults.*{candidate}*.tsv",
                ]
            elif source_type == "stager":
                patterns = [
                    f"getAdjustedPValues.{candidate}.tsv",
                    f"getAdjustedPValues.*{candidate}*.tsv",
                ]
            for pattern in patterns:
                matches = sorted(root.rglob(pattern))
                if matches:
                    return str(matches[0])
        return None

    def _comparison_name_candidates(self, comparison, source_type: str) -> list[str]:
        candidates: list[str] = []
        experiment = comparison.experiment_group or comparison.source_experiment_group
        control = comparison.control_group or comparison.source_control_group
        for name in (
            comparison.source_name(source_type if source_type != "stager" else "dtu"),
            comparison.comparison_id,
            f"{experiment}_vs_{control}" if experiment and control else None,
            f"{control}_vs_{experiment}" if experiment and control else None,
            f"{experiment}-{control}" if experiment and control else None,
            f"{control}-{experiment}" if experiment and control else None,
        ):
            if name and name not in candidates:
                candidates.append(name)
                if "_vs_" in name and name.replace("_vs_", "-") not in candidates:
                    candidates.append(name.replace("_vs_", "-"))
                if "-vs-" in name and name.replace("-vs-", "-") not in candidates:
                    candidates.append(name.replace("-vs-", "-"))
                if " vs " in name and name.replace(" vs ", "-") not in candidates:
                    candidates.append(name.replace(" vs ", "-"))
                if "-" in name and "_vs_" not in name and "-vs-" not in name and " vs " not in name:
                    dashed_to_vs = name.replace("-", "_vs_")
                    if dashed_to_vs not in candidates:
                        candidates.append(dashed_to_vs)
        return candidates

    def _write_isoform_design_template(self, data_dir: Path, manifest: pd.DataFrame) -> Path:
        design_file = data_dir / "isoform_design.tsv"
        frame = manifest.copy()
        if "condition" not in frame.columns:
            frame["condition"] = frame["sample_id"].astype(str)
        if "batch" not in frame.columns:
            frame["batch"] = ""
        frame = frame.rename(columns={"sample_id": "sampleID"})
        output = frame[["sampleID", "condition", "batch"]].copy()
        output.to_csv(design_file, sep="\t", index=False)
        return design_file

    def _autodetect_embedded_tools(self) -> None:
        if self.current_project is None:
            return
        candidate_roots = [
            self.current_project.project_root / "software",
            bundled_software_root(),
        ]
        for software_root in candidate_roots:
            if "jutils" not in self.current_project.tool_paths and (software_root / "Jutils").exists():
                self.current_project.tool_paths["jutils"] = str(software_root / "Jutils")
            if "rmats2sashimiplot" not in self.current_project.tool_paths and (software_root / "rmats2sashimiplot").exists():
                self.current_project.tool_paths["rmats2sashimiplot"] = str(software_root / "rmats2sashimiplot")
            if "suppa" not in self.current_project.tool_paths and (software_root / "SUPPA").exists():
                self.current_project.tool_paths["suppa"] = str(software_root / "SUPPA")
        self._refresh_tool_adapters()

    def _safe_output_name(self, value: str) -> str:
        cleaned = "".join(ch if ch.isalnum() or ch in ("-", "_", ".") else "_" for ch in value).strip("_")
        return cleaned or "comparison"

    def _command_to_text(self, command) -> str:
        if isinstance(command, (list, tuple)):
            return " ".join(str(part) for part in command)
        return str(command)

    def _preview_kind_for_path(self, path: Path) -> str:
        suffix = path.suffix.lower()
        if suffix in {".png", ".jpg", ".jpeg", ".svg"}:
            return "image"
        if suffix in {".tsv", ".csv", ".txt"}:
            return "table"
        if suffix in {".xlsx", ".xls"}:
            return "excel"
        if suffix == ".html":
            return "html"
        if suffix == ".pdf":
            return "pdf"
        return "binary"

    def _selected_comparisons(self):
        if self.current_project is None:
            return []
        if not self.current_project.selected_comparison_ids:
            return []
        selected_ids = set(self.current_project.selected_comparison_ids)
        lookup = {comparison.comparison_id: comparison for comparison in self.current_project.available_comparisons}
        ordered_ids = [item for item in self._ordered_available_comparison_ids() if item in selected_ids]
        return [lookup[item] for item in ordered_ids if item in lookup]

    def _ordered_available_comparison_ids(self) -> list[str]:
        if self.current_project is None:
            return []
        available_ids = [item.comparison_id for item in self.current_project.available_comparisons]
        ordered = [item for item in self.current_project.comparison_order_ids if item in available_ids]
        for comparison_id in available_ids:
            if comparison_id not in ordered:
                ordered.append(comparison_id)
        return ordered

    def _apply_comparison_order_to_available(self) -> None:
        if self.current_project is None:
            return
        ordered_ids = self._ordered_available_comparison_ids()
        lookup = {item.comparison_id: item for item in self.current_project.available_comparisons}
        self.current_project.available_comparisons = [lookup[item] for item in ordered_ids if item in lookup]

    def _ensure_default_visualization_groups(self) -> None:
        if self.current_project is None or self.current_project.visualization_groups:
            return
        selected = self._selected_comparisons()
        if not selected:
            return
        default_ids = [item.comparison_id for item in selected[: min(4, len(selected))]]
        self.current_project.visualization_groups = [
            VisualizationGroupDefinition(
                group_id="viz_1",
                comparison_ids=default_ids,
                display_name=self._default_visualization_group_name(default_ids),
            )
        ]

    def _active_visualization_groups(self) -> list[VisualizationGroupDefinition]:
        if self.current_project is None:
            return []
        valid_ids = {item.comparison_id for item in self._selected_comparisons()}
        return [
            group
            for group in self.current_project.visualization_groups
            if group.enabled and any(item in valid_ids for item in group.comparison_ids)
        ]

    def _ensure_default_comparison_pairs(self) -> None:
        if self.current_project is None or self.current_project.comparison_pairs:
            return
        selected = self._selected_comparisons()
        if len(selected) < 2:
            return
        first, second = selected[:2]
        self.current_project.comparison_pairs = [
            ComparisonPairDefinition(
                pair_id="pair_1",
                comparison_a=first.comparison_id,
                comparison_b=second.comparison_id,
                experiment_group=first.experiment_group,
                control_group=second.experiment_group,
            )
        ]
        self._autofill_comparison_pair(self.current_project.comparison_pairs[0])

    def _active_comparison_pairs(self) -> list[ComparisonPairDefinition]:
        if self.current_project is None:
            return []
        valid_ids = {item.comparison_id for item in self._selected_comparisons()}
        return [
            pair
            for pair in self.current_project.comparison_pairs
            if pair.enabled and pair.comparison_a in valid_ids and pair.comparison_b in valid_ids
        ]

    def _autofill_comparison_pair(
        self,
        pair: ComparisonPairDefinition,
        *,
        update_display: bool = True,
        update_experiment: bool = True,
        update_control: bool = True,
    ) -> None:
        if self.current_project is None:
            return
        if update_display:
            pair.display_name = self._default_pair_display_name(pair)
        if update_experiment:
            pair.experiment_group = self._default_pair_experiment_group(pair)
        if update_control:
            pair.control_group = self._default_pair_control_group(pair)

    def _pair_comparison_lookup(self) -> dict[str, object]:
        if self.current_project is None:
            return {}
        return {item.comparison_id: item for item in self.current_project.available_comparisons}

    def _pair_endpoints(self, pair: ComparisonPairDefinition):
        lookup = self._pair_comparison_lookup()
        return lookup.get(pair.comparison_a or ""), lookup.get(pair.comparison_b or "")

    def _resolve_comparison_reference(self, reference: str | None) -> str | None:
        if self.current_project is None or not reference:
            return None
        reference = str(reference)
        for comparison in self.current_project.available_comparisons:
            if reference in {
                comparison.comparison_id,
                comparison.display_resolved_name,
                comparison.resolved_name,
                comparison.display_name or "",
            }:
                return comparison.comparison_id
        return reference or None

    def _default_visualization_group_name(self, comparison_ids: list[str]) -> str | None:
        if self.current_project is None or not comparison_ids:
            return None
        lookup = {item.comparison_id: item for item in self.current_project.available_comparisons}
        labels = [lookup[item].display_resolved_name for item in comparison_ids if item in lookup]
        if not labels:
            return None
        return " + ".join(labels)

    def _default_pair_display_name(self, pair: ComparisonPairDefinition) -> str | None:
        left, right = self._pair_endpoints(pair)
        if left is not None and right is not None:
            return f"{left.display_resolved_name} compared with {right.display_resolved_name}"
        return None

    def _default_pair_experiment_group(self, pair: ComparisonPairDefinition) -> str | None:
        left, _ = self._pair_endpoints(pair)
        if left is not None:
            return left.experiment_group or left.source_experiment_group
        return None

    def _default_pair_control_group(self, pair: ComparisonPairDefinition) -> str | None:
        _, right = self._pair_endpoints(pair)
        if right is not None:
            return right.experiment_group or right.source_experiment_group
        return None

    @staticmethod
    def _filter_frame_by_comparison(frame: pd.DataFrame, comparison_id: str | None) -> pd.DataFrame:
        if frame is None or frame.empty:
            return pd.DataFrame()
        if not comparison_id or "comparison_id" not in frame.columns:
            return frame.copy()
        return frame.loc[frame["comparison_id"].astype(str) == str(comparison_id)].copy()

    def _candidate_custom_gene_rows(
        self,
        frame: pd.DataFrame,
        comparison_id: str | None,
        genes: list[str],
    ) -> pd.DataFrame:
        if frame.empty or not genes:
            return pd.DataFrame(columns=frame.columns if frame is not None else [])
        normalized = {item.strip().lower() for item in genes if item.strip()}
        if not normalized:
            return pd.DataFrame(columns=frame.columns)
        matched = frame.loc[
            frame["gene_symbol"].fillna("").astype(str).str.lower().isin(normalized)
            | frame["gene_id"].fillna("").astype(str).str.lower().isin(normalized)
        ].copy()
        missing = [
            item
            for item in genes
            if item.strip()
            and item.strip().lower() not in set(matched["gene_symbol"].fillna("").astype(str).str.lower())
            and item.strip().lower() not in set(matched["gene_id"].fillna("").astype(str).str.lower())
        ]
        if not missing:
            return matched
        template = frame.iloc[0].copy()
        rows = []
        for gene in missing:
            row = template.copy()
            row["comparison_id"] = comparison_id or row.get("comparison_id", "")
            row["gene_symbol"] = gene
            row["gene_id"] = gene
            row["candidate_tier"] = "Unclassified"
            row["evidence_class"] = "manual_shortlist_only"
            row["evidence_count"] = 0
            row["DEG_significant"] = False
            row["rMATS_significant"] = False
            row["DEXSeq_significant"] = False
            row["DTU_significant"] = False
            row["DE_FDR"] = pd.NA
            row["best_rMATS_FDR"] = pd.NA
            row["best_DEXSeq_qvalue"] = pd.NA
            row["best_DTU_qvalue"] = pd.NA
            row["standardized_log2FC"] = pd.NA
            row["dominant_rMATS_standardized_dPSI"] = pd.NA
            row["max_abs_dPSI"] = pd.NA
            row["n_rMATS_significant_events"] = 0
            row["dominant_rMATS_event_type"] = ""
            row["best_rMATS_event_id"] = ""
            row["direction_class"] = "not_significant"
            row["shortlist_gene"] = True
            row["blacklist_gene"] = False
            row["blacklist_exclusion_reason"] = ""
            row["candidate_reason"] = "Manual shortlist gene; no significant evidence row was found in this comparison."
            rows.append(row)
        placeholder = pd.DataFrame(rows)
        return pd.concat([matched, placeholder], axis=0, ignore_index=True)

    def _build_candidate_gene_table(
        self,
        integration: pd.DataFrame,
        mechanism: pd.DataFrame,
    ) -> pd.DataFrame:
        if self.current_project is None or integration is None or integration.empty:
            return pd.DataFrame()

        frame = integration.copy()
        frame["normalized_gene_id"] = self._normalized_gene_id_series(frame)
        support = mechanism.copy() if mechanism is not None else pd.DataFrame()
        if not support.empty:
            support = support.copy()
            support["normalized_gene_id"] = support.get("gene_id", pd.Series(index=support.index, dtype="object")).map(normalize_gene_id)
            join_columns = [column for column in ("comparison_id", "normalized_gene_id") if column in support.columns]
            if len(join_columns) == 2:
                frame = frame.merge(
                    support[
                        [
                            "comparison_id",
                            "normalized_gene_id",
                            "gene_id",
                            "gene_symbol",
                            "support_DEXSeq",
                            "support_DTU",
                            "support_SUPPA",
                            "best_DEXSeq_qvalue",
                            "n_DEXSeq_significant_exons",
                            "best_DTU_qvalue",
                            "n_DTU_significant_isoforms",
                            "n_support_methods",
                        ]
                    ],
                    on=["comparison_id", "normalized_gene_id"],
                    how="left",
                    suffixes=("", "_support"),
                )

        for column in ("support_DEXSeq", "support_DTU", "support_SUPPA", "n_support_methods"):
            if column not in frame.columns:
                frame[column] = 0
            frame[column] = pd.to_numeric(frame[column], errors="coerce").fillna(0)
        for column in ("best_DEXSeq_qvalue", "best_DTU_qvalue"):
            if column not in frame.columns:
                frame[column] = pd.NA
            frame[column] = pd.to_numeric(frame[column], errors="coerce")
        for column in ("n_DEXSeq_significant_exons", "n_DTU_significant_isoforms"):
            if column not in frame.columns:
                frame[column] = 0
            frame[column] = pd.to_numeric(frame[column], errors="coerce").fillna(0).astype(int)

        frame["standardized_log2FC"] = pd.to_numeric(frame.get("log2FC"), errors="coerce")
        frame["DE_FDR"] = pd.to_numeric(frame.get("deg_padj"), errors="coerce")
        frame["standardized_dPSI"] = pd.to_numeric(frame.get("representative_dPSI"), errors="coerce")
        frame["rMATS_FDR"] = pd.to_numeric(frame.get("representative_event_FDR"), errors="coerce")
        frame["abs_log2FC"] = frame["standardized_log2FC"].abs()
        frame["max_abs_dPSI"] = frame["standardized_dPSI"].abs()

        frame["DEG_significant"] = frame.get("deg_sig", False).fillna(False).astype(bool)
        frame["rMATS_significant"] = frame.get("splicing_sig", False).fillna(False).astype(bool)
        frame["DEXSeq_significant"] = frame["support_DEXSeq"] > 0
        frame["DTU_significant"] = frame["support_DTU"] > 0
        frame["n_rMATS_significant_events"] = pd.to_numeric(frame.get("n_sig_events_all"), errors="coerce").fillna(0).astype(int)
        frame["best_rMATS_FDR"] = frame["rMATS_FDR"]
        frame["dominant_rMATS_event_type"] = frame.get("representative_event_type")
        frame["best_rMATS_event_id"] = frame.get("representative_event_id")
        frame["dominant_rMATS_standardized_dPSI"] = frame["standardized_dPSI"]
        if "gene_symbol_support" in frame.columns:
            frame["gene_symbol_support"] = frame["gene_symbol_support"].where(
                frame["gene_symbol_support"].notna() & frame["gene_symbol_support"].astype(str).str.strip().ne("")
            )

        comparison_lookup = {item.comparison_id: item for item in self.current_project.available_comparisons}
        annotation_lookup = self._annotation_gene_symbol_lookup()
        normalized_gene_ids = frame["normalized_gene_id"]
        frame["group1"] = frame["comparison_id"].map(
            lambda value: self._comparison_group_labels(comparison_lookup.get(str(value)))[0]
        )
        frame["group2"] = frame["comparison_id"].map(
            lambda value: self._comparison_group_labels(comparison_lookup.get(str(value)))[1]
        )
        frame["rmats_group1"] = frame["comparison_id"].map(
            lambda value: self._comparison_splicing_direction_labels(comparison_lookup.get(str(value)))[0]
        )
        frame["rmats_group2"] = frame["comparison_id"].map(
            lambda value: self._comparison_splicing_direction_labels(comparison_lookup.get(str(value)))[1]
        )
        frame["deg_group1"] = frame["comparison_id"].map(
            lambda value: self._comparison_deg_direction_labels(comparison_lookup.get(str(value)))[0]
        )
        frame["deg_group2"] = frame["comparison_id"].map(
            lambda value: self._comparison_deg_direction_labels(comparison_lookup.get(str(value)))[1]
        )
        frame["log2fc_direction_label"] = frame.apply(
            lambda row: f"log2FC > 0 means transcript expression is higher in {row['deg_group1']} than {row['deg_group2']}"
            if pd.notna(row["deg_group1"]) and pd.notna(row["deg_group2"])
            else "log2FC > 0 follows the configured transcript numerator group",
            axis=1,
        )
        frame["dpsi_direction_label"] = frame.apply(
            lambda row: f"dPSI > 0 means event usage is higher in {row['rmats_group1']} than {row['rmats_group2']}"
            if pd.notna(row["rmats_group1"]) and pd.notna(row["rmats_group2"])
            else "dPSI > 0 follows the configured splicing group1",
            axis=1,
        )
        frame["dexseq_direction_label"] = frame.apply(
            lambda row: f"DEXSeq usage direction follows {row['rmats_group1']} vs {row['rmats_group2']}"
            if pd.notna(row["rmats_group1"]) and pd.notna(row["rmats_group2"])
            else "DEXSeq usage direction follows the configured splicing group1",
            axis=1,
        )
        frame["dtu_direction_label"] = frame.apply(
            lambda row: f"DTU isoform usage direction follows {row['rmats_group1']} vs {row['rmats_group2']}"
            if pd.notna(row["rmats_group1"]) and pd.notna(row["rmats_group2"])
            else "DTU isoform usage direction follows the configured splicing group1",
            axis=1,
        )
        frame["transcript_direction_flipped"] = frame["comparison_id"].map(
            lambda value: bool(comparison_lookup.get(str(value)).reverse_deg) if str(value) in comparison_lookup else False
        )
        frame["splicing_direction_flipped"] = frame["comparison_id"].map(
            lambda value: bool(comparison_lookup.get(str(value)).reverse_splicing) if str(value) in comparison_lookup else False
        )

        frame["evidence_count"] = (
            frame["DEG_significant"].astype(int)
            + frame["rMATS_significant"].astype(int)
            + frame["DEXSeq_significant"].astype(int)
            + frame["DTU_significant"].astype(int)
        )
        frame["evidence_class"] = frame.apply(self._candidate_evidence_class, axis=1)
        frame["direction_class"] = frame.apply(self._candidate_direction_class, axis=1)
        frame["candidate_tier"] = frame.apply(self._candidate_tier, axis=1)
        frame["candidate_reason"] = frame["candidate_tier"].map(self._candidate_reason)
        gene_symbol_base = frame.get("geneSymbol")
        if gene_symbol_base is None:
            gene_symbol_base = pd.Series(pd.NA, index=frame.index, dtype="object")
        frame["gene_symbol"] = gene_symbol_base
        if annotation_lookup:
            mapped_symbols = normalized_gene_ids.map(annotation_lookup)
            frame["gene_symbol"] = mapped_symbols.combine_first(frame["gene_symbol"])
        if "gene_symbol_support" in frame.columns:
            frame["gene_symbol"] = frame["gene_symbol"].combine_first(frame["gene_symbol_support"])
        frame["gene_symbol"] = frame["gene_symbol"].combine_first(frame.get("gene_id"))
        frame["shortlist_gene"] = frame["gene_symbol"].astype(str).isin(self._shortlist_gene_set())
        frame["blacklist_gene"] = frame["gene_symbol"].astype(str).isin(set(self.blacklist_genes_for_display()))
        frame["blacklist_exclusion_reason"] = frame["blacklist_gene"].map(lambda value: "manual blacklist" if value else "")

        frame = frame.sort_values(
            [
                "comparison_id",
                "candidate_tier",
                "evidence_count",
                "best_rMATS_FDR",
                "max_abs_dPSI",
                "DE_FDR",
                "abs_log2FC",
                "n_rMATS_significant_events",
                "gene_symbol",
            ],
            ascending=[True, True, False, True, False, True, False, False, True],
            na_position="last",
        ).reset_index(drop=True)
        frame["rank"] = frame.groupby("comparison_id").cumcount() + 1
        return frame[
            [
                "comparison_id",
                "comparison_name",
                "group1",
                "group2",
                "rank",
                "gene_id",
                "gene_symbol",
                "candidate_tier",
                "evidence_class",
                "evidence_count",
                "DEG_significant",
                "DE_FDR",
                "standardized_log2FC",
                "abs_log2FC",
                "log2fc_direction_label",
                "transcript_direction_flipped",
                "rMATS_significant",
                "n_rMATS_significant_events",
                "best_rMATS_FDR",
                "max_abs_dPSI",
                "dominant_rMATS_event_type",
                "best_rMATS_event_id",
                "dominant_rMATS_standardized_dPSI",
                "dpsi_direction_label",
                "dexseq_direction_label",
                "dtu_direction_label",
                "splicing_direction_flipped",
                "DEXSeq_significant",
                "best_DEXSeq_qvalue",
                "n_DEXSeq_significant_exons",
                "DTU_significant",
                "best_DTU_qvalue",
                "n_DTU_significant_isoforms",
                "direction_class",
                "shortlist_gene",
                "blacklist_gene",
                "blacklist_exclusion_reason",
                "candidate_reason",
            ]
        ].copy()

    def _build_cross_comparison_candidate_matrix(self, frame: pd.DataFrame) -> pd.DataFrame:
        if frame is None or frame.empty:
            return pd.DataFrame()
        working = frame.copy()
        if "blacklist_gene" in working.columns:
            working = working.loc[~working["blacklist_gene"].fillna(False)].copy()
        working["is_candidate"] = working["candidate_tier"].astype(str) != "Unclassified"
        working = working.loc[working["is_candidate"]].copy()
        if working.empty:
            return pd.DataFrame()
        rows: list[dict[str, object]] = []
        total_selected = len(working["comparison_id"].dropna().astype(str).unique().tolist())
        comparison_order = self._configured_comparison_order()
        for gene_id, subset in working.groupby("gene_id", dropna=False):
            row: dict[str, object] = {
                "gene_id": gene_id,
                "gene_symbol": subset["gene_symbol"].dropna().astype(str).iloc[0] if "gene_symbol" in subset.columns else gene_id,
            }
            candidate_comparisons: list[str] = []
            for _, item in subset.iterrows():
                comparison_id = str(item["comparison_id"])
                is_candidate = bool(item["is_candidate"])
                row[f"is_candidate_{comparison_id}"] = is_candidate
                row[f"tier_{comparison_id}"] = item["candidate_tier"]
                row[f"rank_{comparison_id}"] = item["rank"]
                row[f"evidence_class_{comparison_id}"] = item["evidence_class"]
                row[f"evidence_count_{comparison_id}"] = item["evidence_count"]
                row[f"direction_class_{comparison_id}"] = item["direction_class"]
                row[f"DE_FDR_{comparison_id}"] = item["DE_FDR"]
                row[f"standardized_log2FC_{comparison_id}"] = item["standardized_log2FC"]
                row[f"best_rMATS_FDR_{comparison_id}"] = item["best_rMATS_FDR"]
                row[f"dominant_rMATS_standardized_dPSI_{comparison_id}"] = item["dominant_rMATS_standardized_dPSI"]
                row[f"max_abs_dPSI_{comparison_id}"] = item["max_abs_dPSI"]
                row[f"dominant_rMATS_event_type_{comparison_id}"] = item["dominant_rMATS_event_type"]
                if is_candidate:
                    candidate_comparisons.append(comparison_id)
            candidate_pattern, pattern_reason, has_direction_reversal, has_tier_change, has_evidence_class_change = self._cross_comparison_pattern(
                subset,
                candidate_comparisons,
                total_selected,
                comparison_order,
            )
            row["candidate_pattern"] = candidate_pattern
            row["pattern_reason"] = pattern_reason
            row["candidate_comparisons"] = ";".join(candidate_comparisons)
            row["n_comparisons_candidate"] = len(candidate_comparisons)
            row["best_overall_tier"] = subset["candidate_tier"].astype(str).sort_values().iloc[0]
            row["best_overall_rank"] = pd.to_numeric(subset["rank"], errors="coerce").min()
            row["has_direction_reversal"] = has_direction_reversal
            row["has_tier_change"] = has_tier_change
            row["has_evidence_class_change"] = has_evidence_class_change
            rows.append(row)
        return pd.DataFrame(rows)

    @staticmethod
    def _empty_cross_as_pattern_tables() -> dict[str, pd.DataFrame]:
        base_columns = [
            "gene",
            "event_type",
            "event_id",
            "event_matching_key",
            "comparison_A",
            "comparison_B",
            "experiment_A_display_name",
            "experiment_B_display_name",
            "dPSI_A",
            "dPSI_B",
            "FDR_A",
            "FDR_B",
            "significant_A",
            "significant_B",
            "direction_A",
            "direction_B",
            "delta_dPSI",
            "abs_delta_dPSI",
            "same_direction",
            "opposite_direction",
            "same_direction_strength_delta",
            "event_class",
            "abs_dPSI_significance_cutoff_used",
            "large_delta_dPSI_cutoff_used",
        ]
        empty = pd.DataFrame(columns=base_columns)
        return {
            "same_direction_large_delta": empty.copy(),
            "opposite_direction": empty.copy(),
            "all_shared": empty.copy(),
            "condition_specific": empty.copy(),
        }

    def _cross_as_pattern_source_frame(self, comparison_id: str | None, abs_dpsi_cutoff: float) -> pd.DataFrame:
        if not comparison_id or self.current_project is None:
            return pd.DataFrame()
        comparison = next(
            (item for item in self.current_project.available_comparisons if item.comparison_id == comparison_id),
            None,
        )
        if comparison is None:
            return pd.DataFrame()
        events = self._sashimi_event_catalog_frame(comparison_id)
        if events is None or events.empty:
            return pd.DataFrame()

        working = events.copy()
        working["gene"] = working.get("gene_symbol", pd.Series(pd.NA, index=working.index)).fillna(working.get("gene_id"))
        working["event_id"] = working.get("event_id", pd.Series("", index=working.index)).fillna("").astype(str)
        working["coordinates"] = working.get("coordinates", pd.Series("", index=working.index)).fillna("").astype(str)
        working["event_id"] = working["event_id"].where(working["event_id"].str.strip().ne(""), working["coordinates"])
        working["normalized_gene_id"] = self._normalized_gene_id_series(working)
        working["event_matching_key"] = working.apply(
            lambda row: "|".join(
                [
                    str(row.get("normalized_gene_id") or row.get("gene") or ""),
                    str(row.get("event_type") or ""),
                    str(row.get("event_id") or row.get("coordinates") or ""),
                ]
            ),
            axis=1,
        )
        working["dPSI"] = pd.to_numeric(working.get("dPSI"), errors="coerce")
        working["FDR"] = pd.to_numeric(working.get("FDR"), errors="coerce")
        left, right = self._comparison_splicing_direction_labels(comparison)
        working["comparison_id"] = comparison.comparison_id
        working["comparison_display_name"] = comparison.display_resolved_name
        working["direction"] = working["dPSI"].map(
            lambda value: (
                f"inclusion higher in {left} than {right}"
                if pd.notna(value) and float(value) > 0
                else f"inclusion lower in {left} than {right}"
                if pd.notna(value) and float(value) < 0
                else "no directional shift"
            )
        )
        working["significant"] = (
            working["FDR"].le(0.05)
            & working["dPSI"].abs().ge(abs_dpsi_cutoff)
        )
        working = working.sort_values(["significant", "FDR"], ascending=[False, True], na_position="last")
        working = working.drop_duplicates(subset=["event_matching_key"], keep="first").reset_index(drop=True)
        return working[
            [
                "comparison_id",
                "comparison_display_name",
                "gene",
                "gene_symbol",
                "gene_id",
                "event_type",
                "event_id",
                "event_matching_key",
                "coordinates",
                "dPSI",
                "FDR",
                "significant",
                "direction",
            ]
        ].copy()

    def _build_cross_comparison_significant_as_patterns(
        self,
        pair: ComparisonPairDefinition,
        *,
        abs_dpsi_cutoff: float,
        large_delta_dpsi_cutoff: float,
    ) -> dict[str, pd.DataFrame]:
        comparison_a = self._cross_as_pattern_source_frame(pair.comparison_a, abs_dpsi_cutoff)
        comparison_b = self._cross_as_pattern_source_frame(pair.comparison_b, abs_dpsi_cutoff)
        if comparison_a.empty and comparison_b.empty:
            return self._empty_cross_as_pattern_tables()

        merged = comparison_a.merge(
            comparison_b,
            on="event_matching_key",
            how="outer",
            suffixes=("_A", "_B"),
        )
        if merged.empty:
            return self._empty_cross_as_pattern_tables()

        merged["gene"] = (
            merged.get("gene_A", pd.Series(pd.NA, index=merged.index))
            .combine_first(merged.get("gene_B", pd.Series(pd.NA, index=merged.index)))
            .fillna("")
            .astype(str)
        )
        merged["event_type"] = (
            merged.get("event_type_A", pd.Series(pd.NA, index=merged.index))
            .combine_first(merged.get("event_type_B", pd.Series(pd.NA, index=merged.index)))
            .fillna("")
            .astype(str)
        )
        merged["event_id"] = (
            merged.get("event_id_A", pd.Series(pd.NA, index=merged.index))
            .combine_first(merged.get("event_id_B", pd.Series(pd.NA, index=merged.index)))
            .fillna("")
            .astype(str)
        )
        merged["comparison_A"] = pair.comparison_a or ""
        merged["comparison_B"] = pair.comparison_b or ""
        merged["experiment_A_display_name"] = comparison_a.get("comparison_display_name", pd.Series(dtype="object")).iloc[0] if not comparison_a.empty else str(pair.comparison_a or "")
        merged["experiment_B_display_name"] = comparison_b.get("comparison_display_name", pd.Series(dtype="object")).iloc[0] if not comparison_b.empty else str(pair.comparison_b or "")
        merged["dPSI_A"] = pd.to_numeric(merged.get("dPSI_A"), errors="coerce")
        merged["dPSI_B"] = pd.to_numeric(merged.get("dPSI_B"), errors="coerce")
        merged["FDR_A"] = pd.to_numeric(merged.get("FDR_A"), errors="coerce")
        merged["FDR_B"] = pd.to_numeric(merged.get("FDR_B"), errors="coerce")
        merged["significant_A"] = merged.get("significant_A", pd.Series(False, index=merged.index)).fillna(False).astype(bool)
        merged["significant_B"] = merged.get("significant_B", pd.Series(False, index=merged.index)).fillna(False).astype(bool)
        merged["direction_A"] = merged.get("direction_A", pd.Series("", index=merged.index)).fillna("").astype(str)
        merged["direction_B"] = merged.get("direction_B", pd.Series("", index=merged.index)).fillna("").astype(str)
        merged["delta_dPSI"] = merged["dPSI_A"] - merged["dPSI_B"]
        merged["abs_delta_dPSI"] = merged["delta_dPSI"].abs()
        sign_a = np.sign(merged["dPSI_A"].fillna(0))
        sign_b = np.sign(merged["dPSI_B"].fillna(0))
        merged["same_direction"] = (sign_a == sign_b) & sign_a.ne(0) & sign_b.ne(0)
        merged["opposite_direction"] = sign_a.ne(sign_b) & sign_a.ne(0) & sign_b.ne(0)
        merged["same_direction_strength_delta"] = (merged["dPSI_A"].abs() - merged["dPSI_B"].abs()).abs()
        merged["abs_dPSI_significance_cutoff_used"] = abs_dpsi_cutoff
        merged["large_delta_dPSI_cutoff_used"] = large_delta_dpsi_cutoff

        shared = merged.loc[merged["significant_A"] & merged["significant_B"]].copy()
        shared["event_class"] = "shared_significant"
        same_direction_large_delta = shared.loc[
            shared["same_direction"] & shared["same_direction_strength_delta"].ge(large_delta_dpsi_cutoff)
        ].copy()
        same_direction_large_delta["event_class"] = "same_direction_large_delta"
        opposite_direction = shared.loc[shared["opposite_direction"]].copy()
        opposite_direction["event_class"] = "opposite_direction"

        condition_specific = merged.loc[
            (merged["significant_A"] & ~merged["significant_B"])
            | (~merged["significant_A"] & merged["significant_B"])
        ].copy()
        condition_specific["event_class"] = np.where(
            condition_specific["significant_A"],
            "comparison_A_only",
            "comparison_B_only",
        )

        output_columns = [
            "gene",
            "event_type",
            "event_id",
            "event_matching_key",
            "comparison_A",
            "comparison_B",
            "experiment_A_display_name",
            "experiment_B_display_name",
            "dPSI_A",
            "dPSI_B",
            "FDR_A",
            "FDR_B",
            "significant_A",
            "significant_B",
            "direction_A",
            "direction_B",
            "delta_dPSI",
            "abs_delta_dPSI",
            "same_direction",
            "opposite_direction",
            "same_direction_strength_delta",
            "event_class",
            "abs_dPSI_significance_cutoff_used",
            "large_delta_dPSI_cutoff_used",
        ]

        return {
            "same_direction_large_delta": same_direction_large_delta.reindex(columns=output_columns).reset_index(drop=True),
            "opposite_direction": opposite_direction.reindex(columns=output_columns).reset_index(drop=True),
            "all_shared": shared.reindex(columns=output_columns).reset_index(drop=True),
            "condition_specific": condition_specific.reindex(columns=output_columns).reset_index(drop=True),
        }

    def _cross_comparison_pattern(
        self,
        subset: pd.DataFrame,
        candidate_comparisons: list[str],
        total_selected: int,
        comparison_order: list[str],
    ) -> tuple[str, str, bool, bool, bool]:
        if not candidate_comparisons:
            return "comparison_specific", "Candidate present in only one comparison row.", False, False, False
        evidence_classes = set(subset["evidence_class"].dropna().astype(str).tolist()) if "evidence_class" in subset.columns else set()
        direction_classes = set(subset["direction_class"].dropna().astype(str).tolist()) if "direction_class" in subset.columns else set()
        has_direction_reversal = self._has_direction_reversal(direction_classes)
        has_tier_change = len(set(subset["candidate_tier"].dropna().astype(str).tolist())) > 1 if "candidate_tier" in subset.columns else False
        has_evidence_class_change = len(evidence_classes) > 1
        if len(candidate_comparisons) == total_selected and total_selected > 0:
            return "shared_all", "Candidate is present in all selected comparisons.", has_direction_reversal, has_tier_change, has_evidence_class_change
        if len(candidate_comparisons) > 1:
            if has_direction_reversal:
                return "direction_reversed", "Standardized transcript or splicing direction reverses across comparisons.", has_direction_reversal, has_tier_change, has_evidence_class_change
            if has_tier_change:
                return "tier_changed", "Candidate tier changes across comparisons.", has_direction_reversal, has_tier_change, has_evidence_class_change
            if has_evidence_class_change:
                return "evidence_class_changed", "Evidence class changes across comparisons.", has_direction_reversal, has_tier_change, has_evidence_class_change
            transition_pattern = self._cross_comparison_transition_pattern(candidate_comparisons, comparison_order)
            if transition_pattern == "gained":
                return "gained", "Candidate becomes present later in the defined comparison order.", has_direction_reversal, has_tier_change, has_evidence_class_change
            if transition_pattern == "lost":
                return "lost", "Candidate is present early and absent later in the defined comparison order.", has_direction_reversal, has_tier_change, has_evidence_class_change
            if transition_pattern == "order_not_defined":
                return "shared_subset", "order_not_defined: configure comparison order to distinguish gained vs lost patterns.", has_direction_reversal, has_tier_change, has_evidence_class_change
            return "shared_subset", "Candidate is shared across a subset of comparisons.", has_direction_reversal, has_tier_change, has_evidence_class_change
        if evidence_classes == {"DEG"} or evidence_classes == {"DEG_only"}:
            return "DEG_only_in_one_condition", "Candidate is DEG-only in one comparison.", has_direction_reversal, has_tier_change, has_evidence_class_change
        if all(label.startswith("rMATS") or label == "rMATS" for label in evidence_classes if label):
            return "splicing_only_in_one_condition", "Candidate is splicing-only in one comparison.", has_direction_reversal, has_tier_change, has_evidence_class_change
        return "comparison_specific", "Candidate is specific to one comparison.", has_direction_reversal, has_tier_change, has_evidence_class_change

    @staticmethod
    def _has_direction_reversal(direction_classes: set[str]) -> bool:
        if not direction_classes:
            return False
        transcript_up = any("transcript_up" in label for label in direction_classes)
        transcript_down = any("transcript_down" in label for label in direction_classes)
        splicing_up = any(label.endswith("splicing_up") for label in direction_classes)
        splicing_down = any(label.endswith("splicing_down") for label in direction_classes)
        return (transcript_up and transcript_down) or (splicing_up and splicing_down)

    def _configured_comparison_order(self) -> list[str]:
        if self.current_project is None:
            return []
        selected = {item.comparison_id for item in self._selected_comparisons()}
        configured = [item for item in self.current_project.comparison_order_ids if item in selected]
        return configured

    @staticmethod
    def _cross_comparison_transition_pattern(candidate_comparisons: list[str], comparison_order: list[str]) -> str:
        if len(candidate_comparisons) <= 1:
            return "comparison_specific"
        if not comparison_order:
            return "order_not_defined"
        candidate_set = set(candidate_comparisons)
        ordered = list(comparison_order)
        if not ordered:
            return "order_not_defined"
        presence = [comparison_id in candidate_set for comparison_id in ordered]
        gains = any((not previous) and current for previous, current in zip(presence, presence[1:], strict=False))
        losses = any(previous and (not current) for previous, current in zip(presence, presence[1:], strict=False))
        if gains and not losses:
            return "gained"
        if losses and not gains:
            return "lost"
        return "shared_subset"

    def _comparison_group_labels(self, comparison) -> tuple[str, str]:
        if comparison is None:
            return "group1", "group2"
        return (
            comparison.experiment_group or comparison.source_experiment_group or "group1",
            comparison.control_group or comparison.source_control_group or "group2",
        )

    def _comparison_splicing_direction_labels(self, comparison) -> tuple[str, str]:
        if comparison is None:
            return "group1", "group2"
        left, right = comparison.final_as_groups
        return (
            left or comparison.experiment_group or comparison.source_experiment_group or "group1",
            right or comparison.control_group or comparison.source_control_group or "group2",
        )

    def _comparison_deg_direction_labels(self, comparison) -> tuple[str, str]:
        if comparison is None:
            return "group1", "group2"
        left, right = comparison.final_deg_groups
        return (
            left or comparison.experiment_group or comparison.source_experiment_group or "group1",
            right or comparison.control_group or comparison.source_control_group or "group2",
        )

    def _candidate_evidence_class(self, row: pd.Series) -> str:
        labels: list[str] = []
        if bool(row.get("DEG_significant")):
            labels.append("DEG")
        if bool(row.get("rMATS_significant")):
            labels.append("rMATS")
        if bool(row.get("DEXSeq_significant")):
            labels.append("DEXSeq")
        if bool(row.get("DTU_significant")):
            labels.append("DTU")
        return "_".join(labels) if labels else "unclassified"

    def _candidate_direction_class(self, row: pd.Series) -> str:
        deg_sig = bool(row.get("DEG_significant"))
        rmats_sig = bool(row.get("rMATS_significant"))
        log2fc = pd.to_numeric(row.get("standardized_log2FC"), errors="coerce")
        dpsi = pd.to_numeric(row.get("dominant_rMATS_standardized_dPSI"), errors="coerce")
        if deg_sig and rmats_sig and pd.notna(log2fc) and pd.notna(dpsi):
            transcript_up = float(log2fc) > 0
            splicing_up = float(dpsi) > 0
            if transcript_up and splicing_up:
                return "transcript_up_splicing_up"
            if transcript_up and not splicing_up:
                return "transcript_up_splicing_down"
            if not transcript_up and splicing_up:
                return "transcript_down_splicing_up"
            return "transcript_down_splicing_down"
        if deg_sig and not rmats_sig:
            return "DEG_only"
        if rmats_sig and not deg_sig and (bool(row.get("DEXSeq_significant")) or bool(row.get("DTU_significant"))):
            return "rMATS_DEXSeq_or_DTU_only"
        if rmats_sig and not deg_sig:
            return "rMATS_only"
        return "unclassified"

    def _candidate_tier(self, row: pd.Series) -> str:
        deg_sig = bool(row.get("DEG_significant"))
        rmats_sig = bool(row.get("rMATS_significant"))
        dex_sig = bool(row.get("DEXSeq_significant"))
        dtu_sig = bool(row.get("DTU_significant"))
        support = dex_sig or dtu_sig
        if deg_sig and rmats_sig and support:
            return "Tier 1"
        if deg_sig and rmats_sig and not support:
            return "Tier 2"
        if (not deg_sig) and rmats_sig and support:
            return "Tier 3"
        if (not deg_sig) and rmats_sig and not support:
            return "Tier 4"
        if deg_sig and (not rmats_sig) and (not support):
            return "Tier 5"
        return "Unclassified"

    def _candidate_reason(self, tier: str) -> str:
        mapping = {
            "Tier 1": "Tier 1: DEG significant, rMATS significant, and DEXSeq/DTU support. Strong multi-layer candidate.",
            "Tier 2": "Tier 2: DEG significant and rMATS significant. Core expression-splicing candidate.",
            "Tier 3": "Tier 3: rMATS significant with DEXSeq/DTU support but no DEG. Splicing-primary candidate.",
            "Tier 4": "Tier 4: rMATS significant only. Event-level splicing candidate.",
            "Tier 5": "Tier 5: DEG significant only. Expression-only candidate.",
        }
        return mapping.get(str(tier), "Unclassified candidate.")

    def _shortlist_gene_set(self) -> set[str]:
        values: set[str] = set()
        if self.current_project is not None:
            values.update(item.strip() for item in self.current_project.shortlist_genes if item and item.strip())
        shortlist = self.run_state.results.cards_shortlist
        if shortlist is not None and not shortlist.empty:
            for column in ("geneSymbol", "GeneID", "gene_symbol", "gene_id"):
                if column in shortlist.columns:
                    values.update(
                        str(value).strip()
                        for value in shortlist[column].dropna().astype(str).tolist()
                        if str(value).strip()
                    )
        return values

    def _shortlist_output_tables(self) -> dict[str, pd.DataFrame]:
        tables: dict[str, pd.DataFrame] = {
            "ranked_candidates.tsv": self.run_state.results.ranked_candidates,
        }
        shortlist_frames = [
            frame
            for frame in (self.run_state.results.shortlist_dp, self.run_state.results.shortlist_ko)
            if frame is not None and not frame.empty
        ]
        for index, frame in enumerate(shortlist_frames, start=1):
            display_group = (
                str(frame["display_group"].iloc[0])
                if "display_group" in frame.columns and not frame.empty
                else f"shortlist_{index}"
            )
            tables[f"{self._safe_output_name(display_group)}.tsv"] = frame
        if self.run_state.results.cards_shortlist is not None and not self.run_state.results.cards_shortlist.empty:
            tables["combined_shortlist.tsv"] = self.run_state.results.cards_shortlist
        return tables

    def _candidate_output_tables(self) -> dict[str, pd.DataFrame]:
        frame = self.run_state.results.candidate_gene_table
        if frame is None or frame.empty:
            return {}
        tables: dict[str, pd.DataFrame] = {
            "gene_level_integrated_candidates.tsv": frame,
        }
        for tier in sorted(frame["candidate_tier"].dropna().astype(str).unique().tolist()):
            subset = frame.loc[frame["candidate_tier"].astype(str) == tier].copy()
            if subset.empty:
                continue
            tables[f"{self._safe_output_name(tier.lower())}_candidates.tsv"] = subset
        matrix = self.run_state.results.cross_comparison_candidate_matrix
        if matrix is not None and not matrix.empty:
            tables["cross_comparison_candidate_matrix.tsv"] = matrix
        return tables

    def _comparison_config_payload(self) -> dict[str, object]:
        if self.current_project is None:
            return {}
        payload: dict[str, object] = {}
        for comparison in self.current_project.available_comparisons:
            group1, group2 = self._comparison_group_labels(comparison)
            numerator_group = comparison.experiment_group or comparison.source_experiment_group or group1
            denominator_group = comparison.control_group or comparison.source_control_group or group2
            output_prefix = comparison.output_prefix or self._safe_output_name(comparison.comparison_id)
            payload[comparison.comparison_id] = {
                "display_name": comparison.display_resolved_name,
                "biological_question": comparison.biological_question or "",
                "group1": group1,
                "group2": group2,
                "transcript_contrast": {
                    "numerator_group": numerator_group,
                    "denominator_group": denominator_group,
                    "log2fc_positive_means": f"Transcript expression higher in {numerator_group} than {denominator_group}",
                },
                "splicing_contrast": {
                    "group1": group1,
                    "group2": group2,
                    "dpsi_positive_means": f"Event usage higher in {group1} than {group2}",
                },
                "expected_direction_notes": {
                    "notes": comparison.expected_direction_notes or "",
                    "known_positive_control_genes": list(comparison.known_positive_control_genes),
                    "known_negative_control_genes": list(comparison.known_negative_control_genes),
                },
                "output_prefix": output_prefix,
            }
        return payload

    def _comparison_parameters_payload(self, comparison) -> dict[str, object]:
        group1, group2 = self._comparison_group_labels(comparison)
        numerator_group = comparison.experiment_group or comparison.source_experiment_group or group1
        denominator_group = comparison.control_group or comparison.source_control_group or group2
        return {
            "comparison_id": comparison.comparison_id,
            "display_name": comparison.display_resolved_name,
            "biological_question": comparison.biological_question or "",
            "group1": group1,
            "group2": group2,
            "transcript_contrast": {
                "numerator_group": numerator_group,
                "denominator_group": denominator_group,
            },
            "splicing_contrast": {
                "group1": group1,
                "group2": group2,
            },
            "log2fc_direction_label": f"log2FC > 0 means transcript expression is higher in {numerator_group} than {denominator_group}",
            "dpsi_direction_label": f"dPSI > 0 means event usage is higher in {group1} than {group2}",
            "transcript_direction_flipped": bool(comparison.reverse_deg),
            "splicing_direction_flipped": bool(comparison.reverse_splicing),
            "thresholds": {
                "DE_FDR": self.current_thresholds.deg_padj,
                "abs_standardized_log2FC": self.current_thresholds.deg_log2fc,
                "rMATS_FDR": self.current_thresholds.splicing_fdr,
                "abs_standardized_dPSI": self.current_thresholds.splicing_dpsi,
                "DEXSeq_qvalue": self.current_thresholds.dexseq_qvalue,
                "DTU_qvalue": self.current_thresholds.dtu_qvalue,
                "program_delta_dpsi": self.current_thresholds.program_delta_dpsi,
            },
            "input_paths": dict(self.current_project.input_paths) if self.current_project is not None else {},
            "output_prefix": comparison.output_prefix or self._safe_output_name(comparison.comparison_id),
            "shortlist": list(self._shortlist_gene_set()),
            "blacklist": list(self.blacklist_genes_for_display()),
        }

    def _write_json(self, path: Path, payload: dict[str, object]) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2, ensure_ascii=False)

    def _write_text(self, path: Path, text: str) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text, encoding="utf-8")

    def _write_per_comparison_candidate_outputs(self) -> None:
        if self.current_project is None or self.current_project.output_root is None:
            return
        frame = self.run_state.results.candidate_gene_table
        if frame is None or frame.empty:
            return
        base_dir = self.current_project.output_root / "02_candidate_gene_screening"
        for comparison in self._selected_comparisons():
            subset = frame.loc[frame["comparison_id"].astype(str) == comparison.comparison_id].copy()
            if subset.empty:
                continue
            comp_dir = base_dir / (comparison.output_prefix or self._safe_output_name(comparison.comparison_id))
            comp_dir.mkdir(parents=True, exist_ok=True)
            subset.to_csv(comp_dir / "gene_level_integrated_candidates.tsv", sep="\t", index=False)
            for tier in ("Tier 1", "Tier 2", "Tier 3", "Tier 4", "Tier 5"):
                tier_subset = subset.loc[subset["candidate_tier"].astype(str) == tier].copy()
                tier_subset.to_csv(comp_dir / f"{self._safe_output_name(tier.lower())}_candidates.tsv", sep="\t", index=False)
            subset.loc[subset["DEG_significant"]].to_csv(comp_dir / "DEG_significant_genes.tsv", sep="\t", index=False)
            subset.loc[subset["rMATS_significant"]].to_csv(comp_dir / "rMATS_gene_summary.tsv", sep="\t", index=False)
            subset.loc[subset["DEXSeq_significant"]].to_csv(comp_dir / "DEXSeq_gene_summary.tsv", sep="\t", index=False)
            subset.loc[subset["DTU_significant"]].to_csv(comp_dir / "DTU_gene_summary.tsv", sep="\t", index=False)
            subset.loc[subset["shortlist_gene"]].to_csv(comp_dir / "shortlist.tsv", sep="\t", index=False)
            subset.loc[subset["blacklist_gene"]].to_csv(comp_dir / "blacklist.tsv", sep="\t", index=False)
            summary = pd.DataFrame(
                [
                    {
                        "comparison_id": comparison.comparison_id,
                        "display_name": comparison.display_resolved_name,
                        "n_genes": len(subset),
                        "n_tier1": int((subset["candidate_tier"] == "Tier 1").sum()),
                        "n_tier2": int((subset["candidate_tier"] == "Tier 2").sum()),
                        "n_tier3": int((subset["candidate_tier"] == "Tier 3").sum()),
                        "n_tier4": int((subset["candidate_tier"] == "Tier 4").sum()),
                        "n_tier5": int((subset["candidate_tier"] == "Tier 5").sum()),
                        "n_shortlist": int(subset["shortlist_gene"].fillna(False).sum()),
                        "n_blacklist": int(subset["blacklist_gene"].fillna(False).sum()),
                    }
                ]
            )
            summary.to_csv(comp_dir / "summary.tsv", sep="\t", index=False)
            self._write_json(comp_dir / "parameters.json", self._comparison_parameters_payload(comparison))
            self._write_text(
                comp_dir / "run.log",
                f"Generated candidate screening outputs for {comparison.comparison_id}\nRows: {len(subset)}\n",
            )

    def _write_cross_comparison_outputs(self) -> None:
        if self.current_project is None or self.current_project.output_root is None:
            return
        matrix = self.run_state.results.cross_comparison_candidate_matrix
        if matrix is None or matrix.empty:
            return
        base_dir = self.current_project.output_root / "03_cross_comparison_candidate_comparison"
        base_dir.mkdir(parents=True, exist_ok=True)
        matrix.to_csv(base_dir / "cross_comparison_candidate_matrix.tsv", sep="\t", index=False)
        summary = pd.DataFrame(
            [
                {
                    "n_genes": len(matrix),
                    "n_shared_all": int((matrix["candidate_pattern"].astype(str) == "shared_all").sum()) if "candidate_pattern" in matrix.columns else 0,
                    "n_multi_comparison": int(pd.to_numeric(matrix.get("n_comparisons_candidate"), errors="coerce").fillna(0).ge(2).sum()) if "n_comparisons_candidate" in matrix.columns else 0,
                }
            ]
        )
        summary.to_csv(base_dir / "cross_comparison_summary.tsv", sep="\t", index=False)
        self._write_json(
            base_dir / "parameters.json",
            {
                "selected_comparisons": [item.comparison_id for item in self._selected_comparisons()],
                "thresholds": {
                    "DE_FDR": self.current_thresholds.deg_padj,
                    "abs_standardized_log2FC": self.current_thresholds.deg_log2fc,
                    "rMATS_FDR": self.current_thresholds.splicing_fdr,
                    "abs_standardized_dPSI": self.current_thresholds.splicing_dpsi,
                    "DEXSeq_qvalue": self.current_thresholds.dexseq_qvalue,
                    "DTU_qvalue": self.current_thresholds.dtu_qvalue,
                    "program_delta_dpsi": self.current_thresholds.program_delta_dpsi,
                },
                "shortlist": list(self._shortlist_gene_set()),
                "blacklist": list(self.blacklist_genes_for_display()),
            },
        )
        self._write_text(base_dir / "run.log", f"Generated cross-comparison candidate matrix\nRows: {len(matrix)}\n")

    def _filter_sashimi_manifest(self, manifest: pd.DataFrame, shortlist: pd.DataFrame) -> pd.DataFrame:
        if manifest.empty or shortlist.empty:
            return manifest

        filtered = manifest.copy()
        selected = self._selected_comparisons()
        allowed_comparison_names = {
            comparison.source_name("sashimi")
            for comparison in selected
        } | {
            comparison.resolved_name
            for comparison in selected
        } | {
            comparison.comparison_id
            for comparison in selected
        }
        if "comparison_id" in filtered.columns and allowed_comparison_names:
            filtered = filtered[
                filtered["comparison_id"].astype(str).isin(allowed_comparison_names)
            ].copy()

        if "event_uid" in filtered.columns and "event_uid" in shortlist.columns:
            event_ids = {
                str(value)
                for value in shortlist["event_uid"].dropna().astype(str).tolist()
                if str(value).strip()
            }
            if event_ids:
                filtered = filtered[filtered["event_uid"].astype(str).isin(event_ids)].copy()

        if filtered.empty:
            filtered = manifest.copy()
            gene_candidates: set[str] = set()
            for column in ("GeneID", "gene_id", "geneSymbol", "gene_symbol"):
                if column in shortlist.columns:
                    gene_candidates.update(
                        str(value)
                        for value in shortlist[column].dropna().astype(str).tolist()
                        if str(value).strip()
                    )
            gene_columns = [column for column in ("GeneID", "geneSymbol") if column in filtered.columns]
            if gene_candidates and gene_columns:
                mask = pd.Series(False, index=filtered.index)
                for column in gene_columns:
                    mask = mask | filtered[column].astype(str).isin(gene_candidates)
                filtered = filtered[mask].copy()

        return filtered.reset_index(drop=True)

    def get_or_build_sashimi_manifest(self, allow_generate: bool = True) -> pd.DataFrame:
        if self.current_project is None:
            raise RuntimeError("No project loaded.")
        self._log_heavy("get_or_build_sashimi_manifest triggered")
        manifest = self.run_state.results.sashimi_manifest
        if manifest is not None and not manifest.empty:
            self._mark_module_state("sashimi_manifest", "cache_hit", "Loaded cached sashimi manifest.", from_cache=True)
            return manifest.copy()
        if not allow_generate:
            return pd.DataFrame()
        shortlist = self._candidate_selection_sashimi_source()
        if shortlist.empty:
            self._mark_module_state(
                "sashimi_manifest",
                "not_run",
                "Missing input: event follow-up source not available. Please run 5.3.2 / 5.3.4 first.",
            )
            return pd.DataFrame()
        raw_manifest = self.load_first_matching_table("rmats2sashimi_manifest.tsv")
        filtered = self._filter_sashimi_manifest(raw_manifest, shortlist)
        if filtered.empty:
            filtered = self._build_generated_sashimi_manifest(shortlist)
        if filtered is None or filtered.empty:
            filtered = pd.DataFrame()
        else:
            filtered = filtered.reset_index(drop=True)
        self.run_state.results.sashimi_manifest = filtered
        self._mark_module_state(
            "sashimi_manifest",
            "finished",
            "Generated sashimi manifest from candidate event follow-up inputs.",
            from_cache=False,
            output_folder=str(self.current_project.output_root / "06_sashimi") if self.current_project and self.current_project.output_root else "",
        )
        return filtered.copy()

    def _candidate_selection_sashimi_source(self) -> pd.DataFrame:
        if self.current_project is None:
            return pd.DataFrame()
        frames: list[pd.DataFrame] = []
        for comparison in self._selected_comparisons():
            selected, _meta = self.candidate_event_followup_frame(
                comparison.comparison_id,
                source="selection",
                top_n=20,
                allow_generate=False,
            )
            if selected.empty:
                continue
            subset = selected.copy()
            subset = subset.loc[subset["best_rMATS_event_id"].fillna("").astype(str).str.strip().ne("")]
            if subset.empty:
                continue
            subset["event_uid"] = subset["best_rMATS_event_id"].astype(str)
            subset["event_type"] = subset["dominant_rMATS_event_type"].astype(str)
            subset["GeneID"] = subset["gene_id"]
            subset["geneSymbol"] = subset["gene_symbol"]
            subset["display_group"] = comparison.display_resolved_name
            subset["class_interpretation"] = subset["direction_class"].astype(str)
            frames.append(subset)
        if not frames:
            return pd.DataFrame()
        return pd.concat(frames, ignore_index=True)

    def _build_generated_sashimi_manifest(self, shortlist: pd.DataFrame) -> pd.DataFrame:
        if self.current_project is None or self.current_project.output_root is None or shortlist.empty:
            return pd.DataFrame()

        generated_root = self.current_project.output_root / "06_sashimi" / "data" / "generated_manifest"
        event_root = generated_root / "event_files"
        bam_root = generated_root / "bam_lists"
        event_root.mkdir(parents=True, exist_ok=True)
        bam_root.mkdir(parents=True, exist_ok=True)

        selected = self._selected_comparisons() or self.current_project.available_comparisons
        comparison_lookup = {item.comparison_id: item for item in selected}
        mode = self.current_project.selected_rmats_mode or "JC"
        table_cache: dict[tuple[str, str, str], pd.DataFrame] = {}
        rows: list[dict[str, object]] = []
        seen: set[tuple[str, str]] = set()

        for _, row in shortlist.iterrows():
            event_uid = str(row.get("event_uid", "")).strip()
            if not event_uid:
                continue
            event_type = (
                str(row.get("event_type", "")).strip()
                or str(row.get("event_type_A", "")).strip()
                or str(row.get("event_type_B", "")).strip()
                or event_uid.split("|", 1)[0]
            )
            if not event_type:
                continue
            for comparison_ref in self._sashimi_target_comparisons(row):
                comparison_id = self._resolve_comparison_reference(comparison_ref)
                if not comparison_id:
                    continue
                comparison = comparison_lookup.get(comparison_id)
                if comparison is None:
                    continue
                key = (event_uid, comparison_id)
                if key in seen:
                    continue
                event_file = self._write_generated_sashimi_event_file(
                    comparison,
                    event_uid,
                    event_type,
                    mode,
                    event_root,
                    table_cache,
                )
                if event_file is None:
                    continue
                bam_inputs = self._prepare_generated_sashimi_bam_lists(comparison, bam_root)
                if bam_inputs is None:
                    continue
                b1_txt, b2_txt, label1, label2 = bam_inputs
                gene_symbol = (
                    str(row.get("geneSymbol", "")).strip()
                    or str(row.get("gene_symbol", "")).strip()
                    or str(row.get("gene_symbol_A", "")).strip()
                    or str(row.get("gene_symbol_B", "")).strip()
                )
                gene_id = (
                    str(row.get("GeneID", "")).strip()
                    or str(row.get("gene_id", "")).strip()
                    or str(row.get("gene_id_A", "")).strip()
                    or str(row.get("gene_id_B", "")).strip()
                )
                display_group = str(row.get("display_group", "")).strip() or comparison.display_resolved_name
                outdir = Path(self._safe_output_name(display_group)) / self._safe_output_name(
                    f"{comparison.display_resolved_name}_{gene_symbol or gene_id or event_type}_{event_type}"
                )[:120]
                rows.append(
                    {
                        "display_group": display_group,
                        "comparison_id": comparison.comparison_id,
                        "comparison_name": comparison.display_resolved_name,
                        "geneSymbol": gene_symbol,
                        "GeneID": gene_id,
                        "event_uid": event_uid,
                        "event_type": event_type,
                        "class_interpretation": str(row.get("class_interpretation", row.get("class_label", row.get("class", "")))),
                        "label1": label1,
                        "label2": label2,
                        "b1_txt": str(b1_txt),
                        "b2_txt": str(b2_txt),
                        "event_file": str(event_file),
                        "source_rmats_file": str(
                            self._resolve_rmats_event_table_path(comparison, event_type, mode) or ""
                        ),
                        "outdir": str(outdir),
                    }
                )
                seen.add(key)

        return pd.DataFrame(rows)

    def _sashimi_target_comparisons(self, row: pd.Series) -> list[str]:
        class_name = str(row.get("class", "")).strip()
        candidates: list[str] = []
        comparison_a = str(row.get("comparison_A", "")).strip()
        comparison_b = str(row.get("comparison_B", "")).strip()
        is_sig_a = bool(row.get("is_sig_A")) if pd.notna(row.get("is_sig_A")) else False
        is_sig_b = bool(row.get("is_sig_B")) if pd.notna(row.get("is_sig_B")) else False

        if class_name == "A_only" and comparison_a:
            candidates.append(comparison_a)
        elif class_name == "B_only" and comparison_b:
            candidates.append(comparison_b)
        elif class_name in {"opposite_direction", "shared_same_direction", "same_direction_large_delta", "shared_zero_direction", "shared_unresolved"}:
            if comparison_a:
                candidates.append(comparison_a)
            if comparison_b:
                candidates.append(comparison_b)
        else:
            if is_sig_a and comparison_a:
                candidates.append(comparison_a)
            if is_sig_b and comparison_b:
                candidates.append(comparison_b)
            if not candidates:
                if comparison_a:
                    candidates.append(comparison_a)
                if comparison_b:
                    candidates.append(comparison_b)

        deduped: list[str] = []
        seen: set[str] = set()
        for item in candidates:
            if item and item not in seen:
                deduped.append(item)
                seen.add(item)
        return deduped

    def _write_generated_sashimi_event_file(
        self,
        comparison,
        event_uid: str,
        event_type: str,
        mode: str,
        event_root: Path,
        table_cache: dict[tuple[str, str, str], pd.DataFrame],
    ) -> Path | None:
        source_path = self._resolve_rmats_event_table_path(comparison, event_type, mode)
        if source_path is None:
            return None
        cache_key = (comparison.comparison_id, event_type, mode)
        if cache_key not in table_cache:
            frame = pd.read_csv(source_path, sep="\t", low_memory=False)
            if not frame.empty:
                frame = frame.copy()
                frame["event_uid"] = frame.apply(lambda item: event_uid_from_row(item, event_type), axis=1)
            table_cache[cache_key] = frame
        frame = table_cache[cache_key]
        if frame.empty or "event_uid" not in frame.columns:
            return None
        match = frame.loc[frame["event_uid"].astype(str) == event_uid].copy()
        if match.empty:
            return None
        match = match.drop(columns=["event_uid"], errors="ignore").head(1)
        target_dir = event_root / self._safe_output_name(comparison.comparison_id)
        target_dir.mkdir(parents=True, exist_ok=True)
        event_slug = self._safe_output_name(f"{event_type}_{event_uid}")[:140]
        target_path = target_dir / f"{event_slug}.{mode}.tsv"
        match.to_csv(target_path, sep="\t", index=False)
        return target_path

    def _resolve_rmats_event_table_path(self, comparison, event_type: str, mode: str) -> Path | None:
        if comparison.rmats_path is None:
            return None
        candidates = [
            comparison.rmats_path / "rmats_post" / f"{event_type}.MATS.{mode}.txt",
            comparison.rmats_path / "rmats_post" / f"{event_type}.MATS.JC.txt",
            comparison.rmats_path / "rmats_post" / f"{event_type}.MATS.JCEC.txt",
        ]
        for candidate in candidates:
            if candidate.exists():
                return candidate
        return None

    def _prepare_generated_sashimi_bam_lists(self, comparison, bam_root: Path) -> tuple[Path, Path, str, str] | None:
        analysis_exp = (comparison.experiment_group or comparison.source_experiment_group or "").strip()
        analysis_ctrl = (comparison.control_group or comparison.source_control_group or "").strip()
        source_exp, source_ctrl = comparison.source_groups("rmats")
        source_exp = (source_exp or "").strip()
        source_ctrl = (source_ctrl or "").strip()

        label1 = analysis_exp or source_exp or "group_1"
        label2 = analysis_ctrl or source_ctrl or "group_2"
        b1_paths = self._resolve_sashimi_bam_group_paths(label1)
        b2_paths = self._resolve_sashimi_bam_group_paths(label2)

        if not b1_paths and source_exp and source_exp != label1:
            b1_paths = self._resolve_sashimi_bam_group_paths(source_exp)
            if b1_paths:
                label1 = source_exp
        if not b2_paths and source_ctrl and source_ctrl != label2:
            b2_paths = self._resolve_sashimi_bam_group_paths(source_ctrl)
            if b2_paths:
                label2 = source_ctrl
        if not b1_paths or not b2_paths:
            return None

        safe_comp = self._safe_output_name(comparison.comparison_id)
        b1_txt = bam_root / f"{safe_comp}_{self._safe_output_name(label1)}_b1.txt"
        b2_txt = bam_root / f"{safe_comp}_{self._safe_output_name(label2)}_b2.txt"
        b1_txt.write_text("\n".join(str(path) for path in b1_paths) + "\n", encoding="utf-8")
        b2_txt.write_text("\n".join(str(path) for path in b2_paths) + "\n", encoding="utf-8")
        return b1_txt, b2_txt, label1, label2

    def _resolve_sashimi_bam_group_paths(self, group_name: str) -> list[Path]:
        if self.current_project is None or not group_name:
            return []
        bam_list = self._find_sashimi_bam_list_file(group_name)
        if bam_list is None or not bam_list.exists():
            return []

        entries: list[str] = []
        for raw_line in bam_list.read_text(encoding="utf-8").splitlines():
            for token in raw_line.split(","):
                token = token.strip()
                if token:
                    entries.append(token)

        roots: list[Path] = []
        configured = self.current_project.tool_paths.get("bam_root")
        if configured:
            configured_path = Path(configured)
            if configured_path.exists():
                roots.append(configured_path)
        if self.current_project.rmats_root is not None and self.current_project.rmats_root.parent.exists():
            roots.append(self.current_project.rmats_root.parent)
        if len(bam_list.parents) >= 3 and bam_list.parents[2].exists():
            roots.append(bam_list.parents[2])
        roots.append(self.current_project.project_root)

        resolved: list[Path] = []
        seen: set[Path] = set()
        for entry in entries:
            entry_path = Path(entry)
            candidate: Path | None = None
            if entry_path.is_absolute() and entry_path.exists():
                candidate = entry_path
            else:
                for root in roots:
                    direct = root / entry
                    if direct.exists():
                        candidate = direct
                        break
                if candidate is None:
                    for root in roots:
                        matches = list(root.rglob(entry_path.name))
                        if matches:
                            candidate = matches[0]
                            break
            if candidate is not None and candidate not in seen:
                resolved.append(candidate)
                seen.add(candidate)
        return resolved

    def _find_sashimi_bam_list_file(self, group_name: str) -> Path | None:
        if self.current_project is None or self.current_project.rmats_root is None:
            return None
        exact = self.current_project.rmats_root / "bamlist" / f"{group_name}_bamlist.txt"
        if exact.exists():
            return exact
        matches = list(self.current_project.rmats_root.rglob(f"{group_name}*bamlist*.txt"))
        if matches:
            return matches[0]
        return None

    def _write_analysis_outputs(self) -> None:
        if self.current_project is None or self.current_project.output_root is None:
            return
        self.current_project.output_root.mkdir(parents=True, exist_ok=True)
        comparison_config_path = self.current_project.output_root / "comparison_config.yaml"
        with comparison_config_path.open("w", encoding="utf-8") as handle:
            yaml.safe_dump(self._comparison_config_payload(), handle, sort_keys=False, allow_unicode=True)

        self._write_step_output(
            "01_splicing_landscape",
            {"splicing_landscape.tsv": self.run_state.results.splicing_landscape},
        )
        self._write_step_output(
            "02_program_comparison",
            {
                "program_events.tsv": self.run_state.results.program_events,
                "program_summary.tsv": self.run_state.results.program_summary,
            },
        )
        self._write_step_output(
            "03_mechanism_support",
            {"mechanism_support.tsv": self.run_state.results.mechanism_support},
        )
        self._write_step_output(
            "04_tx_splicing_integration",
            {
                "tx_splicing_gene_table.tsv": self.run_state.results.tx_splicing_gene_table,
                "tx_splicing_summary.tsv": self.run_state.results.tx_splicing_summary,
            },
        )
        self._write_step_output(
            "04a_candidate_gene_screening",
            self._candidate_output_tables(),
        )
        self._write_per_comparison_candidate_outputs()
        self._write_cross_comparison_outputs()
        self._write_step_output(
            "04b_shortlist",
            self._shortlist_output_tables(),
        )
        self._write_step_output(
            "05_cards",
            {
                "cards_shortlist.tsv": self.run_state.results.cards_shortlist,
                "cards_expression_support.tsv": self.run_state.results.cards_expression_support,
            },
        )
        self._write_step_output(
            "06_sashimi",
            {
                "rmats2sashimi_manifest.tsv": self.run_state.results.sashimi_manifest,
            },
        )
        self._write_step_output(
            "07_isoform_switch",
            {
                "isoform_quant_manifest.tsv": self.run_state.results.isoform_manifest,
            },
        )
        self._write_step_output(
            "08_jutils",
            {
                "jutils_manifest.tsv": self.run_state.results.jutils_manifest,
            },
        )

        self._save_splicing_landscape_plot()
        self._save_program_summary_plot()
        self._save_tx_splicing_plot()
        self._save_shortlist_plot()
        self._save_cards_summary_plot()
        self._save_sashimi_summary_plot()
        self._save_isoform_summary_plot()
        self._save_jutils_summary_plot()

    def _write_analysis_outputs_for_modules(self, selected_modules: set[str]) -> None:
        if self.current_project is None or self.current_project.output_root is None:
            return
        self.current_project.output_root.mkdir(parents=True, exist_ok=True)
        comparison_config_path = self.current_project.output_root / "comparison_config.yaml"
        with comparison_config_path.open("w", encoding="utf-8") as handle:
            yaml.safe_dump(self._comparison_config_payload(), handle, sort_keys=False, allow_unicode=True)

        if "splicing_landscape" in selected_modules:
            self._write_step_output(
                "01_splicing_landscape",
                {"splicing_landscape.tsv": self.run_state.results.splicing_landscape},
            )
            self._save_splicing_landscape_plot()

        if "program_comparison" in selected_modules:
            self._write_step_output(
                "02_program_comparison",
                {
                    "program_events.tsv": self.run_state.results.program_events,
                    "program_summary.tsv": self.run_state.results.program_summary,
                },
            )
            self._save_program_summary_plot()

        if "mechanism_support" in selected_modules:
            self._write_step_output(
                "03_mechanism_support",
                {"mechanism_support.tsv": self.run_state.results.mechanism_support},
            )

        if "tx_splicing_integration" in selected_modules:
            self._write_step_output(
                "04_tx_splicing_integration",
                {
                    "tx_splicing_gene_table.tsv": self.run_state.results.tx_splicing_gene_table,
                    "tx_splicing_summary.tsv": self.run_state.results.tx_splicing_summary,
                },
            )
            self._save_tx_splicing_plot()

    def _write_step_output(self, step_name: str, tables: dict[str, pd.DataFrame]) -> None:
        if self.current_project is None or self.current_project.output_root is None:
            return
        step_dir = self.current_project.output_root / step_name
        data_dir = step_dir / "data"
        figures_dir = step_dir / "figures"
        data_dir.mkdir(parents=True, exist_ok=True)
        figures_dir.mkdir(parents=True, exist_ok=True)
        for filename, frame in tables.items():
            if frame is not None and not frame.empty:
                frame.to_csv(data_dir / filename, sep="\t", index=False)

    def _save_splicing_landscape_plot(self) -> None:
        if self.current_project is None or self.current_project.output_root is None:
            return
        frame = self.run_state.results.splicing_landscape
        if frame.empty:
            return
        figures_dir = self.current_project.output_root / "01_splicing_landscape" / "figures"
        figures_dir.mkdir(parents=True, exist_ok=True)
        pivot = frame.pivot_table(
            index="comparison_name",
            columns="event_type",
            values="n_sig",
            aggfunc="first",
        ).fillna(0)
        ax = pivot.plot(kind="bar", stacked=True, figsize=(10, 6))
        ax.set_ylabel("Significant events")
        ax.set_title("Splicing landscape by comparison")
        plt.tight_layout()
        plt.savefig(figures_dir / "splicing_landscape_stacked.png", dpi=300)
        plt.close()

    def _save_program_summary_plot(self) -> None:
        if self.current_project is None or self.current_project.output_root is None:
            return
        frame = self.run_state.results.program_summary
        if frame.empty:
            return
        figures_dir = self.current_project.output_root / "02_program_comparison" / "figures"
        figures_dir.mkdir(parents=True, exist_ok=True)
        ax = frame.plot(kind="bar", x="class_label", y="n_events", figsize=(9, 5), legend=False)
        ax.set_ylabel("Events")
        ax.set_title("Program comparison summary")
        plt.xticks(rotation=25, ha="right")
        plt.tight_layout()
        plt.savefig(figures_dir / "program_summary.png", dpi=300)
        plt.close()

    def _save_tx_splicing_plot(self) -> None:
        if self.current_project is None or self.current_project.output_root is None:
            return
        frame = self.run_state.results.tx_splicing_gene_table
        if frame.empty:
            return
        figures_dir = self.current_project.output_root / "04_tx_splicing_integration" / "figures"
        figures_dir.mkdir(parents=True, exist_ok=True)
        plot_frame = frame.copy()
        plot_frame["abs_log2FC"] = pd.to_numeric(
            plot_frame.get("standardized_log2FC", plot_frame.get("log2FC")), errors="coerce"
        ).abs()
        plot_frame["abs_dPSI"] = pd.to_numeric(
            plot_frame.get("standardized_dPSI", plot_frame.get("representative_dPSI")), errors="coerce"
        ).abs()
        ax = plot_frame.plot.scatter(
            x="abs_log2FC",
            y="abs_dPSI",
            figsize=(7, 6),
            alpha=0.6,
        )
        ax.set_xlabel("|standardized log2FC|")
        ax.set_ylabel("|standardized dPSI|")
        ax.set_title("Transcriptome x splicing integration")
        plt.tight_layout()
        plt.savefig(figures_dir / "tx_splicing_scatter.png", dpi=300)
        plt.close()

    def _save_cards_summary_plot(self) -> None:
        if self.current_project is None or self.current_project.output_root is None:
            return
        frame = self.run_state.results.cards_shortlist
        if frame.empty or "display_group" not in frame.columns:
            return
        figures_dir = self.current_project.output_root / "05_cards" / "figures"
        figures_dir.mkdir(parents=True, exist_ok=True)
        counts = frame["display_group"].value_counts().reset_index()
        counts.columns = ["display_group", "n_rows"]
        ax = counts.plot(kind="bar", x="display_group", y="n_rows", figsize=(8, 5), legend=False)
        ax.set_ylabel("Rows")
        ax.set_title("Cards shortlist groups")
        plt.xticks(rotation=20, ha="right")
        plt.tight_layout()
        plt.savefig(figures_dir / "cards_groups.png", dpi=300)
        plt.close()

    def _save_shortlist_plot(self) -> None:
        if self.current_project is None or self.current_project.output_root is None:
            return
        frame = self.run_state.results.cards_shortlist
        if frame.empty or "display_group" not in frame.columns:
            return
        figures_dir = self.current_project.output_root / "04b_shortlist" / "figures"
        figures_dir.mkdir(parents=True, exist_ok=True)
        counts = frame["display_group"].value_counts().reset_index()
        counts.columns = ["display_group", "n_genes"]
        ax = counts.plot(kind="bar", x="display_group", y="n_genes", figsize=(8, 5), legend=False)
        ax.set_ylabel("Genes")
        ax.set_title("Shortlist groups")
        plt.xticks(rotation=20, ha="right")
        plt.tight_layout()
        plt.savefig(figures_dir / "shortlist_groups.png", dpi=300)
        plt.close()

    def _save_sashimi_summary_plot(self) -> None:
        if self.current_project is None or self.current_project.output_root is None:
            return
        frame = self.run_state.results.sashimi_manifest
        if frame.empty or "comparison_id" not in frame.columns:
            return
        figures_dir = self.current_project.output_root / "06_sashimi" / "figures"
        figures_dir.mkdir(parents=True, exist_ok=True)
        counts = frame["comparison_id"].value_counts().reset_index()
        counts.columns = ["comparison_id", "n_rows"]
        ax = counts.plot(kind="bar", x="comparison_id", y="n_rows", figsize=(8, 5), legend=False)
        ax.set_ylabel("Events")
        ax.set_title("Sashimi manifest by comparison")
        plt.xticks(rotation=20, ha="right")
        plt.tight_layout()
        plt.savefig(figures_dir / "sashimi_manifest_counts.png", dpi=300)
        plt.close()

    def _save_isoform_summary_plot(self) -> None:
        if self.current_project is None or self.current_project.output_root is None:
            return
        frame = self.run_state.results.isoform_manifest
        if frame.empty:
            return
        figures_dir = self.current_project.output_root / "07_isoform_switch" / "figures"
        figures_dir.mkdir(parents=True, exist_ok=True)
        counts = frame.copy()
        counts["source_dir"] = counts["quant_sf"].map(lambda value: str(Path(value).parent))
        counts = counts["source_dir"].value_counts().reset_index()
        counts.columns = ["source_dir", "n_files"]
        ax = counts.plot(kind="bar", x="source_dir", y="n_files", figsize=(10, 5), legend=False)
        ax.set_ylabel("quant.sf files")
        ax.set_title("Isoform manifest sources")
        plt.xticks(rotation=30, ha="right")
        plt.tight_layout()
        plt.savefig(figures_dir / "isoform_sources.png", dpi=300)
        plt.close()

    def _save_jutils_summary_plot(self) -> None:
        if self.current_project is None or self.current_project.output_root is None:
            return
        frame = self.run_state.results.jutils_manifest
        if frame.empty or "comparison_name" not in frame.columns:
            return
        figures_dir = self.current_project.output_root / "08_jutils" / "figures"
        figures_dir.mkdir(parents=True, exist_ok=True)
        counts = frame["comparison_name"].value_counts().reset_index()
        counts.columns = ["comparison_name", "n_rows"]
        ax = counts.plot(kind="bar", x="comparison_name", y="n_rows", figsize=(8, 5), legend=False)
        ax.set_ylabel("Entries")
        ax.set_title("Jutils manifest preview")
        plt.xticks(rotation=20, ha="right")
        plt.tight_layout()
        plt.savefig(figures_dir / "jutils_manifest_counts.png", dpi=300)
        plt.close()

    def _load_default_thresholds(self) -> dict:
        config_path = bundled_configs_root() / "default_thresholds.yaml"
        with config_path.open("r", encoding="utf-8") as handle:
            return yaml.safe_load(handle)
