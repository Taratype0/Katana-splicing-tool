from __future__ import annotations

import json
import hashlib
from pathlib import Path
from datetime import datetime
from dataclasses import asdict
import gzip
import math
import re
import shutil
import urllib.error
import urllib.request
from itertools import combinations

import yaml
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import importlib.util
import numpy as np

from src.adapters.jutils_adapter import JutilsAdapter
from src.adapters.isoformswitch_runner import IsoformSwitchRunner
from src.adapters.cards_adapter import CardsAdapter
from src.adapters.sashimi_adapter import SashimiAdapter
from src.adapters.native_sashimi_adapter import NativeSashimiAdapter, _parse_gtf_attrs
from src.analysis.common import EVENT_TYPES, event_uid_from_row, load_annotation_gene_symbol_lookup, normalize_gene_id
from src.analysis.mechanism_support import MechanismSupportAnalyzer
from src.analysis.native_splicing_patterns import (
    psi_pattern_clusters,
    replicate_group_separation,
    row_zscore_matrix,
    sample_similarity_tables,
)
from src.analysis.program_comparison import SplicingProgramComparator
from src.analysis.splicing_landscape import SplicingLandscapeAnalyzer
from src.analysis.stable_as import apply_stable_as_filter, pairwise_direction_consistency
from src.analysis.tx_splicing_integration import TxSplicingIntegrator
from src.runtime_paths import bundled_configs_root, bundled_software_root
from src.services.artifact_registry import AnalysisArtifactRegistry
from src.services.cross_comparison_service import CrossComparisonService
from src.services.jutils_benchmark_service import JutilsBenchmarkService
from src.services.output_directory_service import OutputDirectoryService, OutputMigrationResult
from src.services.project_report_service import ProjectReportService
from src.services.shortlist_service import ShortlistService
from src.adapters.isoformswitch_adapter import IsoformSwitchAdapter
from src.domain.models.analysis import AnalysisRequest, AnalysisResults, AnalysisRunState, ThresholdConfig
from src.domain.models.project import ComparisonPairDefinition, IsoformSampleDefinition, ProjectConfig, VisualizationGroupDefinition
from src.io.scanner import ProjectScanner


class ProjectService:
    NATIVE_PATTERN_SCHEMA_VERSION = "3"
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
        self.native_sashimi_adapter = NativeSashimiAdapter()
        self.isoform_runner: IsoformSwitchRunner | None = None
        self.shortlist_service = ShortlistService()
        self.output_directories = OutputDirectoryService(self.application_source_root())
        self.jutils_benchmark = JutilsBenchmarkService()
        self.default_thresholds = self._load_default_thresholds()
        self.run_state = AnalysisRunState()
        self.current_thresholds = self._threshold_config_from_mapping(self.default_thresholds["presets"]["default"])
        self._last_jutils_report: dict[str, object] = {
            "status": "idle",
            "command_lines": [],
            "output_directory": "",
            "generated_files": [],
            "stdout": "",
            "stderr": "",
            "error_log": "",
        }
        self._last_native_jutils_report: dict[str, object] = {
            "status": "idle",
            "mode": "pattern",
            "output_directory": "",
            "generated_files": [],
            "stdout": "",
            "stderr": "",
            "error_log": "",
        }
        self._last_sashimi_failures = pd.DataFrame()
        self._sashimi_event_catalog_cache: dict[str, pd.DataFrame] = {}
        self._cross_as_pattern_cache: dict[tuple[str, float, float], dict[str, pd.DataFrame]] = {}
        self._cross_as_interaction_cache: dict[tuple[str, float, float, int], dict[str, pd.DataFrame]] = {}
        self._cross_as_interaction_batch_cache: dict[tuple[tuple[str, ...], float, float, int], dict[str, pd.DataFrame]] = {}
        self._candidate_catalog_cache: dict[tuple[str, int, tuple[str, ...]], pd.DataFrame] = {}
        self._significant_as_event_summary_cache: dict[str, pd.DataFrame] = {}
        self._rbp_support_cache: dict[tuple[str, str, int, str], dict[str, pd.DataFrame]] = {}
        self._rbp_cross_support_cache: dict[tuple[str, int, str, tuple[str, ...]], dict[str, pd.DataFrame]] = {}
        self._isoform_support_cache: dict[tuple[str, str, int], dict[str, pd.DataFrame]] = {}
        self._isoform_cross_support_cache: dict[tuple[str, int, tuple[str, ...]], dict[str, pd.DataFrame]] = {}
        self._native_jutils_cache: dict[tuple[str, tuple[str, ...], float, float, str], dict[str, pd.DataFrame]] = {}
        self._gtf_transcript_gene_cache: dict[str, pd.DataFrame] = {}
        self._deg_expression_cache: dict[str, pd.DataFrame] = {}
        self._comparison_expression_support_cache: dict[str, pd.DataFrame] = {}
        self._expression_support_cache = pd.DataFrame()
        self._expression_support_loaded = False
        self._candidate_cache_fingerprint_in_memory = ""
        self._artifact_registry: AnalysisArtifactRegistry | None = None
        self._module_states = self._new_module_states()

    def _clear_candidate_selection_caches(self) -> None:
        self._candidate_catalog_cache = {}
        self._rbp_support_cache = {}
        self._rbp_cross_support_cache = {}
        self._isoform_support_cache = {}
        self._isoform_cross_support_cache = {}
        self._native_jutils_cache = {}

    def _emit_progress(self, callback, message: str) -> None:
        if callback is not None:
            callback(message)

    def _log_heavy(self, message: str) -> None:
        print(f"[katana-heavy] {message}")

    def _log_cache(self, message: str) -> None:
        print(f"[katana-cache] {message}")

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
            "cross_comparison_as_interactions": {"status": "not_run", "message": "Not run."},
            "sashimi_manifest": {"status": "not_run", "message": "Not run."},
            "sashimi_plot": {"status": "not_run", "message": "Not run."},
            "jutils_heatmap": {"status": "not_run", "message": "Not run."},
            "jutils_pca": {"status": "not_run", "message": "Not run."},
            "rbp_support": {"status": "not_run", "message": "Not run."},
            "rbp_cross_support": {"status": "not_run", "message": "Not run."},
            "isoform": {"status": "not_run", "message": "Not run."},
            "isoform_support": {"status": "not_run", "message": "Not run."},
            "isoform_cross_support": {"status": "not_run", "message": "Not run."},
            "deg_explorer": {"status": "not_run", "message": "Not run."},
            "functional_enrichment": {"status": "not_run", "message": "Not run."},
        }

    def _reset_module_states(self) -> None:
        self._module_states = self._new_module_states()

    @staticmethod
    def _safe_bool(value) -> bool:
        if value is None:
            return False
        try:
            if pd.isna(value):
                return False
        except (TypeError, ValueError):
            pass
        if isinstance(value, str):
            return value.strip().lower() in {"1", "true", "yes", "y", "t"}
        return bool(value)

    @staticmethod
    def _safe_text(value, default: str = "") -> str:
        if value is None:
            return default
        try:
            if pd.isna(value):
                return default
        except (TypeError, ValueError):
            pass
        text = str(value).strip()
        return text if text else default

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
        self._cross_as_interaction_cache = {}
        self._candidate_catalog_cache = {}
        self._significant_as_event_summary_cache = {}
        self._rbp_support_cache = {}
        self._rbp_cross_support_cache = {}
        self._isoform_support_cache = {}
        self._isoform_cross_support_cache = {}
        self._native_jutils_cache = {}
        self._gtf_transcript_gene_cache = {}
        self._deg_expression_cache = {}
        self._comparison_expression_support_cache = {}
        self._expression_support_cache = pd.DataFrame()
        self._expression_support_loaded = False
        self._candidate_cache_fingerprint_in_memory = ""
        self._artifact_registry = None
        self.current_project = self.scanner.scan(
            root_dir,
            input_paths=input_paths,
            force_auto_scan=force_auto_scan,
            progress_callback=progress_callback,
        )
        self._emit_progress(progress_callback, "load comparison config")
        self._emit_progress(progress_callback, "load pairing direction")
        self._refresh_tool_adapters()
        self._normalize_project_selection_state()
        if not self.current_project.selected_analysis_modules:
            self.current_project.selected_analysis_modules = list(self.ANALYSIS_MODULES)
        self._ensure_default_comparison_pairs()
        self._ensure_default_visualization_groups()
        saved_thresholds = self.current_project.tool_paths.get("saved_thresholds")
        if isinstance(saved_thresholds, dict) and saved_thresholds:
            self.current_thresholds = self._threshold_config_from_mapping(saved_thresholds)
        self._emit_progress(progress_callback, "load cached candidate selection state")
        self._autodetect_embedded_tools()
        self._emit_progress(progress_callback, "finish load project")
        return self.current_project

    def artifact_registry(self) -> AnalysisArtifactRegistry:
        resolved_root = self.require_output_root()
        if self._artifact_registry is None or self._artifact_registry.output_root != resolved_root:
            self._artifact_registry = AnalysisArtifactRegistry(
                resolved_root,
                self.current_project.project_name,
            )
        return self._artifact_registry

    @staticmethod
    def application_source_root() -> Path:
        return Path(__file__).resolve().parents[2]

    def output_root_validation_error(self, output_root: str | Path | None = None) -> str:
        if self.current_project is None:
            return "No project loaded."
        candidate = output_root if output_root is not None else self.current_project.output_root
        return self.output_directories.validation_error(candidate)

    def require_output_root(self, *, create: bool = True) -> Path:
        if self.current_project is None:
            raise RuntimeError("No project loaded.")
        return self.output_directories.require(self.current_project.output_root, create=create)

    def default_export_path(self, filename: str) -> Path:
        export_root = self.require_output_root() / "exports"
        export_root.mkdir(parents=True, exist_ok=True)
        return export_root / filename

    def build_project_report(self) -> dict[str, Path]:
        if self.current_project is None:
            raise RuntimeError("No project loaded.")
        project = self.current_project
        output_root = self.require_output_root()
        report_service = ProjectReportService(output_root, project.project_name)
        selected_comparisons = self.selected_comparisons_for_display()
        outputs = report_service.build(
            registry=self.artifact_registry(),
            project_summary={
                "project_name": project.project_name,
                "project_root": str(project.project_root),
                "output_root": str(output_root),
                "selected_comparisons": len(selected_comparisons),
                "comparison_ids": ", ".join(item.comparison_id for item in selected_comparisons),
                "pairing_confirmed": bool(project.pairing_confirmed),
                "comparison_sets_confirmed": bool(project.comparison_sets_confirmed),
                "visualization_groups_confirmed": bool(project.visualization_groups_confirmed),
            },
            module_states={key: dict(value) for key, value in self._module_states.items()},
        )
        source_ids = [
            record.artifact_id
            for record in self.artifact_registry().find()
            if record.module_id != "project_report"
        ]
        for artifact_key, path in outputs.items():
            self.artifact_registry().register(
                module_id="project_report",
                artifact_key=artifact_key,
                kind="report" if path.suffix.lower() == ".html" else "table" if path.suffix.lower() == ".tsv" else "parameters",
                entity_type="project",
                path=path,
                source_artifact_ids=source_ids,
                consumer_modules=["project_report"],
                description="Cached project-wide provenance and artifact inventory report.",
            )
        return outputs

    def _normalize_project_selection_state(self) -> None:
        if self.current_project is None:
            return
        available_ids = [item.comparison_id for item in self.current_project.available_comparisons]
        available_lookup = set(available_ids)

        selected = [
            item
            for item in dict.fromkeys(self.current_project.selected_comparison_ids)
            if item in available_lookup
        ]
        if not selected:
            selected = available_ids[:4]
        self.current_project.selected_comparison_ids = selected

        ordered = [
            item
            for item in dict.fromkeys(self.current_project.comparison_order_ids)
            if item in available_lookup
        ]
        for comparison_id in available_ids:
            if comparison_id not in ordered:
                ordered.append(comparison_id)
        self.current_project.comparison_order_ids = ordered

        program_ids = [
            item
            for item in dict.fromkeys(self.current_project.program_comparison_ids)
            if item in set(selected)
        ]
        self.current_project.program_comparison_ids = program_ids[:2]

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
        self._cross_as_interaction_cache = {}
        self._cross_as_interaction_batch_cache = {}
        self._candidate_catalog_cache = {}
        self._significant_as_event_summary_cache = {}
        self._candidate_cache_fingerprint_in_memory = ""
        self._deg_expression_cache = {}
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
        self._cross_as_interaction_cache = {}
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
        return self._dedupe_gene_tokens(self.current_project.candidate_selection_genes.get(comparison_id, []))

    def candidate_blacklist_genes_for_display(self, comparison_id: str) -> list[str]:
        if self.current_project is None or not comparison_id:
            return []
        return self._dedupe_gene_tokens(self.current_project.candidate_blacklist_genes.get(comparison_id, []))

    @staticmethod
    def _dedupe_gene_tokens(values: list[str] | tuple[str, ...] | set[str] | None) -> list[str]:
        ordered: list[str] = []
        seen: set[str] = set()
        for value in values or []:
            text = str(value or "").strip()
            if not text:
                continue
            key = text.lower()
            if key in seen:
                continue
            seen.add(key)
            ordered.append(text)
        return ordered

    def publish_analysis_selection(
        self,
        *,
        module_id: str,
        artifact_key: str,
        entity_type: str,
        frame: pd.DataFrame,
        scope: str,
        scope_id: str,
        comparison_ids: list[str] | tuple[str, ...] = (),
        id_columns: list[str] | tuple[str, ...] = (),
        source_artifact_ids: list[str] | tuple[str, ...] = (),
        consumer_modules: list[str] | tuple[str, ...] = (),
        selection_rule: str = "",
        parameters: dict[str, object] | None = None,
    ):
        return self.artifact_registry().publish_selection(
            module_id=module_id,
            artifact_key=artifact_key,
            entity_type=entity_type,
            frame=frame,
            scope=scope,
            scope_id=scope_id,
            comparison_ids=comparison_ids,
            id_columns=id_columns,
            source_artifact_ids=source_artifact_ids,
            consumer_modules=consumer_modules,
            selection_rule=selection_rule,
            parameters=parameters,
        )

    def publish_candidate_selection_artifacts(self, comparison_id: str) -> None:
        if self.current_project is None or not comparison_id:
            return
        selected_genes = self.candidate_selection_genes_for_display(comparison_id)
        blacklist_genes = self.candidate_blacklist_genes_for_display(comparison_id)
        catalog = self.candidate_selection_catalog(comparison_id, allow_generate=False)
        selected = self._candidate_custom_gene_rows(catalog, comparison_id, selected_genes) if not catalog.empty else pd.DataFrame()
        if selected.empty:
            selected = pd.DataFrame({"gene_symbol": selected_genes})
        if "comparison_id" not in selected.columns:
            selected.insert(0, "comparison_id", comparison_id)
        blacklist = pd.DataFrame(
            {
                "comparison_id": [comparison_id] * len(blacklist_genes),
                "gene_symbol": blacklist_genes,
            }
        )
        self.publish_analysis_selection(
            module_id="5.3.2",
            artifact_key="candidate_gene_selection",
            entity_type="gene",
            frame=selected,
            scope="comparison",
            scope_id=comparison_id,
            comparison_ids=[comparison_id],
            id_columns=["comparison_id", "gene_symbol", "gene_id"],
            consumer_modules=["5.3.3", "5.3.4", "5.3.5", "5.5", "7.1", "8.1"],
            selection_rule="Comparison-specific candidate selection after blacklist exclusion.",
            parameters={"blacklist_count": len(blacklist_genes)},
        )
        self.publish_analysis_selection(
            module_id="5.3.2",
            artifact_key="candidate_gene_blacklist",
            entity_type="gene",
            frame=blacklist,
            scope="comparison",
            scope_id=comparison_id,
            comparison_ids=[comparison_id],
            id_columns=["comparison_id", "gene_symbol"],
            consumer_modules=["5.3.3", "5.3.4", "5.3.5", "5.5", "7.1", "8.1"],
            selection_rule="Genes excluded from downstream candidate visualizations.",
        )

    def downstream_gene_selection(
        self,
        comparison_id: str,
        *,
        consumer_module: str,
    ) -> tuple[pd.DataFrame, dict[str, object]]:
        frame, artifact = self.artifact_registry().latest_selection(
            entity_type="gene",
            scope_id=comparison_id,
            module_id="5.3.2",
            artifact_key="candidate_gene_selection",
        )
        if artifact is not None and not frame.empty:
            return frame, {
                "artifact_id": artifact.artifact_id,
                "source_module": artifact.module_id,
                "selection_rule": artifact.selection_rule,
                "fingerprint": artifact.fingerprint,
                "consumer_module": consumer_module,
                "cache_status": "artifact_registry",
            }
        genes = self.candidate_selection_genes_for_display(comparison_id)
        fallback = pd.DataFrame(
            {
                "comparison_id": [comparison_id] * len(genes),
                "gene_symbol": genes,
            }
        )
        return fallback, {
            "artifact_id": "",
            "source_module": "project_config",
            "selection_rule": "Legacy comparison-specific candidate selection.",
            "fingerprint": "",
            "consumer_module": consumer_module,
            "cache_status": "project_config_fallback",
        }

    def publish_sashimi_event_selection(
        self,
        frame: pd.DataFrame,
        *,
        comparison_id: str,
        selection_rule: str,
    ) -> None:
        if frame is None or frame.empty or not comparison_id:
            return
        selection = frame.copy()
        if "comparison_id" not in selection.columns:
            selection.insert(0, "comparison_id", comparison_id)
        if "event_id" not in selection.columns and "event_uid" in selection.columns:
            selection["event_id"] = selection["event_uid"]
        if "event_matching_key" not in selection.columns:
            gene = selection.get("gene_symbol", selection.get("geneSymbol", pd.Series("", index=selection.index)))
            event_type = selection.get("event_type", pd.Series("", index=selection.index))
            event_id = selection.get("event_id", pd.Series("", index=selection.index))
            selection["event_matching_key"] = (
                gene.fillna("").astype(str)
                + "|"
                + event_type.fillna("").astype(str)
                + "|"
                + event_id.fillna("").astype(str)
            )
        self.publish_analysis_selection(
            module_id="5.5.2",
            artifact_key="selected_sashimi_events",
            entity_type="event",
            frame=selection,
            scope="comparison",
            scope_id=comparison_id,
            comparison_ids=[comparison_id],
            id_columns=["comparison_id", "event_matching_key", "event_id"],
            consumer_modules=["5.5.2", "5.5.3", "9.CoSplice"],
            selection_rule=selection_rule,
        )

    @staticmethod
    def cosplice_module_catalog(output_dir: str | Path) -> pd.DataFrame:
        root = Path(output_dir)
        summary_path = root / "04c_consensus_modules" / "consensus_module_summary.tsv"
        assignments_path = root / "04c_consensus_modules" / "consensus_event_module_assignments.tsv"
        if not assignments_path.exists():
            assignments_path = root / "04_network" / "observed" / "event_module_assignments.tsv"
        summary = pd.DataFrame()
        if summary_path.exists():
            try:
                summary = pd.read_csv(summary_path, sep="\t", low_memory=False)
            except (OSError, pd.errors.ParserError, UnicodeDecodeError):
                summary = pd.DataFrame()
        assignments = pd.DataFrame()
        if assignments_path.exists():
            try:
                assignments = pd.read_csv(assignments_path, sep="\t", low_memory=False)
            except (OSError, pd.errors.ParserError, UnicodeDecodeError):
                assignments = pd.DataFrame()
        if summary.empty and not assignments.empty and "module" in assignments.columns:
            aggregations: dict[str, tuple[str, str]] = {
                "n_events": ("module", "size"),
            }
            for output_name, source_name in (
                ("mean_module_membership", "module_membership"),
                ("mean_network_connectivity", "network_connectivity"),
                ("mean_bootstrap_support", "mean_incident_bootstrap_support"),
            ):
                if source_name in assignments.columns:
                    aggregations[output_name] = (source_name, "mean")
            summary = assignments.groupby("module", dropna=False).agg(**aggregations).reset_index()
        module_column = next(
            (column for column in ("module", "consensus_module", "module_id") if column in summary.columns),
            None,
        )
        if module_column is None:
            return pd.DataFrame()
        if module_column != "module":
            summary = summary.rename(columns={module_column: "module"})
        summary["module"] = summary["module"].fillna("").astype(str)
        summary = summary.loc[summary["module"].ne("") & summary["module"].ne("grey")].copy()
        return summary.drop_duplicates("module", keep="first").sort_values("module").reset_index(drop=True)

    def publish_cosplice_module_selection(
        self,
        output_dir: str | Path,
        module_ids: list[str] | tuple[str, ...],
    ) -> pd.DataFrame:
        clean_ids = list(dict.fromkeys(str(item).strip() for item in module_ids if str(item).strip()))
        if not clean_ids:
            raise RuntimeError("Select one or more CoSplice modules before publishing the handoff.")
        catalog = self.cosplice_module_catalog(output_dir)
        if catalog.empty:
            raise RuntimeError("No cached CoSplice module catalog is available.")
        selected = catalog.loc[catalog["module"].isin(clean_ids)].copy()
        if selected.empty:
            raise RuntimeError("The selected CoSplice modules are not present in the cached module catalog.")
        self.publish_analysis_selection(
            module_id="9.CoSplice",
            artifact_key="selected_cosplice_modules",
            entity_type="module",
            frame=selected,
            scope="project",
            scope_id="project",
            id_columns=["module"],
            consumer_modules=["9.CoSplice.interpretation", "7.RBP", "8.Isoform", "project_report"],
            selection_rule="User-selected cached CoSplice modules; publishing does not rerun the network.",
            parameters={"cosplice_output_dir": str(Path(output_dir).resolve())},
        )
        return selected

    def add_candidate_selection_gene(self, comparison_id: str, gene: str) -> None:
        self.add_candidate_selection_genes(comparison_id, [gene])

    def add_candidate_selection_genes(self, comparison_id: str, genes: list[str] | tuple[str, ...] | set[str]) -> None:
        if self.current_project is None:
            raise RuntimeError("No project loaded.")
        if not comparison_id:
            return
        blocked = {
            item.strip().lower()
            for item in self.current_project.candidate_blacklist_genes.get(comparison_id, [])
            if item and item.strip()
        }
        current = self._dedupe_gene_tokens(self.current_project.candidate_selection_genes.get(comparison_id, []))
        current_keys = {item.lower() for item in current}
        changed = False
        for gene in genes or []:
            value = str(gene or "").strip()
            if not value:
                continue
            value_key = value.lower()
            if value_key in blocked or value_key in current_keys:
                continue
            current.append(value)
            current_keys.add(value_key)
            changed = True
        if changed:
            self.current_project.candidate_selection_genes[comparison_id] = self._dedupe_gene_tokens(current)
            self._clear_candidate_selection_caches()
            self.save_project_config()
            self.publish_candidate_selection_artifacts(comparison_id)

    def remove_candidate_selection_gene(self, comparison_id: str, gene: str) -> None:
        if self.current_project is None or not comparison_id:
            raise RuntimeError("No project loaded.")
        gene_key = str(gene or "").strip().lower()
        self.current_project.candidate_selection_genes[comparison_id] = [
            item
            for item in self.current_project.candidate_selection_genes.get(comparison_id, [])
            if str(item or "").strip().lower() != gene_key
        ]
        self._clear_candidate_selection_caches()
        self.save_project_config()
        self.publish_candidate_selection_artifacts(comparison_id)

    def add_candidate_blacklist_gene(self, comparison_id: str, gene: str) -> None:
        self.add_candidate_blacklist_genes(comparison_id, [gene])

    def add_candidate_blacklist_genes(self, comparison_id: str, genes: list[str] | tuple[str, ...] | set[str]) -> None:
        if self.current_project is None:
            raise RuntimeError("No project loaded.")
        if not comparison_id:
            return
        current = self._dedupe_gene_tokens(self.current_project.candidate_blacklist_genes.get(comparison_id, []))
        current_keys = {item.lower() for item in current}
        blocked_keys: set[str] = set()
        changed = False
        for gene in genes or []:
            value = str(gene or "").strip()
            if not value:
                continue
            value_key = value.lower()
            blocked_keys.add(value_key)
            if value_key in current_keys:
                continue
            current.append(value)
            current_keys.add(value_key)
            changed = True
        if not blocked_keys:
            return
        selected = [
            item for item in self.current_project.candidate_selection_genes.get(comparison_id, [])
            if str(item or "").strip().lower() not in blocked_keys
        ]
        selected_deduped = self._dedupe_gene_tokens(selected)
        if selected_deduped != self.current_project.candidate_selection_genes.get(comparison_id, []):
            changed = True
        if changed:
            self.current_project.candidate_blacklist_genes[comparison_id] = self._dedupe_gene_tokens(current)
            self.current_project.candidate_selection_genes[comparison_id] = selected_deduped
            self._clear_candidate_selection_caches()
            self.save_project_config()
            self.publish_candidate_selection_artifacts(comparison_id)

    def remove_candidate_blacklist_gene(self, comparison_id: str, gene: str) -> None:
        if self.current_project is None or not comparison_id:
            raise RuntimeError("No project loaded.")
        gene_key = str(gene or "").strip().lower()
        self.current_project.candidate_blacklist_genes[comparison_id] = [
            item
            for item in self.current_project.candidate_blacklist_genes.get(comparison_id, [])
            if str(item or "").strip().lower() != gene_key
        ]
        self._clear_candidate_selection_caches()
        self.save_project_config()
        self.publish_candidate_selection_artifacts(comparison_id)

    def reset_candidate_selection(self, comparison_id: str, *, top_n: int = 20) -> None:
        if self.current_project is None or not comparison_id:
            raise RuntimeError("No project loaded.")
        self.current_project.candidate_selection_genes[comparison_id] = self._dedupe_gene_tokens(self._default_candidate_selection_genes(
            comparison_id,
            top_n=top_n,
        ))
        self._clear_candidate_selection_caches()
        self.save_project_config()
        self.publish_candidate_selection_artifacts(comparison_id)

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
        before = asdict(self.current_thresholds)
        for key, value in values.items():
            if hasattr(self.current_thresholds, key):
                setattr(self.current_thresholds, key, value)
        if asdict(self.current_thresholds) != before:
            self.invalidate_analysis_results(
                "Analysis caches are stale because thresholds changed. Re-run only the module whose output you need."
            )
            self._candidate_cache_fingerprint_in_memory = ""
            self._clear_candidate_selection_caches()
        self.save_project_config()

    def apply_threshold_preset(self, preset_name: str) -> None:
        preset = self.default_thresholds["presets"].get(preset_name)
        if preset is None:
            raise KeyError(f"Unknown preset: {preset_name}")
        updated = self._threshold_config_from_mapping(preset)
        if asdict(updated) != asdict(self.current_thresholds):
            self.current_thresholds = updated
            self.invalidate_analysis_results(
                f"Analysis caches are stale because the '{preset_name}' threshold preset was selected."
            )
            self._candidate_cache_fingerprint_in_memory = ""
            self._clear_candidate_selection_caches()
        else:
            self.current_thresholds = updated
        self.save_project_config()

    def preview_splicing_landscape(self):
        if self.current_project is None:
            raise RuntimeError("No project loaded.")
        if self.run_state.status == "completed" and not self.run_state.results.splicing_landscape.empty:
            return self.run_state.results.splicing_landscape
        request = self.build_default_request()
        request.rmats_mode = self.current_project.selected_rmats_mode
        return self.splicing_landscape.run(self.current_project, request)

    def deg_explorer_file_map(self, comparison_id: str) -> dict[str, Path]:
        if self.current_project is None or self.current_project.output_root is None:
            return {}
        comparison = next(
            (item for item in self.current_project.available_comparisons if item.comparison_id == comparison_id),
            None,
        )
        if comparison is None:
            return {}
        base = self.current_project.output_root / "05_deg_explorer" / (
            comparison.output_prefix or self._safe_output_name(comparison.comparison_id)
        )
        return {
            "all": base / "deg_all.tsv",
            "significant": base / "deg_significant.tsv",
            "up": base / "deg_up.tsv",
            "down": base / "deg_down.tsv",
            "summary": base / "summary.tsv",
            "parameters": base / "parameters.json",
            "volcano": base / "volcano.png",
            "ma": base / "ma_plot.png",
        }

    def preview_deg_explorer(
        self,
        comparison_id: str,
        *,
        allow_generate: bool = False,
        force_rebuild: bool = False,
    ) -> pd.DataFrame:
        if self.current_project is None:
            raise RuntimeError("No project loaded.")
        files = self.deg_explorer_file_map(comparison_id)
        parameters = files.get("parameters")
        expected = self._deg_explorer_signature(comparison_id)
        if not force_rebuild and parameters is not None and parameters.exists() and files.get("all", Path()).exists():
            try:
                metadata = json.loads(parameters.read_text(encoding="utf-8"))
                if metadata.get("signature") == expected:
                    frame = pd.read_csv(files["all"], sep="\t", low_memory=False)
                    self._mark_module_state(
                        "deg_explorer",
                        "cache_hit",
                        f"Loaded cached DEG Explorer result for {comparison_id}.",
                        from_cache=True,
                        output_folder=str(files["all"].parent),
                    )
                    return frame
            except (OSError, ValueError, json.JSONDecodeError):
                pass
        if not allow_generate:
            self._mark_module_state("deg_explorer", "not_run", "No cached DEG Explorer result. Click Run DEG Explorer.")
            return pd.DataFrame()
        comparison = next(
            (item for item in self.current_project.available_comparisons if item.comparison_id == comparison_id),
            None,
        )
        if comparison is None or comparison.deg_path is None or not comparison.deg_path.exists():
            raise RuntimeError(f"Missing DEG input for comparison: {comparison_id}")
        self._mark_module_state("deg_explorer", "running", f"Running DEG Explorer for {comparison_id}.")
        self._deg_expression_cache.pop(comparison_id, None)
        frame = self._deg_expression_frame(comparison_id)
        if frame.empty:
            raise RuntimeError(f"DEG input contains no usable rows for comparison: {comparison_id}")
        frame = frame.copy()
        frame.insert(0, "comparison_id", comparison.comparison_id)
        frame.insert(1, "display_name", comparison.display_resolved_name)
        deg_group1, deg_group2 = self._comparison_deg_direction_labels(comparison)
        frame["DEG_direction"] = f"{deg_group1} vs {deg_group2}"
        frame["direction_class"] = "not_significant"
        frame.loc[frame["DEG_significant"] & frame["standardized_log2FC"].gt(0), "direction_class"] = "up"
        frame.loc[frame["DEG_significant"] & frame["standardized_log2FC"].lt(0), "direction_class"] = "down"
        output_dir = files["all"].parent
        output_dir.mkdir(parents=True, exist_ok=True)
        frame.to_csv(files["all"], sep="\t", index=False)
        significant = frame.loc[frame["DEG_significant"].fillna(False)].copy()
        significant.to_csv(files["significant"], sep="\t", index=False)
        significant.loc[significant["standardized_log2FC"].gt(0)].to_csv(files["up"], sep="\t", index=False)
        significant.loc[significant["standardized_log2FC"].lt(0)].to_csv(files["down"], sep="\t", index=False)
        summary = pd.DataFrame(
            [{
                "comparison_id": comparison.comparison_id,
                "display_name": comparison.display_resolved_name,
                "n_tested": len(frame),
                "n_significant": len(significant),
                "n_up": int((significant["standardized_log2FC"] > 0).sum()),
                "n_down": int((significant["standardized_log2FC"] < 0).sum()),
            }]
        )
        summary.to_csv(files["summary"], sep="\t", index=False)
        self._save_deg_explorer_plots(frame, files, comparison.display_resolved_name)
        self._write_json(
            files["parameters"],
            {
                "signature": expected,
                "comparison_id": comparison.comparison_id,
                "display_name": comparison.display_resolved_name,
                "direction": f"log2FC > 0 means expression is higher in {deg_group1} than {deg_group2}",
                "thresholds": {"DE_FDR": self.current_thresholds.deg_padj, "abs_log2FC": self.current_thresholds.deg_log2fc},
                "generated_files": [str(path) for path in files.values()],
            },
        )
        self._mark_module_state(
            "deg_explorer",
            "finished",
            f"Finished DEG Explorer for {comparison_id}.",
            from_cache=False,
            output_folder=str(output_dir),
        )
        return frame

    def _deg_explorer_signature(self, comparison_id: str) -> dict[str, object]:
        if self.current_project is None:
            return {}
        comparison = next(
            (item for item in self.current_project.available_comparisons if item.comparison_id == comparison_id),
            None,
        )
        if comparison is None:
            return {}
        return {
            "cache_version": 1,
            "comparison_id": comparison.comparison_id,
            "display_name": comparison.display_name or "",
            "reverse_deg": bool(comparison.reverse_deg),
            "thresholds": {"deg_padj": self.current_thresholds.deg_padj, "deg_log2fc": self.current_thresholds.deg_log2fc},
            "deg": self.tx_splicing_integration._file_signature(comparison.deg_path, label="DEG"),
        }

    @staticmethod
    def _save_deg_explorer_plots(frame: pd.DataFrame, files: dict[str, Path], title: str) -> None:
        x = pd.to_numeric(frame["standardized_log2FC"], errors="coerce")
        fdr = pd.to_numeric(frame["DE_FDR"], errors="coerce").clip(lower=1e-300)
        significant = frame["DEG_significant"].fillna(False).astype(bool)
        colors = np.where(significant & x.gt(0), "#d95f59", np.where(significant & x.lt(0), "#247f79", "#aab5bd"))
        fig, axis = plt.subplots(figsize=(10, 7))
        axis.scatter(x, -np.log10(fdr), c=colors, s=12, alpha=0.72, linewidths=0)
        axis.axvline(0, color="#596872", linewidth=0.8)
        axis.set_xlabel("standardized log2FC")
        axis.set_ylabel("-log10 adjusted p-value")
        axis.set_title(f"DEG Volcano: {title}", fontweight="bold")
        fig.tight_layout()
        fig.savefig(files["volcano"], dpi=180, bbox_inches="tight")
        plt.close(fig)
        base_mean = pd.to_numeric(frame["baseMean"], errors="coerce").clip(lower=0)
        fig, axis = plt.subplots(figsize=(10, 7))
        axis.scatter(np.log10(base_mean + 1), x, c=colors, s=12, alpha=0.72, linewidths=0)
        axis.axhline(0, color="#596872", linewidth=0.8)
        axis.set_xlabel("log10(baseMean + 1)")
        axis.set_ylabel("standardized log2FC")
        axis.set_title(f"DEG MA Plot: {title}", fontweight="bold")
        fig.tight_layout()
        fig.savefig(files["ma"], dpi=180, bbox_inches="tight")
        plt.close(fig)

    def enrichment_file_map(self, comparison_id: str) -> dict[str, Path]:
        if self.current_project is None or self.current_project.output_root is None:
            return {}
        comparison = next(
            (item for item in self.current_project.available_comparisons if item.comparison_id == comparison_id),
            None,
        )
        if comparison is None:
            return {}
        base = self.current_project.output_root / "05b_functional_enrichment" / (
            comparison.output_prefix or self._safe_output_name(comparison.comparison_id)
        )
        return {
            "results": base / "go_kegg_enrichment.tsv",
            "genes": base / "input_genes.tsv",
            "parameters": base / "parameters.json",
            "plot": base / "go_kegg_enrichment.png",
        }

    def run_functional_enrichment(
        self,
        comparison_id: str,
        *,
        source: str,
        organism: str,
        custom_genes: list[str] | None = None,
        force_rebuild: bool = False,
    ) -> pd.DataFrame:
        if self.current_project is None:
            raise RuntimeError("No project loaded.")
        genes = self._enrichment_input_genes(comparison_id, source, custom_genes or [])
        if not genes:
            raise RuntimeError(
                f"Missing input genes for '{source}'. Run DEG Explorer or Candidate Gene Selection first, or provide a custom gene list."
            )
        files = self.enrichment_file_map(comparison_id)
        signature = {
            "cache_version": 1,
            "comparison_id": comparison_id,
            "source": source,
            "organism": organism,
            "genes_sha256": hashlib.sha256("\n".join(genes).encode("utf-8")).hexdigest(),
            "sources": ["GO:BP", "GO:MF", "GO:CC", "KEGG"],
        }
        if not force_rebuild and files.get("parameters") and files["parameters"].exists() and files["results"].exists():
            try:
                metadata = json.loads(files["parameters"].read_text(encoding="utf-8"))
                if metadata.get("signature") == signature:
                    result = pd.read_csv(files["results"], sep="\t", low_memory=False)
                    self._mark_module_state("functional_enrichment", "cache_hit", "Loaded cached GO/KEGG enrichment.", from_cache=True)
                    return result
            except (OSError, ValueError, json.JSONDecodeError):
                pass
        self._mark_module_state("functional_enrichment", "running", f"Running GO/KEGG enrichment for {comparison_id}.")
        payload = json.dumps(
            {
                "organism": organism,
                "query": genes,
                "sources": ["GO:BP", "GO:MF", "GO:CC", "KEGG"],
                "user_threshold": 0.05,
                "significance_threshold_method": "fdr",
                "no_evidences": False,
            }
        ).encode("utf-8")
        request = urllib.request.Request(
            "https://biit.cs.ut.ee/gprofiler/api/gost/profile/",
            data=payload,
            headers={"Content-Type": "application/json", "User-Agent": "Katana-Splicing-Tool/1.0"},
            method="POST",
        )
        try:
            with urllib.request.urlopen(request, timeout=60) as response:
                data = json.loads(response.read().decode("utf-8"))
        except (urllib.error.URLError, TimeoutError, json.JSONDecodeError) as exc:
            self._mark_module_state("functional_enrichment", "failed", f"GO/KEGG request failed: {exc}")
            raise RuntimeError(f"GO/KEGG enrichment request failed: {exc}") from exc
        rows = []
        for item in data.get("result", []):
            intersections = item.get("intersections", []) or []
            flattened_intersections = []
            for value in intersections:
                if isinstance(value, (list, tuple, set)):
                    flattened_intersections.extend(str(gene) for gene in value if gene)
                elif value:
                    flattened_intersections.append(str(value))
            rows.append(
                {
                    "comparison_id": comparison_id,
                    "source": item.get("source"),
                    "term_id": item.get("native"),
                    "term_name": item.get("name"),
                    "FDR": item.get("p_value"),
                    "term_size": item.get("term_size"),
                    "query_size": item.get("query_size"),
                    "intersection_size": item.get("intersection_size"),
                    "intersection_genes": ";".join(dict.fromkeys(flattened_intersections)),
                    "input_gene_source": source,
                    "organism": organism,
                }
            )
        result = pd.DataFrame(rows)
        output_dir = files["results"].parent
        output_dir.mkdir(parents=True, exist_ok=True)
        result.to_csv(files["results"], sep="\t", index=False)
        pd.DataFrame({"comparison_id": comparison_id, "gene_symbol": genes, "input_gene_source": source}).to_csv(
            files["genes"], sep="\t", index=False
        )
        self._write_json(
            files["parameters"],
            {
                "signature": signature,
                "generated_at": datetime.now().isoformat(timespec="seconds"),
                "generated_files": [str(path) for path in files.values()],
                "note": "Over-representation analysis via the g:Profiler API; enrichment is associative, not causal.",
            },
        )
        self._save_enrichment_plot(result, files["plot"], comparison_id)
        self._mark_module_state(
            "functional_enrichment",
            "finished",
            f"Finished GO/KEGG enrichment for {comparison_id}; {len(result)} terms passed the API threshold.",
            output_folder=str(output_dir),
            from_cache=False,
        )
        return result

    def load_cached_functional_enrichment(self, comparison_id: str) -> pd.DataFrame:
        files = self.enrichment_file_map(comparison_id)
        path = files.get("results")
        if path is None or not path.exists():
            return pd.DataFrame()
        try:
            result = pd.read_csv(path, sep="\t", low_memory=False)
        except Exception:
            return pd.DataFrame()
        self._mark_module_state(
            "functional_enrichment",
            "cache_hit",
            f"Loaded cached GO/KEGG enrichment for {comparison_id}.",
            output_folder=str(path.parent),
            from_cache=True,
        )
        return result

    def _enrichment_input_genes(self, comparison_id: str, source: str, custom_genes: list[str]) -> list[str]:
        if source == "custom":
            return self._dedupe_gene_tokens(custom_genes)
        if source == "candidate_selection":
            selection, _meta = self.downstream_gene_selection(
                comparison_id,
                consumer_module="5.8",
            )
            if selection.empty:
                return []
            values = selection.get(
                "gene_symbol",
                selection.get("gene_id", pd.Series(dtype="object")),
            ).dropna().astype(str).tolist()
            return self._dedupe_gene_tokens(values)
        files = self.deg_explorer_file_map(comparison_id)
        key = {"deg_significant": "significant", "deg_up": "up", "deg_down": "down"}.get(source)
        path = files.get(key or "")
        if path is None or not path.exists():
            return []
        try:
            frame = pd.read_csv(path, sep="\t", low_memory=False, usecols=lambda column: column in {"gene_symbol", "gene_id"})
        except Exception:
            return []
        values = frame.get("gene_symbol", frame.get("gene_id", pd.Series(dtype="object"))).dropna().astype(str).tolist()
        return self._dedupe_gene_tokens(values)

    @staticmethod
    def _save_enrichment_plot(frame: pd.DataFrame, path: Path, comparison_id: str) -> None:
        if frame is None or frame.empty:
            return
        show = frame.copy()
        show["FDR"] = pd.to_numeric(show["FDR"], errors="coerce")
        show = show.dropna(subset=["FDR"]).sort_values("FDR").head(20).sort_values("FDR", ascending=False)
        if show.empty:
            return
        labels = show["term_name"].astype(str).str.slice(0, 58)
        values = -np.log10(show["FDR"].clip(lower=1e-300))
        colors = show["source"].map({"GO:BP": "#247f79", "GO:MF": "#d9903d", "GO:CC": "#4c78a8", "KEGG": "#d95f59"}).fillna("#81909a")
        fig, axis = plt.subplots(figsize=(11, 8))
        axis.barh(labels, values, color=colors)
        axis.set_xlabel("-log10 FDR")
        axis.set_title(f"GO / KEGG enrichment: {comparison_id}", fontweight="bold")
        fig.tight_layout()
        fig.savefig(path, dpi=180, bbox_inches="tight")
        plt.close(fig)

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
        enriched = self._enrich_tx_splicing_table(gene_table)
        self.run_state.results.tx_splicing_gene_table = enriched.copy()
        self.run_state.results.tx_splicing_summary = summary.copy()
        self._mark_module_state(
            "tx_splicing_integration",
            "finished",
            "Transcript + splicing integration loaded/generated for the selected comparisons.",
            output_folder=str(self.current_project.output_root / "04_tx_splicing_integration") if self.current_project.output_root else "",
            from_cache=True,
        )
        return enriched, summary

    def preview_candidate_gene_screening(self, allow_generate: bool = True, *, force_rebuild: bool = False) -> pd.DataFrame:
        if self.current_project is None:
            raise RuntimeError("No project loaded.")
        self._log_cache("candidate ranking request: checking memory/disk cache")
        if not force_rebuild and self.run_state.results.candidate_gene_table is not None and not self.run_state.results.candidate_gene_table.empty:
            cached = self._normalize_candidate_gene_table(self.run_state.results.candidate_gene_table)
            self.run_state.results.candidate_gene_table = cached.copy()
            if not self._candidate_support_table_is_stale(cached):
                self._log_cache("candidate ranking cache hit: in-memory table")
                self._mark_module_state(
                    "candidate_ranking",
                    "cache_hit",
                    "Loaded candidate ranking from in-memory cache.",
                    from_cache=True,
                )
                return cached

        if not force_rebuild:
            cached_disk = self._load_cached_candidate_gene_table()
            if cached_disk is not None and not cached_disk.empty:
                can_rebuild = self.run_state.results.tx_splicing_gene_table is not None and not self.run_state.results.tx_splicing_gene_table.empty
                if not allow_generate or not self._candidate_support_table_is_stale(cached_disk) or not can_rebuild:
                    self.run_state.results.candidate_gene_table = cached_disk.copy()
                    self._log_cache("candidate ranking cache hit: cached TSV output")
                    self._mark_module_state(
                        "candidate_ranking",
                        "cache_hit",
                        "Loaded candidate ranking from cached TSV output.",
                        from_cache=True,
                        output_folder=str(self.current_project.output_root / "04a_candidate_gene_screening") if self.current_project.output_root else "",
                    )
                    return cached_disk.copy()

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
            mechanism = pd.DataFrame()
            self._mark_module_state(
                "mechanism_support",
                "not_run",
                "Mechanism support was not auto-generated by 5.3. Run the mechanism/support module explicitly if DEXSeq/DTU support columns are required.",
            )
        self._log_heavy("candidate ranking rebuild started from cached transcript + splicing integration")
        rebuilt = self._normalize_candidate_gene_table(self._build_candidate_gene_table(integration, mechanism)).copy()
        self.run_state.results.candidate_gene_table = rebuilt.copy()
        self._clear_candidate_selection_caches()
        self._write_step_output(
            "04a_candidate_gene_screening",
            self._candidate_output_tables(),
        )
        self._write_per_comparison_candidate_outputs()
        self._write_candidate_cache_meta()
        self._mark_module_state(
            "candidate_ranking",
            "finished",
            "Candidate ranking built from cached transcript + splicing integration and current support tables, then written to cache/output.",
            from_cache=False,
            output_folder=str(self.current_project.output_root / "02_candidate_gene_screening") if self.current_project.output_root else "",
        )
        return rebuilt

    def _candidate_gene_table_cache_paths(self) -> list[tuple[Path, str | None]]:
        if self.current_project is None or self.current_project.output_root is None:
            return []
        root = self.current_project.output_root
        paths: list[tuple[Path, str | None]] = [
            (root / "04a_candidate_gene_screening" / "data" / "gene_level_integrated_candidates.tsv", None),
            (root / "04a_candidate_gene_screening" / "gene_level_integrated_candidates.tsv", None),
        ]
        for comparison in self.selected_comparisons_for_display():
            path = self.candidate_screening_file_map(comparison.comparison_id).get("gene_level_integrated_candidates")
            if path is not None:
                paths.append((path, comparison.comparison_id))
        return paths

    def _load_cached_candidate_gene_table(self) -> pd.DataFrame:
        expected = self._candidate_cache_signature()
        meta_path = self._candidate_cache_meta_path()
        if meta_path is None or not meta_path.exists():
            self._log_cache("candidate ranking cache miss: parameter fingerprint metadata is missing")
            return pd.DataFrame()
        try:
            cached_meta = json.loads(meta_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            self._log_cache("candidate ranking cache miss: parameter fingerprint metadata is unreadable")
            return pd.DataFrame()
        if cached_meta.get("signature") != expected:
            self._log_cache("candidate ranking cache stale: comparison, direction, threshold, or input fingerprint changed")
            return pd.DataFrame()
        frames: list[pd.DataFrame] = []
        for path, comparison_id in self._candidate_gene_table_cache_paths():
            if not path.exists() or path.stat().st_size <= 0:
                continue
            try:
                frame = pd.read_csv(path, sep="\t", low_memory=False)
            except Exception:
                continue
            if frame.empty:
                continue
            if comparison_id and "comparison_id" not in frame.columns:
                frame["comparison_id"] = comparison_id
            if comparison_id and "comparison_id" in frame.columns:
                frame["comparison_id"] = frame["comparison_id"].fillna(comparison_id).astype(str)
            frames.append(frame)
            if comparison_id is None:
                break
        if not frames:
            return pd.DataFrame()
        self._candidate_cache_fingerprint_in_memory = str(cached_meta.get("fingerprint", ""))
        return self._normalize_candidate_gene_table(pd.concat(frames, axis=0, ignore_index=True, sort=False))

    def _candidate_cache_meta_path(self) -> Path | None:
        if self.current_project is None or self.current_project.output_root is None:
            return None
        return self.current_project.output_root / "04a_candidate_gene_screening" / "candidate_cache_meta.json"

    def _candidate_cache_signature(self) -> dict[str, object]:
        if self.current_project is None:
            return {}
        request = self.build_default_request()
        comparisons = []
        for comparison in self._selected_comparisons():
            comparisons.append(
                {
                    "identity": {
                        "comparison_id": comparison.comparison_id,
                        "display_name": comparison.display_name or "",
                        "experiment_group": comparison.experiment_group or "",
                        "control_group": comparison.control_group or "",
                        "reverse_deg": bool(comparison.reverse_deg),
                        "reverse_splicing": bool(comparison.reverse_splicing),
                        "rmats_name": comparison.rmats_name or "",
                        "deg_name": comparison.deg_name or "",
                        "dexseq_name": comparison.dexseq_name or "",
                        "dtu_name": comparison.dtu_name or "",
                    },
                    "integration": self.tx_splicing_integration._comparison_cache_signature(comparison, request),
                }
            )
        return {
            "cache_version": 4,
            "rmats_mode": self.current_project.selected_rmats_mode,
            "thresholds": asdict(self.current_thresholds),
            "comparison_order": [item.comparison_id for item in self._selected_comparisons()],
            "comparisons": comparisons,
            "input_paths": {key: str(value) for key, value in sorted(self.current_project.input_paths.items())},
        }

    def _write_candidate_cache_meta(self) -> None:
        meta_path = self._candidate_cache_meta_path()
        if meta_path is None:
            return
        signature = self._candidate_cache_signature()
        encoded = json.dumps(signature, sort_keys=True, separators=(",", ":")).encode("utf-8")
        fingerprint = hashlib.sha256(encoded).hexdigest()
        self._write_json(
            meta_path,
            {
                "fingerprint": fingerprint,
                "created_at": datetime.now().isoformat(timespec="seconds"),
                "signature": signature,
            },
        )
        self._candidate_cache_fingerprint_in_memory = fingerprint

    def _normalize_candidate_gene_table(self, frame: pd.DataFrame | None) -> pd.DataFrame:
        if frame is None or frame.empty:
            return pd.DataFrame()
        working = frame.copy()
        for column in [
            "rank",
            "evidence_count",
            "DE_FDR",
            "standardized_log2FC",
            "abs_log2FC",
            "n_rMATS_significant_events",
            "best_rMATS_FDR",
            "max_abs_dPSI",
            "dominant_rMATS_standardized_dPSI",
            "best_DEXSeq_qvalue",
            "n_DEXSeq_significant_exons",
            "best_DTU_qvalue",
            "n_DTU_significant_isoforms",
            "cross_n_comparisons_candidate",
            "cross_best_overall_rank",
            "ranking_score",
        ]:
            if column in working.columns:
                working[column] = pd.to_numeric(working[column], errors="coerce")
        for column in ["DEG_significant", "rMATS_significant", "DEXSeq_significant", "DTU_significant", "shortlist_gene", "blacklist_gene"]:
            if column in working.columns:
                working[column] = working[column].map(self._safe_bool)
        if "comparison_id" not in working.columns:
            working["comparison_id"] = ""
        if "gene_symbol" not in working.columns:
            working["gene_symbol"] = working.get("geneSymbol", working.get("gene_id", pd.Series("", index=working.index)))
        if "gene_id" not in working.columns:
            working["gene_id"] = working.get("normalized_gene_id", working["gene_symbol"])
        if "normalized_gene_id" not in working.columns:
            working["normalized_gene_id"] = working["gene_id"].map(normalize_gene_id)
        if "candidate_unique_key" not in working.columns:
            gene_key = working["gene_symbol"].fillna("").astype(str).str.strip()
            fallback = working["normalized_gene_id"].fillna("").astype(str).str.strip()
            working["candidate_unique_key"] = gene_key.where(gene_key.ne(""), fallback)
        for column in ["significant_as_event_types", "significant_as_event_ids", "significant_as_event_list", "as_event_summary"]:
            if column not in working.columns:
                working[column] = ""
            working[column] = working[column].fillna("").astype(str)
        sort_columns = [column for column in ["comparison_id", "rank", "gene_symbol", "gene_id"] if column in working.columns]
        if sort_columns:
            working = working.sort_values(sort_columns, ascending=[True] * len(sort_columns), na_position="last")
        working = working.drop_duplicates(subset=["comparison_id", "candidate_unique_key"], keep="first").reset_index(drop=True)
        if "rank" in working.columns:
            working["rank"] = working.groupby("comparison_id", dropna=False).cumcount() + 1
        rank_values = pd.to_numeric(
            working.get("rank", pd.Series(range(1, len(working) + 1), index=working.index)),
            errors="coerce",
        )
        if "ranking_score" not in working.columns or pd.to_numeric(working["ranking_score"], errors="coerce").isna().all():
            max_rank = rank_values.groupby(working["comparison_id"].fillna("").astype(str), dropna=False).transform("max")
            with np.errstate(divide="ignore", invalid="ignore"):
                score = ((max_rank - rank_values + 1) / max_rank) * 100.0
            working["ranking_score"] = score.round(3)
        else:
            working["ranking_score"] = pd.to_numeric(working["ranking_score"], errors="coerce")
        return working

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
                "Missing input: no enabled comparison pair is selected for 5.4.3.",
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
                "Loaded cached 5.4.3 cross-comparison significant AS event dPSI pattern analysis.",
                from_cache=True,
                output_folder=str(self.cross_comparison_pattern_output_dir(pair.pair_id) or ""),
            )
            return {key: value.copy() for key, value in cached.items()}

        if not allow_generate:
            self._mark_module_state(
                "cross_comparison_as_patterns",
                "not_run",
                "Not run. Click 'Apply / Recalculate' to build the 5.4.3 cross-comparison significant AS event dPSI pattern tables for the selected comparison pair.",
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
            "Built 5.4.3 cross-comparison significant AS event dPSI pattern tables.",
            from_cache=False,
            output_folder=str(self.cross_comparison_pattern_output_dir(pair.pair_id) or ""),
        )
        self._write_cross_comparison_pattern_outputs(pair, tables, abs_cutoff, delta_cutoff)
        return {key: value.copy() for key, value in tables.items()}

    def preview_cross_comparison_as_interactions(
        self,
        pair_id: str | None,
        *,
        abs_dpsi_cutoff: float = 0.10,
        large_delta_dpsi_cutoff: float = 0.10,
        min_total_reads: int = 20,
        allow_generate: bool = False,
        force_rebuild: bool = False,
    ) -> dict[str, pd.DataFrame]:
        if self.current_project is None:
            raise RuntimeError("No project loaded.")
        pair = next((item for item in self._active_comparison_pairs() if item.pair_id == pair_id), None)
        if pair is None:
            self._mark_module_state(
                "cross_comparison_as_interactions",
                "not_run",
                "Missing input: no enabled comparison pair is selected for 5.4.3.",
            )
            return self._empty_cross_as_interaction_tables()

        abs_cutoff = max(float(abs_dpsi_cutoff), 0.0)
        delta_cutoff = max(float(large_delta_dpsi_cutoff), 0.0)
        min_reads = max(int(min_total_reads), 0)
        cache_key = (pair.pair_id, round(abs_cutoff, 6), round(delta_cutoff, 6), min_reads)
        cached = None if force_rebuild else self._cross_as_interaction_cache.get(cache_key)
        if cached is not None:
            self._mark_module_state(
                "cross_comparison_as_interactions",
                "cache_hit",
                "Loaded cached 5.4.3 cross-comparison AS interaction modeling tables.",
                from_cache=True,
                output_folder=str(self.cross_comparison_interaction_output_dir(pair.pair_id) or ""),
            )
            return {key: value.copy() for key, value in cached.items()}

        if not allow_generate:
            self._mark_module_state(
                "cross_comparison_as_interactions",
                "not_run",
                "Not run. Click 'Apply / Recalculate' to build 5.4.3 for the selected comparison pair.",
            )
            return self._empty_cross_as_interaction_tables()

        tables = self._build_cross_comparison_as_interactions(
            pair,
            abs_dpsi_cutoff=abs_cutoff,
            large_delta_dpsi_cutoff=delta_cutoff,
            min_total_reads=min_reads,
        )
        self._cross_as_interaction_cache[cache_key] = {key: value.copy() for key, value in tables.items()}
        self._mark_module_state(
            "cross_comparison_as_interactions",
            "finished",
            "Built 5.4.3 cross-comparison AS interaction modeling tables.",
            from_cache=False,
            output_folder=str(self.cross_comparison_interaction_output_dir(pair.pair_id) or ""),
        )
        self._write_cross_comparison_interaction_outputs(pair, tables, abs_cutoff, delta_cutoff, min_reads)
        return {key: value.copy() for key, value in tables.items()}

    def preview_cross_comparison_as_interaction_batch(
        self,
        *,
        abs_dpsi_cutoff: float = 0.10,
        large_delta_dpsi_cutoff: float = 0.10,
        min_total_reads: int = 20,
        allow_generate: bool = False,
        force_rebuild: bool = False,
    ) -> dict[str, pd.DataFrame]:
        if self.current_project is None:
            raise RuntimeError("No project loaded.")
        pairs = self._interaction_batch_pairs()
        if not pairs:
            self._mark_module_state(
                "cross_comparison_as_interactions",
                "not_run",
                "Missing input: at least two selected comparisons are required for all-pair 5.4.3 batch summary.",
            )
            return self._empty_cross_as_interaction_batch_tables()

        abs_cutoff = max(float(abs_dpsi_cutoff), 0.0)
        delta_cutoff = max(float(large_delta_dpsi_cutoff), 0.0)
        min_reads = max(int(min_total_reads), 0)
        pair_signature = tuple(pair.pair_id for pair in pairs)
        cache_key = (pair_signature, round(abs_cutoff, 6), round(delta_cutoff, 6), min_reads)
        cached = None if force_rebuild else self._cross_as_interaction_batch_cache.get(cache_key)
        if cached is not None:
            self._mark_module_state(
                "cross_comparison_as_interactions",
                "cache_hit",
                "Loaded cached all-pair 5.4.3 AS interaction batch summary.",
                from_cache=True,
                output_folder=str(self.cross_comparison_interaction_batch_output_dir() or ""),
            )
            return {key: value.copy() for key, value in cached.items()}

        if not allow_generate:
            self._mark_module_state(
                "cross_comparison_as_interactions",
                "not_run",
                "Not run. Select an all-pair 5.4.3 view and click 'Apply / Recalculate' to build the batch summary.",
            )
            return self._empty_cross_as_interaction_batch_tables()

        summary_rows: list[dict[str, object]] = []
        top_frames: list[pd.DataFrame] = []
        exclusive_frames: list[pd.DataFrame] = []
        benchmark_frames: list[pd.DataFrame] = []
        diagnostic_frames: list[pd.DataFrame] = []
        for order, pair in enumerate(pairs, start=1):
            tables = self._build_cross_comparison_as_interactions(
                pair,
                abs_dpsi_cutoff=abs_cutoff,
                large_delta_dpsi_cutoff=delta_cutoff,
                min_total_reads=min_reads,
            )
            pair_cache_key = (pair.pair_id, round(abs_cutoff, 6), round(delta_cutoff, 6), min_reads)
            self._cross_as_interaction_cache[pair_cache_key] = {key: value.copy() for key, value in tables.items()}
            self._write_cross_comparison_interaction_outputs(pair, tables, abs_cutoff, delta_cutoff, min_reads)

            shared = tables.get("shared_interaction", pd.DataFrame()).copy()
            exclusive = tables.get("exclusive", pd.DataFrame()).copy()
            benchmark = tables.get("benchmark", pd.DataFrame()).copy()
            diagnostics = tables.get("diagnostics", pd.DataFrame()).copy()
            if not diagnostics.empty:
                diagnostics.insert(0, "batch_pair_order", order)
                diagnostic_frames.append(diagnostics)
            if not shared.empty:
                shared.insert(0, "batch_pair_order", order)
                shared.insert(1, "pair_id", pair.pair_id)
                shared.insert(2, "pair_name", pair.resolved_name)
                top = shared.sort_values(
                    ["interaction_priority_score", "abs_delta_dPSI", "shrunken_interaction_pvalue"],
                    ascending=[False, False, True],
                    na_position="last",
                ).head(100)
                top_frames.append(top)
            if not exclusive.empty:
                exclusive.insert(0, "batch_pair_order", order)
                exclusive.insert(1, "pair_id", pair.pair_id)
                exclusive.insert(2, "pair_name", pair.resolved_name)
                exclusive_frames.append(exclusive)
            if not benchmark.empty:
                benchmark.insert(0, "batch_pair_order", order)
                benchmark.insert(1, "pair_id", pair.pair_id)
                benchmark.insert(2, "pair_name", pair.resolved_name)
                benchmark_frames.append(benchmark)

            high_confidence = (
                shared.get("interaction_class", pd.Series(dtype="object"))
                .astype(str)
                .eq("high_confidence_interaction_candidate")
                .sum()
                if not shared.empty else 0
            )
            moderate = (
                shared.get("interaction_class", pd.Series(dtype="object"))
                .astype(str)
                .eq("moderate_interaction_candidate")
                .sum()
                if not shared.empty else 0
            )
            top_event = shared.sort_values("interaction_priority_score", ascending=False, na_position="last").head(1) if not shared.empty else pd.DataFrame()
            summary_rows.append(
                {
                    "batch_pair_order": order,
                    "pair_id": pair.pair_id,
                    "pair_name": pair.resolved_name,
                    "comparison_A": pair.comparison_a,
                    "comparison_B": pair.comparison_b,
                    "n_shared_significant": len(shared),
                    "n_high_confidence_interaction": int(high_confidence),
                    "n_moderate_interaction": int(moderate),
                    "n_opposite_direction": len(tables.get("opposite_direction", pd.DataFrame())),
                    "n_same_direction_large_delta": len(tables.get("large_delta", pd.DataFrame())),
                    "n_exclusive": len(exclusive),
                    "n_exclusive_A": int(exclusive.get("exclusive_side", pd.Series(dtype="object")).astype(str).eq("comparison_A_only").sum()) if not exclusive.empty else 0,
                    "n_exclusive_B": int(exclusive.get("exclusive_side", pd.Series(dtype="object")).astype(str).eq("comparison_B_only").sum()) if not exclusive.empty else 0,
                    "n_exclusive_supported": int(
                        exclusive.get("non_significant_side_read_support", pd.Series(dtype="object"))
                        .astype(str)
                        .str.startswith("supported_")
                        .sum()
                    ) if not exclusive.empty else 0,
                    "n_exclusive_A_supported": int(
                        (
                            exclusive.get("exclusive_side", pd.Series(dtype="object")).astype(str).eq("comparison_A_only")
                            & exclusive.get("non_significant_side_read_support", pd.Series(dtype="object")).astype(str).str.startswith("supported_")
                        ).sum()
                    ) if not exclusive.empty else 0,
                    "n_exclusive_B_supported": int(
                        (
                            exclusive.get("exclusive_side", pd.Series(dtype="object")).astype(str).eq("comparison_B_only")
                            & exclusive.get("non_significant_side_read_support", pd.Series(dtype="object")).astype(str).str.startswith("supported_")
                        ).sum()
                    ) if not exclusive.empty else 0,
                    "n_exclusive_low_power": int(
                        exclusive.get("non_significant_side_read_support", pd.Series(dtype="object"))
                        .astype(str)
                        .str.startswith("low_power")
                        .sum()
                    ) if not exclusive.empty else 0,
                    "n_exclusive_counterpart_not_found": int(
                        exclusive.get("non_significant_side_read_support", pd.Series(dtype="object"))
                        .astype(str)
                        .str.startswith("counterpart_event_not_found")
                        .sum()
                    ) if not exclusive.empty else 0,
                    "top_gene": str(top_event.iloc[0].get("gene", "")) if not top_event.empty else "",
                    "top_event_type": str(top_event.iloc[0].get("event_type", "")) if not top_event.empty else "",
                    "top_abs_delta_dPSI": top_event.iloc[0].get("abs_delta_dPSI", pd.NA) if not top_event.empty else pd.NA,
                    "top_interaction_priority_score": top_event.iloc[0].get("interaction_priority_score", pd.NA) if not top_event.empty else pd.NA,
                    "delta_top20_model_top20_overlap": diagnostics.iloc[0].get("delta_top20_model_top20_overlap", pd.NA) if not diagnostics.empty else pd.NA,
                    "delta_vs_model_spearman_rank_correlation": diagnostics.iloc[0].get("delta_vs_model_spearman_rank_correlation", pd.NA) if not diagnostics.empty else pd.NA,
                    "abs_dPSI_significance_cutoff_used": abs_cutoff,
                    "large_delta_dPSI_cutoff_used": delta_cutoff,
                    "min_total_reads_used": min_reads,
                }
            )

        tables = {
            "batch_summary": pd.DataFrame(summary_rows),
            "batch_top_interactions": pd.concat(top_frames, axis=0, ignore_index=True) if top_frames else pd.DataFrame(),
            "batch_exclusive": pd.concat(exclusive_frames, axis=0, ignore_index=True) if exclusive_frames else pd.DataFrame(),
            "batch_benchmark": pd.concat(benchmark_frames, axis=0, ignore_index=True) if benchmark_frames else pd.DataFrame(),
            "batch_diagnostics": pd.concat(diagnostic_frames, axis=0, ignore_index=True) if diagnostic_frames else pd.DataFrame(),
        }
        self._cross_as_interaction_batch_cache[cache_key] = {key: value.copy() for key, value in tables.items()}
        self._write_cross_comparison_interaction_batch_outputs(tables, abs_cutoff, delta_cutoff, min_reads)
        self._mark_module_state(
            "cross_comparison_as_interactions",
            "finished",
            "Built all-pair 5.4.3 AS interaction batch summary.",
            from_cache=False,
            output_folder=str(self.cross_comparison_interaction_batch_output_dir() or ""),
        )
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
        self.current_project.output_root = self.output_directories.require(output_root)
        self._artifact_registry = None
        self.save_project_config()

    def can_migrate_current_output_root(self) -> bool:
        return bool(
            self.current_project is not None
            and self.output_directories.is_recognized_output_root(self.current_project.output_root)
        )

    def migrate_output_root(
        self,
        output_root: str | Path,
        *,
        remove_source: bool = True,
    ) -> OutputMigrationResult:
        if self.current_project is None:
            raise RuntimeError("No project loaded.")
        source = self.current_project.output_root
        if source is None:
            raise RuntimeError("The current project has no existing output directory to migrate.")
        result = self.output_directories.migrate(
            source,
            output_root,
            remove_source=remove_source,
        )
        self.current_project.output_root = result.destination
        self._artifact_registry = None
        self.save_project_config()
        return result

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

    def run_native_sashimi_pipeline(self, manifest: pd.DataFrame | None = None) -> str:
        if self.current_project is None:
            raise RuntimeError("No project loaded.")
        output_root = self.require_output_root()
        active_manifest = manifest.copy() if manifest is not None else self.run_state.results.sashimi_manifest.copy()
        if active_manifest.empty:
            raise RuntimeError("No selected sashimi manifest is available. Please select events and generate the manifest first.")
        gtf_path = self._resolve_gtf_path()
        samtools_path = self.current_project.tool_paths.get("samtools") or "samtools"
        self.native_sashimi_adapter = NativeSashimiAdapter(samtools_path)
        logs, failures = self.native_sashimi_adapter.run_manifest(active_manifest, output_root, gtf_path=gtf_path)
        self._last_sashimi_failures = failures.copy()
        if failures.empty:
            self._mark_module_state(
                "sashimi_plot",
                "finished",
                "Finished Katana native sashimi for the selected events.",
                output_folder=str(output_root / "06_sashimi"),
                from_cache=False,
            )
        else:
            self._mark_module_state(
                "sashimi_plot",
                "failed",
                f"Katana native sashimi finished with {len(failures)} failed event(s).",
                output_folder=str(output_root / "06_sashimi"),
                from_cache=False,
            )
            failed_dir = output_root / "06_sashimi"
            failed_dir.mkdir(parents=True, exist_ok=True)
            failures.to_csv(failed_dir / "failed_native_sashimi.tsv", sep="\t", index=False)
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
        fasta_file = self._resolve_fasta_path("transcript") or self._resolve_fasta_path("any")
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
        configured_samtools = self.current_project.tool_paths.get("samtools") if self.current_project else ""
        samtools_available = bool(configured_samtools and Path(configured_samtools).exists()) or shutil.which("samtools") is not None
        lines.append(f"pysam: {'available' if pysam_available else 'missing'}")
        lines.append(f"samtools: {configured_samtools or shutil.which('samtools') or 'missing'}")
        lines.append(f"Katana native sashimi: {'ready' if samtools_available else 'not ready (samtools missing)'}")
        lines.append(f"Rscript: {self.current_project.tool_paths.get('rscript', 'not set') if self.current_project else 'not set'}")
        lines.append(f"GTF: {self.current_project.tool_paths.get('gtf', 'not set') if self.current_project else 'not set'}")
        lines.append(f"Annotation TSV: {self.current_project.tool_paths.get('annotation_tsv', 'not set') if self.current_project else 'not set'}")
        transcript_fasta = self._resolve_fasta_path("transcript") if self.current_project else None
        genome_fasta = self._resolve_fasta_path("genome") if self.current_project else None
        lines.append(f"Transcript FASTA: {transcript_fasta or 'not set'}")
        lines.append(f"Genome FASTA: {genome_fasta or 'not set'}")
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
        output_root = self.require_output_root()
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
        if (
            self.current_project is not None
            and report.get("status") in {"idle", "not_run", ""}
            and not report.get("generated_files")
        ):
            output_root = self.require_output_root()
            base = output_root / "08_jutils"
            if base.exists():
                browser = self.build_jutils_output_browser()
                if not browser.empty:
                    report["status"] = "cache_hit"
                    report["output_directory"] = str(base)
                    report["generated_files"] = browser.to_dict("records")
        return report

    def build_jutils_output_browser(self) -> pd.DataFrame:
        if self.current_project is None:
            raise RuntimeError("No project loaded.")
        output_root = self.require_output_root()
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

    def native_jutils_report(self, mode: str = "pattern") -> dict[str, object]:
        report = dict(self._last_native_jutils_report)
        if str(report.get("mode", "")) != str(mode):
            report = {
                "status": "idle",
                "mode": mode,
                "output_directory": "",
                "generated_files": [],
                "stdout": "",
                "stderr": "",
                "error_log": "",
            }
        report.setdefault("status", "idle")
        report.setdefault("mode", mode)
        report.setdefault("output_directory", "")
        report.setdefault("generated_files", [])
        report.setdefault("stdout", "")
        report.setdefault("stderr", "")
        report.setdefault("error_log", "")
        if (
            self.current_project is not None
            and report.get("status") in {"idle", "not_run", ""}
            and not report.get("generated_files")
        ):
            browser = self.build_native_jutils_output_browser(mode)
            if not browser.empty:
                output_root = self.require_output_root()
                base = output_root / "06_native_splicing_patterns" / ("candidate_gene_heatmap" if mode == "candidate" else "pattern_explorer")
                report["status"] = "cache_hit"
                report["mode"] = mode
                report["output_directory"] = str(base)
                report["generated_files"] = browser.to_dict("records")
        return report

    def native_pattern_summary(self, mode: str = "pattern") -> dict[str, object]:
        base = self._native_jutils_output_dir(mode)
        summary: dict[str, object] = {
            "status": "not_run",
            "comparison_count": 0,
            "event_count": 0,
            "sample_count": 0,
            "pattern_cluster_count": 0,
            "pair_count": 0,
            "median_direction_similarity": float("nan"),
            "median_weighted_cosine_similarity": float("nan"),
            "low_confidence_pair_count": 0,
            "output_directory": str(base),
        }
        parameters_path = base / "parameters.json"
        if not parameters_path.exists() or not self._native_jutils_cache_schema_matches(base):
            return summary
        try:
            parameters = json.loads(parameters_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            return summary
        summary["status"] = "cache_hit"
        summary["comparison_count"] = len(parameters.get("comparison_ids", []) or [])

        cached_summary = self._read_table_or_empty(base / "native_pattern_summary.tsv")
        if not cached_summary.empty:
            row = cached_summary.iloc[0].to_dict()
            for key in summary:
                if key in row and key not in {"status", "output_directory"}:
                    summary[key] = row[key]
            return summary

        def read_table(name: str) -> pd.DataFrame:
            return self._read_table_or_empty(base / name)

        events = read_table("unified_as_event_table.tsv")
        psi = read_table("sample_level_psi_matrix.tsv")
        clusters = read_table("psi_pattern_cluster_summary.tsv")
        similarity = read_table("pairwise_as_similarity.tsv")
        summary["event_count"] = int(len(events))
        summary["sample_count"] = max(int(psi.shape[1] - 1), 0) if not psi.empty else 0
        summary["pattern_cluster_count"] = int(len(clusters))
        summary["pair_count"] = int(len(similarity))
        if not similarity.empty:
            direction = pd.to_numeric(similarity.get("direction_similarity"), errors="coerce")
            cosine = pd.to_numeric(similarity.get("weighted_cosine_similarity"), errors="coerce")
            summary["median_direction_similarity"] = float(direction.median()) if direction.notna().any() else float("nan")
            summary["median_weighted_cosine_similarity"] = float(cosine.median()) if cosine.notna().any() else float("nan")
            confidence = similarity.get("similarity_confidence", pd.Series("", index=similarity.index)).astype(str)
            summary["low_confidence_pair_count"] = int(confidence.str.contains("low|no_shared", case=False, regex=True).sum())
        return summary

    def _register_native_pattern_outputs(
        self,
        *,
        mode: str,
        generated: list[Path],
        comparison_ids: tuple[str, ...],
        parameters: dict[str, object],
    ) -> None:
        module_id = "5.3.5" if mode == "candidate" else "5.6"
        source_artifact_ids: list[str] = []
        if mode == "candidate":
            for comparison_id in comparison_ids:
                selection = self.artifact_registry().find(
                    module_id="5.3.2",
                    artifact_key="candidate_gene_selection",
                    scope_id=comparison_id,
                )
                if selection:
                    source_artifact_ids.append(selection[0].artifact_id)
        for path in generated:
            suffix = path.suffix.lower()
            kind = (
                "plot"
                if suffix in {".png", ".jpg", ".jpeg", ".svg", ".pdf"}
                else "table"
                if suffix in {".tsv", ".csv", ".xlsx"}
                else "parameters"
                if suffix == ".json"
                else "report"
                if suffix in {".html", ".htm"}
                else "other"
            )
            entity_type = (
                "event"
                if "event" in path.name.lower() or "psi" in path.name.lower() or "dpsi" in path.name.lower()
                else "analysis"
            )
            self.artifact_registry().register(
                module_id=module_id,
                artifact_key=path.relative_to(self._native_jutils_output_dir(mode)).as_posix(),
                kind=kind,
                entity_type=entity_type,
                path=path,
                scope="comparison_set",
                scope_id="|".join(comparison_ids) or "project",
                comparison_ids=comparison_ids,
                source_artifact_ids=source_artifact_ids,
                consumer_modules=["5.6.4", "9"] if mode == "pattern" else ["5.3.4", "5.5"],
                description=(
                    "Candidate-selected splicing pattern artifact."
                    if mode == "candidate"
                    else "Global selected-comparison splicing pattern artifact."
                ),
                parameters=parameters,
            )

    def build_native_jutils_output_browser(self, mode: str = "pattern") -> pd.DataFrame:
        if self.current_project is None:
            raise RuntimeError("No project loaded.")
        output_root = self.require_output_root()
        base = output_root / "06_native_splicing_patterns" / ("candidate_gene_heatmap" if mode == "candidate" else "pattern_explorer")
        columns = [
            "result_group",
            "preview_priority",
            "comparison_id",
            "relative_path",
            "absolute_path",
            "file_type",
            "file_size",
            "last_modified",
            "preview_kind",
            "previewable",
        ]
        if not base.exists():
            return pd.DataFrame(columns=columns)
        if not self._native_jutils_cache_schema_matches(base):
            return pd.DataFrame(columns=columns)
        rows: list[dict[str, object]] = []
        for path in sorted(base.rglob("*")):
            if not path.is_file():
                continue
            rel = path.relative_to(base)
            stat = path.stat()
            suffix = path.suffix.lower()
            preview_kind = self._preview_kind_for_path(path)
            result_group, preview_priority = self._native_jutils_browser_priority(str(rel).replace("\\", "/"), mode)
            comparison_scope = self._native_jutils_browser_scope(str(rel).replace("\\", "/"))
            rows.append(
                {
                    "result_group": result_group,
                    "preview_priority": preview_priority,
                    "comparison_id": comparison_scope,
                    "relative_path": str(rel),
                    "absolute_path": str(path),
                    "file_type": suffix.lstrip(".") or "file",
                    "file_size": stat.st_size,
                    "last_modified": datetime.fromtimestamp(stat.st_mtime).strftime("%Y-%m-%d %H:%M:%S"),
                    "preview_kind": preview_kind,
                    "previewable": preview_kind != "binary",
                }
            )
        if not rows:
            return pd.DataFrame(columns=columns)
        frame = pd.DataFrame(rows, columns=columns)
        return frame.sort_values(["preview_priority", "relative_path"], ascending=[True, True]).reset_index(drop=True)

    @staticmethod
    def _native_jutils_browser_scope(relative_path: str) -> str:
        normalized = relative_path.replace("\\", "/")
        parts = normalized.split("/")
        if len(parts) >= 3 and parts[0] == "plots" and parts[1] == "per_comparison":
            stem = Path(parts[-1]).stem
            for suffix in (
                "_sample_psi_heatmap_clustered",
                "_sample_psi_heatmap",
                "_sample_psi_pca",
            ):
                if stem.endswith(suffix):
                    stem = stem[: -len(suffix)]
                    break
            return f"per_comparison:{stem}"
        return "cross_comparison"

    @staticmethod
    def _native_jutils_browser_priority(relative_path: str, mode: str) -> tuple[str, int]:
        path = relative_path.lower()
        if mode == "candidate":
            priorities = [
                ("plots/native_candidate_gene_as_summary_heatmap.png", "1. Candidate gene AS summary", 0),
                ("candidate_gene_as_summary.tsv", "1. Candidate gene AS summary", 1),
                ("candidate_gene_as_matrix.tsv", "1. Candidate gene AS summary", 2),
                ("plots/per_comparison/", "2. Per-comparison sample PSI", 20),
                ("plots/native_as_dpsi_pattern_heatmap.png", "3. Candidate event dPSI", 30),
                ("dpsi_matrix.tsv", "3. Candidate event dPSI", 31),
            ]
        else:
            priorities = [
                ("plots/per_comparison/", "1. Per-comparison PSI heatmap / PCA", 0),
                ("per_comparison_sample_psi_pca.tsv", "1. Per-comparison PSI heatmap / PCA", 1),
                ("plots/native_sample_level_psi_heatmap.png", "2. Absolute sample PSI heatmap", 20),
                ("sample_level_psi_matrix.tsv", "2. Absolute sample PSI heatmap", 21),
                ("plots/native_sample_level_psi_zscore_heatmap.png", "3. Z-score PSI heatmap / pattern clusters", 30),
                ("plots/native_sample_psi_pattern_clusters.png", "3. Z-score PSI heatmap / pattern clusters", 31),
                ("psi_pattern_clusters.tsv", "3. Z-score PSI heatmap / pattern clusters", 32),
                ("psi_pattern_cluster_summary.tsv", "3. Z-score PSI heatmap / pattern clusters", 33),
                ("plots/native_sample_psi_pca.png", "4. Sample PSI PCA / replicate QC", 40),
                ("sample_psi_pca.tsv", "4. Sample PSI PCA / replicate QC", 41),
                ("plots/native_replicate_group_separation.png", "4. Sample PSI PCA / replicate QC", 42),
                ("replicate_group_separation.tsv", "4. Sample PSI PCA / replicate QC", 43),
                ("plots/native_sample_psi_similarity_heatmap.png", "5. Sample PSI similarity + denominator", 50),
                ("plots/native_sample_psi_shared_event_counts.png", "5. Sample PSI similarity + denominator", 51),
                ("sample_psi_similarity.tsv", "5. Sample PSI similarity + denominator", 52),
                ("sample_psi_spearman_matrix.tsv", "5. Sample PSI similarity + denominator", 53),
                ("sample_psi_shared_event_count_matrix.tsv", "5. Sample PSI similarity + denominator", 54),
                ("plots/native_as_similarity_heatmap.png", "6. Comparison significant-event overlap", 60),
                ("plots/native_as_direction_similarity_heatmap.png", "6. Comparison significant-event overlap", 61),
                ("pairwise_as_similarity.tsv", "6. Comparison significant-event overlap", 62),
                ("plots/native_as_dpsi_pattern_heatmap.png", "7. Cross-comparison dPSI audit", 70),
                ("dpsi_matrix.tsv", "7. Cross-comparison dPSI audit", 71),
            ]
        for marker, group, priority in priorities:
            if path.endswith(marker) or marker.endswith("/") and marker in path:
                return group, priority
        if path.endswith("parameters.json"):
            return "9. Parameters / reproducibility", 900
        if path.endswith(".tsv"):
            return "8. Source tables", 800
        return "8. Other outputs", 850

    def run_native_jutils_pattern_explorer(
        self,
        *,
        fdr_cutoff: float = 0.05,
        abs_dpsi_cutoff: float = 0.05,
        top_events: int = 100,
        selection_mode: str = "supervised",
    ) -> dict[str, pd.DataFrame]:
        if self.current_project is None:
            raise RuntimeError("No project loaded.")
        comparisons = self._selected_comparisons()
        comparison_ids = tuple(comparison.comparison_id for comparison in comparisons)
        selection_mode = "unsupervised" if str(selection_mode).lower().startswith("unsupervised") else "supervised"
        cache_key = (
            self.NATIVE_PATTERN_SCHEMA_VERSION,
            "pattern",
            comparison_ids,
            float(fdr_cutoff),
            float(abs_dpsi_cutoff),
            str(top_events),
            selection_mode,
        )
        self._log_heavy("native Jutils pattern explorer rebuild started")
        tables = self._build_native_jutils_tables(
            comparisons,
            fdr_cutoff=fdr_cutoff,
            abs_dpsi_cutoff=abs_dpsi_cutoff,
            top_events=top_events,
            selection_mode=selection_mode,
            gene_filter=None,
        )
        tables["native_pattern_summary"] = self._native_pattern_summary_table(
            tables,
            comparison_ids=comparison_ids,
            fdr_cutoff=fdr_cutoff,
            abs_dpsi_cutoff=abs_dpsi_cutoff,
            top_events=top_events,
            selection_mode=selection_mode,
        )
        self._native_jutils_cache[cache_key] = {key: value.copy() for key, value in tables.items()}
        base = self._native_jutils_output_dir("pattern")
        output_parameters = {
            "mode": "pattern_explorer",
            "analysis_schema_version": self.NATIVE_PATTERN_SCHEMA_VERSION,
            "description": "Native selected-comparison AS pattern overview. Support-prioritization only; not a causal model.",
            "fdr_cutoff": fdr_cutoff,
            "abs_dpsi_cutoff": abs_dpsi_cutoff,
            "top_events": top_events,
            "selection_mode": selection_mode,
            "comparison_ids": list(comparison_ids),
            "heatmap_row_order": "hierarchical row clustering with weighted linkage and cityblock distance when scipy is available; deterministic signal ordering fallback",
            "heatmap_column_order": "biological sample order for the absolute PSI matrix; clustered views may reorder samples by PSI similarity",
            "pca_method": "SVD PCA on the event x biological-sample PSI matrix after feature standardization and mean imputation",
            "cluster_method": "weighted-linkage hierarchical clustering on row-z-scored sample/replicate PSI trajectories",
            "similarity_rule": "sample PSI Spearman correlation is always reported with the shared measured-event count; significant-event overlap and direction agreement are separate comparison-level metrics",
        }
        generated = self._write_native_jutils_outputs(
            base,
            tables,
            output_parameters,
        )
        self._register_native_pattern_outputs(
            mode="pattern",
            generated=generated,
            comparison_ids=comparison_ids,
            parameters=output_parameters,
        )
        browser = self.build_native_jutils_output_browser("pattern")
        self._last_native_jutils_report = {
            "status": "finished",
            "mode": "pattern",
            "output_directory": str(base),
            "generated_files": browser.to_dict("records"),
            "stdout": "\n".join(str(path) for path in generated),
            "stderr": "",
            "error_log": "",
        }
        return tables

    def run_native_jutils_candidate_gene_heatmap(
        self,
        *,
        gene_text: str = "",
        fdr_cutoff: float = 0.05,
        abs_dpsi_cutoff: float = 0.05,
        top_events: int = 100,
        selection_mode: str = "supervised",
    ) -> dict[str, pd.DataFrame]:
        if self.current_project is None:
            raise RuntimeError("No project loaded.")
        comparisons = self._selected_comparisons()
        gene_map = self._native_jutils_gene_filter_map(gene_text, comparisons)
        genes = set(gene_map)
        if not genes:
            raise RuntimeError("No candidate genes provided or selected. Enter gene symbols/IDs or run 5.3.2 candidate selection first.")
        comparison_ids = tuple(comparison.comparison_id for comparison in comparisons)
        selection_mode = "unsupervised" if str(selection_mode).lower().startswith("unsupervised") else "supervised"
        cache_key = (
            self.NATIVE_PATTERN_SCHEMA_VERSION,
            "candidate",
            comparison_ids,
            float(fdr_cutoff),
            float(abs_dpsi_cutoff),
            str(top_events),
            selection_mode,
            "|".join(sorted(genes)),
        )
        self._log_heavy("native Jutils candidate gene AS heatmap rebuild started")
        tables = self._build_native_jutils_tables(
            comparisons,
            fdr_cutoff=fdr_cutoff,
            abs_dpsi_cutoff=abs_dpsi_cutoff,
            top_events=top_events,
            selection_mode=selection_mode,
            gene_filter=genes,
            gene_filter_display=gene_map,
        )
        tables["native_pattern_summary"] = self._native_pattern_summary_table(
            tables,
            comparison_ids=comparison_ids,
            fdr_cutoff=fdr_cutoff,
            abs_dpsi_cutoff=abs_dpsi_cutoff,
            top_events=top_events,
            selection_mode=selection_mode,
        )
        self._native_jutils_cache[cache_key] = {key: value.copy() for key, value in tables.items()}
        base = self._native_jutils_output_dir("candidate")
        output_parameters = {
            "mode": "candidate_gene_as_heatmap",
            "analysis_schema_version": self.NATIVE_PATTERN_SCHEMA_VERSION,
            "description": "Native candidate-gene AS heatmap. Uses selected genes only; not part of global AS pattern overview.",
            "fdr_cutoff": fdr_cutoff,
            "abs_dpsi_cutoff": abs_dpsi_cutoff,
            "top_events": top_events,
            "selection_mode": selection_mode,
            "genes": list(gene_map.values()),
            "comparison_ids": list(comparison_ids),
            "heatmap_row_order": "hierarchical row clustering with weighted linkage and cityblock distance when scipy is available; deterministic signal ordering fallback",
            "heatmap_column_order": "biological sample order for the absolute PSI matrix; clustered views may reorder samples by PSI similarity",
            "pca_method": "SVD PCA on the event x biological-sample PSI matrix after feature standardization and mean imputation",
            "cluster_method": "weighted-linkage hierarchical clustering on row-z-scored sample/replicate PSI trajectories",
            "similarity_rule": "sample PSI Spearman correlation is always reported with the shared measured-event count; significant-event overlap and direction agreement are separate comparison-level metrics",
        }
        generated = self._write_native_jutils_outputs(
            base,
            tables,
            output_parameters,
        )
        self._register_native_pattern_outputs(
            mode="candidate",
            generated=generated,
            comparison_ids=comparison_ids,
            parameters=output_parameters,
        )
        browser = self.build_native_jutils_output_browser("candidate")
        self._last_native_jutils_report = {
            "status": "finished",
            "mode": "candidate",
            "output_directory": str(base),
            "generated_files": browser.to_dict("records"),
            "stdout": "\n".join(str(path) for path in generated),
            "stderr": "",
            "error_log": "",
        }
        return tables

    def _native_jutils_output_dir(self, mode: str) -> Path:
        if self.current_project is None:
            raise RuntimeError("No project loaded.")
        output_root = self.require_output_root()
        return output_root / "06_native_splicing_patterns" / ("candidate_gene_heatmap" if mode == "candidate" else "pattern_explorer")

    @staticmethod
    def _native_pattern_summary_table(
        tables: dict[str, pd.DataFrame],
        *,
        comparison_ids: tuple[str, ...],
        fdr_cutoff: float,
        abs_dpsi_cutoff: float,
        top_events: int,
        selection_mode: str,
    ) -> pd.DataFrame:
        events = tables.get("unified_as_event_table", pd.DataFrame())
        psi = tables.get("sample_level_psi_matrix", pd.DataFrame())
        clusters = tables.get("psi_pattern_cluster_summary", pd.DataFrame())
        similarity = tables.get("pairwise_as_similarity", pd.DataFrame())
        direction = pd.to_numeric(
            similarity.get(
                "direction_similarity",
                pd.Series(float("nan"), index=similarity.index, dtype="float64"),
            ),
            errors="coerce",
        )
        cosine = pd.to_numeric(
            similarity.get(
                "weighted_cosine_similarity",
                pd.Series(float("nan"), index=similarity.index, dtype="float64"),
            ),
            errors="coerce",
        )
        confidence = similarity.get(
            "similarity_confidence",
            pd.Series("", index=similarity.index, dtype="object"),
        ).astype(str)
        measured_pairs = int(len(similarity))
        low_confidence_pairs = int(confidence.str.contains("low|no_shared", case=False, regex=True).sum())
        interpretation = (
            "No comparison-pair similarity could be estimated."
            if measured_pairs == 0
            else (
                f"{low_confidence_pairs} of {measured_pairs} comparison pairs have low or insufficient shared-event "
                "support; interpret similarity together with its measured-event denominator."
            )
        )
        return pd.DataFrame(
            [
                {
                    "comparison_count": len(comparison_ids),
                    "comparison_ids": "|".join(comparison_ids),
                    "event_count": int(len(events)),
                    "sample_count": max(int(psi.shape[1] - 1), 0) if not psi.empty else 0,
                    "pattern_cluster_count": int(len(clusters)),
                    "pair_count": measured_pairs,
                    "median_direction_similarity": (
                        float(direction.median()) if direction.notna().any() else float("nan")
                    ),
                    "median_weighted_cosine_similarity": (
                        float(cosine.median()) if cosine.notna().any() else float("nan")
                    ),
                    "low_confidence_pair_count": low_confidence_pairs,
                    "fdr_cutoff_used": float(fdr_cutoff),
                    "abs_dpsi_cutoff_used": float(abs_dpsi_cutoff),
                    "top_events_per_comparison": int(top_events),
                    "selection_mode": str(selection_mode),
                    "interpretation": interpretation,
                    "scientific_scope": (
                        "Descriptive sample-PSI pattern, replicate QC and comparison-similarity analysis; "
                        "not a causal regulator model."
                    ),
                }
            ]
        )

    def _native_jutils_cache_schema_matches(self, base: Path) -> bool:
        parameters_path = base / "parameters.json"
        if not parameters_path.exists():
            return False
        try:
            payload = json.loads(parameters_path.read_text(encoding="utf-8"))
        except (OSError, ValueError, TypeError):
            return False
        return str(payload.get("analysis_schema_version", "")) == self.NATIVE_PATTERN_SCHEMA_VERSION

    def _native_jutils_gene_filter_map(self, gene_text: str, comparisons: list) -> dict[str, str]:
        tokens: dict[str, str] = {}
        for token in re.split(r"[\s,;]+", gene_text or ""):
            value = token.strip()
            if not value:
                continue
            tokens.setdefault(value.casefold(), value)
        if tokens:
            return tokens
        selected: dict[str, str] = {}
        for comparison in comparisons:
            selection, _meta = self.downstream_gene_selection(
                comparison.comparison_id,
                consumer_module="5.3.5",
            )
            values = selection.get(
                "gene_symbol",
                selection.get("gene_id", pd.Series(dtype="object")),
            ).dropna().astype(str).tolist()
            for token in values:
                value = str(token or "").strip()
                if value:
                    selected.setdefault(value.casefold(), value)
        return selected

    def _native_jutils_gene_filter(self, gene_text: str, comparisons: list) -> set[str]:
        return set(self._native_jutils_gene_filter_map(gene_text, comparisons))

    def _build_native_jutils_tables(
        self,
        comparisons: list,
        *,
        fdr_cutoff: float,
        abs_dpsi_cutoff: float,
        top_events: int,
        selection_mode: str,
        gene_filter: set[str] | None,
        gene_filter_display: dict[str, str] | None = None,
    ) -> dict[str, pd.DataFrame]:
        event_rows: list[dict[str, object]] = []
        psi_rows: list[dict[str, object]] = []
        annotation_lookup = self._annotation_gene_symbol_lookup()
        sample_ids_by_group = self._native_jutils_sample_ids_by_group()
        mode = self.current_project.selected_rmats_mode if self.current_project else "JC"
        for comparison in comparisons:
            group1_label, group2_label = self._comparison_splicing_direction_labels(comparison)
            display_name = comparison.display_resolved_name
            for event_type in EVENT_TYPES:
                source_path = self._resolve_rmats_event_table_path(comparison, event_type, mode or "JC")
                if source_path is None or not source_path.exists():
                    continue
                try:
                    frame = pd.read_csv(source_path, sep="\t", low_memory=False)
                except Exception:
                    continue
                if frame.empty:
                    continue
                gene_ids = self._normalized_gene_id_series(frame)
                if "geneSymbol" in frame.columns:
                    gene_symbols = frame["geneSymbol"].astype("object")
                elif "gene_symbol" in frame.columns:
                    gene_symbols = frame["gene_symbol"].astype("object")
                else:
                    gene_symbols = gene_ids.map(annotation_lookup)
                gene_symbols = gene_symbols.fillna(gene_ids)
                if gene_filter is not None:
                    symbol_keys = gene_symbols.astype(str).str.casefold()
                    gene_id_keys = gene_ids.astype(str).str.casefold()
                    keep_mask = symbol_keys.isin(gene_filter) | gene_id_keys.isin(gene_filter)
                    if not bool(keep_mask.any()):
                        continue
                    frame = frame.loc[keep_mask].copy()
                    gene_ids = gene_ids.loc[frame.index]
                    gene_symbols = gene_symbols.loc[frame.index]
                fdr_values = pd.to_numeric(frame["FDR"] if "FDR" in frame.columns else frame.get("fdr"), errors="coerce")
                dpsi_values = pd.to_numeric(
                    frame["IncLevelDifference"] if "IncLevelDifference" in frame.columns else frame.get("dPSI"),
                    errors="coerce",
                )
                if comparison.reverse_splicing:
                    dpsi_values = -dpsi_values
                supervised = str(selection_mode).lower() != "unsupervised"
                if supervised:
                    significant_mask = fdr_values.notna() & dpsi_values.notna() & (fdr_values <= fdr_cutoff) & (dpsi_values.abs() >= abs_dpsi_cutoff)
                    if not bool(significant_mask.any()):
                        continue
                    frame = frame.loc[significant_mask].copy()
                    gene_ids = gene_ids.loc[frame.index]
                    gene_symbols = gene_symbols.loc[frame.index]
                    fdr_values = fdr_values.loc[frame.index]
                    dpsi_values = dpsi_values.loc[frame.index]
                event_ids = frame.apply(lambda row: event_uid_from_row(row, event_type), axis=1)
                for idx, row in frame.iterrows():
                    gene_id = str(gene_ids.loc[idx] or "").strip()
                    gene_symbol = str(gene_symbols.loc[idx] or gene_id).strip()
                    if gene_filter is not None:
                        if gene_symbol.casefold() not in gene_filter and gene_id.casefold() not in gene_filter:
                            continue
                    fdr = fdr_values.loc[idx]
                    dpsi = dpsi_values.loc[idx]
                    significant = bool(pd.notna(fdr) and pd.notna(dpsi) and float(fdr) <= fdr_cutoff and abs(float(dpsi)) >= abs_dpsi_cutoff)
                    if supervised and not significant:
                        continue
                    event_id = str(event_ids.loc[idx])
                    matching_key = f"{gene_symbol}|{event_type}|{event_id}"
                    psi_1 = self._parse_float_list(row.get("IncLevel1"))
                    psi_2 = self._parse_float_list(row.get("IncLevel2"))
                    ijc_1 = self._parse_float_list(row.get("IJC_SAMPLE_1"))
                    ijc_2 = self._parse_float_list(row.get("IJC_SAMPLE_2"))
                    sjc_1 = self._parse_float_list(row.get("SJC_SAMPLE_1"))
                    sjc_2 = self._parse_float_list(row.get("SJC_SAMPLE_2"))
                    if comparison.reverse_splicing:
                        psi_1, psi_2 = psi_2, psi_1
                        ijc_1, ijc_2 = ijc_2, ijc_1
                        sjc_1, sjc_2 = sjc_2, sjc_1
                    psi_1_mean = float(np.nanmean(psi_1)) if psi_1 else float("nan")
                    psi_2_mean = float(np.nanmean(psi_2)) if psi_2 else float("nan")
                    event_rows.append(
                        {
                            "comparison_id": comparison.comparison_id,
                            "comparison_display_name": display_name,
                            "gene_symbol": gene_symbol,
                            "gene_id": gene_id,
                            "event_type": event_type,
                            "event_id": event_id,
                            "event_matching_key": matching_key,
                            "FDR": fdr,
                            "dPSI": dpsi,
                            "abs_dPSI": abs(float(dpsi)) if pd.notna(dpsi) else float("nan"),
                            "significant": significant,
                            "direction": "positive" if pd.notna(dpsi) and float(dpsi) > 0 else "negative" if pd.notna(dpsi) and float(dpsi) < 0 else "zero",
                            "psi_group1_mean": psi_1_mean,
                            "psi_group2_mean": psi_2_mean,
                            "group1_label": group1_label,
                            "group2_label": group2_label,
                            "coordinates": self._rmats_event_coordinates(row, event_type),
                            "splicing_direction_flipped": bool(comparison.reverse_splicing),
                            "fdr_cutoff_used": fdr_cutoff,
                            "abs_dpsi_cutoff_used": abs_dpsi_cutoff,
                            "native_jutils_selection_mode": selection_mode,
                        }
                    )
                    max_len = max(len(psi_1), len(psi_2), len(ijc_1), len(ijc_2), len(sjc_1), len(sjc_2))
                    for sample_idx in range(max_len):
                        if sample_idx < len(psi_1):
                            psi_rows.append(
                                {
                                    "comparison_id": comparison.comparison_id,
                                    "comparison_display_name": display_name,
                                    "event_matching_key": matching_key,
                                    "event_id": event_id,
                                    "gene_symbol": gene_symbol,
                                    "event_type": event_type,
                                    "sample_id": f"{comparison.comparison_id}:{group1_label}:rep{sample_idx + 1}",
                                    "biological_sample_id": self._native_jutils_biological_sample_id(
                                        sample_ids_by_group,
                                        group1_label,
                                        sample_idx,
                                    ),
                                    "sample_group": group1_label,
                                    "group_side": "canonical_group1",
                                    "PSI": psi_1[sample_idx],
                                    "IJC": ijc_1[sample_idx] if sample_idx < len(ijc_1) else pd.NA,
                                    "SJC": sjc_1[sample_idx] if sample_idx < len(sjc_1) else pd.NA,
                                }
                            )
                        if sample_idx < len(psi_2):
                            psi_rows.append(
                                {
                                    "comparison_id": comparison.comparison_id,
                                    "comparison_display_name": display_name,
                                    "event_matching_key": matching_key,
                                    "event_id": event_id,
                                    "gene_symbol": gene_symbol,
                                    "event_type": event_type,
                                    "sample_id": f"{comparison.comparison_id}:{group2_label}:rep{sample_idx + 1}",
                                    "biological_sample_id": self._native_jutils_biological_sample_id(
                                        sample_ids_by_group,
                                        group2_label,
                                        sample_idx,
                                    ),
                                    "sample_group": group2_label,
                                    "group_side": "canonical_group2",
                                    "PSI": psi_2[sample_idx],
                                    "IJC": ijc_2[sample_idx] if sample_idx < len(ijc_2) else pd.NA,
                                    "SJC": sjc_2[sample_idx] if sample_idx < len(sjc_2) else pd.NA,
                                }
                            )

        events = pd.DataFrame(event_rows)
        psi_long = pd.DataFrame(psi_rows)
        if events.empty:
            tables = self._empty_native_jutils_tables()
            if gene_filter:
                summary, matrix = self._native_jutils_candidate_gene_summary(
                    events,
                    gene_filter_display or {key: key for key in sorted(gene_filter)},
                    comparisons,
                    fdr_cutoff=fdr_cutoff,
                    abs_dpsi_cutoff=abs_dpsi_cutoff,
                )
                tables["candidate_gene_as_summary"] = summary
                tables["candidate_gene_as_matrix"] = matrix
            return tables
        events = events.sort_values(["significant", "FDR", "abs_dPSI", "event_type", "gene_symbol"], ascending=[False, True, False, True, True], na_position="last")
        top_n = max(int(top_events), 1)
        if str(selection_mode).lower() == "unsupervised" and not psi_long.empty:
            scores = (
                psi_long.assign(PSI_numeric=pd.to_numeric(psi_long["PSI"], errors="coerce"))
                .groupby(["comparison_id", "event_matching_key"], dropna=False)["PSI_numeric"]
                .apply(lambda values: float((values - values.mean()).abs().mean()) if values.notna().any() else float("nan"))
                .reset_index(name="native_jutils_selection_score")
            )
            events = events.merge(scores, on=["comparison_id", "event_matching_key"], how="left")
            events = events.sort_values(
                ["native_jutils_selection_score", "abs_dPSI", "FDR", "event_type", "gene_symbol"],
                ascending=[False, False, True, True, True],
                na_position="last",
            )
        else:
            events["native_jutils_selection_score"] = pd.to_numeric(events["abs_dPSI"], errors="coerce")
        events = events.groupby("comparison_id", sort=False, group_keys=False).head(top_n).copy()
        keep_pairs = set(zip(events["comparison_id"].astype(str), events["event_matching_key"].astype(str)))
        if not psi_long.empty:
            psi_long = psi_long.loc[
                [
                    (str(row.comparison_id), str(row.event_matching_key)) in keep_pairs
                    for row in psi_long[["comparison_id", "event_matching_key"]].itertuples(index=False)
                ]
            ].copy()
        dpsi_matrix = self._native_jutils_dpsi_matrix(events)
        psi_matrix = self._native_jutils_psi_matrix(psi_long)
        psi_zscore_matrix = row_zscore_matrix(psi_matrix)
        similarity = self._native_jutils_similarity(events, comparisons)
        replicate_qc = self._native_jutils_replicate_qc(events)
        sample_pca = self._native_jutils_sample_pca(psi_matrix, psi_long)
        per_comparison_sample_pca = self._native_jutils_per_comparison_sample_pca(psi_long)
        sample_similarity, sample_similarity_matrix, sample_overlap_matrix = sample_similarity_tables(psi_matrix)
        replicate_separation = replicate_group_separation(psi_long, sample_similarity)
        pattern_clusters, pattern_cluster_summary = psi_pattern_clusters(psi_matrix, events)
        tables = {
            "unified_as_event_table": events.reset_index(drop=True),
            "sample_level_psi_long": psi_long.reset_index(drop=True),
            "dpsi_matrix": dpsi_matrix,
            "sample_level_psi_matrix": psi_matrix,
            "sample_level_psi_zscore_matrix": psi_zscore_matrix,
            "pairwise_as_similarity": similarity,
            "sample_psi_similarity": sample_similarity,
            "sample_psi_spearman_matrix": sample_similarity_matrix,
            "sample_psi_shared_event_count_matrix": sample_overlap_matrix,
            "sample_psi_pca": sample_pca,
            "per_comparison_sample_psi_pca": per_comparison_sample_pca,
            "replicate_psi_qc": replicate_qc,
            "replicate_group_separation": replicate_separation,
            "psi_pattern_clusters": pattern_clusters,
            "psi_pattern_cluster_summary": pattern_cluster_summary,
        }
        if gene_filter:
            summary, matrix = self._native_jutils_candidate_gene_summary(
                events,
                gene_filter_display or {key: key for key in sorted(gene_filter)},
                comparisons,
                fdr_cutoff=fdr_cutoff,
                abs_dpsi_cutoff=abs_dpsi_cutoff,
            )
            tables["candidate_gene_as_summary"] = summary
            tables["candidate_gene_as_matrix"] = matrix
        return tables

    def _empty_native_jutils_tables(self) -> dict[str, pd.DataFrame]:
        return {
            "unified_as_event_table": pd.DataFrame(),
            "sample_level_psi_long": pd.DataFrame(),
            "dpsi_matrix": pd.DataFrame(),
            "sample_level_psi_matrix": pd.DataFrame(),
            "sample_level_psi_zscore_matrix": pd.DataFrame(),
            "pairwise_as_similarity": pd.DataFrame(),
            "sample_psi_similarity": pd.DataFrame(),
            "sample_psi_spearman_matrix": pd.DataFrame(),
            "sample_psi_shared_event_count_matrix": pd.DataFrame(),
            "sample_psi_pca": pd.DataFrame(),
            "per_comparison_sample_psi_pca": pd.DataFrame(),
            "replicate_psi_qc": pd.DataFrame(),
            "replicate_group_separation": pd.DataFrame(),
            "psi_pattern_clusters": pd.DataFrame(),
            "psi_pattern_cluster_summary": pd.DataFrame(),
            "candidate_gene_as_summary": pd.DataFrame(),
            "candidate_gene_as_matrix": pd.DataFrame(),
        }

    def _native_jutils_dpsi_matrix(self, events: pd.DataFrame) -> pd.DataFrame:
        matrix = events.pivot_table(index="event_matching_key", columns="comparison_display_name", values="dPSI", aggfunc="first")
        matrix = matrix.reindex(index=events.drop_duplicates("event_matching_key")["event_matching_key"].tolist())
        return matrix.reset_index()

    def _native_jutils_psi_matrix(self, psi_long: pd.DataFrame) -> pd.DataFrame:
        if psi_long.empty:
            return pd.DataFrame()
        sample_column = "biological_sample_id" if "biological_sample_id" in psi_long.columns else "sample_id"
        matrix = psi_long.pivot_table(index="event_matching_key", columns=sample_column, values="PSI", aggfunc="mean")
        return matrix.reset_index()

    def _native_jutils_sample_ids_by_group(self) -> dict[str, list[str]]:
        if self.current_project is None:
            return {}
        grouped: dict[str, list[str]] = {}
        for sample in self.current_project.isoform_samples:
            if not sample.include:
                continue
            group = str(sample.sample_group or sample.condition or "").strip()
            sample_id = str(sample.sample_id or "").strip()
            if not group or not sample_id:
                continue
            grouped.setdefault(group.casefold(), []).append(sample_id)

        def natural_key(value: str) -> list[object]:
            return [int(token) if token.isdigit() else token.casefold() for token in re.split(r"(\d+)", value)]

        return {key: sorted(set(values), key=natural_key) for key, values in grouped.items()}

    @staticmethod
    def _native_jutils_biological_sample_id(
        sample_ids_by_group: dict[str, list[str]],
        group_label: str,
        sample_index: int,
    ) -> str:
        candidates = sample_ids_by_group.get(str(group_label).strip().casefold(), [])
        if sample_index < len(candidates):
            return candidates[sample_index]
        return f"{group_label}:rep{sample_index + 1}"

    def _native_jutils_candidate_gene_summary(
        self,
        events: pd.DataFrame,
        gene_filter_display: dict[str, str],
        comparisons: list,
        *,
        fdr_cutoff: float,
        abs_dpsi_cutoff: float,
    ) -> tuple[pd.DataFrame, pd.DataFrame]:
        rows: list[dict[str, object]] = []
        event_lookup = events.copy()
        if not event_lookup.empty:
            event_lookup["_gene_symbol_key"] = event_lookup["gene_symbol"].astype(str).str.casefold()
            event_lookup["_gene_id_key"] = event_lookup["gene_id"].astype(str).str.casefold()
            event_lookup["FDR_numeric"] = pd.to_numeric(event_lookup["FDR"], errors="coerce")
            event_lookup["dPSI_numeric"] = pd.to_numeric(event_lookup["dPSI"], errors="coerce")
            event_lookup["abs_dPSI_numeric"] = event_lookup["dPSI_numeric"].abs()
        for query_key, display_gene in gene_filter_display.items():
            for comparison in comparisons:
                display_name = comparison.display_resolved_name
                subset = pd.DataFrame()
                if not event_lookup.empty:
                    subset = event_lookup.loc[
                        (event_lookup["comparison_id"].astype(str) == comparison.comparison_id)
                        & (
                            (event_lookup["_gene_symbol_key"] == query_key)
                            | (event_lookup["_gene_id_key"] == query_key)
                        )
                    ].copy()
                significant = subset.loc[subset["significant"].fillna(False)].copy() if not subset.empty else pd.DataFrame()
                if significant.empty:
                    rows.append(
                        {
                            "gene_symbol": display_gene,
                            "gene_query": display_gene,
                            "comparison_id": comparison.comparison_id,
                            "comparison_display_name": display_name,
                            "n_significant_as_events": 0,
                            "event_types": "",
                            "dominant_event_type": "",
                            "dominant_event_id": "",
                            "best_FDR": pd.NA,
                            "best_dPSI": pd.NA,
                            "best_abs_dPSI": pd.NA,
                            "direction": "no_significant_AS",
                            "plotted_in_event_heatmap": False,
                            "selection_note": "selected gene retained; no significant AS event under active cutoffs",
                            "fdr_cutoff_used": fdr_cutoff,
                            "abs_dpsi_cutoff_used": abs_dpsi_cutoff,
                        }
                    )
                    continue
                significant = significant.sort_values(
                    ["abs_dPSI_numeric", "FDR_numeric", "event_type"],
                    ascending=[False, True, True],
                    na_position="last",
                )
                dominant = significant.iloc[0]
                event_types = ";".join(
                    significant["event_type"].dropna().astype(str).drop_duplicates().tolist()
                )
                rows.append(
                    {
                        "gene_symbol": str(dominant.get("gene_symbol") or display_gene),
                        "gene_query": display_gene,
                        "comparison_id": comparison.comparison_id,
                        "comparison_display_name": display_name,
                        "n_significant_as_events": int(len(significant)),
                        "event_types": event_types,
                        "dominant_event_type": dominant.get("event_type", ""),
                        "dominant_event_id": dominant.get("event_id", ""),
                        "best_FDR": dominant.get("FDR", pd.NA),
                        "best_dPSI": dominant.get("dPSI", pd.NA),
                        "best_abs_dPSI": dominant.get("abs_dPSI", pd.NA),
                        "direction": dominant.get("direction", ""),
                        "plotted_in_event_heatmap": True,
                        "selection_note": "dominant event selected by max abs(dPSI) among significant events",
                        "fdr_cutoff_used": fdr_cutoff,
                        "abs_dpsi_cutoff_used": abs_dpsi_cutoff,
                    }
                )
        summary = pd.DataFrame(rows)
        if summary.empty:
            return summary, pd.DataFrame()
        matrix = summary.pivot_table(
            index="gene_symbol",
            columns="comparison_display_name",
            values="best_dPSI",
            aggfunc="first",
            dropna=False,
        )
        ordered_genes = []
        seen: set[str] = set()
        for gene in summary["gene_symbol"].astype(str).tolist():
            if gene not in seen:
                seen.add(gene)
                ordered_genes.append(gene)
        matrix = matrix.reindex(index=ordered_genes)
        return summary, matrix.reset_index()

    def _native_jutils_similarity(self, events: pd.DataFrame, comparisons: list) -> pd.DataFrame:
        rows: list[dict[str, object]] = []
        by_comparison = {
            comparison.display_resolved_name: events.loc[
                (events["comparison_id"].astype(str) == comparison.comparison_id) & events["significant"].fillna(False)
            ].copy()
            for comparison in comparisons
        }
        for comp_a, comp_b in combinations(comparisons, 2):
            name_a = comp_a.display_resolved_name
            name_b = comp_b.display_resolved_name
            frame_a = by_comparison.get(name_a, pd.DataFrame())
            frame_b = by_comparison.get(name_b, pd.DataFrame())
            map_a = frame_a.drop_duplicates("event_matching_key").set_index("event_matching_key") if not frame_a.empty else pd.DataFrame()
            map_b = frame_b.drop_duplicates("event_matching_key").set_index("event_matching_key") if not frame_b.empty else pd.DataFrame()
            keys_a = set(map_a.index.astype(str)) if not map_a.empty else set()
            keys_b = set(map_b.index.astype(str)) if not map_b.empty else set()
            shared = keys_a & keys_b
            union = keys_a | keys_b
            same_direction = 0
            opposite_direction = 0
            dpsi_a: list[float] = []
            dpsi_b: list[float] = []
            for key in shared:
                value_a = float(map_a.loc[key, "dPSI"])
                value_b = float(map_b.loc[key, "dPSI"])
                dpsi_a.append(value_a)
                dpsi_b.append(value_b)
                if np.sign(value_a) == np.sign(value_b):
                    same_direction += 1
                else:
                    opposite_direction += 1
            signed_a = {f"{key}|{int(np.sign(float(map_a.loc[key, 'dPSI'])))}" for key in keys_a if pd.notna(map_a.loc[key, "dPSI"])}
            signed_b = {f"{key}|{int(np.sign(float(map_b.loc[key, 'dPSI'])))}" for key in keys_b if pd.notna(map_b.loc[key, "dPSI"])}
            weighted_cosine = float("nan")
            if union:
                vector_a = []
                vector_b = []
                for key in sorted(union):
                    vector_a.append(float(map_a.loc[key, "dPSI"]) if key in keys_a else 0.0)
                    vector_b.append(float(map_b.loc[key, "dPSI"]) if key in keys_b else 0.0)
                a_arr = np.asarray(vector_a, dtype=float)
                b_arr = np.asarray(vector_b, dtype=float)
                denom = float(np.linalg.norm(a_arr) * np.linalg.norm(b_arr))
                weighted_cosine = float(np.dot(a_arr, b_arr) / denom) if denom > 0 else float("nan")
            rows.append(
                {
                    "comparison_A": comp_a.comparison_id,
                    "comparison_B": comp_b.comparison_id,
                    "comparison_A_display_name": name_a,
                    "comparison_B_display_name": name_b,
                    "n_significant_A": len(keys_a),
                    "n_significant_B": len(keys_b),
                    "shared_significant_events": len(shared),
                    "union_significant_events": len(union),
                    "shared_event_support_fraction": (
                        len(shared) / min(len(keys_a), len(keys_b))
                        if min(len(keys_a), len(keys_b)) > 0
                        else float("nan")
                    ),
                    "jaccard_overlap": len(shared) / len(union) if union else float("nan"),
                    "signed_jaccard_overlap": len(signed_a & signed_b) / len(signed_a | signed_b) if (signed_a | signed_b) else float("nan"),
                    "direction_similarity": same_direction / len(shared) if shared else float("nan"),
                    "opposite_direction_fraction": opposite_direction / len(shared) if shared else float("nan"),
                    "spearman_dPSI_correlation": pd.Series(dpsi_a).corr(pd.Series(dpsi_b), method="spearman") if len(dpsi_a) >= 2 else float("nan"),
                    "weighted_cosine_similarity": weighted_cosine,
                    "similarity_evidence_weight": min(1.0, len(shared) / 20.0),
                    "confidence_weighted_cosine_similarity": (
                        weighted_cosine * min(1.0, len(shared) / 20.0)
                        if pd.notna(weighted_cosine)
                        else float("nan")
                    ),
                    "similarity_confidence": (
                        "high_shared_support"
                        if len(shared) >= 50
                        else "moderate_shared_support"
                        if len(shared) >= 10
                        else "low_shared_support_interpret_cautiously"
                        if len(shared) > 0
                        else "no_shared_events"
                    ),
                    "similarity_interpretation": (
                        f"{len(shared)} shared significant event(s) out of {len(union)} union event(s). "
                        "Interpret weighted_cosine_similarity together with shared_event_support_fraction; "
                        "confidence_weighted_cosine_similarity downweights pairs with few shared events."
                    ),
                    "specificity_A": (len(keys_a - keys_b) / len(keys_a)) if keys_a else float("nan"),
                    "specificity_B": (len(keys_b - keys_a) / len(keys_b)) if keys_b else float("nan"),
                }
            )
        return pd.DataFrame(rows)

    def _native_jutils_sample_pca(self, psi_matrix: pd.DataFrame, psi_long: pd.DataFrame) -> pd.DataFrame:
        if psi_matrix.empty or psi_matrix.shape[1] <= 2:
            return pd.DataFrame()
        numeric = psi_matrix.drop(columns=[psi_matrix.columns[0]]).apply(pd.to_numeric, errors="coerce")
        scores, explained = self._pca_scores(numeric.transpose())
        if scores.empty:
            return pd.DataFrame()
        output = scores.rename_axis("sample_id").reset_index()
        if not psi_long.empty and "biological_sample_id" in psi_long.columns:
            metadata = (
                psi_long.groupby("biological_sample_id", dropna=False)
                .agg(
                    sample_group=("sample_group", lambda values: "; ".join(dict.fromkeys(values.astype(str)))),
                    comparison_memberships=("comparison_id", lambda values: "; ".join(dict.fromkeys(values.astype(str)))),
                )
            )
            output["sample_group"] = output["sample_id"].map(metadata["sample_group"])
            output["comparison_memberships"] = output["sample_id"].map(metadata["comparison_memberships"])
        else:
            output["sample_group"] = output["sample_id"].astype(str).str.split(":").str[0]
            output["comparison_memberships"] = ""
        output["PC1_variance_explained"] = explained[0] if len(explained) > 0 else float("nan")
        output["PC2_variance_explained"] = explained[1] if len(explained) > 1 else float("nan")
        return output

    def _native_jutils_per_comparison_sample_pca(self, psi_long: pd.DataFrame) -> pd.DataFrame:
        if psi_long.empty or not {"comparison_id", "event_matching_key", "sample_id", "PSI"}.issubset(psi_long.columns):
            return pd.DataFrame()
        rows: list[pd.DataFrame] = []
        for comparison_id, subset in psi_long.groupby("comparison_id", sort=False):
            pivot = subset.pivot_table(
                index="event_matching_key",
                columns="sample_id",
                values="PSI",
                aggfunc="mean",
            )
            if pivot.shape[0] < 2 or pivot.shape[1] < 2:
                continue
            scores, explained = self._pca_scores(pivot.transpose())
            if scores.empty:
                continue
            output = scores.rename_axis("sample_id").reset_index()
            output["comparison_id"] = str(comparison_id)
            output["comparison_display_name"] = (
                str(subset["comparison_display_name"].dropna().iloc[0])
                if "comparison_display_name" in subset.columns and subset["comparison_display_name"].notna().any()
                else str(comparison_id)
            )
            output["sample_group"] = output["sample_id"].astype(str).str.split(":").str[1]
            output["PC1_variance_explained"] = explained[0] if len(explained) > 0 else float("nan")
            output["PC2_variance_explained"] = explained[1] if len(explained) > 1 else float("nan")
            output["pca_scope"] = "per_comparison_sample_PSI"
            rows.append(output)
        return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()

    def _pca_scores(self, matrix: pd.DataFrame) -> tuple[pd.DataFrame, list[float]]:
        if matrix.empty:
            return pd.DataFrame(), []
        numeric = matrix.copy()
        numeric = numeric.apply(pd.to_numeric, errors="coerce")
        numeric = numeric.loc[:, numeric.notna().any(axis=0)]
        numeric = numeric.fillna(numeric.mean(axis=0)).fillna(0.0)
        if numeric.shape[0] < 2 or numeric.shape[1] < 2:
            return pd.DataFrame(), []
        centered = numeric - numeric.mean(axis=0)
        std = centered.std(axis=0).replace(0, 1.0).fillna(1.0)
        scaled = centered / std
        try:
            u, s, _vt = np.linalg.svd(scaled.to_numpy(dtype=float), full_matrices=False)
        except np.linalg.LinAlgError:
            return pd.DataFrame(), []
        scores = u[:, :2] * s[:2]
        explained = (s**2) / max(float(np.sum(s**2)), 1e-12)
        output = pd.DataFrame(scores, index=numeric.index, columns=["PC1", "PC2"][: scores.shape[1]])
        if "PC2" not in output.columns:
            output["PC2"] = 0.0
        return output, explained[:2].tolist()

    def _native_jutils_replicate_qc(self, events: pd.DataFrame) -> pd.DataFrame:
        if events.empty:
            return pd.DataFrame()
        working = events.copy()
        working["psi_delta_from_group_means"] = pd.to_numeric(working["psi_group1_mean"], errors="coerce") - pd.to_numeric(working["psi_group2_mean"], errors="coerce")
        working["psi_group_mean_disagrees_with_dPSI"] = np.sign(working["psi_delta_from_group_means"].fillna(0)) != np.sign(pd.to_numeric(working["dPSI"], errors="coerce").fillna(0))
        return working[
            [
                "comparison_id",
                "comparison_display_name",
                "event_matching_key",
                "gene_symbol",
                "event_type",
                "FDR",
                "dPSI",
                "psi_group1_mean",
                "psi_group2_mean",
                "psi_delta_from_group_means",
                "psi_group_mean_disagrees_with_dPSI",
            ]
        ].copy()

    def _write_native_jutils_outputs(self, base: Path, tables: dict[str, pd.DataFrame], parameters: dict[str, object]) -> list[Path]:
        if base.exists():
            shutil.rmtree(base)
        base.mkdir(parents=True, exist_ok=True)
        plots_dir = base / "plots"
        plots_dir.mkdir(parents=True, exist_ok=True)
        generated: list[Path] = []
        for name, frame in tables.items():
            path = base / f"{name}.tsv"
            frame.to_csv(path, sep="\t", index=False)
            generated.append(path)
        plot_files = self._plot_native_jutils_tables(tables, plots_dir)
        benchmark_files = self._write_native_vs_legacy_jutils_benchmark(base, tables)
        parameters = dict(parameters)
        parameters["generated_plots"] = [str(path) for path in plot_files]
        parameters["generated_tables"] = [str(base / f"{name}.tsv") for name in tables]
        parameters["generated_benchmarks"] = [str(path) for path in benchmark_files]
        self._write_json(base / "parameters.json", parameters)
        generated.append(base / "parameters.json")
        generated.extend(plot_files)
        generated.extend(benchmark_files)
        return generated

    def _write_native_vs_legacy_jutils_benchmark(self, base: Path, tables: dict[str, pd.DataFrame]) -> list[Path]:
        if self.current_project is None:
            return []
        return self.jutils_benchmark.write(base, self.require_output_root(), tables)

    def _plot_native_jutils_tables(self, tables: dict[str, pd.DataFrame], plots_dir: Path) -> list[Path]:
        generated: list[Path] = []
        candidate_summary = tables.get("candidate_gene_as_summary", pd.DataFrame())
        if not candidate_summary.empty:
            generated.append(
                self._plot_candidate_gene_as_heatmap(
                    candidate_summary,
                    plots_dir / "native_candidate_gene_as_summary_heatmap.png",
                )
            )
        generated.extend(self._plot_native_per_comparison_psi_heatmaps(tables, plots_dir / "per_comparison"))
        generated.extend(self._plot_native_per_comparison_psi_pcas(tables, plots_dir / "per_comparison"))
        dpsi_matrix = tables.get("dpsi_matrix", pd.DataFrame())
        if not dpsi_matrix.empty:
            generated.append(
                self._plot_matrix_heatmap(
                    dpsi_matrix,
                    "Native AS dPSI pattern heatmap",
                    plots_dir / "native_as_dpsi_pattern_heatmap.png",
                    center_zero=True,
                    cluster_rows=True,
                )
            )
        psi_matrix = tables.get("sample_level_psi_matrix", pd.DataFrame())
        if not psi_matrix.empty:
            generated.append(
                self._plot_matrix_heatmap(
                    psi_matrix,
                    "Native sample-level PSI heatmap",
                    plots_dir / "native_sample_level_psi_heatmap.png",
                    center_zero=False,
                    cluster_rows=True,
                )
            )
        psi_zscore_matrix = tables.get("sample_level_psi_zscore_matrix", pd.DataFrame())
        if not psi_zscore_matrix.empty:
            generated.append(
                self._plot_matrix_heatmap(
                    psi_zscore_matrix,
                    "Sample-level PSI heatmap (row Z-score)",
                    plots_dir / "native_sample_level_psi_zscore_heatmap.png",
                    center_zero=True,
                    cluster_rows=True,
                    cluster_columns=True,
                )
            )
        similarity = tables.get("pairwise_as_similarity", pd.DataFrame())
        if not similarity.empty:
            generated.append(self._plot_similarity_heatmap(similarity, "weighted_cosine_similarity", "AS weighted cosine similarity", plots_dir / "native_as_similarity_heatmap.png"))
            generated.append(self._plot_similarity_heatmap(similarity, "direction_similarity", "AS direction agreement", plots_dir / "native_as_direction_similarity_heatmap.png"))
        sample_similarity = tables.get("sample_psi_spearman_matrix", pd.DataFrame())
        sample_overlap = tables.get("sample_psi_shared_event_count_matrix", pd.DataFrame())
        if not sample_similarity.empty:
            generated.append(
                self._plot_square_matrix_heatmap(
                    sample_similarity,
                    "Sample PSI similarity (Spearman)",
                    plots_dir / "native_sample_psi_similarity_heatmap.png",
                    value_label="Spearman rho",
                    minimum=-1.0,
                    maximum=1.0,
                    cmap="RdBu_r",
                )
            )
        if not sample_overlap.empty:
            maximum = pd.to_numeric(
                sample_overlap.drop(columns=[sample_overlap.columns[0]]).stack(), errors="coerce"
            ).max()
            generated.append(
                self._plot_square_matrix_heatmap(
                    sample_overlap,
                    "Shared measured AS events per sample pair",
                    plots_dir / "native_sample_psi_shared_event_counts.png",
                    value_label="shared events",
                    minimum=0.0,
                    maximum=max(float(maximum) if pd.notna(maximum) else 1.0, 1.0),
                    cmap="YlGnBu",
                    integer_labels=True,
                )
            )
        sample_pca = tables.get("sample_psi_pca", pd.DataFrame())
        if not sample_pca.empty:
            generated.append(
                self._plot_native_pca(
                    sample_pca,
                    "sample_id",
                    "Native sample-level PSI PCA",
                    plots_dir / "native_sample_psi_pca.png",
                    color_column="sample_group",
                )
            )
        pattern_clusters = tables.get("psi_pattern_clusters", pd.DataFrame())
        if not pattern_clusters.empty:
            generated.append(
                self._plot_psi_pattern_cluster_heatmap(
                    pattern_clusters,
                    plots_dir / "native_sample_psi_pattern_clusters.png",
                )
            )
        replicate_separation = tables.get("replicate_group_separation", pd.DataFrame())
        if not replicate_separation.empty:
            generated.append(
                self._plot_replicate_group_separation(
                    replicate_separation,
                    plots_dir / "native_replicate_group_separation.png",
                )
            )
        return [path for path in generated if path.exists()]

    def _plot_square_matrix_heatmap(
        self,
        frame: pd.DataFrame,
        title: str,
        path: Path,
        *,
        value_label: str,
        minimum: float,
        maximum: float,
        cmap: str,
        integer_labels: bool = False,
    ) -> Path:
        if frame.empty or frame.shape[1] <= 1:
            return path
        labels = frame.iloc[:, 0].astype(str).tolist()
        numeric = frame.iloc[:, 1:].apply(pd.to_numeric, errors="coerce")
        numeric.columns = numeric.columns.astype(str)
        numeric = numeric.reindex(columns=labels)
        data = numeric.to_numpy(dtype=float)
        size = max(6.5, min(16, 3.5 + len(labels) * 0.75))
        fig, ax = plt.subplots(figsize=(size, size))
        image = ax.imshow(np.ma.masked_invalid(data), cmap=cmap, vmin=minimum, vmax=maximum)
        ax.set_title(title, fontweight="bold", pad=14)
        ax.set_xticks(range(len(labels)))
        ax.set_xticklabels(labels, rotation=90, ha="center", fontsize=8)
        ax.set_yticks(range(len(labels)))
        ax.set_yticklabels(labels, fontsize=8)
        if len(labels) <= 18:
            for row_index in range(len(labels)):
                for column_index in range(len(labels)):
                    value = data[row_index, column_index]
                    if np.isfinite(value):
                        label = f"{int(round(value))}" if integer_labels else f"{value:.2f}"
                        ax.text(column_index, row_index, label, ha="center", va="center", fontsize=7)
        fig.colorbar(image, ax=ax, fraction=0.04, pad=0.03, label=value_label)
        fig.tight_layout()
        path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(path, dpi=220, bbox_inches="tight")
        plt.close(fig)
        return path

    def _plot_psi_pattern_cluster_heatmap(self, frame: pd.DataFrame, path: Path) -> Path:
        if frame.empty or "psi_pattern_cluster" not in frame.columns:
            return path
        metadata = {
            "event_matching_key",
            "gene_symbol",
            "event_type",
            "event_id",
            "psi_pattern_cluster",
        }
        value_columns = [
            column
            for column in frame.columns
            if column not in metadata and pd.api.types.is_numeric_dtype(frame[column])
        ]
        if not value_columns:
            return path
        working = frame.sort_values(["psi_pattern_cluster", "event_type", "gene_symbol"], na_position="last").head(120)
        data = working[value_columns].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
        fig, ax = plt.subplots(
            figsize=(max(9, min(18, 3 + len(value_columns) * 0.75)), max(6, min(18, 2 + len(working) * 0.16)))
        )
        image = ax.imshow(np.ma.masked_invalid(data), aspect="auto", cmap="RdBu_r", vmin=-2.5, vmax=2.5)
        ax.set_title("Sample-level PSI trajectory clusters", fontweight="bold", pad=14)
        ax.set_xticks(range(len(value_columns)))
        ax.set_xticklabels(value_columns, rotation=90, fontsize=7)
        labels = [
            f"C{row.psi_pattern_cluster} | {str(getattr(row, 'gene_symbol', ''))[:18]} | {getattr(row, 'event_type', '')}"
            for row in working.itertuples(index=False)
        ]
        ax.set_yticks(range(len(labels)))
        ax.set_yticklabels(labels, fontsize=6)
        fig.colorbar(image, ax=ax, fraction=0.025, pad=0.02, label="row Z-score PSI")
        fig.tight_layout()
        path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(path, dpi=220, bbox_inches="tight")
        plt.close(fig)
        return path

    def _plot_replicate_group_separation(self, frame: pd.DataFrame, path: Path) -> Path:
        working = frame.copy()
        if working.empty or "comparison_id" not in working.columns:
            return path
        within = pd.to_numeric(working["median_within_group_spearman"], errors="coerce")
        between = pd.to_numeric(working["median_between_group_spearman"], errors="coerce")
        x = np.arange(len(working))
        fig, ax = plt.subplots(figsize=(max(8, len(working) * 1.35), 5.5))
        width = 0.36
        ax.bar(x - width / 2, within, width, label="within-group replicates", color="#20a39e")
        ax.bar(x + width / 2, between, width, label="between groups", color="#f19c79")
        ax.axhline(0, color="#6d7b86", linewidth=0.8)
        ax.set_ylim(-1, 1)
        ax.set_ylabel("Median Spearman PSI similarity")
        ax.set_title("Replicate/group PSI similarity QC", fontweight="bold")
        ax.set_xticks(x)
        ax.set_xticklabels(working["comparison_id"].astype(str), rotation=35, ha="right")
        ax.legend(frameon=False, loc="upper right")
        fig.tight_layout()
        path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(path, dpi=220, bbox_inches="tight")
        plt.close(fig)
        return path

    def _plot_candidate_gene_as_heatmap(self, summary: pd.DataFrame, path: Path) -> Path:
        required = {"gene_symbol", "comparison_display_name", "best_dPSI", "n_significant_as_events"}
        if summary.empty or not required.issubset(summary.columns):
            return path
        working = summary.copy()
        working["best_dPSI_numeric"] = pd.to_numeric(working["best_dPSI"], errors="coerce")
        working["n_significant_numeric"] = pd.to_numeric(working["n_significant_as_events"], errors="coerce").fillna(0)
        pivot = working.pivot_table(
            index="gene_symbol",
            columns="comparison_display_name",
            values="best_dPSI_numeric",
            aggfunc="first",
            dropna=False,
        )
        counts = working.pivot_table(
            index="gene_symbol",
            columns="comparison_display_name",
            values="n_significant_numeric",
            aggfunc="max",
            fill_value=0,
            dropna=False,
        ).reindex(index=pivot.index, columns=pivot.columns, fill_value=0)
        order = (
            pd.DataFrame(
                {
                    "n_total_significant": counts.sum(axis=1),
                    "max_abs_dPSI": pivot.abs().max(axis=1).fillna(0),
                    "gene": pivot.index.astype(str),
                },
                index=pivot.index,
            )
            .sort_values(["n_total_significant", "max_abs_dPSI", "gene"], ascending=[False, False, True])
            .index
        )
        pivot = pivot.loc[order].head(80)
        counts = counts.loc[pivot.index]
        data = pivot.to_numpy(dtype=float)
        masked = np.ma.masked_invalid(data)
        vmax = float(np.nanmax(np.abs(data))) if np.isfinite(data).any() else 0.2
        vmax = max(vmax, 0.2)
        cmap = plt.get_cmap("RdBu_r").copy()
        cmap.set_bad("#d8d8d8")
        width = max(8, min(18, 2.8 + pivot.shape[1] * 1.4))
        height = max(5.5, min(18, 2.2 + pivot.shape[0] * 0.34))
        fig, ax = plt.subplots(figsize=(width, height))
        image = ax.imshow(masked, aspect="auto", cmap=cmap, vmin=-vmax, vmax=vmax)
        ax.set_title("Candidate gene AS summary heatmap\nselected genes retained; gray = no significant AS under active cutoffs")
        ax.set_xticks(range(pivot.shape[1]))
        ax.set_xticklabels([str(col) for col in pivot.columns], rotation=35, ha="right", fontsize=8)
        ax.set_yticks(range(pivot.shape[0]))
        ax.set_yticklabels([str(item)[:45] for item in pivot.index], fontsize=8)
        for row_idx in range(pivot.shape[0]):
            for col_idx in range(pivot.shape[1]):
                count_value = int(counts.iloc[row_idx, col_idx]) if pd.notna(counts.iloc[row_idx, col_idx]) else 0
                if count_value > 0:
                    ax.text(col_idx, row_idx, f"n={count_value}", ha="center", va="center", fontsize=6, color="black")
        fig.colorbar(image, ax=ax, fraction=0.025, pad=0.02, label="dominant significant AS dPSI")
        fig.tight_layout()
        path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(path, dpi=220)
        plt.close(fig)
        return path

    def _plot_native_per_comparison_psi_pcas(self, tables: dict[str, pd.DataFrame], plots_dir: Path) -> list[Path]:
        frame = tables.get("per_comparison_sample_psi_pca", pd.DataFrame())
        if frame.empty or "comparison_id" not in frame.columns:
            return []
        generated: list[Path] = []
        for comparison_id, subset in frame.groupby("comparison_id", sort=False):
            display_name = (
                str(subset["comparison_display_name"].dropna().iloc[0])
                if "comparison_display_name" in subset.columns and subset["comparison_display_name"].notna().any()
                else str(comparison_id)
            )
            path = plots_dir / f"{self._safe_output_name(display_name)}_sample_psi_pca.png"
            generated.append(
                self._plot_native_pca(
                    subset,
                    "sample_id",
                    f"Native sample PSI PCA: {display_name}",
                    path,
                    color_column="sample_group",
                )
            )
        return generated

    def _plot_native_per_comparison_psi_heatmaps(self, tables: dict[str, pd.DataFrame], plots_dir: Path) -> list[Path]:
        psi_long = tables.get("sample_level_psi_long", pd.DataFrame())
        events = tables.get("unified_as_event_table", pd.DataFrame())
        if psi_long.empty or "comparison_id" not in psi_long.columns:
            return []
        generated: list[Path] = []
        event_meta = pd.DataFrame()
        if not events.empty and {"comparison_id", "event_matching_key"}.issubset(events.columns):
            meta_cols = [
                column
                for column in [
                    "comparison_id",
                    "comparison_display_name",
                    "event_matching_key",
                    "gene_symbol",
                    "event_type",
                    "FDR",
                    "dPSI",
                    "abs_dPSI",
                    "native_jutils_selection_score",
                ]
                if column in events.columns
            ]
            event_meta = events[meta_cols].drop_duplicates(["comparison_id", "event_matching_key"])
        for comparison_id, subset in psi_long.groupby("comparison_id", sort=False):
            pivot = subset.pivot_table(
                index="event_matching_key",
                columns="sample_id",
                values="PSI",
                aggfunc="mean",
            )
            if pivot.empty:
                continue
            if not event_meta.empty:
                meta = event_meta.loc[event_meta["comparison_id"].astype(str) == str(comparison_id)].copy()
                if not meta.empty:
                    sort_columns = [column for column in ["native_jutils_selection_score", "abs_dPSI", "FDR"] if column in meta.columns]
                    sort_ascending = [False if column != "FDR" else True for column in sort_columns]
                    order_source = meta.sort_values(sort_columns, ascending=sort_ascending, na_position="last") if sort_columns else meta
                    order = order_source["event_matching_key"].astype(str).tolist()
                    pivot = pivot.loc[[key for key in order if key in pivot.index]]
                    label_lookup = {
                        str(row.event_matching_key): f"{str(row.gene_symbol)[:18]} | {str(row.event_type)}"
                        for row in meta[["event_matching_key", "gene_symbol", "event_type"]].itertuples(index=False)
                    }
                    pivot.index = [label_lookup.get(str(index), str(index)[:50]) for index in pivot.index]
            display_name = str(subset["comparison_display_name"].dropna().iloc[0]) if "comparison_display_name" in subset.columns and subset["comparison_display_name"].notna().any() else str(comparison_id)
            path = plots_dir / f"{self._safe_output_name(display_name)}_sample_psi_heatmap.png"
            generated.append(
                self._plot_matrix_heatmap(
                    pivot.reset_index().rename(columns={"index": "event"}),
                    f"Native sample PSI heatmap: {display_name}",
                    path,
                    center_zero=False,
                    cluster_rows=True,
                )
            )
            clustered_path = plots_dir / f"{self._safe_output_name(display_name)}_sample_psi_heatmap_clustered.png"
            generated.append(
                self._plot_matrix_heatmap(
                    pivot.reset_index().rename(columns={"index": "event"}),
                    f"Native clustered sample PSI heatmap: {display_name}",
                    clustered_path,
                    center_zero=False,
                    cluster_rows=True,
                    cluster_columns=True,
                )
            )
        return generated

    def _plot_matrix_heatmap(
        self,
        frame: pd.DataFrame,
        title: str,
        path: Path,
        *,
        center_zero: bool,
        cluster_rows: bool = True,
        cluster_columns: bool = False,
    ) -> Path:
        working = frame.copy()
        if working.empty or working.shape[1] <= 1:
            return path
        index_col = working.columns[0]
        working = working.set_index(index_col)
        numeric = working.apply(pd.to_numeric, errors="coerce")
        row_score = numeric.abs().max(axis=1).fillna(0)
        numeric = numeric.loc[row_score.sort_values(ascending=False).index].head(120)
        if cluster_rows and len(numeric) > 2:
            numeric = self._cluster_order_matrix_rows(numeric)
        if cluster_columns and numeric.shape[1] > 2:
            numeric = self._cluster_order_matrix_rows(numeric.T).T
        numeric = numeric.head(80)
        data = numeric.to_numpy(dtype=float)
        masked = np.ma.masked_invalid(data)
        width = max(8, min(18, 2 + numeric.shape[1] * 1.2))
        height = max(5, min(18, 2 + numeric.shape[0] * 0.22))
        fig, ax = plt.subplots(figsize=(width, height))
        if center_zero:
            vmax = float(np.nanmax(np.abs(data))) if np.isfinite(data).any() else 1.0
            vmax = max(vmax, 0.2)
            image = ax.imshow(masked, aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
        else:
            image = ax.imshow(masked, aspect="auto", cmap="magma", vmin=0, vmax=1)
        ax.set_title(title)
        ax.set_xticks(range(numeric.shape[1]))
        ax.set_xticklabels([str(col) for col in numeric.columns], rotation=45, ha="right", fontsize=8)
        y_labels = [str(item)[:70] for item in numeric.index]
        ax.set_yticks(range(numeric.shape[0]))
        ax.set_yticklabels(y_labels, fontsize=6)
        fig.colorbar(image, ax=ax, fraction=0.025, pad=0.02)
        fig.tight_layout()
        path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(path, dpi=220)
        plt.close(fig)
        return path

    def _cluster_order_matrix_rows(self, numeric: pd.DataFrame) -> pd.DataFrame:
        working = numeric.apply(pd.to_numeric, errors="coerce")
        filled = working.fillna(working.mean(axis=0)).fillna(0.0)
        if len(filled) <= 2:
            return working
        try:
            from scipy.cluster.hierarchy import leaves_list, linkage
            from scipy.spatial.distance import pdist

            distances = pdist(filled.to_numpy(dtype=float), metric="cityblock")
            if len(distances) and np.isfinite(distances).all():
                ordered_index = filled.index[leaves_list(linkage(distances, method="weighted"))]
                return working.loc[ordered_index]
        except Exception:
            pass
        # Stable fallback when SciPy is unavailable or the matrix is degenerate.
        order_frame = pd.DataFrame(
            {
                "row": filled.index.astype(str),
                "mean_abs_signal": filled.abs().mean(axis=1).to_numpy(dtype=float),
                "max_abs_signal": filled.abs().max(axis=1).to_numpy(dtype=float),
            },
            index=filled.index,
        )
        ordered_index = order_frame.sort_values(["max_abs_signal", "mean_abs_signal", "row"], ascending=[False, False, True]).index
        return working.loc[ordered_index]

    def _plot_similarity_heatmap(self, frame: pd.DataFrame, metric: str, title: str, path: Path) -> Path:
        labels = sorted(set(frame["comparison_A_display_name"].astype(str)) | set(frame["comparison_B_display_name"].astype(str)))
        matrix = pd.DataFrame(np.nan, index=labels, columns=labels, dtype=float)
        for label in labels:
            matrix.loc[label, label] = 1.0
        for _, row in frame.iterrows():
            left = str(row["comparison_A_display_name"])
            right = str(row["comparison_B_display_name"])
            value = pd.to_numeric(row.get(metric), errors="coerce")
            matrix.loc[left, right] = value
            matrix.loc[right, left] = value
        fig, ax = plt.subplots(figsize=(max(6, len(labels) * 1.2), max(5, len(labels) * 1.0)))
        image = ax.imshow(np.ma.masked_invalid(matrix.to_numpy(dtype=float)), cmap="viridis", vmin=-1 if "cosine" in metric else 0, vmax=1)
        ax.set_title(title)
        ax.set_xticks(range(len(labels)))
        ax.set_xticklabels(labels, rotation=45, ha="right")
        ax.set_yticks(range(len(labels)))
        ax.set_yticklabels(labels)
        for row_idx in range(len(labels)):
            for col_idx in range(len(labels)):
                value = matrix.iloc[row_idx, col_idx]
                if pd.notna(value):
                    ax.text(col_idx, row_idx, f"{float(value):.2f}", ha="center", va="center", color="white" if float(value) < 0.45 else "black", fontsize=8)
        fig.colorbar(image, ax=ax, fraction=0.04, pad=0.03)
        fig.tight_layout()
        path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(path, dpi=220)
        plt.close(fig)
        return path

    def _plot_native_pca(self, frame: pd.DataFrame, label_column: str, title: str, path: Path, *, color_column: str | None = None) -> Path:
        if frame.empty or "PC1" not in frame.columns or "PC2" not in frame.columns:
            return path
        fig, ax = plt.subplots(figsize=(8, 6))
        label_offsets: dict[object, tuple[int, int]] = {}
        if color_column and color_column in frame.columns:
            groups = frame[color_column].astype(str).fillna("")
            unique_groups = sorted(groups.unique().tolist())
            cmap = plt.get_cmap("tab10")
            for idx, group in enumerate(unique_groups):
                subset = frame.loc[groups == group]
                ax.scatter(subset["PC1"], subset["PC2"], s=70, label=group, color=cmap(idx % 10), edgecolor="black", linewidth=0.4)
                ordered_indices = subset.sort_values(["PC2", "PC1"]).index.tolist()
                offsets = [-12, 0, 12, -20, 20]
                for rank, row_index in enumerate(ordered_indices):
                    label_offsets[row_index] = (7, offsets[rank % len(offsets)])
        else:
            ax.scatter(frame["PC1"], frame["PC2"], s=85, color="#4C78A8", edgecolor="black", linewidth=0.5)
        for row_index, row in frame.iterrows():
            label = str(row.get(label_column, ""))[:38]
            offset = label_offsets.get(row_index, (7, 0))
            ax.annotate(
                label,
                (float(row["PC1"]), float(row["PC2"])),
                xytext=offset,
                textcoords="offset points",
                fontsize=8,
                va="center",
                arrowprops={"arrowstyle": "-", "color": "#7b8794", "linewidth": 0.45},
            )
        pc1_var = pd.to_numeric(frame.get("PC1_variance_explained"), errors="coerce").dropna()
        pc2_var = pd.to_numeric(frame.get("PC2_variance_explained"), errors="coerce").dropna()
        x_label = f"PC1 ({float(pc1_var.iloc[0]) * 100:.1f}% var)" if not pc1_var.empty else "PC1"
        y_label = f"PC2 ({float(pc2_var.iloc[0]) * 100:.1f}% var)" if not pc2_var.empty else "PC2"
        ax.set_title(title)
        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)
        ax.axhline(0, color="#999999", linewidth=0.8)
        ax.axvline(0, color="#999999", linewidth=0.8)
        if color_column and color_column in frame.columns:
            ax.legend(title=color_column, bbox_to_anchor=(1.02, 1), loc="upper left")
        fig.tight_layout()
        path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(path, dpi=220)
        plt.close(fig)
        return path

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

    def rbp_support_output_dir(self, comparison_id: str | None) -> Path | None:
        if self.current_project is None or self.current_project.output_root is None or not comparison_id:
            return None
        return self.current_project.output_root / "07_rbp_support" / self._safe_output_name(comparison_id)

    def rbp_support_file_map(self, comparison_id: str | None) -> dict[str, Path]:
        base = self.rbp_support_output_dir(comparison_id)
        if base is None:
            return {}
        return {
            "summary": base / "rbp_enrichment_summary.tsv",
            "matched_genes": base / "rbp_candidate_gene_matches.tsv",
            "parameters": base / "parameters.json",
        }

    def rbp_cross_support_output_dir(self) -> Path | None:
        if self.current_project is None or self.current_project.output_root is None:
            return None
        return self.current_project.output_root / "07_rbp_support" / "cross_comparison"

    def rbp_cross_support_plot_dir(self) -> Path | None:
        base = self.rbp_cross_support_output_dir()
        return base / "plots" if base is not None else None

    def rbp_cross_support_file_map(self) -> dict[str, Path]:
        base = self.rbp_cross_support_output_dir()
        if base is None:
            return {}
        return {
            "enrichment_matrix": base / "cross_rbp_enrichment_matrix.tsv",
            "matched_gene_matrix": base / "cross_rbp_matched_gene_matrix.tsv",
            "comparison_summary": base / "cross_rbp_comparison_summary.tsv",
            "parameters": base / "parameters.json",
        }

    def isoform_support_output_dir(self, comparison_id: str | None) -> Path | None:
        if self.current_project is None or self.current_project.output_root is None or not comparison_id:
            return None
        return self.current_project.output_root / "08_isoform_candidate_support" / self._safe_output_name(comparison_id)

    def isoform_support_file_map(self, comparison_id: str | None) -> dict[str, Path]:
        base = self.isoform_support_output_dir(comparison_id)
        if base is None:
            return {}
        return {
            "gene_summary": base / "candidate_isoform_gene_summary.tsv",
            "transcript_table": base / "candidate_isoform_transcript_support.tsv",
            "parameters": base / "parameters.json",
        }

    def isoform_cross_support_output_dir(self) -> Path | None:
        if self.current_project is None or self.current_project.output_root is None:
            return None
        return self.current_project.output_root / "08_isoform_candidate_support" / "cross_comparison"

    def isoform_cross_support_plot_dir(self) -> Path | None:
        base = self.isoform_cross_support_output_dir()
        return base / "plots" if base is not None else None

    def isoform_cross_support_file_map(self) -> dict[str, Path]:
        base = self.isoform_cross_support_output_dir()
        if base is None:
            return {}
        return {
            "gene_support_matrix": base / "cross_isoform_gene_support_matrix.tsv",
            "transcript_support_matrix": base / "cross_isoform_transcript_support_matrix.tsv",
            "dominant_isoform_switches": base / "cross_isoform_dominant_isoform_switches.tsv",
            "parameters": base / "parameters.json",
        }

    def _cached_parameters_match(self, path: Path | None, expected: dict[str, object]) -> bool:
        if path is None or not path.exists():
            return False
        try:
            cached = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            return False
        for key, value in expected.items():
            if str(cached.get(key, "")) != str(value):
                return False
        return True

    def cross_comparison_output_dir(self) -> Path | None:
        if self.current_project is None or self.current_project.output_root is None:
            return None
        return CrossComparisonService(self.current_project.output_root).base_dir()

    def cross_comparison_file_map(self) -> dict[str, Path]:
        if self.current_project is None or self.current_project.output_root is None:
            return {}
        return CrossComparisonService(self.current_project.output_root).candidate_files()

    def cross_comparison_pattern_output_dir(self, pair_id: str | None) -> Path | None:
        if self.current_project is None or self.current_project.output_root is None:
            return None
        return CrossComparisonService(self.current_project.output_root).pattern_dir(pair_id)

    def cross_comparison_pattern_file_map(self, pair_id: str | None) -> dict[str, Path]:
        if self.current_project is None or self.current_project.output_root is None:
            return {}
        return CrossComparisonService(self.current_project.output_root).pattern_files(pair_id)

    def cross_comparison_interaction_output_dir(self, pair_id: str | None) -> Path | None:
        if self.current_project is None or self.current_project.output_root is None:
            return None
        return CrossComparisonService(self.current_project.output_root).interaction_dir(pair_id)

    def cross_comparison_interaction_file_map(self, pair_id: str | None) -> dict[str, Path]:
        if self.current_project is None or self.current_project.output_root is None:
            return {}
        return CrossComparisonService(self.current_project.output_root).interaction_files(pair_id)

    def cross_comparison_interaction_batch_output_dir(self) -> Path | None:
        if self.current_project is None or self.current_project.output_root is None:
            return None
        return CrossComparisonService(self.current_project.output_root).interaction_batch_dir()

    def cross_comparison_interaction_batch_file_map(self) -> dict[str, Path]:
        if self.current_project is None or self.current_project.output_root is None:
            return {}
        return CrossComparisonService(self.current_project.output_root).interaction_batch_files()

    def _write_cross_comparison_pattern_outputs(
        self,
        pair: ComparisonPairDefinition,
        tables: dict[str, pd.DataFrame],
        abs_cutoff: float,
        delta_cutoff: float,
    ) -> None:
        if self.current_project is None or self.current_project.output_root is None:
            return
        generated = CrossComparisonService(self.current_project.output_root).write_pattern_outputs(
            pair=pair,
            tables=tables,
            abs_cutoff=abs_cutoff,
            delta_cutoff=delta_cutoff,
        )
        self._register_cross_comparison_outputs(
            generated,
            module_id="5.4.legacy_pattern",
            scope_id=pair.pair_id,
            comparison_ids=(pair.comparison_a, pair.comparison_b),
            parameters={"abs_dpsi_cutoff": abs_cutoff, "large_delta_dpsi_cutoff": delta_cutoff},
        )

    def _write_cross_comparison_interaction_outputs(
        self,
        pair: ComparisonPairDefinition,
        tables: dict[str, pd.DataFrame],
        abs_cutoff: float,
        delta_cutoff: float,
        min_total_reads: int,
    ) -> None:
        if self.current_project is None or self.current_project.output_root is None:
            return
        generated = CrossComparisonService(self.current_project.output_root).write_interaction_outputs(
            pair=pair,
            tables=tables,
            abs_cutoff=abs_cutoff,
            delta_cutoff=delta_cutoff,
            min_total_reads=min_total_reads,
        )
        self._register_cross_comparison_outputs(
            generated,
            module_id="5.4.3",
            scope_id=pair.pair_id,
            comparison_ids=(pair.comparison_a, pair.comparison_b),
            parameters={
                "abs_dpsi_cutoff": abs_cutoff,
                "large_delta_dpsi_cutoff": delta_cutoff,
                "min_total_reads": min_total_reads,
            },
        )

    def _write_cross_comparison_interaction_batch_outputs(
        self,
        tables: dict[str, pd.DataFrame],
        abs_cutoff: float,
        delta_cutoff: float,
        min_total_reads: int,
    ) -> None:
        if self.current_project is None or self.current_project.output_root is None:
            return
        generated = CrossComparisonService(self.current_project.output_root).write_interaction_batch_outputs(
            tables=tables,
            abs_cutoff=abs_cutoff,
            delta_cutoff=delta_cutoff,
            min_total_reads=min_total_reads,
        )
        comparison_ids = tuple(comparison.comparison_id for comparison in self._selected_comparisons())
        self._register_cross_comparison_outputs(
            generated,
            module_id="5.4.3",
            scope_id="all_pairs",
            comparison_ids=comparison_ids,
            parameters={
                "abs_dpsi_cutoff": abs_cutoff,
                "large_delta_dpsi_cutoff": delta_cutoff,
                "min_total_reads": min_total_reads,
                "batch": True,
            },
        )

    def _register_cross_comparison_outputs(
        self,
        generated: list[Path],
        *,
        module_id: str,
        scope_id: str,
        comparison_ids: tuple[str, ...],
        parameters: dict[str, object],
    ) -> None:
        if not generated:
            return
        for path in generated:
            suffix = path.suffix.lower()
            kind = "table" if suffix in {".tsv", ".csv", ".xlsx"} else "parameters" if suffix == ".json" else "other"
            self.artifact_registry().register(
                module_id=module_id,
                artifact_key=f"{scope_id}/{path.name}",
                kind=kind,
                entity_type="event",
                path=path,
                scope="comparison_pair" if scope_id != "all_pairs" else "comparison_set",
                scope_id=scope_id,
                comparison_ids=comparison_ids,
                consumer_modules=["5.5", "7.2", "8.2", "9"],
                description="Cross-comparison AS event output registered for downstream reuse and project reporting.",
                parameters=parameters,
            )

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
            f"stable_AS mean(IJC+SJC)>={self.current_thresholds.stable_min_junction_reads:g}, "
            f"stable_AS reps/group>={self.current_thresholds.stable_min_replicates_per_group:g}, "
            f"stable_AS direction_consistency>={self.current_thresholds.stable_min_pairwise_consistency:g}, "
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
            lambda row: self._gene_lookup_key(
                self._safe_text(row.get("normalized_gene_id"), self._safe_text(row.get("gene_id"))),
                row.get("gene_symbol"),
            ),
            axis=1,
        )
        self._expression_support_cache = working
        return working.copy()

    def _comparison_expression_support_frame(self, comparison_id: str | None) -> pd.DataFrame:
        if self.current_project is None or not comparison_id:
            return pd.DataFrame()
        cached = self._comparison_expression_support_cache.get(comparison_id)
        if cached is not None:
            return cached.copy()
        comparison = next(
            (item for item in self.current_project.available_comparisons if item.comparison_id == comparison_id),
            None,
        )
        counts_path = self.current_project.counts_path
        if comparison is None or counts_path is None or not counts_path.exists():
            self._comparison_expression_support_cache[comparison_id] = pd.DataFrame()
            return pd.DataFrame()
        try:
            frame = pd.read_csv(counts_path, sep="\t", low_memory=False)
        except Exception:
            self._comparison_expression_support_cache[comparison_id] = pd.DataFrame()
            return pd.DataFrame()
        if frame.empty:
            self._comparison_expression_support_cache[comparison_id] = pd.DataFrame()
            return pd.DataFrame()

        experiment_group, control_group = self._comparison_deg_direction_labels(comparison)
        if not experiment_group or not control_group:
            self._comparison_expression_support_cache[comparison_id] = pd.DataFrame()
            return pd.DataFrame()

        def _match_group_columns(group_label: str) -> list[str]:
            import re
            label = str(group_label or "").strip()
            if not label:
                return []
            pattern = re.compile(rf"^{re.escape(label)}(?:$|[_\-.])", flags=re.IGNORECASE)
            return [column for column in frame.columns if column != "gene_id" and pattern.search(str(column))]

        exp_columns = _match_group_columns(experiment_group)
        ctrl_columns = _match_group_columns(control_group)
        if not exp_columns or not ctrl_columns:
            self._comparison_expression_support_cache[comparison_id] = pd.DataFrame()
            return pd.DataFrame()

        working = frame.copy()
        working["gene_id"] = self._normalized_gene_id_series(working)
        annotation_lookup = self._annotation_gene_symbol_lookup()
        working["gene_symbol"] = (
            working["gene_id"].map(annotation_lookup)
            if annotation_lookup
            else pd.Series(pd.NA, index=working.index, dtype="object")
        )
        working["gene_symbol"] = working["gene_symbol"].combine_first(working["gene_id"])
        exp_values = working[exp_columns].apply(pd.to_numeric, errors="coerce")
        ctrl_values = working[ctrl_columns].apply(pd.to_numeric, errors="coerce")
        working["expr_experiment"] = exp_values.mean(axis=1, skipna=True)
        working["expr_control"] = ctrl_values.mean(axis=1, skipna=True)
        working["gene_card_key"] = working.apply(
            lambda row: self._gene_lookup_key(row.get("gene_id"), row.get("gene_symbol")),
            axis=1,
        )

        rows: list[dict[str, object]] = []
        for _, row in working.dropna(subset=["gene_card_key"]).iterrows():
            rows.append(
                {
                    "comparison_id": comparison_id,
                    "comparison_display_name": comparison.display_resolved_name,
                    "gene_card_key": row.get("gene_card_key"),
                    "gene_id": row.get("gene_id"),
                    "gene_symbol": row.get("gene_symbol"),
                    "group": experiment_group,
                    "group_role": "experiment",
                    "expr": row.get("expr_experiment"),
                    "sample_count": len(exp_columns),
                    "expression_source": "counts_path_group_mean",
                }
            )
            rows.append(
                {
                    "comparison_id": comparison_id,
                    "comparison_display_name": comparison.display_resolved_name,
                    "gene_card_key": row.get("gene_card_key"),
                    "gene_id": row.get("gene_id"),
                    "gene_symbol": row.get("gene_symbol"),
                    "group": control_group,
                    "group_role": "control",
                    "expr": row.get("expr_control"),
                    "sample_count": len(ctrl_columns),
                    "expression_source": "counts_path_group_mean",
                }
            )
        result = pd.DataFrame(rows)
        self._comparison_expression_support_cache[comparison_id] = result
        return result.copy()

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
            working = self._apply_stable_as_columns(working, dpsi_column="dPSI", reverse_splicing=bool(comparison.reverse_splicing))
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
            working["significant_event"] = working["stable_as_event"].fillna(False).map(self._safe_bool)
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
                        "rmats_sig_before_stability",
                        "stable_as_event",
                        "mean_junction_reads_group1",
                        "mean_junction_reads_group2",
                        "replicate_count_group1",
                        "replicate_count_group2",
                        "supported_replicate_count_group1",
                        "supported_replicate_count_group2",
                        "replicate_direction_consistency",
                        "stable_support_status",
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
                "gene_symbol": self._safe_text(row.get("gene_symbol"), self._safe_text(row.get("gene_id"))),
                "gene_id": self._safe_text(row.get("gene_id")),
                "label": (
                    f"{self._safe_text(row.get('gene_symbol'), self._safe_text(row.get('gene_id')))} [{self._safe_text(row.get('gene_id'))}]"
                    if pd.notna(row.get("gene_id")) and str(row.get("gene_id")).strip()
                    else self._safe_text(row.get("gene_symbol"), self._safe_text(row.get("gene_id")))
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
        gene_id = self._safe_text(candidate_row.get("gene_id"), self._safe_text(fallback_gene_id))
        gene_symbol = self._safe_text(candidate_row.get("gene_symbol"), self._safe_text(fallback_gene_symbol, gene_id))

        deg_frame = self._deg_expression_frame(comparison_id)
        deg_row = pd.Series(dtype="object")
        if not deg_frame.empty:
            subset = deg_frame.loc[deg_frame["gene_card_key"].astype(str) == current_key].copy()
            if not subset.empty:
                deg_row = subset.iloc[0]

        expr_support = self._comparison_expression_support_frame(comparison_id)
        if expr_support.empty:
            expr_support = self._expression_support_frame()
        expr_subset = pd.DataFrame()
        if not expr_support.empty:
            expr_subset = expr_support.loc[expr_support["gene_card_key"].astype(str) == current_key].copy()
            if "comparison_id" in expr_subset.columns and comparison_id:
                expr_subset = expr_subset.loc[expr_subset["comparison_id"].astype(str) == str(comparison_id)].copy()
            if expr_subset.empty:
                gene_symbol_key = str(gene_symbol or "").strip().casefold()
                gene_id_key = normalize_gene_id(gene_id)
                fallback_mask = pd.Series(False, index=expr_support.index)
                if gene_symbol_key and "gene_symbol" in expr_support.columns:
                    fallback_mask |= expr_support["gene_symbol"].fillna("").astype(str).str.casefold().eq(gene_symbol_key)
                if gene_id_key and "gene_id" in expr_support.columns:
                    fallback_mask |= self._normalized_gene_id_series(expr_support).fillna("").astype(str).eq(gene_id_key)
                expr_subset = expr_support.loc[fallback_mask].copy()
                if "comparison_id" in expr_subset.columns and comparison_id:
                    expr_subset = expr_subset.loc[expr_subset["comparison_id"].astype(str) == str(comparison_id)].copy()
        comparison_expression_rows: list[dict[str, object]] = []
        comparison_dpsi_rows: list[dict[str, object]] = []
        event_frames_by_comparison: dict[str, tuple[pd.DataFrame, str]] = {}
        for comparison in selected_comparisons:
            comparison_expr = self._comparison_expression_support_frame(comparison.comparison_id)
            if comparison_expr.empty:
                subset = pd.DataFrame()
            else:
                subset = comparison_expr.loc[comparison_expr["gene_card_key"].astype(str) == current_key].copy()
                if subset.empty:
                    gene_symbol_key = str(gene_symbol or "").strip().casefold()
                    gene_id_key = normalize_gene_id(gene_id)
                    fallback_mask = pd.Series(False, index=comparison_expr.index)
                    if gene_symbol_key and "gene_symbol" in comparison_expr.columns:
                        fallback_mask |= comparison_expr["gene_symbol"].fillna("").astype(str).str.casefold().eq(gene_symbol_key)
                    if gene_id_key and "gene_id" in comparison_expr.columns:
                        fallback_mask |= self._normalized_gene_id_series(comparison_expr).fillna("").astype(str).eq(gene_id_key)
                    subset = comparison_expr.loc[fallback_mask].copy()
            if not subset.empty:
                for _, expr_row in subset.iterrows():
                    expr_value = pd.to_numeric(expr_row.get("expr"), errors="coerce")
                    if pd.isna(expr_value):
                        continue
                    comparison_expression_rows.append(
                        {
                            "comparison_id": comparison.comparison_id,
                            "comparison_display_name": comparison.display_resolved_name,
                            "group": expr_row.get("group"),
                            "group_role": expr_row.get("group_role", ""),
                            "expr": float(expr_value),
                            "sample_count": expr_row.get("sample_count"),
                            "expression_source": expr_row.get("expression_source", "counts_path_group_mean"),
                            "current_comparison": bool(comparison.comparison_id == comparison_id),
                        }
                    )

            comparison_events, _comparison_rule = self._candidate_card_event_frame(comparison.comparison_id, current_key)
            event_frames_by_comparison[comparison.comparison_id] = (comparison_events.copy(), _comparison_rule)
            significant_events = (
                comparison_events.loc[comparison_events["significant_event"].fillna(False)].copy()
                if not comparison_events.empty and "significant_event" in comparison_events.columns
                else pd.DataFrame()
            )
            dominant_for_comparison = pd.Series(dtype="object")
            if not significant_events.empty:
                significant_events = significant_events.assign(
                    abs_dPSI=pd.to_numeric(significant_events.get("dPSI"), errors="coerce").abs()
                ).sort_values(
                    by=["abs_dPSI", "FDR", "event_type", "event_id"],
                    ascending=[False, True, True, True],
                    na_position="last",
                )
                dominant_for_comparison = significant_events.iloc[0]
            comparison_dpsi_rows.append(
                {
                    "comparison_id": comparison.comparison_id,
                    "comparison_display_name": comparison.display_resolved_name,
                    "current_comparison": bool(comparison.comparison_id == comparison_id),
                    "n_significant_as_events": int(len(significant_events)),
                    "dominant_event_type": dominant_for_comparison.get("event_type", pd.NA),
                    "dominant_event_id": dominant_for_comparison.get("event_id", pd.NA),
                    "dominant_dPSI": dominant_for_comparison.get("dPSI", pd.NA),
                    "dominant_FDR": dominant_for_comparison.get("FDR", pd.NA),
                    "dominant_inclusion_direction": dominant_for_comparison.get("inclusion_direction", pd.NA),
                    "dominant_psi_experiment": dominant_for_comparison.get("psi_experiment", pd.NA),
                    "dominant_psi_control": dominant_for_comparison.get("psi_control", pd.NA),
                }
            )
        if str(comparison_id or "") in event_frames_by_comparison:
            event_frame, dominant_rule = event_frames_by_comparison[str(comparison_id or "")]
        else:
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
            "comparison_expression": pd.DataFrame(comparison_expression_rows),
            "comparison_dpsi": pd.DataFrame(comparison_dpsi_rows),
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
        filtered = self._dedupe_candidate_frame(filtered)
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

        selected = self._dedupe_candidate_frame(selected)
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
        cache_key = (
            str(comparison_id or ""),
            int(len(frame)) if frame is not None else 0,
            tuple(str(column) for column in (frame.columns.tolist() if frame is not None else [])),
        )
        cached = self._candidate_catalog_cache.get(cache_key)
        if cached is not None:
            return cached.copy()
        filtered = self._filter_frame_by_comparison(frame, comparison_id)
        filtered = self._dedupe_candidate_frame(filtered)
        if filtered.empty:
            return pd.DataFrame(columns=["comparison_id", "rank", "ranking_score", "gene_symbol", "gene_id", "candidate_tier", "evidence_class", "direction_class"])
        filtered = filtered.sort_values(["rank", "gene_symbol"], ascending=[True, True], na_position="last").copy()
        columns = [
            column for column in [
                "comparison_id",
                "rank",
                "ranking_score",
                "gene_symbol",
                "gene_id",
                "candidate_tier",
                "evidence_class",
                "direction_class",
                "n_rMATS_significant_events",
                "significant_as_event_types",
                "significant_as_event_ids",
                "significant_as_event_list",
                "as_event_summary",
                "cross_candidate_pattern",
                "cross_pattern_reason",
                "cross_candidate_comparisons",
                "cross_n_comparisons_candidate",
                "cross_best_overall_rank",
                "cross_comparison_effect_summary",
                "shared_candidate",
                "specific_candidate",
                "gained_candidate",
                "lost_candidate",
                "dominant_rMATS_event_type",
                "best_rMATS_event_id",
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
        result = filtered[columns].copy().reset_index(drop=True)
        self._candidate_catalog_cache[cache_key] = result.copy()
        return result

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
        frame = self._dedupe_candidate_frame(frame)
        if frame.empty:
            return []
        working = frame.sort_values(["rank", "gene_symbol"], ascending=[True, True], na_position="last").head(top_n).copy()
        genes: list[str] = []
        seen: set[str] = set()
        for _, row in working.iterrows():
            gene = self._safe_text(row.get("gene_symbol"), self._safe_text(row.get("gene_id"))).strip()
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
        self.current_project.candidate_selection_genes[comparison_id] = self._dedupe_gene_tokens(
            self.current_project.candidate_selection_genes.get(comparison_id, [])
        )
        self._clear_candidate_selection_caches()
        self.save_project_config()

    def _candidate_selected_gene_rows(self, filtered: pd.DataFrame, comparison_id: str | None, *, top_n: int = 20, allow_generate: bool = True) -> pd.DataFrame:
        if self.current_project is None or not comparison_id:
            return filtered.sort_values(["rank", "gene_symbol"], ascending=[True, True], na_position="last").head(top_n).copy()
        self._ensure_default_candidate_selection(comparison_id, top_n=top_n, allow_generate=allow_generate)
        selected_genes = self._dedupe_gene_tokens(self.current_project.candidate_selection_genes.get(comparison_id, []))
        if not selected_genes:
            return pd.DataFrame(columns=filtered.columns)
        selected = self._candidate_custom_gene_rows(filtered, comparison_id, selected_genes)
        return selected.sort_values(["rank", "gene_symbol"], ascending=[True, True], na_position="last").reset_index(drop=True)

    def _candidate_support_table_is_stale(self, frame: pd.DataFrame) -> bool:
        if frame is None or frame.empty:
            return False
        required_columns = {
            "DEXSeq_significant",
            "DTU_significant",
            "significant_as_event_ids",
            "significant_as_event_list",
        }
        if not required_columns.issubset(set(frame.columns)):
            return True
        has_any_candidate_support = bool(
            frame["DEXSeq_significant"].map(self._safe_bool).any()
            or frame["DTU_significant"].map(self._safe_bool).any()
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

    def preview_rbp_support(
        self,
        comparison_id: str | None,
        *,
        source: str = "selection",
        top_n: int = 50,
        custom_gene_set_path: str | None = None,
        allow_generate: bool = False,
        force_rebuild: bool = False,
    ) -> dict[str, pd.DataFrame]:
        if self.current_project is None:
            raise RuntimeError("No project loaded.")
        comparison = self._comparison_by_id(comparison_id)
        if comparison is None:
            return self._empty_rbp_support_tables("No comparison selected.")

        file_map = self.rbp_support_file_map(comparison.comparison_id)
        custom_key = str(Path(custom_gene_set_path).resolve()) if custom_gene_set_path else ""
        cache_key = (comparison.comparison_id, source, int(top_n), custom_key)
        if not force_rebuild and cache_key in self._rbp_support_cache:
            self._mark_module_state("rbp_support", "cache_hit", "Loaded RBP support from in-memory cache.")
            return {key: value.copy() for key, value in self._rbp_support_cache[cache_key].items()}
        expected_parameters = {
            "comparison_id": comparison.comparison_id,
            "candidate_source": source,
            "top_n": int(top_n),
            "custom_gene_set_path": custom_gene_set_path or "",
        }
        if (
            not force_rebuild
            and file_map.get("summary")
            and file_map["summary"].exists()
            and self._cached_parameters_match(file_map.get("parameters"), expected_parameters)
        ):
            tables = {
                "summary": self._read_table_or_empty(file_map.get("summary")),
                "matched_genes": self._read_table_or_empty(file_map.get("matched_genes")),
            }
            self._rbp_support_cache[cache_key] = {key: value.copy() for key, value in tables.items()}
            self._mark_module_state(
                "rbp_support",
                "cache_hit",
                "Loaded cached RBP support tables.",
                output_folder=str(self.rbp_support_output_dir(comparison.comparison_id) or ""),
                from_cache=True,
            )
            return tables
        if not allow_generate:
            self._mark_module_state("rbp_support", "not_run", "No cached RBP support result. Click Run RBP support.")
            return self._empty_rbp_support_tables("No cached RBP support result.")

        selected, selection_meta = self._candidate_rows_for_support(comparison.comparison_id, source=source, top_n=top_n)
        universe = self.candidate_selection_catalog(comparison.comparison_id, allow_generate=False)
        universe = self._dedupe_candidate_frame(universe)
        if selected.empty or universe.empty:
            tables = self._empty_rbp_support_tables("Missing input: candidate ranking/selection is not available for this comparison.")
            self._mark_module_state("rbp_support", "failed", tables["summary"].iloc[0]["note"] if not tables["summary"].empty else "Missing candidate input.")
            return tables

        gene_sets = self._load_rbp_gene_sets(custom_gene_set_path)
        if gene_sets.empty:
            tables = self._empty_rbp_support_tables("No RBP gene set available. Provide a custom TSV/CSV or use the built-in seed set.")
            self._mark_module_state("rbp_support", "failed", tables["summary"].iloc[0]["note"])
            return tables

        universe_tokens = self._gene_token_set(universe)
        selected_tokens = self._gene_token_set(selected)
        selected_lookup = self._gene_token_to_rows(selected)
        summary_rows: list[dict[str, object]] = []
        match_rows: list[dict[str, object]] = []
        for (gene_set, source_name), subset in gene_sets.groupby(["gene_set", "source"], dropna=False):
            set_tokens = {
                token
                for token in subset["gene_token"].dropna().astype(str)
                if token and token in universe_tokens
            }
            overlap = sorted(selected_tokens & set_tokens)
            pvalue = self._hypergeom_enrichment_pvalue(
                universe_size=len(universe_tokens),
                set_size=len(set_tokens),
                selected_size=len(selected_tokens),
                overlap_size=len(overlap),
            )
            overlap_genes = sorted(
                {
                    self._safe_text(row.get("gene_symbol"), self._safe_text(row.get("gene_id")))
                    for token in overlap
                    for row in selected_lookup.get(token, [])
                }
            )
            summary_rows.append(
                {
                    "comparison_id": comparison.comparison_id,
                    "display_name": comparison.display_resolved_name,
                    "rbp_gene_set": self._safe_text(gene_set, "RBP gene set"),
                    "rbp_source": self._safe_text(source_name, "custom/built-in"),
                    "candidate_source": source,
                    "selection_rule": selection_meta.get("selection_rule", ""),
                    "universe_gene_count": len(universe_tokens),
                    "selected_gene_count": len(selected_tokens),
                    "gene_set_size_in_universe": len(set_tokens),
                    "overlap_count": len(overlap),
                    "overlap_genes": "; ".join(overlap_genes),
                    "pvalue": pvalue,
                    "note": "Hypergeometric enrichment over the current candidate-ranking universe; use as prioritization, not causal proof.",
                }
            )
            for token in overlap:
                for row in selected_lookup.get(token, []):
                    match = dict(row)
                    match.update(
                        {
                            "comparison_id": comparison.comparison_id,
                            "display_name": comparison.display_resolved_name,
                            "rbp_gene_set": self._safe_text(gene_set, "RBP gene set"),
                            "rbp_source": self._safe_text(source_name, "custom/built-in"),
                        }
                    )
                    match_rows.append(match)

        summary = pd.DataFrame(summary_rows)
        if "pvalue" in summary.columns:
            summary["FDR"] = self._bh_adjust(summary["pvalue"])
            summary = summary.sort_values(["FDR", "pvalue", "overlap_count"], ascending=[True, True, False], na_position="last").reset_index(drop=True)
        matched = pd.DataFrame(match_rows)
        if not matched.empty:
            front = [column for column in ["comparison_id", "display_name", "rbp_gene_set", "gene_symbol", "gene_id", "rank", "candidate_tier", "evidence_class", "direction_class"] if column in matched.columns]
            rest = [column for column in matched.columns if column not in front]
            matched = matched[front + rest].reset_index(drop=True)
        tables = {"summary": summary, "matched_genes": matched}
        self._write_rbp_support_outputs(comparison.comparison_id, tables, source, top_n, custom_gene_set_path)
        self._rbp_support_cache[cache_key] = {key: value.copy() for key, value in tables.items()}
        self._mark_module_state(
            "rbp_support",
            "finished",
            "RBP support/enrichment completed for the selected comparison.",
            output_folder=str(self.rbp_support_output_dir(comparison.comparison_id) or ""),
            from_cache=False,
        )
        return tables

    def preview_rbp_cross_comparison_support(
        self,
        *,
        source: str = "selection",
        top_n: int = 50,
        custom_gene_set_path: str | None = None,
        allow_generate: bool = False,
        force_rebuild: bool = False,
    ) -> dict[str, pd.DataFrame]:
        if self.current_project is None:
            raise RuntimeError("No project loaded.")
        comparisons = self._selected_comparisons()
        if not comparisons:
            self._mark_module_state("rbp_cross_support", "not_run", "Missing input: no selected comparisons.")
            return self._empty_rbp_cross_support_tables("No selected comparisons.")

        selected_comparison_ids = [comparison.comparison_id for comparison in comparisons]
        custom_key = str(Path(custom_gene_set_path).resolve()) if custom_gene_set_path else ""
        cache_key = (source, int(top_n), custom_key, tuple(selected_comparison_ids))
        if not force_rebuild and cache_key in self._rbp_cross_support_cache:
            self._mark_module_state(
                "rbp_cross_support",
                "cache_hit",
                "Loaded cross-comparison RBP support from in-memory cache.",
                from_cache=True,
                output_folder=str(self.rbp_cross_support_output_dir() or ""),
            )
            return {key: value.copy() for key, value in self._rbp_cross_support_cache[cache_key].items()}

        file_map = self.rbp_cross_support_file_map()
        expected_parameters = {
            "candidate_source": source,
            "top_n": int(top_n),
            "custom_gene_set_path": custom_gene_set_path or "",
            "selected_comparisons": selected_comparison_ids,
        }
        if (
            not force_rebuild
            and file_map.get("enrichment_matrix")
            and file_map["enrichment_matrix"].exists()
            and self._cached_parameters_match(file_map.get("parameters"), expected_parameters)
        ):
            tables = {
                "enrichment_matrix": self._read_table_or_empty(file_map.get("enrichment_matrix")),
                "matched_gene_matrix": self._read_table_or_empty(file_map.get("matched_gene_matrix")),
                "comparison_summary": self._read_table_or_empty(file_map.get("comparison_summary")),
            }
            self._rbp_cross_support_cache[cache_key] = {key: value.copy() for key, value in tables.items()}
            self._mark_module_state(
                "rbp_cross_support",
                "cache_hit",
                "Loaded cached cross-comparison RBP support tables.",
                from_cache=True,
                output_folder=str(self.rbp_cross_support_output_dir() or ""),
            )
            return tables

        if not allow_generate:
            self._mark_module_state(
                "rbp_cross_support",
                "not_run",
                "No cached cross-comparison RBP support result. Click Run Cross-comparison RBP support.",
            )
            return self._empty_rbp_cross_support_tables("No cached cross-comparison RBP support result.")

        summary_frames: list[pd.DataFrame] = []
        matched_frames: list[pd.DataFrame] = []
        missing_comparisons: list[str] = []
        for comparison in comparisons:
            tables = self.preview_rbp_support(
                comparison.comparison_id,
                source=source,
                top_n=top_n,
                custom_gene_set_path=custom_gene_set_path,
                allow_generate=False,
                force_rebuild=False,
            )
            summary = tables.get("summary", pd.DataFrame()).copy()
            matched = tables.get("matched_genes", pd.DataFrame()).copy()
            if summary.empty and matched.empty:
                missing_comparisons.append(comparison.display_resolved_name)
                continue
            if not summary.empty and "rbp_gene_set" in summary.columns:
                summary = summary.loc[summary["rbp_gene_set"].astype(str).str.strip().ne("")].copy()
                summary_frames.append(summary)
            if not matched.empty:
                matched_frames.append(matched.copy())

        summary_all = pd.concat(summary_frames, axis=0, ignore_index=True) if summary_frames else pd.DataFrame()
        matched_all = pd.concat(matched_frames, axis=0, ignore_index=True) if matched_frames else pd.DataFrame()
        if summary_all.empty:
            tables = self._empty_rbp_cross_support_tables("No per-comparison RBP support rows were available.")
        else:
            tables = {
                "enrichment_matrix": self._build_rbp_cross_enrichment_matrix(summary_all, comparisons),
                "matched_gene_matrix": self._build_rbp_cross_matched_gene_matrix(matched_all, comparisons),
                "comparison_summary": self._build_rbp_cross_comparison_summary(summary_all, comparisons),
            }

        self._write_rbp_cross_support_outputs(tables, source, top_n, custom_gene_set_path)
        self._rbp_cross_support_cache[cache_key] = {key: value.copy() for key, value in tables.items()}
        missing_note = (
            f" Missing per-comparison RBP support cache for: {', '.join(missing_comparisons)}."
            if missing_comparisons
            else ""
        )
        self._mark_module_state(
            "rbp_cross_support",
            "finished",
            "Cross-comparison RBP support completed from cached per-comparison RBP support tables. "
            "This is a support-pattern view, not causal RBP inference."
            + missing_note,
            from_cache=False,
            output_folder=str(self.rbp_cross_support_output_dir() or ""),
        )
        return {key: value.copy() for key, value in tables.items()}

    def preview_isoform_candidate_support(
        self,
        comparison_id: str | None,
        *,
        source: str = "selection",
        top_n: int = 50,
        allow_generate: bool = False,
        force_rebuild: bool = False,
    ) -> dict[str, pd.DataFrame]:
        if self.current_project is None:
            raise RuntimeError("No project loaded.")
        comparison = self._comparison_by_id(comparison_id)
        if comparison is None:
            return self._empty_isoform_support_tables("No comparison selected.")
        cache_key = (comparison.comparison_id, source, int(top_n))
        file_map = self.isoform_support_file_map(comparison.comparison_id)
        if not force_rebuild and cache_key in self._isoform_support_cache:
            self._mark_module_state("isoform_support", "cache_hit", "Loaded isoform support from in-memory cache.")
            return {key: value.copy() for key, value in self._isoform_support_cache[cache_key].items()}
        expected_parameters = {
            "comparison_id": comparison.comparison_id,
            "candidate_source": source,
            "top_n": int(top_n),
        }
        if (
            not force_rebuild
            and file_map.get("gene_summary")
            and file_map["gene_summary"].exists()
            and self._cached_parameters_match(file_map.get("parameters"), expected_parameters)
        ):
            tables = {
                "gene_summary": self._read_table_or_empty(file_map.get("gene_summary")),
                "transcript_table": self._read_table_or_empty(file_map.get("transcript_table")),
            }
            self._isoform_support_cache[cache_key] = {key: value.copy() for key, value in tables.items()}
            self._mark_module_state(
                "isoform_support",
                "cache_hit",
                "Loaded cached isoform candidate support tables.",
                output_folder=str(self.isoform_support_output_dir(comparison.comparison_id) or ""),
                from_cache=True,
            )
            return tables
        if not allow_generate:
            self._mark_module_state("isoform_support", "not_run", "No cached isoform candidate support result. Click Run Isoform support.")
            return self._empty_isoform_support_tables("No cached isoform candidate support result.")

        selected, selection_meta = self._candidate_rows_for_support(comparison.comparison_id, source=source, top_n=top_n)
        if selected.empty:
            tables = self._empty_isoform_support_tables("Missing input: candidate ranking/selection is not available for this comparison.")
            self._mark_module_state("isoform_support", "failed", tables["gene_summary"].iloc[0]["isoform_support_status"])
            return tables

        gtf_frame = self._gtf_transcript_gene_frame()
        quant_samples = self._comparison_isoform_samples(comparison)
        dtu_support = self._dtu_isoform_support_frame(comparison)
        if gtf_frame.empty:
            gene_summary = self._selected_isoform_gene_summary_without_quant(selected, comparison, "Missing GTF transcript-to-gene annotation.")
            transcript_table = pd.DataFrame(columns=self._isoform_transcript_columns())
            tables = {"gene_summary": gene_summary, "transcript_table": transcript_table}
            self._write_isoform_support_outputs(comparison.comparison_id, tables, source, top_n)
            self._mark_module_state("isoform_support", "failed", "Missing GTF transcript-to-gene annotation.")
            return tables
        if not quant_samples:
            gene_summary = self._selected_isoform_gene_summary_without_quant(selected, comparison, "No matching quant.sf files for this comparison.")
            transcript_table = pd.DataFrame(columns=self._isoform_transcript_columns())
            tables = {"gene_summary": gene_summary, "transcript_table": transcript_table}
            self._write_isoform_support_outputs(comparison.comparison_id, tables, source, top_n)
            self._mark_module_state("isoform_support", "failed", "No matching quant.sf files for this comparison.")
            return tables

        transcript_table = self._build_isoform_transcript_support(selected, comparison, quant_samples, gtf_frame, dtu_support)
        gene_summary = self._build_isoform_gene_summary(selected, comparison, transcript_table, selection_meta)
        tables = {"gene_summary": gene_summary, "transcript_table": transcript_table}
        self._write_isoform_support_outputs(comparison.comparison_id, tables, source, top_n)
        self._isoform_support_cache[cache_key] = {key: value.copy() for key, value in tables.items()}
        self._mark_module_state(
            "isoform_support",
            "finished",
            "Isoform candidate support completed from quant.sf + GTF + DTU/candidate evidence.",
            output_folder=str(self.isoform_support_output_dir(comparison.comparison_id) or ""),
            from_cache=False,
        )
        return tables

    def preview_isoform_cross_comparison_support(
        self,
        *,
        source: str = "selection",
        top_n: int = 50,
        allow_generate: bool = False,
        force_rebuild: bool = False,
    ) -> dict[str, pd.DataFrame]:
        if self.current_project is None:
            raise RuntimeError("No project loaded.")
        comparisons = self._selected_comparisons()
        if not comparisons:
            self._mark_module_state("isoform_cross_support", "not_run", "Missing input: no selected comparisons.")
            return self._empty_isoform_cross_support_tables("No selected comparisons.")

        selected_comparison_ids = [comparison.comparison_id for comparison in comparisons]
        cache_key = (source, int(top_n), tuple(selected_comparison_ids))
        if not force_rebuild and cache_key in self._isoform_cross_support_cache:
            self._mark_module_state(
                "isoform_cross_support",
                "cache_hit",
                "Loaded cross-comparison isoform support from in-memory cache.",
                from_cache=True,
                output_folder=str(self.isoform_cross_support_output_dir() or ""),
            )
            return {key: value.copy() for key, value in self._isoform_cross_support_cache[cache_key].items()}

        file_map = self.isoform_cross_support_file_map()
        expected_parameters = {
            "candidate_source": source,
            "top_n": int(top_n),
            "selected_comparisons": selected_comparison_ids,
        }
        if (
            not force_rebuild
            and file_map.get("gene_support_matrix")
            and file_map["gene_support_matrix"].exists()
            and self._cached_parameters_match(file_map.get("parameters"), expected_parameters)
        ):
            tables = {
                "gene_support_matrix": self._read_table_or_empty(file_map.get("gene_support_matrix")),
                "transcript_support_matrix": self._read_table_or_empty(file_map.get("transcript_support_matrix")),
                "dominant_isoform_switches": self._read_table_or_empty(file_map.get("dominant_isoform_switches")),
            }
            self._isoform_cross_support_cache[cache_key] = {key: value.copy() for key, value in tables.items()}
            self._mark_module_state(
                "isoform_cross_support",
                "cache_hit",
                "Loaded cached cross-comparison isoform support tables.",
                from_cache=True,
                output_folder=str(self.isoform_cross_support_output_dir() or ""),
            )
            return tables

        if not allow_generate:
            self._mark_module_state(
                "isoform_cross_support",
                "not_run",
                "No cached cross-comparison isoform support result. Click Run Cross-comparison Isoform support.",
            )
            return self._empty_isoform_cross_support_tables("No cached cross-comparison isoform support result.")

        gene_frames: list[pd.DataFrame] = []
        transcript_frames: list[pd.DataFrame] = []
        missing_comparisons: list[str] = []
        for comparison in comparisons:
            tables = self.preview_isoform_candidate_support(
                comparison.comparison_id,
                source=source,
                top_n=top_n,
                allow_generate=False,
                force_rebuild=False,
            )
            gene_summary = tables.get("gene_summary", pd.DataFrame()).copy()
            transcript_table = tables.get("transcript_table", pd.DataFrame()).copy()
            if gene_summary.empty and transcript_table.empty:
                missing_comparisons.append(comparison.display_resolved_name)
                continue
            if not gene_summary.empty and "comparison_id" in gene_summary.columns:
                gene_summary = gene_summary.loc[gene_summary["comparison_id"].astype(str).eq(comparison.comparison_id)].copy()
                gene_frames.append(gene_summary)
            if not transcript_table.empty and "comparison_id" in transcript_table.columns:
                transcript_table = transcript_table.loc[transcript_table["comparison_id"].astype(str).eq(comparison.comparison_id)].copy()
                transcript_frames.append(transcript_table)

        gene_all = pd.concat(gene_frames, axis=0, ignore_index=True) if gene_frames else pd.DataFrame()
        transcript_all = pd.concat(transcript_frames, axis=0, ignore_index=True) if transcript_frames else pd.DataFrame()
        if gene_all.empty:
            tables = self._empty_isoform_cross_support_tables("No per-comparison isoform support rows were available.")
        else:
            gene_matrix = self._build_isoform_cross_gene_matrix(gene_all, comparisons)
            transcript_matrix = self._build_isoform_cross_transcript_matrix(transcript_all, comparisons)
            switches = self._build_isoform_cross_dominant_switches(gene_matrix)
            tables = {
                "gene_support_matrix": gene_matrix,
                "transcript_support_matrix": transcript_matrix,
                "dominant_isoform_switches": switches,
            }

        self._write_isoform_cross_support_outputs(tables, source, top_n)
        self._isoform_cross_support_cache[cache_key] = {key: value.copy() for key, value in tables.items()}
        missing_note = (
            f" Missing per-comparison isoform support cache for: {', '.join(missing_comparisons)}."
            if missing_comparisons
            else ""
        )
        self._mark_module_state(
            "isoform_cross_support",
            "finished",
            "Cross-comparison isoform support completed from cached per-comparison isoform support tables. "
            "This compares support patterns and dominant isoforms across comparisons."
            + missing_note,
            from_cache=False,
            output_folder=str(self.isoform_cross_support_output_dir() or ""),
        )
        return {key: value.copy() for key, value in tables.items()}

    def rbp_cross_support_plot_files(self) -> pd.DataFrame:
        return self._cross_support_plot_file_frame(
            self.rbp_cross_support_output_dir(),
            self.rbp_cross_support_plot_dir(),
            expected_selected_comparisons=[comparison.comparison_id for comparison in self._selected_comparisons()],
        )

    def isoform_cross_support_plot_files(self) -> pd.DataFrame:
        return self._cross_support_plot_file_frame(
            self.isoform_cross_support_output_dir(),
            self.isoform_cross_support_plot_dir(),
            expected_selected_comparisons=[comparison.comparison_id for comparison in self._selected_comparisons()],
        )

    def _cross_support_plot_file_frame(
        self,
        output_dir: Path | None,
        plot_dir: Path | None,
        *,
        expected_selected_comparisons: list[str] | None = None,
    ) -> pd.DataFrame:
        columns = ["plot_name", "absolute_path", "file_size", "last_modified", "preview_kind"]
        if output_dir is None or plot_dir is None or not plot_dir.exists():
            return pd.DataFrame(columns=columns)
        parameters_path = output_dir / "parameters.json"
        preferred: list[Path] = []
        try:
            payload = json.loads(parameters_path.read_text(encoding="utf-8")) if parameters_path.exists() else {}
            cached_selected = [str(item) for item in payload.get("selected_comparisons", [])]
            expected = [str(item) for item in (expected_selected_comparisons or [])]
            if expected and cached_selected and cached_selected != expected:
                return pd.DataFrame(columns=columns)
            for item in payload.get("generated_plot_files", []):
                path = Path(str(item))
                if path.exists() and path.suffix.lower() == ".png":
                    preferred.append(path)
        except (OSError, json.JSONDecodeError):
            preferred = []
        paths = preferred or sorted(plot_dir.glob("*.png"))
        rows: list[dict[str, object]] = []
        for path in paths:
            try:
                stat = path.stat()
            except OSError:
                continue
            rows.append(
                {
                    "plot_name": path.name,
                    "absolute_path": str(path),
                    "file_size": stat.st_size,
                    "last_modified": datetime.fromtimestamp(stat.st_mtime).isoformat(timespec="seconds"),
                    "preview_kind": "image",
                }
            )
        return pd.DataFrame(rows, columns=columns)

    def _comparison_by_id(self, comparison_id: str | None):
        if self.current_project is None or not comparison_id:
            return None
        return next((item for item in self.current_project.available_comparisons if item.comparison_id == comparison_id), None)

    @staticmethod
    def _read_table_or_empty(path: Path | None) -> pd.DataFrame:
        if path is None or not path.exists() or path.stat().st_size <= 0:
            return pd.DataFrame()
        try:
            suffix = path.suffix.lower()
            sep = "," if suffix == ".csv" else "\t"
            return pd.read_csv(path, sep=sep, low_memory=False)
        except Exception:
            return pd.DataFrame()

    def _candidate_rows_for_support(self, comparison_id: str, *, source: str, top_n: int) -> tuple[pd.DataFrame, dict[str, object]]:
        if source == "selection":
            artifact_rows, artifact_meta = self.downstream_gene_selection(
                comparison_id,
                consumer_module="7/8",
            )
            if not artifact_rows.empty:
                catalog = self.candidate_selection_catalog(comparison_id, allow_generate=False)
                if not catalog.empty:
                    gene_keys = {
                        str(value).strip().casefold()
                        for column in ("gene_symbol", "gene_id")
                        if column in artifact_rows.columns
                        for value in artifact_rows[column].dropna().astype(str)
                        if str(value).strip()
                    }
                    mask = pd.Series(False, index=catalog.index)
                    for column in ("gene_symbol", "gene_id"):
                        if column in catalog.columns:
                            mask |= catalog[column].fillna("").astype(str).str.casefold().isin(gene_keys)
                    matched = self._dedupe_candidate_frame(catalog.loc[mask].copy())
                    if not matched.empty:
                        artifact_meta.update(
                            {
                                "gene_count": int(len(matched)),
                                "source": "artifact_registry_selection",
                            }
                        )
                        return matched, artifact_meta
                artifact_meta.update(
                    {
                        "gene_count": int(len(artifact_rows)),
                        "source": "artifact_registry_selection",
                    }
                )
                return self._dedupe_candidate_frame(artifact_rows), artifact_meta
        selected, meta = self.candidate_selection_frame(
            comparison_id,
            source=source,
            top_n=top_n,
            allow_generate=False,
        )
        if not selected.empty:
            return self._dedupe_candidate_frame(selected), meta
        catalog = self.candidate_selection_catalog(comparison_id, allow_generate=False)
        if catalog.empty:
            return pd.DataFrame(), {
                "comparison_id": comparison_id,
                "source": source,
                "selection_rule": "No candidate selection and no cached candidate ranking table.",
                "gene_count": 0,
            }
        fallback = self._dedupe_candidate_frame(catalog).sort_values(
            ["rank", "gene_symbol"],
            ascending=[True, True],
            na_position="last",
        ).head(top_n).copy()
        return fallback, {
            "comparison_id": comparison_id,
            "source": "cached_ranking_top_n",
            "selection_rule": f"No saved candidate selection; using top {top_n} unique genes from cached 5.3.1 ranking without writing back to 5.3.2.",
            "gene_count": int(len(fallback)),
            "top_n": int(top_n),
        }

    def _empty_rbp_support_tables(self, note: str) -> dict[str, pd.DataFrame]:
        summary = pd.DataFrame(
            [
                {
                    "comparison_id": "",
                    "display_name": "",
                    "rbp_gene_set": "",
                    "rbp_source": "",
                    "candidate_source": "",
                    "selection_rule": "",
                    "universe_gene_count": 0,
                    "selected_gene_count": 0,
                    "gene_set_size_in_universe": 0,
                    "overlap_count": 0,
                    "overlap_genes": "",
                    "pvalue": pd.NA,
                    "FDR": pd.NA,
                    "note": note,
                }
            ]
        )
        return {"summary": summary, "matched_genes": pd.DataFrame()}

    def _empty_isoform_support_tables(self, note: str) -> dict[str, pd.DataFrame]:
        gene_summary = pd.DataFrame(
            [
                {
                    "comparison_id": "",
                    "display_name": "",
                    "gene_symbol": "",
                    "gene_id": "",
                    "candidate_tier": "",
                    "rank": pd.NA,
                    "n_quantified_isoforms": 0,
                    "n_DTU_significant_isoforms": 0,
                    "best_DTU_qvalue": pd.NA,
                    "isoform_support_status": note,
                }
            ]
        )
        return {"gene_summary": gene_summary, "transcript_table": pd.DataFrame(columns=self._isoform_transcript_columns())}

    def _empty_rbp_cross_support_tables(self, note: str) -> dict[str, pd.DataFrame]:
        summary = pd.DataFrame([{"support_status": note}])
        return {
            "enrichment_matrix": summary.copy(),
            "matched_gene_matrix": pd.DataFrame([{"support_status": note}]),
            "comparison_summary": summary.copy(),
        }

    def _empty_isoform_cross_support_tables(self, note: str) -> dict[str, pd.DataFrame]:
        summary = pd.DataFrame([{"support_status": note}])
        return {
            "gene_support_matrix": summary.copy(),
            "transcript_support_matrix": pd.DataFrame([{"support_status": note}]),
            "dominant_isoform_switches": pd.DataFrame([{"support_status": note}]),
        }

    @staticmethod
    def _support_presence_pattern(present_comparisons: list[str], total_comparisons: int) -> tuple[str, str]:
        n_present = len(present_comparisons)
        if n_present <= 0:
            return "no_support", "No support detected in the selected comparisons."
        if n_present == total_comparisons and total_comparisons > 0:
            return "shared_all", "Support is present in all selected comparisons."
        if n_present == 1:
            return "comparison_specific", "Support is present in only one selected comparison."
        return "shared_subset", "Support is shared by a subset of selected comparisons."

    def _comparison_display_lookup(self, comparisons: list) -> dict[str, str]:
        return {comparison.comparison_id: comparison.display_resolved_name for comparison in comparisons}

    def _build_rbp_cross_enrichment_matrix(self, summary_all: pd.DataFrame, comparisons: list) -> pd.DataFrame:
        comparison_ids = [comparison.comparison_id for comparison in comparisons]
        display_lookup = self._comparison_display_lookup(comparisons)
        rows: list[dict[str, object]] = []
        group_cols = ["rbp_gene_set", "rbp_source"]
        for keys, subset in summary_all.groupby(group_cols, dropna=False):
            rbp_gene_set, rbp_source = keys
            row: dict[str, object] = {
                "rbp_gene_set": self._safe_text(rbp_gene_set),
                "rbp_source": self._safe_text(rbp_source),
            }
            present: list[str] = []
            overlap_union: set[str] = set()
            for comparison_id in comparison_ids:
                comp = subset.loc[subset["comparison_id"].astype(str).eq(comparison_id)].copy()
                overlap_count = int(pd.to_numeric(comp.get("overlap_count"), errors="coerce").fillna(0).max()) if not comp.empty else 0
                fdr = pd.to_numeric(comp.get("FDR"), errors="coerce").min() if not comp.empty and "FDR" in comp.columns else float("nan")
                genes_text = self._safe_text(comp.get("overlap_genes", pd.Series([""])).dropna().head(1).iloc[0], "") if not comp.empty and "overlap_genes" in comp.columns and not comp["overlap_genes"].dropna().empty else ""
                if overlap_count > 0:
                    present.append(comparison_id)
                    overlap_union.update(gene.strip() for gene in genes_text.split(";") if gene.strip())
                row[f"{comparison_id}_display_name"] = display_lookup.get(comparison_id, comparison_id)
                row[f"{comparison_id}_overlap_count"] = overlap_count
                row[f"{comparison_id}_FDR"] = fdr
                row[f"{comparison_id}_overlap_genes"] = genes_text
            pattern, reason = self._support_presence_pattern(present, len(comparison_ids))
            row.update(
                {
                    "support_pattern": pattern,
                    "pattern_reason": reason,
                    "support_comparisons": ";".join(present),
                    "n_comparisons_with_support": len(present),
                    "n_selected_comparisons": len(comparison_ids),
                    "best_FDR": pd.to_numeric(subset.get("FDR"), errors="coerce").min() if "FDR" in subset.columns else float("nan"),
                    "max_overlap_count": int(pd.to_numeric(subset.get("overlap_count"), errors="coerce").fillna(0).max()) if "overlap_count" in subset.columns else 0,
                    "overlap_gene_union": "; ".join(sorted(overlap_union)),
                    "interpretation": "RBP support pattern across comparisons; prioritization only, not causal RBP regulation.",
                }
            )
            rows.append(row)
        if not rows:
            return pd.DataFrame()
        return pd.DataFrame(rows).sort_values(["best_FDR", "max_overlap_count", "rbp_gene_set"], ascending=[True, False, True], na_position="last").reset_index(drop=True)

    def _build_rbp_cross_matched_gene_matrix(self, matched_all: pd.DataFrame, comparisons: list) -> pd.DataFrame:
        if matched_all.empty:
            return pd.DataFrame()
        comparison_ids = [comparison.comparison_id for comparison in comparisons]
        display_lookup = self._comparison_display_lookup(comparisons)
        working = matched_all.copy()
        working["gene_key"] = working.apply(lambda row: self._gene_lookup_key(row.get("gene_id"), row.get("gene_symbol")), axis=1)
        rows: list[dict[str, object]] = []
        group_cols = ["rbp_gene_set", "gene_key", "gene_symbol", "gene_id"]
        for keys, subset in working.groupby(group_cols, dropna=False):
            rbp_gene_set, _gene_key, gene_symbol, gene_id = keys
            row: dict[str, object] = {
                "rbp_gene_set": self._safe_text(rbp_gene_set),
                "gene_symbol": self._safe_text(gene_symbol),
                "gene_id": self._safe_text(gene_id),
            }
            present: list[str] = []
            ranks = []
            for comparison_id in comparison_ids:
                comp = subset.loc[subset["comparison_id"].astype(str).eq(comparison_id)].copy()
                if not comp.empty:
                    present.append(comparison_id)
                    ranks.extend(pd.to_numeric(comp.get("rank"), errors="coerce").dropna().tolist())
                row[f"{comparison_id}_display_name"] = display_lookup.get(comparison_id, comparison_id)
                row[f"{comparison_id}_present"] = not comp.empty
                row[f"{comparison_id}_rank"] = pd.to_numeric(comp.get("rank"), errors="coerce").min() if not comp.empty and "rank" in comp.columns else pd.NA
                row[f"{comparison_id}_candidate_tier"] = self._safe_text(comp.get("candidate_tier", pd.Series([""])).dropna().head(1).iloc[0], "") if not comp.empty and "candidate_tier" in comp.columns and not comp["candidate_tier"].dropna().empty else ""
                row[f"{comparison_id}_evidence_class"] = self._safe_text(comp.get("evidence_class", pd.Series([""])).dropna().head(1).iloc[0], "") if not comp.empty and "evidence_class" in comp.columns and not comp["evidence_class"].dropna().empty else ""
                row[f"{comparison_id}_direction_class"] = self._safe_text(comp.get("direction_class", pd.Series([""])).dropna().head(1).iloc[0], "") if not comp.empty and "direction_class" in comp.columns and not comp["direction_class"].dropna().empty else ""
            pattern, reason = self._support_presence_pattern(present, len(comparison_ids))
            row.update(
                {
                    "support_pattern": pattern,
                    "pattern_reason": reason,
                    "support_comparisons": ";".join(present),
                    "n_comparisons_with_support": len(present),
                    "best_rank": min(ranks) if ranks else pd.NA,
                    "interpretation": "RBP/splicing-factor candidate membership across comparisons; overlap does not prove causal regulation.",
                }
            )
            rows.append(row)
        return pd.DataFrame(rows).sort_values(["n_comparisons_with_support", "best_rank", "gene_symbol"], ascending=[False, True, True], na_position="last").reset_index(drop=True)

    def _build_rbp_cross_comparison_summary(self, summary_all: pd.DataFrame, comparisons: list) -> pd.DataFrame:
        rows: list[dict[str, object]] = []
        display_lookup = self._comparison_display_lookup(comparisons)
        for comparison in comparisons:
            subset = summary_all.loc[summary_all["comparison_id"].astype(str).eq(comparison.comparison_id)].copy()
            overlap_counts = pd.to_numeric(subset.get("overlap_count"), errors="coerce").fillna(0) if not subset.empty and "overlap_count" in subset.columns else pd.Series(dtype=float)
            fdr = pd.to_numeric(subset.get("FDR"), errors="coerce") if not subset.empty and "FDR" in subset.columns else pd.Series(dtype=float)
            rows.append(
                {
                    "comparison_id": comparison.comparison_id,
                    "display_name": display_lookup.get(comparison.comparison_id, comparison.comparison_id),
                    "n_rbp_gene_sets_tested": int(len(subset)),
                    "n_rbp_gene_sets_with_overlap": int((overlap_counts > 0).sum()) if not overlap_counts.empty else 0,
                    "best_FDR": fdr.min() if not fdr.empty else pd.NA,
                    "max_overlap_count": int(overlap_counts.max()) if not overlap_counts.empty else 0,
                    "interpretation": "Per-comparison RBP support summary used by the cross-comparison view.",
                }
            )
        return pd.DataFrame(rows)

    def _build_isoform_cross_gene_matrix(self, gene_all: pd.DataFrame, comparisons: list) -> pd.DataFrame:
        comparison_ids = [comparison.comparison_id for comparison in comparisons]
        display_lookup = self._comparison_display_lookup(comparisons)
        working = gene_all.copy()
        working["gene_key"] = working.apply(lambda row: self._gene_lookup_key(row.get("gene_id"), row.get("gene_symbol")), axis=1)
        rows: list[dict[str, object]] = []
        for (_gene_key, gene_symbol, gene_id), subset in working.groupby(["gene_key", "gene_symbol", "gene_id"], dropna=False):
            row: dict[str, object] = {
                "gene_symbol": self._safe_text(gene_symbol),
                "gene_id": self._safe_text(gene_id),
            }
            present: list[str] = []
            dtu_present: list[str] = []
            dominant_g1: dict[str, str] = {}
            dominant_g2: dict[str, str] = {}
            abs_deltas: list[float] = []
            ranks: list[float] = []
            for comparison_id in comparison_ids:
                comp = subset.loc[subset["comparison_id"].astype(str).eq(comparison_id)].copy()
                if not comp.empty:
                    present.append(comparison_id)
                    if comp.get("DTU_significant", pd.Series([False])).map(self._safe_bool).any():
                        dtu_present.append(comparison_id)
                    rank_values = pd.to_numeric(comp.get("rank"), errors="coerce").dropna()
                    ranks.extend(rank_values.tolist())
                    delta = pd.to_numeric(comp.get("delta_total_TPM_group1_minus_group2"), errors="coerce").dropna()
                    abs_deltas.extend(delta.abs().tolist())
                    dom1 = self._safe_text(comp.get("dominant_transcript_group1", pd.Series([""])).dropna().head(1).iloc[0], "") if "dominant_transcript_group1" in comp.columns and not comp["dominant_transcript_group1"].dropna().empty else ""
                    dom2 = self._safe_text(comp.get("dominant_transcript_group2", pd.Series([""])).dropna().head(1).iloc[0], "") if "dominant_transcript_group2" in comp.columns and not comp["dominant_transcript_group2"].dropna().empty else ""
                    if dom1:
                        dominant_g1[comparison_id] = dom1
                    if dom2:
                        dominant_g2[comparison_id] = dom2
                row[f"{comparison_id}_display_name"] = display_lookup.get(comparison_id, comparison_id)
                row[f"{comparison_id}_has_isoform_support"] = not comp.empty
                row[f"{comparison_id}_rank"] = pd.to_numeric(comp.get("rank"), errors="coerce").min() if not comp.empty and "rank" in comp.columns else pd.NA
                row[f"{comparison_id}_candidate_tier"] = self._safe_text(comp.get("candidate_tier", pd.Series([""])).dropna().head(1).iloc[0], "") if not comp.empty and "candidate_tier" in comp.columns and not comp["candidate_tier"].dropna().empty else ""
                row[f"{comparison_id}_evidence_class"] = self._safe_text(comp.get("evidence_class", pd.Series([""])).dropna().head(1).iloc[0], "") if not comp.empty and "evidence_class" in comp.columns and not comp["evidence_class"].dropna().empty else ""
                row[f"{comparison_id}_DTU_significant"] = comp.get("DTU_significant", pd.Series([False])).map(self._safe_bool).any() if not comp.empty else False
                row[f"{comparison_id}_best_DTU_qvalue"] = pd.to_numeric(comp.get("best_DTU_qvalue"), errors="coerce").min() if not comp.empty and "best_DTU_qvalue" in comp.columns else pd.NA
                row[f"{comparison_id}_total_TPM_group1"] = pd.to_numeric(comp.get("total_TPM_group1"), errors="coerce").sum() if not comp.empty and "total_TPM_group1" in comp.columns else pd.NA
                row[f"{comparison_id}_total_TPM_group2"] = pd.to_numeric(comp.get("total_TPM_group2"), errors="coerce").sum() if not comp.empty and "total_TPM_group2" in comp.columns else pd.NA
                row[f"{comparison_id}_delta_total_TPM_group1_minus_group2"] = pd.to_numeric(comp.get("delta_total_TPM_group1_minus_group2"), errors="coerce").sum() if not comp.empty and "delta_total_TPM_group1_minus_group2" in comp.columns else pd.NA
                row[f"{comparison_id}_dominant_transcript_group1"] = dominant_g1.get(comparison_id, "")
                row[f"{comparison_id}_dominant_transcript_group2"] = dominant_g2.get(comparison_id, "")
            pattern, reason = self._support_presence_pattern(present, len(comparison_ids))
            row.update(
                {
                    "isoform_support_pattern": pattern,
                    "pattern_reason": reason,
                    "support_comparisons": ";".join(present),
                    "n_comparisons_with_isoform_support": len(present),
                    "DTU_support_comparisons": ";".join(dtu_present),
                    "n_comparisons_with_DTU_support": len(dtu_present),
                    "best_rank": min(ranks) if ranks else pd.NA,
                    "max_abs_delta_total_TPM": max(abs_deltas) if abs_deltas else pd.NA,
                    "dominant_group1_by_comparison": "; ".join(f"{cid}:{tx}" for cid, tx in dominant_g1.items()),
                    "dominant_group2_by_comparison": "; ".join(f"{cid}:{tx}" for cid, tx in dominant_g2.items()),
                    "dominant_group1_n_distinct": len(set(dominant_g1.values())),
                    "dominant_group2_n_distinct": len(set(dominant_g2.values())),
                    "interpretation": "Isoform support pattern across comparisons based on quant.sf + GTF + optional DTU qvalues.",
                }
            )
            rows.append(row)
        return pd.DataFrame(rows).sort_values(["n_comparisons_with_DTU_support", "n_comparisons_with_isoform_support", "best_rank"], ascending=[False, False, True], na_position="last").reset_index(drop=True)

    def _build_isoform_cross_transcript_matrix(self, transcript_all: pd.DataFrame, comparisons: list) -> pd.DataFrame:
        if transcript_all.empty:
            return pd.DataFrame()
        comparison_ids = [comparison.comparison_id for comparison in comparisons]
        display_lookup = self._comparison_display_lookup(comparisons)
        rows: list[dict[str, object]] = []
        group_cols = ["gene_id", "gene_symbol", "transcript_id", "transcript_name", "transcript_type"]
        for keys, subset in transcript_all.groupby(group_cols, dropna=False):
            gene_id, gene_symbol, transcript_id, transcript_name, transcript_type = keys
            row: dict[str, object] = {
                "gene_symbol": self._safe_text(gene_symbol),
                "gene_id": self._safe_text(gene_id),
                "transcript_id": self._safe_text(transcript_id),
                "transcript_name": self._safe_text(transcript_name),
                "transcript_type": self._safe_text(transcript_type),
            }
            present: list[str] = []
            dtu_present: list[str] = []
            deltas: list[float] = []
            for comparison_id in comparison_ids:
                comp = subset.loc[subset["comparison_id"].astype(str).eq(comparison_id)].copy()
                if not comp.empty:
                    present.append(comparison_id)
                    if comp.get("DTU_significant", pd.Series([False])).map(self._safe_bool).any():
                        dtu_present.append(comparison_id)
                    deltas.extend(pd.to_numeric(comp.get("delta_TPM_group1_minus_group2"), errors="coerce").dropna().tolist())
                row[f"{comparison_id}_display_name"] = display_lookup.get(comparison_id, comparison_id)
                row[f"{comparison_id}_mean_TPM_group1"] = pd.to_numeric(comp.get("mean_TPM_group1"), errors="coerce").mean() if not comp.empty and "mean_TPM_group1" in comp.columns else pd.NA
                row[f"{comparison_id}_mean_TPM_group2"] = pd.to_numeric(comp.get("mean_TPM_group2"), errors="coerce").mean() if not comp.empty and "mean_TPM_group2" in comp.columns else pd.NA
                row[f"{comparison_id}_delta_TPM_group1_minus_group2"] = pd.to_numeric(comp.get("delta_TPM_group1_minus_group2"), errors="coerce").mean() if not comp.empty and "delta_TPM_group1_minus_group2" in comp.columns else pd.NA
                row[f"{comparison_id}_DTU_significant"] = comp.get("DTU_significant", pd.Series([False])).map(self._safe_bool).any() if not comp.empty else False
                row[f"{comparison_id}_DTU_qvalue"] = pd.to_numeric(comp.get("DTU_qvalue"), errors="coerce").min() if not comp.empty and "DTU_qvalue" in comp.columns else pd.NA
            pattern, reason = self._support_presence_pattern(present, len(comparison_ids))
            row.update(
                {
                    "transcript_support_pattern": pattern,
                    "pattern_reason": reason,
                    "support_comparisons": ";".join(present),
                    "n_comparisons_with_transcript_support": len(present),
                    "DTU_support_comparisons": ";".join(dtu_present),
                    "n_comparisons_with_DTU_support": len(dtu_present),
                    "max_abs_delta_TPM": max([abs(value) for value in deltas]) if deltas else pd.NA,
                    "interpretation": "Transcript-level isoform support across comparisons.",
                }
            )
            rows.append(row)
        return pd.DataFrame(rows).sort_values(["n_comparisons_with_DTU_support", "n_comparisons_with_transcript_support", "max_abs_delta_TPM"], ascending=[False, False, False], na_position="last").reset_index(drop=True)

    def _build_isoform_cross_dominant_switches(self, gene_matrix: pd.DataFrame) -> pd.DataFrame:
        if gene_matrix.empty or "gene_symbol" not in gene_matrix.columns:
            return pd.DataFrame()
        rows: list[dict[str, object]] = []
        for _, row in gene_matrix.iterrows():
            n_g1 = int(pd.to_numeric(pd.Series([row.get("dominant_group1_n_distinct")]), errors="coerce").fillna(0).iloc[0])
            n_g2 = int(pd.to_numeric(pd.Series([row.get("dominant_group2_n_distinct")]), errors="coerce").fillna(0).iloc[0])
            if n_g1 <= 1 and n_g2 <= 1:
                continue
            rows.append(
                {
                    "gene_symbol": row.get("gene_symbol", ""),
                    "gene_id": row.get("gene_id", ""),
                    "switch_type": "dominant_group1_and_group2" if n_g1 > 1 and n_g2 > 1 else "dominant_group1" if n_g1 > 1 else "dominant_group2",
                    "dominant_group1_n_distinct": n_g1,
                    "dominant_group2_n_distinct": n_g2,
                    "dominant_group1_by_comparison": row.get("dominant_group1_by_comparison", ""),
                    "dominant_group2_by_comparison": row.get("dominant_group2_by_comparison", ""),
                    "support_comparisons": row.get("support_comparisons", ""),
                    "n_comparisons_with_isoform_support": row.get("n_comparisons_with_isoform_support", 0),
                    "DTU_support_comparisons": row.get("DTU_support_comparisons", ""),
                    "best_rank": row.get("best_rank", pd.NA),
                    "interpretation": "Dominant isoform differs across comparisons; follow up with transcript-level table and sashimi/event context.",
                }
            )
        if not rows:
            return pd.DataFrame([{"support_status": "No dominant isoform switches found across selected comparisons."}])
        return pd.DataFrame(rows).sort_values(["n_comparisons_with_isoform_support", "best_rank", "gene_symbol"], ascending=[False, True, True], na_position="last").reset_index(drop=True)

    def _write_rbp_cross_support_plots(self, tables: dict[str, pd.DataFrame]) -> list[str]:
        plot_dir = self.rbp_cross_support_plot_dir()
        if plot_dir is None:
            return []
        plot_dir.mkdir(parents=True, exist_ok=True)
        self._clear_png_outputs(plot_dir)
        comparisons = self._selected_comparisons()
        if not comparisons:
            return []
        generated: list[Path] = []
        enrichment = tables.get("enrichment_matrix", pd.DataFrame()).copy()
        matched = tables.get("matched_gene_matrix", pd.DataFrame()).copy()
        if not enrichment.empty:
            generated.append(self._plot_rbp_cross_enrichment(enrichment, comparisons, plot_dir / "cross_rbp_enrichment_support.png"))
        if not matched.empty:
            generated.append(self._plot_rbp_cross_matched_matrix(matched, comparisons, plot_dir / "cross_rbp_matched_gene_matrix.png"))
            for comp_a, comp_b in combinations(comparisons, 2):
                filename = f"pairwise_rbp_rank_support_{self._safe_output_name(comp_a.comparison_id)}_vs_{self._safe_output_name(comp_b.comparison_id)}.png"
                generated.append(self._plot_pairwise_rbp_rank_support(matched, comp_a, comp_b, plot_dir / filename))
        if len(comparisons) < 2:
            generated.append(self._plot_no_data(plot_dir / "pairwise_rbp_rank_support_not_available.png", "Pairwise RBP rank support", "Need at least two comparisons for pairwise plots."))
        return [str(path) for path in generated if path.exists()]

    def _write_isoform_cross_support_plots(self, tables: dict[str, pd.DataFrame]) -> list[str]:
        plot_dir = self.isoform_cross_support_plot_dir()
        if plot_dir is None:
            return []
        plot_dir.mkdir(parents=True, exist_ok=True)
        self._clear_png_outputs(plot_dir)
        comparisons = self._selected_comparisons()
        if not comparisons:
            return []
        generated: list[Path] = []
        gene_matrix = tables.get("gene_support_matrix", pd.DataFrame()).copy()
        switches = tables.get("dominant_isoform_switches", pd.DataFrame()).copy()
        if not gene_matrix.empty:
            generated.append(self._plot_isoform_cross_gene_heatmap(gene_matrix, comparisons, plot_dir / "cross_isoform_gene_support_heatmap.png"))
            for comp_a, comp_b in combinations(comparisons, 2):
                filename = f"pairwise_isoform_delta_tpm_{self._safe_output_name(comp_a.comparison_id)}_vs_{self._safe_output_name(comp_b.comparison_id)}.png"
                generated.append(self._plot_pairwise_isoform_delta_tpm(gene_matrix, comp_a, comp_b, plot_dir / filename))
        if not switches.empty:
            generated.append(self._plot_isoform_dominant_switches(switches, plot_dir / "cross_isoform_dominant_isoform_switches.png"))
        else:
            generated.append(self._plot_no_data(plot_dir / "cross_isoform_dominant_isoform_switches.png", "Dominant isoform switch candidates", "No dominant isoform switches found."))
        if len(comparisons) < 2:
            generated.append(self._plot_no_data(plot_dir / "pairwise_isoform_delta_tpm_not_available.png", "Pairwise isoform delta TPM", "Need at least two comparisons for pairwise plots."))
        return [str(path) for path in generated if path.exists()]

    @staticmethod
    def _clear_png_outputs(plot_dir: Path) -> None:
        for path in plot_dir.glob("*.png"):
            try:
                path.unlink()
            except OSError:
                pass

    def _plot_no_data(self, path: Path, title: str, message: str) -> Path:
        fig, ax = plt.subplots(figsize=(8, 4.5))
        ax.text(0.5, 0.55, title, ha="center", va="center", fontsize=15, fontweight="bold")
        ax.text(0.5, 0.42, message, ha="center", va="center", fontsize=11, color="#4b5563")
        ax.axis("off")
        fig.tight_layout()
        fig.savefig(path, dpi=180)
        plt.close(fig)
        return path

    def _plot_display_name(self, frame: pd.DataFrame, comparison_id: str, fallback: str) -> str:
        column = f"{comparison_id}_display_name"
        if column in frame.columns:
            values = frame[column].dropna().astype(str)
            values = values.loc[values.str.strip().ne("")]
            if not values.empty:
                return values.iloc[0]
        return fallback

    def _plot_rbp_cross_enrichment(self, enrichment: pd.DataFrame, comparisons: list, path: Path) -> Path:
        if enrichment.empty:
            return self._plot_no_data(path, "Cross-comparison RBP support", "No RBP enrichment rows available.")
        row = enrichment.sort_values(["best_FDR", "max_overlap_count"], ascending=[True, False], na_position="last").iloc[0]
        rows: list[dict[str, object]] = []
        for comparison in comparisons:
            cid = comparison.comparison_id
            rows.append(
                {
                    "display_name": self._plot_display_name(enrichment, cid, comparison.display_resolved_name),
                    "overlap_count": pd.to_numeric(pd.Series([row.get(f"{cid}_overlap_count", 0)]), errors="coerce").fillna(0).iloc[0],
                    "FDR": pd.to_numeric(pd.Series([row.get(f"{cid}_FDR", np.nan)]), errors="coerce").iloc[0],
                }
            )
        frame = pd.DataFrame(rows)
        frame["neg_log10_FDR"] = -np.log10(pd.to_numeric(frame["FDR"], errors="coerce").clip(lower=1e-300))
        fig, ax1 = plt.subplots(figsize=(10.5, 6.2))
        x = np.arange(len(frame))
        bars = ax1.bar(x, frame["overlap_count"], color="#2f6f9f", width=0.62, label="RBP overlap genes")
        ax1.set_ylabel("Overlap candidate genes", fontsize=12, fontweight="bold")
        ax1.set_xticks(x)
        ax1.set_xticklabels(frame["display_name"], rotation=28, ha="right")
        gene_set = self._safe_text(row.get("rbp_gene_set", "RBP gene set"), "RBP gene set")
        ax1.set_title(f"Cross-comparison RBP support\n{gene_set} overlap across selected comparisons", fontsize=15, fontweight="bold", pad=18)
        for rect, value in zip(bars, frame["overlap_count"]):
            ax1.text(rect.get_x() + rect.get_width() / 2, rect.get_height() + 0.2, f"{int(value)}", ha="center", va="bottom", fontsize=10, fontweight="bold")
        ax2 = ax1.twinx()
        ax2.plot(x, frame["neg_log10_FDR"], color="#d95f02", marker="o", linewidth=2.5, label="-log10(FDR)")
        ax2.set_ylabel("-log10(FDR)", fontsize=12, fontweight="bold")
        ax1.grid(axis="y", alpha=0.25)
        lines, labels = ax1.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        ax1.legend(lines + lines2, labels + labels2, loc="upper right", frameon=False)
        fig.text(0.02, 0.02, "Interpretation: support-prioritization only. Overlap/FDR does not prove causal RBP regulation.", fontsize=9, color="#4b5563")
        fig.tight_layout(rect=[0.03, 0.06, 0.98, 0.95])
        fig.savefig(path, dpi=180)
        plt.close(fig)
        return path

    def _plot_rbp_cross_matched_matrix(self, matched: pd.DataFrame, comparisons: list, path: Path) -> Path:
        present_cols = [f"{comparison.comparison_id}_present" for comparison in comparisons if f"{comparison.comparison_id}_present" in matched.columns]
        if matched.empty or not present_cols:
            return self._plot_no_data(path, "Cross-comparison RBP matched gene matrix", "No matched RBP candidate genes available.")
        working = matched.copy()
        working["_support_n"] = working[present_cols].apply(lambda values: sum(str(value).lower() == "true" for value in values), axis=1)
        working["_best_rank_sort"] = pd.to_numeric(working.get("best_rank", np.nan), errors="coerce")
        show = working.sort_values(["_support_n", "_best_rank_sort", "gene_symbol"], ascending=[False, True, True], na_position="last").head(30).copy()
        matrix: list[list[int]] = []
        annotations: list[list[str]] = []
        labels: list[str] = []
        for comparison in comparisons:
            cid = comparison.comparison_id
            present_col = f"{cid}_present"
            if present_col not in show.columns:
                continue
            labels.append(self._plot_display_name(show, cid, comparison.display_resolved_name))
            values: list[int] = []
            texts: list[str] = []
            rank_col = f"{cid}_rank"
            for _, row in show.iterrows():
                present = str(row.get(present_col, "")).lower() == "true"
                values.append(1 if present else 0)
                rank = pd.to_numeric(pd.Series([row.get(rank_col, np.nan)]), errors="coerce").iloc[0] if rank_col in show.columns else np.nan
                texts.append(f"rank {int(rank)}" if present and pd.notna(rank) else "")
            matrix.append(values)
            annotations.append(texts)
        if not matrix:
            return self._plot_no_data(path, "Cross-comparison RBP matched gene matrix", "No comparison columns available.")
        mat = np.array(matrix).T
        ann = np.array(annotations).T
        fig, ax = plt.subplots(figsize=(11.5, max(5.5, 0.34 * len(show) + 2)))
        ax.imshow(mat, aspect="auto", cmap=ListedColormap(["#f3f4f6", "#197278"]), vmin=0, vmax=1)
        ax.set_yticks(np.arange(len(show)))
        ax.set_yticklabels(show["gene_symbol"].astype(str), fontsize=9)
        ax.set_xticks(np.arange(len(labels)))
        ax.set_xticklabels(labels, rotation=30, ha="right", fontsize=10)
        ax.set_title("Cross-comparison RBP matched gene matrix\nRows = RBP/splicing-factor candidate genes; columns = selected comparisons", fontsize=14, fontweight="bold", pad=16)
        for i in range(mat.shape[0]):
            for j in range(mat.shape[1]):
                if mat[i, j] == 1:
                    ax.text(j, i, ann[i, j], ha="center", va="center", fontsize=8, color="white")
        ax.set_xlabel("Comparison")
        ax.set_ylabel("Matched RBP/splicing-factor candidate gene")
        for spine in ax.spines.values():
            spine.set_visible(False)
        fig.tight_layout(rect=[0.03, 0.04, 0.98, 0.95])
        fig.savefig(path, dpi=180)
        plt.close(fig)
        return path

    def _plot_pairwise_rbp_rank_support(self, matched: pd.DataFrame, comp_a, comp_b, path: Path) -> Path:
        cid_a = comp_a.comparison_id
        cid_b = comp_b.comparison_id
        present_a = f"{cid_a}_present"
        present_b = f"{cid_b}_present"
        if matched.empty or present_a not in matched.columns or present_b not in matched.columns:
            return self._plot_no_data(path, "Pairwise RBP rank support", "No matched RBP candidate genes available for this pair.")
        frame = matched.copy()
        label_a = self._plot_display_name(frame, cid_a, comp_a.display_resolved_name)
        label_b = self._plot_display_name(frame, cid_b, comp_b.display_resolved_name)
        frame["_present_a"] = frame[present_a].astype(str).str.lower().eq("true")
        frame["_present_b"] = frame[present_b].astype(str).str.lower().eq("true")
        frame = frame.loc[frame["_present_a"] | frame["_present_b"]].copy()
        if frame.empty:
            return self._plot_no_data(path, "Pairwise RBP rank support", "No RBP candidate support in either comparison.")
        frame["_rank_a"] = pd.to_numeric(frame.get(f"{cid_a}_rank", np.nan), errors="coerce")
        frame["_rank_b"] = pd.to_numeric(frame.get(f"{cid_b}_rank", np.nan), errors="coerce")
        max_rank = np.nanmax(frame[["_rank_a", "_rank_b"]].to_numpy(dtype=float))
        if not np.isfinite(max_rank):
            max_rank = 1000
        missing_position = max_rank + max_rank * 0.12 + 10
        frame["_x"] = frame["_rank_a"].fillna(missing_position)
        frame["_y"] = frame["_rank_b"].fillna(missing_position)
        classes = []
        colors = []
        for _, row in frame.iterrows():
            if bool(row["_present_a"]) and bool(row["_present_b"]):
                classes.append("shared")
                colors.append("#197278")
            elif bool(row["_present_a"]):
                classes.append(f"{label_a} only")
                colors.append("#2f6f9f")
            else:
                classes.append(f"{label_b} only")
                colors.append("#d95f02")
        frame["_support_class"] = classes
        frame["_color"] = colors
        fig, ax = plt.subplots(figsize=(8.6, 7.2))
        for support_class, subset in frame.groupby("_support_class"):
            ax.scatter(subset["_x"], subset["_y"], s=80, label=support_class, c=subset["_color"], edgecolor="white", linewidth=0.8, alpha=0.9)
        label_frame = frame.sort_values(["_present_a", "_present_b", "best_rank"], ascending=[False, False, True], na_position="last").head(14)
        for _, row in label_frame.iterrows():
            ax.text(row["_x"] + 5, row["_y"] + 5, str(row.get("gene_symbol", "")), fontsize=8)
        ax.axvline(missing_position - 5, color="#9ca3af", linestyle="--", linewidth=1)
        ax.axhline(missing_position - 5, color="#9ca3af", linestyle="--", linewidth=1)
        ax.invert_xaxis()
        ax.invert_yaxis()
        ax.set_xlabel(f"{label_a} candidate rank (lower is better)")
        ax.set_ylabel(f"{label_b} candidate rank (lower is better)")
        ax.set_title(f"Pairwise RBP candidate support\n{label_a} vs {label_b}", fontsize=14, fontweight="bold", pad=14)
        ax.legend(frameon=False, loc="upper left")
        ax.grid(alpha=0.25)
        fig.text(0.02, 0.02, "Dashed edge = absent in that comparison. This compares candidate support patterns, not causal RBP regulation.", fontsize=8.5, color="#4b5563")
        fig.tight_layout(rect=[0.04, 0.06, 0.98, 0.95])
        fig.savefig(path, dpi=180)
        plt.close(fig)
        return path

    def _plot_isoform_cross_gene_heatmap(self, gene_matrix: pd.DataFrame, comparisons: list, path: Path) -> Path:
        support_cols = [f"{comparison.comparison_id}_has_isoform_support" for comparison in comparisons if f"{comparison.comparison_id}_has_isoform_support" in gene_matrix.columns]
        if gene_matrix.empty or not support_cols:
            return self._plot_no_data(path, "Cross-comparison isoform support", "No isoform support rows available.")
        working = gene_matrix.copy()
        working["_support_n"] = working[support_cols].apply(lambda values: sum(str(value).lower() == "true" for value in values), axis=1)
        working["_best_rank_sort"] = pd.to_numeric(working.get("best_rank", np.nan), errors="coerce")
        show = working.sort_values(["_support_n", "_best_rank_sort", "max_abs_delta_total_TPM"], ascending=[False, True, False], na_position="last").head(30).copy()
        labels: list[str] = []
        value_matrix: list[list[float]] = []
        text_matrix: list[list[str]] = []
        for comparison in comparisons:
            cid = comparison.comparison_id
            support_col = f"{cid}_has_isoform_support"
            if support_col not in show.columns:
                continue
            labels.append(self._plot_display_name(show, cid, comparison.display_resolved_name))
            delta_col = f"{cid}_delta_total_TPM_group1_minus_group2"
            dtu_col = f"{cid}_DTU_significant"
            values: list[float] = []
            texts: list[str] = []
            for _, row in show.iterrows():
                supported = str(row.get(support_col, "")).lower() == "true"
                delta = pd.to_numeric(pd.Series([row.get(delta_col, np.nan)]), errors="coerce").iloc[0]
                values.append(delta if supported and pd.notna(delta) else np.nan)
                dtu = self._safe_bool(row.get(dtu_col, False))
                texts.append((f"{delta:.1f}\nDTU" if dtu else f"{delta:.1f}") if supported and pd.notna(delta) else "")
            value_matrix.append(values)
            text_matrix.append(texts)
        values = np.array(value_matrix, dtype=float).T
        texts = np.array(text_matrix).T
        max_abs = np.nanmax(np.abs(values)) if np.isfinite(values).any() else 1
        fig, ax = plt.subplots(figsize=(12, max(6, 0.31 * len(show) + 2)))
        image = ax.imshow(np.ma.masked_invalid(values), aspect="auto", cmap="RdBu_r", vmin=-max_abs, vmax=max_abs)
        ax.set_yticks(np.arange(len(show)))
        ax.set_yticklabels(show["gene_symbol"].astype(str), fontsize=8.5)
        ax.set_xticks(np.arange(len(labels)))
        ax.set_xticklabels(labels, rotation=30, ha="right", fontsize=10)
        ax.set_title("Cross-comparison isoform support\nColor = delta total TPM (group1 - group2) for supported candidate genes", fontsize=14, fontweight="bold", pad=16)
        for i in range(values.shape[0]):
            for j in range(values.shape[1]):
                if texts[i, j]:
                    ax.text(j, i, texts[i, j], ha="center", va="center", fontsize=6.8, color="black")
        for spine in ax.spines.values():
            spine.set_visible(False)
        colorbar = fig.colorbar(image, ax=ax, shrink=0.8, pad=0.02)
        colorbar.set_label("delta total TPM: group1 - group2")
        fig.text(0.02, 0.02, "Interpretation: isoform support pattern from quant.sf + GTF + optional DTU q-values; follow-up support, not standalone causality.", fontsize=8.5, color="#4b5563")
        fig.tight_layout(rect=[0.03, 0.06, 0.98, 0.95])
        fig.savefig(path, dpi=180)
        plt.close(fig)
        return path

    def _plot_isoform_dominant_switches(self, switches: pd.DataFrame, path: Path) -> Path:
        if switches.empty or "gene_symbol" not in switches.columns:
            return self._plot_no_data(path, "Cross-comparison dominant isoform switch candidates", "No dominant isoform switches found.")
        show = switches.loc[~switches.get("support_status", pd.Series([""] * len(switches))).astype(str).str.contains("No dominant", case=False, na=False)].head(15).copy()
        if show.empty:
            return self._plot_no_data(path, "Cross-comparison dominant isoform switch candidates", "No dominant isoform switches found.")
        fig, ax = plt.subplots(figsize=(10.8, max(4.8, 0.48 * len(show) + 2)))
        y = np.arange(len(show))
        counts = pd.to_numeric(show.get("n_comparisons_with_isoform_support", 0), errors="coerce").fillna(0)
        ax.barh(y, counts, color="#7a5195")
        ax.set_yticks(y)
        ax.set_yticklabels(show["gene_symbol"].astype(str), fontsize=10)
        ax.invert_yaxis()
        ax.set_xlabel("Comparisons with isoform support")
        ax.set_title("Cross-comparison dominant isoform switch candidates\nGenes where dominant transcript differs across selected comparisons", fontsize=14, fontweight="bold", pad=16)
        for index, (_, row) in enumerate(show.iterrows()):
            label = f"{row.get('switch_type', '')} | best rank {row.get('best_rank', '')}"
            ax.text(counts.iloc[index] + 0.04, index, label, va="center", fontsize=9)
        ax.grid(axis="x", alpha=0.25)
        fig.tight_layout(rect=[0.03, 0.04, 0.98, 0.95])
        fig.savefig(path, dpi=180)
        plt.close(fig)
        return path

    def _plot_pairwise_isoform_delta_tpm(self, gene_matrix: pd.DataFrame, comp_a, comp_b, path: Path) -> Path:
        cid_a = comp_a.comparison_id
        cid_b = comp_b.comparison_id
        support_a = f"{cid_a}_has_isoform_support"
        support_b = f"{cid_b}_has_isoform_support"
        if gene_matrix.empty or support_a not in gene_matrix.columns or support_b not in gene_matrix.columns:
            return self._plot_no_data(path, "Pairwise isoform delta TPM", "No isoform support rows available for this pair.")
        frame = gene_matrix.copy()
        label_a = self._plot_display_name(frame, cid_a, comp_a.display_resolved_name)
        label_b = self._plot_display_name(frame, cid_b, comp_b.display_resolved_name)
        frame["_support_a"] = frame[support_a].astype(str).str.lower().eq("true")
        frame["_support_b"] = frame[support_b].astype(str).str.lower().eq("true")
        frame = frame.loc[frame["_support_a"] | frame["_support_b"]].copy()
        if frame.empty:
            return self._plot_no_data(path, "Pairwise isoform delta TPM", "No isoform support in either comparison.")
        frame["_delta_a"] = pd.to_numeric(frame.get(f"{cid_a}_delta_total_TPM_group1_minus_group2", np.nan), errors="coerce")
        frame["_delta_b"] = pd.to_numeric(frame.get(f"{cid_b}_delta_total_TPM_group1_minus_group2", np.nan), errors="coerce")
        frame["_max_abs_delta"] = frame[["_delta_a", "_delta_b"]].abs().max(axis=1)
        show = frame.sort_values(["_support_a", "_support_b", "_max_abs_delta"], ascending=[False, False, False], na_position="last").head(18).copy()
        y = np.arange(len(show))
        height = 0.38
        fig, ax = plt.subplots(figsize=(11, max(5.5, 0.36 * len(show) + 2)))
        ax.barh(y - height / 2, show["_delta_a"].fillna(0), height=height, color="#2f6f9f", label=label_a)
        ax.barh(y + height / 2, show["_delta_b"].fillna(0), height=height, color="#d95f02", label=label_b)
        ax.axvline(0, color="#6b7280", linewidth=1)
        ax.set_yticks(y)
        ax.set_yticklabels(show["gene_symbol"].astype(str), fontsize=9)
        ax.invert_yaxis()
        ax.set_xlabel("delta total TPM (group1 - group2)")
        ax.set_title(f"Pairwise isoform support comparison\n{label_a} vs {label_b}", fontsize=14, fontweight="bold", pad=14)
        ax.legend(frameon=False, loc="upper right")
        ax.grid(axis="x", alpha=0.25)
        fig.text(0.02, 0.02, "Bars compare isoform-support expression deltas across two comparisons. This is support-pattern prioritization, not a causal isoform model.", fontsize=8.5, color="#4b5563")
        fig.tight_layout(rect=[0.04, 0.06, 0.98, 0.95])
        fig.savefig(path, dpi=180)
        plt.close(fig)
        return path

    @staticmethod
    def _candidate_gene_token_from_row(row: pd.Series) -> str:
        symbol = str(row.get("gene_symbol", "") or "").strip()
        if symbol:
            return f"sym:{symbol.casefold()}"
        gene_id = normalize_gene_id(row.get("gene_id"))
        if gene_id:
            return f"id:{gene_id}"
        return ""

    def _gene_token_set(self, frame: pd.DataFrame) -> set[str]:
        if frame is None or frame.empty:
            return set()
        tokens = {self._candidate_gene_token_from_row(row) for _, row in frame.iterrows()}
        return {item for item in tokens if item}

    def _gene_token_to_rows(self, frame: pd.DataFrame) -> dict[str, list[pd.Series]]:
        lookup: dict[str, list[pd.Series]] = {}
        if frame is None or frame.empty:
            return lookup
        for _, row in frame.iterrows():
            token = self._candidate_gene_token_from_row(row)
            if token:
                lookup.setdefault(token, []).append(row)
        return lookup

    def _load_rbp_gene_sets(self, custom_gene_set_path: str | None = None) -> pd.DataFrame:
        rows: list[dict[str, object]] = []
        builtin_symbols = [
            "Srsf1", "Srsf2", "Srsf3", "Srsf4", "Srsf5", "Srsf6", "Srsf7", "Srsf9", "Srsf10",
            "Hnrnpa1", "Hnrnpa2b1", "Hnrnpc", "Hnrnpf", "Hnrnph1", "Hnrnpk", "Hnrnpm", "Hnrnpu",
            "Rbfox1", "Rbfox2", "Rbfox3", "Mbnl1", "Mbnl2", "Mbnl3", "Celf1", "Celf2", "Ptbp1", "Ptbp2",
            "Qk", "Nova1", "Nova2", "Tardbp", "Fus", "Sfpq", "Tra2a", "Tra2b", "Elavl1", "Elavl2", "Elavl3", "Elavl4",
        ]
        for symbol in builtin_symbols:
            rows.append(
                {
                    "gene_set": "katana_builtin_splicing_rbp_seed",
                    "source": "built-in generic RBP/splicing-factor seed",
                    "gene_symbol": symbol,
                    "gene_id": "",
                }
            )
        path = Path(custom_gene_set_path) if custom_gene_set_path else None
        if path is not None and path.exists():
            try:
                sep = "," if path.suffix.lower() == ".csv" else "\t"
                custom = pd.read_csv(path, sep=sep, low_memory=False)
                gene_col = next((col for col in ["gene_symbol", "gene", "symbol", "Gene", "GeneSymbol"] if col in custom.columns), None)
                gene_id_col = next((col for col in ["gene_id", "GeneID", "ensembl_gene_id"] if col in custom.columns), None)
                set_col = next((col for col in ["gene_set", "set", "category", "RBP", "rbp_gene_set"] if col in custom.columns), None)
                source_col = next((col for col in ["source", "database"] if col in custom.columns), None)
                if gene_col or gene_id_col:
                    for _, row in custom.iterrows():
                        rows.append(
                            {
                                "gene_set": self._safe_text(row.get(set_col), "custom_rbp_gene_set") if set_col else "custom_rbp_gene_set",
                                "source": self._safe_text(row.get(source_col), str(path)) if source_col else str(path),
                                "gene_symbol": self._safe_text(row.get(gene_col), "") if gene_col else "",
                                "gene_id": self._safe_text(row.get(gene_id_col), "") if gene_id_col else "",
                            }
                        )
            except Exception as exc:
                rows.append(
                    {
                        "gene_set": "custom_rbp_gene_set",
                        "source": str(path),
                        "gene_symbol": "",
                        "gene_id": "",
                        "load_error": str(exc),
                    }
                )
        frame = pd.DataFrame(rows)
        if frame.empty:
            return pd.DataFrame(columns=["gene_set", "source", "gene_symbol", "gene_id", "gene_token"])
        frame["gene_token"] = frame.apply(
            lambda row: (
                f"sym:{str(row.get('gene_symbol', '')).strip().casefold()}"
                if str(row.get("gene_symbol", "") or "").strip()
                else f"id:{normalize_gene_id(row.get('gene_id'))}"
                if normalize_gene_id(row.get("gene_id"))
                else ""
            ),
            axis=1,
        )
        frame = frame.loc[frame["gene_token"].astype(str).ne("")].drop_duplicates(
            subset=["gene_set", "source", "gene_token"],
            keep="first",
        )
        return frame.reset_index(drop=True)

    @staticmethod
    def _hypergeom_enrichment_pvalue(*, universe_size: int, set_size: int, selected_size: int, overlap_size: int) -> float:
        if universe_size <= 0 or set_size <= 0 or selected_size <= 0 or overlap_size <= 0:
            return float("nan")
        try:
            from scipy.stats import hypergeom

            return float(hypergeom.sf(overlap_size - 1, universe_size, set_size, selected_size))
        except Exception:
            return float("nan")

    def _write_rbp_support_outputs(
        self,
        comparison_id: str,
        tables: dict[str, pd.DataFrame],
        source: str,
        top_n: int,
        custom_gene_set_path: str | None,
    ) -> None:
        base = self.rbp_support_output_dir(comparison_id)
        if base is None:
            return
        base.mkdir(parents=True, exist_ok=True)
        file_map = self.rbp_support_file_map(comparison_id)
        tables.get("summary", pd.DataFrame()).to_csv(file_map["summary"], sep="\t", index=False)
        tables.get("matched_genes", pd.DataFrame()).to_csv(file_map["matched_genes"], sep="\t", index=False)
        self._write_json(
            file_map["parameters"],
            {
                "module": "RBP support",
                "comparison_id": comparison_id,
                "candidate_source": source,
                "top_n": int(top_n),
                "custom_gene_set_path": custom_gene_set_path or "",
                "method": "Hypergeometric overlap/enrichment over the current candidate-ranking universe.",
                "interpretation": "Prioritization only; overlap does not prove causal RBP regulation.",
            },
        )

    def _write_rbp_cross_support_outputs(
        self,
        tables: dict[str, pd.DataFrame],
        source: str,
        top_n: int,
        custom_gene_set_path: str | None,
    ) -> None:
        base = self.rbp_cross_support_output_dir()
        if base is None:
            return
        base.mkdir(parents=True, exist_ok=True)
        file_map = self.rbp_cross_support_file_map()
        tables.get("enrichment_matrix", pd.DataFrame()).to_csv(file_map["enrichment_matrix"], sep="\t", index=False)
        tables.get("matched_gene_matrix", pd.DataFrame()).to_csv(file_map["matched_gene_matrix"], sep="\t", index=False)
        tables.get("comparison_summary", pd.DataFrame()).to_csv(file_map["comparison_summary"], sep="\t", index=False)
        generated_plot_files = self._write_rbp_cross_support_plots(tables)
        self._write_json(
            file_map["parameters"],
            {
                "module": "Cross-comparison RBP support",
                "candidate_source": source,
                "top_n": int(top_n),
                "custom_gene_set_path": custom_gene_set_path or "",
                "selected_comparisons": [comparison.comparison_id for comparison in self._selected_comparisons()],
                "method": "Compares per-comparison RBP/splicing-factor overlap and enrichment support patterns.",
                "interpretation": "Support-prioritization only; this does not prove causal RBP regulation.",
                "generated_plot_files": generated_plot_files,
            },
        )

    @staticmethod
    def _normalize_transcript_id(value: object) -> str:
        text = str(value or "").strip()
        if not text:
            return ""
        return re.sub(r"\.\d+$", "", text)

    def _gtf_transcript_gene_frame(self) -> pd.DataFrame:
        gtf_path = self._resolve_gtf_path()
        if gtf_path is None or not gtf_path.exists():
            return pd.DataFrame()
        cache_key = str(gtf_path.resolve())
        cached = self._gtf_transcript_gene_cache.get(cache_key)
        if cached is not None:
            return cached.copy()
        rows: list[dict[str, object]] = []
        opener = gzip.open if gtf_path.suffix.lower() == ".gz" else open
        try:
            with opener(gtf_path, "rt", encoding="utf-8", errors="replace") as handle:
                for line in handle:
                    if not line or line.startswith("#"):
                        continue
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) < 9 or parts[2] not in {"transcript", "exon"}:
                        continue
                    attrs = _parse_gtf_attrs(parts[8])
                    transcript_id = attrs.get("transcript_id") or attrs.get("transcript")
                    gene_id = attrs.get("gene_id")
                    if not transcript_id or not gene_id:
                        continue
                    rows.append(
                        {
                            "transcript_id": transcript_id,
                            "transcript_key": self._normalize_transcript_id(transcript_id),
                            "gene_id": normalize_gene_id(gene_id),
                            "gene_symbol": attrs.get("gene_name") or attrs.get("gene_symbol") or normalize_gene_id(gene_id),
                            "transcript_name": attrs.get("transcript_name", ""),
                            "transcript_type": attrs.get("transcript_type") or attrs.get("transcript_biotype") or attrs.get("gene_type") or "",
                        }
                    )
        except Exception:
            return pd.DataFrame()
        result = pd.DataFrame(rows)
        if not result.empty:
            result = result.drop_duplicates(subset=["transcript_key"], keep="first").reset_index(drop=True)
        self._gtf_transcript_gene_cache[cache_key] = result.copy()
        return result

    def _comparison_isoform_samples(self, comparison) -> list[IsoformSampleDefinition]:
        if self.current_project is None:
            return []
        group1, group2 = self._comparison_splicing_direction_labels(comparison)
        groups = {group1, group2}
        samples = [
            sample
            for sample in self.current_project.isoform_samples
            if sample.include
            and sample.quant_sf.exists()
            and ((sample.sample_group in groups if sample.sample_group else False) or (sample.condition in groups if sample.condition else False))
        ]
        if samples:
            return samples
        return [
            sample
            for sample in self.current_project.isoform_samples
            if sample.include and sample.quant_sf.exists() and sample.comparison_id == comparison.comparison_id
        ]

    def _sample_group_label(self, sample: IsoformSampleDefinition) -> str:
        return sample.sample_group or sample.condition or sample.comparison_id or sample.sample_id

    def _dtu_isoform_support_frame(self, comparison) -> pd.DataFrame:
        paths = [
            self._resolve_comparison_support_file(comparison, "dtu"),
            self._resolve_comparison_support_file(comparison, "stager"),
        ]
        rows: list[pd.DataFrame] = []
        for raw_path in paths:
            path = Path(raw_path) if raw_path else None
            if path is None or not path.exists():
                continue
            frame = self._read_table_or_empty(path)
            if frame.empty:
                continue
            working = frame.copy()
            q_col = next((col for col in ["padj", "qvalue", "q_value", "adj_pvalue", "adjusted_pvalue"] if col in working.columns), None)
            if q_col is None:
                continue
            tx_col = next((col for col in ["featureID", "transcript_id", "tx_id", "isoform_id", "transcript"] if col in working.columns), None)
            gene_col = next((col for col in ["groupID", "gene_id", "GeneID", "gene"] if col in working.columns), None)
            working["DTU_qvalue"] = pd.to_numeric(working[q_col], errors="coerce")
            working["gene_id"] = working[gene_col].map(normalize_gene_id) if gene_col else pd.NA
            working["transcript_id"] = working[tx_col].astype(str) if tx_col else ""
            working["transcript_key"] = working["transcript_id"].map(self._normalize_transcript_id)
            working["DTU_source_file"] = str(path)
            rows.append(working[["gene_id", "transcript_id", "transcript_key", "DTU_qvalue", "DTU_source_file"]].copy())
        if not rows:
            return pd.DataFrame(columns=["gene_id", "transcript_id", "transcript_key", "DTU_qvalue", "DTU_source_file"])
        result = pd.concat(rows, axis=0, ignore_index=True)
        result = result.sort_values("DTU_qvalue", na_position="last").drop_duplicates(
            subset=["gene_id", "transcript_key"],
            keep="first",
        )
        return result.reset_index(drop=True)

    def _isoform_transcript_columns(self) -> list[str]:
        return [
            "comparison_id",
            "display_name",
            "gene_symbol",
            "gene_id",
            "transcript_id",
            "transcript_name",
            "transcript_type",
            "group1",
            "group2",
            "mean_TPM_group1",
            "mean_TPM_group2",
            "delta_TPM_group1_minus_group2",
            "log2_TPM_ratio_group1_vs_group2",
            "mean_NumReads_group1",
            "mean_NumReads_group2",
            "n_group1_samples",
            "n_group2_samples",
            "replicate_TPM_group1",
            "replicate_TPM_group2",
            "DTU_qvalue",
            "DTU_significant",
            "DTU_source_file",
            "isoform_support_status",
        ]

    def _selected_isoform_gene_summary_without_quant(self, selected: pd.DataFrame, comparison, status: str) -> pd.DataFrame:
        context = self.comparison_context(comparison.comparison_id)
        group1, group2 = self._comparison_splicing_direction_labels(comparison)
        rows = []
        for _, row in selected.iterrows():
            rows.append(
                {
                    "comparison_id": comparison.comparison_id,
                    "display_name": context.get("display_name", comparison.display_resolved_name),
                    "gene_symbol": self._safe_text(row.get("gene_symbol"), self._safe_text(row.get("gene_id"))),
                    "gene_id": self._safe_text(row.get("gene_id"), ""),
                    "candidate_tier": self._safe_text(row.get("candidate_tier"), ""),
                    "rank": row.get("rank", pd.NA),
                    "DTU_significant": bool(row.get("DTU_significant", False)),
                    "best_DTU_qvalue": row.get("best_DTU_qvalue", pd.NA),
                    "n_DTU_significant_isoforms": row.get("n_DTU_significant_isoforms", 0),
                    "n_quantified_isoforms": 0,
                    "group1": group1,
                    "group2": group2,
                    "total_TPM_group1": pd.NA,
                    "total_TPM_group2": pd.NA,
                    "dominant_transcript_group1": "",
                    "dominant_transcript_group2": "",
                    "isoform_support_status": status,
                }
            )
        return pd.DataFrame(rows)

    def _build_isoform_transcript_support(
        self,
        selected: pd.DataFrame,
        comparison,
        samples: list[IsoformSampleDefinition],
        gtf_frame: pd.DataFrame,
        dtu_support: pd.DataFrame,
    ) -> pd.DataFrame:
        selected_tokens = self._gene_token_set(selected)
        selected_gene_ids = {
            normalize_gene_id(row.get("gene_id"))
            for _, row in selected.iterrows()
            if normalize_gene_id(row.get("gene_id"))
        }
        selected_symbols = {
            str(row.get("gene_symbol", "") or "").strip().casefold()
            for _, row in selected.iterrows()
            if str(row.get("gene_symbol", "") or "").strip()
        }
        quant_rows: list[pd.DataFrame] = []
        for sample in samples:
            try:
                frame = pd.read_csv(
                    sample.quant_sf,
                    sep="\t",
                    low_memory=False,
                    usecols=lambda col: col in {"Name", "TPM", "NumReads", "Length", "EffectiveLength"},
                )
            except Exception:
                continue
            if frame.empty or "Name" not in frame.columns:
                continue
            working = frame.copy()
            working["transcript_key"] = working["Name"].map(self._normalize_transcript_id)
            working["TPM"] = pd.to_numeric(working.get("TPM"), errors="coerce")
            working["NumReads"] = pd.to_numeric(working.get("NumReads"), errors="coerce")
            working["sample_id"] = sample.sample_id
            working["sample_group"] = self._sample_group_label(sample)
            quant_rows.append(working)
        if not quant_rows:
            return pd.DataFrame(columns=self._isoform_transcript_columns())
        quant = pd.concat(quant_rows, axis=0, ignore_index=True)
        quant = quant.merge(gtf_frame, on="transcript_key", how="left")
        quant["gene_key"] = quant.apply(lambda row: self._gene_lookup_key(row.get("gene_id"), row.get("gene_symbol")), axis=1)
        quant = quant.loc[
            quant["gene_key"].isin(selected_tokens)
            | quant["gene_id"].map(normalize_gene_id).isin(selected_gene_ids)
            | quant["gene_symbol"].fillna("").astype(str).str.casefold().isin(selected_symbols)
        ].copy()
        if quant.empty:
            return pd.DataFrame(columns=self._isoform_transcript_columns())
        group1, group2 = self._comparison_splicing_direction_labels(comparison)
        quant = quant.loc[quant["sample_group"].isin({group1, group2})].copy()
        if quant.empty:
            return pd.DataFrame(columns=self._isoform_transcript_columns())
        if dtu_support is not None and not dtu_support.empty:
            transcript_dtu = dtu_support.loc[dtu_support["transcript_key"].astype(str).ne("")].copy()
            transcript_dtu = transcript_dtu.sort_values("DTU_qvalue", na_position="last").drop_duplicates("transcript_key", keep="first")
            gene_dtu = dtu_support.loc[dtu_support["gene_id"].astype(str).ne("")].copy()
            gene_dtu = gene_dtu.sort_values("DTU_qvalue", na_position="last").drop_duplicates("gene_id", keep="first")
            if not transcript_dtu.empty:
                quant = quant.merge(
                    transcript_dtu[["transcript_key", "DTU_qvalue", "DTU_source_file"]],
                    on="transcript_key",
                    how="left",
                )
            else:
                quant["DTU_qvalue"] = pd.NA
                quant["DTU_source_file"] = ""
            if not gene_dtu.empty:
                quant = quant.merge(
                    gene_dtu[["gene_id", "DTU_qvalue", "DTU_source_file"]].rename(
                        columns={"DTU_qvalue": "DTU_gene_qvalue", "DTU_source_file": "DTU_gene_source_file"}
                    ),
                    on="gene_id",
                    how="left",
                )
                quant["DTU_qvalue"] = quant["DTU_qvalue"].combine_first(quant["DTU_gene_qvalue"])
                quant["DTU_source_file"] = quant["DTU_source_file"].combine_first(quant["DTU_gene_source_file"])
                quant = quant.drop(columns=["DTU_gene_qvalue", "DTU_gene_source_file"], errors="ignore")
        else:
            quant["DTU_qvalue"] = pd.NA
            quant["DTU_source_file"] = ""
        rows: list[dict[str, object]] = []
        group_cols = ["gene_id", "gene_symbol", "transcript_key", "transcript_id", "transcript_name", "transcript_type"]
        for keys, subset in quant.groupby(group_cols, dropna=False):
            gene_id, gene_symbol, transcript_key, transcript_id, transcript_name, transcript_type = keys
            by_group = subset.groupby("sample_group", dropna=False)
            g1 = by_group.get_group(group1) if group1 in by_group.groups else pd.DataFrame()
            g2 = by_group.get_group(group2) if group2 in by_group.groups else pd.DataFrame()
            mean_tpm_g1 = float(g1["TPM"].mean()) if not g1.empty else float("nan")
            mean_tpm_g2 = float(g2["TPM"].mean()) if not g2.empty else float("nan")
            mean_reads_g1 = float(g1["NumReads"].mean()) if not g1.empty else float("nan")
            mean_reads_g2 = float(g2["NumReads"].mean()) if not g2.empty else float("nan")
            qvalue = pd.to_numeric(subset.get("DTU_qvalue"), errors="coerce").min()
            rows.append(
                {
                    "comparison_id": comparison.comparison_id,
                    "display_name": comparison.display_resolved_name,
                    "gene_symbol": gene_symbol,
                    "gene_id": gene_id,
                    "transcript_id": transcript_id or transcript_key,
                    "transcript_name": transcript_name,
                    "transcript_type": transcript_type,
                    "group1": group1,
                    "group2": group2,
                    "mean_TPM_group1": mean_tpm_g1,
                    "mean_TPM_group2": mean_tpm_g2,
                    "delta_TPM_group1_minus_group2": mean_tpm_g1 - mean_tpm_g2 if np.isfinite(mean_tpm_g1) and np.isfinite(mean_tpm_g2) else float("nan"),
                    "log2_TPM_ratio_group1_vs_group2": math.log2((mean_tpm_g1 + 0.1) / (mean_tpm_g2 + 0.1)) if np.isfinite(mean_tpm_g1) and np.isfinite(mean_tpm_g2) else float("nan"),
                    "mean_NumReads_group1": mean_reads_g1,
                    "mean_NumReads_group2": mean_reads_g2,
                    "n_group1_samples": int(len(g1)),
                    "n_group2_samples": int(len(g2)),
                    "replicate_TPM_group1": "",
                    "replicate_TPM_group2": "",
                    "DTU_qvalue": qvalue,
                    "DTU_significant": bool(pd.notna(qvalue) and qvalue <= self.current_thresholds.dtu_qvalue),
                    "DTU_source_file": self._safe_text(subset.get("DTU_source_file", pd.Series([""])).dropna().head(1).iloc[0] if "DTU_source_file" in subset.columns and not subset["DTU_source_file"].dropna().empty else ""),
                    "isoform_support_status": "quantified isoform; DTU qvalue shown if available",
                }
            )
            rows[-1]["replicate_TPM_group1"] = "; ".join(
                f"{row.sample_id}:{float(row.TPM):.4g}" for row in g1[["sample_id", "TPM"]].itertuples(index=False) if pd.notna(row.TPM)
            )
            rows[-1]["replicate_TPM_group2"] = "; ".join(
                f"{row.sample_id}:{float(row.TPM):.4g}" for row in g2[["sample_id", "TPM"]].itertuples(index=False) if pd.notna(row.TPM)
            )
        result = pd.DataFrame(rows)
        if result.empty:
            return pd.DataFrame(columns=self._isoform_transcript_columns())
        return result.sort_values(["gene_symbol", "DTU_significant", "DTU_qvalue", "mean_TPM_group1"], ascending=[True, False, True, False], na_position="last").reset_index(drop=True)

    def _build_isoform_gene_summary(
        self,
        selected: pd.DataFrame,
        comparison,
        transcript_table: pd.DataFrame,
        selection_meta: dict[str, object],
    ) -> pd.DataFrame:
        if transcript_table.empty:
            return self._selected_isoform_gene_summary_without_quant(selected, comparison, "No quantified transcript matched selected candidate genes.")
        selected_lookup = {}
        for _, row in selected.iterrows():
            for key in {
                self._candidate_gene_token_from_row(row),
                self._gene_lookup_key(row.get("gene_id"), row.get("gene_symbol")),
                normalize_gene_id(row.get("gene_id")),
                str(row.get("gene_symbol", "") or "").strip().casefold(),
            }:
                if key:
                    selected_lookup[key] = row
        rows: list[dict[str, object]] = []
        for (gene_id, gene_symbol), subset in transcript_table.groupby(["gene_id", "gene_symbol"], dropna=False):
            lookup_keys = [
                self._gene_lookup_key(gene_id, gene_symbol),
                normalize_gene_id(gene_id),
                str(gene_symbol or "").strip().casefold(),
            ]
            candidate_row = pd.Series(dtype="object")
            for lookup_key in lookup_keys:
                if lookup_key and lookup_key in selected_lookup:
                    candidate_row = selected_lookup[lookup_key]
                    break
            group1 = self._safe_text(subset["group1"].dropna().head(1).iloc[0], "") if "group1" in subset.columns and not subset["group1"].dropna().empty else ""
            group2 = self._safe_text(subset["group2"].dropna().head(1).iloc[0], "") if "group2" in subset.columns and not subset["group2"].dropna().empty else ""
            total_g1 = pd.to_numeric(subset["mean_TPM_group1"], errors="coerce").sum()
            total_g2 = pd.to_numeric(subset["mean_TPM_group2"], errors="coerce").sum()
            dominant_g1 = subset.sort_values("mean_TPM_group1", ascending=False, na_position="last").head(1)
            dominant_g2 = subset.sort_values("mean_TPM_group2", ascending=False, na_position="last").head(1)
            best_q = pd.to_numeric(subset["DTU_qvalue"], errors="coerce").min()
            rows.append(
                {
                    "comparison_id": comparison.comparison_id,
                    "display_name": comparison.display_resolved_name,
                    "gene_symbol": gene_symbol,
                    "gene_id": gene_id,
                    "candidate_tier": self._safe_text(candidate_row.get("candidate_tier"), ""),
                    "rank": candidate_row.get("rank", pd.NA),
                    "evidence_class": self._safe_text(candidate_row.get("evidence_class"), ""),
                    "selection_rule": selection_meta.get("selection_rule", ""),
                    "n_quantified_isoforms": int(len(subset)),
                    "group1": group1,
                    "group2": group2,
                    "total_TPM_group1": float(total_g1),
                    "total_TPM_group2": float(total_g2),
                    "delta_total_TPM_group1_minus_group2": float(total_g1 - total_g2),
                    "dominant_transcript_group1": self._safe_text(dominant_g1.iloc[0].get("transcript_id"), "") if not dominant_g1.empty else "",
                    "dominant_transcript_group2": self._safe_text(dominant_g2.iloc[0].get("transcript_id"), "") if not dominant_g2.empty else "",
                    "best_DTU_qvalue": best_q,
                    "DTU_significant": bool(pd.notna(best_q) and best_q <= self.current_thresholds.dtu_qvalue),
                    "n_DTU_significant_isoforms": int(subset["DTU_significant"].map(self._safe_bool).sum()),
                    "isoform_support_status": "quant.sf support available; DTU qvalue shown if mapped",
                }
            )
        return pd.DataFrame(rows).sort_values(["DTU_significant", "best_DTU_qvalue", "rank"], ascending=[False, True, True], na_position="last").reset_index(drop=True)

    def _write_isoform_support_outputs(self, comparison_id: str, tables: dict[str, pd.DataFrame], source: str, top_n: int) -> None:
        base = self.isoform_support_output_dir(comparison_id)
        if base is None:
            return
        base.mkdir(parents=True, exist_ok=True)
        file_map = self.isoform_support_file_map(comparison_id)
        tables.get("gene_summary", pd.DataFrame()).to_csv(file_map["gene_summary"], sep="\t", index=False)
        tables.get("transcript_table", pd.DataFrame()).to_csv(file_map["transcript_table"], sep="\t", index=False)
        self._write_json(
            file_map["parameters"],
            {
                "module": "Isoform candidate support",
                "comparison_id": comparison_id,
                "candidate_source": source,
                "top_n": int(top_n),
                "inputs": "candidate ranking/selection + quant.sf + GTF + optional DTU qvalues",
                "interpretation": "Isoform usage support layer for candidate follow-up; DEG/rMATS remain the core candidate definition.",
            },
        )

    def _write_isoform_cross_support_outputs(self, tables: dict[str, pd.DataFrame], source: str, top_n: int) -> None:
        base = self.isoform_cross_support_output_dir()
        if base is None:
            return
        base.mkdir(parents=True, exist_ok=True)
        file_map = self.isoform_cross_support_file_map()
        tables.get("gene_support_matrix", pd.DataFrame()).to_csv(file_map["gene_support_matrix"], sep="\t", index=False)
        tables.get("transcript_support_matrix", pd.DataFrame()).to_csv(file_map["transcript_support_matrix"], sep="\t", index=False)
        tables.get("dominant_isoform_switches", pd.DataFrame()).to_csv(file_map["dominant_isoform_switches"], sep="\t", index=False)
        generated_plot_files = self._write_isoform_cross_support_plots(tables)
        self._write_json(
            file_map["parameters"],
            {
                "module": "Cross-comparison isoform support",
                "candidate_source": source,
                "top_n": int(top_n),
                "selected_comparisons": [comparison.comparison_id for comparison in self._selected_comparisons()],
                "inputs": "per-comparison candidate isoform support tables from quant.sf + GTF + optional DTU qvalues",
                "interpretation": "Compares isoform-support patterns and dominant transcript changes across comparisons.",
                "generated_plot_files": generated_plot_files,
            },
        )

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
        force_refresh: bool = False,
    ) -> pd.DataFrame:
        _ = allow_generate  # Explicitly unused: opening the page must not generate downstream inputs.
        return self._sashimi_event_catalog_frame(comparison_id, force_rebuild=force_refresh)

    def refresh_sashimi_event_catalog(self, comparison_id: str | None) -> pd.DataFrame:
        return self._sashimi_event_catalog_frame(comparison_id, force_rebuild=True)

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
                "rmats_sig_before_stability",
                "stable_as_event",
                "mean_junction_reads_group1",
                "mean_junction_reads_group2",
                "replicate_count_group1",
                "replicate_count_group2",
                "supported_replicate_count_group1",
                "supported_replicate_count_group2",
                "replicate_direction_consistency",
                "stable_support_status",
                "splicing_direction_flipped",
            ]
        )

    def _sashimi_event_catalog_file_map(self, comparison_id: str | None) -> dict[str, Path]:
        if self.current_project is None or self.current_project.output_root is None or not comparison_id:
            return {}
        base = self.current_project.output_root / "06_sashimi" / "event_catalog"
        safe_id = self._safe_output_name(str(comparison_id))
        return {
            "catalog": base / f"{safe_id}_sashimi_event_catalog.tsv",
            "parameters": base / f"{safe_id}_sashimi_event_catalog.parameters.json",
        }

    def _sashimi_event_catalog_parameters(self, comparison) -> dict[str, object]:
        return {
            "comparison_id": comparison.comparison_id,
            "rmats_mode": self.current_project.selected_rmats_mode if self.current_project else "JC",
            "rmats_name": comparison.rmats_name or "",
            "rmats_path": str(comparison.rmats_path or ""),
            "reverse_splicing": bool(comparison.reverse_splicing),
            "splicing_fdr": self.current_thresholds.splicing_fdr,
            "splicing_dpsi": self.current_thresholds.splicing_dpsi,
            "stable_min_junction_reads": self.current_thresholds.stable_min_junction_reads,
            "stable_min_replicates_per_group": self.current_thresholds.stable_min_replicates_per_group,
            "stable_min_pairwise_consistency": self.current_thresholds.stable_min_pairwise_consistency,
            "catalog_version": 3,
        }

    def _sashimi_event_catalog_frame(self, comparison_id: str | None, *, force_rebuild: bool = False) -> pd.DataFrame:
        if self.current_project is None or not comparison_id:
            return self._empty_sashimi_events_frame()
        comparison = next(
            (item for item in self.current_project.available_comparisons if item.comparison_id == comparison_id),
            None,
        )
        if comparison is None or comparison.rmats_path is None:
            return self._empty_sashimi_events_frame()
        mode = self.current_project.selected_rmats_mode or "JC"
        cache_key = "|".join(
            [
                comparison.comparison_id,
                str(mode),
                str(int(bool(comparison.reverse_splicing))),
                f"{float(self.current_thresholds.splicing_fdr):.6g}",
                f"{float(self.current_thresholds.splicing_dpsi):.6g}",
                f"{float(self.current_thresholds.stable_min_junction_reads):.6g}",
                str(int(self.current_thresholds.stable_min_replicates_per_group)),
                f"{float(self.current_thresholds.stable_min_pairwise_consistency):.6g}",
            ]
        )
        cached = None if force_rebuild else self._sashimi_event_catalog_cache.get(cache_key)
        if cached is not None:
            self._log_cache(f"sashimi event catalog cache hit: memory for {comparison.comparison_id}")
            return cached.copy()

        file_map = self._sashimi_event_catalog_file_map(comparison.comparison_id)
        expected_parameters = self._sashimi_event_catalog_parameters(comparison)
        catalog_path = file_map.get("catalog")
        if (
            not force_rebuild
            and catalog_path is not None
            and catalog_path.exists()
            and self._cached_parameters_match(file_map.get("parameters"), expected_parameters)
        ):
            try:
                disk_cached = pd.read_csv(catalog_path, sep="\t", low_memory=False)
            except Exception:
                disk_cached = pd.DataFrame()
            if list(disk_cached.columns):
                self._log_cache(f"sashimi event catalog cache hit: disk for {comparison.comparison_id}")
                self._sashimi_event_catalog_cache[cache_key] = disk_cached.copy()
                return disk_cached

        self._log_heavy(f"sashimi event catalog rebuild started for {comparison.comparison_id}")

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
            working["normalized_gene_id"] = normalized_gene_ids

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
            working["splicing_direction_flipped"] = bool(comparison.reverse_splicing)
            working["FDR"] = pd.to_numeric(working["FDR"] if "FDR" in working.columns else working.get("fdr"), errors="coerce")
            working = self._apply_stable_as_columns(working, dpsi_column="dPSI", reverse_splicing=bool(comparison.reverse_splicing))
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
                        "normalized_gene_id",
                        "event_type",
                        "event_id",
                        "event_matching_key",
                        "dPSI",
                        "FDR",
                        "direction",
                        "coordinates",
                        "rmats_sig_before_stability",
                        "stable_as_event",
                        "mean_junction_reads_group1",
                        "mean_junction_reads_group2",
                        "replicate_count_group1",
                        "replicate_count_group2",
                        "supported_replicate_count_group1",
                        "supported_replicate_count_group2",
                        "replicate_direction_consistency",
                        "stable_support_status",
                        "splicing_direction_flipped",
                    ]
                ].copy()
            )

        if not rows:
            empty = self._empty_sashimi_events_frame()
            self._sashimi_event_catalog_cache[cache_key] = empty
            self._write_sashimi_event_catalog_cache(empty, file_map, expected_parameters)
            return empty.copy()

        combined = pd.concat(rows, ignore_index=True)
        combined = combined.drop_duplicates(subset=["comparison_id", "event_matching_key"]).reset_index(drop=True)
        combined = combined.sort_values(
            by=["FDR", "event_type", "gene_symbol", "event_id"],
            ascending=[True, True, True, True],
            na_position="last",
        ).reset_index(drop=True)
        self._sashimi_event_catalog_cache[cache_key] = combined.copy()
        self._write_sashimi_event_catalog_cache(combined, file_map, expected_parameters)
        return combined

    def _write_sashimi_event_catalog_cache(
        self,
        frame: pd.DataFrame,
        file_map: dict[str, Path],
        parameters: dict[str, object],
    ) -> None:
        catalog_path = file_map.get("catalog")
        parameters_path = file_map.get("parameters")
        if catalog_path is None or parameters_path is None:
            return
        catalog_path.parent.mkdir(parents=True, exist_ok=True)
        frame.to_csv(catalog_path, sep="\t", index=False)
        self._write_json(parameters_path, parameters)

    def _significant_as_event_summary_frame(self, comparison_id: str | None) -> pd.DataFrame:
        cache_key = "|".join(
            [
                str(comparison_id or ""),
                f"{float(self.current_thresholds.splicing_fdr):.6g}",
                f"{float(self.current_thresholds.splicing_dpsi):.6g}",
                f"{float(self.current_thresholds.stable_min_junction_reads):.6g}",
                str(int(self.current_thresholds.stable_min_replicates_per_group)),
                f"{float(self.current_thresholds.stable_min_pairwise_consistency):.6g}",
            ]
        )
        cached = self._significant_as_event_summary_cache.get(cache_key)
        if cached is not None:
            return cached.copy()
        events = self._sashimi_event_catalog_frame(comparison_id)
        if events.empty:
            empty = pd.DataFrame(
                columns=[
                    "comparison_id",
                    "normalized_gene_id",
                    "significant_as_event_ids",
                    "significant_as_event_list",
                    "significant_as_stability_summary",
                    "stable_as_min_supported_replicates_group1",
                    "stable_as_min_supported_replicates_group2",
                    "stable_as_min_direction_consistency",
                ]
            )
            self._significant_as_event_summary_cache[cache_key] = empty.copy()
            return empty
        working = events.copy()
        working["FDR"] = pd.to_numeric(working.get("FDR"), errors="coerce")
        working["dPSI"] = pd.to_numeric(working.get("dPSI"), errors="coerce")
        if "stable_as_event" in working.columns:
            working["significant_event"] = working["stable_as_event"].fillna(False).map(self._safe_bool)
        else:
            working["significant_event"] = (
                working["FDR"].le(self.current_thresholds.splicing_fdr)
                & working["dPSI"].abs().ge(self.current_thresholds.splicing_dpsi)
            )
        working = working.loc[working["significant_event"].fillna(False)].copy()
        if working.empty:
            empty = pd.DataFrame(
                columns=[
                    "comparison_id",
                    "normalized_gene_id",
                    "significant_as_event_ids",
                    "significant_as_event_list",
                    "significant_as_stability_summary",
                    "stable_as_min_supported_replicates_group1",
                    "stable_as_min_supported_replicates_group2",
                    "stable_as_min_direction_consistency",
                ]
            )
            self._significant_as_event_summary_cache[cache_key] = empty.copy()
            return empty
        working["event_label"] = working.apply(
            lambda row: f"{row.get('event_type', '')}:{row.get('event_id', '')}",
            axis=1,
        )
        rows: list[dict[str, object]] = []
        for (comp_id, gene_id), subset in working.groupby(["comparison_id", "normalized_gene_id"], dropna=False):
            event_ids = []
            event_labels = []
            stability_labels = []
            seen_ids: set[str] = set()
            seen_labels: set[str] = set()
            ordered_subset = subset.sort_values(
                ["FDR", "event_type", "event_id"],
                ascending=[True, True, True],
                na_position="last",
            )
            for _, row in ordered_subset.iterrows():
                event_id = self._safe_text(row.get("event_id", ""))
                event_label = self._safe_text(row.get("event_label", ""))
                if event_id and event_id not in seen_ids:
                    seen_ids.add(event_id)
                    event_ids.append(event_id)
                if event_label and event_label not in seen_labels:
                    seen_labels.add(event_label)
                    event_labels.append(event_label)
                    supported_1 = pd.to_numeric(row.get("supported_replicate_count_group1"), errors="coerce")
                    supported_2 = pd.to_numeric(row.get("supported_replicate_count_group2"), errors="coerce")
                    total_1 = pd.to_numeric(row.get("replicate_count_group1"), errors="coerce")
                    total_2 = pd.to_numeric(row.get("replicate_count_group2"), errors="coerce")
                    reads_1 = pd.to_numeric(row.get("mean_junction_reads_group1"), errors="coerce")
                    reads_2 = pd.to_numeric(row.get("mean_junction_reads_group2"), errors="coerce")
                    consistency = pd.to_numeric(row.get("replicate_direction_consistency"), errors="coerce")
                    stability_labels.append(
                        f"{event_label} | supported replicates "
                        f"{int(supported_1) if pd.notna(supported_1) else 'NA'}/{int(total_1) if pd.notna(total_1) else 'NA'} vs "
                        f"{int(supported_2) if pd.notna(supported_2) else 'NA'}/{int(total_2) if pd.notna(total_2) else 'NA'} | "
                        f"mean IJC+SJC {float(reads_1):.1f}/{float(reads_2):.1f} | "
                        f"direction consistency {float(consistency):.2f}"
                        if pd.notna(reads_1) and pd.notna(reads_2) and pd.notna(consistency)
                        else f"{event_label} | replicate stability metrics unavailable"
                    )
            supported_group1 = pd.to_numeric(
                ordered_subset.get("supported_replicate_count_group1"), errors="coerce"
            )
            supported_group2 = pd.to_numeric(
                ordered_subset.get("supported_replicate_count_group2"), errors="coerce"
            )
            consistency_values = pd.to_numeric(
                ordered_subset.get("replicate_direction_consistency"), errors="coerce"
            )
            rows.append(
                {
                    "comparison_id": self._safe_text(comp_id),
                    "normalized_gene_id": self._safe_text(gene_id),
                    "significant_as_event_ids": "; ".join(event_ids),
                    "significant_as_event_list": "; ".join(event_labels),
                    "significant_as_stability_summary": "; ".join(stability_labels),
                    "stable_as_min_supported_replicates_group1": supported_group1.min(skipna=True),
                    "stable_as_min_supported_replicates_group2": supported_group2.min(skipna=True),
                    "stable_as_min_direction_consistency": consistency_values.min(skipna=True),
                }
            )
        result = pd.DataFrame(rows)
        self._significant_as_event_summary_cache[cache_key] = result.copy()
        return result

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
        manifest = self._enrich_sashimi_manifest_with_event_metadata(manifest, source)
        self.run_state.results.sashimi_manifest = manifest.copy()
        selected_comparison_id = str(source["comparison_id"].iloc[0]) if "comparison_id" in source.columns else str(comparison_id or "")
        self.publish_sashimi_event_selection(
            source,
            comparison_id=selected_comparison_id,
            selection_rule="User-selected rMATS events from the 5.5 event catalog.",
        )
        self._mark_module_state(
            "sashimi_manifest",
            "finished" if not manifest.empty else "failed",
            "Generated sashimi manifest for the selected events." if not manifest.empty else "No sashimi manifest rows could be generated for the selected events.",
            output_folder=str(self.current_project.output_root / "06_sashimi") if self.current_project.output_root else "",
            from_cache=False,
        )
        return manifest

    def build_custom_native_sashimi_manifest(self, comparison_id: str | None, custom_event: dict[str, object]) -> pd.DataFrame:
        """Build a manifest row for a user-defined SE-like skipped-block event.

        Coordinates are generic 1-based inclusive genomic intervals. This is
        intentionally separate from rMATS event discovery so custom junctions
        that are absent from rMATS tables can still be visualized.
        """
        if self.current_project is None:
            raise RuntimeError("No project loaded.")
        comparison = next(
            (item for item in self.current_project.available_comparisons if item.comparison_id == comparison_id),
            None,
        )
        if comparison is None:
            raise RuntimeError("Select a valid comparison before building a custom sashimi event.")

        event_type = str(custom_event.get("event_type", "SE") or "SE").strip().upper()
        if event_type != "SE":
            raise RuntimeError("Custom native events currently support SE-like skipped-block coordinates.")

        def text_value(key: str, default: str = "") -> str:
            value = custom_event.get(key, default)
            if value is None or (isinstance(value, float) and pd.isna(value)):
                return default
            text = str(value).strip()
            return text if text else default

        def int_value(key: str) -> int:
            raw = text_value(key)
            if not raw:
                raise RuntimeError(f"Missing custom event coordinate: {key}")
            return int(float(raw.replace(",", "")))

        gene_symbol = text_value("gene_symbol")
        gene_id = text_value("gene_id")
        chrom = text_value("chrom")
        strand = text_value("strand", "+")
        if not gene_symbol and not gene_id:
            raise RuntimeError("Custom event requires gene symbol or gene ID.")
        if not chrom:
            raise RuntimeError("Custom event requires chromosome.")
        if strand not in {"+", "-"}:
            raise RuntimeError("Custom event strand must be + or -.")

        upstream_start = int_value("upstream_start")
        upstream_end = int_value("upstream_end")
        event_start = int_value("event_start")
        event_end = int_value("event_end")
        downstream_start = int_value("downstream_start")
        downstream_end = int_value("downstream_end")
        if min(upstream_start, upstream_end, event_start, event_end, downstream_start, downstream_end) <= 0:
            raise RuntimeError("Custom event coordinates must be positive 1-based genomic positions.")

        event_id = text_value("event_id") or self._safe_output_name(
            f"custom_{gene_symbol or gene_id}_{chrom}_{event_start}_{event_end}_{event_type}"
        )
        output_root = self.require_output_root()
        event_root = output_root / "06_sashimi" / "custom_events" / self._safe_output_name(comparison.comparison_id)
        bam_root = output_root / "06_sashimi" / "generated_bam_lists"
        event_root.mkdir(parents=True, exist_ok=True)
        bam_root.mkdir(parents=True, exist_ok=True)

        event_row = {
            "ID": event_id,
            "GeneID": gene_id,
            "geneSymbol": gene_symbol,
            "chr": chrom,
            "strand": strand,
            "exonStart_0base": min(event_start, event_end) - 1,
            "exonEnd": max(event_start, event_end),
            "upstreamES": min(upstream_start, upstream_end) - 1,
            "upstreamEE": max(upstream_start, upstream_end),
            "downstreamES": min(downstream_start, downstream_end) - 1,
            "downstreamEE": max(downstream_start, downstream_end),
            "IncLevelDifference": text_value("dpsi"),
            "FDR": text_value("fdr"),
            "IncLevel1": text_value("psi1"),
            "IncLevel2": text_value("psi2"),
        }
        event_file = event_root / f"{self._safe_output_name(event_id)[:140]}.SE.tsv"
        pd.DataFrame([event_row]).to_csv(event_file, sep="\t", index=False)

        bam_inputs = self._prepare_generated_sashimi_bam_lists(comparison, bam_root)
        if bam_inputs is None:
            raise RuntimeError("Could not resolve BAM lists for this comparison. Check sample group names and BAM inputs.")
        b1_txt, b2_txt, label1, label2 = bam_inputs
        display_group = comparison.display_resolved_name
        outdir = Path(self._safe_output_name(display_group)) / self._safe_output_name(
            f"{display_group}_{gene_symbol or gene_id}_{event_type}_{event_id}"
        )[:120]
        manifest = pd.DataFrame(
            [
                {
                    "display_group": display_group,
                    "comparison_id": comparison.comparison_id,
                    "comparison_name": display_group,
                    "geneSymbol": gene_symbol,
                    "GeneID": gene_id,
                    "event_uid": event_id,
                    "event_type": event_type,
                    "class_interpretation": "user_defined_custom_event",
                    "dPSI": text_value("dpsi"),
                    "FDR": text_value("fdr"),
                    "direction": "custom event uses the selected comparison canonical AS direction",
                    "coordinates": (
                        f"{event_type} | chr={chrom} | strand={strand} | "
                        f"upstream={min(upstream_start, upstream_end)}-{max(upstream_start, upstream_end)} | "
                        f"event={min(event_start, event_end)}-{max(event_start, event_end)} | "
                        f"downstream={min(downstream_start, downstream_end)}-{max(downstream_start, downstream_end)}"
                    ),
                    "event_matching_key": f"{gene_symbol or gene_id}|{event_type}|{event_id}",
                    "splicing_direction_flipped": False,
                    "label1": label1,
                    "label2": label2,
                    "b1_txt": str(b1_txt),
                    "b2_txt": str(b2_txt),
                    "event_file": str(event_file),
                    "source_rmats_file": "",
                    "outdir": str(outdir),
                }
            ]
        )
        self.run_state.results.sashimi_manifest = manifest.copy()
        self.publish_sashimi_event_selection(
            manifest,
            comparison_id=comparison.comparison_id,
            selection_rule="User-defined native sashimi event coordinates.",
        )
        self._mark_module_state(
            "sashimi_manifest",
            "finished",
            "Generated custom native sashimi manifest for the selected comparison.",
            output_folder=str(output_root / "06_sashimi"),
            from_cache=False,
        )
        return manifest

    def build_custom_native_sashimi_manifest_for_rows(self, comparison_id: str | None, custom_events: pd.DataFrame) -> pd.DataFrame:
        if custom_events is None or custom_events.empty:
            raise RuntimeError("No custom events selected.")
        frames: list[pd.DataFrame] = []
        for _, row in custom_events.iterrows():
            frames.append(self.build_custom_native_sashimi_manifest(comparison_id, row.to_dict()))
        manifest = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
        self.run_state.results.sashimi_manifest = manifest.copy()
        self._mark_module_state(
            "sashimi_manifest",
            "finished" if not manifest.empty else "failed",
            f"Generated custom native sashimi manifest for {len(manifest)} selected event(s)."
            if not manifest.empty
            else "No custom native sashimi manifest rows could be generated.",
            output_folder=str((self.current_project.output_root if self.current_project else Path(".")) / "06_sashimi")
            if self.current_project and self.current_project.output_root
            else "",
            from_cache=False,
        )
        return manifest

    def detect_significant_continuous_skipped_exon_junctions(
        self,
        comparison_id: str | None,
        *,
        gene_query: str,
        fdr_cutoff: float = 0.05,
        abs_dpsi_cutoff: float = 0.10,
    ) -> pd.DataFrame:
        """Detect rMATS-significant SE events and annotate skipped exon labels.

        This is a targeted gene-level detector. It avoids project-load or
        page-load scanning and only runs after the user clicks Detect. Significance
        comes only from rMATS FDR and rMATS dPSI; BAMs are not used here.
        """
        if self.current_project is None:
            raise RuntimeError("No project loaded.")
        query = str(gene_query or "").strip()
        if not query:
            raise RuntimeError("Enter a gene symbol or gene ID before detecting continuous skipped-exon junctions.")
        comparison = next(
            (item for item in self.current_project.available_comparisons if item.comparison_id == comparison_id),
            None,
        )
        if comparison is None:
            raise RuntimeError("Select a valid comparison before detecting continuous skipped-exon junctions.")
        gtf_path = self._resolve_gtf_path()
        if gtf_path is None or not gtf_path.exists():
            raise RuntimeError("Genome GTF is required for continuous skipped-exon detection.")
        return self._detect_rmats_significant_continuous_skips(
            comparison,
            gtf_path,
            query,
            fdr_cutoff=float(fdr_cutoff),
            abs_dpsi_cutoff=float(abs_dpsi_cutoff),
        )

    def _detect_rmats_significant_continuous_skips(
        self,
        comparison,
        gtf_path: Path,
        gene_query: str,
        *,
        fdr_cutoff: float,
        abs_dpsi_cutoff: float,
    ) -> pd.DataFrame:
        mode = self.current_project.selected_rmats_mode if self.current_project else "JC"
        source_path = self._resolve_rmats_event_table_path(comparison, "SE", mode or "JC")
        if source_path is None or not source_path.exists():
            return self._empty_continuous_skip_frame()
        frame = pd.read_csv(source_path, sep="\t", low_memory=False)
        if frame.empty:
            return self._empty_continuous_skip_frame()
        working = frame.copy()
        if "GeneID" not in working.columns:
            working["GeneID"] = ""
        if "geneSymbol" not in working.columns:
            working["geneSymbol"] = ""
        query = str(gene_query or "").strip()
        query_core = normalize_gene_id(query)
        gene_mask = (
            working["geneSymbol"].fillna("").astype(str).str.strip().eq(query)
            | working["GeneID"].fillna("").astype(str).str.strip().eq(query)
            | working["GeneID"].fillna("").astype(str).map(normalize_gene_id).eq(query_core)
        )
        working = working.loc[gene_mask].copy()
        if working.empty:
            return self._empty_continuous_skip_frame()
        working["FDR_numeric"] = pd.to_numeric(working.get("FDR"), errors="coerce")
        dpsi = pd.to_numeric(working.get("IncLevelDifference"), errors="coerce")
        if comparison.reverse_splicing:
            dpsi = -dpsi
        working["standardized_dPSI"] = dpsi
        working = working.loc[
            working["FDR_numeric"].le(fdr_cutoff) & working["standardized_dPSI"].abs().ge(abs_dpsi_cutoff)
        ].copy()
        if working.empty:
            return self._empty_continuous_skip_frame()

        transcript_cache: dict[str, list[dict[str, object]]] = {}
        rows: list[dict[str, object]] = []
        for _, row in working.iterrows():
            gene_symbol = str(row.get("geneSymbol", "")).strip().strip('"')
            gene_id = str(row.get("GeneID", "")).strip().strip('"')
            cache_key = gene_id or gene_symbol
            if cache_key not in transcript_cache:
                transcript_cache[cache_key] = self._gtf_transcript_exons_for_gene(gtf_path, gene_id or gene_symbol)
            transcripts = transcript_cache[cache_key]
            event_start = int(float(row.get("exonStart_0base"))) + 1
            event_end = int(float(row.get("exonEnd")))
            chrom = str(row.get("chr", "")).strip()
            strand = str(row.get("strand", "")).strip()
            skipped_exons, transcript_id, transcript_name, n_skipped = self._label_event_block_from_transcripts(
                transcripts,
                chrom,
                strand,
                event_start,
                event_end,
            )
            if not skipped_exons:
                continue
            upstream_start = int(float(row.get("upstreamES"))) + 1
            upstream_end = int(float(row.get("upstreamEE")))
            downstream_start = int(float(row.get("downstreamES"))) + 1
            downstream_end = int(float(row.get("downstreamEE")))
            event_id = event_uid_from_row(row, "SE")
            psi1 = str(row.get("IncLevel1", ""))
            psi2 = str(row.get("IncLevel2", ""))
            if comparison.reverse_splicing:
                psi1, psi2 = psi2, psi1
            dpsi_value = float(row.get("standardized_dPSI"))
            fdr_value = float(row.get("FDR_numeric"))
            rows.append(
                {
                    "source": "continuous_skip_detector",
                    "comparison_id": comparison.comparison_id,
                    "comparison_display_name": comparison.display_resolved_name,
                    "gene_symbol": gene_symbol,
                    "gene_id": gene_id,
                    "event_type": "SE",
                    "event_id": event_id,
                    "event_matching_key": f"{gene_id}|{transcript_id}|{event_id}",
                    "dPSI": dpsi_value,
                    "FDR": fdr_value,
                    "qvalue": fdr_value,
                    "direction": "rMATS significant SE; dPSI follows canonical AS direction",
                    "coordinates": (
                        f"rMATS significant continuous skipped exons {skipped_exons}; "
                        f"{chrom}:{event_start}-{event_end}"
                    ),
                    "splicing_direction_flipped": bool(comparison.reverse_splicing),
                    "transcript_id": transcript_id,
                    "transcript_name": transcript_name,
                    "chrom": chrom,
                    "strand": strand,
                    "skipped_exons": skipped_exons,
                    "n_skipped_exons": n_skipped,
                    "PSI_group1": np.nanmean([float(item) for item in str(psi1).split(",") if item not in {"", "NA", "nan"}])
                    if str(psi1).strip()
                    else np.nan,
                    "PSI_group2": np.nanmean([float(item) for item in str(psi2).split(",") if item not in {"", "NA", "nan"}])
                    if str(psi2).strip()
                    else np.nan,
                    "psi1": psi1,
                    "psi2": psi2,
                    "fdr": fdr_value,
                    "dpsi": dpsi_value,
                    "upstream_start": upstream_start,
                    "upstream_end": upstream_end,
                    "event_start": event_start,
                    "event_end": event_end,
                    "downstream_start": downstream_start,
                    "downstream_end": downstream_end,
                }
            )
        result = pd.DataFrame(rows)
        if result.empty:
            return self._empty_continuous_skip_frame()
        output_root = self.require_output_root()
        outdir = output_root / "06_sashimi" / "continuous_skip_detector" / self._safe_output_name(comparison.comparison_id)
        outdir.mkdir(parents=True, exist_ok=True)
        result.to_csv(outdir / "significant_continuous_skipped_exon_junctions.tsv", sep="\t", index=False)
        return result.sort_values(["FDR", "n_skipped_exons", "gene_symbol"], ascending=[True, False, True]).reset_index(drop=True)

    def _label_event_block_from_transcripts(
        self,
        transcripts: list[dict[str, object]],
        chrom: str,
        strand: str,
        event_start: int,
        event_end: int,
    ) -> tuple[str, str, str, int]:
        best: tuple[int, str, str, str] = (0, "", "", "")
        for transcript in transcripts:
            if str(transcript.get("chrom", "")) != chrom:
                continue
            if strand and str(transcript.get("strand", "")) != strand:
                continue
            labels: list[str] = []
            for exon in transcript.get("exons", []):
                exon_start = int(exon["start"])
                exon_end = int(exon["end"])
                overlap = max(0, min(event_end, exon_end) - max(event_start, exon_start) + 1)
                if overlap <= 0:
                    continue
                exon_fraction = overlap / max(1, exon_end - exon_start + 1)
                if exon_fraction >= 0.50:
                    labels.append(str(exon["label"]))
            formatted = self._format_detected_exon_range(labels)
            score = len(set(labels))
            if formatted and score > best[0]:
                best = (
                    score,
                    formatted,
                    str(transcript.get("transcript_id", "")),
                    str(transcript.get("transcript_name", "")),
                )
        return best[1], best[2], best[3], best[0]

    @staticmethod
    def _empty_continuous_skip_frame() -> pd.DataFrame:
        return pd.DataFrame(
            columns=[
                "source",
                "comparison_id",
                "comparison_display_name",
                "gene_symbol",
                "gene_id",
                "event_type",
                "event_id",
                "dPSI",
                "FDR",
                "qvalue",
                "direction",
                "coordinates",
                "event_matching_key",
                "splicing_direction_flipped",
                "transcript_id",
                "transcript_name",
                "chrom",
                "strand",
                "skipped_exons",
                "n_skipped_exons",
                "PSI_group1",
                "PSI_group2",
                "psi1",
                "psi2",
                "fdr",
                "dpsi",
                "upstream_start",
                "upstream_end",
                "event_start",
                "event_end",
                "downstream_start",
                "downstream_end",
            ]
        )

    def _gtf_transcript_exons_for_gene(self, gtf_path: Path, gene_query: str) -> list[dict[str, object]]:
        opener = gzip.open if gtf_path.suffix == ".gz" else open
        query = gene_query.strip()
        query_core = normalize_gene_id(query)
        transcripts: dict[str, dict[str, object]] = {}
        with opener(gtf_path, "rt", encoding="utf-8", errors="ignore") as handle:
            for line in handle:
                if line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 9 or parts[2] != "exon":
                    continue
                attrs = _parse_gtf_attrs(parts[8])
                gene_symbol = str(attrs.get("gene_name", "")).strip()
                gene_id = str(attrs.get("gene_id", "")).strip()
                if query and query not in {gene_symbol, gene_id, normalize_gene_id(gene_id), query_core}:
                    continue
                transcript_id = str(attrs.get("transcript_id", "") or attrs.get("transcript_name", "") or "unknown_transcript")
                transcript = transcripts.setdefault(
                    transcript_id,
                    {
                        "gene_symbol": gene_symbol or query,
                        "gene_id": gene_id,
                        "transcript_id": transcript_id,
                        "transcript_name": str(attrs.get("transcript_name", transcript_id)),
                        "chrom": parts[0],
                        "strand": parts[6],
                        "exons": [],
                    },
                )
                exon_number = str(attrs.get("exon_number", "")).strip()
                label = f"E{exon_number}" if exon_number else str(attrs.get("exon_id", "exon"))
                transcript["exons"].append({"start": int(parts[3]), "end": int(parts[4]), "label": label})
        return [item for item in transcripts.values() if len(item.get("exons", [])) >= 3]

    @staticmethod
    def _format_detected_exon_range(labels: list[str]) -> str:
        numbers: list[int] = []
        for label in labels:
            match = re.fullmatch(r"E(\d+)", str(label))
            if match:
                numbers.append(int(match.group(1)))
        if numbers and len(numbers) == len(labels):
            numbers = sorted(set(numbers))
            if len(numbers) == 1:
                return f"E{numbers[0]}"
            if numbers == list(range(numbers[0], numbers[-1] + 1)):
                return f"E{numbers[0]}-E{numbers[-1]}"
        return "/".join(str(label) for label in labels)

    @staticmethod
    def _enrich_sashimi_manifest_with_event_metadata(manifest: pd.DataFrame, events: pd.DataFrame) -> pd.DataFrame:
        if manifest.empty or events.empty:
            return manifest
        enriched = manifest.copy()
        event_key = "event_uid" if "event_uid" in enriched.columns else "event_id" if "event_id" in enriched.columns else None
        if event_key is None or "event_id" not in events.columns:
            return enriched
        metadata_columns = [
            "event_id",
            "dPSI",
            "FDR",
            "direction",
            "coordinates",
            "event_matching_key",
            "comparison_display_name",
            "splicing_direction_flipped",
        ]
        available = [column for column in metadata_columns if column in events.columns]
        if not available:
            return enriched
        metadata = events[available].drop_duplicates(subset=["event_id"], keep="first").copy()
        enriched = enriched.merge(metadata, how="left", left_on=event_key, right_on="event_id", suffixes=("", "_selected"))
        if "event_id_selected" in enriched.columns:
            enriched = enriched.drop(columns=["event_id_selected"])
        return enriched

    def last_sashimi_failures(self) -> pd.DataFrame:
        return self._last_sashimi_failures.copy()

    def build_sashimi_output_browser(self, comparison_id: str | None = None) -> pd.DataFrame:
        if self.current_project is None:
            raise RuntimeError("No project loaded.")
        output_root = self.require_output_root()
        base = output_root / "06_sashimi" / "plots"
        columns = [
            "gene",
            "event_type",
            "event_id",
            "comparison_display_name",
            "dPSI",
            "FDR",
            "comparison_id",
            "relative_path",
            "absolute_path",
            "file_type",
            "file_size",
            "last_modified",
            "preview_kind",
            "previewable",
        ]
        if not base.exists():
            return pd.DataFrame(columns=columns)

        safe_display_to_comparison = {
            self._safe_output_name(comparison.display_resolved_name): comparison.comparison_id
            for comparison in self.current_project.available_comparisons
        }
        manifest_meta = self._sashimi_manifest_output_metadata()
        rows: list[dict[str, object]] = []
        for path in sorted(base.rglob("*")):
            if not path.is_file():
                continue
            rel = path.relative_to(base)
            inferred_comparison = ""
            if rel.parts:
                inferred_comparison = safe_display_to_comparison.get(rel.parts[0], "")
            metadata = self._metadata_for_sashimi_output(rel, manifest_meta)
            if not metadata:
                metadata = self._metadata_from_sashimi_output_path(rel)
            if not inferred_comparison:
                inferred_comparison = str(metadata.get("comparison_id", ""))
            if comparison_id and inferred_comparison and inferred_comparison != comparison_id:
                continue
            stat = path.stat()
            preview_kind = self._preview_kind_for_path(path)
            rows.append(
                {
                    "gene": metadata.get("gene", ""),
                    "event_type": metadata.get("event_type", ""),
                    "event_id": metadata.get("event_id", ""),
                    "comparison_display_name": metadata.get("comparison_display_name", ""),
                    "dPSI": metadata.get("dPSI", ""),
                    "FDR": metadata.get("FDR", ""),
                    "comparison_id": inferred_comparison,
                    "relative_path": str(rel),
                    "absolute_path": str(path),
                    "file_type": path.suffix.lower().lstrip(".") or "file",
                    "file_size": stat.st_size,
                    "last_modified": datetime.fromtimestamp(stat.st_mtime).isoformat(timespec="seconds"),
                    "preview_kind": preview_kind,
                    "previewable": preview_kind in {"image", "table", "excel", "html"},
                }
            )
        frame = pd.DataFrame(rows, columns=columns)
        if frame.empty:
            return frame
        if comparison_id:
            exact = frame.loc[frame["comparison_id"].astype(str) == str(comparison_id)].copy()
            if not exact.empty:
                frame = exact
        return frame.sort_values(["comparison_id", "gene", "event_type", "relative_path"], na_position="last").reset_index(drop=True)

    def _sashimi_manifest_output_metadata(self) -> list[tuple[str, dict[str, object]]]:
        manifest = self.run_state.results.sashimi_manifest
        if manifest is None or manifest.empty or "outdir" not in manifest.columns:
            return []
        rows: list[tuple[str, dict[str, object]]] = []
        for _, row in manifest.iterrows():
            outdir = self._safe_text(row.get("outdir", "")).replace("\\", "/").strip("/")
            if not outdir:
                continue
            metadata = {
                "gene": (
                    self._safe_text(row.get("geneSymbol"))
                    or self._safe_text(row.get("gene_symbol"))
                    or self._safe_text(row.get("GeneID"))
                    or self._safe_text(row.get("gene_id"))
                ),
                "event_type": self._safe_text(row.get("event_type", "")),
                "event_id": self._safe_text(row.get("event_uid")) or self._safe_text(row.get("event_id")),
                "comparison_id": self._safe_text(row.get("comparison_id", "")),
                "comparison_display_name": self._safe_text(row.get("comparison_name")) or self._safe_text(row.get("display_group")),
                "dPSI": row.get("dPSI", ""),
                "FDR": row.get("FDR", ""),
            }
            rows.append((outdir, metadata))
        rows.sort(key=lambda item: len(item[0]), reverse=True)
        return rows

    @staticmethod
    def _metadata_for_sashimi_output(relative_path: Path, metadata_rows: list[tuple[str, dict[str, object]]]) -> dict[str, object]:
        rel = relative_path.as_posix().strip("/")
        for prefix, metadata in metadata_rows:
            if rel == prefix or rel.startswith(prefix.rstrip("/") + "/"):
                return metadata
        return {}

    @staticmethod
    def _metadata_from_sashimi_output_path(relative_path: Path) -> dict[str, object]:
        parts = list(relative_path.parts)
        if len(parts) < 2:
            return {}
        folder = str(parts[1])
        tokens = [token for token in re.split(r"[_\s]+", folder) if token]
        event_type = next((token for token in reversed(tokens) if token in EVENT_TYPES), "")
        gene = ""
        if event_type and event_type in tokens:
            event_index = len(tokens) - 1 - list(reversed(tokens)).index(event_type)
            if event_index > 0:
                gene = tokens[event_index - 1]
        return {
            "gene": gene,
            "event_type": event_type,
            "event_id": folder,
            "comparison_display_name": str(parts[0]),
        }

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
                "stable_min_junction_reads": self.current_thresholds.stable_min_junction_reads,
                "stable_min_replicates_per_group": self.current_thresholds.stable_min_replicates_per_group,
                "stable_min_pairwise_consistency": self.current_thresholds.stable_min_pairwise_consistency,
                "program_delta_dpsi": self.current_thresholds.program_delta_dpsi,
                "deg_padj": self.current_thresholds.deg_padj,
                "deg_log2fc": self.current_thresholds.deg_log2fc,
                "dexseq_qvalue": self.current_thresholds.dexseq_qvalue,
                "dtu_qvalue": self.current_thresholds.dtu_qvalue,
            },
            "tool_paths": {
                key: value
                for key, value in self.current_project.tool_paths.items()
                if key in {"rscript", "gtf", "annotation_tsv", "fasta", "genome_fasta", "transcript_fasta", "jutils", "rmats2sashimiplot", "suppa", "bam_root", "samtools"}
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

    def _resolve_gtf_path(self) -> Path | None:
        if self.current_project is None:
            return None
        configured = self.current_project.tool_paths.get("gtf")
        if configured and Path(configured).exists():
            return Path(configured)
        roots: list[Path] = []
        for root in [self.current_project.project_root, self.current_project.output_root]:
            if root:
                roots.extend([root, *root.parents])
        for value in self.current_project.input_paths.values():
            try:
                path = Path(str(value))
            except TypeError:
                continue
            if path.exists():
                roots.extend([path if path.is_dir() else path.parent, *(path if path.is_dir() else path.parent).parents])
        seen: set[Path] = set()
        for root in roots:
            if root in seen or not root.exists():
                continue
            seen.add(root)
            direct = list(root.glob("*.gtf")) + list(root.glob("*.gtf.gz"))
            recursive = list(root.glob("*gencode*.gtf")) + list(root.glob("*gencode*.gtf.gz")) + list(root.glob("*annotation*.gtf")) + list(root.glob("*annotation*.gtf.gz"))
            for candidate in direct + recursive:
                if _looks_like_genome_gtf(candidate):
                    self.current_project.tool_paths["gtf"] = str(candidate)
                    return candidate
        return None

    def _resolve_fasta_path(self, kind: str = "any") -> Path | None:
        if self.current_project is None:
            return None
        kind = (kind or "any").lower()
        configured_keys = {
            "genome": ("genome_fasta", "fasta"),
            "transcript": ("transcript_fasta", "fasta"),
            "any": ("transcript_fasta", "genome_fasta", "fasta"),
        }.get(kind, ("fasta",))
        for key in configured_keys:
            configured = self.current_project.tool_paths.get(key)
            if configured and Path(configured).exists():
                return Path(configured)

        roots: list[Path] = []
        for root in [self.current_project.project_root, self.current_project.output_root]:
            if root:
                roots.extend([root, *root.parents])
        for value in self.current_project.input_paths.values():
            try:
                path = Path(str(value))
            except TypeError:
                continue
            if path.exists():
                base = path if path.is_dir() else path.parent
                roots.extend([base, *base.parents])

        seen: set[Path] = set()
        candidates: list[Path] = []
        for root in roots:
            if root in seen or not root.exists():
                continue
            seen.add(root)
            candidates.extend(
                path
                for pattern in ("*.fa", "*.fasta", "*.fa.gz", "*.fasta.gz")
                for path in root.glob(pattern)
                if path.is_file()
            )
        if not candidates:
            return None

        def score(path: Path) -> tuple[int, int, str]:
            name = path.name.lower()
            if kind == "transcript":
                preferred = int("transcript" not in name and "cdna" not in name)
                avoided = int("genome" in name or "dna.primary" in name or "primary_assembly" in name)
            elif kind == "genome":
                preferred = int("genome" not in name and "dna.primary" not in name and "primary_assembly" not in name)
                avoided = int("transcript" in name or "cdna" in name)
            else:
                preferred = 0
                avoided = 0
            return (preferred, avoided, len(path.parts), str(path))

        chosen = sorted(candidates, key=score)[0]
        if kind == "transcript":
            self.current_project.tool_paths["transcript_fasta"] = str(chosen)
            self.current_project.tool_paths.setdefault("fasta", str(chosen))
        elif kind == "genome":
            self.current_project.tool_paths["genome_fasta"] = str(chosen)
        else:
            self.current_project.tool_paths["fasta"] = str(chosen)
        return chosen

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
        if suffix in {".tsv", ".csv", ".txt", ".md"}:
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

    def _interaction_batch_pairs(self) -> list[ComparisonPairDefinition]:
        if self.current_project is None:
            return []
        selected = self.selected_comparisons_for_display()
        if len(selected) < 2:
            return []
        comparison_lookup = {comparison.comparison_id: comparison for comparison in selected}
        pairs: list[ComparisonPairDefinition] = []
        seen: set[frozenset[str]] = set()
        for pair in self._active_comparison_pairs():
            if not pair.comparison_a or not pair.comparison_b:
                continue
            key = frozenset([pair.comparison_a, pair.comparison_b])
            if len(key) != 2 or key in seen:
                continue
            pairs.append(pair)
            seen.add(key)
        for left, right in combinations(selected, 2):
            key = frozenset([left.comparison_id, right.comparison_id])
            if key in seen:
                continue
            pair_id = f"allpair_{self._safe_output_name(left.comparison_id)}__vs__{self._safe_output_name(right.comparison_id)}"
            pairs.append(
                ComparisonPairDefinition(
                    pair_id=pair_id,
                    comparison_a=left.comparison_id,
                    comparison_b=right.comparison_id,
                    display_name=f"{comparison_lookup[left.comparison_id].display_resolved_name} compared with {comparison_lookup[right.comparison_id].display_resolved_name}",
                    experiment_group=left.display_resolved_name,
                    control_group=right.display_resolved_name,
                )
            )
            seen.add(key)
        return pairs

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

    @staticmethod
    def _dedupe_candidate_frame(frame: pd.DataFrame) -> pd.DataFrame:
        if frame is None or frame.empty:
            return pd.DataFrame()
        working = frame.copy()
        if "candidate_unique_key" not in working.columns:
            gene_key = working.get("gene_symbol", pd.Series("", index=working.index, dtype="object")).fillna("").astype(str).str.strip()
            if "gene_id" in working.columns:
                fallback = working["gene_id"].fillna("").astype(str).str.strip()
                gene_key = gene_key.where(gene_key.ne(""), fallback)
            working["candidate_unique_key"] = gene_key
        subset = ["candidate_unique_key"]
        if "comparison_id" in working.columns:
            subset = ["comparison_id", "candidate_unique_key"]
        sort_columns = [column for column in ["comparison_id", "rank", "gene_symbol", "gene_id"] if column in working.columns]
        if sort_columns:
            working = working.sort_values(sort_columns, ascending=[True] * len(sort_columns), na_position="last")
        return working.drop_duplicates(subset=subset, keep="first").reset_index(drop=True)

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
            return self._dedupe_candidate_frame(matched)
        template = frame.iloc[0].copy()
        rows = []
        for gene in missing:
            row = template.copy()
            row["comparison_id"] = self._safe_text(comparison_id) or self._safe_text(row.get("comparison_id", ""))
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
        return self._dedupe_candidate_frame(pd.concat([matched, placeholder], axis=0, ignore_index=True))

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
            support = support.sort_values(
                by=["comparison_id", "normalized_gene_id", "best_DEXSeq_qvalue", "best_DTU_qvalue", "n_support_methods"],
                ascending=[True, True, True, True, False],
                na_position="last",
            )
            support = support.drop_duplicates(subset=["comparison_id", "normalized_gene_id"], keep="first")
            join_columns = [column for column in ("comparison_id", "normalized_gene_id") if column in support.columns]
            if len(join_columns) == 2:
                existing_keys = set(
                    zip(
                        frame.get("comparison_id", pd.Series("", index=frame.index)).fillna("").astype(str),
                        frame["normalized_gene_id"].fillna("").astype(str),
                        strict=False,
                    )
                )
                support_only = support.loc[
                    ~support.apply(
                        lambda row: (
                            str(row.get("comparison_id", "")),
                            str(row.get("normalized_gene_id", "")),
                        )
                        in existing_keys,
                        axis=1,
                    )
                ].copy()
                if not support_only.empty:
                    comparison_name_lookup = {
                        item.comparison_id: item.display_resolved_name
                        for item in self.current_project.available_comparisons
                    }
                    support_rows = pd.DataFrame(
                        {
                            "comparison_id": support_only["comparison_id"].astype(str),
                            "comparison_name": support_only["comparison_id"].astype(str).map(comparison_name_lookup),
                            "gene_id": support_only.get("gene_id", support_only["normalized_gene_id"]),
                            "geneSymbol": support_only.get("gene_symbol", support_only["normalized_gene_id"]),
                            "normalized_gene_id": support_only["normalized_gene_id"],
                            "log2FC": pd.NA,
                            "deg_padj": pd.NA,
                            "deg_sig": False,
                            "representative_event_id": pd.NA,
                            "representative_event_type": pd.NA,
                            "representative_dPSI": pd.NA,
                            "representative_event_FDR": pd.NA,
                            "splicing_sig": False,
                            "combined_class": "support-only",
                            "direction_class": "support_only",
                            "agreement_class": "support_only",
                            "n_sig_events_all": 0,
                            "sig_event_types": "",
                        }
                    )
                    frame = pd.concat([frame, support_rows], axis=0, ignore_index=True, sort=False)
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

        deg_sig_source = frame["deg_sig"] if "deg_sig" in frame.columns else pd.Series(False, index=frame.index)
        splicing_sig_source = frame["splicing_sig"] if "splicing_sig" in frame.columns else pd.Series(False, index=frame.index)
        frame["DEG_significant"] = deg_sig_source.map(self._safe_bool)
        frame["rMATS_significant"] = splicing_sig_source.map(self._safe_bool)
        frame["stable_AS_confirmed"] = frame["rMATS_significant"]
        frame["DEXSeq_significant"] = frame["support_DEXSeq"] > 0
        frame["DTU_significant"] = frame["support_DTU"] > 0
        frame["n_rMATS_significant_events"] = pd.to_numeric(frame.get("n_sig_events_all"), errors="coerce").fillna(0).astype(int)
        frame["best_rMATS_FDR"] = frame["rMATS_FDR"]
        frame["dominant_rMATS_event_type"] = frame.get("representative_event_type")
        frame["best_rMATS_event_id"] = frame.get("representative_event_id")
        frame["dominant_rMATS_standardized_dPSI"] = frame["standardized_dPSI"]
        frame["dominant_rMATS_mean_junction_reads_group1"] = pd.to_numeric(
            frame.get("representative_mean_junction_reads_group1"), errors="coerce"
        )
        frame["dominant_rMATS_mean_junction_reads_group2"] = pd.to_numeric(
            frame.get("representative_mean_junction_reads_group2"), errors="coerce"
        )
        frame["dominant_rMATS_replicate_count_group1"] = pd.to_numeric(
            frame.get("representative_replicate_count_group1"), errors="coerce"
        )
        frame["dominant_rMATS_replicate_count_group2"] = pd.to_numeric(
            frame.get("representative_replicate_count_group2"), errors="coerce"
        )
        frame["dominant_rMATS_supported_replicates_group1"] = pd.to_numeric(
            frame.get("representative_supported_replicate_count_group1"), errors="coerce"
        )
        frame["dominant_rMATS_supported_replicates_group2"] = pd.to_numeric(
            frame.get("representative_supported_replicate_count_group2"), errors="coerce"
        )
        frame["dominant_rMATS_replicate_direction_consistency"] = pd.to_numeric(
            frame.get("representative_replicate_direction_consistency"), errors="coerce"
        )
        frame["dominant_rMATS_stable_support_status"] = frame.get(
            "representative_stable_support_status", pd.Series(pd.NA, index=frame.index)
        )
        frame["significant_as_event_types"] = frame.get("sig_event_types", pd.Series("", index=frame.index)).fillna("").astype(str)
        frame["as_event_summary"] = frame.apply(
            lambda row: (
                f"{int(pd.to_numeric(row.get('n_rMATS_significant_events', 0), errors='coerce') if pd.notna(pd.to_numeric(row.get('n_rMATS_significant_events', 0), errors='coerce')) else 0)} significant AS event(s)"
                + (f" [{self._safe_text(row.get('significant_as_event_types', '')).strip()}]" if self._safe_text(row.get("significant_as_event_types", "")).strip() else "")
            ),
            axis=1,
        )
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
        frame["candidate_reason"] = frame.apply(self._candidate_reason, axis=1)
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
        frame["candidate_unique_key"] = frame["gene_symbol"].fillna("").astype(str).str.strip()
        frame.loc[frame["candidate_unique_key"].eq(""), "candidate_unique_key"] = frame["normalized_gene_id"].fillna("").astype(str)
        summary_frames: list[pd.DataFrame] = []
        for comparison_id in frame["comparison_id"].dropna().astype(str).unique().tolist():
            summary = self._significant_as_event_summary_frame(comparison_id)
            if summary is not None and not summary.empty:
                summary_frames.append(summary)
        if summary_frames:
            event_summary = pd.concat(summary_frames, axis=0, ignore_index=True)
            event_summary = event_summary.drop_duplicates(subset=["comparison_id", "normalized_gene_id"], keep="first")
            frame = frame.merge(
                event_summary,
                on=["comparison_id", "normalized_gene_id"],
                how="left",
            )
        if "significant_as_event_ids" not in frame.columns:
            frame["significant_as_event_ids"] = ""
        if "significant_as_event_list" not in frame.columns:
            frame["significant_as_event_list"] = ""
        if "significant_as_stability_summary" not in frame.columns:
            frame["significant_as_stability_summary"] = ""
        for column in (
            "stable_as_min_supported_replicates_group1",
            "stable_as_min_supported_replicates_group2",
            "stable_as_min_direction_consistency",
        ):
            if column not in frame.columns:
                frame[column] = pd.NA
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
        frame = frame.drop_duplicates(subset=["comparison_id", "candidate_unique_key"], keep="first").reset_index(drop=True)
        frame["rank"] = frame.groupby("comparison_id").cumcount() + 1
        rank_values = pd.to_numeric(frame["rank"], errors="coerce")
        max_rank = rank_values.groupby(frame["comparison_id"].fillna("").astype(str), dropna=False).transform("max")
        with np.errstate(divide="ignore", invalid="ignore"):
            frame["ranking_score"] = (((max_rank - rank_values + 1) / max_rank) * 100.0).round(3)
        frame = self._attach_cross_candidate_pattern_columns(frame)
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
                "stable_AS_confirmed",
                "n_rMATS_significant_events",
                "significant_as_event_types",
                "significant_as_event_ids",
                "significant_as_event_list",
                "significant_as_stability_summary",
                "stable_as_min_supported_replicates_group1",
                "stable_as_min_supported_replicates_group2",
                "stable_as_min_direction_consistency",
                "as_event_summary",
                "cross_candidate_pattern",
                "cross_pattern_reason",
                "cross_candidate_comparisons",
                "cross_n_comparisons_candidate",
                "cross_best_overall_rank",
                "cross_comparison_effect_summary",
                "shared_candidate",
                "specific_candidate",
                "gained_candidate",
                "lost_candidate",
                "best_rMATS_FDR",
                "max_abs_dPSI",
                "dominant_rMATS_event_type",
                "best_rMATS_event_id",
                "dominant_rMATS_standardized_dPSI",
                "dominant_rMATS_mean_junction_reads_group1",
                "dominant_rMATS_mean_junction_reads_group2",
                "dominant_rMATS_replicate_count_group1",
                "dominant_rMATS_replicate_count_group2",
                "dominant_rMATS_supported_replicates_group1",
                "dominant_rMATS_supported_replicates_group2",
                "dominant_rMATS_replicate_direction_consistency",
                "dominant_rMATS_stable_support_status",
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
                "ranking_score",
            ]
        ].copy()

    def _attach_cross_candidate_pattern_columns(self, frame: pd.DataFrame) -> pd.DataFrame:
        pattern_columns = [
            "cross_candidate_pattern",
            "cross_pattern_reason",
            "cross_candidate_comparisons",
            "cross_n_comparisons_candidate",
            "cross_best_overall_rank",
            "cross_comparison_effect_summary",
            "shared_candidate",
            "specific_candidate",
            "gained_candidate",
            "lost_candidate",
        ]
        if frame is None or frame.empty:
            return frame
        working = frame.copy()
        for column in pattern_columns:
            if column not in working.columns:
                working[column] = ""
        matrix = self._build_cross_comparison_candidate_matrix(working)
        if matrix is None or matrix.empty or "gene_id" not in matrix.columns:
            return working
        matrix_columns = [
            "gene_id",
            "candidate_pattern",
            "pattern_reason",
            "candidate_comparisons",
            "n_comparisons_candidate",
            "best_overall_rank",
        ]
        matrix_subset = matrix[[column for column in matrix_columns if column in matrix.columns]].drop_duplicates(
            subset=["gene_id"],
            keep="first",
        )
        rename_map = {
            "candidate_pattern": "cross_candidate_pattern",
            "pattern_reason": "cross_pattern_reason",
            "candidate_comparisons": "cross_candidate_comparisons",
            "n_comparisons_candidate": "cross_n_comparisons_candidate",
            "best_overall_rank": "cross_best_overall_rank",
        }
        matrix_subset = matrix_subset.rename(columns=rename_map)
        effect_summary = self._cross_candidate_effect_summary(working)
        if not effect_summary.empty:
            matrix_subset = matrix_subset.merge(effect_summary, on="gene_id", how="left")
        working = working.drop(columns=[column for column in rename_map.values() if column in working.columns], errors="ignore")
        if "cross_comparison_effect_summary" in working.columns:
            working = working.drop(columns=["cross_comparison_effect_summary"], errors="ignore")
        working = working.merge(matrix_subset, on="gene_id", how="left")
        for column in pattern_columns:
            if column not in working.columns:
                working[column] = ""
            working[column] = working[column].fillna("")
        pattern = working["cross_candidate_pattern"].fillna("").astype(str)
        working["shared_candidate"] = pattern.isin({"shared_all", "shared_subset"})
        working["specific_candidate"] = pattern.isin(
            {
                "comparison_specific",
                "DEG_only_in_one_condition",
                "splicing_only_in_one_condition",
            }
        )
        working["gained_candidate"] = pattern.eq("gained")
        working["lost_candidate"] = pattern.eq("lost")
        return working

    def _cross_candidate_effect_summary(self, frame: pd.DataFrame) -> pd.DataFrame:
        if frame is None or frame.empty or "gene_id" not in frame.columns:
            return pd.DataFrame(columns=["gene_id", "cross_comparison_effect_summary"])
        display_lookup = {
            item.comparison_id: item.display_resolved_name
            for item in self.current_project.available_comparisons
        } if self.current_project is not None else {}
        order = self._configured_comparison_order()
        observed = frame.get("comparison_id", pd.Series("", index=frame.index)).dropna().astype(str).unique().tolist()
        ordered_ids = [item for item in order if item in observed]
        ordered_ids.extend(item for item in observed if item not in ordered_ids)
        order_index = {comparison_id: index for index, comparison_id in enumerate(ordered_ids)}

        def _fmt_number(value) -> str:
            numeric = pd.to_numeric(value, errors="coerce")
            if pd.isna(numeric):
                return "NA"
            return f"{float(numeric):.3g}"

        rows: list[dict[str, object]] = []
        for gene_id, subset in frame.groupby("gene_id", dropna=False):
            local = subset.copy()
            local["_comparison_order"] = local.get("comparison_id", pd.Series("", index=local.index)).astype(str).map(
                lambda value: order_index.get(value, len(order_index))
            )
            local = local.sort_values(["_comparison_order", "rank"], ascending=[True, True], na_position="last")
            parts: list[str] = []
            for _, row in local.iterrows():
                comparison_id = self._safe_text(row.get("comparison_id"))
                comparison_name = display_lookup.get(comparison_id, comparison_id)
                event_type = self._safe_text(row.get("dominant_rMATS_event_type")).strip()
                event_id = self._safe_text(row.get("best_rMATS_event_id")).strip()
                event_label = f"{event_type}:{event_id}" if event_type or event_id else "event=NA"
                parts.append(
                    "; ".join(
                        [
                            f"{comparison_name}",
                            f"rank={_fmt_number(row.get('rank'))}",
                            f"tier={row.get('candidate_tier', '')}",
                            f"evidence={row.get('evidence_class', '')}",
                            f"dPSI={_fmt_number(row.get('dominant_rMATS_standardized_dPSI'))}",
                            f"FDR={_fmt_number(row.get('best_rMATS_FDR'))}",
                            event_label,
                        ]
                    )
                )
            rows.append(
                {
                    "gene_id": gene_id,
                    "cross_comparison_effect_summary": " | ".join(parts),
                }
            )
        return pd.DataFrame(rows)

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
        selected_ids = [item.comparison_id for item in self._selected_comparisons()]
        total_selected = len(selected_ids) or len(working["comparison_id"].dropna().astype(str).unique().tolist())
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

    @staticmethod
    def _empty_cross_as_interaction_tables() -> dict[str, pd.DataFrame]:
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
            "interaction_log_odds",
            "interaction_se",
            "interaction_z",
            "interaction_pvalue",
            "interaction_FDR",
            "interaction_model",
            "interaction_model_status",
            "interaction_overdispersion_scale",
            "interaction_model_n_observations",
            "interaction_model_df_resid",
            "interaction_model_pearson_chi2",
            "interaction_model_deviance",
            "interaction_model_total_reads",
            "quasi_binomial_dispersion",
            "quasi_binomial_overdispersed",
            "beta_binomial_rho_estimate",
            "beta_binomial_interaction_log_odds",
            "beta_binomial_likelihood_ratio",
            "beta_binomial_pvalue",
            "beta_binomial_FDR",
            "beta_binomial_rho",
            "beta_binomial_full_loglik",
            "beta_binomial_null_loglik",
            "beta_binomial_model_status",
            "beta_binomial_warning",
            "model_valid",
            "model_diagnostic_class",
            "model_warning",
            "eb_prior_variance",
            "shrinkage_weight",
            "shrunken_interaction_log_odds",
            "shrunken_interaction_z",
            "shrunken_interaction_pvalue",
            "shrunken_interaction_FDR",
            "interaction_priority_score",
            "simple_delta_score",
            "simple_delta_priority_score",
            "delta_dPSI_rank",
            "interaction_model_rank",
            "delta_vs_model_rank_shift",
            "simple_delta_top20",
            "model_top20",
            "benchmark_class",
            "benchmark_score_delta",
            "interaction_class",
            "read_support_class",
            "replicate_consistency_score",
            "mean_total_reads_A",
            "mean_total_reads_B",
            "replicate_count_A",
            "replicate_count_B",
            "exclusive_side",
            "counterpart_event_found",
            "significant_side_mean_total_reads",
            "non_significant_side_mean_total_reads",
            "significant_side_replicate_count",
            "non_significant_side_replicate_count",
            "non_significant_side_dPSI",
            "non_significant_side_FDR",
            "non_significant_side_read_support",
            "event_class",
            "unified_result_type",
            "unified_display_class",
            "support_class",
            "abs_dPSI_significance_cutoff_used",
            "large_delta_dPSI_cutoff_used",
            "min_total_reads_used",
        ]
        empty = pd.DataFrame(columns=base_columns)
        return {
            "unified_event_table": empty.copy(),
            "shared_interaction": empty.copy(),
            "opposite_direction": empty.copy(),
            "large_delta": empty.copy(),
            "exclusive": empty.copy(),
            "benchmark": empty.copy(),
            "diagnostics": pd.DataFrame(
                columns=[
                    "pair_id",
                    "pair_name",
                    "comparison_A",
                    "comparison_B",
                    "n_events_A",
                    "n_events_B",
                    "n_shared_significant",
                    "n_opposite_direction",
                    "n_large_delta",
                    "n_exclusive_supported",
                    "n_exclusive_low_power",
                    "n_exclusive_counterpart_not_found",
                    "n_exclusive_not_shared_significant",
                    "n_model_converged",
                    "n_model_fallback",
                    "n_model_failed",
                    "median_quasi_binomial_dispersion",
                    "max_quasi_binomial_dispersion",
                    "median_beta_binomial_rho_estimate",
                    "median_beta_binomial_rho",
                    "n_beta_binomial_converged",
                    "n_beta_binomial_failed",
                    "n_beta_binomial_FDR_0_05",
                    "n_overdispersed_events",
                    "delta_top20_model_top20_overlap",
                    "delta_top20_model_top20_overlap_fraction",
                    "delta_vs_model_spearman_rank_correlation",
                    "n_delta_only_top20",
                    "n_model_only_top20",
                    "abs_dPSI_significance_cutoff_used",
                    "large_delta_dPSI_cutoff_used",
                    "min_total_reads_used",
                ]
            ),
        }

    @staticmethod
    def _empty_cross_as_interaction_batch_tables() -> dict[str, pd.DataFrame]:
        return {
            "batch_summary": pd.DataFrame(
                columns=[
                    "batch_pair_order",
                    "pair_id",
                    "pair_name",
                    "comparison_A",
                    "comparison_B",
                    "n_shared_significant",
                    "n_high_confidence_interaction",
                    "n_moderate_interaction",
                    "n_opposite_direction",
                    "n_same_direction_large_delta",
                    "n_exclusive",
                    "n_exclusive_A",
                    "n_exclusive_B",
                    "n_exclusive_supported",
                    "n_exclusive_A_supported",
                    "n_exclusive_B_supported",
                    "n_exclusive_low_power",
                    "n_exclusive_counterpart_not_found",
                    "top_gene",
                    "top_event_type",
                    "top_abs_delta_dPSI",
                    "top_interaction_priority_score",
                    "delta_top20_model_top20_overlap",
                    "delta_vs_model_spearman_rank_correlation",
                    "abs_dPSI_significance_cutoff_used",
                    "large_delta_dPSI_cutoff_used",
                    "min_total_reads_used",
                ]
            ),
            "batch_top_interactions": ProjectService._empty_cross_as_interaction_tables()["shared_interaction"].copy(),
            "batch_exclusive": ProjectService._empty_cross_as_interaction_tables()["exclusive"].copy(),
            "batch_benchmark": ProjectService._empty_cross_as_interaction_tables()["benchmark"].copy(),
            "batch_diagnostics": ProjectService._empty_cross_as_interaction_tables()["diagnostics"].copy(),
        }

    @staticmethod
    def _parse_float_list(value: object) -> list[float]:
        if value is None:
            return []
        try:
            if pd.isna(value):
                return []
        except (TypeError, ValueError):
            pass
        values: list[float] = []
        for token in str(value).replace(";", ",").split(","):
            text = token.strip()
            if not text or text.upper() in {"NA", "NAN", "NONE"}:
                continue
            try:
                values.append(float(text))
            except ValueError:
                continue
        return values

    @staticmethod
    def _paired_read_totals(ijc_values: list[float], sjc_values: list[float]) -> list[float]:
        n_items = max(len(ijc_values), len(sjc_values))
        totals: list[float] = []
        for idx in range(n_items):
            ijc = ijc_values[idx] if idx < len(ijc_values) and np.isfinite(ijc_values[idx]) else 0.0
            sjc = sjc_values[idx] if idx < len(sjc_values) and np.isfinite(sjc_values[idx]) else 0.0
            totals.append(float(ijc) + float(sjc))
        return totals

    def _apply_stable_as_columns(self, frame: pd.DataFrame, *, dpsi_column: str = "dPSI", reverse_splicing: bool = False) -> pd.DataFrame:
        working = frame.copy()
        fdr = pd.to_numeric(working.get("FDR"), errors="coerce")
        standardized_dpsi = pd.to_numeric(working.get(dpsi_column), errors="coerce")
        return apply_stable_as_filter(
            working,
            fdr=fdr,
            standardized_dpsi=standardized_dpsi,
            reverse_splicing=reverse_splicing,
            thresholds=self.current_thresholds,
        )

    @staticmethod
    def _safe_logit(value: float) -> float:
        if not np.isfinite(value):
            return float("nan")
        clipped = min(max(float(value), 1e-4), 1.0 - 1e-4)
        return float(np.log(clipped / (1.0 - clipped)))

    @staticmethod
    def _normal_two_sided_pvalue(z_value: float) -> float:
        if not np.isfinite(z_value):
            return float("nan")
        return float(math.erfc(abs(float(z_value)) / math.sqrt(2.0)))

    @staticmethod
    def _bh_adjust(pvalues: pd.Series) -> pd.Series:
        numeric = pd.to_numeric(pvalues, errors="coerce")
        adjusted = pd.Series(np.nan, index=numeric.index, dtype="float64")
        valid = numeric.dropna()
        if valid.empty:
            return adjusted
        ordered = valid.sort_values()
        m = float(len(ordered))
        previous = 1.0
        for rank in range(len(ordered), 0, -1):
            item_index = ordered.index[rank - 1]
            pvalue = float(ordered.iloc[rank - 1])
            estimate = min(previous, pvalue * m / float(rank))
            adjusted.loc[item_index] = min(max(estimate, 0.0), 1.0)
            previous = estimate
        return adjusted

    def _interaction_group_stats(self, psi_values: list[float], ijc_values: list[float], sjc_values: list[float]) -> dict[str, float]:
        count_psi: list[float] = []
        totals: list[float] = []
        max_len = max(len(ijc_values), len(sjc_values))
        for idx in range(max_len):
            ijc = ijc_values[idx] if idx < len(ijc_values) and np.isfinite(ijc_values[idx]) else 0.0
            sjc = sjc_values[idx] if idx < len(sjc_values) and np.isfinite(sjc_values[idx]) else 0.0
            total = ijc + sjc
            totals.append(total)
            count_psi.append((ijc + 0.5) / (total + 1.0) if total > 0 else float("nan"))
        usable_psi = [float(value) for value in psi_values if np.isfinite(value)]
        if not usable_psi:
            usable_psi = [float(value) for value in count_psi if np.isfinite(value)]
        logits = [self._safe_logit(value) for value in usable_psi if np.isfinite(value)]
        n = len(logits)
        mean_logit = float(np.nanmean(logits)) if logits else float("nan")
        var_logit = float(np.nanvar(logits, ddof=1)) if len(logits) > 1 else float("nan")
        mean_psi = float(np.nanmean(usable_psi)) if usable_psi else float("nan")
        finite_totals = [float(value) for value in totals if np.isfinite(value)]
        mean_total = float(np.nanmean(finite_totals)) if finite_totals else float("nan")
        total_sum = float(np.nansum(finite_totals)) if finite_totals else 0.0
        return {
            "mean_psi": mean_psi,
            "mean_logit": mean_logit,
            "var_logit": var_logit,
            "n": float(n),
            "mean_total_reads": mean_total,
            "total_reads_sum": total_sum,
        }

    def _replicate_direction_consistency(self, psi_group1: list[float], psi_group2: list[float], dpsi_value: float) -> float:
        return pairwise_direction_consistency(psi_group1, psi_group2, dpsi_value)

    def _fit_binomial_glm_interaction(self, row: pd.Series) -> dict[str, object]:
        """Fit IJC/SJC ~ comparison + group + comparison:group.

        interaction_log_odds is reported as comparison A group effect minus
        comparison B group effect so its sign matches delta_dPSI = dPSI_A -
        dPSI_B. The GLM coefficient itself is B-minus-A, so we flip it.
        """

        def as_list(key: str) -> list[float]:
            value = row.get(key, [])
            if isinstance(value, list):
                return [float(item) for item in value if np.isfinite(item)]
            return self._parse_float_list(value)

        records: list[dict[str, float]] = []
        specs = [
            ("A", 0.0, 1.0, as_list("ijc_group1_values_A"), as_list("sjc_group1_values_A")),
            ("A", 0.0, 0.0, as_list("ijc_group2_values_A"), as_list("sjc_group2_values_A")),
            ("B", 1.0, 1.0, as_list("ijc_group1_values_B"), as_list("sjc_group1_values_B")),
            ("B", 1.0, 0.0, as_list("ijc_group2_values_B"), as_list("sjc_group2_values_B")),
        ]
        for _label, comparison_b, group1, ijc_values, sjc_values in specs:
            for ijc, sjc in zip(ijc_values, sjc_values, strict=False):
                if not np.isfinite(ijc) or not np.isfinite(sjc):
                    continue
                total = float(ijc) + float(sjc)
                if total <= 0:
                    continue
                records.append(
                    {
                        "ijc": max(float(ijc), 0.0),
                        "sjc": max(float(sjc), 0.0),
                        "comparison_b": comparison_b,
                        "group1": group1,
                        "interaction": comparison_b * group1,
                    }
                )
        if len(records) < 4:
            return {
                "interaction_log_odds": float("nan"),
                "interaction_se": float("nan"),
                "interaction_z": float("nan"),
                "interaction_pvalue": float("nan"),
                "interaction_model": "binomial_glm",
                "interaction_model_status": "insufficient_replicate_count",
                "interaction_overdispersion_scale": float("nan"),
                "interaction_model_n_observations": len(records),
                "interaction_model_df_resid": float("nan"),
                "interaction_model_pearson_chi2": float("nan"),
                "interaction_model_deviance": float("nan"),
                "interaction_model_total_reads": float("nan"),
                "quasi_binomial_dispersion": float("nan"),
                "quasi_binomial_overdispersed": False,
                "beta_binomial_rho_estimate": float("nan"),
                "model_valid": False,
                "model_diagnostic_class": "insufficient_replicates",
                "model_warning": "Fewer than four replicate group observations with IJC/SJC support.",
            }
        try:
            import statsmodels.api as sm

            frame = pd.DataFrame(records)
            endog = frame[["ijc", "sjc"]].to_numpy(dtype=float)
            exog = frame[["comparison_b", "group1", "interaction"]].to_numpy(dtype=float)
            exog = sm.add_constant(exog, has_constant="add")
            model = sm.GLM(endog, exog, family=sm.families.Binomial())
            result = model.fit(maxiter=100, disp=0)
            raw_coef = float(result.params[3])
            raw_se = float(result.bse[3])
            scale = 1.0
            pearson_chi2 = float(getattr(result, "pearson_chi2", float("nan")))
            df_resid = float(getattr(result, "df_resid", float("nan")))
            deviance = float(getattr(result, "deviance", float("nan")))
            try:
                if np.isfinite(df_resid) and df_resid > 0 and np.isfinite(pearson_chi2):
                    scale = max(1.0, pearson_chi2 / df_resid)
            except Exception:
                scale = 1.0
            adjusted_se = raw_se * math.sqrt(scale) if np.isfinite(raw_se) and raw_se > 0 else float("nan")
            interaction_log_odds = -raw_coef
            interaction_z = (-raw_coef / adjusted_se) if np.isfinite(adjusted_se) and adjusted_se > 0 else float("nan")
            interaction_pvalue = self._normal_two_sided_pvalue(interaction_z)
            total_reads = float((frame["ijc"] + frame["sjc"]).sum())
            mean_total = float((frame["ijc"] + frame["sjc"]).mean()) if not frame.empty else float("nan")
            rho = (scale - 1.0) / max(mean_total - 1.0, 1.0) if np.isfinite(mean_total) and mean_total > 1 else 0.0
            rho = float(np.clip(rho, 0.0, 0.99)) if np.isfinite(rho) else float("nan")
            overdispersed = bool(np.isfinite(scale) and scale >= 1.5)
            model_valid = bool(np.isfinite(interaction_pvalue) and np.isfinite(adjusted_se) and adjusted_se > 0)
            return {
                "interaction_log_odds": interaction_log_odds,
                "interaction_se": adjusted_se,
                "interaction_z": interaction_z,
                "interaction_pvalue": interaction_pvalue,
                "interaction_model": "binomial_glm_ijc_sjc_with_pearson_overdispersion",
                "interaction_model_status": "converged" if bool(getattr(result, "converged", True)) else "not_converged",
                "interaction_overdispersion_scale": scale,
                "interaction_model_n_observations": len(records),
                "interaction_model_df_resid": df_resid,
                "interaction_model_pearson_chi2": pearson_chi2,
                "interaction_model_deviance": deviance,
                "interaction_model_total_reads": total_reads,
                "quasi_binomial_dispersion": scale,
                "quasi_binomial_overdispersed": overdispersed,
                "beta_binomial_rho_estimate": rho,
                "model_valid": model_valid,
                "model_diagnostic_class": "quasi_binomial_overdispersed" if overdispersed else "quasi_binomial_ok",
                "model_warning": "" if model_valid else "Model converged but interaction coefficient or SE was not finite.",
            }
        except Exception as exc:
            return {
                "interaction_log_odds": float("nan"),
                "interaction_se": float("nan"),
                "interaction_z": float("nan"),
                "interaction_pvalue": float("nan"),
                "interaction_model": "binomial_glm_ijc_sjc_with_pearson_overdispersion",
                "interaction_model_status": f"failed: {exc}",
                "interaction_overdispersion_scale": float("nan"),
                "interaction_model_n_observations": len(records),
                "interaction_model_df_resid": float("nan"),
                "interaction_model_pearson_chi2": float("nan"),
                "interaction_model_deviance": float("nan"),
                "interaction_model_total_reads": float("nan"),
                "quasi_binomial_dispersion": float("nan"),
                "quasi_binomial_overdispersed": False,
                "beta_binomial_rho_estimate": float("nan"),
                "model_valid": False,
                "model_diagnostic_class": "model_failed",
                "model_warning": str(exc),
            }

    def _fit_beta_binomial_interaction(self, row: pd.Series) -> dict[str, object]:
        """Fit a beta-binomial interaction model on replicate IJC/SJC counts.

        The full model is inclusion ~ comparison + group + comparison:group.
        The null model removes the interaction term. The likelihood-ratio
        p-value tests whether the AS effect differs between comparison A and B.
        """

        def as_list(key: str) -> list[float]:
            value = row.get(key, [])
            if isinstance(value, list):
                return [float(item) for item in value if np.isfinite(item)]
            return self._parse_float_list(value)

        records: list[tuple[float, float, float, float, float]] = []
        specs = [
            (0.0, 1.0, as_list("ijc_group1_values_A"), as_list("sjc_group1_values_A")),
            (0.0, 0.0, as_list("ijc_group2_values_A"), as_list("sjc_group2_values_A")),
            (1.0, 1.0, as_list("ijc_group1_values_B"), as_list("sjc_group1_values_B")),
            (1.0, 0.0, as_list("ijc_group2_values_B"), as_list("sjc_group2_values_B")),
        ]
        for comparison_b, group1, ijc_values, sjc_values in specs:
            for ijc, sjc in zip(ijc_values, sjc_values, strict=False):
                if not np.isfinite(ijc) or not np.isfinite(sjc):
                    continue
                success = max(float(ijc), 0.0)
                failure = max(float(sjc), 0.0)
                total = success + failure
                if total <= 0:
                    continue
                records.append((success, total, comparison_b, group1, comparison_b * group1))
        if len(records) < 4:
            return {
                "beta_binomial_interaction_log_odds": float("nan"),
                "beta_binomial_likelihood_ratio": float("nan"),
                "beta_binomial_pvalue": float("nan"),
                "beta_binomial_rho": float("nan"),
                "beta_binomial_full_loglik": float("nan"),
                "beta_binomial_null_loglik": float("nan"),
                "beta_binomial_model_status": "insufficient_replicate_count",
                "beta_binomial_warning": "Fewer than four replicate group observations with IJC/SJC support.",
            }
        try:
            from scipy.optimize import minimize
            from scipy.special import betaln, expit, gammaln
            from scipy.stats import chi2

            frame = pd.DataFrame(records, columns=["success", "total", "comparison_b", "group1", "interaction"])
            y = frame["success"].to_numpy(dtype=float)
            n = frame["total"].to_numpy(dtype=float)
            x_full = np.column_stack(
                [
                    np.ones(len(frame), dtype=float),
                    frame["comparison_b"].to_numpy(dtype=float),
                    frame["group1"].to_numpy(dtype=float),
                    frame["interaction"].to_numpy(dtype=float),
                ]
            )
            x_null = x_full[:, :3]

            def beta_binomial_loglik(params: np.ndarray, design: np.ndarray) -> float:
                beta = params[:-1]
                rho = float(expit(params[-1]) * 0.95)
                rho = min(max(rho, 1e-6), 0.95)
                mu = np.clip(expit(design @ beta), 1e-6, 1.0 - 1e-6)
                concentration = (1.0 - rho) / rho
                alpha = np.clip(mu * concentration, 1e-6, 1e8)
                beta_param = np.clip((1.0 - mu) * concentration, 1e-6, 1e8)
                log_comb = gammaln(n + 1.0) - gammaln(y + 1.0) - gammaln(n - y + 1.0)
                value = log_comb + betaln(y + alpha, n - y + beta_param) - betaln(alpha, beta_param)
                return float(np.sum(value))

            def negative_loglik(params: np.ndarray, design: np.ndarray) -> float:
                value = beta_binomial_loglik(params, design)
                if not np.isfinite(value):
                    return 1e100
                return -value

            observed_rate = np.clip((y.sum() + 0.5) / (n.sum() + 1.0), 1e-4, 1.0 - 1e-4)
            intercept = float(np.log(observed_rate / (1.0 - observed_rate)))
            rho_start = float(pd.to_numeric(row.get("beta_binomial_rho_estimate", np.nan), errors="coerce"))
            if not np.isfinite(rho_start) or rho_start <= 0:
                rho_start = 0.02
            rho_start = min(max(rho_start, 1e-5), 0.80)
            rho_param_start = float(np.log((rho_start / 0.95) / max(1.0 - (rho_start / 0.95), 1e-6)))
            full_start = np.array([intercept, 0.0, 0.0, 0.0, rho_param_start], dtype=float)
            null_start = np.array([intercept, 0.0, 0.0, rho_param_start], dtype=float)
            bounds_full = [(-20.0, 20.0), (-20.0, 20.0), (-20.0, 20.0), (-20.0, 20.0), (-12.0, 3.0)]
            bounds_null = [(-20.0, 20.0), (-20.0, 20.0), (-20.0, 20.0), (-12.0, 3.0)]
            full_fit = minimize(
                negative_loglik,
                full_start,
                args=(x_full,),
                method="L-BFGS-B",
                bounds=bounds_full,
                options={"maxiter": 120, "ftol": 1e-7},
            )
            null_fit = minimize(
                negative_loglik,
                null_start,
                args=(x_null,),
                method="L-BFGS-B",
                bounds=bounds_null,
                options={"maxiter": 120, "ftol": 1e-7},
            )
            full_ll = beta_binomial_loglik(full_fit.x, x_full)
            null_ll = beta_binomial_loglik(null_fit.x, x_null)
            lr_stat = max(0.0, 2.0 * (full_ll - null_ll))
            pvalue = float(chi2.sf(lr_stat, df=1))
            rho = float(expit(full_fit.x[-1]) * 0.95)
            warning_parts = []
            if not full_fit.success:
                warning_parts.append(f"full model: {full_fit.message}")
            if not null_fit.success:
                warning_parts.append(f"null model: {null_fit.message}")
            return {
                "beta_binomial_interaction_log_odds": -float(full_fit.x[3]),
                "beta_binomial_likelihood_ratio": lr_stat,
                "beta_binomial_pvalue": pvalue,
                "beta_binomial_rho": rho,
                "beta_binomial_full_loglik": full_ll,
                "beta_binomial_null_loglik": null_ll,
                "beta_binomial_model_status": "converged" if full_fit.success and null_fit.success else "fit_warning",
                "beta_binomial_warning": "; ".join(warning_parts),
            }
        except Exception as exc:
            return {
                "beta_binomial_interaction_log_odds": float("nan"),
                "beta_binomial_likelihood_ratio": float("nan"),
                "beta_binomial_pvalue": float("nan"),
                "beta_binomial_rho": float("nan"),
                "beta_binomial_full_loglik": float("nan"),
                "beta_binomial_null_loglik": float("nan"),
                "beta_binomial_model_status": "failed",
                "beta_binomial_warning": str(exc),
            }

    def _cross_as_interaction_source_frame(
        self,
        comparison_id: str | None,
        abs_dpsi_cutoff: float,
        min_total_reads: int,
        *,
        significant_only: bool = True,
        candidate_keys: set[str] | None = None,
    ) -> pd.DataFrame:
        if not comparison_id or self.current_project is None:
            return pd.DataFrame()
        comparison = next(
            (item for item in self.current_project.available_comparisons if item.comparison_id == comparison_id),
            None,
        )
        if comparison is None:
            return pd.DataFrame()
        mode = self.current_project.selected_rmats_mode or "JC"
        annotation_lookup = self._annotation_gene_symbol_lookup()
        left_label, right_label = self._comparison_splicing_direction_labels(comparison)
        rows: list[dict[str, object]] = []
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
            fdr_values = pd.to_numeric(working.get("FDR"), errors="coerce")
            dpsi_values = pd.to_numeric(working.get("IncLevelDifference"), errors="coerce")
            if comparison.reverse_splicing:
                dpsi_values = -dpsi_values
            working["FDR"] = fdr_values
            working["dPSI"] = dpsi_values
            working = self._apply_stable_as_columns(working, dpsi_column="dPSI", reverse_splicing=bool(comparison.reverse_splicing))
            fdr_values = pd.to_numeric(working.get("FDR"), errors="coerce")
            dpsi_values = pd.to_numeric(working.get("dPSI"), errors="coerce")
            normalized_gene_ids = self._normalized_gene_id_series(working)
            if "geneSymbol" in working.columns:
                gene_symbols = working["geneSymbol"].astype("object")
            elif "gene_symbol" in working.columns:
                gene_symbols = working["gene_symbol"].astype("object")
            else:
                gene_symbols = pd.Series(pd.NA, index=working.index, dtype="object")
            missing_symbols = gene_symbols.isna() | gene_symbols.astype(str).str.strip().eq("")
            if missing_symbols.any() and annotation_lookup:
                gene_symbols.loc[missing_symbols] = normalized_gene_ids.loc[missing_symbols].map(annotation_lookup)
            gene_symbols = gene_symbols.combine_first(normalized_gene_ids)

            event_ids = working.apply(lambda row: event_uid_from_row(row, event_type), axis=1).astype(str)
            coordinates = working.apply(lambda row: self._rmats_event_coordinates(row, event_type), axis=1).astype(str)
            matching_ids = event_ids.where(event_ids.str.strip().ne(""), coordinates)
            key_genes = normalized_gene_ids.fillna("").astype(str)
            key_genes = key_genes.where(key_genes.str.strip().ne(""), gene_symbols.fillna("").astype(str))
            event_keys = (
                key_genes
                + "|"
                + event_type
                + "|"
                + matching_ids.fillna("").astype(str)
            )
            significant_mask = working["stable_as_event"].fillna(False).map(self._safe_bool) & dpsi_values.abs().ge(abs_dpsi_cutoff)
            if candidate_keys is not None:
                keep_mask = event_keys.isin(candidate_keys)
            elif significant_only:
                keep_mask = significant_mask
            else:
                keep_mask = pd.Series(True, index=working.index)
            working = working.loc[keep_mask].copy()
            if working.empty:
                continue
            fdr_values = fdr_values.loc[working.index]
            dpsi_values = dpsi_values.loc[working.index]
            gene_symbols = gene_symbols.loc[working.index]
            normalized_gene_ids = normalized_gene_ids.loc[working.index]
            event_ids = event_ids.loc[working.index]
            coordinates = coordinates.loc[working.index]
            event_keys = event_keys.loc[working.index]
            significant_mask = significant_mask.loc[working.index]
            stable_as_event = working["stable_as_event"].loc[working.index]
            rmats_sig_before_stability = working["rmats_sig_before_stability"].loc[working.index]

            for idx, row in working.iterrows():
                gene_id = str(normalized_gene_ids.loc[idx] or "").strip()
                gene_symbol = str(gene_symbols.loc[idx] or gene_id or "").strip()
                event_id = event_ids.loc[idx]
                coordinate_text = coordinates.loc[idx]
                event_matching_key = event_keys.loc[idx]
                fdr = fdr_values.loc[idx]
                dpsi = dpsi_values.loc[idx]

                psi_1 = self._parse_float_list(row.get("IncLevel1"))
                psi_2 = self._parse_float_list(row.get("IncLevel2"))
                ijc_1 = self._parse_float_list(row.get("IJC_SAMPLE_1"))
                ijc_2 = self._parse_float_list(row.get("IJC_SAMPLE_2"))
                sjc_1 = self._parse_float_list(row.get("SJC_SAMPLE_1"))
                sjc_2 = self._parse_float_list(row.get("SJC_SAMPLE_2"))
                if comparison.reverse_splicing:
                    psi_1, psi_2 = psi_2, psi_1
                    ijc_1, ijc_2 = ijc_2, ijc_1
                    sjc_1, sjc_2 = sjc_2, sjc_1

                stats_1 = self._interaction_group_stats(psi_1, ijc_1, sjc_1)
                stats_2 = self._interaction_group_stats(psi_2, ijc_2, sjc_2)
                effect_logit = stats_1["mean_logit"] - stats_2["mean_logit"]
                var_terms = []
                if stats_1["n"] > 1 and np.isfinite(stats_1["var_logit"]):
                    var_terms.append(stats_1["var_logit"] / max(stats_1["n"], 1.0))
                if stats_2["n"] > 1 and np.isfinite(stats_2["var_logit"]):
                    var_terms.append(stats_2["var_logit"] / max(stats_2["n"], 1.0))
                effect_se = float(np.sqrt(np.nansum(var_terms))) if var_terms else float("nan")
                if not np.isfinite(effect_se) or effect_se <= 0:
                    effect_se = float(
                        np.sqrt(
                            1.0 / max(stats_1["total_reads_sum"], 1.0)
                            + 1.0 / max(stats_2["total_reads_sum"], 1.0)
                        )
                    )
                mean_total_reads = float(np.nanmean([stats_1["mean_total_reads"], stats_2["mean_total_reads"]]))
                replicate_count = int(max(stats_1["n"], stats_2["n"]))
                consistency = self._replicate_direction_consistency(psi_1, psi_2, float(dpsi) if pd.notna(dpsi) else float("nan"))
                significant = bool(significant_mask.loc[idx])
                support_ok = bool(np.isfinite(mean_total_reads) and mean_total_reads >= min_total_reads and replicate_count >= 2)
                rows.append(
                    {
                        "comparison_id": comparison.comparison_id,
                        "comparison_display_name": comparison.display_resolved_name,
                        "gene": gene_symbol or gene_id,
                        "gene_symbol": gene_symbol,
                        "gene_id": gene_id,
                        "event_type": event_type,
                        "event_id": event_id,
                        "event_matching_key": event_matching_key,
                        "coordinates": coordinate_text,
                        "dPSI": float(dpsi) if pd.notna(dpsi) else float("nan"),
                        "FDR": float(fdr) if pd.notna(fdr) else float("nan"),
                        "significant": significant,
                        "rmats_sig_before_stability": bool(rmats_sig_before_stability.loc[idx]),
                        "stable_as_event": bool(stable_as_event.loc[idx]),
                        "direction": (
                            f"inclusion higher in {left_label} than {right_label}"
                            if pd.notna(dpsi) and float(dpsi) > 0
                            else f"inclusion lower in {left_label} than {right_label}"
                            if pd.notna(dpsi) and float(dpsi) < 0
                            else "no directional shift"
                        ),
                        "group1_label": left_label,
                        "group2_label": right_label,
                        "psi_group1_mean": stats_1["mean_psi"],
                        "psi_group2_mean": stats_2["mean_psi"],
                        "psi_group1_values": psi_1,
                        "psi_group2_values": psi_2,
                        "ijc_group1_values": ijc_1,
                        "ijc_group2_values": ijc_2,
                        "sjc_group1_values": sjc_1,
                        "sjc_group2_values": sjc_2,
                        "mean_total_reads": mean_total_reads,
                        "replicate_count": replicate_count,
                        "effect_logit": effect_logit,
                        "effect_se": effect_se,
                        "replicate_consistency_score": consistency,
                        "read_support_ok": support_ok,
                        "abs_dPSI_significance_cutoff_used": abs_dpsi_cutoff,
                        "min_total_reads_used": min_total_reads,
                    }
                )
        if not rows:
            return pd.DataFrame()
        output = pd.DataFrame(rows)
        output = output.sort_values(["significant", "FDR", "event_type", "gene"], ascending=[False, True, True, True], na_position="last")
        return output.drop_duplicates(subset=["event_matching_key"], keep="first").reset_index(drop=True)

    def _build_cross_comparison_as_interactions(
        self,
        pair: ComparisonPairDefinition,
        *,
        abs_dpsi_cutoff: float,
        large_delta_dpsi_cutoff: float,
        min_total_reads: int,
    ) -> dict[str, pd.DataFrame]:
        comparison_a = self._cross_as_interaction_source_frame(pair.comparison_a, abs_dpsi_cutoff, min_total_reads)
        comparison_b = self._cross_as_interaction_source_frame(pair.comparison_b, abs_dpsi_cutoff, min_total_reads)
        if comparison_a.empty and comparison_b.empty:
            return self._empty_cross_as_interaction_tables()
        n_significant_events_a = len(comparison_a)
        n_significant_events_b = len(comparison_b)

        keys_a = set(comparison_a.get("event_matching_key", pd.Series(dtype="object")).dropna().astype(str))
        keys_b = set(comparison_b.get("event_matching_key", pd.Series(dtype="object")).dropna().astype(str))
        a_only_keys = keys_a - keys_b
        b_only_keys = keys_b - keys_a
        if a_only_keys and pair.comparison_b:
            counterpart_b = self._cross_as_interaction_source_frame(
                pair.comparison_b,
                abs_dpsi_cutoff,
                min_total_reads,
                significant_only=False,
                candidate_keys=a_only_keys,
            )
            if not counterpart_b.empty:
                comparison_b = (
                    pd.concat([comparison_b, counterpart_b], axis=0, ignore_index=True)
                    .sort_values(["significant", "FDR"], ascending=[False, True], na_position="last")
                    .drop_duplicates(subset=["event_matching_key"], keep="first")
                    .reset_index(drop=True)
                )
        if b_only_keys and pair.comparison_a:
            counterpart_a = self._cross_as_interaction_source_frame(
                pair.comparison_a,
                abs_dpsi_cutoff,
                min_total_reads,
                significant_only=False,
                candidate_keys=b_only_keys,
            )
            if not counterpart_a.empty:
                comparison_a = (
                    pd.concat([comparison_a, counterpart_a], axis=0, ignore_index=True)
                    .sort_values(["significant", "FDR"], ascending=[False, True], na_position="last")
                    .drop_duplicates(subset=["event_matching_key"], keep="first")
                    .reset_index(drop=True)
                )

        merged = comparison_a.merge(comparison_b, on="event_matching_key", how="outer", suffixes=("_A", "_B"))
        if merged.empty:
            return self._empty_cross_as_interaction_tables()

        def coalesce(left: str, right: str, default: object = "") -> pd.Series:
            return (
                merged.get(left, pd.Series(pd.NA, index=merged.index))
                .combine_first(merged.get(right, pd.Series(pd.NA, index=merged.index)))
                .fillna(default)
            )

        merged["gene"] = coalesce("gene_A", "gene_B", "").astype(str)
        merged["event_type"] = coalesce("event_type_A", "event_type_B", "").astype(str)
        merged["event_id"] = coalesce("event_id_A", "event_id_B", "").astype(str)
        merged["comparison_A"] = pair.comparison_a or ""
        merged["comparison_B"] = pair.comparison_b or ""
        merged["experiment_A_display_name"] = (
            comparison_a.get("comparison_display_name", pd.Series(dtype="object")).iloc[0] if not comparison_a.empty else str(pair.comparison_a or "")
        )
        merged["experiment_B_display_name"] = (
            comparison_b.get("comparison_display_name", pd.Series(dtype="object")).iloc[0] if not comparison_b.empty else str(pair.comparison_b or "")
        )
        for column in ("dPSI_A", "dPSI_B", "FDR_A", "FDR_B", "effect_logit_A", "effect_logit_B", "effect_se_A", "effect_se_B", "mean_total_reads_A", "mean_total_reads_B", "replicate_consistency_score_A", "replicate_consistency_score_B"):
            merged[column] = pd.to_numeric(merged.get(column), errors="coerce")
        merged["significant_A"] = merged.get("significant_A", pd.Series(False, index=merged.index)).map(
            lambda value: bool(value) if pd.notna(value) else False
        )
        merged["significant_B"] = merged.get("significant_B", pd.Series(False, index=merged.index)).map(
            lambda value: bool(value) if pd.notna(value) else False
        )
        merged["direction_A"] = merged.get("direction_A", pd.Series("", index=merged.index)).fillna("").astype(str)
        merged["direction_B"] = merged.get("direction_B", pd.Series("", index=merged.index)).fillna("").astype(str)
        merged["replicate_count_A"] = pd.to_numeric(merged.get("replicate_count_A"), errors="coerce").fillna(0).astype(int)
        merged["replicate_count_B"] = pd.to_numeric(merged.get("replicate_count_B"), errors="coerce").fillna(0).astype(int)
        merged["delta_dPSI"] = merged["dPSI_A"] - merged["dPSI_B"]
        merged["abs_delta_dPSI"] = merged["delta_dPSI"].abs()
        sign_a = np.sign(merged["dPSI_A"].fillna(0))
        sign_b = np.sign(merged["dPSI_B"].fillna(0))
        merged["same_direction"] = (sign_a == sign_b) & sign_a.ne(0) & sign_b.ne(0)
        merged["opposite_direction"] = sign_a.ne(sign_b) & sign_a.ne(0) & sign_b.ne(0)
        merged["same_direction_strength_delta"] = (merged["dPSI_A"].abs() - merged["dPSI_B"].abs()).abs()
        merged["fallback_interaction_log_odds"] = merged["effect_logit_A"] - merged["effect_logit_B"]
        merged["fallback_interaction_se"] = np.sqrt(
            pd.to_numeric(merged["effect_se_A"], errors="coerce").pow(2)
            + pd.to_numeric(merged["effect_se_B"], errors="coerce").pow(2)
        )
        merged["fallback_interaction_z"] = merged["fallback_interaction_log_odds"] / merged["fallback_interaction_se"].replace(0, np.nan)
        merged["fallback_interaction_pvalue"] = merged["fallback_interaction_z"].map(self._normal_two_sided_pvalue)
        merged["interaction_log_odds"] = merged["fallback_interaction_log_odds"]
        merged["interaction_se"] = merged["fallback_interaction_se"]
        merged["interaction_z"] = merged["fallback_interaction_z"]
        merged["interaction_pvalue"] = merged["fallback_interaction_pvalue"]
        merged["interaction_model"] = "fallback_logit_difference"
        merged["interaction_model_status"] = "fallback_before_binomial_glm"
        merged["interaction_overdispersion_scale"] = np.nan
        merged["replicate_consistency_score"] = pd.concat(
            [
                merged["replicate_consistency_score_A"],
                merged["replicate_consistency_score_B"],
            ],
            axis=1,
        ).min(axis=1, skipna=True)
        merged["read_support_class"] = np.where(
            merged["mean_total_reads_A"].ge(min_total_reads) & merged["mean_total_reads_B"].ge(min_total_reads),
            "adequate_read_support",
            "low_power_or_missing_read_support",
        )
        merged["abs_dPSI_significance_cutoff_used"] = abs_dpsi_cutoff
        merged["large_delta_dPSI_cutoff_used"] = large_delta_dpsi_cutoff
        merged["min_total_reads_used"] = min_total_reads

        shared = merged.loc[merged["significant_A"] & merged["significant_B"]].copy()
        if not shared.empty:
            glm_results = shared.apply(lambda row: pd.Series(self._fit_binomial_glm_interaction(row)), axis=1)
            for column in glm_results.columns:
                shared[column] = glm_results[column]
            fallback_mask = shared["interaction_pvalue"].isna()
            if fallback_mask.any():
                shared.loc[fallback_mask, "interaction_log_odds"] = shared.loc[fallback_mask, "fallback_interaction_log_odds"]
                shared.loc[fallback_mask, "interaction_se"] = shared.loc[fallback_mask, "fallback_interaction_se"]
                shared.loc[fallback_mask, "interaction_z"] = shared.loc[fallback_mask, "fallback_interaction_z"]
                shared.loc[fallback_mask, "interaction_pvalue"] = shared.loc[fallback_mask, "fallback_interaction_pvalue"]
                shared.loc[fallback_mask, "interaction_model"] = "fallback_logit_difference"
                shared.loc[fallback_mask, "interaction_model_status"] = (
                    "binomial_glm_unavailable; fallback to replicate logit difference"
                )
                shared.loc[fallback_mask, "interaction_overdispersion_scale"] = np.nan
                shared.loc[fallback_mask, "quasi_binomial_dispersion"] = np.nan
                shared.loc[fallback_mask, "quasi_binomial_overdispersed"] = False
                shared.loc[fallback_mask, "beta_binomial_rho_estimate"] = np.nan
                shared.loc[fallback_mask, "model_valid"] = False
                shared.loc[fallback_mask, "model_diagnostic_class"] = "fallback_logit_difference"
                shared.loc[fallback_mask, "model_warning"] = "Binomial GLM unavailable; using replicate logit-difference fallback."
            beta_results = shared.apply(lambda row: pd.Series(self._fit_beta_binomial_interaction(row)), axis=1)
            for column in beta_results.columns:
                shared[column] = beta_results[column]
        shared["interaction_FDR"] = self._bh_adjust(shared["interaction_pvalue"]) if not shared.empty else pd.Series(dtype="float64")
        shared["beta_binomial_FDR"] = self._bh_adjust(shared["beta_binomial_pvalue"]) if not shared.empty and "beta_binomial_pvalue" in shared.columns else pd.Series(dtype="float64")
        if not shared.empty:
            effect = pd.to_numeric(shared["interaction_log_odds"], errors="coerce")
            se = pd.to_numeric(shared["interaction_se"], errors="coerce")
            valid = effect.notna() & se.notna() & se.gt(0)
            if valid.sum() > 1:
                observed_variance = float(effect.loc[valid].var(ddof=1))
                # A few near-separated GLMs can have enormous SEs; the median gives a
                # robust empirical-Bayes prior scale without letting those dominate.
                sampling_variance = float(se.loc[valid].pow(2).median())
                prior_variance = max(observed_variance - sampling_variance, 0.0)
            else:
                prior_variance = 0.0
            if prior_variance > 0:
                shrinkage_weight = prior_variance / (prior_variance + se.pow(2))
            else:
                shrinkage_weight = pd.Series(0.0, index=shared.index, dtype="float64")
            shared["eb_prior_variance"] = prior_variance
            shared["shrinkage_weight"] = shrinkage_weight.clip(lower=0.0, upper=1.0)
            shared["shrunken_interaction_log_odds"] = effect * shared["shrinkage_weight"]
            shared["shrunken_interaction_z"] = shared["shrunken_interaction_log_odds"] / se.replace(0, np.nan)
            shared["shrunken_interaction_pvalue"] = shared["shrunken_interaction_z"].map(self._normal_two_sided_pvalue)
            shared["shrunken_interaction_FDR"] = self._bh_adjust(shared["shrunken_interaction_pvalue"])
            effect_component = (shared["abs_delta_dPSI"] / max(large_delta_dpsi_cutoff, 1e-6)).clip(lower=0.0, upper=3.0) / 3.0
            confidence_component = (
                -np.log10(pd.to_numeric(shared["shrunken_interaction_pvalue"], errors="coerce").clip(lower=1e-300))
            ).clip(lower=0.0, upper=5.0) / 5.0
            consistency_component = pd.to_numeric(shared["replicate_consistency_score"], errors="coerce").fillna(0.0).clip(lower=0.0, upper=1.0)
            read_support_component = (
                pd.concat(
                    [
                        pd.to_numeric(shared["mean_total_reads_A"], errors="coerce"),
                        pd.to_numeric(shared["mean_total_reads_B"], errors="coerce"),
                    ],
                    axis=1,
                )
                .min(axis=1, skipna=True)
                .fillna(0.0)
                .div(max(float(min_total_reads), 1.0) * 5.0)
                .clip(lower=0.0, upper=1.0)
            )
            shared["interaction_priority_score"] = (
                100.0
                * (
                    0.35 * effect_component
                    + 0.25 * confidence_component
                    + 0.20 * consistency_component
                    + 0.20 * read_support_component
                )
            ).round(3)
            simple_delta_score = pd.to_numeric(shared["abs_delta_dPSI"], errors="coerce")
            shared["simple_delta_score"] = simple_delta_score
            shared["simple_delta_priority_score"] = (
                100.0
                * (simple_delta_score / max(float(large_delta_dpsi_cutoff), 1e-6)).clip(lower=0.0, upper=3.0)
                / 3.0
            ).round(3)
            shared["delta_dPSI_rank"] = simple_delta_score.rank(method="min", ascending=False, na_option="bottom").astype("Int64")
            shared["interaction_model_rank"] = pd.to_numeric(
                shared["interaction_priority_score"], errors="coerce"
            ).rank(method="min", ascending=False, na_option="bottom").astype("Int64")
            shared["delta_vs_model_rank_shift"] = (
                pd.to_numeric(shared["delta_dPSI_rank"], errors="coerce")
                - pd.to_numeric(shared["interaction_model_rank"], errors="coerce")
            )
            top_n = min(20, len(shared))
            shared["simple_delta_top20"] = pd.to_numeric(shared["delta_dPSI_rank"], errors="coerce").le(top_n)
            shared["model_top20"] = pd.to_numeric(shared["interaction_model_rank"], errors="coerce").le(top_n)
            shared["benchmark_score_delta"] = shared["interaction_priority_score"] - shared["simple_delta_priority_score"]
            shared["benchmark_class"] = np.select(
                [
                    shared["simple_delta_top20"] & shared["model_top20"],
                    shared["simple_delta_top20"] & ~shared["model_top20"],
                    ~shared["simple_delta_top20"] & shared["model_top20"],
                ],
                [
                    "top_by_both_delta_and_model",
                    "delta_only_top20",
                    "model_only_top20",
                ],
                default="not_top20",
            )
        else:
            for column in (
                "eb_prior_variance",
                "shrinkage_weight",
                "shrunken_interaction_log_odds",
                "shrunken_interaction_z",
                "shrunken_interaction_pvalue",
                "shrunken_interaction_FDR",
                "interaction_priority_score",
                "simple_delta_score",
                "simple_delta_priority_score",
                "delta_dPSI_rank",
                "interaction_model_rank",
                "delta_vs_model_rank_shift",
                "benchmark_score_delta",
            ):
                shared[column] = pd.Series(dtype="float64")
            for column in ("simple_delta_top20", "model_top20"):
                shared[column] = pd.Series(dtype="bool")
            shared["benchmark_class"] = pd.Series(dtype="object")
        shared["interaction_class"] = np.where(
            shared["shrunken_interaction_FDR"].le(0.05) & shared["abs_delta_dPSI"].ge(large_delta_dpsi_cutoff) & shared["read_support_class"].eq("adequate_read_support"),
            "high_confidence_interaction_candidate",
            np.where(
                shared["shrunken_interaction_pvalue"].le(0.05) & shared["abs_delta_dPSI"].ge(large_delta_dpsi_cutoff),
                "moderate_interaction_candidate",
                "pattern_followup_only",
            ),
        )
        shared["event_class"] = "shared_significant_interaction_tested"
        opposite = shared.loc[shared["opposite_direction"]].copy()
        opposite["event_class"] = "opposite_direction_shared_significant"
        large_delta = shared.loc[
            shared["same_direction"] & shared["same_direction_strength_delta"].ge(large_delta_dpsi_cutoff)
        ].copy()
        large_delta["event_class"] = "same_direction_large_dpsi_delta"
        benchmark = shared.sort_values(
            ["benchmark_class", "interaction_model_rank", "delta_dPSI_rank"],
            ascending=[True, True, True],
            na_position="last",
        ).copy()

        exclusive = merged.loc[
            (merged["significant_A"] & ~merged["significant_B"])
            | (~merged["significant_A"] & merged["significant_B"])
        ].copy()
        if not exclusive.empty:
            exclusive["interaction_FDR"] = np.nan
            exclusive["interaction_class"] = "not_interaction_tested_exclusive_event"
            exclusive["exclusive_side"] = np.where(exclusive["significant_A"], "comparison_A_only", "comparison_B_only")
            exclusive["counterpart_event_found"] = np.where(
                exclusive["significant_A"],
                exclusive["dPSI_B"].notna(),
                exclusive["dPSI_A"].notna(),
            )
            exclusive["significant_side_mean_total_reads"] = np.where(
                exclusive["significant_A"],
                exclusive["mean_total_reads_A"],
                exclusive["mean_total_reads_B"],
            )
            exclusive["non_significant_side_mean_total_reads"] = np.where(
                exclusive["significant_A"],
                exclusive["mean_total_reads_B"],
                exclusive["mean_total_reads_A"],
            )
            exclusive["significant_side_replicate_count"] = np.where(
                exclusive["significant_A"],
                exclusive["replicate_count_A"],
                exclusive["replicate_count_B"],
            )
            exclusive["non_significant_side_replicate_count"] = np.where(
                exclusive["significant_A"],
                exclusive["replicate_count_B"],
                exclusive["replicate_count_A"],
            )
            exclusive["non_significant_side_dPSI"] = np.where(
                exclusive["significant_A"],
                exclusive["dPSI_B"],
                exclusive["dPSI_A"],
            )
            exclusive["non_significant_side_FDR"] = np.where(
                exclusive["significant_A"],
                exclusive["FDR_B"],
                exclusive["FDR_A"],
            )
            exclusive["non_significant_side_read_support"] = np.select(
                [
                    exclusive["significant_A"] & ~exclusive["counterpart_event_found"],
                    ~exclusive["significant_A"] & ~exclusive["counterpart_event_found"],
                    exclusive["significant_A"] & exclusive["non_significant_side_mean_total_reads"].astype(float).ge(min_total_reads),
                    ~exclusive["significant_A"] & exclusive["non_significant_side_mean_total_reads"].astype(float).ge(min_total_reads),
                    exclusive["significant_A"],
                    ~exclusive["significant_A"],
                ],
                [
                    "counterpart_event_not_found_in_B",
                    "counterpart_event_not_found_in_A",
                    "supported_non_significant_in_B",
                    "supported_non_significant_in_A",
                    "low_power_or_absent_in_B",
                    "low_power_or_absent_in_A",
                ],
                default="not_classified",
            )
            exclusive["event_class"] = exclusive["exclusive_side"] + "|" + exclusive["non_significant_side_read_support"]

        output_columns = self._empty_cross_as_interaction_tables()["shared_interaction"].columns.tolist()
        model_status = shared.get("interaction_model_status", pd.Series(dtype="object")).fillna("").astype(str) if not shared.empty else pd.Series(dtype="object")
        diagnostic_class = shared.get("model_diagnostic_class", pd.Series(dtype="object")).fillna("").astype(str) if not shared.empty else pd.Series(dtype="object")
        top20_overlap = int((shared.get("simple_delta_top20", pd.Series(dtype=bool)) & shared.get("model_top20", pd.Series(dtype=bool))).sum()) if not shared.empty else 0
        top20_union = int((shared.get("simple_delta_top20", pd.Series(dtype=bool)) | shared.get("model_top20", pd.Series(dtype=bool))).sum()) if not shared.empty else 0
        rank_corr = float("nan")
        if not shared.empty and len(shared) > 1:
            rank_frame = shared[["delta_dPSI_rank", "interaction_model_rank"]].apply(pd.to_numeric, errors="coerce").dropna()
            if len(rank_frame) > 1:
                rank_corr = float(rank_frame["delta_dPSI_rank"].corr(rank_frame["interaction_model_rank"], method="spearman"))
        diagnostics = pd.DataFrame(
            [
                {
                    "pair_id": pair.pair_id,
                    "pair_name": pair.resolved_name,
                    "comparison_A": pair.comparison_a,
                    "comparison_B": pair.comparison_b,
                    "n_events_A": n_significant_events_a,
                    "n_events_B": n_significant_events_b,
                    "n_events_A_with_counterparts": len(comparison_a),
                    "n_events_B_with_counterparts": len(comparison_b),
                    "n_shared_significant": len(shared),
                    "n_opposite_direction": len(opposite),
                    "n_large_delta": len(large_delta),
                    "n_exclusive_supported": int(exclusive.get("non_significant_side_read_support", pd.Series(dtype="object")).astype(str).str.startswith("supported_").sum()) if not exclusive.empty else 0,
                    "n_exclusive_low_power": int(exclusive.get("non_significant_side_read_support", pd.Series(dtype="object")).astype(str).str.startswith("low_power").sum()) if not exclusive.empty else 0,
                    "n_exclusive_counterpart_not_found": int(exclusive.get("non_significant_side_read_support", pd.Series(dtype="object")).astype(str).str.startswith("counterpart_event_not_found").sum()) if not exclusive.empty else 0,
                    "n_exclusive_not_shared_significant": len(exclusive),
                    "n_model_converged": int(model_status.eq("converged").sum()),
                    "n_model_fallback": int(diagnostic_class.eq("fallback_logit_difference").sum()),
                    "n_model_failed": int(diagnostic_class.eq("model_failed").sum()),
                    "median_quasi_binomial_dispersion": float(pd.to_numeric(shared.get("quasi_binomial_dispersion", pd.Series(dtype="float64")), errors="coerce").median()) if not shared.empty else float("nan"),
                    "max_quasi_binomial_dispersion": float(pd.to_numeric(shared.get("quasi_binomial_dispersion", pd.Series(dtype="float64")), errors="coerce").max()) if not shared.empty else float("nan"),
                    "median_beta_binomial_rho_estimate": float(pd.to_numeric(shared.get("beta_binomial_rho_estimate", pd.Series(dtype="float64")), errors="coerce").median()) if not shared.empty else float("nan"),
                    "median_beta_binomial_rho": float(pd.to_numeric(shared.get("beta_binomial_rho", pd.Series(dtype="float64")), errors="coerce").median()) if not shared.empty else float("nan"),
                    "n_beta_binomial_converged": int(shared.get("beta_binomial_model_status", pd.Series(dtype="object")).astype(str).eq("converged").sum()) if not shared.empty else 0,
                    "n_beta_binomial_failed": int(shared.get("beta_binomial_model_status", pd.Series(dtype="object")).astype(str).isin(["failed", "fit_warning"]).sum()) if not shared.empty else 0,
                    "n_beta_binomial_FDR_0_05": int(pd.to_numeric(shared.get("beta_binomial_FDR", pd.Series(dtype="float64")), errors="coerce").le(0.05).sum()) if not shared.empty else 0,
                    "n_overdispersed_events": int(shared.get("quasi_binomial_overdispersed", pd.Series(dtype=bool)).map(self._safe_bool).sum()) if not shared.empty else 0,
                    "delta_top20_model_top20_overlap": top20_overlap,
                    "delta_top20_model_top20_overlap_fraction": float(top20_overlap / top20_union) if top20_union else float("nan"),
                    "delta_vs_model_spearman_rank_correlation": rank_corr,
                    "n_delta_only_top20": int(shared.get("benchmark_class", pd.Series(dtype="object")).astype(str).eq("delta_only_top20").sum()) if not shared.empty else 0,
                    "n_model_only_top20": int(shared.get("benchmark_class", pd.Series(dtype="object")).astype(str).eq("model_only_top20").sum()) if not shared.empty else 0,
                    "abs_dPSI_significance_cutoff_used": abs_dpsi_cutoff,
                    "large_delta_dPSI_cutoff_used": large_delta_dpsi_cutoff,
                    "min_total_reads_used": min_total_reads,
                }
            ]
        )
        unified_parts: list[pd.DataFrame] = []
        if not shared.empty:
            shared_unified = shared.copy()
            shared_unified["unified_result_type"] = "shared_significant_modeled"
            shared_unified["unified_display_class"] = "Shared significant modeled"
            shared_unified.loc[shared_unified.get("opposite_direction", pd.Series(False, index=shared_unified.index)).map(self._safe_bool), "unified_display_class"] = "Opposite direction"
            shared_unified.loc[
                shared_unified.get("same_direction", pd.Series(False, index=shared_unified.index)).map(self._safe_bool)
                & pd.to_numeric(shared_unified.get("same_direction_strength_delta", pd.Series(np.nan, index=shared_unified.index)), errors="coerce").ge(large_delta_dpsi_cutoff),
                "unified_display_class",
            ] = "Same direction + large dPSI delta"
            shared_unified["support_class"] = shared_unified.get("read_support_class", "shared_significant")
            unified_parts.append(shared_unified)
        if not exclusive.empty:
            exclusive_unified = exclusive.copy()
            side = exclusive_unified.get("exclusive_side", pd.Series("", index=exclusive_unified.index)).fillna("").astype(str)
            support = exclusive_unified.get("non_significant_side_read_support", pd.Series("", index=exclusive_unified.index)).fillna("").astype(str)
            exclusive_unified["unified_result_type"] = np.where(
                side.eq("comparison_A_only"),
                "A_specific_significant",
                np.where(side.eq("comparison_B_only"), "B_specific_significant", "exclusive_significant"),
            )
            exclusive_unified["unified_display_class"] = np.where(
                side.eq("comparison_A_only"),
                "A-specific significant",
                np.where(side.eq("comparison_B_only"), "B-specific significant", "Specific significant"),
            )
            exclusive_unified["support_class"] = support
            unified_parts.append(exclusive_unified)
        if unified_parts:
            unified = pd.concat(unified_parts, ignore_index=True, sort=False)
            unified["_sort_priority"] = unified["unified_result_type"].map(
                {
                    "shared_significant_modeled": 0,
                    "A_specific_significant": 1,
                    "B_specific_significant": 2,
                    "exclusive_significant": 3,
                }
            ).fillna(9)
            unified = unified.sort_values(
                ["_sort_priority", "interaction_priority_score", "abs_delta_dPSI", "gene"],
                ascending=[True, False, False, True],
                na_position="last",
            ).drop(columns=["_sort_priority"], errors="ignore")
        else:
            unified = pd.DataFrame(columns=output_columns)
        return {
            "unified_event_table": unified.reindex(columns=output_columns).reset_index(drop=True),
            "shared_interaction": shared.reindex(columns=output_columns).reset_index(drop=True),
            "opposite_direction": opposite.reindex(columns=output_columns).reset_index(drop=True),
            "large_delta": large_delta.reindex(columns=output_columns).reset_index(drop=True),
            "exclusive": exclusive.reindex(columns=output_columns).reset_index(drop=True),
            "benchmark": benchmark.reindex(columns=output_columns).reset_index(drop=True),
            "diagnostics": diagnostics.reset_index(drop=True),
        }

    def _resolve_genome_fasta_path(self) -> Path | None:
        if self.current_project is None:
            return None
        for key in ("genome_fasta", "genome", "fasta", "primary_assembly_fasta"):
            value = self.current_project.input_paths.get(key) if self.current_project.input_paths else None
            if value and Path(value).exists():
                return Path(value)
        roots = [self.current_project.project_root, *self.current_project.project_root.parents]
        names = (
            "GRCm39.primary_assembly.genome.fa",
            "Mus_musculus.GRCm39.dna.primary_assembly.fa",
            "genome.fa",
        )
        for root in roots[:5]:
            for name in names:
                candidate = root / name
                if candidate.exists():
                    return candidate
        return None

    @staticmethod
    def _load_fasta_index(fasta_path: Path) -> dict[str, tuple[int, int, int, int]]:
        index_path = Path(str(fasta_path) + ".fai")
        if not index_path.exists():
            return {}
        index: dict[str, tuple[int, int, int, int]] = {}
        with open(index_path, "r", encoding="utf-8") as handle:
            for line in handle:
                parts = line.rstrip("\n").split("\t")
                if len(parts) >= 5:
                    index[parts[0]] = (int(parts[1]), int(parts[2]), int(parts[3]), int(parts[4]))
        return index

    @staticmethod
    def _fetch_fasta_sequence(fasta_path: Path, index: dict[str, tuple[int, int, int, int]], chrom: str, start: int, end: int) -> str:
        if chrom not in index or start > end:
            return ""
        length, offset, line_bases, line_width = index[chrom]
        start = max(1, int(start))
        end = min(int(end), int(length))
        if start > end:
            return ""
        chunks: list[str] = []
        with open(fasta_path, "rb") as handle:
            pos = start
            while pos <= end:
                zero_based = pos - 1
                line_index = zero_based // line_bases
                line_offset = zero_based % line_bases
                read_bases = min(end - pos + 1, line_bases - line_offset)
                handle.seek(offset + line_index * line_width + line_offset)
                chunks.append(handle.read(read_bases).decode("ascii", errors="ignore"))
                pos += read_bases
        return "".join(chunks).upper()

    @staticmethod
    def _reverse_complement(seq: str) -> str:
        return seq.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1].upper()

    @staticmethod
    def _parse_event_region(coordinates: object) -> tuple[str, str, int, int] | None:
        text = str(coordinates or "")
        chrom_match = re.search(r"chr=([^|\s]+)", text)
        strand_match = re.search(r"strand=([+-])", text)
        numbers = [int(item) for item in re.findall(r"=(\d+)", text)]
        if not chrom_match or not numbers:
            return None
        return chrom_match.group(1), strand_match.group(1) if strand_match else "+", min(numbers) + 1, max(numbers)

    @staticmethod
    def _kmer_counts(seq: str, k: int) -> dict[str, int]:
        counts: dict[str, int] = {}
        if len(seq) < k:
            return counts
        for idx in range(0, len(seq) - k + 1):
            kmer = seq[idx : idx + k].upper()
            if len(kmer) == k and set(kmer) <= {"A", "C", "G", "T"}:
                counts[kmer] = counts.get(kmer, 0) + 1
        return counts

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
                    self._safe_text(row.get("normalized_gene_id")) or self._safe_text(row.get("gene")),
                    self._safe_text(row.get("event_type")),
                    self._safe_text(row.get("event_id")) or self._safe_text(row.get("coordinates")),
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
        stable_source = (
            working["stable_as_event"].fillna(False).map(self._safe_bool)
            if "stable_as_event" in working.columns
            else working["FDR"].le(self.current_thresholds.splicing_fdr)
        )
        working["significant"] = stable_source & working["dPSI"].abs().ge(abs_dpsi_cutoff)
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
                "mean_junction_reads_group1",
                "mean_junction_reads_group2",
                "replicate_count_group1",
                "replicate_count_group2",
                "supported_replicate_count_group1",
                "supported_replicate_count_group2",
                "replicate_direction_consistency",
                "stable_support_status",
            ]
        ].copy()

    def _build_cross_comparison_significant_as_patterns(
        self,
        pair: ComparisonPairDefinition,
        *,
        abs_dpsi_cutoff: float,
        large_delta_dpsi_cutoff: float,
    ) -> dict[str, pd.DataFrame]:
        import numpy as np

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
        merged["significant_A"] = merged.get(
            "significant_A",
            pd.Series(False, index=merged.index, dtype=bool),
        ).map(self._safe_bool)
        merged["significant_B"] = merged.get(
            "significant_B",
            pd.Series(False, index=merged.index, dtype=bool),
        ).map(self._safe_bool)
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
        if self._safe_bool(row.get("DEG_significant")):
            labels.append("DEG")
        if self._safe_bool(row.get("rMATS_significant")):
            labels.append("rMATS")
        if self._safe_bool(row.get("DEXSeq_significant")):
            labels.append("DEXSeq")
        if self._safe_bool(row.get("DTU_significant")):
            labels.append("DTU")
        return "_".join(labels) if labels else "unclassified"

    def _candidate_direction_class(self, row: pd.Series) -> str:
        deg_sig = self._safe_bool(row.get("DEG_significant"))
        rmats_sig = self._safe_bool(row.get("rMATS_significant"))
        dex_sig = self._safe_bool(row.get("DEXSeq_significant"))
        dtu_sig = self._safe_bool(row.get("DTU_significant"))
        support = dex_sig or dtu_sig
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
        if deg_sig and not rmats_sig and support:
            return "DEG_DEXSeq_or_DTU_support"
        if deg_sig and not rmats_sig:
            return "DEG_only"
        if rmats_sig and not deg_sig and support:
            return "rMATS_DEXSeq_or_DTU_only"
        if rmats_sig and not deg_sig:
            return "rMATS_only"
        if (not deg_sig) and (not rmats_sig) and dex_sig and dtu_sig:
            return "DEXSeq_DTU_support_only"
        if (not deg_sig) and (not rmats_sig) and dex_sig:
            return "DEXSeq_support_only"
        if (not deg_sig) and (not rmats_sig) and dtu_sig:
            return "DTU_support_only"
        return "unclassified"

    def _candidate_tier(self, row: pd.Series) -> str:
        deg_sig = self._safe_bool(row.get("DEG_significant"))
        rmats_sig = self._safe_bool(row.get("rMATS_significant"))
        dex_sig = self._safe_bool(row.get("DEXSeq_significant"))
        dtu_sig = self._safe_bool(row.get("DTU_significant"))
        support = dex_sig or dtu_sig
        if deg_sig and rmats_sig and support:
            return "Tier 1"
        if deg_sig and rmats_sig and not support:
            return "Tier 2"
        if ((not deg_sig) and rmats_sig and support) or (deg_sig and (not rmats_sig) and support):
            return "Tier 3"
        if ((not deg_sig) and rmats_sig and not support) or ((not deg_sig) and (not rmats_sig) and support):
            return "Tier 4"
        if deg_sig and (not rmats_sig) and (not support):
            return "Tier 5"
        return "Unclassified"

    def _candidate_reason(self, row: pd.Series) -> str:
        tier = self._safe_text(row.get("candidate_tier", ""))
        evidence = self._safe_text(row.get("evidence_class", ""), "unclassified")
        if tier == "Tier 1":
            return "Tier 1: DEG significant, rMATS significant, and DEXSeq/DTU support. Strong multi-layer candidate."
        if tier == "Tier 2":
            return "Tier 2: DEG significant and rMATS significant. Core expression-splicing candidate."
        if tier == "Tier 3":
            return f"Tier 3: support-backed single-primary candidate ({evidence}). Includes rMATS+support or DEG+support genes."
        if tier == "Tier 4":
            return f"Tier 4: single-layer splicing/support candidate ({evidence}). Includes rMATS-only or DEXSeq/DTU-only support genes."
        if tier == "Tier 5":
            return "Tier 5: DEG significant only. Expression-only candidate."
        return f"Unclassified candidate ({evidence})."

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
            deg_group1, deg_group2 = self._comparison_deg_direction_labels(comparison)
            as_group1, as_group2 = self._comparison_splicing_direction_labels(comparison)
            output_prefix = comparison.output_prefix or self._safe_output_name(comparison.comparison_id)
            payload[comparison.comparison_id] = {
                "display_name": comparison.display_resolved_name,
                "biological_question": comparison.biological_question or "",
                "group1": group1,
                "group2": group2,
                "transcript_contrast": {
                    "numerator_group": deg_group1,
                    "denominator_group": deg_group2,
                    "log2fc_positive_means": f"Transcript expression higher in {deg_group1} than {deg_group2}",
                },
                "splicing_contrast": {
                    "group1": as_group1,
                    "group2": as_group2,
                    "dpsi_positive_means": f"Event usage higher in {as_group1} than {as_group2}",
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
        deg_group1, deg_group2 = self._comparison_deg_direction_labels(comparison)
        as_group1, as_group2 = self._comparison_splicing_direction_labels(comparison)
        return {
            "comparison_id": comparison.comparison_id,
            "display_name": comparison.display_resolved_name,
            "biological_question": comparison.biological_question or "",
            "group1": group1,
            "group2": group2,
            "transcript_contrast": {
                "numerator_group": deg_group1,
                "denominator_group": deg_group2,
            },
            "splicing_contrast": {
                "group1": as_group1,
                "group2": as_group2,
            },
            "log2fc_direction_label": f"log2FC > 0 means transcript expression is higher in {deg_group1} than {deg_group2}",
            "dpsi_direction_label": f"dPSI > 0 means event usage is higher in {as_group1} than {as_group2}",
            "transcript_direction_flipped": bool(comparison.reverse_deg),
            "splicing_direction_flipped": bool(comparison.reverse_splicing),
            "thresholds": {
                "DE_FDR": self.current_thresholds.deg_padj,
                "abs_standardized_log2FC": self.current_thresholds.deg_log2fc,
                "rMATS_FDR": self.current_thresholds.splicing_fdr,
                "abs_standardized_dPSI": self.current_thresholds.splicing_dpsi,
                "stable_AS_min_IJC_plus_SJC_per_replicate": self.current_thresholds.stable_min_junction_reads,
                "stable_AS_min_replicates_per_group": self.current_thresholds.stable_min_replicates_per_group,
                "stable_AS_direction_consistency": self.current_thresholds.stable_min_pairwise_consistency,
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
                    "stable_AS_min_IJC_plus_SJC_per_replicate": self.current_thresholds.stable_min_junction_reads,
                    "stable_AS_min_replicates_per_group": self.current_thresholds.stable_min_replicates_per_group,
                    "stable_AS_direction_consistency": self.current_thresholds.stable_min_pairwise_consistency,
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
                display_group = self._safe_text(row.get("display_group"), comparison.display_resolved_name)
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
                        "dPSI": row.get("dPSI", row.get("standardized_dPSI", "")),
                        "FDR": row.get("FDR", row.get("best_rMATS_FDR", "")),
                        "direction": row.get("direction", ""),
                        "coordinates": row.get("coordinates", ""),
                        "event_matching_key": row.get("event_matching_key", event_uid),
                        "splicing_direction_flipped": bool(comparison.reverse_splicing),
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
        is_sig_a = self._safe_bool(row.get("is_sig_A"))
        is_sig_b = self._safe_bool(row.get("is_sig_B"))

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
        if not deduped:
            fallback = self._safe_text(row.get("comparison_id")) or self._safe_text(row.get("comparison"))
            if fallback:
                deduped.append(fallback)
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
        final_exp, final_ctrl = comparison.final_as_groups
        analysis_exp = (final_exp or comparison.experiment_group or comparison.source_experiment_group or "").strip()
        analysis_ctrl = (final_ctrl or comparison.control_group or comparison.source_control_group or "").strip()
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

    @staticmethod
    def _threshold_config_from_mapping(values: dict | None) -> ThresholdConfig:
        allowed = set(ThresholdConfig.__dataclass_fields__.keys())
        cleaned = {key: value for key, value in (values or {}).items() if key in allowed}
        return ThresholdConfig(**cleaned)


def _looks_like_genome_gtf(path: Path) -> bool:
    name = path.name.lower()
    if not (name.endswith(".gtf") or name.endswith(".gtf.gz")):
        return False
    if name.startswith("events_") or "per_local_event" in str(path).lower():
        return False
    return any(token in name for token in ("gencode", "annotation", "genes"))
