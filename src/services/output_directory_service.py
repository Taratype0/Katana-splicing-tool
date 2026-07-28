from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import shutil


_OUTPUT_MARKERS = {
    "00_project_report",
    "01_splicing_landscape",
    "02_candidate_gene_screening",
    "03_cross_comparison_candidate_comparison",
    "04_tx_splicing_integration",
    "05_deg_explorer",
    "06_native_splicing_patterns",
    "06_sashimi",
    "07_rbp_support",
    "08_isoform_candidate_support",
    "09_cosplice_network",
}


@dataclass(frozen=True)
class OutputMigrationResult:
    source: Path
    destination: Path
    moved_entries: int
    source_removed: bool


class OutputDirectoryService:
    """Validate and migrate Katana output roots without touching input data."""

    def __init__(self, application_source_root: Path) -> None:
        self.application_source_root = application_source_root.expanduser().resolve()

    def validation_error(self, output_root: str | Path | None) -> str:
        if output_root is None or not str(output_root).strip():
            return "Output directory is not configured."
        resolved = Path(output_root).expanduser().resolve()
        if resolved == self.application_source_root or self.application_source_root in resolved.parents:
            return (
                "The output directory cannot be inside the Katana source directory. "
                "Choose a separate project results directory."
            )
        return ""

    def require(self, output_root: str | Path | None, *, create: bool = True) -> Path:
        error = self.validation_error(output_root)
        if error:
            raise RuntimeError(error)
        resolved = Path(output_root).expanduser().resolve()
        if create:
            resolved.mkdir(parents=True, exist_ok=True)
        return resolved

    def is_recognized_output_root(self, output_root: str | Path | None) -> bool:
        if output_root is None:
            return False
        source = Path(output_root).expanduser().resolve()
        if not source.is_dir() or source == self.application_source_root:
            return False
        if source.name.casefold() in {"katana_output", "katana-results", "katana_results"}:
            return True
        try:
            child_names = {child.name for child in source.iterdir() if child.is_dir()}
        except OSError:
            return False
        return bool(child_names & _OUTPUT_MARKERS)

    def migrate(
        self,
        source_root: str | Path,
        destination_root: str | Path,
        *,
        remove_source: bool = True,
    ) -> OutputMigrationResult:
        source = Path(source_root).expanduser().resolve()
        destination = self.require(destination_root)
        if source == destination:
            return OutputMigrationResult(source, destination, 0, False)
        if not self.is_recognized_output_root(source):
            raise RuntimeError(
                "The current directory does not look like a Katana output root, so it was not moved. "
                "Use Choose Output Directory to change the path without migrating input files."
            )
        if source in destination.parents or destination in source.parents:
            raise RuntimeError("The new output directory must not be nested inside the current output directory.")

        moved_entries = 0
        for child in source.iterdir():
            target = destination / child.name
            if child.is_dir():
                shutil.copytree(child, target, dirs_exist_ok=True)
            else:
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(child, target)
            moved_entries += 1

        source_removed = False
        if remove_source:
            shutil.rmtree(source)
            source_removed = True
        return OutputMigrationResult(source, destination, moved_entries, source_removed)
