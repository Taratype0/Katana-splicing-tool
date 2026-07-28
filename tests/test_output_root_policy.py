from __future__ import annotations

from pathlib import Path

import pytest

from src.domain.models.project import ProjectConfig
from src.io.scanner import ProjectScanner
from src.services.project_service import ProjectService


def test_new_scan_does_not_infer_output_directory(tmp_path: Path) -> None:
    project = ProjectScanner().scan(tmp_path)
    assert project.output_root is None


def test_output_root_inside_source_tree_is_rejected(tmp_path: Path) -> None:
    service = ProjectService()
    service.current_project = ProjectConfig(
        project_name="test",
        project_root=tmp_path,
        output_root=None,
    )
    with pytest.raises(RuntimeError, match="cannot be inside"):
        service.update_output_root(service.application_source_root() / "generated")


def test_selected_external_output_root_is_created_and_used(tmp_path: Path) -> None:
    project_root = tmp_path / "inputs"
    output_root = tmp_path / "results"
    project_root.mkdir()
    service = ProjectService()
    service.current_project = ProjectConfig(
        project_name="test",
        project_root=project_root,
        config_path=project_root / ".katana_project.yaml",
    )

    service.update_output_root(output_root)

    assert service.require_output_root() == output_root.resolve()
    assert service.default_export_path("table.tsv") == output_root.resolve() / "exports" / "table.tsv"
    assert output_root.exists()


def test_recognized_output_tree_can_be_migrated(tmp_path: Path) -> None:
    project_root = tmp_path / "inputs"
    old_output = tmp_path / "katana_output"
    new_output = tmp_path / "migrated_results"
    project_root.mkdir()
    (old_output / "06_native_splicing_patterns").mkdir(parents=True)
    (old_output / "06_native_splicing_patterns" / "result.tsv").write_text("x\n1\n", encoding="utf-8")
    service = ProjectService()
    service.current_project = ProjectConfig(
        project_name="test",
        project_root=project_root,
        config_path=project_root / ".katana_project.yaml",
        output_root=old_output,
    )

    result = service.migrate_output_root(new_output)

    assert result.moved_entries == 1
    assert result.source_removed is True
    assert not old_output.exists()
    assert (new_output / "06_native_splicing_patterns" / "result.tsv").exists()
    assert service.current_project.output_root == new_output.resolve()


def test_unrecognized_directory_is_never_migrated(tmp_path: Path) -> None:
    project_root = tmp_path / "inputs"
    project_root.mkdir()
    (project_root / "raw.fastq").write_text("raw", encoding="utf-8")
    service = ProjectService()
    service.current_project = ProjectConfig(
        project_name="test",
        project_root=project_root,
        config_path=project_root / ".katana_project.yaml",
        output_root=project_root,
    )

    with pytest.raises(RuntimeError, match="does not look like"):
        service.migrate_output_root(tmp_path / "results")

    assert (project_root / "raw.fastq").exists()
