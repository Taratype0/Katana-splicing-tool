"""R discovery uses explicit configuration and installed tools, never a workstation path."""
from pathlib import Path

import pytest

from katana_cosplice import runtime_discovery


@pytest.fixture
def isolated_discovery(monkeypatch, tmp_path):
    for name in ("KATANA_RSCRIPT", "CONDA_PREFIX", "CONDA_EXE", "R_HOME", "ProgramFiles", "ProgramFiles(x86)"):
        monkeypatch.delenv(name, raising=False)
    monkeypatch.setattr(runtime_discovery.shutil, "which", lambda name: None)
    monkeypatch.setattr(runtime_discovery.sys, "prefix", str(tmp_path / "python"))
    monkeypatch.setattr(Path, "home", classmethod(lambda cls: tmp_path / "home"))
    return tmp_path


def executable(path):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.touch()
    return path.resolve()


def test_explicit_setting_precedes_path_and_deduplicates(monkeypatch, isolated_discovery):
    root = isolated_discovery
    configured = executable(root / "chosen installation" / "Rscript")
    available = executable(root / "path installation" / "Rscript")
    monkeypatch.setenv("KATANA_RSCRIPT", str(configured))
    monkeypatch.setattr(runtime_discovery.shutil, "which", lambda name: str(available) if name == "Rscript" else None)
    assert runtime_discovery.rscript_candidates() == [configured, available]
    monkeypatch.setenv("KATANA_RSCRIPT", str(available))
    assert runtime_discovery.rscript_candidates() == [available]


def test_active_environment_finds_posix_and_windows_layouts(monkeypatch, isolated_discovery):
    prefix = isolated_discovery / "active environment"
    unix = executable(prefix / "bin" / "Rscript")
    windows = executable(prefix / "Scripts" / "Rscript.exe")
    sibling = executable(prefix / "envs" / "analysis" / "Lib" / "R" / "bin" / "Rscript.exe")
    monkeypatch.setenv("CONDA_PREFIX", str(prefix))
    assert runtime_discovery.rscript_candidates() == [unix, windows, sibling]


def test_r_home_is_supported(monkeypatch, isolated_discovery):
    home = isolated_discovery / "R installation"
    selected = executable(home / "bin" / "Rscript")
    monkeypatch.setenv("R_HOME", str(home))
    assert runtime_discovery.rscript_candidates() == [selected]


def test_windows_standard_installations_derive_from_environment(monkeypatch, isolated_discovery):
    installed = isolated_discovery / "program files"
    selected = executable(installed / "R" / "R-4.4.0" / "bin" / "Rscript.exe")
    monkeypatch.setattr(runtime_discovery.sys, "platform", "win32")
    monkeypatch.setenv("ProgramFiles", str(installed))
    assert runtime_discovery.rscript_candidates() == [selected]


def test_missing_r_does_not_fabricate_an_executable(monkeypatch, isolated_discovery):
    monkeypatch.setenv("KATANA_RSCRIPT", str(isolated_discovery / "missing"))
    assert runtime_discovery.rscript_candidates() == []


def test_existing_discovery_entrypoints_share_the_portable_candidates(monkeypatch, isolated_discovery):
    selected = executable(isolated_discovery / "selected Rscript")
    monkeypatch.setattr(runtime_discovery, "rscript_candidates", lambda: [selected])
    from katana_cosplice.benchmark import _rscript_candidates
    from katana_next.full_fixture import _find_rscript
    from katana_next.migration import _discover_rscript
    assert _rscript_candidates() == [selected]
    assert _find_rscript() == selected
    assert _discover_rscript(None) == str(selected)
    explicit = executable(isolated_discovery / "configured Rscript")
    assert _discover_rscript(str(explicit)) == str(explicit)
