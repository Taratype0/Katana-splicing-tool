"""Installed-runtime paths are portable; no native fit or server is started."""
from pathlib import Path
import sys

import pandas as pd
import pytest

from katana_next import launcher


@pytest.mark.parametrize("platform, variable, relative", [
    ("win32", "LOCALAPPDATA", Path("Katana/runtime")),
    ("linux", "XDG_DATA_HOME", Path("katana/runtime")),
])
def test_launcher_uses_absolute_user_data_environment(tmp_path, monkeypatch, platform, variable, relative):
    user_data = tmp_path / "user-data"
    monkeypatch.setattr(launcher.sys, "platform", platform)
    monkeypatch.setenv(variable, str(user_data))
    assert launcher.default_runtime_root() == user_data / relative
    assert not user_data.exists()  # Resolving defaults performs no writes.


@pytest.mark.parametrize("platform, variable, relative", [
    ("win32", "LOCALAPPDATA", Path("AppData/Local/Katana/runtime")),
    ("linux", "XDG_DATA_HOME", Path(".local/share/katana/runtime")),
    ("darwin", "XDG_DATA_HOME", Path("Library/Application Support/Katana/runtime")),
])
def test_launcher_fallback_is_user_home_not_package(tmp_path, monkeypatch, platform, variable, relative):
    monkeypatch.setattr(launcher.sys, "platform", platform)
    monkeypatch.setattr(launcher.Path, "home", classmethod(lambda cls: tmp_path / "home"))
    monkeypatch.setenv(variable, "relative-path-is-not-a-user-data-root")
    assert launcher.default_runtime_root() == tmp_path / "home" / relative


@pytest.mark.parametrize("explicit", [False, True])
def test_launcher_main_writes_only_selected_user_root(tmp_path, monkeypatch, explicit):
    from katana_next import desktop
    default = tmp_path / "user-data/runtime"
    selected = tmp_path / "explicit-runtime" if explicit else default
    monkeypatch.setattr(launcher, "default_runtime_root", lambda: default)
    monkeypatch.setattr(launcher, "__file__", str(tmp_path / "read-only-install/katana_next/launcher.py"))
    monkeypatch.setattr(launcher.faulthandler, "enable", lambda **kw: None)
    monkeypatch.setattr(launcher, "_launch_runtime", None)
    monkeypatch.setattr(desktop, "desktop_runtime_root", lambda root: root)
    seen = []
    monkeypatch.setattr(launcher, "ensure_service", lambda port, root: seen.append(root) or "http://127.0.0.1:8765/")
    monkeypatch.setattr(launcher.webbrowser, "open", lambda *a, **k: pytest.fail("No browser requested"))
    monkeypatch.setattr(launcher.sys, "argv", ["katana-open", "--no-open"] +
                        (["--runtime-root", str(selected)] if explicit else []))
    launcher.main()
    assert seen == [selected.resolve()]
    assert (selected / "launcher-events.jsonl").is_file()
    assert not (tmp_path / "read-only-install").exists()
    if explicit:
        assert not default.exists()


def test_standalone_wgcna_resolves_bare_rscript_before_input_validation(monkeypatch):
    from katana_cosplice import r_wgcna
    monkeypatch.setattr(r_wgcna.shutil, "which", lambda name: sys.executable if name == "Rscript" else None)
    monkeypatch.setattr(r_wgcna.subprocess, "run", lambda *a, **k: pytest.fail("No native calculation allowed"))
    with pytest.raises(ValueError, match="nonempty matrix"):
        r_wgcna.exact_wgcna_labels(pd.DataFrame(), soft_power=6, min_module_size=20,
                                  merge_cut_height=.25, deep_split=2,
                                  pam_respects_dendro=False, rscript="Rscript")


@pytest.mark.parametrize("parameter_selected", [False, True])
def test_flow_network_persists_resolved_rscript_before_execution(tmp_path, monkeypatch, parameter_selected):
    from katana_cosplice import config as config_module
    from katana_next.flow import stages
    from katana_next.workflows import programs
    config = config_module.CoSpliceConfig(project_name="test", rmats_root=tmp_path,
                                          sample_metadata=tmp_path / "samples.tsv", output_dir=tmp_path / "old")
    config.network.wgcna_rscript = "Rscript" if not parameter_selected else None
    monkeypatch.setattr(config_module, "load_config", lambda path: config)
    monkeypatch.setattr(config_module.CoSpliceConfig, "validate", lambda self: None)
    monkeypatch.delenv("KATANA_CONTAINER", raising=False)
    monkeypatch.setattr(stages.shutil, "which", lambda name: sys.executable if name == "Rscript" else None)
    class StopBeforeScience(Exception):
        pass
    def capture(actual, destination):
        assert actual.network.wgcna_rscript == str(Path(sys.executable).absolute())
        raise StopBeforeScience
    monkeypatch.setattr(programs, "write_production_configuration", capture)
    with pytest.raises(StopBeforeScience):
        stages._network({"adapter": "network.fit", "inputs": {"config": "declared.yaml"},
                         "parameters": {"rscript": "Rscript"} if parameter_selected else {}}, tmp_path / "stage")
    assert not (tmp_path / "stage").exists()
