import importlib.util
from pathlib import Path

import pytest


def test_empty_process_log_is_valid_but_empty_data_is_not(tmp_path):
    source = Path(__file__).resolve().parents[2] / "scripts/accept_rbp_upgrade.py"
    spec = importlib.util.spec_from_file_location("rbp_acceptance", source)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    log = tmp_path / "fimo.log"
    log.touch()
    module.validate_output(log)
    data = tmp_path / "enrichment.parquet"
    data.touch()
    with pytest.raises(ValueError, match="empty output"):
        module.validate_output(data)
    with pytest.raises(ValueError, match="Missing"):
        module.validate_output(tmp_path / "missing.log")
