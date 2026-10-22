from pathlib import Path

import pandas as pd
import pytest

from katana_cosplice.config import CoSpliceConfig, NetworkConfig, load_config
from katana_cosplice.network import fit_network
from katana_cosplice.pipeline import CoSplicePipeline


def config(tmp_path: Path, minimum=15):
    return CoSpliceConfig(
        project_name="small_sample_regression",
        rmats_root=tmp_path / "rmats",
        sample_metadata=tmp_path / "samples.tsv",
        output_dir=tmp_path / "output",
        network=NetworkConfig(minimum_network_samples=minimum),
    )


def test_sample_gate_default_and_explicit_override():
    # Too few events stops fitting after the sample gate, without running R.
    matrix = pd.DataFrame(0.5, index=["e1", "e2"], columns=[f"s{i}" for i in range(14)])
    qc = pd.DataFrame({"stability_score": [1.0, 1.0]}, index=matrix.index)
    with pytest.raises(ValueError, match="At least 15"):
        fit_network(matrix, qc, NetworkConfig())
    with pytest.warns(UserWarning, match="Exploratory network with 14 samples"):
        with pytest.raises(ValueError, match="Too few stable events"):
            fit_network(matrix, qc, NetworkConfig(minimum_network_samples=14))


@pytest.mark.parametrize("minimum", [True, 2, -1, 14.5, "14"])
def test_invalid_sample_gate_rejected_at_config_validation(tmp_path, minimum):
    with pytest.raises(ValueError, match="minimum_network_samples"):
        config(tmp_path, minimum).validate()


def test_explicit_sample_gate_survives_yaml_roundtrip(tmp_path):
    original = config(tmp_path, 14)
    path = tmp_path / "config.yaml"
    original.write(path)
    restored = load_config(path)
    assert restored.network.minimum_network_samples == 14
    assert restored.network.bootstrap_iterations == original.network.bootstrap_iterations
    assert restored.network.max_events == original.network.max_events


@pytest.mark.parametrize("with_fdr", [False, True])
def test_report_handles_unconfigured_factorial_without_false_zero(tmp_path, monkeypatch, with_fdr):
    import katana_cosplice.pipeline as module

    # Isolate HTML report behavior from PDF rendering and browser dependencies.
    monkeypatch.setattr(module, "build_module_browser", lambda *args: args[-1].write_text("browser"))
    monkeypatch.setattr(module, "build_validation_pdf", lambda *args: args[-1].write_bytes(b"pdf stub"))
    pipeline = CoSplicePipeline(config(tmp_path))
    network = pipeline.output / "04_network/observed"
    network.mkdir(parents=True)
    pd.DataFrame({"module": ["M01"], "largest_module_fraction": [0.2]}).to_csv(
        network / "network_diagnostics.tsv", sep="\t", index=False
    )
    association = pipeline.output / "05_association/observed"
    association.mkdir(parents=True)
    table = pd.DataFrame({"FDR": [0.01, 0.20, None]}) if with_fdr else pd.DataFrame({"status": ["not_configured"]})
    table.to_csv(association / "factorial_module_association.tsv", sep="\t", index=False)
    report = pipeline.report().read_text(encoding="utf-8")
    if with_fdr:
        assert "1 module-term associations pass FDR 0.05" in report
    else:
        assert "Factorial association was not assessed" in report
        assert "0 module-term associations pass" not in report
