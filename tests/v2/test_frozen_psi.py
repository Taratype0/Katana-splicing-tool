from dataclasses import asdict
from types import SimpleNamespace
import json
import numpy as np
import pandas as pd
import pytest

from katana_cosplice.config import NetworkConfig
from katana_cosplice.network import _program_scores
from katana_next.science.psi_network import PsiFilterConfig, psi_filter_audit
from katana_next.science.frozen_psi import freeze_psi_modules, project_psi_modules, save_frozen_model, load_frozen_model


def fixture():
    rng = np.random.default_rng(17)
    ids = [f"train{i}" for i in range(30)]
    matrix = pd.DataFrame(rng.uniform(.01, .99, (45, 30)), index=[f"E{i}" for i in range(45)], columns=ids)
    matrix.iloc[0, 0] = np.nan
    means = matrix.mean(axis=1)
    complete = matrix.T.fillna(means).T
    labels = np.array(["M01"]*20+["M02"]*20+["grey"]*5)
    baseline = SimpleNamespace(selected_matrix=complete, assignments=pd.DataFrame({"event_key":matrix.index, "module":labels}),
        program_scores=_program_scores(complete, labels))
    samples = pd.DataFrame({"sample_id":ids,"biological_sample_id":ids,"condition":"A"})
    traits = pd.DataFrame({"trait":3+2*baseline.program_scores.M01}, index=ids)
    model = freeze_psi_modules(matrix, means, baseline, samples, feature_space_id="hg19|test", network=NetworkConfig(), filters=PsiFilterConfig(), traits=traits)
    return matrix, baseline, samples, model


def test_all_qc_failures_and_old_default_selection_are_preserved():
    rng = np.random.default_rng(1)
    x = pd.DataFrame(rng.uniform(size=(50, 40)))
    x.iloc[0] = np.nan
    x.iloc[1] = .05
    x.iloc[2, :30] = np.nan
    q = psi_filter_audit(x)
    old = x.notna().mean(axis=1).ge(.8) & x.var(axis=1).ge(.0025) & (x.max(axis=1)-x.min(axis=1)).ge(.5) & (x.ge(.2)&x.le(.8)).mean(axis=1).ge(.1)
    assert q.passed.equals(old)
    assert q.failure_count.iloc[0] == 4
    assert q.filter_reasons.iloc[1] == "variance;range;alternative_fraction"
    assert q.failed_valid_fraction.iloc[2]
    assert asdict(PsiFilterConfig()) == dict(min_valid_fraction=.8, min_variance=.0025, min_range=.5, min_alternative_fraction=.1)


def test_frozen_projection_reproduces_scores_and_training_only_trait_coefficients(tmp_path):
    matrix, baseline, samples, model = fixture()
    # Identical numeric observations with distinct biological IDs permit a mathematical replay.
    new = matrix.rename(columns=lambda x:"new_"+x)
    meta = samples.assign(sample_id=new.columns, biological_sample_id=new.columns)
    result = project_psi_modules(model, new, meta, feature_space_id="hg19|test")
    assert np.allclose(result["scores"], baseline.program_scores, atol=1e-12)
    assert model["trait_models"][0]["module"] == "M01"
    assert model["trait_models"][0]["slope"] == pytest.approx(2.)
    path = tmp_path/"frozen.json"
    save_frozen_model(model, path)
    loaded = load_frozen_model(path)
    assert loaded == model
    old = path.read_bytes()
    # Reordering and addition of an irrelevant event/sample do not change existing scores.
    permuted = new.iloc[::-1, ::-1].copy()
    permuted["another"] = .75
    permuted.loc["unused"] = .9
    newmeta = pd.DataFrame({"sample_id":permuted.columns,"biological_sample_id":permuted.columns})
    second = project_psi_modules(loaded, permuted, newmeta, feature_space_id="hg19|test")
    assert np.allclose(second["scores"].loc[new.columns], result["scores"])
    assert path.read_bytes() == old


def test_missing_event_audit_coverage_and_no_zero_score_substitution():
    matrix, baseline, samples, model = fixture()
    new = matrix.rename(columns=lambda x:"new_"+x).drop(index=["E0"])
    new.loc[[f"E{i}" for i in range(1, 10)], new.columns[0]] = np.nan
    meta = samples.assign(sample_id=new.columns, biological_sample_id=new.columns)
    result = project_psi_modules(model, new, meta, feature_space_id="hg19|test")
    assert np.isnan(result["scores"].M01.iloc[0])
    assert np.isfinite(result["scores"].M02.iloc[0])
    assert result["event_audit"].query("event_key == 'E0'").projection_status.iloc[0].startswith("absent_event")
    assert result["sample_coverage"].query("module == 'M01'").score_status.iloc[0] == "insufficient_measured_events"


def test_independence_identity_and_model_corruption_rejected(tmp_path):
    matrix, baseline, samples, model = fixture()
    with pytest.raises(ValueError, match="overlaps"):
        project_psi_modules(model, matrix, samples, feature_space_id="hg19|test")
    new = matrix.rename(columns=lambda x:"new_"+x)
    meta = samples.assign(sample_id=new.columns, biological_sample_id=new.columns)
    with pytest.raises(ValueError, match="feature_space"):
        project_psi_modules(model, new, meta, feature_space_id="hg38|other")
    meta.loc[1,"biological_sample_id"] = meta.loc[0,"biological_sample_id"]
    with pytest.raises(ValueError, match="technical"):
        project_psi_modules(model, new, meta, feature_space_id="hg19|test")
    path = tmp_path/"model.json"
    save_frozen_model(model, path)
    model["modules"][0]["centers"][0] += .01
    path.write_text(json.dumps(model))
    with pytest.raises(ValueError, match="modified"):
        load_frozen_model(path)


def test_public_entrypoints_registered_with_original_network_defaults():
    from katana_next.workflows.registry import WorkflowRegistry
    registry = WorkflowRegistry()
    values = registry.get("programs.psi_network").normalize_parameters({})
    assert (values["signed"], values["method"], values["soft_power"]) == (True,"bicor",6)
    assert registry.get("programs.psi_projection").definition.workspace_id == "programs"


def test_r_backend_rejects_silent_method_substitution():
    from katana_cosplice.r_wgcna import exact_wgcna_labels
    with pytest.raises(ValueError, match="unsupported"):
        exact_wgcna_labels(pd.DataFrame(), soft_power=6, min_module_size=20, merge_cut_height=.25,
                          deep_split=2, pam_respects_dendro=False, cor_type="spearman")


def test_cli_tsv_roundtrip_uses_public_projection_entrypoint(tmp_path, monkeypatch):
    from katana_next.psi_cli import main
    matrix, baseline, samples, model = fixture()
    new = matrix.rename(columns=lambda x:"new_"+x)
    meta = samples.assign(sample_id=new.columns, biological_sample_id=new.columns)
    model_path = tmp_path/"model.json"
    save_frozen_model(model, model_path)
    new.to_csv(tmp_path/"psi.tsv", sep="\t")
    meta.to_csv(tmp_path/"samples.tsv", sep="\t", index=False)
    monkeypatch.setattr("sys.argv", ["psi_cli", "project", "--psi", str(tmp_path/"psi.tsv"),
        "--samples", str(tmp_path/"samples.tsv"), "--model", str(model_path),
        "--out", str(tmp_path/"result"), "--feature-space-id", "hg19|test"])
    main()
    result = pd.read_csv(tmp_path/"result/scores.tsv", sep="\t", index_col=0)
    assert np.allclose(result, baseline.program_scores, atol=1e-12)
    assert json.loads((tmp_path/"result/COMPLETE.json").read_text())["module_refit"] is False
    assert len(pd.read_csv(tmp_path/"result/event_audit.tsv", sep="\t")) == 40


def test_release_launcher_does_not_reuse_an_original_version_server(monkeypatch):
    import io
    from katana_next import __version__
    from katana_next.launcher import service_ready
    health = {"status":"ok","version":"2.0.0b2"}
    def response(url, timeout):
        payload = health if url.endswith("health") else [{"workspace_id":k} for k in ("project","differential","programs","reports")]
        return io.BytesIO(json.dumps(payload).encode())
    monkeypatch.setattr("urllib.request.urlopen", response)
    assert not service_ready("http://127.0.0.1:8765/")
    health["version"] = __version__
    assert service_ready("http://127.0.0.1:18765/")
