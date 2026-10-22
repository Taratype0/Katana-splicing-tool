import json
from pathlib import Path

from katana_next.workflows.registry import WorkflowRegistry

APP = Path(__file__).resolve().parents[2]
RELEASE = APP.parent
BASE = RELEASE.parent / "cosplice-usable-20261018"


def read(path):
    return json.loads(path.read_text(encoding="utf-8"))


def test_registry_retains_all_frozen_routes_and_only_adds_parallel_entry():
    original = read(BASE / "DEFAULTS_LOCK.json")["workflows"]
    current = {w.definition.workflow_id: w for w in WorkflowRegistry().all()}
    assert set(current) == set(original) | {"prioritize.motiftrace_parallel"}
    assert len(current) == 54
    for key, old in original.items():
        defaults = {p.key: p.default for p in current[key].definition.parameters}
        for name, value in old["defaults"].items():
            assert defaults[name] == value, (key, name)
    assert current["programs.hub_events"].version == "1.0.1"
    assert current["programs.djexpress"].version == original["programs.djexpress"]["version"]


def test_network_and_hub_implementations_and_export_controls_are_retained():
    import hashlib
    protected = [
        "katana_cosplice/config.py", "katana_cosplice/network.py", "katana_cosplice/consensus.py",
        "katana_next/science/hub_events.py", "katana_next/workflows/hub_events.py",
        "katana_next/workflows/djexpress.py", "katana_next/workflows/psi_network.py",
        "apps/web/src/components/ModuleDashboard.tsx", "apps/web/src/components/moduleDownloads.ts",
    ]
    for name in protected:
        assert hashlib.sha256((APP / name).read_bytes()).digest() == hashlib.sha256((BASE / "app" / name).read_bytes()).digest(), name
    old = (BASE / "app/katana_cosplice/pipeline.py").read_text(encoding="utf-8")
    expected = old.replace(" cosplice report</title>", " Katana-Splice-Network report</title>").replace(
        "<body><h1>cosplice:", "<body><h1>Katana-Splice-Network:")
    assert (APP / "katana_cosplice/pipeline.py").read_text(encoding="utf-8") == expected


def test_official_motiftrace_science_is_exact_and_original_registry_is_only_merged():
    import hashlib
    original = read(RELEASE / "evidence/motiftrace_original/motiftrace.official.json")
    for name, identity in original["runtime_identity"].items():
        if name == "katana_next/workflows/registry.py":
            continue
        assert hashlib.sha256((APP / name).read_bytes()).hexdigest() == identity["sha256"], name


def test_parallel_and_three_layer_routes_exposed_without_new_default_analysis():
    registry = WorkflowRegistry()
    main = registry.get("prioritize.rbp_mechanism")
    assert main.normalize_parameters({})["sequence_scope"] == "boundary_association"
    assert main.parent_artifacts(None, {"sequence_scope": "completed_parallel_models"}) == []
    for key in ("prioritize.rbp_support", "prioritize.motiftrace_consequences"):
        assert registry.get(key).definition.to_dict()["presentation_role"] == "primary"
    assert registry.get("prioritize.motiftrace_parallel").definition.required_asset_kinds == ("motiftrace_parallel_results",)
