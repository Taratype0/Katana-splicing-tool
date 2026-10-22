"""Public, data-free guard for the accepted MotifTrace scientific release.

Digests originate from the accepted 1.0.0-frozen.20261017 source archive.
Changing one requires a separately versioned scientific release, rather than
silently changing calculations while adding packaging or presentation adapters.
No local study files, paths, sequences, or completed-result bundles are needed.
"""
import hashlib
import importlib
import json
from pathlib import Path
import re

import pandas as pd
import pytest


FROZEN_NUMERICAL_SHA256 = {
    "motiftrace_events": "4b85600fa8741133150e50ab35d7976035311ce9d47a9552c9bb391957c90f7e",
    "motiftrace_clusters": "c4f2937ad5404e6a35990ee7e351bab89ebd36f839f49e001966d00b2d5001e2",
    "motiftrace_consequences": "7cf122e48e87caf482da542758e625ce69bc23be3d2f16d0af1e03e25f82ccf1",
    "motiftrace_usage": "e26db6e87d6ce6c5ca14d405aa249c8dbaecf96d7e064f97e73a1dc58f3c39e0",
    "motiftrace_parallel": "73b4aa75ce81c128ebbb1a4dd4be0e04aed37e48ab1e7c8352e6ae19c54320a8",
    "motiftrace_windows": "fec10019bef1e9ed396c6568d77af52f4dedfeb3f3b7b259ba2e1174473ddeec",
    "motiftrace_custom": "56bf46ba62ea90f77d5d0d75cc71756ff81e149c9b768191213425c79fd8dd8c",
    "rbp_sequence_results": "565b32b6e4b84323f8c95df8853d61d32bcaabf3102e34cb4402e0407a52e19e",
    "rbp_short_extension": "e3aef52b9c81aa1593e034d43476dde9bbad97dc35335d30fe6e6d70ea2c3b6a",
}


def test_integration_record_binds_existing_versions_parameters_and_sources():
    from katana_motiftrace import mode_catalog
    from katana_next.workflows.registry import WorkflowRegistry
    root = Path(__file__).resolve().parents[2]
    record = json.loads((root / 'docs/motiftrace/integration.v1.json').read_text())
    assert record['scientific_release'] == mode_catalog()['scientific_release']
    assert not record['scientific_calculation_changes']
    assert not record['scientific_parameter_changes']
    registry = WorkflowRegistry()
    for mode in record['modes']:
        workflow = registry.get(mode['workflow_id'])
        assert workflow.version == mode['workflow_version']
        params = [p.to_dict() for p in workflow.definition.parameters]
        digest = hashlib.sha256(json.dumps(params, sort_keys=True, separators=(',', ':')).encode()).hexdigest()
        assert digest == mode['parameter_definitions_sha256']
    for name, digest in record['interface_source_sha256'].items():
        assert hashlib.sha256((root / name).read_bytes()).hexdigest() == digest
    for name, values in record['scientific_sources'].items():
        assert hashlib.sha256((root / name).read_bytes()).hexdigest() == values['current_sha256']


@pytest.mark.parametrize("module,expected", FROZEN_NUMERICAL_SHA256.items())
def test_packaging_preserves_frozen_numerical_source(module, expected):
    loaded = importlib.import_module("katana_next.science." + module)
    assert hashlib.sha256(Path(loaded.__file__).read_bytes()).hexdigest() == expected


def test_parallel_records_preserve_native_extreme_P_NA_and_direction():
    from katana_next.science.motiftrace_parallel import records

    original = pd.DataFrame({
        "native_P": [5e-310, float("nan"), 0.05],
        "native_FDR": [1e-250, float("nan"), 0.20],
        "dpsi": [-0.125, 0.0, 0.25],
        "start0": [0, 10, 20],
        "end0": [5, 15, 25],
    })
    values = json.loads(json.dumps(records(original), allow_nan=False))
    assert values[0]["native_P"] == 5e-310 > 0
    assert values[0]["native_FDR"] == 1e-250
    assert values[1]["native_P"] is None
    assert values[1]["native_FDR"] is None
    assert [row["dpsi"] for row in values] == [-0.125, 0.0, 0.25]
    assert [(row["start0"], row["end0"]) for row in values] == [(0, 5), (10, 15), (20, 25)]


def test_parallel_presentation_keeps_distinct_evidence_units(tmp_path):
    from katana_next.science.motiftrace_parallel import render_parallel

    payload = {
        "family": {"models": 1, "windows": 1, "units_per_case": 1},
        "new_inference": False,
        "population": {"MotifTrace_conditional_family_P": 1e-200,
                       "MAPP_empirical_P_if_estimable": None},
        "loci": [{"event_id": "synthetic", "RNA_start0": 0, "RNA_end0": 5,
                  "blocks": [[100, 103], [200, 202]], "strand": "-",
                  "evidence": "descriptive_reference_sequence", "site_P": None}],
        "limitations": ["Population statistics do not become site statistics."],
    }
    report = render_parallel(payload, tmp_path / "parallel.html")
    text = Path(report).read_text(encoding="utf-8")
    embedded = json.loads(re.search(r'<script[^>]*type="application/json"[^>]*>(.*?)</script>', text, re.S).group(1))
    assert embedded["population"] == payload["population"]
    assert embedded["loci"]["$rows"] == [[payload["loci"][0][k] for k in embedded["loci"]["$columns"]]]
    assert embedded["new_inference"] is False
    assert embedded["limitations"] == payload["limitations"]
