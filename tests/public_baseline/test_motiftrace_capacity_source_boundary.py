"""Prove the three source edits remove declared capacity guards only.

The archived source bytes keep their original accepted digests. An explicit,
exhaustive AST edit list describes the allowed execution changes; everything
else, including scores, matching, windows and P/FDR/FWER, must remain identical.
"""
import ast
import copy
import hashlib
import json
from pathlib import Path
import textwrap

import pytest


ROOT = Path(__file__).resolve().parents[2]
SNAPSHOT = ROOT / "tests/fixtures/frozen_motiftrace_capacity_20261017"
BASELINE_RECORD_SHA256 = "d907f6dbb06b323a9d39092be373b40409124de837c09d3763d61de7b9eac7fa"
CAPACITY_SOURCES = {
    "motiftrace_events": "4b85600fa8741133150e50ab35d7976035311ce9d47a9552c9bb391957c90f7e",
    "motiftrace_custom": "56bf46ba62ea90f77d5d0d75cc71756ff81e149c9b768191213425c79fd8dd8c",
    "motiftrace_clusters": "c4f2937ad5404e6a35990ee7e351bab89ebd36f839f49e001966d00b2d5001e2",
}


def _dump(node):
    return ast.dump(node, include_attributes=False)


def _statement(text):
    tree = ast.parse(textwrap.dedent(text).strip())
    assert len(tree.body) == 1
    return tree.body[0]


class _DocumentaryText(ast.NodeTransformer):
    """Only module/function docstrings are outside the numerical AST contract."""
    def _visit_body(self, node):
        self.generic_visit(node)
        if (node.body and isinstance(node.body[0], ast.Expr) and
                isinstance(node.body[0].value, ast.Constant) and
                isinstance(node.body[0].value.value, str)):
            node.body = node.body[1:]
        return node

    visit_Module = _visit_body
    visit_FunctionDef = _visit_body
    visit_AsyncFunctionDef = _visit_body


def _replace_statement(tree, function, old, new=None):
    """Replace exactly one complete declared statement in its declared scope."""
    wanted = _dump(_statement(old))
    replacement = _statement(new) if new is not None else None
    matches = []

    class Edit(ast.NodeTransformer):
        def __init__(self):
            self.scope = None

        def visit_FunctionDef(self, node):
            previous, self.scope = self.scope, node.name
            result = self.generic_visit(node)
            self.scope = previous
            return result

        def generic_visit(self, node):
            if self.scope == function and _dump(node) == wanted:
                matches.append(node)
                return copy.deepcopy(replacement)
            return super().generic_visit(node)

    tree = Edit().visit(tree)
    assert len(matches) == 1, (function, old, len(matches))
    return tree


def _compatibility_default(tree, function, name, old_default):
    functions = [node for node in tree.body if isinstance(node, ast.FunctionDef) and node.name == function]
    assert len(functions) == 1
    args = functions[0].args
    indices = [i for i, arg in enumerate(args.kwonlyargs) if arg.arg == name]
    assert len(indices) == 1
    index = indices[0]
    assert _dump(args.kw_defaults[index]) == _dump(ast.Constant(old_default))
    args.kw_defaults[index] = ast.Constant(None)


def _expected_capacity_only_tree(module, source):
    tree = ast.parse(source)
    if module == "motiftrace_events":
        _compatibility_default(tree, "extract_event_rna", "max_total_bases", 20_000_000)
        tree = _replace_statement(tree, "extract_event_rna",
            'fasta, rows, total = IndexedFasta(Path(fasta_path)), [], 0',
            'fasta, rows = IndexedFasta(Path(fasta_path)), []')
        tree = _replace_statement(tree, "extract_event_rna", 'total += sum(b-a for a, b in blocks)')
        tree = _replace_statement(tree, "extract_event_rna", '''
            if total > max_total_bases:
                raise ValueError("Explicit local RNA base budget exceeded; nothing is silently truncated")
        ''')
        tree = _replace_statement(tree, "validate_models", '''
            if not models or len(models) > 256:
                raise ValueError("Require 1-256 explicitly selected models; no silent library truncation")
        ''', '''
            if not models:
                raise ValueError("Require at least one explicitly selected model; no silent library truncation")
        ''')
        tree = _replace_statement(tree, "trace_instances", '''
            if rna.loc[rna.layer.eq("mature_local_form"), "length"].sum()*len(models) > 25_000_000:
                raise ValueError("Explicit motif-by-local-RNA budget exceeded; no silent model/event truncation")
        ''')
    elif module == "motiftrace_custom":
        tree = _replace_statement(tree, None, 'MAX_MODELS = 256', 'MAX_MODELS = None')
        tree = _replace_statement(tree, "_file", '''
            if raw.stat().st_size > 16 * 1024 ** 2:
                raise ValueError("Custom library exceeds the explicit 16MiB local input budget")
        ''')
        tree = _replace_statement(tree, "load_custom_library", '''
            if not set(FIELDS).issubset(table) or table.empty or len(table) > 1024:
                raise ValueError("Require 1-1024 source rows with motif_id/RBP/species/evidence/source")
        ''', '''
            if not set(FIELDS).issubset(table) or table.empty:
                raise ValueError("Require source rows with motif_id/RBP/species/evidence/source")
        ''')
        tree = _replace_statement(tree, "load_custom_library", '''
            if len(models) > MAX_MODELS:
                raise ValueError(f"Custom job exceeds the explicit {MAX_MODELS}-model budget; never silently truncate a library")
        ''')
        tree = _replace_statement(tree, "add_custom_results", '''
            if len(library["models"]) * (14 + 2 * sum(span for span, _ in windows.values())) > 1_000_000:
                raise ValueError("Custom analysis exceeds one million planned tests per contrast; preserve inputs and use an explicitly budgeted job")
        ''')
    elif module == "motiftrace_clusters":
        _compatibility_default(tree, "build_paired_datasets", "max_tests", 100_000)
        _compatibility_default(tree, "build_paired_datasets", "max_matrix_cells", 8_000_000)
        tree = _replace_statement(tree, "paired_cluster_family", '''
            if not 99 <= permutations <= 9999 or minimum_pairs < 2:
                raise ValueError("Require 99-9999 fixed permutations and at least two independent pairs")
        ''', '''
            if (isinstance(permutations, bool) or not isinstance(permutations, (int, np.integer))
                    or permutations < 99 or minimum_pairs < 2):
                raise ValueError("Require at least 99 integer permutations and at least two independent pairs")
        ''')
        tree = _replace_statement(tree, "paired_cluster_family", '''
            if permutations*len(pair_ids) > 8_000_000:
                raise ValueError("Explicit whole-pair permutation memory budget exceeded")
        ''')
        tree = _replace_statement(tree, "build_paired_datasets", '''
            if tests > max_tests:
                raise ValueError("Explicit new-test budget exceeded; no models/windows are silently dropped")
        ''')
        tree = _replace_statement(tree, "build_paired_datasets", '''
            if cells > max_matrix_cells:
                raise ValueError("Explicit paired-matrix memory budget exceeded; no gene pairs silently dropped")
        ''')
    else:
        raise AssertionError("Unapproved changed source: " + module)
    return _DocumentaryText().visit(tree)


def historical_integration_record():
    path = SNAPSHOT / "integration.before-capacity-removal.json"
    raw = path.read_bytes()
    assert hashlib.sha256(raw).hexdigest() == BASELINE_RECORD_SHA256
    return json.loads(raw)


def assert_capacity_only_change(module, current_source=None):
    assert module in CAPACITY_SOURCES
    snapshot = SNAPSHOT / (module + ".py")
    raw = snapshot.read_bytes()
    assert hashlib.sha256(raw).hexdigest() == CAPACITY_SOURCES[module]
    expected = _expected_capacity_only_tree(module, raw.decode("utf-8"))
    if current_source is None:
        current_source = (ROOT / "katana_next/science" / (module + ".py")).read_text("utf-8")
    from performance_source_boundary import restore_execution_definitions
    actual = _DocumentaryText().visit(ast.parse(restore_execution_definitions(module, current_source)))
    assert _dump(actual) == _dump(expected), "Scientific AST changed outside declared capacity-only boundary: " + module


@pytest.mark.parametrize("module", CAPACITY_SOURCES)
def test_exact_frozen_source_allows_only_declared_capacity_edits(module):
    assert_capacity_only_change(module)


def test_contract_preserves_old_source_records_and_binds_three_archived_snapshots():
    record = json.loads((ROOT / "docs/motiftrace/integration.v1.json").read_text("utf-8"))
    historical = historical_integration_record()
    assert set(record["scientific_sources"]) == set(historical["scientific_sources"])
    assert record["integration_version"] == "1.0.0.post7"
    assert not record["scientific_calculation_changes"]
    assert not record["scientific_parameter_changes"]
    assert not record["completed_studies_recomputed"]
    for name, old in historical["scientific_sources"].items():
        item = record["scientific_sources"][name]
        assert item["frozen_sha256"] == old["frozen_sha256"]
        module = Path(name).stem
        if module in CAPACITY_SOURCES:
            assert old["current_sha256"] == CAPACITY_SOURCES[module]
            assert item["change_kind"] == ("execution_capacity_guards_removed" if module == "motiftrace_custom" else "execution_scoring_reuse_and_streaming")
            assert item["same"] is False
            assert item["frozen_snapshot"] == (SNAPSHOT / (module + ".py")).relative_to(ROOT).as_posix()
            assert_capacity_only_change(module)
        elif name == "katana_next/science/motiftrace_consequences.py":
            from performance_source_boundary import restore_execution_definitions
            actual = restore_execution_definitions(module, (ROOT / name).read_text("utf8"))
            reference = (ROOT / "tests/fixtures/performance_reference_20261026" / (module + ".py")).read_text("utf8")
            assert hashlib.sha256(reference.encode()).hexdigest() == old["current_sha256"] or hashlib.sha256((ROOT / "tests/fixtures/performance_reference_20261026" / (module + ".py")).read_bytes()).hexdigest() == old["current_sha256"]
            assert _dump(ast.parse(actual)) == _dump(ast.parse(reference))
        else:
            # Do not legitimize additional numerical or presentation changes by
            # merely updating their current hashes in the new register.
            assert item["current_sha256"] == old["current_sha256"]
            assert hashlib.sha256((ROOT / name).read_bytes()).hexdigest() == old["current_sha256"]


@pytest.mark.parametrize("module,old,new", [
    ("motiftrace_events", "candidate[starts] = rel >= .8", "candidate[starts] = rel >= .7"),
    ("motiftrace_events", "background=(.25, .25, .25, .25)", "background=(.4, .2, .2, .2)"),
    ("motiftrace_custom", '"minimum_per_group": 20', '"minimum_per_group": 19'),
    ("motiftrace_clusters", "CLUSTER_Z = float(norm.isf(.01))", "CLUSTER_Z = float(norm.isf(.02))"),
    ("motiftrace_clusters", "q[valid] = bh_adjust(p[valid])", "q[valid] = p[valid]"),
    ("motiftrace_clusters", "minimum_pairs=20", "minimum_pairs=19"),
])
def test_boundary_rejects_scoring_or_statistical_changes(module, old, new):
    current = (ROOT / "katana_next/science" / (module + ".py")).read_text("utf-8")
    assert current.count(old) == 1
    with pytest.raises(AssertionError, match="Scientific AST changed"):
        assert_capacity_only_change(module, current.replace(old, new, 1))


@pytest.mark.parametrize("module", CAPACITY_SOURCES)
def test_boundary_rejects_unrelated_new_code(module):
    current = (ROOT / "katana_next/science" / (module + ".py")).read_text("utf-8")
    with pytest.raises(AssertionError, match="Scientific AST changed"):
        assert_capacity_only_change(module, current + "\nUNRELATED_SCIENTIFIC_CHANGE = True\n")
