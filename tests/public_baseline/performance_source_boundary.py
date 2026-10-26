"""Explicit execution-only AST revision, backed by numerical reference tests."""
import ast
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SNAPSHOT = ROOT / 'tests/fixtures/performance_reference_20261026'


def restore_execution_definitions(module, source):
    record = json.loads((SNAPSHOT / 'execution-definitions.json').read_text())
    approved = record['definitions'].get(module)
    if approved is None:
        return source
    old_path = SNAPSHOT / (module + '.py')
    sources = json.loads((SNAPSHOT / 'sources.json').read_text())['sources']
    assert hashlib.sha256(old_path.read_bytes()).hexdigest() == sources['katana_next/science/' + module + '.py']
    old = ast.parse(old_path.read_text('utf8'))
    before = {n.name:n for n in old.body if hasattr(n,'name')}
    current = ast.parse(source)
    found = set()
    body = []
    for node in current.body:
        name = getattr(node,'name',None)
        if name in approved:
            digest = hashlib.sha256(ast.dump(node,include_attributes=False).encode()).hexdigest()
            assert digest == approved[name], 'Scientific AST changed outside declared execution revision: ' + str(name)
            found.add(name)
            if name in before:
                body.append(before[name])
        else:
            body.append(node)
    assert found == set(approved), 'Scientific AST changed: execution definition missing'
    current.body = body
    return ast.unparse(current)
