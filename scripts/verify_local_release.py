"""Check the installed local application's own source/build ledger."""
import hashlib
import json
from pathlib import Path

root = Path(__file__).resolve().parents[1]
ledger = json.loads((root / 'release-integrity.json').read_text('utf-8'))
for name, expected in ledger['files'].items():
    path = (root / name).resolve()
    if not path.is_relative_to(root) or not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest() != expected:
        raise SystemExit('Release file changed or missing: ' + name)
print('Katana local release verified:', ledger['release'], len(ledger['files']), 'files')
