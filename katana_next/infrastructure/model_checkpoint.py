"""Atomic model-prefix checkpoints; partial fits are never published as results."""
import hashlib
import json
from pathlib import Path


class ModelCheckpoint:
    def __init__(self, path: Path, fingerprint: str):
        self.path, self.fingerprint = path, fingerprint

    def load(self):
        if not self.path.exists():
            return 0, []
        record = json.loads(self.path.read_text(encoding="utf-8"))
        payload = record["payload"]
        digest = hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest()
        if record.get("sha256") != digest or payload.get("fingerprint") != self.fingerprint:
            raise ValueError("Model checkpoint fingerprint or integrity mismatch; use a new run directory")
        return int(payload["completed_groups"]), payload["rows"]

    def save(self, completed, rows):
        payload = {"fingerprint": self.fingerprint, "completed_groups": completed, "rows": rows}
        record = {"payload": payload, "sha256": hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest()}
        self.path.parent.mkdir(parents=True, exist_ok=True)
        temporary = self.path.with_suffix(".partial")
        temporary.write_text(json.dumps(record), encoding="utf-8")
        temporary.replace(self.path)


def model_groups(grouped, rows, checkpoint=None, report=None):
    completed, prior = checkpoint.load() if checkpoint else (0, [])
    rows.extend(prior)
    total = len(grouped)
    if completed > total:
        raise ValueError("Model checkpoint contains more groups than the current input")
    for index, value in enumerate(grouped):
        if index < completed:
            continue
        if checkpoint and index % 100 == 0:
            checkpoint.save(index, rows)
        if report and index % 10 == 0:
            report(index, total)
        yield value
    if checkpoint:
        checkpoint.save(total, rows)
    if report:
        report(total, total)
