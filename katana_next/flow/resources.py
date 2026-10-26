"""Execution-only resource overrides that can change during managed resume."""
from __future__ import annotations
import re

FIELDS = {"cpus", "memory", "time", "queue", "account", "maxRetries"}


def _directives(values):
    if not isinstance(values, dict) or set(values) - FIELDS:
        raise ValueError("Resource fields: " + ", ".join(sorted(FIELDS)))
    lines = []
    for key, value in values.items():
        if key in {"cpus", "maxRetries"}:
            if isinstance(value, bool) or not isinstance(value, int) or value < (1 if key == "cpus" else 0):
                raise ValueError(f"{key} must be an integer >= {1 if key == 'cpus' else 0}")
            lines.append(f"{key} = {value}")
        elif key in {"memory", "time"}:
            units = r"(?:KB|MB|GB|TB)" if key == "memory" else r"(?:ms|s|m|h|d)"
            if not isinstance(value, str) or not re.fullmatch(r"[1-9][0-9]*(?:\.[0-9]+)?\s*" + units, value):
                raise ValueError(f"{key} requires a positive amount and unit")
            # Resource retry escalation remains explicit and bounded by the
            # site's scheduler; no queue is inferred from requested memory.
            lines.append(f"{key} = '{value}'")
        else:
            if not isinstance(value, str) or not re.fullmatch(r"[A-Za-z0-9_.-]+", value):
                raise ValueError(f"Invalid scheduler {key}")
            lines.append(f"clusterOptions = '--account={value}'" if key == "account" else f"queue = '{value}'")
    return lines


def render_resources(value):
    if not isinstance(value, dict) or set(value) - {"defaults", "labels", "processes"}:
        raise ValueError("Resource configuration accepts defaults, labels and processes only; scientific parameters belong in flow_config")
    lines = ["// Generated execution-only settings. No scientific parameter overrides.", "process {"]
    if not isinstance(value.get("defaults", {}), dict) or set(value.get("defaults", {})) - {"queue", "account"}:
        raise ValueError("defaults accepts queue/account; set CPU, memory and time under exact process names")
    lines += ["    " + line for line in _directives(value.get("defaults", {}))]
    labels = value.get("labels", {})
    if not isinstance(labels, dict):
        raise ValueError("labels must map Nextflow process labels to resource settings")
    for name, options in labels.items():
        if not isinstance(name, str) or not re.fullmatch(r"[a-z][a-z0-9_]*", name):
            raise ValueError("Use an exact lowercase Nextflow process label")
        directives = _directives(options)
        if not directives:
            continue
        lines.append(f"    withLabel: '{name}' {{")
        lines += ["        " + line for line in directives]
        lines.append("    }")
    processes = value.get("processes", {})
    if not isinstance(processes, dict):
        raise ValueError("processes must map process names to resource settings")
    for name, options in processes.items():
        if not isinstance(name, str) or not re.fullmatch(r"[A-Z][A-Z0-9_]*", name):
            raise ValueError("Use the exact uppercase Nextflow process alias")
        directives = _directives(options)
        if not directives:
            continue
        lines.append(f"    withName: '(^|.*:){name}' {{")
        lines += ["        " + line for line in directives]
        lines.append("    }")
    if len(lines) == 2:
        return lines[0] + "\n"
    return "\n".join([*lines, "}", ""])
