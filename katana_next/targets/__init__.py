"""Portable selected-gene/event requests and independent evidence reports.

Importing this module does not open a project, launch a service or start jobs.
The same request can be executed by a local CLI or a server Python service.
"""
from __future__ import annotations


def export_targets(*args, **kwargs):
    from .bundle import export_targets as export
    return export(*args, **kwargs)


def plan_targets(*args, **kwargs):
    from .runner import plan_targets as plan
    return plan(*args, **kwargs)


def run_targets(*args, **kwargs):
    from .runner import run_targets as run
    return run(*args, **kwargs)


def inspect_results(*args, **kwargs):
    from .runner import inspect_results as inspect
    return inspect(*args, **kwargs)


__all__ = ["export_targets", "plan_targets", "run_targets", "inspect_results"]
