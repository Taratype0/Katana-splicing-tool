"""Project-bound Python API using Katana's existing execution services."""
from __future__ import annotations

from os import PathLike
from typing import TYPE_CHECKING

from .catalog import PACKAGE_VERSION, SCIENTIFIC_RELEASE, mode_catalog, resolve_mode, resolve_parameters

if TYPE_CHECKING:
    from katana_next.client import KatanaProject


class MotifTrace:
    """Open a project manifest or borrow an already-open ``KatanaProject``.

    Nothing runs on construction, inventory, planning or result listing. Analysis
    is an explicit ``run`` call. Closing a borrowed client leaves its owner open.
    """
    def __init__(self, project: str | PathLike | KatanaProject):
        self._owns_project = isinstance(project, (str, PathLike))
        if self._owns_project:
            from katana_next.client import KatanaProject
            self.project = KatanaProject(project)
        else:
            if not all(hasattr(project, name) for name in ("plan", "run", "artifacts")):
                raise TypeError("Provide a project manifest path or an open KatanaProject")
            self.project = project
        self._closed = False

    def _require_open(self):
        if self._closed:
            raise RuntimeError("MotifTrace client is closed")

    def inventory(self) -> dict:
        self._require_open()
        payload = mode_catalog()
        definitions = {entry["workflow_id"]: entry for entry in self.project.inventory()}
        for item in payload["modes"]:
            item["workflow"] = definitions[item["workflow_id"]]
        return payload

    @staticmethod
    def _parameters(mode, parameters, selection_id, contrast_ids):
        values = dict(parameters or {}) if isinstance(parameters, dict) or parameters is None else parameters
        spec, values = resolve_parameters(mode, values)
        for key, explicit in (("selection_id", selection_id), ("contrast_ids", contrast_ids)):
            if explicit is not None:
                if key in values and values[key] != explicit:
                    raise ValueError(f"Conflicting {key} in parameters and keyword argument")
                values[key] = explicit
        return spec, values

    def plan(self, mode: str, parameters: dict | None = None, *, selection_id=None, contrast_ids=None) -> dict:
        self._require_open()
        spec, values = self._parameters(mode, parameters, selection_id, contrast_ids)
        plan = dict(self.project.plan(spec["workflow_id"], values))
        plan["motiftrace"] = {"mode_id": mode, "package_version": PACKAGE_VERSION,
                             "scientific_release": SCIENTIFIC_RELEASE, "operation": spec["operation"]}
        return plan

    def run(self, mode: str, parameters: dict | None = None, *, selection_id=None, contrast_ids=None,
            force=False, timeout=None):
        self._require_open()
        spec, values = self._parameters(mode, parameters, selection_id, contrast_ids)
        # The shared service performs source/selection/input validation and exact
        # cache matching. No private runner or scientific fallback lives here.
        return self.project.run(spec["workflow_id"], values, force=force, timeout=timeout)

    def artifacts(self, mode: str | None = None, *, run_fingerprint: str | None = None):
        self._require_open()
        spec = resolve_mode(mode) if mode else None
        ids = [spec["workflow_id"]] if spec else list(dict.fromkeys(
            item["workflow_id"] for item in mode_catalog()["modes"]))
        return [item for workflow_id in ids for item in self.project.artifacts(workflow_id)
                if (run_fingerprint is None or item.run_fingerprint == run_fingerprint)
                and (not spec or all(item.parameters.get(key) == value
                                     for key, value in spec["fixed_parameters"].items()))]

    def select(self, source_artifact_id: str, row_ids, *, name="MotifTrace events"):
        """Save an exact event-row selection through the shared project service."""
        self._require_open()
        from katana_next.domain.models import SelectionSet
        if isinstance(row_ids, (str, bytes)):
            raise ValueError("row_ids must be a sequence of exact source row IDs")
        rows = tuple(str(value) for value in row_ids)
        if not rows or len(set(rows)) != len(rows):
            raise ValueError("Select nonempty, unique source row IDs")
        source = self.project.context.catalog.get_artifact(source_artifact_id)
        if source is None or source.project_id != self.project.project_id:
            raise ValueError("The source artifact must belong to this project")
        selection = SelectionSet(project_id=self.project.project_id, name=name,
            source_artifact_id=source_artifact_id, entity_type=source.entity_type, entity_ids=rows)
        return self.project.projects.save_selection(selection)

    def close(self):
        if not self._closed:
            if self._owns_project:
                self.project.close()
            self._closed = True

    def __enter__(self):
        self._require_open()
        return self

    def __exit__(self, *_):
        self.close()
