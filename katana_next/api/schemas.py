from __future__ import annotations

from typing import Any

from pydantic import BaseModel, ConfigDict, Field


class ApiModel(BaseModel):
    model_config = ConfigDict(extra="forbid")


class OpenProjectRequest(ApiModel):
    manifest_path: str


class ConfigureGeneFilterRequest(ApiModel):
    # Defaults apply only to a newly requested isolated scope. Opening an existing
    # manifest never injects these settings or changes its statistical history.
    exclude_mitochondrial: bool = True
    exclude_pseudogenes: bool = True
    exclude_rrna: bool = False
    protein_coding_only: bool = False
    excluded_genes: list[str] = Field(default_factory=list)
    min_expression: float = Field(default=0, ge=0, allow_inf_nan=False)
    min_samples: int = Field(default=1, ge=1)
    annotation_path: str = ""
    expression_matrix_path: str = ""


class MigrateProjectRequest(ApiModel):
    legacy_config_path: str
    output_root: str
    manifest_path: str | None = None


class ConfigureResourcesRequest(ApiModel):
    manifest_path: str
    resource_bundle_path: str
    output_root: str


class RunWorkflowRequest(ApiModel):
    parameters: dict[str, Any] = Field(default_factory=dict)
    force: bool = False


class CreateSelectionRequest(ApiModel):
    name: str
    entity_type: str
    source_artifact_id: str
    entity_ids: list[str]
    filters: dict[str, Any] = Field(default_factory=dict)
