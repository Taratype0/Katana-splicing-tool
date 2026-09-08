from __future__ import annotations

from typing import Any

from pydantic import BaseModel, ConfigDict, Field


class ApiModel(BaseModel):
    model_config = ConfigDict(extra="forbid")


class OpenProjectRequest(ApiModel):
    manifest_path: str


class RunWorkflowRequest(ApiModel):
    parameters: dict[str, Any] = Field(default_factory=dict)
    force: bool = False


class CreateSelectionRequest(ApiModel):
    name: str
    entity_type: str
    source_artifact_id: str
    entity_ids: list[str]
    filters: dict[str, Any] = Field(default_factory=dict)
