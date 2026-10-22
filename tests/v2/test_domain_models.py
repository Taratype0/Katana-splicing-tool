from __future__ import annotations

from pathlib import Path

import pytest
from pydantic import ValidationError

from katana_next.domain.models import (
    AnnotationReference,
    Contrast,
    ContrastArm,
    ProjectManifest,
    Sample,
)
from katana_next.infrastructure.fingerprints import (
    fingerprint_file,
    workflow_fingerprint,
)


def test_contrast_rejects_overlapping_arms(tmp_path: Path) -> None:
    with pytest.raises(ValidationError, match="overlap"):
        Contrast(
            contrast_id="mut-vs-wt",
            label="MUT vs WT",
            numerator=ContrastArm(label="MUT", sample_ids=("s1",)),
            denominator=ContrastArm(label="WT", sample_ids=("s1",)),
        )


def test_manifest_rejects_unknown_contrast_sample(tmp_path: Path) -> None:
    with pytest.raises(ValidationError, match="unknown samples"):
        ProjectManifest(
            name="demo",
            input_root=tmp_path,
            output_root=tmp_path / "out",
            annotation=AnnotationReference(genome_build="GRCm39"),
            samples=(
                Sample(
                    sample_id="wt1",
                    biological_sample_id="mouse-1",
                    condition="WT",
                ),
            ),
            contrasts=(
                Contrast(
                    contrast_id="mut-vs-wt",
                    label="MUT vs WT",
                    numerator=ContrastArm(label="MUT", sample_ids=("mut1",)),
                    denominator=ContrastArm(label="WT", sample_ids=("wt1",)),
                ),
            ),
        )


def test_workflow_fingerprint_changes_with_input_and_parameters(tmp_path: Path) -> None:
    source = tmp_path / "input.tsv"
    source.write_text("event\nA\n", encoding="utf-8")
    first_input = fingerprint_file(source, force_full=True)
    first = workflow_fingerprint(
        workflow_id="demo",
        workflow_version="1",
        project_fingerprint="project",
        inputs=[first_input],
        parameters={"threshold": 0.1},
    )
    changed_parameter = workflow_fingerprint(
        workflow_id="demo",
        workflow_version="1",
        project_fingerprint="project",
        inputs=[first_input],
        parameters={"threshold": 0.2},
    )
    source.write_text("event\nB\n", encoding="utf-8")
    changed_input = workflow_fingerprint(
        workflow_id="demo",
        workflow_version="1",
        project_fingerprint="project",
        inputs=[fingerprint_file(source, force_full=True)],
        parameters={"threshold": 0.1},
    )

    assert len(first) == 64
    assert first != changed_parameter
    assert first != changed_input
