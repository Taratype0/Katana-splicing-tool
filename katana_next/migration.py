from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import sys
from pathlib import Path

import yaml

from katana_next.domain.models import (
    AnnotationReference,
    Contrast,
    ContrastArm,
    InputAsset,
    ProjectManifest,
    Sample,
)
from katana_next.infrastructure.manifest_store import ManifestStore


EVENT_TYPES = ("SE", "RI", "A3SS", "A5SS", "MXE")


def migrate_legacy_project(
    legacy_config: str | Path,
    output_root: str | Path,
    destination: str | Path | None = None,
) -> Path:
    source = Path(legacy_config).expanduser().resolve()
    if not source.is_file():
        raise FileNotFoundError(source)
    raw = yaml.safe_load(source.read_text(encoding="utf-8")) or {}
    project_root = source.parent
    target_root = Path(output_root).expanduser().resolve()
    target_manifest = (
        Path(destination).expanduser().resolve()
        if destination
        else target_root.parent / "katana.project.json"
    )
    target_root.mkdir(parents=True, exist_ok=True)

    sample_rows = [row for row in raw.get("isoform_samples", []) if row.get("include", True)]
    groups: dict[str, list[str]] = {}
    samples: list[Sample] = []
    for row in sample_rows:
        sample_id = str(row.get("sample_id", "")).strip()
        condition = str(row.get("condition") or row.get("sample_group") or "").strip()
        if not sample_id or not condition:
            continue
        groups.setdefault(condition, []).append(sample_id)
        assets: dict[str, Path] = {}
        quant = _existing(row.get("quant_sf"))
        if quant:
            assets["quant_sf"] = quant
        bam = _find_first(
            project_root,
            (f"{sample_id}_sorted.bam", f"{sample_id}.bam"),
            excluded=("katana_output",),
        )
        if bam:
            assets["bam"] = bam
            bai = next(
                (
                    candidate
                    for candidate in (
                        Path(str(bam) + ".bai"),
                        bam.with_suffix(".bai"),
                    )
                    if candidate.is_file()
                ),
                None,
            )
            if bai:
                assets["bai"] = bai
        samples.append(
            Sample(
                sample_id=sample_id,
                biological_sample_id=sample_id,
                condition=condition,
                batch=str(row.get("batch") or "") or None,
                assets=assets,
            )
        )
    if not samples:
        raise ValueError("Legacy project has no usable isoform_samples/sample metadata.")

    contrasts: list[Contrast] = []
    assets: list[InputAsset] = []
    rmats_mode = str(raw.get("selected_rmats_mode", "JC")).upper()
    comparisons = [item for item in raw.get("comparisons", []) if item.get("comparison_id")]
    all_files = [
        path
        for path in project_root.rglob("*")
        if path.is_file() and "katana_output" not in {part.lower() for part in path.parts}
    ]
    for item in comparisons:
        contrast_id = str(item["comparison_id"])
        numerator = str(item.get("experiment_group") or item.get("source_experiment_group") or "").strip()
        denominator = str(item.get("control_group") or item.get("source_control_group") or "").strip()
        if numerator not in groups or denominator not in groups:
            continue
        label = str(item.get("display_name") or f"{numerator} vs {denominator}")
        contrasts.append(
            Contrast(
                contrast_id=contrast_id,
                label=label,
                numerator=ContrastArm(label=numerator, sample_ids=tuple(groups[numerator])),
                denominator=ContrastArm(label=denominator, sample_ids=tuple(groups[denominator])),
            )
        )
        rmats_name = str(item.get("rmats_name") or contrast_id)
        reverse_splicing = bool(item.get("reverse_splicing", item.get("reverse_direction", False)))
        for event_type in EVENT_TYPES:
            filename = f"{event_type}.MATS.{rmats_mode}.txt"
            match = _best_comparison_file(all_files, filename, rmats_name)
            if match:
                assets.append(
                    InputAsset(
                        asset_id=f"rmats::{contrast_id}::{event_type}",
                        kind="rmats_result",
                        path=match,
                        metadata={
                            "contrast_id": contrast_id,
                            "event_type": event_type,
                            "reverse_direction": reverse_splicing,
                            "numerator_sample_ids": groups[numerator],
                            "denominator_sample_ids": groups[denominator],
                            "source_comparison": rmats_name,
                            "mode": rmats_mode,
                        },
                    )
                )
        deg_name = str(item.get("deg_name") or contrast_id)
        deg = next(
            (
                path
                for path in all_files
                if path.name.lower() == f"{deg_name}.deseq2.results.tsv".lower()
            ),
            None,
        )
        if deg:
            assets.append(
                InputAsset(
                    asset_id=f"deg::{contrast_id}",
                    kind="deg_result",
                    path=deg,
                    metadata={
                        "contrast_id": contrast_id,
                        "reverse_direction": bool(item.get("reverse_deg", False)),
                    },
                )
            )
        dtu_name = str(item.get("dtu_name") or contrast_id)
        dtu = _best_dtu_file(all_files, dtu_name)
        if dtu:
            assets.append(
                InputAsset(
                    asset_id=f"dtu::{contrast_id}",
                    kind="dtu_result",
                    path=dtu,
                    metadata={"contrast_id": contrast_id},
                )
            )
    if not contrasts:
        raise ValueError("No legacy comparison could be mapped to declared sample groups.")
    if not any(asset.kind == "rmats_result" for asset in assets):
        raise ValueError("No rMATS event tables were discovered for the migrated contrasts.")

    cosplice = _existing(
        project_root
        / "katana_output"
        / "09_cosplice_dynamic_uncertainty_validation"
        / "cosplice_config.yaml"
    )
    if cosplice:
        assets.append(
            InputAsset(
                asset_id="cosplice::legacy-config",
                kind="cosplice_config",
                path=cosplice,
                required=False,
            )
        )
    tool_paths = raw.get("tool_paths") or {}
    gtf = _existing(tool_paths.get("gtf"))
    genome_fasta = _find_reference(project_root, "*primary_assembly.genome.fa")
    transcript_fasta = _find_reference(project_root, "*transcripts.fa")
    manifest = ProjectManifest(
        project_id=_slug(project_root.name),
        name=f"{project_root.name} / Katana Next",
        input_root=project_root,
        output_root=target_root,
        annotation=AnnotationReference(
            genome_build=str(raw.get("genome_build") or "GRCm39"),
            gtf=gtf,
            genome_fasta=genome_fasta,
            transcript_fasta=transcript_fasta,
        ),
        assets=tuple(assets),
        samples=tuple(samples),
        contrasts=tuple(contrasts),
        metadata={
            "migrated_from": str(source),
            "direction_contract": "all effects are numerator minus denominator; no runtime direction selector",
            "samtools_path": str((raw.get("tool_paths") or {}).get("samtools", "samtools")),
            "rscript_path": _discover_rscript(tool_paths.get("rscript")),
        },
    )
    saved = ManifestStore().save(manifest, target_manifest)
    audit = {
        "source": str(source),
        "manifest": str(saved),
        "output_root": str(target_root),
        "samples": len(samples),
        "contrasts": len(contrasts),
        "assets_by_kind": {
            kind: sum(asset.kind == kind for asset in assets)
            for kind in sorted({asset.kind for asset in assets})
        },
    }
    saved.with_name("migration_audit.json").write_text(
        json.dumps(audit, indent=2), encoding="utf-8"
    )
    return saved


def _existing(value: object) -> Path | None:
    if not value:
        return None
    path = Path(str(value)).expanduser().resolve()
    return path if path.is_file() else None


def _find_first(
    root: Path, names: tuple[str, ...], *, excluded: tuple[str, ...] = ()
) -> Path | None:
    lowered = {name.lower() for name in names}
    return next(
        (
            path
            for path in root.rglob("*")
            if path.is_file()
            and path.name.lower() in lowered
            and not any(part.lower() in excluded for part in path.parts)
        ),
        None,
    )


def _best_comparison_file(files: list[Path], filename: str, comparison: str) -> Path | None:
    normalized = _slug(comparison)
    candidates = [path for path in files if path.name.lower() == filename.lower()]
    return next(
        (path for path in candidates if normalized in _slug(str(path.parent))),
        None,
    )


def _best_dtu_file(files: list[Path], comparison: str) -> Path | None:
    normalized = _slug(comparison)
    candidates = [
        path
        for path in files
        if normalized in _slug(str(path.parent))
        and path.suffix.lower() == ".tsv"
        and (
            path.name.lower().startswith("getadjustedpvalues")
            or path.name.lower().startswith("dexseqresults")
        )
    ]
    candidates.sort(key=lambda path: (not path.name.lower().startswith("getadjustedpvalues"), len(str(path))))
    return candidates[0] if candidates else None


def _find_reference(project_root: Path, pattern: str) -> Path | None:
    for root in (project_root, project_root.parent, project_root.parent.parent):
        match = next((path for path in root.glob(pattern) if path.is_file()), None)
        if match:
            return match.resolve()
    return None


def _discover_rscript(configured: object) -> str:
    configured_path = _existing(configured)
    if configured_path:
        return str(configured_path)
    candidates: list[Path] = []
    on_path = shutil.which("Rscript")
    if on_path:
        candidates.append(Path(on_path))
    prefixes = [Path(sys.prefix)]
    conda_prefix = os.environ.get("CONDA_PREFIX")
    if conda_prefix:
        prefixes.append(Path(conda_prefix))
    conda_executable = os.environ.get("CONDA_EXE")
    if conda_executable:
        prefixes.append(Path(conda_executable).resolve().parent.parent)
    prefixes.extend(
        [
            Path("D:/anaconda3"),
            Path("C:/ProgramData/anaconda3"),
            Path.home() / "anaconda3",
            Path.home() / "miniconda3",
        ]
    )
    for prefix in dict.fromkeys(path.resolve() for path in prefixes):
        candidates.extend(
            [prefix / "Scripts" / "Rscript.exe", prefix / "Lib" / "R" / "bin" / "Rscript.exe"]
        )
        env_root = prefix / "envs"
        if env_root.is_dir():
            candidates.extend(sorted(env_root.glob("*/Scripts/Rscript.exe")))
            candidates.extend(sorted(env_root.glob("*/Lib/R/bin/Rscript.exe")))
    return str(next((path.resolve() for path in candidates if path.is_file()), ""))


def _slug(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "-", value.lower()).strip("-")


def main() -> None:
    parser = argparse.ArgumentParser(description="Migrate a Katana 1.x YAML project to Katana Next.")
    parser.add_argument("legacy_config")
    parser.add_argument("--output-root", required=True)
    parser.add_argument("--manifest")
    arguments = parser.parse_args()
    print(migrate_legacy_project(arguments.legacy_config, arguments.output_root, arguments.manifest))


if __name__ == "__main__":
    main()
