from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from katana_cosplice.config import CoSpliceConfig
from katana_next.domain.models import (
    AnnotationReference,
    Contrast,
    ContrastArm,
    InputAsset,
    ProjectManifest,
    Sample,
)
from katana_next.infrastructure.manifest_store import ManifestStore


CONDITIONS = ("Baseline", "Perturbation", "Rescue")
EVENT_TYPES = ("SE", "RI", "A3SS", "A5SS", "MXE")


def create_full_fixture(root: str | Path) -> Path:
    base = Path(root).expanduser().resolve()
    inputs = base / "inputs"
    output = base / "results"
    rmats_root = inputs / "rmats"
    quant_root = inputs / "quant"
    for path in (inputs, output, rmats_root, quant_root):
        path.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(1701)
    samples = _samples(inputs, quant_root)
    contrasts = (
        _contrast("perturbation-vs-baseline", "Perturbation", "Baseline", samples),
        _contrast("rescue-vs-baseline", "Rescue", "Baseline", samples),
        _contrast("rescue-vs-perturbation", "Rescue", "Perturbation", samples),
    )
    event_records = _event_records(80)
    assets: list[InputAsset] = []
    bamlist_dir = rmats_root / "bamlist"
    bamlist_dir.mkdir(parents=True, exist_ok=True)
    for condition in CONDITIONS:
        condition_samples = [sample.sample_id for sample in samples if sample.condition == condition]
        (bamlist_dir / f"{condition}_bamlist.txt").write_text(
            ",".join(condition_samples), encoding="utf-8"
        )
    for contrast in contrasts:
        contrast_dir = (
            rmats_root
            / f"{contrast.numerator.label}-{contrast.denominator.label}"
            / "rmats_post"
        )
        contrast_dir.mkdir(parents=True, exist_ok=True)
        numerator = [sample.sample_id for sample in samples if sample.sample_id in contrast.numerator.sample_ids]
        denominator = [sample.sample_id for sample in samples if sample.sample_id in contrast.denominator.sample_ids]
        for event_type in EVENT_TYPES:
            records = [row for row in event_records if row["event_type"] == event_type]
            frame = _rmats_frame(records, contrast, rng)
            path = contrast_dir / f"{event_type}.MATS.JC.txt"
            frame.to_csv(path, sep="\t", index=False)
            assets.append(
                InputAsset(
                    asset_id=f"rmats::{contrast.contrast_id}::{event_type}",
                    kind="rmats_result",
                    path=path,
                    metadata={
                        "contrast_id": contrast.contrast_id,
                        "event_type": event_type,
                        "mode": "JC",
                        "numerator_sample_ids": numerator,
                        "denominator_sample_ids": denominator,
                    },
                )
            )
        deg_path = inputs / f"{contrast.contrast_id}.deg.tsv"
        _deg_table(contrast, deg_path, rng)
        assets.append(InputAsset(asset_id=f"deg::{contrast.contrast_id}", kind="deg_result", path=deg_path, metadata={"contrast_id": contrast.contrast_id}))
        dtu_path = inputs / f"{contrast.contrast_id}.dtu.tsv"
        _dtu_table(contrast, dtu_path, rng)
        assets.append(InputAsset(asset_id=f"dtu::{contrast.contrast_id}", kind="dtu_result", path=dtu_path, metadata={"contrast_id": contrast.contrast_id}))

    gtf = inputs / "fixture.annotation.gtf"
    transcript_fasta = inputs / "fixture.transcripts.fa"
    genome_fasta = inputs / "fixture.genome.fa"
    _references(gtf, transcript_fasta, genome_fasta)
    pathways = inputs / "fixture.pathways.gmt"
    rbp_sets = inputs / "fixture.rbp.gmt"
    pathways.write_text(
        "RNA_splicing\tfixture\tGENE1\tGENE2\tGENE3\tGENE4\tGENE5\tGENE6\n"
        "Chromatin_regulation\tfixture\tGENE7\tGENE8\tGENE9\tGENE10\tGENE11\n"
        "Stress_response\tfixture\tGENE12\tGENE13\tGENE14\tGENE15\n",
        encoding="utf-8",
    )
    rbp_sets.write_text(
        "SRSF7_targets\tfixture\tGENE1\tGENE3\tGENE5\tGENE8\n"
        "PTBP1_targets\tfixture\tGENE2\tGENE4\tGENE6\tGENE10\n"
        "MBNL_targets\tfixture\tGENE7\tGENE9\tGENE11\tGENE13\n",
        encoding="utf-8",
    )
    assets.extend(
        [
            InputAsset(asset_id="genesets", kind="gene_sets_gmt", path=pathways),
            InputAsset(asset_id="rbp-sets", kind="rbp_gene_sets", path=rbp_sets),
        ]
    )
    sc_metadata, sc_counts = _single_cell_fixture(inputs, rng)
    assets.extend(
        [
            InputAsset(asset_id="single-cell-metadata", kind="single_cell_metadata", path=sc_metadata),
            InputAsset(asset_id="single-cell-counts", kind="single_cell_junction_counts", path=sc_counts),
        ]
    )
    cosplice_metadata = inputs / "cosplice_samples.csv"
    pd.DataFrame(
        {
            "sample_id": [sample.sample_id for sample in samples],
            "condition": [sample.condition for sample in samples],
        }
    ).to_csv(cosplice_metadata, index=False)
    cosplice_config = CoSpliceConfig(
        project_name="Katana Next full fixture",
        rmats_root=rmats_root,
        sample_metadata=cosplice_metadata,
        output_dir=output / "engines" / "cosplice",
        quant_root=quant_root,
        gtf=gtf,
        genome_fasta=genome_fasta,
        transcript_fasta=transcript_fasta,
        genome_build="fixture-v1",
        annotation_version="fixture-v1",
        rmats_mode="JC",
    )
    cosplice_config.network.min_module_size = 8
    cosplice_config.network.tom_dense_max_events = 500
    cosplice_config.network.max_events = 80
    cosplice_config.network.bootstrap_iterations = 10
    cosplice_config.consensus.bootstrap_iterations = 10
    cosplice_config.consensus.min_module_size = 8
    cosplice_config.validation.bootstrap_iterations = 5
    cosplice_path = inputs / "cosplice.yaml"
    cosplice_config.write(cosplice_path)
    assets.append(InputAsset(asset_id="cosplice-config", kind="cosplice_config", path=cosplice_path))
    manifest = ProjectManifest(
        project_id="katana-next-full-fixture",
        name="Katana Next complete acceptance fixture",
        input_root=inputs,
        output_root=output,
        annotation=AnnotationReference(
            genome_build="fixture-v1",
            gtf=gtf,
            genome_fasta=genome_fasta,
            transcript_fasta=transcript_fasta,
        ),
        assets=tuple(assets),
        samples=tuple(samples),
        contrasts=contrasts,
        metadata={
            "purpose": "Software and scientific-contract acceptance only; not a biological result.",
            "direction_contract": "all effects are numerator minus denominator",
            "rscript_path": str(_find_rscript() or ""),
        },
    )
    return ManifestStore().save(manifest, base / "katana.project.json")


def _samples(inputs: Path, quant_root: Path) -> list[Sample]:
    samples = []
    for condition in CONDITIONS:
        for replicate in range(1, 6):
            sample_id = f"{condition.lower()}-{replicate}"
            quant = quant_root / sample_id / "quant.sf"
            quant.parent.mkdir(parents=True, exist_ok=True)
            _quant_table(condition, replicate, quant)
            samples.append(
                Sample(
                    sample_id=sample_id,
                    biological_sample_id=sample_id,
                    condition=condition,
                    batch=f"batch-{replicate}",
                    assets={"quant_sf": quant},
                )
            )
    return samples


def _contrast(identifier: str, numerator: str, denominator: str, samples: list[Sample]) -> Contrast:
    return Contrast(
        contrast_id=identifier,
        label=f"{numerator} vs {denominator}",
        numerator=ContrastArm(label=numerator, sample_ids=tuple(sample.sample_id for sample in samples if sample.condition == numerator)),
        denominator=ContrastArm(label=denominator, sample_ids=tuple(sample.sample_id for sample in samples if sample.condition == denominator)),
    )


def _event_records(count: int) -> list[dict[str, object]]:
    rows = []
    for index in range(count):
        start = 10_000 + index * 500
        rows.append(
            {
                "event_type": EVENT_TYPES[index % len(EVENT_TYPES)],
                "ID": index,
                "GeneID": f"ENSMUSG{index % 20 + 1:011d}.1",
                "geneSymbol": f"GENE{index % 20 + 1}",
                "chr": f"chr{index % 4 + 1}",
                "strand": "+" if index % 2 == 0 else "-",
                "exonStart_0base": start,
                "exonEnd": start + 90,
                "upstreamES": start - 300,
                "upstreamEE": start - 200,
                "downstreamES": start + 300,
                "downstreamEE": start + 400,
                "riExonStart_0base": start - 200,
                "riExonEnd": start + 400,
                "longExonStart_0base": start,
                "longExonEnd": start + 140,
                "shortES": start + 45,
                "shortEE": start + 140,
                "flankingES": start - 260,
                "flankingEE": start - 170,
                "1stExonStart_0base": start,
                "1stExonEnd": start + 80,
                "2ndExonStart_0base": start + 150,
                "2ndExonEnd": start + 235,
                "module": index % 4,
            }
        )
    return rows


def _condition_psi(module: int, condition: str, index: int) -> float:
    patterns = {
        0: {"Baseline": 0.20, "Perturbation": 0.70, "Rescue": 0.58},
        1: {"Baseline": 0.72, "Perturbation": 0.30, "Rescue": 0.68},
        2: {"Baseline": 0.35, "Perturbation": 0.55, "Rescue": 0.75},
        3: {"Baseline": 0.52, "Perturbation": 0.50, "Rescue": 0.54},
    }
    return float(np.clip(patterns[module][condition] + (index % 3 - 1) * 0.015, 0.02, 0.98))


def _rmats_frame(records: list[dict[str, object]], contrast: Contrast, rng: np.random.Generator) -> pd.DataFrame:
    rows = []
    replicate_count = len(contrast.numerator.sample_ids)
    for record in records:
        index = int(record["ID"])
        numerator = np.clip([_condition_psi(int(record["module"]), contrast.numerator.label, index) + rng.normal(0, 0.012) for _ in range(replicate_count)], 0.01, 0.99)
        denominator = np.clip([_condition_psi(int(record["module"]), contrast.denominator.label, index) + rng.normal(0, 0.012) for _ in range(replicate_count)], 0.01, 0.99)
        dpsi = float(np.mean(numerator) - np.mean(denominator))
        total = np.linspace(85, 105, replicate_count, dtype=int) + index % 7
        row = {key: value for key, value in record.items() if key not in {"event_type", "module"}}
        row.update(
            {
                "FDR": 0.001 + (index % 8) * 0.002 if abs(dpsi) >= 0.10 else 0.45,
                "IncLevelDifference": dpsi,
                "IncLevel1": ",".join(f"{value:.4f}" for value in numerator),
                "IncLevel2": ",".join(f"{value:.4f}" for value in denominator),
                "IJC_SAMPLE_1": ",".join(str(int(round(value * depth))) for value, depth in zip(numerator, total, strict=True)),
                "SJC_SAMPLE_1": ",".join(str(int(round((1 - value) * depth))) for value, depth in zip(numerator, total, strict=True)),
                "IJC_SAMPLE_2": ",".join(str(int(round(value * depth))) for value, depth in zip(denominator, total, strict=True)),
                "SJC_SAMPLE_2": ",".join(str(int(round((1 - value) * depth))) for value, depth in zip(denominator, total, strict=True)),
            }
        )
        rows.append(row)
    return pd.DataFrame(rows)


def _deg_table(contrast: Contrast, path: Path, rng: np.random.Generator) -> None:
    rows = []
    for index in range(1, 31):
        effect = (1.8 if index <= 5 else -1.5 if index <= 9 else rng.normal(0, 0.25))
        rows.append({"gene_id": f"ENSMUSG{index:011d}.1", "gene_name": f"GENE{index}", "log2FoldChange": effect, "padj": 0.004 * index if index <= 9 else 0.5})
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)


def _dtu_table(contrast: Contrast, path: Path, rng: np.random.Generator) -> None:
    rows = []
    for index in range(1, 31):
        rows.append({"groupID": f"ENSMUSG{index:011d}.1", "featureID": f"ENSMUST{index * 2:011d}.1", "adj_pval": 0.006 * index if index <= 7 else 0.55})
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)


def _quant_table(condition: str, replicate: int, path: Path) -> None:
    rows = []
    condition_shift = {"Baseline": 0.0, "Perturbation": 3.0, "Rescue": 1.0}[condition]
    for gene in range(1, 31):
        for isoform in (1, 2):
            transcript = (gene - 1) * 2 + isoform
            base = 7.0 + gene % 6
            shift = condition_shift if isoform == 1 and gene <= 8 else -condition_shift / 2 if isoform == 2 and gene <= 8 else 0
            tpm = max(0.05, base + shift + (replicate - 2) * 0.2)
            rows.append({"Name": f"ENSMUST{transcript:011d}.1", "Length": 1000 + transcript, "EffectiveLength": 850 + transcript, "TPM": tpm, "NumReads": tpm * 12})
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)


def _references(gtf: Path, transcript_fasta: Path, genome_fasta: Path) -> None:
    gtf_lines = []
    fasta_lines = []
    for gene in range(1, 31):
        for isoform in (1, 2):
            transcript = (gene - 1) * 2 + isoform
            start = 1000 + gene * 200
            attributes = f'gene_id "ENSMUSG{gene:011d}.1"; transcript_id "ENSMUST{transcript:011d}.1"; gene_name "GENE{gene}";'
            gtf_lines.append(f"chr1\tfixture\ttranscript\t{start}\t{start + 99}\t.\t+\t.\t{attributes}")
            gtf_lines.append(f"chr1\tfixture\texon\t{start}\t{start + 99}\t.\t+\t.\t{attributes}")
            fasta_lines.append(f">ENSMUST{transcript:011d}.1\n" + "ACGT" * 25)
    gtf.write_text("\n".join(gtf_lines) + "\n", encoding="utf-8")
    transcript_fasta.write_text("\n".join(fasta_lines) + "\n", encoding="utf-8")
    genome_fasta.write_text(">chr1\n" + "ACGT" * 20_000 + "\n", encoding="utf-8")


def _single_cell_fixture(inputs: Path, rng: np.random.Generator) -> tuple[Path, Path]:
    metadata_rows = []
    count_rows = []
    for condition in ("Baseline", "Perturbation"):
        for donor_index in range(1, 5):
            donor = f"{condition[:1]}D{donor_index}"
            for cell_type in ("Neuron", "Glia"):
                for cell_index in range(1, 13):
                    cell = f"{donor}-{cell_type[:1]}-{cell_index}"
                    metadata_rows.append({"cell_id": cell, "sample_id": donor, "donor_id": donor, "cell_type": cell_type, "condition": condition, "batch": f"batch-{donor_index % 2 + 1}", "assay": "droplet_3prime"})
                    for event_index in range(1, 21):
                        baseline = 0.25 + (event_index % 5) * 0.10
                        effect = 0.22 if condition == "Perturbation" and event_index <= 8 and cell_type == "Neuron" else -0.16 if condition == "Perturbation" and 9 <= event_index <= 13 and cell_type == "Glia" else 0.0
                        depth = int(rng.poisson(3.0))
                        inclusion = int(rng.binomial(depth, np.clip(baseline + effect, 0.03, 0.97)))
                        count_rows.append({"cell_id": cell, "event_id": f"SC_EVENT_{event_index}", "gene_id": f"ENSMUSG{event_index:011d}", "gene_symbol": f"GENE{event_index}", "event_type": "junction_choice", "annotated": True, "inclusion": inclusion, "exclusion": depth - inclusion})
    metadata = inputs / "single_cell_metadata.tsv"
    counts = inputs / "single_cell_junction_counts.tsv"
    pd.DataFrame(metadata_rows).to_csv(metadata, sep="\t", index=False)
    pd.DataFrame(count_rows).to_csv(counts, sep="\t", index=False)
    return metadata, counts


def _find_rscript() -> Path | None:
    candidates = (
        Path("D:/anaconda3/envs/katana-splicing-tool/Scripts/Rscript.exe"),
        Path.home() / "miniconda3/envs/katana-splicing-tool/Scripts/Rscript.exe",
        Path.home() / "anaconda3/envs/katana-splicing-tool/Scripts/Rscript.exe",
    )
    return next((path.resolve() for path in candidates if path.is_file()), None)


def main() -> None:
    parser = argparse.ArgumentParser(description="Create a complete Katana Next acceptance fixture.")
    parser.add_argument("root")
    arguments = parser.parse_args()
    print(create_full_fixture(arguments.root))


if __name__ == "__main__":
    main()
