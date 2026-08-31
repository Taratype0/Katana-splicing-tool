from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import pandas as pd


EVENT_METADATA_COLUMNS = 10


def _metadata_column_name(value: str) -> str:
    return "_".join(value.strip().split()).replace("/", "_")


def prepare_tcga_spliceseq(
    source: Path,
    output_dir: Path,
    *,
    cohort: str,
    event_type: str = "ES",
) -> dict[str, object]:
    source = source.resolve()
    output_dir = output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    with source.open(encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        sample_columns = [
            column
            for column in header[EVENT_METADATA_COLUMNS:]
            if column.startswith("TCGA_")
        ]
        sample_positions = {
            column: header.index(column)
            for column in sample_columns
        }
        event_rows: list[dict[str, object]] = []
        matrix_rows: list[list[object]] = []
        event_keys: list[str] = []
        clinical_rows: dict[str, dict[str, str]] = {}
        for row in reader:
            if len(row) < EVENT_METADATA_COLUMNS:
                continue
            symbol = row[0].strip()
            as_id = row[1].strip()
            splice_type = row[2].strip().upper()
            if not symbol or not as_id or not splice_type:
                clinical_name = _metadata_column_name(symbol)
                if clinical_name:
                    clinical_rows[clinical_name] = {
                        sample: row[position] if position < len(row) else ""
                        for sample, position in sample_positions.items()
                    }
                continue
            if splice_type != event_type.upper():
                continue
            event_key = f"{symbol}|{splice_type}|{float(as_id):.1f}"
            event_keys.append(event_key)
            matrix_rows.append(
                [
                    row[position] if position < len(row) else ""
                    for position in sample_positions.values()
                ]
            )
            event_rows.append(
                {
                    "event_key": event_key,
                    "gene_symbol": symbol,
                    "as_id": as_id,
                    "event_type": splice_type,
                    "exons": row[3],
                    "from_exon": row[4],
                    "to_exon": row[5],
                    "novel_splice": row[6],
                    "pct_with_values": row[7],
                    "psi_range": row[8],
                    "std_psi": row[9],
                }
            )

    matrix = pd.DataFrame(matrix_rows, index=event_keys, columns=sample_columns)
    matrix = matrix.apply(pd.to_numeric, errors="coerce")
    matrix = matrix.loc[~matrix.index.duplicated(keep="first")]
    tumor_columns = [
        column
        for column in matrix.columns
        if not column.casefold().endswith("_norm")
    ]
    tumor = matrix.loc[:, tumor_columns]
    metadata = pd.DataFrame(event_rows).drop_duplicates("event_key", keep="first")
    sample_metadata = pd.DataFrame(
        {
            "sample_id": sample_columns,
            "cohort": cohort,
            "phenotype": [
                "normal" if sample.casefold().endswith("_norm") else "tumor"
                for sample in sample_columns
            ],
            "sample_type": [
                "adjacent_normal" if sample.casefold().endswith("_norm") else "primary_tumor"
                for sample in sample_columns
            ],
        }
    )
    for clinical_name, values in clinical_rows.items():
        sample_metadata[clinical_name] = sample_metadata["sample_id"].map(values)

    event_type = event_type.upper()
    matrix_path = output_dir / f"{cohort}_{event_type}_PSI_matrix.tsv"
    tumor_path = output_dir / f"{cohort}_{event_type}_PSI_tumor_only.tsv"
    metadata_path = output_dir / f"{cohort}_{event_type}_event_metadata.tsv"
    sample_metadata_path = output_dir / f"{cohort}_sample_metadata.tsv"
    parameters_path = (
        output_dir / "parameters.json"
        if event_type == "ES"
        else output_dir / f"parameters_{event_type}.json"
    )
    matrix.to_csv(matrix_path, sep="\t")
    tumor.to_csv(tumor_path, sep="\t")
    metadata.to_csv(metadata_path, sep="\t", index=False)
    sample_metadata.to_csv(sample_metadata_path, sep="\t", index=False)

    parameters = {
        "source": str(source),
        "cohort": cohort,
        "event_type": event_type,
        "event_key": "gene_symbol|event_type|TCGASpliceSeq_AS_ID",
        "n_events": int(len(matrix)),
        "n_all_samples": int(matrix.shape[1]),
        "n_tumor_samples": int(tumor.shape[1]),
        "n_normal_samples": int(matrix.shape[1] - tumor.shape[1]),
        "normal_sample_rule": "sample column suffix _Norm",
        "phenotype_definition": (
            "phenotype is derived only from the TCGA SpliceSeq _Norm suffix; "
            "downloaded clinical fields are retained as separate columns"
        ),
        "clinical_fields": sorted(clinical_rows),
        "generated_files": [
            str(matrix_path),
            str(tumor_path),
            str(metadata_path),
            str(sample_metadata_path),
        ],
    }
    parameters_path.write_text(json.dumps(parameters, indent=2), encoding="utf-8")
    return parameters


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Convert an official TCGA SpliceSeq download into event-by-sample PSI matrices."
    )
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--cohort", required=True)
    parser.add_argument("--event-type", default="ES")
    args = parser.parse_args()
    parameters = prepare_tcga_spliceseq(
        args.source,
        args.output,
        cohort=args.cohort,
        event_type=args.event_type,
    )
    print(json.dumps(parameters, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
