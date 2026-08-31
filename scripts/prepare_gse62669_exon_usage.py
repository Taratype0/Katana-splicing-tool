from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


def prepare_exon_usage(
    source: Path,
    output_matrix: Path,
    output_metadata: Path,
    output_audit: Path,
    *,
    min_exons_per_gene: int = 3,
    min_gene_total: float = 10.0,
    min_feature_sd: float = 1e-8,
) -> None:
    counts = pd.read_csv(source, sep="\t", index_col=0, low_memory=False)
    counts = counts.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    counts.index = counts.index.astype(str)
    genes = counts.index.to_series().str.rsplit("_", n=1).str[0]
    exon_counts = genes.value_counts()
    retained_genes = exon_counts.loc[exon_counts.ge(min_exons_per_gene)].index
    counts = counts.loc[genes.isin(retained_genes)]
    genes = genes.loc[counts.index]

    usage_parts: list[pd.DataFrame] = []
    audit_rows: list[dict[str, object]] = []
    for gene_id, gene_counts in counts.groupby(genes, sort=False):
        totals = gene_counts.sum(axis=0)
        usable = totals.ge(min_gene_total)
        if int(usable.sum()) < 3:
            continue
        usage = gene_counts.div(totals.where(totals.gt(0)), axis=1)
        usage = usage.loc[usage.std(axis=1).ge(min_feature_sd)]
        if usage.empty:
            continue
        usage.index = [f"{gene_id}|{event_id}" for event_id in usage.index]
        usage_parts.append(usage)
        audit_rows.append(
            {
                "gene_id": gene_id,
                "n_exons": len(gene_counts),
                "n_variable_exons": len(usage),
                "n_samples_with_gene_total_at_least_cutoff": int(usable.sum()),
            }
        )
    if not usage_parts:
        raise ValueError("No exon-usage features passed the configured filters.")
    matrix = pd.concat(usage_parts)
    matrix.index.name = "event_key"
    output_matrix.parent.mkdir(parents=True, exist_ok=True)
    matrix.to_csv(output_matrix, sep="\t")
    pd.DataFrame(
        {
            "sample_id": matrix.columns,
            "condition": "HSCC_brain",
            "cohort": "GSE62669",
        }
    ).to_csv(output_metadata, sep="\t", index=False)
    pd.DataFrame(audit_rows).to_csv(output_audit, sep="\t", index=False)
    parameters = {
        "source": str(source.resolve()),
        "output_matrix": str(output_matrix.resolve()),
        "definition": (
            "Each feature is one annotated exon count divided by the sum of annotated exon counts "
            "for the same Ensembl gene and sample."
        ),
        "scientific_boundary": (
            "This is an exon-usage fraction proxy derived from the paper's normalized exon-count "
            "matrix. It is not rMATS PSI and is not the original gene-level Mantel statistic."
        ),
        "min_exons_per_gene": min_exons_per_gene,
        "min_gene_total": min_gene_total,
        "min_feature_sd": min_feature_sd,
        "n_features": len(matrix),
        "n_samples": matrix.shape[1],
        "n_genes": int(matrix.index.to_series().str.split("|", regex=False).str[0].nunique()),
    }
    output_matrix.with_suffix(".parameters.json").write_text(
        json.dumps(parameters, indent=2), encoding="utf-8"
    )


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Prepare a transparent GSE62669 exon-usage proxy matrix."
    )
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output-matrix", type=Path, required=True)
    parser.add_argument("--output-metadata", type=Path, required=True)
    parser.add_argument("--output-audit", type=Path, required=True)
    parser.add_argument("--min-exons-per-gene", type=int, default=3)
    parser.add_argument("--min-gene-total", type=float, default=10.0)
    parser.add_argument(
        "--min-feature-sd",
        type=float,
        default=1e-8,
        help="Minimum across-sample exon-usage standard deviation; use 0 to retain invariant exons.",
    )
    args = parser.parse_args()
    prepare_exon_usage(
        args.source,
        args.output_matrix,
        args.output_metadata,
        args.output_audit,
        min_exons_per_gene=args.min_exons_per_gene,
        min_gene_total=args.min_gene_total,
        min_feature_sd=args.min_feature_sd,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
