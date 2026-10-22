"""A declared downstream bulk-gene universe, independent of upstream model fitting."""
from __future__ import annotations

from types import SimpleNamespace
from pathlib import Path

import pandas as pd

from src.analysis.gene_filter import project_gene_filter, gene_filter_signature
from src.domain.models.gene_filter import GeneFilterConfig


def filter_config(manifest) -> GeneFilterConfig:
    return GeneFilterConfig.from_mapping(manifest.metadata.get("gene_filter", {}).get("config"))


def filter_adapter(manifest):
    settings = manifest.metadata.get("gene_filter", {})
    tools = {}
    annotation = settings.get("annotation_tsv") or next((str(a.path) for a in manifest.assets if a.kind in {"gene_annotation", "annotation_tsv"}), "")
    if annotation:
        tools["annotation_tsv"] = annotation
    gtf = settings.get("gtf") or manifest.annotation.gtf
    if gtf:
        tools["gtf"] = str(gtf)
    expression = settings.get("expression_matrix") or next((str(a.path) for a in manifest.assets if a.kind in {"normalized_counts", "gene_counts", "counts_matrix"}), "")
    return SimpleNamespace(gene_filter=filter_config(manifest), tool_paths=tools,
                           counts_path=Path(expression) if expression else None)


def gene_filter_details(manifest) -> dict:
    adapter = filter_adapter(manifest)
    if not adapter.gene_filter.active:
        return {"active": False, "config": adapter.gene_filter.to_dict()}
    engine = project_gene_filter(adapter)
    audit = engine.annotation_audit()
    excluded = audit.loc[audit.excluded]
    return {"active": True, "config": adapter.gene_filter.to_dict(),
            "identity": gene_filter_signature(adapter),
            "annotation_genes": len(audit), "excluded_annotation_genes": len(excluded),
            "reason_counts": excluded.reasons.str.split(";").explode().value_counts().to_dict(),
            "unknown_biotype_annotation_genes": int(audit.gene_type.eq("").sum()),
            "unknown_chromosome_annotation_genes": int(audit.chr.eq("").sum()),
            "expression_sample_columns": engine.expression_columns,
            "statistics_policy": "Original upstream effects, P values and FDR retained; no upstream refit or subset FDR recalculation.",
            "unknown_annotation_policy": "Retained unless protein-coding-only or expression filtering is enabled; annotation counts are not tested-gene counts."}


def apply_gene_filter(context, frame: pd.DataFrame) -> pd.DataFrame:
    if not filter_config(context.manifest).active or frame.empty:
        return frame
    return project_gene_filter(filter_adapter(context.manifest)).apply(frame)


def validate_gene_filter_identity(manifest) -> None:
    expected = manifest.metadata.get("gene_filter_identity")
    if expected and gene_filter_signature(filter_adapter(manifest)) != expected:
        raise ValueError("Gene-filter annotation or expression source changed. Reconfigure exclusions from the original project.")


def filter_exon_bins(context, frame: pd.DataFrame) -> pd.DataFrame:
    """DEXSeq bin identities encode one or more genes before the last colon."""
    if not filter_config(context.manifest).active or frame.empty:
        return frame
    if "gene_id" in frame or "groupID" in frame:
        return apply_gene_filter(context, frame)
    if "bin_id" not in frame:
        raise ValueError("Gene exclusion requires a gene or DEXSeq bin identity.")
    working = frame.assign(gene_id=frame.bin_id.astype(str).str.rsplit(":", n=1).str[0])
    return apply_gene_filter(context, working).drop(columns="gene_id")


def expression_filter_rules(context, path: Path) -> tuple[Path, str]:
    """Pass only downstream exclusions to R; the full statistical fit stays intact."""
    adapter = filter_adapter(context.manifest)
    engine = project_gene_filter(adapter)
    audit = engine.annotation_audit()
    rules = [{"identifier": row.gene_id, "action": "exclude" if row.excluded else "keep"}
             for row in audit.itertuples()]
    rules.extend({"identifier": gene, "action": "exclude"} for gene in adapter.gene_filter.excluded_genes)
    pd.DataFrame(rules, columns=["identifier", "action"]).to_csv(path, sep="\t", index=False)
    strict = adapter.gene_filter.protein_coding_only or adapter.gene_filter.min_expression > 0
    return path, "strict" if strict else "retain_unknown"


def gene_filter_footnote(manifest) -> str:
    details = gene_filter_details(manifest)
    if not details["active"]:
        return "Gene exclusions: none."
    config = details["config"]
    selected = [label for key, label in (
        ("exclude_mitochondrial", "mitochondrial-encoded genes"),
        ("exclude_pseudogenes", "pseudogenes"),
        ("exclude_rrna", "rRNA"), ("protein_coding_only", "non-protein-coding or unannotated genes")) if config[key]]
    if config["excluded_genes"]:
        selected.append(f"custom list ({len(config['excluded_genes'])})")
    if config["min_expression"] > 0:
        selected.append(f"expression < {config['min_expression']} in the required {config['min_samples']} samples (supplied matrix units)")
    return ("Bulk gene exclusions: " + "; ".join(selected) + ". "
            + f"{details['excluded_annotation_genes']} / {details['annotation_genes']} annotated genes excluded (not tested-gene counts). "
            + "Original upstream effect sizes, P values and FDR are retained; differential models were not refitted. "
            + "PCA/heatmaps and enrichment must be regenerated in this filtered project scope.")
