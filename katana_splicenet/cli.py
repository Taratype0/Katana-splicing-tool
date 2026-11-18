"""Standalone CLI for katana-splicenet (v2.3.2) — pure Python out-of-the-box, nf-core compatible."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import time
import pandas as pd

from katana_splicenet.evaluation import enrich_module_genes_v2, score_modules_v2
from katana_splicenet.network import consensus_core_audit, fit_wgcna_v2
from katana_splicenet.qc_and_selection import select_events_v2, unit_robust_audit_v2

__version__ = "2.3.2"
def write_nfcore_versions_yml(path: Path | str, process_name: str = "KATANA_SPLICENET") -> None:
    """Emit standard nf-core `versions.yml` recording tool and Python versions."""
    import platform
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    yml = (
        f'"{process_name}":\n'
        f'    katana-splicenet: "2.3.2"\n'
        f'    python: "{platform.python_version()}"\n'
    )
    p.write_text(yml, encoding="utf-8")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="katana-splicenet",
        description="Katana Bulk CoSplice v2.3.2: Excess-BioVar debiased signed TOM co-splicing network & module discovery.",
    )
    parser.add_argument("--version", action="version", version=f"katana-splicenet {__version__}")
    parser.add_argument("--usage", type=Path, required=False, help="Relative usage / PSI matrix TSV (rows=events, cols=samples)")
    parser.add_argument("--samples", type=Path, required=False, help="Sample metadata TSV (sample_id, biological_sample_id)")
    parser.add_argument("--support", type=Path, default=None, help="Optional total read-support matrix TSV")
    parser.add_argument("--annotation", type=Path, default=None, help="Optional event-to-gene TSV (event_key, gene_id)")
    parser.add_argument("--gmt", type=Path, default=None, help="Optional MSigDB/GO .gmt pathway file")
    parser.add_argument("-o", "--output-dir", "--out", dest="output_dir", type=Path, required=False, help="Output directory")
    parser.add_argument("--threads", type=int, default=0, help="CPU thread hint for BLAS operations")
    parser.add_argument("--feature-cap", type=int, default=5000)
    parser.add_argument("--min-module-size", type=int, default=20)
    parser.add_argument("--ranking-metric", choices=["excess_bio_variance", "variance"], default="excess_bio_variance")
    parser.add_argument("--power", type=int, default=0, help="Signed TOM soft-threshold power (0 = auto scale-free pick_soft_threshold_v2)")
    parser.add_argument("--backend", choices=["python", "r-wgcna"], default="python")
    parser.add_argument("--subsamples", type=int, default=0)
    parser.add_argument("--emit-versions-yml", type=Path, default=None, help="Write nf-core compatible versions.yml")
    parser.add_argument("--process-name", default="KATANA_SPLICENET")
    args = parser.parse_args(argv)

    if args.emit_versions_yml is not None and args.output_dir is None:
        write_nfcore_versions_yml(args.emit_versions_yml, args.process_name)
        return 0
    if args.usage is None or args.samples is None or args.output_dir is None:
        parser.error("--usage, --samples, and --output-dir are required.")

    t0 = time.monotonic()
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    usage = pd.read_csv(args.usage, sep="\t", index_col=0)
    meta = pd.read_csv(args.samples, sep="\t", dtype={"sample_id": str, "biological_sample_id": str})
    support = pd.read_csv(args.support, sep="\t", index_col=0) if args.support is not None else None

    masked, audit = unit_robust_audit_v2(usage, meta, support)
    selected, sel_audit = select_events_v2(masked, audit, feature_cap=args.feature_cap, ranking_metric=args.ranking_metric)
    sel_audit.to_csv(out_dir / "event_selection_audit_v2.tsv.gz", sep="\t", index_label="event_key")

    assignments = fit_wgcna_v2(
        selected, out_dir / "wgcna_fit", power=args.power, min_module_size=args.min_module_size, backend=args.backend
    )
    stability_metrics: dict[str, float] = {}
    if args.subsamples > 0:
        assignments, stability_metrics = consensus_core_audit(
            usage, meta, support, assignments,
            feature_cap=args.feature_cap, ranking_metric=args.ranking_metric,
            min_module_size=args.min_module_size, n_subsamples=args.subsamples,
        )

    assignments.to_csv(out_dir / "module_assignments_v2.tsv", sep="\t", index=False)
    scores, _, mod_summary, basis = score_modules_v2(selected, assignments)
    scores.to_csv(out_dir / "module_scores_v2.tsv", sep="\t", index_label="sample_id")
    mod_summary.to_csv(out_dir / "module_summary_v2.tsv", sep="\t", index=False)
    basis.to_csv(out_dir / "module_pc1_basis_v2.tsv", sep="\t", index=False)

    if args.annotation is not None:
        ann_df = pd.read_csv(args.annotation, sep="\t", dtype=str).set_index("event_key")
        gene_map = ann_df["gene_id"]
        mod_genes = assignments[assignments["module"].ne(0)].copy()
        mod_genes["gene_id"] = mod_genes["event_key"].map(gene_map)
        mod_genes.to_csv(out_dir / "module_genes_v2.tsv", sep="\t", index=False)
        if args.gmt is not None:
            enr = enrich_module_genes_v2(assignments, gene_map, args.gmt)
            enr.to_csv(out_dir / "module_pathway_enrichment_v2.tsv", sep="\t", index=False)

    receipt = {
        "version": __version__,
        "backend": args.backend,
        "input_events": int(len(usage)),
        "qc_eligible_events": int(audit["eligible"].sum()),
        "selected_events": int(len(selected)),
        "non_grey_modules": int(assignments.loc[assignments["module"].ne(0), "module"].nunique()),
        "assigned_events": int(assignments["module"].ne(0).sum()),
        "stability": stability_metrics,
        "elapsed_seconds": round(time.monotonic() - t0, 2),
    }
    (out_dir / "run_receipt.json").write_text(json.dumps(receipt, indent=2), encoding="utf-8")
    if args.emit_versions_yml is not None:
        write_nfcore_versions_yml(args.emit_versions_yml, args.process_name)
    print(json.dumps(receipt))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
