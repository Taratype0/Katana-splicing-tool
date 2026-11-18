"""Standalone CLI for katana-sc-splicenet (v2.2.2) — Smart-seq3 support-precision co-splicing network & modules."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import pandas as pd

from katana_sc_splicenet.evaluation import enrich_sc_module_genes_v2
from katana_sc_splicenet.network import fit_sc_cosplice_v2
from katana_sc_splicenet.qc_and_selection import (
    attach_split_state_reliability_and_psix_v2,
    audit_sc_events_v2,
    select_sc_events_v2,
)

__version__ = "2.2.2"
def write_nfcore_versions_yml(path: Path | str, process_name: str = "KATANA_SC_COSPLICE") -> None:
    """Emit standard nf-core `versions.yml` recording tool and Python versions."""
    import platform
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    yml = (
        f'"{process_name}":\n'
        f'    katana-sc-splicenet: "2.2.2"\n'
        f'    python: "{platform.python_version()}"\n'
    )
    p.write_text(yml, encoding="utf-8")


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        prog="katana-sc-splicenet",
        description="Katana SC CoSplice v2.2.2: Smart-seq3 Hybrid-Manifold ExcessBioVar & support-precision TOM co-splicing network.",
    )
    p.add_argument("--version", action="version", version=f"katana-sc-splicenet {__version__}")
    p.add_argument("--usage", type=Path, required=False, help="TSV(.gz): cell_id rows x event_id columns observed usage")
    p.add_argument("--metadata", type=Path, required=False, help="TSV: QC-passed cell metadata (cell_id, donor_id, cell_type, condition, qc_pass)")
    p.add_argument("--event-genes", type=Path, required=False, help="TSV: event_id, gene_id[, gene_symbol]")
    p.add_argument("--support", type=Path, default=None, help="Optional TSV(.gz): cell_id x event_id count denominators")
    p.add_argument("--host-expr", type=Path, default=None, help="Optional TSV(.gz): cell_id x gene_id log-expression for Hybrid Manifold selector")
    p.add_argument("--select-budget", type=int, default=0, help="If > 0, select top-K splicing events first using --selector")
    p.add_argument(
        "--selector",
        choices=["hybrid_manifold_excess_bio_var", "reliability_excess_bio_var", "excess_bio_var"],
        default="hybrid_manifold_excess_bio_var",
    )
    p.add_argument("--level", default="cell_state_distribution", choices=["between_donor", "cell_state_distribution", "within_state_cell"])
    p.add_argument("--top-fraction", type=float, default=0.05)
    p.add_argument("--min-module-size", type=int, default=10)
    p.add_argument("--repeats", type=int, default=8)
    p.add_argument("--threads", type=int, default=0, help="CPU thread hint for BLAS matrix kernels")
    p.add_argument("--seed", type=int, default=20261117)
    p.add_argument("--gmt", type=Path, default=None, help="Optional GMT pathway file for module enrichment")
    p.add_argument("-o", "--out", "--output-dir", dest="out", type=Path, required=False, help="Output directory")
    p.add_argument("--emit-versions-yml", type=Path, default=None)
    p.add_argument("--process-name", default="KATANA_SCSPLICENET")
    args = p.parse_args(argv)

    if args.emit_versions_yml is not None and args.out is None:
        write_nfcore_versions_yml(args.emit_versions_yml, args.process_name)
        return 0
    if args.usage is None or args.metadata is None or args.event_genes is None or args.out is None:
        p.error("--usage, --metadata, --event-genes, and --out are required.")

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    usage = pd.read_csv(args.usage, sep="\t", index_col=0)
    metadata = pd.read_csv(args.metadata, sep="\t")
    eg_df = pd.read_csv(args.event_genes, sep="\t")
    support = pd.read_csv(args.support, sep="\t", index_col=0) if args.support else None
    eg_ser = eg_df.set_index("event_id")["gene_id"].astype(str)

    if args.select_budget > 0 and usage.shape[1] > args.select_budget:
        host_df = pd.read_csv(args.host_expr, sep="\t", index_col=0) if args.host_expr else None
        aud = audit_sc_events_v2(usage, eg_ser, support=support)
        aud = attach_split_state_reliability_and_psix_v2(aud, usage, metadata, support=support, host_log_expr=host_df, seed=args.seed)
        sel = select_sc_events_v2(aud, budget=args.select_budget, max_per_gene=2, ranking_metric=args.selector)
        keep_ids = sel["event_id"].astype(str).tolist()
        usage = usage.loc[:, keep_ids]
        if support is not None:
            support = support.loc[:, keep_ids]
        eg_df = eg_df[eg_df["event_id"].astype(str).isin(set(keep_ids))].copy()
        sel.to_csv(out / "sc_cosplice_selected_events_v21.tsv", sep="\t", index=False)

    net = fit_sc_cosplice_v2(
        usage=usage, metadata=metadata, event_genes=eg_df, support=support,
        level=args.level, top_fraction=args.top_fraction, repeats=args.repeats,
        random_seed=args.seed, min_module_size=args.min_module_size, use_r_dynamic_cut=False,
    )
    net.edges.to_csv(out / "sc_cosplice_edges_v21.tsv", sep="\t", index=False)
    net.modules.to_csv(out / "sc_cosplice_modules_v21.tsv", sep="\t", index=False)
    net.event_audit.to_csv(out / "sc_cosplice_event_audit_v21.tsv", sep="\t", index=False)

    if args.gmt is not None:
        sym_map = eg_df.set_index("event_id")["gene_symbol"] if "gene_symbol" in eg_df.columns else None
        enr = enrich_sc_module_genes_v2(net.modules, Path(args.gmt), gene_symbol_map=sym_map)
        enr.to_csv(out / "sc_cosplice_module_enrichment_v21.tsv", sep="\t", index=False)

    diag = {"version": __version__, **net.diagnostics}
    (out / "sc_cosplice_summary_v21.json").write_text(json.dumps(diag, indent=2), encoding="utf-8")
    if args.emit_versions_yml is not None:
        write_nfcore_versions_yml(args.emit_versions_yml, args.process_name)
    print(json.dumps({
        "version": __version__,
        "status": diag["status"],
        "selected_edges": diag["selected_edges"],
        "discovered_modules": diag["discovered_modules"],
    }))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
