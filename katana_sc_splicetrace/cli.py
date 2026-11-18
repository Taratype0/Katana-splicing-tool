"""Standalone CLI for katana-sc-splicetrace (2.4.1) with Extensible Motif DB v3.0, Module Bridge & GTF/UniProt Proteome Viz."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import platform
import time
import pandas as pd

from katana_sc_splicetrace.contrast_discovery import (
    compute_sc_splicing_contrast_v2,
    compute_sc_stratified_meta_contrast_v2,
    discover_sc_rbp_programs_v2,
)
from katana_sc_splicetrace.fast_accel import configure_server_threads, fast_score_mean_log_odds
from katana_sc_splicetrace.module_bridge import bridge_splicenet_modules_to_rbps
from katana_sc_splicetrace.protein_consequences import trace_protein_and_nmd_consequences
from katana_sc_splicetrace.rna_forms import IndexedFasta, extract_event_rna_v2
from katana_sc_splicetrace.sequence_and_consequence import (
    build_sc_directional_motif_tensor_v2,
    extract_sc_event_sequences_v2,
    trace_sc_mature_motif_consequences_v2,
)
from katana_sc_splicetrace.target_and_psi_models import rank_rbp_target_events_v2

__version__ = "2.4.1"


def write_nfcore_versions_yml(path: Path | str, process_name: str = "KATANA_SCSPLICETRACE") -> None:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(f'"{process_name}":\n    katana-sc-splicetrace: "2.4.1"\n    python: "{platform.python_version()}"\n', encoding="utf-8")


def add_extensible_library_args(parser: argparse.ArgumentParser) -> None:
    grp = parser.add_argument_group("Motif Library & Extensible Database Options")
    grp.add_argument("--species", choices=["human", "mouse"], default="human", help="Target species namespace (`human` or `mouse`)")
    grp.add_argument("--extra-motifs", nargs="+", type=Path, default=None)
    grp.add_argument("--no-v3-expansion", action="store_true")
    grp.add_argument("--no-user-packs", action="store_true")
    grp.add_argument("--list-rbps", action="store_true")
    grp.add_argument("--filter-rbp", default=None)
    grp.add_argument("--import-motifs", type=Path, default=None)
    grp.add_argument("--export-pack", type=Path, default=None)
    grp.add_argument("--install-pack", type=Path, default=None)
    grp.add_argument("--pack-name", default="custom_user_pack")
    grp.add_argument("--remove-pack", default=None)


def handle_library_management_commands(args: argparse.Namespace) -> tuple[bool, object | None]:
    from katana_sc_splicetrace.library import (
        export_katana_motif_pack,
        get_user_motif_packs_dir,
        load_bundled_human_rbp_library,
        load_curated_human_rbp_library,
        load_external_motif_entries,
    )

    if args.remove_pack:
        pdir = get_user_motif_packs_dir()
        removed = []
        for cand in (pdir / f"{args.remove_pack}.json.gz", pdir / f"{args.remove_pack}.json"):
            if cand.is_file():
                cand.unlink()
                removed.append(str(cand))
        print(json.dumps({"action": "remove_pack", "removed": removed}))
        return True, None

    if args.import_motifs is not None:
        if args.export_pack is None:
            raise ValueError("--export-pack <out.json.gz> is required with --import-motifs")
        ents = load_external_motif_entries(args.import_motifs, default_source=args.pack_name)
        info = export_katana_motif_pack(ents, args.export_pack, pack_name=args.pack_name)
        print(json.dumps({"action": "import_motifs", **info}, indent=2))
        return True, None

    if args.install_pack is not None:
        pdir = get_user_motif_packs_dir()
        pdir.mkdir(parents=True, exist_ok=True)
        dst = pdir / f"{args.pack_name}.json.gz"
        ents = load_external_motif_entries(args.install_pack, default_source=args.pack_name)
        info = export_katana_motif_pack(ents, dst, pack_name=args.pack_name)
        print(json.dumps({"action": "install_pack", **info}, indent=2))
        return True, None

    if getattr(args, "models_json", None) is not None and getattr(args, "assignments_tsv", None) is not None:
        lib = load_curated_human_rbp_library(
            args.models_json, args.assignments_tsv,
            include_v3_expansion=not args.no_v3_expansion, extra_motif_files=args.extra_motifs,
            include_installed_user_packs=not args.no_user_packs, species=args.species,
        )
    else:
        lib = load_bundled_human_rbp_library(
            include_v3_expansion=not args.no_v3_expansion, extra_motif_files=args.extra_motifs,
            include_installed_user_packs=not args.no_user_packs, species=args.species,
        )

    if args.list_rbps:
        df_sum = lib.to_summary_dataframe()
        if args.filter_rbp:
            q = str(args.filter_rbp).upper()
            df_sum = df_sum[
                df_sum["RBP"].str.upper().str.contains(q)
                | df_sum["mouse_symbol"].str.upper().str.contains(q)
                | df_sum["category"].str.upper().str.contains(q)
            ]
        print(df_sum.to_string(index=False))
        return True, lib

    return False, lib


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        prog="katana-sc-splicetrace",
        description="Katana SC-SpliceTrace 2.4.1: Smart-seq3 RBP discovery, mature-mRNA motif, Module Bridge & GTF/UniProt Proteome Consequences.",
    )
    p.add_argument("--version", action="version", version=f"katana-sc-splicetrace {__version__}")
    p.add_argument("--events", type=Path, required=False)
    p.add_argument("--psi", type=Path, required=False)
    p.add_argument("--metadata", type=Path, required=False)
    p.add_argument("--fasta", type=Path, required=False)
    p.add_argument("--gtf", type=Path, default=None, help="Optional GTF/GTF.gz for automated single-cell ORF/NMD/UniProt consequences & 5-style figures")
    p.add_argument("--splicenet-modules", type=Path, default=None, help="Optional sc_cosplice_modules.tsv to bridge SC-SpliceNet modules -> RBP motifs")
    p.add_argument("--threads", type=int, default=0, help="CPU threads for Numba acceleration (0 = auto)")
    p.add_argument("--support", type=Path, default=None)
    p.add_argument("--transpose-psi", action="store_true")
    p.add_argument("--group-col", default="condition")
    p.add_argument("--numerator", default=None)
    p.add_argument("--denominator", default=None)
    p.add_argument("--donor-col", default="donor_id")
    p.add_argument("--query-rbp", default=None)
    p.add_argument("--models-json", type=Path, default=None)
    p.add_argument("--assignments-tsv", type=Path, default=None)
    p.add_argument("--permutations", type=int, default=1999)
    p.add_argument("--seed", type=int, default=20261117)
    p.add_argument("-o", "--out", "--output-dir", dest="out", type=Path, required=False)
    p.add_argument("--emit-versions-yml", type=Path, default=None)
    p.add_argument("--process-name", default="KATANA_SCSPLICETRACE")
    add_extensible_library_args(p)
    args = p.parse_args(argv)

    handled, lib = handle_library_management_commands(args)
    if handled:
        return 0

    if args.emit_versions_yml is not None and args.out is None:
        write_nfcore_versions_yml(args.emit_versions_yml, args.process_name)
        return 0
    if args.events is None or args.fasta is None or args.out is None:
        p.error("--events, --fasta, and --out (or --list-rbps / --import-motifs / --install-pack) are required.")

    active_threads = configure_server_threads(args.threads)
    t0 = time.monotonic()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    events_df = pd.read_csv(args.events, sep="\t")
    rna_df, cov_df = extract_sc_event_sequences_v2(events_df, Path(args.fasta))
    rbp_tensor, ev_ids, rbps, region_names = build_sc_directional_motif_tensor_v2(events_df, rna_df, lib)

    summary: dict = {"version": __version__, "library": lib.name, "n_rbps": lib.n_rbps, "n_models": lib.n_models, "active_threads": int(active_threads), "status": "complete"}
    contrast_df = None
    conseq_df = pd.DataFrame()

    if args.query_rbp is not None:
        resolved_q = lib.resolve_rbp_symbol(args.query_rbp) or args.query_rbp
        t_df = rank_rbp_target_events_v2(resolved_q, ev_ids, events_df["gene_id"].astype(str).to_numpy(), rbp_tensor, rbps, region_names)
        t_df.to_csv(out / f"sc_motiftrace_targets_{resolved_q}.tsv", sep="\t", index=False)
        summary["ranked_targets_for"] = resolved_q

    if args.psi is not None and args.metadata is not None and args.numerator is not None and args.denominator is not None:
        psi_df = pd.read_csv(args.psi, sep="\t", index_col=0)
        sup_df = pd.read_csv(args.support, sep="\t", index_col=0) if args.support else None
        if args.transpose_psi:
            psi_df = psi_df.T
            if sup_df is not None:
                sup_df = sup_df.T
        md = pd.read_csv(args.metadata, sep="\t")
        eg = events_df.set_index("event_id")["gene_id"]
        u_aligned = psi_df.reindex(columns=events_df["event_id"].astype(str))
        s_aligned = sup_df.reindex(columns=events_df["event_id"].astype(str)) if sup_df is not None else None

        d_col = args.donor_col if args.donor_col in md.columns else ("donor" if "donor" in md.columns else None)
        n_donors = md[d_col].nunique() if d_col else 1
        if 2 <= n_donors < 5:
            contrast_df, _, c_info = compute_sc_stratified_meta_contrast_v2(
                u_aligned, md, eg, support=s_aligned, group_col=args.group_col,
                numerator_group=args.numerator, denominator_group=args.denominator, donor_col=d_col, seed=args.seed,
            )
            donor_deltas = None
        else:
            contrast_df, donor_deltas, c_info = compute_sc_splicing_contrast_v2(
                u_aligned, md, eg, support=s_aligned, group_col=args.group_col,
                numerator_group=args.numerator, denominator_group=args.denominator, donor_col=args.donor_col, seed=args.seed,
            )

        conseq_df = trace_sc_mature_motif_consequences_v2(rna_df, contrast_df, lib)
        res = discover_sc_rbp_programs_v2(
            contrast_df, rbp_tensor, rbps, region_names, cov_df,
            donor_delta_df=donor_deltas, mature_consequences_df=conseq_df,
            n_permutations=args.permutations, seed=args.seed,
        )
        res.rbp_ranking.to_csv(out / "sc_splicetrace_rbp_ranking_v24.tsv", sep="\t", index=False)
        res.region_tests.to_csv(out / "sc_splicetrace_region_tests_v24.tsv", sep="\t", index=False)
        res.event_contrasts.to_csv(out / "sc_splicetrace_event_contrasts_v24.tsv", sep="\t", index=False)
        res.mature_consequences.to_csv(out / "sc_splicetrace_mature_consequences_v24.tsv", sep="\t", index=False)
        summary.update({**c_info, **res.diagnostics})

    if args.splicenet_modules is not None:
        mod_in = pd.read_csv(args.splicenet_modules, sep="\t")
        br_df = bridge_splicenet_modules_to_rbps(
            mod_in, events_df, args.fasta, lib,
            extract_rna_fn=extract_event_rna_v2, score_odds_fn=fast_score_mean_log_odds,
            n_permutations=min(args.permutations, 999), seed=args.seed,
        )
        br_df.to_csv(out / "sc_module_x_rbp_bridge_v24.tsv", sep="\t", index=False)
        summary["bridged_module_rbp_rows"] = int(len(br_df))

    if args.gtf is not None:
        if contrast_df is not None:
            elig_c = contrast_df[contrast_df["eligible"]].copy()
            elig_c["abs_z"] = elig_c["z_shift"].abs()
            shifted = elig_c[elig_c["abs_z"] >= 1.5]
            if len(shifted) < 30:
                shifted = elig_c.sort_values("abs_z", ascending=False).head(250)
            else:
                shifted = shifted.sort_values("abs_z", ascending=False).head(350)
            ev_for_prot = events_df.drop(columns=[c for c in ("dpsi",) if c in events_df.columns]).merge(
                shifted[["event_id", "raw_dpsi", "z_shift"]].rename(columns={"raw_dpsi": "dpsi"}),
                on="event_id",
                how="inner",
            )
        else:
            ev_for_prot = events_df.head(350).copy()
            if "dpsi" not in ev_for_prot.columns:
                ev_for_prot["dpsi"] = 0.0

        p_df = trace_protein_and_nmd_consequences(
            ev_for_prot, args.gtf, args.fasta, indexed_fasta_cls=IndexedFasta, species=args.species
        )
        p_df.to_csv(out / "sc_protein_and_nmd_consequences_v24.tsv", sep="\t", index=False)
        summary["protein_consequences_evaluated"] = int(len(p_df))
        summary["nmd_switch_events"] = int(p_df["nmd_switch_class"].isin(["inclusion_triggers_NMD_poison_exon", "skipping_triggers_NMD_decay"]).sum()) if not p_df.empty else 0

        try:
            from katana_sc_splicetrace.protein_viz import (
                plot_causal_alluvial_flow_v24,
                plot_consequence_landscape_v24,
                plot_dual_track_nmd_architecture_v24,
                plot_lollipop_protein_splicing_map_v24,
                plot_side_by_side_inframe_domains_v24,
            )
            if not p_df.empty:
                fig_dir = out / "figures"
                cid_lbl = f"{args.numerator}_vs_{args.denominator}" if (args.numerator and args.denominator) else "sc_panel"
                m_map = {}
                if not conseq_df.empty and "event_id" in conseq_df.columns:
                    from katana_sc_splicetrace.protein_viz import _clean_rbp_label
                    for eid_m, grp_m in conseq_df.groupby("event_id"):
                        m_map[str(eid_m)] = ", ".join(list(dict.fromkeys(f"{_clean_rbp_label(r.RBP)}({r.hit_RNA})" for r in grp_m.itertuples()))[:3])
                plot_consequence_landscape_v24({cid_lbl: p_df}, cid_lbl, fig_dir / "sc_style1_proteome_landscape", dataset_title="Single-Cell Smart-seq3")
                plot_dual_track_nmd_architecture_v24(p_df, m_map, fig_dir / "sc_style2_dual_track_nmd", cid_lbl)
                plot_side_by_side_inframe_domains_v24(p_df, fig_dir / "sc_style3_inframe_domains", cid_lbl)
                plot_lollipop_protein_splicing_map_v24(p_df, fig_dir / "sc_style4_lollipop_protein_map", cid_lbl)
                plot_causal_alluvial_flow_v24(p_df, conseq_df, fig_dir / "sc_style5_causal_alluvial", cid_lbl)
        except Exception:
            pass

    summary["elapsed_seconds"] = round(time.monotonic() - t0, 2)
    (out / "sc_splicetrace_summary_v24.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    if args.emit_versions_yml is not None:
        write_nfcore_versions_yml(args.emit_versions_yml, args.process_name)
    print(json.dumps(summary))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
