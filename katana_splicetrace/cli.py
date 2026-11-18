"""Standalone CLI for katana-splicetrace (2.4.3) — Package & nf-core unified interface."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import platform
import time
import pandas as pd

from katana_splicetrace.continuous_rescue import run_bulk_continuous_rescue_v22
from katana_splicetrace.fast_accel import configure_server_threads, fast_score_mean_log_odds
from katana_splicetrace.matching import match_events_v2
from katana_splicetrace.module_bridge import bridge_splicenet_modules_to_rbps
from katana_splicetrace.protein_consequences import trace_protein_and_nmd_consequences
from katana_splicetrace.rmats_io import load_rmats_measurable_events
from katana_splicetrace.rna_forms import IndexedFasta, compute_flank_gc_fast, extract_event_rna_v2, trace_motif_instances_v2
from katana_splicetrace.scanner import run_motiftrace_v2_scan

__version__ = "2.4.3"


def write_nfcore_versions_yml(path: Path | str, process_name: str = "KATANA_SPLICETRACE") -> None:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(f'"{process_name}":\n    katana-splicetrace: "{__version__}"\n    python: "{platform.python_version()}"\n', encoding="utf-8")


def add_extensible_library_args(parser: argparse.ArgumentParser) -> None:
    grp = parser.add_argument_group("Motif Library & Extensible Database Options")
    grp.add_argument("--species", choices=["human", "mouse"], default="human", help="RBP & UniProt species namespace (`human` or `mouse`)")
    grp.add_argument("--extra-motifs", nargs="+", type=Path, default=None, help="External motif files (.meme/.motif/.jaspar/.pfm/.tsv/.csv/.json[.gz])")
    grp.add_argument("--no-v3-expansion", action="store_true", help="Disable v3.0 m6A & ENCODE RBNS/eCLIP expansion")
    grp.add_argument("--no-user-packs", action="store_true", help="Ignore user-installed motif packs in ~/.katana/motif_packs")
    grp.add_argument("--list-rbps", action="store_true", help="Print full RBP catalog summary and exit")
    grp.add_argument("--filter-rbp", default=None, help="Substring/symbol filter when using --list-rbps")
    grp.add_argument("--import-motifs", type=Path, default=None, help="Convert external motif file into a portable Katana Motif Pack")
    grp.add_argument("--export-pack", type=Path, default=None, help="Destination .json[.gz] path for --import-motifs")
    grp.add_argument("--install-pack", type=Path, default=None, help="Install external motif pack into ~/.katana/motif_packs")
    grp.add_argument("--pack-name", default="custom_user_pack", help="Pack name for --install-pack or --remove-pack")
    grp.add_argument("--remove-pack", default=None, help="Remove an installed user motif pack by name and exit")


def handle_library_management_commands(args: argparse.Namespace) -> tuple[bool, object | None]:
    from katana_splicetrace.library import (
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


def _read_table(path: Path) -> pd.DataFrame:
    name = path.name.lower()
    if name.endswith(".parquet"):
        return pd.read_parquet(path)
    if name.endswith(".csv") or name.endswith(".csv.gz"):
        return pd.read_csv(path)
    return pd.read_csv(path, sep="\t")


def _run_gtf_and_5plots(
    ev_df: pd.DataFrame,
    conseq_df: pd.DataFrame,
    gtf_path: Path,
    fasta_path: Path,
    out_dir: Path,
    species: str,
    contrast_label: str = "contrast",
) -> dict:
    prot_df = trace_protein_and_nmd_consequences(
        ev_df.head(500), gtf_path, fasta_path, indexed_fasta_cls=IndexedFasta, species=species
    )
    prot_df.to_csv(out_dir / "protein_and_nmd_consequences_v24.tsv", sep="\t", index=False)
    info = {
        "protein_consequences_evaluated": int(len(prot_df)),
        "nmd_switch_events": int(prot_df["nmd_switch_class"].isin(["inclusion_triggers_NMD_poison_exon", "skipping_triggers_NMD_decay"]).sum()) if not prot_df.empty else 0,
        "uniprot_domain_hits": int(prot_df["disrupted_protein_features"].ne("none").sum()) if not prot_df.empty else 0,
    }
    try:
        from katana_splicetrace.protein_viz import (
            _clean_rbp_label,
            plot_causal_alluvial_flow_v24,
            plot_consequence_landscape_v24,
            plot_dual_track_nmd_architecture_v24,
            plot_lollipop_protein_splicing_map_v24,
            plot_side_by_side_inframe_domains_v24,
        )
        if not prot_df.empty:
            fig_dir = out_dir / "figures"
            m_map = {}
            if not conseq_df.empty and "event_id" in conseq_df.columns:
                for eid_m, grp_m in conseq_df.groupby("event_id"):
                    m_map[str(eid_m)] = ", ".join(list(dict.fromkeys(f"{_clean_rbp_label(r.RBP)}({r.hit_RNA})" for r in grp_m.itertuples()))[:3])
            plot_consequence_landscape_v24({contrast_label: prot_df}, contrast_label, fig_dir / "style1_proteome_consequence_landscape")
            plot_dual_track_nmd_architecture_v24(prot_df, m_map, fig_dir / "style2_nmd_dual_track_architecture", contrast_label)
            plot_side_by_side_inframe_domains_v24(prot_df, fig_dir / "style3_inframe_domain_remodeling", contrast_label)
            plot_lollipop_protein_splicing_map_v24(prot_df, fig_dir / "style4_lollipop_protein_map", contrast_label)
            plot_causal_alluvial_flow_v24(prot_df, conseq_df, fig_dir / "style5_causal_rbp_to_proteome_alluvial", contrast_label)
    except Exception:
        pass
    return info


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="katana-splicetrace",
        description="Katana-SpliceTrace v2.4.3: Bulk AS RBP discovery, mature-mRNA motif, SpliceNet bridge & GTF/UniProt ORF/NMD suite.",
    )
    parser.add_argument("--version", action="version", version=f"katana-splicetrace {__version__}")
    parser.add_argument("--rmats", type=Path, default=None, help="rMATS directory, .tar.gz archive, or measurable events table")
    parser.add_argument("--foreground", type=Path, default=None, help="Discrete foreground event table (.tsv[.gz]/.csv/.parquet)")
    parser.add_argument("--background", type=Path, default=None, help="Discrete background pool event table (.tsv[.gz]/.csv/.parquet)")
    parser.add_argument("--events", type=Path, default=None, help="Standalone event coordinate table for --gtf or --splicenet-modules runs")
    parser.add_argument("--fasta", type=Path, required=False, help="Indexed reference genome FASTA (.fa with .fa.fai)")
    parser.add_argument("--gtf", type=Path, default=None, help="Optional annotation GTF (.gtf[.gz]) for ORF, frameshift, 50-nt EJC NMD & UniProt consequences + 5 plots")
    parser.add_argument("--splicenet-modules", type=Path, default=None, help="Optional SpliceNet module_assignments.tsv for Cross-Module RBP motif bridge")
    parser.add_argument("--threads", type=int, default=0, help="CPU threads for Numba acceleration (0 = auto-detect SLURM_CPUS_PER_TASK / all cores)")
    parser.add_argument("-o", "--output-dir", "--out", dest="output_dir", type=Path, required=False, help="Output directory")
    parser.add_argument("--mode", choices=["auto", "discrete", "rescue"], default="auto")
    parser.add_argument("--models-json", type=Path, default=None)
    parser.add_argument("--assignments-tsv", type=Path, default=None)
    parser.add_argument("--flank-gc-caliper", type=float, default=0.08)
    parser.add_argument("--acceptor-core-nt", type=int, default=3)
    parser.add_argument("--permutations", type=int, default=999)
    parser.add_argument("--seed", type=int, default=20261117)
    parser.add_argument("--emit-versions-yml", type=Path, default=None)
    parser.add_argument("--process-name", default="KATANA_SPLICETRACE")
    add_extensible_library_args(parser)
    args = parser.parse_args(argv)

    handled, lib = handle_library_management_commands(args)
    if handled:
        return 0

    if args.emit_versions_yml is not None and args.output_dir is None:
        write_nfcore_versions_yml(args.emit_versions_yml, args.process_name)
        return 0
    if args.fasta is None or args.output_dir is None:
        parser.error("--fasta and --output-dir (or --out) are required for analysis runs.")

    active_threads = configure_server_threads(args.threads)
    t0 = time.monotonic()
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    bridge_pool_df: pd.DataFrame | None = None

    if args.rmats is not None and (args.mode in {"auto", "rescue"} and args.foreground is None):
        meas_df = load_rmats_measurable_events(args.rmats, args.fasta)
        bridge_pool_df = meas_df
        res22 = run_bulk_continuous_rescue_v22(meas_df, args.fasta, lib, n_permutations=args.permutations, seed=args.seed, threads=active_threads)
        res22.rbp_ranking.to_csv(out_dir / "rbp_ranking_v22.tsv", sep="\t", index=False)
        res22.region_tests.to_csv(out_dir / "region_tests_v22.tsv", sep="\t", index=False)
        res22.event_contrasts.to_csv(out_dir / "evaluated_events_v22.tsv", sep="\t", index=False)
        if not res22.mature_consequences.empty:
            res22.mature_consequences.to_csv(out_dir / "motif_consequence_changes_v22.tsv", sep="\t", index=False)
        receipt = {"version": __version__, "library": lib.name, "n_rbps": lib.n_rbps, "n_models": lib.n_models, "mode": "rescue_v24", "active_threads": int(active_threads), **res22.diagnostics}
        if args.gtf is not None:
            fg_prot_ev = res22.event_contrasts[res22.event_contrasts["_is_fg"]]
            receipt.update(_run_gtf_and_5plots(fg_prot_ev, res22.mature_consequences, args.gtf, args.fasta, out_dir, args.species))

    elif args.foreground is not None and args.background is not None:
        fg = _read_table(args.foreground).copy()
        bg = _read_table(args.background).copy()
        comb_all = pd.concat([fg, bg], ignore_index=True).drop_duplicates("event_id")
        bridge_pool_df = comb_all
        if args.flank_gc_caliper > 0:
            if "flank_GC" not in fg.columns:
                fg["flank_GC"] = fg["event_id"].astype(str).map(compute_flank_gc_fast(fg, args.fasta))
            if "flank_GC" not in bg.columns:
                bg["flank_GC"] = bg["event_id"].astype(str).map(compute_flank_gc_fast(bg, args.fasta))
            caliper_val: float | None = args.flank_gc_caliper
        else:
            caliper_val = None
        matched, _, _ = match_events_v2(fg, bg, flank_gc_caliper=caliper_val)
        m_ok = matched[matched["status"].eq("matched")]
        conseq_fg = pd.DataFrame()
        if args.mode == "auto" and len(m_ok) < 250:
            res22 = run_bulk_continuous_rescue_v22(comb_all, args.fasta, lib, n_permutations=args.permutations, seed=args.seed, threads=active_threads)
            res22.rbp_ranking.to_csv(out_dir / "rbp_ranking_all_v2.tsv", sep="\t", index=False)
            res22.region_tests.to_csv(out_dir / "whole_regions_v2.tsv", sep="\t", index=False)
            conseq_fg = res22.mature_consequences
            receipt = {"version": __version__, "library": lib.name, "n_rbps": lib.n_rbps, "mode": "auto_small_sample_rescue", "active_threads": int(active_threads), **res22.diagnostics}
        else:
            matched.to_csv(out_dir / "matching_v2.tsv", sep="\t", index=False)
            needed = set(m_ok["foreground_event"].astype(str)) | set(m_ok["background_event"].astype(str))
            comb_ev = pd.concat([fg[fg["event_id"].isin(needed)], bg[bg["event_id"].isin(needed)]], ignore_index=True)
            rna_df = extract_event_rna_v2(comb_ev, args.fasta, acceptor_core_nt=args.acceptor_core_nt)
            res = run_motiftrace_v2_scan(matched, rna_df, lib, permutations=args.permutations, seed=args.seed)
            res["whole_regions"].to_csv(out_dir / "whole_regions_v2.tsv", sep="\t", index=False)
            res["clusters"].to_csv(out_dir / "clusters_v2.tsv", sep="\t", index=False)
            res["rbp_ranking_all"].to_csv(out_dir / "rbp_ranking_all_v2.tsv", sep="\t", index=False)
            receipt = {"version": __version__, "library": lib.name, "n_rbps": lib.n_rbps, "mode": "discrete_matched_v2", "active_threads": int(active_threads), "matched_events": int(len(m_ok)), "significant_rbps_fwer05": int((res["rbp_ranking_all"]["rbp_P_FWER"] <= 0.05).sum())}
        if args.gtf is not None:
            if "dpsi" not in fg.columns:
                fg["dpsi"] = 0.0
            receipt.update(_run_gtf_and_5plots(fg, conseq_fg, args.gtf, args.fasta, out_dir, args.species))

    elif args.events is not None and (args.gtf is not None or args.splicenet_modules is not None):
        ev_only = _read_table(args.events).copy()
        if "dpsi" not in ev_only.columns:
            ev_only["dpsi"] = 0.0
        bridge_pool_df = ev_only
        receipt = {"version": __version__, "library": lib.name, "mode": "consequences_or_bridge_only", "active_threads": int(active_threads)}
        if args.gtf is not None:
            receipt.update(_run_gtf_and_5plots(ev_only, pd.DataFrame(), args.gtf, args.fasta, out_dir, args.species))
    else:
        parser.error("Provide --rmats, --foreground/--background, --events (with --gtf or --splicenet-modules), or a library flag (--list-rbps).")

    if args.splicenet_modules is not None:
        ev_bridge = _read_table(args.events) if args.events is not None else bridge_pool_df
        if ev_bridge is not None and not ev_bridge.empty:
            mod_in = _read_table(args.splicenet_modules)
            br_df = bridge_splicenet_modules_to_rbps(
                mod_in, ev_bridge, args.fasta, lib,
                extract_rna_fn=extract_event_rna_v2, score_odds_fn=fast_score_mean_log_odds,
                n_permutations=args.permutations, seed=args.seed,
            )
            br_df.to_csv(out_dir / "module_x_rbp_bridge_v24.tsv", sep="\t", index=False)
            receipt["bridged_module_rbp_rows"] = int(len(br_df))

    receipt["elapsed_seconds"] = round(time.monotonic() - t0, 2)
    (out_dir / "run_receipt.json").write_text(json.dumps(receipt, indent=2), encoding="utf-8")
    if args.emit_versions_yml is not None:
        write_nfcore_versions_yml(args.emit_versions_yml, args.process_name)
    print(json.dumps(receipt))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
