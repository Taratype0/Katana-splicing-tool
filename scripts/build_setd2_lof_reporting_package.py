from __future__ import annotations

import argparse
import base64
import csv
import hashlib
import html
import json
import platform
import shutil
import subprocess
import sys
import zipfile
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


DEFAULT_SOURCE = Path(
    r"E:\Project\RNASEQ\Result\Setd2_results\katana_output"
    r"\09_cosplice_network\06i_lof_hypothesis_closure"
)
DEFAULT_OUTPUT = DEFAULT_SOURCE.parent / "06j_lof_reporting_package"
EVENT_TYPES = ("SE", "RI", "A3SS", "A5SS", "MXE")
COMPARISONS = {
    "DP": {
        "MUT": "MUT_DP-WT_DP",
        "CLK3": "WT_DP_CLK3-WT_DP",
        "Pla": "WT_DP_Pla-WT_DP",
    },
    "noDP": {
        "MUT": "MUT_NO_DP-WT_NO_DP",
        "CLK3": "WT_NO_DP_CLK3-WT_NO_DP",
        "Pla": "WT_NO_DP_Pla-WT_NO_DP",
    },
}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _git(repo: Path, *args: str) -> str:
    try:
        completed = subprocess.run(
            ["git", *args],
            cwd=repo,
            check=True,
            capture_output=True,
            text=True,
        )
        return completed.stdout.strip()
    except Exception as exc:  # provenance must not prevent report generation
        return f"unavailable: {exc}"


def _copy_tree_files(source: Path, destination: Path) -> list[tuple[Path, Path]]:
    copied: list[tuple[Path, Path]] = []
    for path in sorted(source.rglob("*")):
        if not path.is_file():
            continue
        relative = path.relative_to(source)
        target = destination / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(path, target)
        copied.append((path, target))
    return copied


def _upstream_inputs(source: Path) -> list[tuple[str, Path]]:
    project_root = source.parents[2]
    network_root = project_root / "katana_output" / "09_cosplice_network"
    exact_root = network_root / "exact_wgcna_full_stable_mut_effect"
    rmats_root = project_root / "Setd2_Splicing_genomebam" / "star_salmon" / "rmats"
    inputs = [
        ("module_assignment", exact_root / "event_module_assignments.tsv"),
        ("module_score", exact_root / "module_program_scores.tsv"),
        ("module_pathway", exact_root / "module_go_kegg_enrichment.tsv"),
        ("stable_psi", network_root / "03_stable_as_qc" / "stable_psi_matrix.tsv"),
        ("sample_metadata", network_root / "00_config" / "factorial_sample_metadata.tsv"),
        ("event_metadata", network_root / "02_matrix" / "event_metadata.tsv"),
        (
            "rbp_expression",
            network_root / "06e_module_interpretation" / "rbp_gene_expression_matrix.tsv",
        ),
        (
            "rbp_motif_hits",
            network_root / "06e_module_interpretation" / "event_rbp_motif_hits.tsv",
        ),
    ]
    comparison_ids = sorted(
        {comparison for mapping in COMPARISONS.values() for comparison in mapping.values()}
    )
    inputs.extend(
        (
            f"rmats_{event_type}_{comparison}",
            rmats_root
            / comparison
            / "rmats_post"
            / f"{event_type}.MATS.JC.txt",
        )
        for comparison in comparison_ids
        for event_type in EVENT_TYPES
    )
    return inputs


def _write_upstream_input_manifest(source: Path, destination: Path) -> None:
    rows = []
    for role, path in _upstream_inputs(source):
        exists = path.is_file()
        stat = path.stat() if exists else None
        rows.append(
            {
                "input_role": role,
                "absolute_path": str(path.resolve()) if exists else str(path),
                "exists": exists,
                "bytes": stat.st_size if stat else "",
                "modified_utc": (
                    datetime.fromtimestamp(stat.st_mtime, tz=timezone.utc).isoformat()
                    if stat
                    else ""
                ),
                "sha256": _sha256(path) if exists else "",
            }
        )
    with destination.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def _fmt_p(value: object) -> str:
    number = float(value)
    return f"{number:.2e}" if number < 0.001 else f"{number:.3f}"


def _fmt(value: object, digits: int = 3) -> str:
    return f"{float(value):.{digits}f}"


def _data_uri(path: Path) -> str:
    mime = "image/png" if path.suffix.lower() == ".png" else "application/octet-stream"
    encoded = base64.b64encode(path.read_bytes()).decode("ascii")
    return f"data:{mime};base64,{encoded}"


def _markdown_table(frame: pd.DataFrame) -> str:
    return frame.to_markdown(index=False)


def _html_table(frame: pd.DataFrame) -> str:
    return frame.to_html(index=False, border=0, classes="data-table", escape=True)


def _event_summary(source: Path) -> pd.DataFrame:
    frame = pd.read_csv(source / "exact_shared_event_lof_summary.tsv", sep="\t")
    out = pd.DataFrame(
        {
            "Context": frame["context"],
            "Reference": frame["inhibitor"],
            "MUT significant": frame["n_MUT_significant_total"].astype(int),
            "Reference significant": frame["n_inhibitor_significant_total"].astype(int),
            "Exact shared": frame["n_shared_significant"].astype(int),
            "MUT covered": frame["MUT_contained_pct"].map(lambda x: f"{x:.1f}%"),
            "Overlap P": frame["overlap_enrichment_pvalue"].map(_fmt_p),
            "Overlap FDR": frame["overlap_enrichment_FDR"].map(_fmt_p),
            "Same direction": frame["shared_same_direction_pct"].map(lambda x: f"{x:.1f}%"),
            "dPSI Spearman": frame["shared_dPSI_spearman"].map(lambda x: f"{x:.3f}"),
            "Median |delta dPSI|": frame["median_shared_abs_delta_dPSI"].map(
                lambda x: f"{x:.3f}"
            ),
        }
    )
    return out


def _core_summary(source: Path) -> pd.DataFrame:
    genes = pd.read_csv(source / "exact_core_gene_membership.tsv", sep="\t")
    events = pd.read_csv(source / "exact_core_event_type_summary.tsv", sep="\t")
    gene_counts = genes.groupby(["context", "core_set"])["gene_id"].nunique()
    event_counts = events.groupby(["context", "core_set"])["n_events"].sum()
    rows = []
    labels = ["CLK3_like", "PlaB_like", "MUT_CLK3_PlaB_triple"]
    for context in ["DP", "noDP"]:
        for core in labels:
            rows.append(
                {
                    "Context": context,
                    "Core": core,
                    "Events": int(event_counts.get((context, core), 0)),
                    "Unique genes": int(gene_counts.get((context, core), 0)),
                }
            )
    return pd.DataFrame(rows)


def _functional_summary(source: Path) -> pd.DataFrame:
    frame = pd.read_csv(source / "exact_core_go_kegg_enrichment.tsv", sep="\t")
    significant = frame[pd.to_numeric(frame["FDR"], errors="coerce") <= 0.05].copy()
    if significant.empty:
        return pd.DataFrame(columns=["Context", "Core", "Source", "Term", "FDR", "Genes"])
    significant = significant.sort_values(["FDR", "intersection_size"]).head(12)
    return pd.DataFrame(
        {
            "Context": significant["context"],
            "Core": significant["core_set"],
            "Source": significant["source"],
            "Term": significant["term_name"],
            "FDR": significant["FDR"].map(_fmt_p),
            "Genes": significant["intersection_size"].astype(int),
        }
    )


def _module_summary(source: Path) -> pd.DataFrame:
    frame = pd.read_csv(source / "exact_module_annotation_summary.tsv", sep="\t")
    frame = frame[frame["module"].isin(["M03", "M06"])].copy()
    return frame[
        [
            "module",
            "interpretation_label",
            "n_events",
            "n_genes",
            "top_hub_genes",
            "top_pathways",
            "top_expression_correlated_RBPs",
        ]
    ].rename(
        columns={
            "module": "Module",
            "interpretation_label": "Interpretation",
            "n_events": "Events",
            "n_genes": "Genes",
            "top_hub_genes": "Top hub genes",
            "top_pathways": "Top pathways",
            "top_expression_correlated_RBPs": "Prioritized RBPs",
        }
    )


def _occlusion_summary(source: Path) -> pd.DataFrame:
    frame = pd.read_csv(source / "exact_module_occlusion_convergence.tsv", sep="\t")
    frame = frame[frame["module"].isin(["M03", "M06"])].copy()
    return pd.DataFrame(
        {
            "Module": frame["module"],
            "Context": frame["DP_status"],
            "Inhibitor": frame["inhibitor"],
            "WT effect": frame["WT_inhibitor_effect"].map(lambda x: f"{x:+.3f}"),
            "MUT effect": frame["MUT_inhibitor_effect"].map(lambda x: f"{x:+.3f}"),
            "Interaction FDR": frame["occlusion_interaction_FDR"].map(_fmt_p),
            "Formal occlusion": frame["formal_occlusion_supported"].map(
                lambda x: "yes" if bool(x) else "no"
            ),
            "Descriptive convergence": frame["descriptive_convergence"].map(
                lambda x: "yes" if bool(x) else "no"
            ),
        }
    )


def _rbp_summary(source: Path) -> pd.DataFrame:
    frame = pd.read_csv(source / "exact_module_rbp_regulators.tsv", sep="\t")
    frame = frame[frame["module"].isin(["M03", "M06"])].copy()
    frame = frame.sort_values(["module", "FDR"]).groupby("module", as_index=False).head(8)
    return pd.DataFrame(
        {
            "Module": frame["module"],
            "RBP": frame["RBP"],
            "Expression-module r": frame["expression_module_correlation"].map(
                lambda x: f"{x:+.3f}"
            ),
            "FDR": frame["FDR"].map(_fmt_p),
            "Module gene": frame["RBP_is_module_gene"].map(lambda x: "yes" if bool(x) else "no"),
            "max |kME|": frame["max_abs_kME"].map(
                lambda x: "" if pd.isna(x) else f"{x:.3f}"
            ),
        }
    )


def _build_report_markdown(
    event: pd.DataFrame,
    cores: pd.DataFrame,
    functional: pd.DataFrame,
    modules: pd.DataFrame,
    occlusion: pd.DataFrame,
    rbp: pd.DataFrame,
    author: str,
) -> str:
    return f"""# Setd2 MUT LOF-compatible splicing evidence report

Prepared by: {author}
Generated: {datetime.now().astimezone().isoformat(timespec='seconds')}

## Executive answer

Setd2 MUT contains a statistically non-random, directionally concordant subset of the CLK3- and PlaB-inhibitor splicing responses. This supports **partial functional convergence compatible with a LOF component**, but it does not prove complete inhibitor phenocopy, a saturated linear pathway, or definitive SETD2 LOF. Combination-module tests show retained inhibitor responses and no systematic formal occlusion.

## Evidence pyramid

1. **Direct event evidence — strongest:** exact shared significant rMATS events are far above chance, with 86.7-96.4% directional agreement and strong dPSI correlation.
2. **Core functional interpretation — limited:** CLK3-like, PlaB-like and triple cores are real event sets, but GO/KEGG reveals no broad shared splicing pathway after FDR correction.
3. **Combination/epistasis — does not prove one saturated pathway:** inhibitor responses remain in MUT; M03 DP+PlaB is enhanced rather than occluded.
4. **Co-splicing modules — supporting context:** M03/M06 are programs discovered across all 36 samples; treatment and DP dominate, while genotype-by-DP is significant.
5. **Core-to-module localization — orthogonal layers:** inhibitor-like cores do not significantly concentrate in M03/M06.
6. **RBP prioritization — hypothesis generation only:** correlations and module membership prioritize candidates but do not establish MUT-specific or causal regulators.

## 1. Exact event overlap

Significant AS rule: `rMATS FDR <= 0.05` and `abs(standardized dPSI) >= 0.10`. Events are matched by exact canonical rMATS identity: gene ID, event type, chromosome, strand and event coordinates.

{_markdown_table(event)}

Interpretation: overlap enrichment FDR ranges from approximately `1.57e-40` to `7.05e-100`. The common subset is strongly concordant, but only 20-34% of the MUT signature is contained in either inhibitor reference.

## 2. CLK3-like, PlaB-like and triple cores

- CLK3-like: exact event significant in MUT-vs-WT and WT+CLK3-vs-WT in the same DP context, with the same standardized dPSI sign.
- PlaB-like: corresponding MUT and PlaB events significant and same-direction.
- Triple core: the exact event is significant in MUT, CLK3 and PlaB effects and all three signs agree.

{_markdown_table(cores)}

## 3. Functional interpretation of the cores

GO/KEGG uses unique genes and a custom background restricted to genes represented by exact AS events testable in the three relevant contrasts.

{_markdown_table(functional)}

Only the DP CLK3-like core passes FDR 0.05, driven by four genes in single-strand-annealing DNA repair terms. PlaB-like and triple cores do not have robust FDR-significant GO/KEGG terms. This is a negative result that prevents over-claiming a unified functional pathway.

## 4. M03 and M06 are all-sample co-splicing programs

{_markdown_table(modules)}

The module membership was estimated across WT/MUT, DP/noDP and Baseline/CLK3/PlaB samples. A module is not present or absent in one condition; every sample has a module score. M03/M06 are described as DP-dependent MUT-associated programs because genotype-by-DP is significant, while the genotype main effect is not.

## 5. MUT + inhibitor combination logic

Formal occlusion tests `[MUT inhibitor - MUT baseline] - [WT inhibitor - WT baseline]`. If MUT had already saturated exactly the same linear pathway, the incremental inhibitor response in MUT should be significantly attenuated.

{_markdown_table(occlusion)}

No M03/M06 condition shows formal occlusion. M03 under DP+PlaB shows significant enhancement instead of attenuation. The combination results therefore support shared output programs plus independent or context-dependent effects, not a proven single saturated pathway.

## 6. Core-to-module localization

The exact inhibitor-like cores are not significantly enriched in M03/M06 after multiple testing. Event-level overlap and all-sample module covariance answer different questions and must be shown as complementary evidence, not collapsed into one claim.

## 7. RBP prioritization

{_markdown_table(rbp)}

These correlations use all 36 samples. They are not MUT-specific by construction. `Module gene=yes` means an AS event from that RBP belongs to the module; it does not mean the RBP is expressed or spliced only in MUT. Use the label **module-associated RBP prioritization**, not **MUT-specific upstream regulators**.

## Claim boundaries

### Supported

- A finite Setd2 MUT AS subset converges with CLK3 and PlaB inhibition in exact event identity, direction and dPSI magnitude.
- M03/M06 are treatment-responsive co-splicing programs with significant genotype-by-DP effects.
- Several RBP genes are associated with module behavior and are suitable validation candidates.

### Not supported

- Complete inhibitor phenocopy.
- Formal proof that MUT is SETD2 LOF from RNA-seq alone.
- A single saturated linear pathway shared by MUT and inhibitor.
- Causal or MUT-specific RBP regulators.

## Recommended validation sequence

1. Confirm SETD2 biochemical/function loss: protein abundance, relevant chromatin mark/activity, or an assay tied to the mutation.
2. Rescue the mutation or compare with an isogenic SETD2 loss model.
3. Validate representative triple-core, CLK3-specific and PlaB-specific AS events by RT-PCR or targeted sequencing.
4. Test top module-associated RBPs only after stratified/partial-correlation and motif/CLIP convergence.
5. Use the combination experiment as an epistasis test, reporting interaction effect and FDR rather than visual similarity alone.
"""


def _build_storyboard(author: str) -> str:
    slides = [
        ("1", "Setd2 MUT shows partial inhibitor-like splicing convergence", "Minimal title slide", "Question: does MUT resemble functional splicing inhibition, and does this support LOF?"),
        ("2", "The RNA-seq result supports a partial LOF-compatible component, not definitive LOF", "Evidence pyramid", "State the answer first. Separate supported claims from unsupported claims."),
        ("3", "The analysis compares matched effects within DP and noDP contexts", "Study-design schematic", "Show MUT-vs-WT, WT+CLK3-vs-WT and WT+PlaB-vs-WT. Define exact event matching and thresholds."),
        ("4", "Exact AS overlap is far above chance in all four comparisons", "Four-row evidence table", "Display shared counts, containment, hypergeometric P/FDR. Do not collapse CLK3 and PlaB."),
        ("5", "Shared events agree in direction and dPSI magnitude", "01_exact_event_lof_evidence.png", "Emphasize 86.7-96.4% same direction and Spearman 0.72-0.92."),
        ("6", "MUT contains distinct CLK3-like, PlaB-like and triple event cores", "Core count matrix", "Define each core explicitly and show DP/noDP event and gene counts."),
        ("7", "Core gene enrichment is weak and does not establish one shared pathway", "GO/KEGG evidence table", "Only DP CLK3-like SSA repair terms pass FDR; show the four-gene caveat."),
        ("8", "M03 and M06 are all-sample programs, not MUT-only modules", "Module definition + compact table", "Explain that all 36 samples receive scores; treatment and DP dominate; genotype-by-DP is significant."),
        ("9", "Inhibitors move M03/M06 strongly, while the MUT difference depends on DP", "02_M03_M06_condition_scores.png", "Interpret within-module shifts only; module eigengene sign is arbitrary."),
        ("10", "Combination tests do not show systematic pathway occlusion", "03_M03_M06_occlusion_convergence.png", "Show formal interaction FDR. Flag M03 DP+PlaB enhancement rather than occlusion."),
        ("11", "The exact inhibitor-like core is not concentrated in M03/M06", "04_shared_core_module_mapping.png", "Event overlap and co-splicing covariance are complementary, not equivalent."),
        ("12", "RBP analysis prioritizes validation candidates but is not causal", "RBP evidence table", "Use all-sample correlations; distinguish module membership, expression correlation and motif evidence."),
        ("13", "The strongest conclusion rests on event identity, direction and effect agreement", "Claim/evidence matrix", "Rank direct event evidence above modules, GO and RBP prioritization."),
        ("14", "Orthogonal validation is required to classify the mutation as LOF", "Validation sequence", "Biochemical loss, rescue/isogenic loss, event RT-PCR, then candidate RBP perturbation."),
        ("15", "Every number and figure is traceable to a packaged source table", "Reproducibility map", "Show package directories, parameters, manifest, checksums and regeneration command."),
    ]
    lines = [
        "# PowerPoint-ready storyboard",
        "",
        f"Prepared for: {author}",
        "",
        "Audience outcome: the lab should understand exactly what the RNA-seq supports about LOF, what it does not prove, and which validation experiments are decision-critical.",
        "",
    ]
    for number, title, visual, notes in slides:
        lines.extend(
            [
                f"## Slide {number} — {title}",
                "",
                f"Primary visual: {visual}",
                "",
                f"Speaker intent: {notes}",
                "",
            ]
        )
    return "\n".join(lines)


def _build_speaker_notes() -> str:
    return """# 中文讲稿建议

## 开场

我们要回答的不是抑制剂是否造成很多剪接变化，而是 Setd2 MUT 的剪接表型是否包含与功能抑制一致的分子成分，以及组合实验是否支持同一路径关系。

## 一句话答案

MUT 与 CLK3/PlaB 抑制存在统计上极显著、方向高度一致的共同 AS 核心，因此与部分 LOF-like response 相容；但 MUT 仍保留明显抑制剂响应，module 层面没有系统性 occlusion，所以不能仅凭 RNA-seq 宣称完整 LOF 或完整 phenocopy。

## 事件层

先强调匹配的是同一个 rMATS event，而不是同一个 gene。四组重合检验 FDR 均远小于 1e-39。共同事件中 86.7-96.4% 的 dPSI 方向一致，效应量相关性也很高。这个层级是最直接的证据。

## 功能层

CLK3-like、PlaB-like 和三方核心必须分开。当前 GO/KEGG 的阳性结果很少，不能把事件重合自动解释成统一剪接通路。下一步更适合做 event-level coding consequence、RBP motif/CLIP 和定向实验验证。

## Module 层

M03/M06 是用全部样本建立的共剪接程序，不是 MUT-only module。它们的 treatment effect 很强，DP effect 和 genotype-by-DP 也显著。Module 证据说明存在协同变化的剪接程序，但不能替代具体共同事件。

## 组合实验

如果 MUT 与抑制剂作用于完全相同且已饱和的线性路径，抑制剂在 MUT 中的附加效应应显著减弱。当前没有系统性 occlusion，M03 的 DP+PlaB 甚至增强。因此更符合部分共享输出加独立或上下文依赖作用。

## RBP

PTBP1、SRSF4、MBNL1/2、SRSF7 是 all-sample module-associated candidates，不是 MUT-specific regulators。相关性、module membership 和 motif 只能用于排验证优先级，不能写成上游因果机制。

## 收尾

汇报最后落到三个决定：优先验证哪些共同事件；用什么正交实验确认 SETD2 功能损失；组合实验是否需要更多重复或更直接的 epistasis readout。
"""


def _build_html_report(
    event: pd.DataFrame,
    cores: pd.DataFrame,
    functional: pd.DataFrame,
    modules: pd.DataFrame,
    occlusion: pd.DataFrame,
    rbp: pd.DataFrame,
    figure_dir: Path,
    author: str,
) -> str:
    figures = []
    for name in [
        "01_exact_event_lof_evidence.png",
        "02_M03_M06_condition_scores.png",
        "03_M03_M06_occlusion_convergence.png",
        "04_shared_core_module_mapping.png",
    ]:
        path = figure_dir / name
        if path.exists():
            figures.append((name, _data_uri(path)))
    figure_html = "".join(
        f'<figure><img src="{uri}" alt="{html.escape(name)}"><figcaption>{html.escape(name)}</figcaption></figure>'
        for name, uri in figures
    )
    return f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Setd2 MUT LOF-compatible splicing evidence</title>
<style>
:root {{ --ink:#111827; --mut:#2f67d8; --clk:#e88910; --pla:#168b75; --panel:#f3f5f7; --rule:#c8ced8; }}
* {{ box-sizing:border-box; }}
body {{ margin:0; color:var(--ink); font-family:Arial,Helvetica,sans-serif; line-height:1.5; background:white; }}
main {{ max-width:1180px; margin:0 auto; padding:52px 62px 80px; }}
h1 {{ font-size:44px; line-height:1.05; margin:0 0 14px; }}
h2 {{ font-size:28px; margin:54px 0 16px; padding-top:12px; border-top:2px solid var(--ink); }}
h3 {{ font-size:20px; margin:28px 0 10px; }}
.subtitle {{ color:#4b5563; font-size:18px; }}
.answer {{ margin:32px 0; padding:24px 28px; background:var(--panel); border-left:8px solid var(--mut); font-size:22px; font-weight:700; }}
.pyramid {{ display:grid; grid-template-columns:repeat(3,1fr); gap:14px; margin:22px 0; }}
.pyramid div {{ padding:18px; border:1px solid var(--rule); min-height:130px; }}
.pyramid b {{ display:block; font-size:18px; margin-bottom:8px; }}
.data-table {{ width:100%; border-collapse:collapse; font-size:13px; margin:16px 0 26px; }}
.data-table th {{ text-align:left; background:#111827; color:white; padding:8px; }}
.data-table td {{ border-bottom:1px solid #d7dce4; padding:7px 8px; vertical-align:top; }}
figure {{ margin:36px 0 54px; }}
figure img {{ display:block; width:100%; height:auto; border:1px solid #e0e4ea; }}
figcaption {{ color:#596273; font-size:12px; margin-top:8px; }}
.limit {{ color:#7f1d1d; background:#fff1f2; padding:18px 22px; }}
.supported {{ color:#064e3b; background:#ecfdf5; padding:18px 22px; }}
code {{ background:#eef0f3; padding:1px 4px; }}
@media print {{ main {{ max-width:none; padding:24px; }} h2 {{ page-break-before:always; }} figure {{ page-break-inside:avoid; }} }}
@media (max-width:800px) {{ main {{ padding:28px 20px; }} .pyramid {{ grid-template-columns:1fr; }} .data-table {{ font-size:11px; }} }}
</style>
</head>
<body><main>
<h1>Setd2 MUT LOF-compatible splicing evidence</h1>
<p class="subtitle">Prepared by {html.escape(author)} · Generated {datetime.now().astimezone().isoformat(timespec='seconds')}</p>
<div class="answer">Setd2 MUT contains a statistically non-random, directionally concordant subset of the CLK3- and PlaB-inhibitor splicing responses. This supports partial functional convergence compatible with a LOF component, but not definitive LOF or complete inhibitor phenocopy.</div>
<section class="pyramid">
<div><b>Strongest evidence</b>Exact shared significant AS events, overlap enrichment, direction and dPSI agreement.</div>
<div><b>Supporting evidence</b>All-sample M03/M06 co-splicing programs and genotype-by-DP associations.</div>
<div><b>Hypothesis generation</b>Core GO/KEGG, module localization and module-associated RBP prioritization.</div>
</section>
<h2>1 · Exact event overlap is far above chance</h2>
<p>Events are exact canonical rMATS matches. Significance requires FDR ≤ 0.05 and |standardized dPSI| ≥ 0.10.</p>
{_html_table(event)}
{figure_html.split('</figure>')[0]}</figure>
<h2>2 · Distinct CLK3-like, PlaB-like and triple cores</h2>
{_html_table(cores)}
<p>These sets are defined separately in DP and noDP. They must not be collapsed into one inhibitor category.</p>
<h2>3 · Functional enrichment is limited</h2>
{_html_table(functional)}
<p class="limit"><b>Boundary:</b> only the DP CLK3-like core has FDR-significant GO terms, driven by four genes. PlaB-like and triple cores do not reveal a robust common pathway.</p>
<h2>4 · M03/M06 are all-sample co-splicing programs</h2>
{_html_table(modules)}
{figures[1][1] if False else ''}
<figure><img src="{figures[1][1] if len(figures)>1 else ''}" alt="M03 M06 condition scores"><figcaption>Within-module condition shifts; module score sign is orientation-dependent.</figcaption></figure>
<h2>5 · Combination tests do not establish pathway saturation</h2>
{_html_table(occlusion)}
<figure><img src="{figures[2][1] if len(figures)>2 else ''}" alt="Occlusion and convergence"><figcaption>Formal inference uses interaction FDR; visual similarity is descriptive.</figcaption></figure>
<h2>6 · Shared cores are not concentrated in M03/M06</h2>
<figure><img src="{figures[3][1] if len(figures)>3 else ''}" alt="Core module mapping"><figcaption>Event overlap and network covariance are complementary evidence layers.</figcaption></figure>
<h2>7 · RBP results prioritize validation candidates</h2>
{_html_table(rbp)}
<p class="limit">These are all-sample module associations, not MUT-specific or causal regulators.</p>
<h2>8 · Claim boundaries and validation</h2>
<div class="supported"><b>Supported:</b> partial event-level convergence with CLK3/PlaB inhibition; strong direction and dPSI agreement; DP-dependent MUT-associated module behavior.</div>
<div class="limit"><b>Not supported:</b> complete inhibitor phenocopy, formal pathway occlusion, definitive SETD2 LOF, or causal MUT-specific RBPs.</div>
<ol><li>Biochemical SETD2 function/mark measurement.</li><li>Rescue or isogenic SETD2 loss comparison.</li><li>RT-PCR validation of triple, CLK3-like and PlaB-like events.</li><li>Only then prioritize RBP perturbation using multi-layer evidence.</li></ol>
</main></body></html>"""


def _write_manifest(output: Path, source_map: dict[Path, Path]) -> None:
    rows = []
    for path in sorted(p for p in output.rglob("*") if p.is_file()):
        if path.name in {"manifest.tsv", "SHA256SUMS.txt"}:
            continue
        relative = path.relative_to(output)
        source = source_map.get(path.resolve())
        rows.append(
            {
                "package_path": relative.as_posix(),
                "role": relative.parts[0] if len(relative.parts) > 1 else "report",
                "source_path": str(source) if source else "generated",
                "bytes": path.stat().st_size,
                "modified_utc": datetime.fromtimestamp(
                    path.stat().st_mtime, tz=timezone.utc
                ).isoformat(),
                "sha256": _sha256(path),
            }
        )
    manifest = output / "manifest.tsv"
    with manifest.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    checksum_lines = [f"{row['sha256']}  {row['package_path']}" for row in rows]
    checksum_lines.append(f"{_sha256(manifest)}  manifest.tsv")
    (output / "SHA256SUMS.txt").write_text("\n".join(checksum_lines) + "\n", encoding="utf-8")


def build(source: Path, output: Path, repo: Path, author: str) -> Path:
    if not source.exists():
        raise FileNotFoundError(f"Closure output not found: {source}")
    required = [
        "exact_shared_event_lof_summary.tsv",
        "exact_core_gene_membership.tsv",
        "exact_core_event_type_summary.tsv",
        "exact_core_go_kegg_enrichment.tsv",
        "exact_module_annotation_summary.tsv",
        "exact_module_occlusion_convergence.tsv",
        "exact_module_rbp_regulators.tsv",
        "parameters.json",
    ]
    missing = [name for name in required if not (source / name).exists()]
    if missing:
        raise FileNotFoundError(f"Missing required closure outputs: {', '.join(missing)}")

    output.mkdir(parents=True, exist_ok=True)
    data_dir = output / "data" / "closure"
    figure_dir = output / "figures"
    provenance_dir = output / "provenance"
    source_dir = output / "source_code"
    for directory in [data_dir, figure_dir, provenance_dir, source_dir]:
        directory.mkdir(parents=True, exist_ok=True)

    source_map: dict[Path, Path] = {}
    for original, copied in _copy_tree_files(source, data_dir):
        source_map[copied.resolve()] = original.resolve()
    for plot in sorted((source / "plots").glob("*.png")):
        target = figure_dir / plot.name
        shutil.copy2(plot, target)
        source_map[target.resolve()] = plot.resolve()

    for code_path in [
        repo / "scripts" / "analyze_setd2_lof_hypothesis.py",
        repo / "katana_cosplice" / "functional_enrichment.py",
        Path(__file__).resolve(),
    ]:
        if code_path.exists():
            target = source_dir / code_path.name
            shutil.copy2(code_path, target)
            source_map[target.resolve()] = code_path.resolve()

    event = _event_summary(source)
    cores = _core_summary(source)
    functional = _functional_summary(source)
    modules = _module_summary(source)
    occlusion = _occlusion_summary(source)
    rbp = _rbp_summary(source)

    reports = {
        "SETD2_LOF_ANALYSIS_REPORT.md": _build_report_markdown(
            event, cores, functional, modules, occlusion, rbp, author
        ),
        "PRESENTATION_STORYBOARD_EN.md": _build_storyboard(author),
        "SPEAKER_NOTES_CN.md": _build_speaker_notes(),
        "SETD2_LOF_EVIDENCE_REPORT.html": _build_html_report(
            event, cores, functional, modules, occlusion, rbp, figure_dir, author
        ),
    }
    for name, content in reports.items():
        (output / name).write_text(content, encoding="utf-8")

    slide_sources = pd.DataFrame(
        [
            (4, "Exact overlap significance", "exact_shared_event_lof_summary.tsv"),
            (5, "Containment, direction and dPSI agreement", "exact_shared_event_lof_summary.tsv; 01_exact_event_lof_evidence.png"),
            (6, "CLK3/PlaB/triple core counts", "exact_core_gene_membership.tsv; exact_core_event_type_summary.tsv"),
            (7, "Core GO/KEGG", "exact_core_go_kegg_enrichment.tsv; exact_core_go_kegg_provider.json"),
            (8, "M03/M06 identity", "exact_module_annotation_summary.tsv; exact_module_factor_contrasts.tsv"),
            (9, "Module condition scores", "exact_module_condition_scores.tsv; 02_M03_M06_condition_scores.png"),
            (10, "Occlusion/convergence", "exact_module_occlusion_convergence.tsv; 03_M03_M06_occlusion_convergence.png"),
            (11, "Core-to-module mapping", "exact_core_module_enrichment.tsv; 04_shared_core_module_mapping.png"),
            (12, "RBP prioritization", "exact_module_rbp_regulators.tsv; exact_module_rbp_motif_enrichment.tsv"),
        ],
        columns=["slide", "claim", "source_files"],
    )
    slide_sources.to_csv(output / "FIGURE_AND_SLIDE_SOURCE_MAP.tsv", sep="\t", index=False)

    provenance = {
        "generated_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "author": author,
        "source_closure": str(source.resolve()),
        "report_package": str(output.resolve()),
        "python": sys.version,
        "platform": platform.platform(),
        "git_commit": _git(repo, "rev-parse", "HEAD"),
        "git_branch": _git(repo, "branch", "--show-current"),
        "git_status": _git(repo, "status", "--short"),
        "regeneration_command": (
            f'python "{Path(__file__).resolve()}" --source "{source.resolve()}" '
            f'--output "{output.resolve()}" --author "{author}"'
        ),
    }
    (provenance_dir / "provenance.json").write_text(
        json.dumps(provenance, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    (provenance_dir / "git_status.txt").write_text(provenance["git_status"] + "\n", encoding="utf-8")
    (provenance_dir / "python_environment.txt").write_text(
        subprocess.run(
            [sys.executable, "-m", "pip", "freeze"],
            capture_output=True,
            text=True,
            check=False,
        ).stdout,
        encoding="utf-8",
    )
    _write_upstream_input_manifest(
        source, provenance_dir / "UPSTREAM_INPUT_MANIFEST.tsv"
    )

    readme = f"""# Setd2 LOF-compatible splicing reporting package

This package contains the report-ready interpretation of the Katana/CoSplice closure analysis.

## Start here

- `SETD2_LOF_EVIDENCE_REPORT.html`: self-contained visual report.
- `SETD2_LOF_ANALYSIS_REPORT.md`: article-style scientific narrative.
- `PRESENTATION_STORYBOARD_EN.md`: 15-slide PowerPoint-ready structure.
- `SPEAKER_NOTES_CN.md`: Chinese presentation notes.
- `FIGURE_AND_SLIDE_SOURCE_MAP.tsv`: maps claims and slides to source files.

## Traceability

- `data/closure/`: complete closure outputs, including all TSV/JSON/plots and the original interpretation report.
- `figures/`: presentation-ready copies of generated figures.
- `source_code/`: analysis and packaging scripts used for this release.
- `provenance/`: environment, Git state, source/output paths and regeneration command.
- `provenance/UPSTREAM_INPUT_MANIFEST.tsv`: every analysis input path, file size, timestamp and SHA-256.
- `manifest.tsv`: file inventory, source path, size, timestamp and SHA-256.
- `SHA256SUMS.txt`: package integrity checksums.

Raw BAM/rMATS inputs are not duplicated. Every rMATS table actually consumed by this closure, together with the stable PSI matrix, sample/event metadata, module tables and RBP inputs, is recorded and hashed in `provenance/UPSTREAM_INPUT_MANIFEST.tsv`.

## Regenerate

```powershell
{provenance['regeneration_command']}
```

## Scientific boundary

The package supports partial inhibitor-like molecular convergence. It does not claim complete phenocopy, formal pathway saturation, definitive SETD2 LOF, or causal MUT-specific RBP regulators.
"""
    (output / "README.md").write_text(readme, encoding="utf-8")

    _write_manifest(output, source_map)

    archive = output.with_suffix(".zip")
    with zipfile.ZipFile(archive, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=6) as zipped:
        for path in sorted(p for p in output.rglob("*") if p.is_file()):
            archive_path = Path(output.name) / path.relative_to(output)
            zipped.write(path, arcname=archive_path.as_posix())
    return archive


def main() -> None:
    parser = argparse.ArgumentParser(description="Build a traceable Setd2 LOF reporting package.")
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--repo", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--author", default="Tara")
    args = parser.parse_args()
    archive = build(args.source, args.output, args.repo, args.author)
    print(f"Report package: {args.output.resolve()}")
    print(f"Archive: {archive.resolve()}")


if __name__ == "__main__":
    main()
