"""Readable, reload-verified SE motif evidence, RNA maps and sequence logos."""
from __future__ import annotations

import hashlib
import html
import json
import re
import shutil
from pathlib import Path
from zipfile import ZIP_DEFLATED, ZipFile

import numpy as np
import pandas as pd

from katana_cosplice.motif_scanner import IndexedFasta, load_motif_database
from katana_next.science.se_rbp import sequence_from_blocks


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda: f.read(8*1024*1024), b""):
            h.update(block)
    return h.hexdigest()


def verify_sites(directory):
    directory = Path(directory)
    frozen = json.loads((directory/"ANALYSIS_FREEZE.json").read_text(encoding="utf-8"))
    fasta_path = frozen.get("genome_fasta") or next(x["path"] for x in frozen["input_files"] if x["path"].endswith(".fa"))
    fasta = IndexedFasta(Path(fasta_path))
    sites = pd.read_parquet(directory/"sites.parquet")
    windows = pd.read_parquet(directory/"sequences.parquet").set_index("window_id")
    unique = sites.drop_duplicates(["scan_token", "window_id", "start_1based", "stop_1based"])
    cache = {}
    with fasta.fasta.open("rb") as stream:
        for row in unique.itertuples(index=False):
            w = windows.loc[row.window_id]
            a, b = int(row.start_1based)-1, int(row.stop_1based)
            if row.matched_rna_sequence != w.sequence[a:b].replace("T", "U"):
                raise ValueError("Matched sequence differs from original scanned string")
            key = (row.chr, row.strand, row.hit_genomic_blocks_json)
            if key not in cache:
                cache[key] = sequence_from_blocks(fasta, stream, row.chr, row.strand,
                                                 json.loads(row.hit_genomic_blocks_json)).replace("T", "U")
            if cache[key] != row.matched_rna_sequence:
                raise ValueError("Actual FASTA block sequence round-trip mismatch")
            if "junction" in row.region and not a < row.junction_offset < b:
                raise ValueError("Noncrossing site incorrectly marked as splice-junction motif")
            if len(row.matched_rna_sequence) != row.width or row.site_pvalue > frozen["policy"]["site_pvalue"]:
                raise ValueError("FIMO width/threshold mismatch")
    result = {"all_unique_sites_checked": len(unique), "genome_substrings_reloaded": len(cache),
              "negative_strand_sites_checked": int(unique.strand.eq("-").sum()),
              "junction_sites_checked": int(unique.region.str.contains("junction").sum()),
              "exact_sequence_and_block_round_trip": True,
              "original_input_hashes_unchanged": all(digest(x["path"]) == x["sha256"] for x in frozen["input_files"])}
    if not result["original_input_hashes_unchanged"]:
        raise ValueError("Original input changed")
    (directory/"SEQUENCE_VALIDATION.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    return result


def _logo(pwm, title, path):
    import matplotlib.pyplot as plt
    from matplotlib.font_manager import FontProperties
    from matplotlib.patches import PathPatch
    from matplotlib.textpath import TextPath
    from matplotlib.transforms import Affine2D
    pwm = np.asarray(pwm, float)
    pwm = pwm/pwm.sum(axis=1, keepdims=True)
    entropy = -(pwm*np.log2(np.maximum(pwm, 1e-300))).sum(axis=1)
    heights = pwm*(2-entropy)[:, None]
    colors = dict(A="#228833", C="#4477AA", G="#CCBB44", U="#CC6677")
    fig, ax = plt.subplots(figsize=(max(4, len(pwm)*.38), 2.4))
    for i, column in enumerate(heights):
        bottom = 0
        for j in np.argsort(column):
            height, base = column[j], "ACGU"[j]
            if height < 1e-7:
                continue
            text = TextPath((0, 0), base, size=1, prop=FontProperties(family="DejaVu Sans", weight="bold"))
            bounds = text.get_extents()
            transform = Affine2D().translate(-bounds.x0, -bounds.y0).scale(.86/bounds.width, height/bounds.height).translate(i+.07, bottom)
            ax.add_patch(PathPatch(text, transform=transform+ax.transData, color=colors[base], lw=0))
            bottom += height
    ax.set(xlim=(0, len(pwm)), ylim=(0, 2.05), ylabel="Information (bits)", xlabel="Motif position")
    ax.set_xticks(np.arange(len(pwm))+.5, np.arange(1, len(pwm)+1))
    ax.set_title(title, fontsize=10)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    fig.savefig(path.with_suffix(".svg"))
    plt.close(fig)


def _map(sites, windows, pairs, rbp, evidence, output):
    import matplotlib.pyplot as plt
    own = sites[sites.RBP.str.casefold().eq(rbp.casefold()) & sites.evidence_class.eq(evidence)]
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.2), sharey=True)
    curves = []
    for ax, region, edges in zip(axes, ["upstream_intron", "cassette_exon", "downstream_intron"],
                                [np.arange(-200, -19, 20), np.linspace(0, 100, 11), np.linspace(6, 200, 11)], strict=True):
        valid = set(windows.loc[windows.region.eq(region), "event_id"])
        for direction, color in [("inclusion_up", "#0072B2"), ("inclusion_down", "#D55E00")]:
            paired = pairs[pairs.status.eq("matched") & pairs.direction.eq(direction)]
            for role, column, style in [("foreground", "foreground_event", "-"), ("background", "background_event", "--")]:
                ids = set(paired[column]) & valid
                hits = own[own.region.eq(region) & own.event_id.isin(ids)]
                values = []
                for a, b in zip(edges[:-1], edges[1:]):
                    h = hits[hits.relative_position.ge(a) & hits.relative_position.lt(b)]
                    fraction = h.event_id.nunique()/len(ids) if ids else np.nan
                    values.append(fraction)
                    curves.append({"RBP": rbp, "evidence_class": evidence, "region": region, "direction": direction,
                                   "role": role, "bin_start": a, "bin_end": b, "event_denominator": len(ids), "fraction_with_site": fraction})
                ax.plot((edges[:-1]+edges[1:])/2, values, style, color=color,
                        label=f"{direction.replace('inclusion_', '')} {role} (n={len(ids)})", lw=1.6)
        ax.set_title(region.replace("_", " "), fontsize=10)
        ax.set_xlabel("Cassette length (%)" if region == "cassette_exon" else "RNA-oriented distance from cassette boundary (nt)", fontsize=8)
        ax.spines[["top", "right"]].set_visible(False)
        ax.grid(alpha=.15)
    axes[0].set_ylabel("Fraction of matched gene events with site")
    axes[-1].legend(fontsize=7, loc="upper left", bbox_to_anchor=(1.01, 1))
    fig.suptitle(f"{rbp} | {evidence} | descriptive positional map, not per-bin significance", fontsize=10)
    fig.tight_layout()
    fig.savefig(output, dpi=160, bbox_inches="tight")
    fig.savefig(output.with_suffix(".svg"), bbox_inches="tight")
    plt.close(fig)
    return curves


def build_report(directory):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    directory = Path(directory)
    (directory/"figures").mkdir(exist_ok=True)
    (directory/"logos").mkdir(exist_ok=True)
    summary = json.loads((directory/"RESULTS.json").read_text(encoding="utf-8"))
    frozen = json.loads((directory/"ANALYSIS_FREEZE.json").read_text(encoding="utf-8"))
    enrichment = pd.read_parquet(directory/"enrichment.parquet")
    catalog = pd.read_parquet(directory/"motif_catalog.parquet")
    sites = pd.read_parquet(directory/"sites.parquet")
    windows = pd.read_parquet(directory/"sequences.parquet")
    pairs = pd.read_parquet(directory/"matching.parquet")
    fasta_text = "".join(f">{r.window_id} {r.region} {r.gene_symbol} {r.role} genomic_strand={r.genomic_strand}\n{r.sequence_rna}\n"
                         for r in windows.itertuples(index=False))
    fasta_export = directory/"all_sequences.rna.fa"
    if fasta_export.exists() and fasta_export.read_text(encoding="ascii") != fasta_text:
        old = directory/"all_sequences.original_scan_export.rna.fa"
        if not old.exists():
            shutil.copy2(fasta_export, old)
    fasta_export.write_text(fasta_text, encoding="ascii")
    # A necessary lower bound on ANY best-word P is min(background)^width.
    # This flags proven lack of screening power, without deleting hypotheses or
    # retrospectively relaxing the site threshold / BH family after seeing P.
    coverage = catalog.copy()
    coverage["site_P_hard_lower_bound_by_width"] = min(summary["fimo"]["background_ACGT"])**coverage.width
    coverage["provably_unattainable_at_fixed_site_P"] = coverage.site_P_hard_lower_bound_by_width.gt(summary["policy"]["site_pvalue"])
    coverage.to_csv(directory/"SITE_THRESHOLD_COVERAGE.tsv", sep="\t", index=False)
    impossible = int(coverage.provably_unattainable_at_fixed_site_P.sum())
    significant = enrichment[enrichment.enrichment_fdr.lt(.05)].copy()
    significant.to_csv(directory/"significant_enrichment_BH005.tsv", sep="\t", index=False)
    foreground = sites[sites.role.eq("foreground")].copy()
    # Enrichment is an RBP/evidence union test, not individual PWM significance.
    fg = foreground.merge(enrichment[["RBP", "evidence_class", "region", "direction", "enrichment_pvalue", "enrichment_fdr", "status"]],
                          on=["RBP", "evidence_class", "region", "direction"], how="left", validate="many_to_one")
    fg.to_csv(directory/"DAS_actual_motif_hits_with_enrichment.tsv", sep="\t", index=False)
    fg[fg.enrichment_fdr.lt(.05)].to_csv(directory/"significant_RBP_actual_hit_sequences.tsv", sep="\t", index=False)
    bed = []
    for i, row in enumerate(foreground.drop_duplicates(["scan_token", "window_id", "start_1based", "stop_1based"]).itertuples(index=False)):
        blocks = sorted(json.loads(row.hit_genomic_blocks_json))
        start, end = blocks[0][0], blocks[-1][1]
        bed.append([row.chr, start, end, f"site{i:07d}:{row.scan_token}:{row.window_id}", 0, row.strand, start, end, "0", len(blocks),
                    ",".join(str(b-a) for a, b in blocks)+",", ",".join(str(a-start) for a, b in blocks)+","])
    pd.DataFrame(bed).to_csv(directory/"DAS_motif_sites.split_blocks.bed", sep="\t", header=False, index=False)
    motifs = []
    for db in frozen["databases"]:
        motifs.extend(load_motif_database(Path(db["path"]), organism=frozen.get("species", "Mus_musculus")))
    by_identity = {(m.rbp, m.motif_id, m.source, m.evidence_class): m for m in motifs}
    motif_lookup = {row.motif_token: by_identity[(row.RBP, row.motif_id, row.source, row.evidence_class)]
                    for row in catalog.itertuples(index=False)}
    ranked = enrichment.sort_values(["enrichment_fdr", "enrichment_pvalue", "RBP", "region", "direction"])
    display = ranked.drop_duplicates(["RBP", "evidence_class"]).head(6)
    logo_tokens = set(t for s in significant.motif_tokens for t in s.split(";"))
    # Also illustrate a prespecified top-six descriptive view when no discovery
    # survives correction; it is never relabelled as significant enrichment.
    for row in display.itertuples(index=False):
        choices = catalog[(catalog.RBP.str.casefold() == row.RBP.casefold()) & catalog.evidence_class.eq(row.evidence_class)]
        observed = foreground[foreground.motif_token.isin(choices.motif_token)].groupby("motif_token").event_id.nunique()
        token = observed.sort_values(ascending=False).index[0] if len(observed) else choices.motif_token.iloc[0]
        logo_tokens.add(token)
    logo_manifest = []
    for token in sorted(logo_tokens):
        motif = motif_lookup[token]
        pwm = np.asarray(motif.pwm)/np.asarray(motif.pwm).sum(axis=1, keepdims=True)
        pd.DataFrame(pwm, columns=list("ACGU")).rename_axis("position_0based").to_csv(directory/"logos"/(token+".pwm.tsv"), sep="\t")
        _logo(pwm, f"{motif.rbp} | {token}\n{motif.evidence_class}", directory/"logos"/(token+".png"))
        logo_manifest.append({"motif_token": token, "RBP": motif.rbp, "motif_id": motif.motif_id,
                              "source": motif.source, "evidence_class": motif.evidence_class,
                              "consensus_not_actual_hit": "".join("ACGU"[x] for x in pwm.argmax(axis=1)), "logo": f"logos/{token}.png"})
    logos = pd.DataFrame(logo_manifest)
    logos.to_csv(directory/"logo_manifest.tsv", sep="\t", index=False)
    curves, map_files = [], []
    for i, row in enumerate(display.itertuples(index=False)):
        path = directory/"figures"/f"position_{i+1}_{re.sub('[^A-Za-z0-9]', '_', row.RBP)}.png"
        curves.extend(_map(sites, windows, pairs, row.RBP, row.evidence_class, path))
        map_files.append((row, path))
    pd.DataFrame(curves).to_csv(directory/"positional_map_source.tsv", sep="\t", index=False)
    top = ranked.head(20).iloc[::-1]
    fig, ax = plt.subplots(figsize=(11, 7))
    y = np.arange(len(top))
    colors = np.where(top.enrichment_fdr.lt(.05), "#0072B2", "#8C929A")
    ax.scatter(top.hit_fraction_difference, y, s=25+30*-np.log10(top.enrichment_pvalue.clip(lower=1e-12)), c=colors)
    labels = [f"{r.RBP} | {r.region} | {r.direction} | {r.evidence_class}" for r in top.itertuples()]
    ax.set_yticks(y, labels, fontsize=7)
    ax.axvline(0, color="black", lw=.7)
    ax.set_xlabel("Foreground minus matched-control hit fraction")
    ax.set_title("SE RBP hypotheses: top 20 ranked, blue only if global BH FDR < 0.05")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(directory/"figures/enrichment_overview.png", dpi=160)
    fig.savefig(directory/"figures/enrichment_overview.svg")
    plt.close(fig)
    validation = verify_sites(directory)
    columns = ["RBP", "evidence_class", "region", "direction", "matched_gene_pairs", "foreground_hit_genes", "background_hit_genes", "hit_fraction_difference", "enrichment_pvalue", "enrichment_fdr"]
    heading = (f"发现 {len(significant)} 个通过统一 BH FDR < 0.05 的 RBP/证据/区域/方向假设。" if len(significant) else
               "本次未发现通过统一 BH FDR < 0.05 的 motif 富集。不能据此宣称 DAS 已造成 RBP 功能改变，也不能证明不存在作用。")
    parts = ["<!doctype html><meta charset='utf-8'><title>Katana SE RBP real evidence</title>",
        "<style>body{font:16px Georgia,serif;max-width:1280px;margin:40px auto;color:#20313c;padding:0 25px}h1,h2{color:#12576b}table{border-collapse:collapse;font:12px sans-serif;max-width:100%}td,th{padding:6px;border:1px solid #cbd4d9}img{max-width:100%}code{overflow-wrap:anywhere}a{color:#006b8c}.note{background:#f0f5f7;padding:18px}.rna{font:13px monospace;overflow-wrap:anywhere}details{padding:8px;border-bottom:1px solid #ddd}</style>",
        f"<h1>Katana：SE 序列影响与 RBP motif</h1><p>{html.escape(summary['comparison'])}；{html.escape(frozen.get('species', 'Mus_musculus'))} / {html.escape(frozen['genome_build'])}；dPSI = {html.escape(frozen.get('numerator', 'numerator'))} − {html.escape(frozen.get('denominator', 'denominator'))} inclusion。</p>",
        f"<div class='note'><b>{heading}</b><p>真实已有 DAS 后处理；没有重算量化、DAS 模型或提交 HPC 作业。motif 命中是序列预测，不是本细胞中的实际结合，更不是因果功能验证。</p></div>",
        f"<p>保留 {summary['selected_SE']} 个 SE / {summary['selected_genes']} 个基因；独立匹配 {sum(summary['matched_by_direction'].values())} 对基因；方向分母 {html.escape(str(summary['matched_by_direction']))}。数据库 {summary['motif_assignments']} 条物种匹配 motif 指派 / {summary['RBPs']} 个 RBP，证据层分开；共 {summary['hypotheses']} 个假设。</p>",
        "<h2>1. 富集结论</h2><p>下表展示前 20 行，完整分母、P/FDR、阴性、匹配损耗均在附件。灰色排名不是显著发现。site P / site q 与富集 FDR 是不同检验。</p>",
        f"<p><b>检出能力限制：</b>固定 site P≤{summary['policy']['site_pvalue']:g} 对部分短 motif 过严；{impossible}/{len(catalog)} 条 motif 指派仅按宽度即可证明达不到阈值，必须记为筛查无检出能力，而不能解释为该 RBP 不存在作用。原 {summary['hypotheses']:,} 行检验家族没有因看到结果而删项或放宽阈值。详见 SITE_THRESHOLD_COVERAGE.tsv。</p>",
        ranked[columns].head(20).to_html(index=False, float_format=lambda x: f"{x:.4g}"),
        "<img src='figures/enrichment_overview.png'>",
        "<h2>2. 为什么这是 DAS 影响的序列</h2><p>扫描被包含/跳过的完整 cassette exon；另扫描两种 inclusion junction 和 exclusion junction，只有真正跨接头的命中才计为接头位点。用 dPSI 方向标注其在 RNA 群体中的预期表示量增加/减少；这不是 DNA motif 突变或结合亲和力改变。</p><p>上游/下游内含子另外分析，RNA 方向统一；屏蔽靠近剪接位点的 20/6 nt。它们只能提示调控关联，不能写成被 SE 直接删掉的序列。</p>",
        "<h2>3. 位置图与 motif logo</h2><p>位置曲线为描述性分箱比例，不进行逐箱显著性宣称。实线 DAS，虚线匹配背景；蓝色 inclusion 增加，橙色 inclusion 减少。logo 来自数据库原 PWM，不能代替下面的真实命中序列。</p>"]
    for row, path in map_files:
        parts.append(f"<h3>{html.escape(row.RBP)} / {html.escape(row.evidence_class)}：最佳行 FDR={row.enrichment_fdr:.4g}</h3><img src='{path.relative_to(directory).as_posix()}'>")
        available = logos[(logos.RBP.str.casefold() == row.RBP.casefold()) & logos.evidence_class.eq(row.evidence_class)]
        for logo in available.head(4).itertuples():
            parts.append(f"<img style='width:420px' src='{logo.logo}'><p>{html.escape(logo.motif_id)}；{html.escape(logo.source)}</p>")
    parts.append("<h2>4. 真实命中序列示例</h2><p>预先固定展示排序：富集 FDR、DAS FDR、event ID、motif token。示例不是额外独立验证；所有命中及阴性假设见完整附件。</p>")
    examples = fg.sort_values(["enrichment_fdr", "DAS_fdr", "event_id", "motif_token"]).drop_duplicates(["event_id", "region"]).head(12)
    for row in examples.itertuples():
        parts.append(f"<details open><summary>{html.escape(row.gene_symbol)} | {html.escape(row.RBP)} | {row.region} | dPSI={row.dpsi:+.3f}, DAS FDR={row.DAS_fdr:.4g}</summary><p>{row.chr} ({row.strand}); RNA-ordered 0-based half-open blocks: <code>{row.hit_genomic_blocks_json}</code></p><p class='rna'>{html.escape(row.rna_context)}</p><p>{row.availability_interpretation}; site P={row.site_pvalue:.3g}, enrichment FDR={row.enrichment_fdr:.4g}；{row.evidence_class}。共享相同 PWM 的 RBP：{html.escape(row.indistinguishable_PWM_RBPs)}</p></details>")
    parts.extend(["<h2>5. 数据库覆盖与边界</h2><p>本次实际扫描 CISBP-RNA 2.00 与 ATtRACT 中的可用 mouse PWM。RBPDB、RNAcompete、HTR-SELEX、RBNS 等与聚合库有重叠，按来源保留而不把数据库条目当独立生物学证据。CISBP 中的 RNA-Bind-n-Seq 指派已纳入，其 mouse 同源推断仍标为 inferred；没有把人源实验伪称为鼠直接实验。mCross / ENCORI 等独立资源仍有未接入缺口；oRNAment 是预测位点，CLIP 是参考细胞实验结合，不等于本样本功能。详见 16 项公开资源清单与 native source audit；不能宣称穷尽所有公开 motif 或所有 RBP。</p>",
        "<h2>6. 附件与可复核性</h2><ul>"])
    for path in ["enrichment.tsv", "significant_enrichment_BH005.tsv", "matching.tsv", "sequence_audit.tsv", "frozen_selected_DAS_SE.tsv", "DAS_actual_motif_hits_with_enrichment.tsv", "significant_RBP_actual_hit_sequences.tsv", "DAS_actual_motif_hits.rna.fa", "all_sequences.rna.fa", "DAS_motif_sites.split_blocks.bed", "logo_manifest.tsv", "motif_catalog.tsv", "positional_map_source.tsv", "ANALYSIS_FREEZE.json", "SEQUENCE_VALIDATION.json", "PUBLIC_DATABASE_COVERAGE.tsv", "AGGREGATED_SOURCE_AUDIT.tsv", "SITE_THRESHOLD_COVERAGE.tsv", "fimo/motifs.meme"]:
        parts.append(f"<li><a href='{path}'>{path}</a></li>")
    parts.extend(["</ul>", f"<p>所有 {validation['all_unique_sites_checked']} 个 unique PWM 位点均重新从 FASTA 分段取回并逐字符核对；其中负链 {validation['negative_strand_sites_checked']}，接头 {validation['junction_sites_checked']}。原输入全部 SHA 保持。</p>",
        "<p>富集背景控制 GC、cassette 长度、读深、denominator PSI；内含子长度/组成、细胞类型及复合蛋白作用仍可能混杂。近端 PWM 不能替代 RNA 结构、eCLIP 或 RBP 扰动验证。CI 为配对比例差的近似区间，小分母以 exact P 为准。</p>",
        "<h2>7. 其他剪接类型</h2><p>A3SS/A5SS 改变外显子延伸片段和接头；RI 改变保留的内含子；MXE 替换外显子及接头。都可以由同版本基因组和事件坐标恢复相应序列，并已有 rMAPS2 等先例。当前实现只接受 SE，不把其他类型混入本次结果；后续需逐类验证路径与坐标。</p>"])
    (directory/"REPORT.html").write_text("\n".join(parts), encoding="utf-8")
    (directory/"CONCLUSION.md").write_text(f"# {summary['comparison']} SE RBP\n\n{heading}\n\nSelected SE: {summary['selected_SE']}; matched directions: {summary['matched_by_direction']}; hypotheses: {summary['hypotheses']}.\n\nNo quantification or DAS refit. Motif association/predicted availability is not functional or causal validation.\n", encoding="utf-8")
    return {"significant_rows": len(significant), "logos": len(logos), "sequence_validation": validation}


def package_report(directory, archive=None):
    directory = Path(directory)
    archive = Path(archive) if archive else directory.parent/(directory.name+"_evidence.zip")
    files = sorted(p for p in directory.rglob("*") if p.is_file())
    manifest = [{"path": p.relative_to(directory).as_posix(), "bytes": p.stat().st_size, "sha256": digest(p)} for p in files if p.name != "EVIDENCE_MANIFEST.json"]
    (directory/"EVIDENCE_MANIFEST.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    with ZipFile(archive, "x", ZIP_DEFLATED) as z:
        for p in sorted(p for p in directory.rglob("*") if p.is_file()):
            z.write(p, p.relative_to(directory))
    with ZipFile(archive) as z:
        if z.testzip() is not None:
            raise ValueError("Evidence ZIP CRC failed")
        for item in manifest:
            if hashlib.sha256(z.read(item["path"])).hexdigest() != item["sha256"]:
                raise ValueError("Evidence ZIP member SHA failed")
    receipt = {"path": str(archive), "bytes": archive.stat().st_size, "sha256": digest(archive), "CRC_and_member_SHA_verified": True, "files": len(manifest)+1}
    archive.with_suffix(".verification.json").write_text(json.dumps(receipt, indent=2), encoding="utf-8")
    return receipt
