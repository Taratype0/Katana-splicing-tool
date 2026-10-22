"""Portable, evidence-separated report of a completed RNA-usage motif case."""
from __future__ import annotations

import html
import json
import re
import shutil
import time
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib.figure import Figure
from matplotlib.patches import Rectangle

from katana_next.science.motiftrace_usage import read_table, sha256, write_json


def _save(figure, output, stem):
    figure.savefig(output/(stem+".png"), dpi=160, bbox_inches="tight")
    figure.savefig(output/(stem+".svg"), bbox_inches="tight")


def _table(frame):
    return frame.to_html(index=False, escape=True, border=0, na_rep="NA / not computed", float_format=lambda x: f"{x:.5g}")


def _locus_bin_overlap(locus, bin_):
    if locus["chr"] != bin_["genomicData.seqnames"] or locus["strand"] != bin_["genomicData.strand"]:
        return 0
    return sum(max(0, min(b, bin_["end0"])-max(a, bin_["start0"])) for a, b in json.loads(locus["genomic_blocks_json"]))


def _select_example(loci, bins):
    candidates = loci[loci.model_kind.eq("consensus") & loci.pair_state.isin(["increased_form_only", "decreased_form_only"])].copy()
    passing = bins[bins.DEXSeq_bin_FDR.le(.05)].to_dict("records")
    def matched_bin(row):
        width = sum(b-a for a, b in json.loads(row["genomic_blocks_json"]))
        matches = [b for b in passing if _locus_bin_overlap(row, b) == width]
        return sorted(matches, key=lambda b: str(b["featureID"]))[0]["featureID"] if matches else ""
    candidates["contained_in_significant_bin"] = [matched_bin(row) for row in candidates.to_dict("records")]
    candidates = candidates[candidates.contained_in_significant_bin.ne("")]
    if candidates.empty:
        raise ValueError("No predeclared coordinate-supported consensus example; preserve this result rather than tune selection")
    candidates["abs_change"] = candidates.delta_fraction_availability.abs()
    return candidates.sort_values(["abs_change", "model_id", "site_key"], ascending=[False, True, True]).iloc[0].to_dict()


def _positions_figure(rna, pair_hits, case, site_key, selected_bin, output, symbol):
    txs = [case["decreased_transcript"], case["increased_transcript"]]
    forms = rna.set_index("transcript_id")
    figure = Figure(figsize=(12.5, 8), layout="constrained")
    axes = figure.subplots(3, 1)
    colors = ["#b86432", "#126a68"]
    for y, tx, color in zip([1, 0], txs, colors):
        row = forms.loc[tx]
        blocks = json.loads(row.blocks_json)
        axes[0].plot([min(a for a, b in blocks), max(b for a, b in blocks)], [y, y], color=color, lw=1)
        for a, b in blocks:
            axes[0].add_patch(Rectangle((a, y-.11), b-a, .22, color=color))
    axes[0].axvspan(selected_bin.start0, selected_bin.end0, color="#d2a149", alpha=.3,
                    label=f"DEXSeq {selected_bin.featureID}, q={selected_bin.DEXSeq_bin_FDR:.3g}")
    axes[0].set_yticks([1, 0], ["Decreased RNA", "Increased RNA"])
    axes[0].set_ylim(-.45, 1.45)
    axes[0].set_xlabel(f"{forms.iloc[0]['chr']} genomic coordinate (0-based); strand {forms.iloc[0].strand}")
    axes[0].legend(frameon=False, fontsize=9)
    axes[0].set_title(f"{symbol}: real annotated isoforms and the tested exon-counting bin", loc="left")
    for axis, tx, color in zip(axes[1:], txs, colors):
        row = forms.loc[tx]
        offset = 0
        for a, b in json.loads(row.blocks_json):
            axis.add_patch(Rectangle((offset, -.12), b-a, .24, facecolor="#e3e8e2", edgecolor="#687b78", linewidth=.5))
            offset += b-a
        for position in json.loads(row.junctions_json):
            axis.axvline(position, color="#aab6b1", lw=.6, ls=":")
        for hit in pair_hits[pair_hits.transcript_id.eq(tx)].to_dict("records"):
            selected = hit["site_key"] == site_key
            axis.add_patch(Rectangle((hit["RNA_start0"], .19 if selected else -.08), hit["RNA_end0"]-hit["RNA_start0"], .2 if selected else .16,
                                    color="#b9432e" if selected else color))
            if selected:
                axis.annotate(hit["hit_RNA"], (hit["RNA_start0"], .42), fontsize=9, color="#9e2e20")
        axis.set_xlim(-max(1, row.length_nt*.01), row.length_nt*1.02)
        axis.set_ylim(-.25, .65)
        axis.set_yticks([])
        axis.set_title(f"{tx}: {row.length_nt:,} nt; {len(pair_hits[pair_hits.transcript_id.eq(tx)])} literal hits of the selected model", fontsize=10, loc="left")
        axis.set_xlabel("Full real spliced reference RNA position, 5' to 3' (nt); dotted = exon joins")
    for axis in axes:
        axis.spines[["top", "right"]].set_visible(False)
    _save(figure, output, "REAL_RNA_FORMS_AND_MOTIF_POSITIONS")


def _sample_figure(usage, sample_loci, selected_key, sample_bins, selected_bin, output, symbol):
    samples = usage.drop_duplicates("sample_id").sort_values(["arm", "sample_id"]).sample_id.tolist()
    arms = usage.drop_duplicates("sample_id").set_index("sample_id").arm
    labels = [s.replace("EGFP_PBS_", "EGFP ").replace("CRE_PBS_", "CRE ") for s in samples]
    colors = ["#b86432" if arms[s] == "denominator" else "#126a68" for s in samples]
    x = np.arange(len(samples))
    figure = Figure(figsize=(13, 9), layout="constrained")
    axes = figure.subplots(2, 2)
    gene = usage.groupby("sample_id").gene_TPM.first().reindex(samples)
    axes[0, 0].bar(x, gene, color=colors)
    axes[0, 0].set_title(f"{symbol}: existing Salmon gene TPM", loc="left")
    axes[0, 0].set_ylabel("Sum of existing transcript TPM")
    fractions = usage.pivot(index="sample_id", columns="transcript_id", values="transcript_fraction").reindex(samples)
    bottom = np.zeros(len(samples))
    palette = ["#1c6f73", "#bd6c32", "#8ca9a1", "#b6c397", "#b58b9b", "#6c82a0", "#dac49a", "#8e7770", "#9ebbd0", "#a69e71", "#c3a5b3", "#9caaa1"]
    for index, tx in enumerate(fractions):
        values = fractions[tx].to_numpy()
        axes[0, 1].bar(x, values, bottom=bottom, color=palette[index % len(palette)], label=tx.replace("ENSMUST", ""))
        bottom += values
    axes[0, 1].set_ylim(0, 1)
    axes[0, 1].set_ylabel("Within-gene TPM fraction")
    axes[0, 1].set_title("All annotated quantified isoforms retained", loc="left")
    axes[0, 1].legend(fontsize=7, ncol=3, frameon=False, loc="upper center", bbox_to_anchor=(.5, -.22))
    site = sample_loci[sample_loci.site_key.eq(selected_key)].set_index("sample_id").reindex(samples)
    axes[1, 0].bar(x, site.RNA_fraction_with_locus, color=colors)
    axes[1, 0].set_ylim(0, 1)
    axes[1, 0].set_ylabel("Fraction of gene RNA carrying the locus")
    axes[1, 0].set_title("Selected real motif locus: usage-weighted availability", loc="left")
    bin_samples = sample_bins[sample_bins.featureID.eq(selected_bin.featureID)].set_index("sample_id").reindex(samples)
    axes[1, 1].bar(x, bin_samples.within_gene_bin_count_fraction, color=colors)
    axes[1, 1].set_ylabel("Fraction of this gene's exon-bin counts")
    axes[1, 1].set_title(f"DEXSeq {selected_bin.featureID}: observed counts, not PSI", loc="left")
    for axis in axes.flat:
        axis.set_xticks(x, labels, rotation=35, ha="right")
        axis.spines[["top", "right"]].set_visible(False)
    figure.suptitle("Existing biological samples only; no re-quantification, no newly inferred motif P/FDR", fontsize=12)
    _save(figure, output, "PER_SAMPLE_USAGE_EXPRESSION_AND_MOTIF_AVAILABILITY")


def _literal_logo(pattern, output):
    tokens = re.findall(r"\[[ACGT]+\]|[ACGT]", pattern)
    figure = Figure(figsize=(max(5, len(tokens)*.6), 2.7), layout="constrained")
    axis = figure.subplots()
    colors = {"A": "#417d46", "C": "#397aa2", "G": "#c8952c", "T": "#b64d43"}
    for i, token in enumerate(tokens):
        letters = token.strip("[]")
        bottom = 0.
        for base in letters:
            height = 1/len(letters)
            axis.bar(i, height, bottom=bottom, color=colors[base], width=.8)
            axis.text(i, bottom+height/2, "U" if base == "T" else base, ha="center", va="center", color="white", fontsize=17, weight="bold")
            bottom += height
    axis.set_xticks(range(len(tokens)), np.arange(1, len(tokens)+1))
    axis.set_ylim(0, 1)
    axis.set_ylabel("Allowed-base proportion")
    axis.set_title("Native consensus display, not an experimentally estimated PWM", fontsize=10)
    axis.spines[["top", "right"]].set_visible(False)
    _save(figure, output, "NATIVE_CONSENSUS_MOTIF_DISPLAY")


def build_report(prepared, output):
    prepared, output = Path(prepared), Path(output)
    if (output/"REPORT_MANIFEST.json").exists():
        history = output/"report_history"/str(time.time_ns())
        history.mkdir(parents=True)
        for p in output.iterdir():
            if p.is_file() and (p.suffix in {".png", ".svg", ".html", ".md"} or p.name in {"REPORT_MANIFEST.json", "ACTUAL_RESULT_REVIEW.json"}):
                shutil.copyfile(p, history/p.name)
    complete = json.loads((output/"USAGE_MOTIF_COMPLETE.json").read_text())
    plan = json.loads((prepared/"PLAN_BEFORE_CASE_SELECTION_AND_SCANNING.json").read_text())
    case = complete["case"]
    usage = pd.read_parquet(prepared/"EXISTING_GENE_TRANSCRIPT_USAGE.parquet")
    rna = pd.read_parquet(prepared/"REFERENCE_TRANSCRIPT_RNA.parquet")
    loci = pd.read_csv(output/"ALL_LOCUS_GAIN_LOSS_RETAINED_AND_COORDINATES.tsv.gz", sep="\t")
    bins = pd.read_csv(output/"ORIGINAL_DEXSEQ_BINS_P_FDR_AND_COUNTS.tsv", sep="\t")
    original_dtu = pd.read_csv(output/"ORIGINAL_STAGER_GENE_AND_TRANSCRIPT_Q.tsv", sep="\t")
    original_deg = pd.read_csv(output/"ORIGINAL_DEG_P_FDR_AND_EFFECT.tsv", sep="\t")
    original_se = pd.read_csv(output/"ORIGINAL_SE_P_FDR_AND_PSI.tsv", sep="\t")
    symbol = str(usage.gene_symbol.iloc[0]) if "gene_symbol" in usage else case["gene_id"]
    example = _select_example(loci, bins)
    selected_bin = bins[bins.featureID.eq(example["contained_in_significant_bin"])].iloc[0]
    model = example["model_id"]
    parts = sorted((output/"model_chunks").glob("*/COMPLETE.json"))
    selected_folder = next(p.parent for p in parts if model in json.loads(p.read_text())["models"])
    hits = pd.read_parquet(selected_folder/"hits.parquet")
    hits = hits[hits.model_id.eq(model)]
    sample_loci = pd.read_parquet(selected_folder/"locus_samples.parquet")
    selected_samples = sample_loci[sample_loci.site_key.eq(example["site_key"])].copy()
    selected_hits = hits[hits.site_key.eq(example["site_key"])].copy()
    sample_bins = pd.read_csv(output/"DEXSEQ_PER_SAMPLE_COUNTS_AND_WITHIN_GENE_FRACTIONS.tsv", sep="\t")
    _positions_figure(rna, hits, case, example["site_key"], selected_bin, output, symbol)
    _sample_figure(usage, sample_loci, example["site_key"], sample_bins, selected_bin, output, symbol)
    catalog = pd.read_csv(output/"FULL_MODEL_CATALOG_AND_SOURCE_EVIDENCE.tsv", sep="\t")
    model_info = catalog[catalog.model_id.eq(model)].iloc[0].to_dict()
    pattern_file = next(r for r in plan["source_inputs"] if Path(r["path"]).name == "NATIVE_CONSENSUS_PATTERNS.json")
    pattern = json.loads(Path(pattern_file["path"]).read_text())[model]
    _literal_logo(pattern, output)
    example.update(native_pattern=pattern, model_source=model_info)
    selected_samples.to_csv(output/"CASE_REAL_MOTIF_LOCUS_PER_SAMPLE_AVAILABILITY.tsv", sep="\t", index=False)
    selected_hits.to_csv(output/"CASE_REAL_MOTIF_HITS_AND_EXACT_POSITIONS.tsv", sep="\t", index=False)
    hits.to_csv(output/"CASE_MODEL_ALL_REAL_HITS_ON_ALL_GENE_ISOFORMS.tsv", sep="\t", index=False)
    evidence_rows = []
    for row in original_se.to_dict("records"):
        blocks = json.loads(example["genomic_blocks_json"])
        n = sum(max(0, min(b, row["exonEnd"])-max(a, row["exonStart_0base"])) for a, b in blocks) if row["chr"] == example["chr"] and row["strand"] == example["strand"] else 0
        evidence_rows.append({"SE_ID": row["ID"], "motif_overlap_target_exon_nt": n,
                              "DAS_P": row["PValue"], "DAS_FDR": row["FDR"], "DAS_dPSI": row["IncLevelDifference"],
                              "DEXSeq_featureID": selected_bin.featureID, "DEXSeq_P": selected_bin.DEXSeq_P, "DEXSeq_bin_FDR": selected_bin.DEXSeq_bin_FDR,
                              "not_a_combined_P_or_independent_replication": True})
    evidence = pd.DataFrame(evidence_rows)
    evidence.to_csv(output/"CASE_MOTIF_DAS_DEXSEQ_COORDINATE_EVIDENCE.tsv", sep="\t", index=False)
    psi_rows = []
    supported = evidence[evidence.motif_overlap_target_exon_nt.gt(0)].sort_values(["motif_overlap_target_exon_nt", "SE_ID"], ascending=[False, True])
    if len(supported):
        source_event = original_se[original_se.ID.eq(supported.iloc[0].SE_ID)].iloc[0]
        for arm, suffix in (("numerator", "1"), ("denominator", "2")):
            ids = plan[arm]["sample_ids"]
            values = str(source_event["IncLevel"+suffix]).split(",")
            ijc = str(source_event["IJC_SAMPLE_"+suffix]).split(",")
            sjc = str(source_event["SJC_SAMPLE_"+suffix]).split(",")
            if len(values) != len(ids) or len(ijc) != len(ids) or len(sjc) != len(ids):
                raise ValueError("Original rMATS replicate order disagrees with declared sample IDs")
            for sample, value, inc, exc in zip(ids, values, ijc, sjc):
                psi_rows.append({"sample_id": sample, "arm": arm, "SE_ID": source_event.ID,
                                 "original_PSI": pd.to_numeric(value, errors="coerce"),
                                 "original_inclusion_junction_count": pd.to_numeric(inc, errors="coerce"),
                                 "original_skipping_junction_count": pd.to_numeric(exc, errors="coerce"),
                                 "DAS_P": source_event.PValue, "DAS_FDR": source_event.FDR})
        psi = pd.DataFrame(psi_rows)
        psi.to_csv(output/"CASE_DAS_ORIGINAL_PSI_AND_COUNTS_PER_REAL_SAMPLE.tsv", sep="\t", index=False)
        ordered = selected_samples.sort_values(["arm", "sample_id"])
        sample_ids = ordered.sample_id.tolist()
        labels = [s.replace("EGFP_PBS_", "EGFP ").replace("CRE_PBS_", "CRE ") for s in sample_ids]
        colors = ["#b86432" if a == "denominator" else "#126a68" for a in ordered.arm]
        figure = Figure(figsize=(12, 4.7), layout="constrained")
        axes = figure.subplots(1, 2)
        axes[0].bar(range(len(sample_ids)), psi.set_index("sample_id").reindex(sample_ids).original_PSI, color=colors)
        axes[0].set_ylim(0, 1)
        axes[0].set_ylabel("Original rMATS junction PSI")
        axes[0].set_title(f"SE {int(source_event.ID)}: dPSI={source_event.IncLevelDifference:+.3f}, FDR={source_event.FDR:.3g}", loc="left", fontsize=10)
        axes[1].bar(range(len(sample_ids)), ordered.transcript_TPM_with_locus, color=colors)
        axes[1].set_ylabel("Sum of existing TPM for isoforms carrying this locus")
        axes[1].set_title("Same motif locus: abundance-weighted sequence opportunity", loc="left", fontsize=10)
        for axis in axes:
            axis.set_xticks(range(len(sample_ids)), labels, rotation=30, ha="right")
            axis.spines[["top", "right"]].set_visible(False)
        _save(figure, output, "ORIGINAL_DAS_PSI_AND_MOTIF_LOCUS_TPM_PER_SAMPLE")
    contexts = []
    with (output/"CASE_REAL_HIT_CONTEXTS_RNA.fasta").open("w", encoding="ascii") as handle:
        for row in selected_hits.to_dict("records"):
            form = rna[rna.transcript_id.eq(row["transcript_id"])].iloc[0]
            start, end = max(0, row["RNA_start0"]-40), min(len(form.sequence_RNA), row["RNA_end0"]+40)
            context = form.sequence_RNA[start:end]
            handle.write(f'>{row["transcript_id"]} RNA_start0={start} RNA_end0={end} motif_start0={row["RNA_start0"]} motif_end0={row["RNA_end0"]}\n{context}\n')
            contexts.append({"transcript_id": row["transcript_id"], "context_RNA_start0": start, "context_RNA_end0": end,
                             "context_RNA": context, "motif_start_within_context0": row["RNA_start0"]-start,
                             "motif_end_within_context0": row["RNA_end0"]-start})
    pd.DataFrame(contexts).to_csv(output/"CASE_REAL_HIT_CONTEXTS_AND_MOTIF_OFFSETS.tsv", sep="\t", index=False)
    write_json(output/"ACTUAL_SEQUENCE_CASE_AND_BOUNDARIES.json", {"symbol": symbol, "case": case, "example": example,
        "CLIP": "not_computed_no_matched_case_CLIP_supplied", "new_motif_P_FDR": "not_computed",
        "not_claimed": ["two significant DTU transcripts", "experimental motif disruption", "sample-specific RNA variants", "validated RBP occupancy"]})
    fractions = usage.groupby(["transcript_id", "arm"]).transcript_fraction.mean().unstack("arm")
    fractions["delta_mean_fraction"] = fractions.numerator-fractions.denominator
    fractions = fractions.reset_index().merge(original_dtu[["transcript_id", "stageR_gene_screen_q", "stageR_transcript_confirmation_q"]], on="transcript_id", how="left")
    fractions.to_csv(output/"ALL_TRANSCRIPTS_MEAN_USAGE_AND_SEPARATE_DTU_Q.tsv", sep="\t", index=False)
    gene_means = usage.groupby(["sample_id", "arm"]).gene_TPM.first().reset_index().groupby("arm").gene_TPM.mean().to_dict()
    site_means = selected_samples.groupby("arm")[["RNA_fraction_with_locus", "transcript_TPM_with_locus"]].mean().to_dict()
    unique_counts = loci.groupby(["model_kind", "pair_state"]).size().rename("model_locus_count").reset_index()
    unique_counts.to_csv(output/"ALL_MODEL_LOCUS_STATES_NOT_INDEPENDENT_RBPS.tsv", sep="\t", index=False)
    paired = fractions[fractions.transcript_id.isin([case["increased_transcript"], case["decreased_transcript"]])]
    contexts_html = []
    for row in contexts:
        a, b = row["motif_start_within_context0"], row["motif_end_within_context0"]
        word = row["context_RNA"]
        marked = html.escape(word[:a])+"<mark>"+html.escape(word[a:b])+"</mark>"+html.escape(word[b:])
        contexts_html.append("<p><strong>"+html.escape(row["transcript_id"])+"</strong> RNA "+str(row["context_RNA_start0"])+"–"+str(row["context_RNA_end0"])+" (0-based, half-open)</p><pre>"+marked+"</pre>")
    attachments = "".join('<li><a href="'+html.escape(p.name)+'">'+html.escape(p.name)+"</a></li>" for p in sorted(output.iterdir()) if p.is_file() and p.suffix not in {".png", ".svg", ".html"})
    intro = f"{symbol}; {plan['comparison_id']}; {plan['genome_build']}. Retains all {len(rna)} isoforms, {usage.sample_id.nunique()} existing biological samples and {complete['models']} motif models."
    text = f"""# Katana MotifTrace: RNA usage and actual motif sequences

{intro}

The case was selected before motif scanning using DTU/DEXSeq/DAS coordinate support and transcript-fraction changes, rather than motif positivity.
Mean fraction changes: increased transcript {case['increased_transcript']}, {case['increased_delta_fraction']:.6f}; decreased transcript {case['decreased_transcript']}, {case['decreased_delta_fraction']:.6f}.
Each transcript retains its own stageR confirmation q. A non-significant increased transcript cannot borrow a gene-level or another transcript's q.

Sequence example selected by the declared rule: {example['RBP']}; model {model}; original pattern {pattern}; actual sequence match {example['hit_RNA']}; {example['chr']} {example['strand']} {example['genomic_blocks_json']}.
Original DEXSeq bin {selected_bin.featureID}: P={selected_bin.DEXSeq_P:.8g}, FDR={selected_bin.DEXSeq_bin_FDR:.8g}.
Weighting uses the complete locus-presence relationship across the 12 RNA forms; mean RNA-fraction availability change is {example['delta_fraction_availability']:.6f}.
Group mean gene TPM={gene_means}. Group mean motif-locus availability={site_means}.

Purpose: connect evidence for one gene to specific RNA forms and coordinates, separating composition from total expression.
Limitations: no new motif-enrichment P/FDR, matched CLIP, affinity or causal experiments. PWM thresholds are for candidate display. Sequences derive from reference GRCm39, not directly measured sample-specific allelic RNA.
A motif may correspond to multiple RBPs; model counts are not independent RBP counts. Original species and database provenance are retained; predictions from human models are not mouse experimental validation.
Original DAS/stageR/DEXSeq/DEG statistics remain separate. No combined P-values, requantification, background changes, or submission/cancellation of other tasks.
"""
    (output/"REAL_CASE_CONCLUSION_EN.md").write_text(text, encoding="utf-8")
    page = '''<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Katana MotifTrace | RNA usage and motif consequences</title><style>body{margin:0;background:linear-gradient(130deg,#e8efe7,#faf3e7);color:#183d3b;font:16px/1.75 Georgia,"Microsoft YaHei",serif}main{max-width:1200px;margin:auto;padding:30px}section{background:#fffdf8;border:1px solid #ced9ce;padding:22px;margin:20px 0}h1{font-size:34px}h2{font-size:23px}table{border-collapse:collapse;font:13px/1.5 sans-serif;display:block;overflow:auto;max-width:100%}td,th{padding:8px;border-bottom:1px solid #d6dfd5;white-space:nowrap}img{max-width:100%;height:auto}a{color:#126765}pre{white-space:pre-wrap;overflow-wrap:anywhere;font:16px/1.8 Consolas,monospace;background:#edf2e7;padding:16px}mark{background:#ffd48a}.boundary{border-left:5px solid #b66832;padding-left:15px}code{overflow-wrap:anywhere}</style><main>'''
    page += "<h1>Katana MotifTrace</h1><h2>From RNA usage to actual motif sequences</h2><p>"+html.escape(intro)+"</p>"
    page += '<section><h2>Results and interpretation boundaries</h2><p>This prepared case includes reloaded quantification, concordance across three sequence references and a scan of all 5712 models. It describes motif-bearing RNA forms and changes in their usage; it is not a new population-enrichment test.</p><p class="boundary">Increased and decreased transcripts retain their own DTU q; gene-level significance is not transferred to transcripts. Motif candidates, DEXSeq overlap and original DAS statistics are not three independent binding validations. Matched CLIP: not computed; no matched input was provided.</p>'+_table(paired)+'<img src="PER_SAMPLE_USAGE_EXPRESSION_AND_MOTIF_AVAILABILITY.png"></section>'
    page += '<section><h2>RNA sequences, motif loci and original DEXSeq bins</h2><p>Displayed model: '+html.escape(str(example["RBP"]))+'; '+html.escape(model)+'; original pattern '+html.escape(pattern)+'; actual sequence match <strong>'+html.escape(example["hit_RNA"])+"</strong>. This example follows the declared selection rule after the full-model run; it is not independent validation.</p><p>Coordinates (0-based, half-open): "+html.escape(example["chr"]+" "+example["strand"]+" "+example["genomic_blocks_json"])+"</p><img src=\"REAL_RNA_FORMS_AND_MOTIF_POSITIONS.png\"><img src=\"NATIVE_CONSENSUS_MOTIF_DISPLAY.png\">"+"".join(contexts_html)+"</section>"
    page += '<section><h2>Separate P/FDR for each evidence layer</h2><h3>stageR transcript confirmation: all isoforms; missing is not zero</h3>'+_table(fractions)+'<h3>Original DEG statistics</h3>'+_table(original_deg)+'<h3>Original DEXSeq bin statistics</h3>'+_table(bins[["featureID", "start0", "end0", "DEXSeq_P", "DEXSeq_bin_FDR"]])+'<h3>Coordinate correspondence among motif, DAS and DEXSeq</h3>'+_table(evidence)+'<p>All original DAS P/FDR and sample PSI values are attached. rMATS local-form PSI and Salmon all-isoform fractions have different denominators and need not agree numerically; DEXSeq bin fractions are not PSI. New motif P/FDR: not computed; neither DEXSeq nor DTU q is substituted.</p>'
    if psi_rows:
        page += '<img src="ORIGINAL_DAS_PSI_AND_MOTIF_LOCUS_TPM_PER_SAMPLE.png">'
    page += '</section>'
    page += '<section><h2>Method status and related work</h2><p>This case adds native DTU/DEXSeq/expression adapters, sample-wise weighting across isoforms, shared-locus deduplication, separation of composition and total abundance, and a coordinate evidence chain. Existing capabilities include five-event local RNA reconstruction, junction motif gain/loss, PSI-weighted opportunity and independent window/region statistics.</p><p>Related frameworks include <a href="https://doi.org/10.1093/nar/gkaa237">rMAPS2</a> for motif/CLIP positional enrichment across five event types, <a href="https://www.nature.com/articles/s41467-024-48046-1">MAPP</a> for motif associations with splice/polyadenylation-site usage, and <a href="https://www.bioconductor.org/packages/release/bioc/html/IsoformSwitchAnalyzeR.html">IsoformSwitchAnalyzeR</a> for isoform switches and predicted consequences. Larger libraries, sliding windows or usage weighting alone do not establish novelty. MotifTrace integrates traceable evidence layers and sequence consequences; independent methodological claims require comparative ablation and validation.</p></section>'
    page += '<section><h2>Complete results and review attachments</h2><p>All models, zero hits, non-significant results and NA values are retained. The 5712 deduplicated models are not 5712 independent RBPs. Original species and evidence classes remain explicit; sequence matches are not observed binding in the same tissue.</p><ul>'+attachments+'</ul></section></main></html>'
    (output/"INDEX.html").write_text(page, encoding="utf-8")
    review = {"gene_symbol": symbol, "case": case, "example_RBP": example["RBP"], "model_id": model,
              "hit_RNA": example["hit_RNA"], "genomic_blocks": json.loads(example["genomic_blocks_json"]),
              "chr": example["chr"], "strand": example["strand"], "site_arm_means": site_means, "gene_TPM_arm_means": gene_means,
              "all_rows": complete["rows"], "all_models": complete["models"]}
    write_json(output/"ACTUAL_RESULT_REVIEW.json", review)
    files = {p.name: {"bytes": p.stat().st_size, "SHA256": sha256(p)} for p in sorted(output.iterdir()) if p.is_file() and p.name != "REPORT_MANIFEST.json"}
    write_json(output/"REPORT_MANIFEST.json", {"files": files, "source_prepared_case": str(prepared), "source_complete_SHA256": sha256(output/"USAGE_MOTIF_COMPLETE.json")})
    return review
