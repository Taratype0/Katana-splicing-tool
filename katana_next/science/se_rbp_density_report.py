"""Evidence exports for frozen density sensitivity, not a native rMAPS run."""
from __future__ import annotations

import hashlib
import html
import json
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from katana_cosplice.motif_scanner import IndexedFasta
from katana_next.science.se_rbp import mapped_blocks, sequence_from_blocks
from katana_next.science.se_rbp_density import consensus_pwm
from katana_next.science.se_rbp_report import _logo


def _best_preference_site(window, motifs, indexes, background):
    best = None
    for i in indexes:
        matrix = motifs[i].pwm
        matrix = matrix/matrix.sum(axis=1, keepdims=True)
        matrix = (20*matrix+.1*background)/20.1
        logs = np.log2(matrix/background)
        for a in range(len(window.sequence)-len(matrix)+1):
            b = a+len(matrix)
            if np.isfinite(window.junction_offset) and not a < window.junction_offset < b:
                continue
            word = window.sequence[a:b]
            if set(word)-set("ACGT"):
                continue
            score = float(sum(logs[j, "ACGT".index(base)] for j, base in enumerate(word)))
            if best is None or score > best[0]:
                best = (score, i, a, b)
    return best


def build_density_report(directory, motifs, catalog, patterns):
    directory = Path(directory)
    result = json.loads((directory/"RESULTS.json").read_text())
    freeze = json.loads((directory/"ANALYSIS_FREEZE.json").read_text())
    tests = pd.read_parquet(directory/"all_enrichment_tests.parquet")
    sig = tests[tests.FDR_global_all_three_cases_and_methods.lt(.05)]
    figures = directory/"figures"
    figures.mkdir()
    logos = directory/"logos"
    logos.mkdir()
    selected = []
    # Fixed descriptive display rule; it never changes which hypotheses are tested.
    for case, group in tests.groupby("case", sort=True):
        selected.extend(group.sort_values(["FDR_global_all_three_cases_and_methods", "pvalue", "catalog_scope", "RBP"]).head(2).index)
    display = tests.loc[sorted(set(selected))]
    display.to_csv(directory/"displayed_hypotheses_not_extra_tests.tsv", sep="\t", index=False)
    requests = pd.concat([sig, display]).drop_duplicates(["case", "catalog_scope", "RBP", "evidence_class", "region", "direction"])
    attachment_rows, plot_rows, image_records = [], [], []
    fasta = IndexedFasta(Path(freeze["input_files"][next(i for i, x in enumerate(freeze["input_files"]) if x["path"].endswith("GRCm39.primary_assembly.genome.fa"))]["path"]))
    sequence_checks = 0
    for case, own_requests in requests.groupby("case", sort=True):
        own = directory/case
        windows = pd.read_parquet(own/"sequences.parquet")
        pairs = pd.read_parquet(own/"matching.parquet")
        scores = np.load(own/"all_sequence_scores.npz")
        for test_index, row in own_requests.iterrows():
            is_pwm = row.method == "PWM_average_odds"
            source = catalog if is_pwm else patterns
            indexes = source.index[source.motif_token.isin(row.motif_tokens.split(";"))].tolist()
            matched = pairs[pairs.status.eq("matched") & pairs.direction.eq(row.direction)]
            ids = set(matched.foreground_event) | set(matched.background_event)
            window_subset = windows[windows.region.eq(row.region) & windows.event_id.isin(ids)]
            values = scores["PWM" if is_pwm else "consensus"][indexes]
            means = np.divide(np.nansum(values, axis=0), np.isfinite(values).sum(axis=0),
                              out=np.full(len(windows), np.nan), where=np.isfinite(values).sum(axis=0) > 0)
            lookup = pd.Series(means, index=windows.window_id)
            for window in window_subset.itertuples(index=False):
                plot_rows.append({"hypothesis_index": int(test_index), "case": case, "RBP": row.RBP,
                    "evidence_class": row.evidence_class, "method": row.method, "region": row.region,
                    "direction": row.direction, "role": window.role, "event_id": window.event_id,
                    "gene_id": window.gene_id, "score": lookup[window.window_id],
                    "global_FDR": row.FDR_global_all_three_cases_and_methods})
                if window.role != "foreground":
                    continue
                own_hits = []
                if is_pwm:
                    best = _best_preference_site(window, motifs, indexes, scores["background"])
                    if best is not None:
                        logodds, i, a, b = best
                        own_hits.append((i, a, b, "representative_max_odds_substring_NOT_significant_binding_site", logodds))
                else:
                    for i in indexes:
                        for match in re.finditer("(?=("+patterns.iloc[i].pattern+"))", window.sequence):
                            a, b = match.start(1), match.end(1)
                            if np.isfinite(window.junction_offset) and not a < window.junction_offset < b:
                                continue
                            own_hits.append((i, a, b, "native_fixed_consensus_match_NOT_observed_binding", np.nan))
                for i, a, b, kind, logodds in own_hits:
                    blocks = mapped_blocks(json.loads(window.blocks_json), window.genomic_strand, a, b)
                    word = window.sequence[a:b]
                    with fasta.fasta.open("rb") as stream:
                        actual = sequence_from_blocks(fasta, stream, window.chr, window.genomic_strand, blocks)
                    if actual != word:
                        raise ValueError("Real RNA hit did not roundtrip to original genome blocks")
                    sequence_checks += 1
                    attachment_rows.append({"hypothesis_index": int(test_index), "case": case,
                        "RBP": row.RBP, "evidence_class": row.evidence_class, "catalog_scope": row.catalog_scope,
                        "region": row.region, "direction": row.direction, "global_FDR": row.FDR_global_all_three_cases_and_methods,
                        "event_id": window.event_id, "gene_id": window.gene_id, "gene_symbol": window.gene_symbol,
                        "DAS_dpsi": window.dpsi, "DAS_FDR": window.DAS_fdr,
                        "motif_token": source.iloc[i].motif_token, "motif_id": source.iloc[i].motif_id,
                        "motif_source": source.iloc[i].source, "chr": window.chr, "genomic_strand": window.genomic_strand,
                        "RNA_start0": a, "RNA_stop_exclusive": b, "genomic_blocks_RNA_order": json.dumps(blocks),
                        "matched_RNA": word.replace("T", "U"),
                        "context_RNA": window.sequence[max(0,a-15):a].lower().replace("t","u")+"["+word.replace("T","U")+"]"+window.sequence[b:b+15].lower().replace("t","u"),
                        "relative_position": (a+b)/2/len(window.sequence)*100 if window.region == "cassette_exon" else
                            window.anchor_start+(a+b)/2 if window.form == "pre_mRNA" else (a+b)/2-window.junction_offset,
                        "representation": kind, "representative_site_log2_odds": logodds,
                        "significant_enrichment_hypothesis": bool(row.FDR_global_all_three_cases_and_methods < .05)})
            if test_index not in display.index:
                continue
            tag = f"hypothesis_{int(test_index):05d}"
            logo = logos/(tag+".svg")
            representative = indexes[0]
            pwm = motifs[representative].pwm if is_pwm else consensus_pwm(patterns.iloc[representative].pattern)
            logo_label = "measured/model PWM" if is_pwm else "schematic native consensus; NOT measured PWM"
            _logo(pwm, row.RBP+" / "+logo_label, logo)
            rows = pd.DataFrame([x for x in plot_rows if x["hypothesis_index"] == test_index])
            fig, axes = plt.subplots(1, 2, figsize=(10.5, 3.5), layout="constrained")
            a = rows[rows.role.eq("foreground")].score.dropna()
            b = rows[rows.role.eq("background")].score.dropna()
            axes[0].boxplot([a, b], tick_labels=["DAS gene", "Matched control"], showfliers=False)
            axes[0].set_ylabel("log2(1+average odds)" if is_pwm else "Covered-base fraction")
            hits = pd.DataFrame([x for x in attachment_rows if x["hypothesis_index"] == test_index])
            if len(hits):
                axes[1].hist(hits.relative_position, bins=25, color="#147d78")
            axes[1].set_ylabel("Descriptive matched substrings")
            axes[1].set_xlabel("Cassette exon position (%)" if row.region == "cassette_exon" else "RNA-relative position (nt)")
            fig.suptitle(f"{case}\n{row.RBP} / {row.region} / {row.direction}; global q={row.FDR_global_all_three_cases_and_methods:.3g}", fontsize=9)
            path = figures/(tag+".svg")
            fig.savefig(path)
            plt.close(fig)
            image_records.append({"hypothesis_index": int(test_index), "case": case, "RBP": row.RBP,
                                  "logo": str(logo.relative_to(directory)), "plot": str(path.relative_to(directory)),
                                  "logo_representation": logo_label,
                                  "all_positions_descriptive_not_windowwise_inference": True})
    attachments = pd.DataFrame(attachment_rows)
    attachments.to_csv(directory/"real_RNA_substrings_and_split_genome_coordinates.tsv", sep="\t", index=False)
    pd.DataFrame(plot_rows).to_csv(directory/"display_and_significant_hypothesis_scores.tsv", sep="\t", index=False)
    with (directory/"real_RNA_substrings.fa").open("w", encoding="ascii") as stream:
        for i, row in enumerate(attachment_rows):
            stream.write(f">hit{i} hypothesis={row['hypothesis_index']} case={row['case']} event={row['event_id']} RBP={row['RBP']}\n{row['matched_RNA']}\n")
    (directory/"FIGURE_MANIFEST.json").write_text(json.dumps(image_records, indent=2), encoding="utf-8")
    (directory/"RNA_COORDINATE_VALIDATION.json").write_text(json.dumps({"actual_substrings_roundtripped": sequence_checks,
        "RNA_forward_and_split_junction_coordinates": True, "all_exported_substrings_exact": True,
        "per_site_significance_or_actual_binding_claimed": False}, indent=2), encoding="utf-8")
    counts = tests.groupby(["case", "catalog_scope"]).agg(hypotheses=("pvalue", "size"),
        global_significant=("FDR_global_all_three_cases_and_methods", lambda x: int(x.lt(.05).sum())),
        within_scope_secondary=("FDR_within_declared_scope", lambda x: int(x.lt(.05).sum())),
        minimum_global_q=("FDR_global_all_three_cases_and_methods", "min"))
    counts.to_csv(directory/"scope_result_summary.tsv", sep="\t")
    images = "".join(f"<h3>{html.escape(x['case']+' / '+x['RBP'])}</h3><img src='{x['plot']}'><img class='logo' src='{x['logo']}'>" for x in image_records)
    comparison = pd.DataFrame([
        ["Core question", "Regulatory motif/CLIP position around AS", "Also reports cassette/junction motif availability after SE"],
        ["Short motifs", "Native consensus density; no requirement for FIMO P=1e-4", "Previously fixed-site PWM hits; now adds continuous PWM and native consensus sensitivity"],
        ["Inference", "Sliding windows, Wilcoxon rank-sum", "Unique gene-matched exact sign tests; whole-region tests, global BH across cases/methods"],
        ["Coverage", "SE/MXE/A5SS/A3SS/RI; optional CLIP", "New sequence-consequence mode currently SE only; not a comprehensive upgrade"],
        ["Evidence", "Motif association, not causal binding proof", "Explicit species/source/evidence labels, split-coordinate actual RNA attachments; still not binding/function proof"],
    ], columns=["Aspect", "rMAPS2", "This Katana analysis"])
    inventory = pd.DataFrame([{k: c[k] for k in ["study", "contrast", "all_SE_rows", "measurement_eligible_SE", "DAS_SE", "DAS_SE_genes", "independent_matched_gene_pairs"]} for c in result["cases"]])
    top = sig if len(sig) else tests.sort_values("pvalue").head(15)
    show = top[["case", "catalog_scope", "RBP", "evidence_class", "region", "direction", "matched_pairs", "mean_score_difference", "pvalue", "FDR_global_all_three_cases_and_methods"]]
    contents = f"""<!doctype html><meta charset='utf-8'><title>SE RBP: matched-density sensitivity and rMAPS2 comparison</title>
<style>body{{font:16px Georgia,serif;max-width:1400px;margin:36px auto;padding:0 25px;color:#173c42;background:#f6f5ef}}table{{border-collapse:collapse;font:13px monospace;background:white;display:block;overflow:auto}}td,th{{padding:7px;border:1px solid #bdcbc9}}img{{max-width:100%;background:white}}.logo{{max-width:650px}}a{{color:#096d73}}h1,h2{{font-family:Georgia,serif}}</style>
<h1>SE RBP motif sensitivity: declared study contrasts</h1>
<p><b>{len(sig)} hypotheses pass the prespecified global BH FDR &lt;0.05</b>, out of {result['tested_hypotheses']} tested.
The frozen global family includes all declared contrasts, catalog/method scopes, inclusion directions and regions.
No quantification, DAS test, sample selection or cutoff was rerun or changed. Old results remain untouched.</p>
<h2>What differs from rMAPS2?</h2>{comparison.to_html(index=False,escape=True)}
<p>This is Katana's separately declared matched-density analysis, not execution or a validated reproduction of rMAPS2 or AME.
The continuous average-odds score borrows AME's scoring idea; the exact paired sign test is different.
Density sensitivity changes the endpoint as well as removing site filtering, so any new signal cannot be attributed solely to the old threshold.</p>
<h2>Actual input completeness and denominators</h2>{inventory.to_html(index=False)}
<p>Original complete SE ID sets and dPSI directions were checked against existing rMATS SE.MATS.JC.txt.
Only existing read-supported stable DAS-SE are foreground. All selected events and their sequences remain in the exports;
one lexicographic event per gene and unique background genes enter inference. Missing/ineligible/unmatched events are audited, not silently dropped.
Each declared contrast retains its original matched pairs and numerator/denominator direction, as recorded in the inventory and analysis freeze.</p>
<h2>Database scope and short-motif limitation</h2><p>Short PWM assignments may be unable to reach a fixed stringent site threshold; eligibility must be read from the recorded resource audit.
This sensitivity scores all valid positions without a site cutoff. It uses {len(catalog)} numeric assignments ({catalog.pwm_sha256.nunique()} exact-unique PWMs),
plus {len(patterns)} native finite consensus/binding-sequence assignments. Evidence classes remain separate.
oRNAment is a cross-species reference, not mouse binding validation. RBPmap native motif labels are used as consensuses, not invented numeric matrices.
RBPDB mouse measurements are retained separately; long/variable probe notation is not converted into guessed motifs.</p>
<p><b>Resource limits:</b> the acquisition and eligibility receipts define which numeric matrices and native consensus records are available;
unavailable archives and unrepresented experimental contexts remain missing. Sequence analysis does not supply sample-matched CLIP or methylation evidence.
Not significant therefore does not prove absence of RBP regulation.</p>
<h2>All method scopes, including negatives</h2>{counts.to_html()}
<p>Within-scope BH is secondary sensitivity only. Primary interpretation uses the global column, not whichever method yields the best P.</p>
<h2>{'Globally significant hypotheses' if len(sig) else 'Lowest P values, descriptive only: none globally significant'}</h2>{show.to_html(index=False,float_format=lambda x:f'{x:.6g}')}
<h2>Biological limits</h2><p>Intronic hits suggest possible regulation of splicing, whereas included cassette/junction hits describe sequence availability.
An exon skip does not mutate the retained DNA motif; it changes RNA isoform representation. Local junction analysis does not reconstruct whole-transcript structure or infer affinity.
Finding a sequence motif or DRACH pattern alone does not measure methylation or occupancy of an RNA-binding protein.
Published RBP mechanisms provide context, not automatically validated results in the current samples.
Matching is observational, not randomized; region-specific composition and transcript context can remain confounded.</p>
<h2>Actual RNA, coordinates, logos and maps</h2><p>Figures use the first two hypotheses per case by frozen global-q/P ordering for display only.
All significant groups are retained in tables. Native consensus substrings are real pattern matches; PWM attachments are representative maximum-odds substrings, NOT individually significant binding sites.
Consensus logos are schematic equal-weight letters, not experimental PWMs. Position histograms are descriptive and not sliding-window significance maps.</p>{images}
<h2>Downloadable evidence</h2><ul>
<li><a href='all_enrichment_tests.tsv'>Every hypothesis, P, global and within-scope FDR</a></li>
<li><a href='significant_global_BH005.tsv'>Global FDR significant hypotheses (header-only if none)</a></li>
<li><a href='real_RNA_substrings_and_split_genome_coordinates.tsv'>Actual RNA substrings, source motifs, DAS effects and split genomic blocks</a></li>
<li><a href='real_RNA_substrings.fa'>RNA substring FASTA</a></li>
<li><a href='PWM_catalog.tsv'>Numeric motif provenance</a> / <a href='native_consensus_catalog.tsv'>Native consensus provenance</a></li>
<li><a href='resource_record_eligibility.tsv'>Included and excluded public records</a></li>
<li><a href='ANALYSIS_FREEZE.json'>Pre-scoring protocol, original SHA hashes, sample identities, methods and limitations</a></li>
<li><a href='RESULTS.json'>Machine-readable complete results and checks</a></li></ul>
<h2>Primary sources</h2><p><a href='https://academic.oup.com/nar/article/48/W1/W300/5819937'>rMAPS2 methods</a>;
<a href='https://meme-suite.org/meme/doc/ame.html'>AME average-odds scoring</a>;
<a href='https://github.com/UofLBioinformatics/rMAPS3'>Author-maintained consensus source (pinned in acquisition receipt)</a>;
<a href='https://rbpdb.ccbr.utoronto.ca/download.php'>RBPDB</a>;
<a href='https://rnabiology.ircm.qc.ca/oRNAment/downloads/'>oRNAment</a>;
<a href='https://rbpmap.technion.ac.il/manual.html'>RBPmap</a>;
<a href='https://pubmed.ncbi.nlm.nih.gov/26876937/'>Example primary study of an RBP splicing mechanism</a>.</p>"""
    (directory/"REPORT.html").write_text(contents, encoding="utf-8")
    manifest = []
    for path in sorted(directory.rglob("*")):
        if path.is_file() and path.name not in {"EVIDENCE_MANIFEST.json", "STATUS.json"}:
            digest = hashlib.sha256(path.read_bytes()).hexdigest()
            manifest.append({"relative_path": str(path.relative_to(directory)), "bytes": path.stat().st_size, "sha256": digest})
    (directory/"EVIDENCE_MANIFEST.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return result
