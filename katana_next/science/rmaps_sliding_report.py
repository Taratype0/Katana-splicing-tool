"""Readable event-level RNA maps and exact native-consensus sequence exports."""
from __future__ import annotations

import csv
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
from katana_next.science.rmaps_sliding_comparison import REGION_SPECS, author_intervals
from katana_next.science.se_rbp import sequence_from_blocks
from katana_next.science.se_rbp_density import consensus_pwm
from katana_next.science.se_rbp_report import _logo


Q = "q_primary_per_arm_both_contrasts_all_positions"
MASK_Q = "q_secondary_canonical_mask_per_arm"
NATIVE = "author_rMAPS3_counts"
COVERAGE = "Katana_overlapping_coverage"


def read_fasta(path):
    strings = []
    with Path(path).open() as stream:
        for header in stream:
            if not header.startswith(">"):
                raise ValueError("Malformed native FASTA")
            strings.append(next(stream).strip())
    return strings


def build_report(root, out, tests, motifs, freeze):
    root, out = Path(root), Path(out)
    summary = json.loads((out/"RESULTS.json").read_text())
    native_and_coverage = tests[tests.arm.isin([NATIVE, COVERAGE])]
    # Same method families are used in both cases, with every position retained.
    comparison_rows = []
    for contrast, case in native_and_coverage.groupby("contrast", observed=True):
        a = case[case.arm.eq(NATIVE)].sort_values(["motif_index", "direction", "region_index", "native_position"])
        b = case[case.arm.eq(COVERAGE)].sort_values(["motif_index", "direction", "region_index", "native_position"])
        if not np.array_equal(a[["motif_index", "direction", "region_index", "native_position"]].to_numpy(),
                              b[["motif_index", "direction", "region_index", "native_position"]].to_numpy()):
            raise ValueError("Cross-tool comparison lost shared hypothesis identity")
        pa, pb = a.raw_P.to_numpy(float), b.raw_P.to_numpy(float)
        sa, sb = a[Q].lt(.05).to_numpy(), b[Q].lt(.05).to_numpy()
        comparison_rows.append(dict(contrast=str(contrast), common_window_hypotheses=len(a),
            significant_both=int((sa & sb).sum()), native_only=int((sa & ~sb).sum()),
            coverage_only=int((~sa & sb).sum()), neither=int((~sa & ~sb).sum()),
            median_absolute_log10_P_difference=float(np.nanmedian(np.abs(np.log10(np.maximum(pa, 1e-300))-np.log10(np.maximum(pb, 1e-300)))))))
    pd.DataFrame(comparison_rows).to_csv(out/"CROSS_TOOL_WINDOW_AGREEMENT.tsv", sep="\t", index=False)
    images = out/"figures"
    images.mkdir()
    logos = out/"logos"
    logos.mkdir()
    selections = []
    for contrast, case in native_and_coverage.groupby("contrast", observed=True):
        display = set(motifs.loc[motifs.RBP.str.casefold().isin(["ptbp1", "srsf3", "srsf10", "nova2"]), "motif_index"])
        top_names = case.sort_values([Q, "raw_P", "motif_index"]).drop_duplicates("RBP").head(3).RBP
        display.update(motifs.loc[motifs.RBP.isin(top_names), "motif_index"])
        for index in sorted(display):
            selections.append((str(contrast), int(index)))
    figure_manifest = []
    region_names = ["UpstreamIntron", "TargetExon_5prime", "TargetExon-3prime", "DownstreamIntron"]
    for contrast, index in selections:
        motif = motifs.iloc[index]
        subset = native_and_coverage[native_and_coverage.contrast.eq(contrast) & native_and_coverage.motif_index.eq(index)]
        fig, axes = plt.subplots(3, 4, figsize=(17, 8.2), layout="constrained")
        for column, region in enumerate(region_names):
            own = subset[subset.region.eq(region)]
            for row, arm in ((0, NATIVE), (1, COVERAGE)):
                arm_data = own[own.arm.eq(arm)]
                for direction, color, label in (("inclusion_up", "#bc4749", "Inclusion up"), ("inclusion_down", "#2877ad", "Inclusion down")):
                    data = arm_data[arm_data.direction.eq(direction)].sort_values("native_position")
                    axes[row, column].plot(data.window_center, data.mean_foreground, color=color, label=label, lw=1.4)
                data = arm_data[arm_data.direction.eq("inclusion_up")].sort_values("native_position")
                axes[row, column].plot(data.window_center, data.mean_background, color="#222222", label="Same background", lw=1.2)
                axes[row, column].set_ylim(bottom=0)
            for arm, style in ((NATIVE, "--"), (COVERAGE, "-")):
                for direction, color in (("inclusion_up", "#bc4749"), ("inclusion_down", "#2877ad")):
                    data = own[own.arm.eq(arm) & own.direction.eq(direction)].sort_values("native_position")
                    y = -np.log10(np.maximum(data.raw_P, 1e-300))
                    axes[2, column].plot(data.window_center, y, color=color, ls=style, lw=1.2,
                        label=("Author count" if arm == NATIVE else "Covered-base")+" / "+direction)
                    sig = data[Q].lt(.05)
                    axes[2, column].scatter(data.loc[sig, "window_center"], y[sig], color=color, s=9)
            axes[0, column].set_title(region, fontsize=10)
            axes[2, column].set_xlabel("50nt window centre relative to exon edge (nt)", fontsize=8)
            if "Intron" in region:
                eligibility = own.drop_duplicates("native_position").sort_values("native_position")
                excluded = eligibility[~eligibility.canonical_mask_eligible]
                if len(excluded):
                    for row in range(3):
                        axes[row, column].axvspan(excluded.window_center.min()-.5, excluded.window_center.max()+.5, color="#888888", alpha=.13)
            for row in range(3):
                axes[row, column].spines[["top", "right"]].set_visible(False)
        axes[0, 0].set_ylabel("Author: mean hit count/window")
        axes[1, 0].set_ylabel("Katana: covered-base fraction")
        axes[2, 0].set_ylabel("-log10(raw P); dots: primary q<0.05")
        axes[0, 0].legend(fontsize=7)
        axes[2, 0].legend(fontsize=6)
        fig.suptitle(f"{contrast} / {motif.RBP} / native pattern {motif.pattern}\n"
                     "All DAS events; identical background, genome, motif and 50nt/1nt grid. Grey: windows excluded only in canonical-mask sensitivity.", fontsize=11)
        name = contrast+f"__m{index:03d}"
        fig.savefig(images/(name+".png"), dpi=130)
        fig.savefig(images/(name+".svg"))
        plt.close(fig)
        _logo(consensus_pwm(motif.pattern), motif.RBP+" / native consensus schematic, not experimental PWM", logos/(f"m{index:03d}.png"))
        figure_manifest.append(dict(contrast=contrast, motif_index=index, RBP=motif.RBP, pattern=motif.pattern,
            plot="figures/"+name+".png", logo=f"logos/m{index:03d}.png",
            selection="All prespecified mechanistic motifs plus first three RBP names by primary-q/raw-P for display only"))
    (out/"FIGURE_MANIFEST.json").write_text(json.dumps(figure_manifest, indent=2), encoding="utf-8")
    export_hits(root, out, native_and_coverage, motifs, freeze)
    rows = pd.DataFrame(summary["results"])
    rows.to_csv(out/"SUMMARY_BY_CONTRAST_AND_ARM.tsv", sep="\t", index=False)
    sig = native_and_coverage[native_and_coverage[Q].lt(.05)]
    top = sig if len(sig) else native_and_coverage.sort_values("raw_P").head(25)
    top = top.sort_values([Q, "raw_P"]).groupby(["contrast", "arm", "RBP", "direction", "region"], observed=True).head(1).head(40)
    top_columns = ["contrast", "arm", "RBP", "pattern", "direction", "region", "native_position",
                   "window_start_RNA_relative", "window_end_RNA_relative", "mean_foreground", "mean_background", "raw_P", Q, MASK_Q]
    top[top_columns].to_csv(out/"DISPLAY_TOP_LOCAL_WINDOWS_NOT_EXTRA_TESTS.tsv", sep="\t", index=False)
    illustrations = "".join(f"<h3>{html.escape(x['contrast']+' / '+x['RBP']+' / '+x['pattern'])}</h3><a href='{x['plot']}'><img src='{x['plot']}'></a><img class='logo' src='{x['logo']}'>"
        for x in figure_manifest if x["RBP"].casefold() == "ptbp1")
    counts = [{"contrast": case["contrast"], "DAS_up_events": case["counts"]["up"]["events"],
        "DAS_down_events": case["counts"]["dn"]["events"], "same_background_events": case["counts"]["bg"]["events"]} for case in freeze["cases"]]
    content = f"""<!doctype html><meta charset='utf-8'><title>Event-level RNA maps: author rMAPS3 vs Katana</title>
<style>body{{font:16px Georgia,serif;max-width:1600px;margin:32px auto;padding:0 25px;background:#f7f6ee;color:#183c40}}table{{display:block;overflow:auto;border-collapse:collapse;background:white;font:12px monospace}}td,th{{padding:7px;border:1px solid #c4cfca}}img{{max-width:100%;background:white}}.logo{{max-width:650px}}a{{color:#08786f}}h1,h2{{font-family:Georgia,serif}}</style>
<h1>Event-level sliding-window comparison</h1>
<p>Native author-run completion and exit codes must be checked against the source execution receipts. This comparison uses the separately labelled author rMAPS3 implementation.
The pinned revision is {freeze['revision']}. This is NOT a claimed native rMAPS2 run, and version number 2.0.1 inside an engine log does not change the repository identity.</p>
<h2>What the inputs and counts mean</h2>{pd.DataFrame(counts).to_html(index=False)}
<p>The unit is an exact SE EVENT, with chromosome, strand and all three exon coordinates. No foreground event was collapsed by gene in these primary RNA maps.
The recorded background entries are non-DAS control EVENTS, not patient pairs or whole-gene sequence alignments.
They are the exact old motif-blind background union, held fixed for every new method. Up and down share that same union in each contrast.</p>
<p>Earlier Katana inference used one lexicographic DAS event per gene plus one matched control event, to reduce same-gene pseudoreplication.
It scanned all selected event sequences but excluded additional events/unmatched genes from that paired test. That was a distinct, more restrictive question, not native rMAPS2.
The secondary old-selection arm now reuses the author's saved count distributions with those prior event IDs and direction-specific controls; it does not rerun any sample measurement.</p>
<h2>Results: all regions and all 1nt positions</h2>{rows.to_html(index=False,float_format=lambda x:f'{x:.5g}')}
<p>Each primary BH family contains all contrasts, inclusion directions, native motif rows and valid regional positions declared in the frozen plan, separately for each arm.
Raw minimum P per RBP is NOT a corrected discovery. Canonical-mask restricted BH and across-all-arm BH are secondary and retain their own labels.
Significant window counts are not independent biological discoveries: adjacent windows, related motifs and events from the same gene are correlated.</p>
<h2>Fairness and real implementation differences</h2><ul>
<li>Shared in both tools: full foreground events, exact background IDs, the same registered reference FASTA and native consensus patterns, intron span200/exon span50, window50/step1 and one-sided Mann-Whitney test. No FIMO/PWM library expansion is mixed into this comparison.</li>
<li>Unmodified author code uses nonoverlapping re.finditer matches and hit counts, with strand/edge indexing conventions. Its raw outputs are retained, not rewritten to match the paper.</li>
<li>Katana's comparison uses all overlapping motif matches and the fraction of the 50 window bases covered. Whole-profile means retain identical inputs but average away position before the same test.</li>
<li>Native code lacks an explicit paper-style 20/6nt splice-signal mask. The prespecified secondary eligibility excludes intronic windows touching those canonical zones identically for both methods, without modifying native P values.</li>
<li>Short/ambiguous windows are NA in the covered-base analysis, not invented zeros. Shared validity requires at least20 foreground and background events. Native raw P remains available outside this common family.</li>
<li>The author's fixed 250nt extraction may extend across a neighbouring exon when an intron is short; every such event is audited. Same-gene dependence is disclosed, not transformed into independent sample replication.</li></ul>
<h2>Independent acceptance</h2><p>{summary['independent_native_P_checks']:,} native P values were independently recomputed from original count distributions; maximum absolute difference {summary['maximum_native_vs_independent_P_difference']:.3g}.
Every author-generated sequence was reloaded from the actual genomic FASTA with ordered event and strand checks. All frozen original input SHA values remain unchanged.
Only local, single-worker/below-normal-priority postprocessing was performed; no SSH, scheduler, quantification, DAS refit or other-task changes.</p>
<h2>Local peaks (all tests remain downloadable)</h2>{top[top_columns].to_html(index=False,float_format=lambda x:f'{x:.5g}')}
<h2>RNA maps</h2><p>Solid red/blue: covered-base up/down; dashed: author count test in the P panel. Dots mark primary-window q&lt;0.05, not per-site binding.
Grey identifies windows removed only from canonical-mask sensitivity. All original author maps and per-event count files remain in the parent native folders.</p>{illustrations}
<p><a href='FIGURE_MANIFEST.json'>Every displayed map and motif logo</a></p>
<h2>Actual motif-to-event sequence correspondence</h2><p><a href='ACTUAL_CONSENSUS_HITS_AND_DAS_COORDINATES.tsv'>Known motif pattern, actual matched RNA, event dPSI/FDR, exact genomic coordinates and source region</a> |
<a href='ACTUAL_CONSENSUS_HITS.rna.fa'>Matched RNA FASTA</a>. Every exported substring is independently reloaded from the genome.
Significant-window matches are labelled by test arm; non-significant mechanistic examples are limited to the first5 lexical events per group and explicitly descriptive.
Consensus logos are schematic allowed bases, NOT PWMs learned from these samples.</p>
<h2>Biological interpretation</h2><p>An intronic motif association is a possible regulator of splicing, not a sequence deleted by exon skipping.
Only cassette-exon motif availability changes directly with inclusion; constitutive flanks and introns are separate regulatory/context layers.
No new splice-junction consequence result is attributed to rMAPS3, which does not test that endpoint here.
A sequence motif does not prove RBP binding or functional effect. YTHDC1 also requires methylation evidence; a DRACH match alone is not m6A or YTHDC1 occupancy.</p>
<h2>Evidence files</h2><ul><li><a href='ALL_WINDOW_AND_PROFILE_TESTS.tsv.gz'>All window/profile raw P, effects, valid event counts and all labelled FDR scopes (gzip)</a></li>
<li><a href='SIGNIFICANT_PRIMARY_WINDOWS.tsv'>All primary significant windows/profiles</a></li>
<li><a href='CROSS_TOOL_WINDOW_AGREEMENT.tsv'>Same-window agreement and disagreements</a></li>
<li><a href='../ANALYSIS_PLAN_AND_FREEZE.json'>Pre-execution event/background/parameter/statistical freeze</a></li>
<li><a href='RESULTS.json'>Actual acceptance results</a></li>
<li><a href='FROZEN_WINDOW_COORDINATES.tsv'>All physical window boundaries and mask eligibility</a></li></ul>
<p>Primary sources: <a href='https://pmc.ncbi.nlm.nih.gov/articles/PMC7319468/'>rMAPS2 paper</a>;
<a href='https://github.com/UofLBioinformatics/rMAPS3/tree/{freeze['revision']}'>Pinned author rMAPS3 source</a>;
<a href='https://meme-suite.org/meme/doc/ame.html'>AME scoring reference for the separate prior whole-region PWM analysis, not this shared-consensus comparison</a>.</p>"""
    (out/"REPORT.html").write_text(content, encoding="utf-8")


def export_hits(root, out, tests, motifs, freeze):
    significant = tests[tests[Q].lt(.05)]
    requests = significant[["contrast", "motif_index", "region", "direction"]].drop_duplicates()
    predefined = []
    for case in freeze["cases"]:
        for motif in motifs[motifs.RBP.str.casefold().isin(["ptbp1", "srsf3", "srsf10", "nova2"])].itertuples(index=False):
            for region in ("TargetExon_5prime", "TargetExon-3prime", "UpstreamIntron", "DownstreamIntron"):
                for direction in ("inclusion_up", "inclusion_down"):
                    predefined.append(dict(contrast=case["contrast"], motif_index=motif.motif_index, region=region, direction=direction))
    requests = pd.concat([requests, pd.DataFrame(predefined)]).drop_duplicates()
    fasta_path = next(x["path"] for x in freeze["input_files"] if x["path"].endswith("GRCm39.primary_assembly.genome.fa"))
    genome = IndexedFasta(Path(fasta_path))
    fields = ["contrast", "RBP", "native_RBP", "pattern", "region", "direction", "event_id", "gene_id", "gene_symbol",
        "DAS_dpsi", "DAS_FDR", "chr", "genomic_strand", "genomic_start0", "genomic_end_exclusive", "RNA_start0", "RNA_stop_exclusive",
        "matched_RNA", "RNA_context", "counted_by_native_nonoverlap_finditer", "export_scope", "significant_arms",
        "overlapping_significant_window_indices", "sequence_consequence_class"]
    count, negative = 0, 0
    with (out/"ACTUAL_CONSENSUS_HITS_AND_DAS_COORDINATES.tsv").open("w", newline="", encoding="utf-8") as handle, \
         (out/"ACTUAL_CONSENSUS_HITS.rna.fa").open("w", encoding="ascii") as fa, genome.fasta.open("rb") as stream:
        writer = csv.DictWriter(handle, fields, delimiter="\t")
        writer.writeheader()
        for request in requests.itertuples(index=False):
            motif = motifs.iloc[int(request.motif_index)]
            subset = significant[significant.contrast.eq(request.contrast) & significant.motif_index.eq(request.motif_index)
                & significant.region.eq(request.region) & significant.direction.eq(request.direction)]
            is_significant = len(subset) > 0
            specification = next(x for x in REGION_SPECS if x[0] == request.region)
            group = "up" if request.direction == "inclusion_up" else "dn"
            case = root/str(request.contrast)
            frame = pd.read_parquet(case/"input"/(group+".events.parquet"))
            sequences = read_fasta(case/"author_native/fasta"/(group+"."+specification[1]+".fasta"))
            indexes = range(len(frame)) if is_significant else range(min(5, len(frame)))
            for i in indexes:
                event = next(frame.iloc[[i]].itertuples(index=False))
                sequence = sequences[i]
                interval = author_intervals(event)[specification[1]]
                native_starts = {match.start() for match in re.finditer(motif.pattern, sequence)}
                for hit in re.finditer("(?=("+motif.pattern+"))", sequence):
                    a, b = hit.start(1), hit.end(1)
                    if is_significant:
                        local_starts = subset.window_start_RNA_relative.to_numpy(float)
                        local_stops = subset.window_end_RNA_relative.to_numpy(float)
                        if specification[2] == "negative":
                            local_starts += len(sequence)
                            local_stops += len(sequence)
                        relevant = (a < local_stops) & (b > local_starts)
                        if not relevant.any():
                            continue
                        windows = ";".join(map(str, sorted(subset.loc[relevant, "native_position"].unique())))
                    else:
                        windows = ""
                    block = (interval[0]+a, interval[0]+b) if event.strand == "+" else (interval[1]-b, interval[1]-a)
                    actual = sequence_from_blocks(genome, stream, event.chr, event.strand, [block])
                    if actual != sequence[a:b]:
                        raise ValueError("Exported consensus substring did not map back to the original genome")
                    count += 1
                    negative += event.strand == "-"
                    writer.writerow(dict(contrast=request.contrast, RBP=motif.RBP, native_RBP=motif.native_RBP,
                        pattern=motif.pattern, region=request.region, direction=request.direction,
                        event_id=event.event_id, gene_id=event.gene_id, gene_symbol=event.gene_symbol,
                        DAS_dpsi=event.dpsi, DAS_FDR=event.fdr, chr=event.chr, genomic_strand=event.strand,
                        genomic_start0=block[0], genomic_end_exclusive=block[1], RNA_start0=a, RNA_stop_exclusive=b,
                        matched_RNA=actual.replace("T", "U"), RNA_context=sequence[max(0,a-15):a].lower().replace("t","u")+"["+actual.replace("T","U")+"]"+sequence[b:b+15].lower().replace("t","u"),
                        counted_by_native_nonoverlap_finditer=a in native_starts,
                        export_scope="all_foreground_hits_overlapping_primary_significant_windows" if is_significant else "first5_lexical_events_descriptive_NOT_significant",
                        significant_arms=";".join(sorted(subset.arm.astype(str).unique())), overlapping_significant_window_indices=windows,
                        sequence_consequence_class="cassette_exon_RNA_availability" if specification[1] == "TargetExon" else "intronic_regulatory_association" if "Intron" in request.region else "constitutive_exon_context"))
                    fa.write(f">hit{count} case={request.contrast} event={event.event_id} RBP={motif.RBP} region={request.region}\n{actual.replace('T','U')}\n")
    (out/"HIT_SEQUENCE_VALIDATION.json").write_text(json.dumps({"exact_real_genome_substrings_checked": count,
        "negative_strand_hits_checked": int(negative), "all_roundtrips_exact": True, "observed_binding_claimed": False}, indent=2), encoding="utf-8")
