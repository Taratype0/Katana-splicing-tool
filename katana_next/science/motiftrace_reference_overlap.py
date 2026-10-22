"""Fixed event-window overlap audit; missing reference coordinates stay unknown."""
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

from katana_next.science.motiftrace_external import PeakIndex, interval_from_rna, overlapping_hit_positions
from katana_next.science.rmaps_sliding_comparison import author_intervals, window_design


def audit_windows(windows, groups, sequences, score_root, references):
    design = window_design()
    summaries, events, hits = [], [], []
    indexes = {rbp: PeakIndex(reference["rows"]) for rbp, reference in references.items()}
    for item in windows.to_dict("records"):
        window = design[design.region_index.eq(item["region_index"]) & design.native_position.eq(item["native_position"])].iloc[0]
        if not np.array_equal([window.window_start_RNA_relative, window.window_end_RNA_relative],
                              [item["window_start_RNA_relative"], item["window_end_RNA_relative"]]):
            raise ValueError("Frozen physical window geometry differs")
        foreground = "up" if item["direction"] == "inclusion_up" else "dn"
        case = str(item["contrast"])
        reference = references.get(item["RBP"])
        index = indexes.get(item["RBP"])
        count = {g: dict(N=0, peak=0, motif_peak=0, genes=set()) for g in ("fg", "bg")}
        with np.load(score_root/f"m{int(item['motif_index']):03d}.npz", allow_pickle=False) as cache:
            masks = {g: np.isfinite(cache[g][:,int(window.name)]) for g in (foreground,"bg")}
        for group, label in ((foreground,"fg"),("bg","bg")):
            if len(masks[group]) != len(groups[group]):
                raise ValueError("Cached opportunity matrix event count differs")
            for ordinal, event in enumerate(groups[group].itertuples(index=False)):
                seq = sequences[group][window.fasta_region][ordinal]
                span = int(design.region_index.eq(item["region_index"]).sum())
                left = int(window.native_position) if window.side == "positive" else len(seq)-span-49+int(window.native_position)
                right = left+50
                good = 0 <= left < right <= len(seq) and set(seq[left:right]) <= set("ACGT")
                if good != bool(masks[group][ordinal]):
                    raise ValueError("Cached complete-window eligibility changed")
                if not good:
                    continue
                block = author_intervals(event)[window.fasta_region]
                start, stop = interval_from_rna(block,event.strand,left,right)
                observed = None if index is None else index.overlap(event.chr,event.strand,start,stop)
                motif_peak = None if index is None else False
                for a,b in overlapping_hit_positions(seq,item["pattern"]):
                    if a >= right or b <= left:
                        continue
                    mstart,mstop = interval_from_rna(block,event.strand,a,b)
                    istart,istop = interval_from_rna(block,event.strand,max(a,left),min(b,right))
                    match = None if index is None else index.overlap(event.chr,event.strand,istart,istop)
                    if match:
                        motif_peak = True
                    if label == "fg":
                        hits.append(dict(contrast=case,window_id=item["window_id"],RBP=item["RBP"],
                            dataset=reference["label"] if reference else None,event_id=event.event_id,gene_id=event.gene_id,
                            gene_symbol=event.gene_symbol,chr=event.chr,strand=event.strand,window_start0=start,window_stop=stop,
                            motif_start0=mstart,motif_stop=mstop,actual_motif_RNA=seq[a:b].replace("T","U"),
                            window_RNA=seq[left:right].replace("T","U"),pattern=item["pattern"],
                            experimental_peak_overlaps_motif_in_window=match))
                count[label]["N"] += 1
                count[label]["genes"].add(event.gene_id)
                if index is not None:
                    count[label]["peak"] += int(observed)
                    count[label]["motif_peak"] += int(motif_peak)
                events.append(dict(contrast=case,window_id=item["window_id"],RBP=item["RBP"],group=label,
                    event_id=event.event_id,gene_id=event.gene_id,gene_symbol=event.gene_symbol,
                    chr=event.chr,strand=event.strand,start=start,end=stop,window_RNA=seq[left:right].replace("T","U"),
                    CLIP_overlap=observed,motif_CLIP_overlap=motif_peak))
        fg,bg = count["fg"],count["bg"]
        if (fg["N"],bg["N"]) != (int(item["n_foreground"]),int(item["n_background"])):
            raise ValueError("Original window foreground/background denominators differ")
        p = float(fisher_exact([[fg["peak"],fg["N"]-fg["peak"]],[bg["peak"],bg["N"]-bg["peak"]]],alternative="greater").pvalue) if index is not None and min(fg["N"],bg["N"])>=20 else np.nan
        row = {key:item[key] for key in ("contrast","window_id","RBP","pattern","direction","region","original_author_P","original_author_q","original_MotifTrace_P","original_MotifTrace_q")}
        row.update(relative_start=int(window.window_start_RNA_relative),relative_stop=int(window.window_end_RNA_relative),
            dataset=reference["label"] if reference else None,quality=reference["quality"] if reference else "NA_no_usable_mouse_window_coordinates",
            foreground_N=fg["N"],background_N=bg["N"],foreground_genes=len(fg["genes"]),background_genes=len(bg["genes"]),
            foreground_CLIP=fg["peak"] if index is not None else np.nan,background_CLIP=bg["peak"] if index is not None else np.nan,
            foreground_motif_CLIP=fg["motif_peak"] if index is not None else np.nan,
            background_motif_CLIP=bg["motif_peak"] if index is not None else np.nan,external_CLIP_raw_P=p)
        summaries.append(row)
    return pd.DataFrame(summaries),pd.DataFrame(events),pd.DataFrame(hits)
