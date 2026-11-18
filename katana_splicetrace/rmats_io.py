"""Direct rMATS directory, .tar.gz archive, TSV, and Parquet loader for katana-splicetrace."""
from __future__ import annotations

from pathlib import Path
import tarfile
import numpy as np
import pandas as pd

EVENT_COORD_COLS = {
    "SE": ["exonStart_0base", "exonEnd", "upstreamES", "upstreamEE", "downstreamES", "downstreamEE"],
    "A3SS": ["longExonStart_0base", "longExonEnd", "shortES", "shortEE", "flankingES", "flankingEE"],
    "A5SS": ["longExonStart_0base", "longExonEnd", "shortES", "shortEE", "flankingES", "flankingEE"],
    "MXE": ["1stExonStart_0base", "1stExonEnd", "2ndExonStart_0base", "2ndExonEnd", "upstreamES", "upstreamEE", "downstreamES", "downstreamEE"],
    "RI": ["riExonStart_0base", "riExonEnd", "upstreamES", "upstreamEE", "downstreamES", "downstreamEE"],
}


def _standardize_rmats_df(d: pd.DataFrame, kind: str, valid_contigs: set[str] | None, comparison_id: str) -> pd.DataFrame:
    if d.empty:
        return pd.DataFrame()
    if valid_contigs is not None:
        d = d[d["chr"].astype(str).isin(valid_contigs) & d["strand"].isin(["+", "-"])].copy()
    else:
        d = d[d["strand"].isin(["+", "-"])].copy()
    if d.empty:
        return pd.DataFrame()
    d["event_type"] = kind
    d["gene_id"] = d["GeneID"].astype(str).str.strip('"').str.split(".").str[0]
    d["gene_symbol"] = d.get("geneSymbol", d["gene_id"]).astype(str).str.strip('"')
    coord_cols = EVENT_COORD_COLS[kind]
    d["event_id"] = (
        kind + ":" + d["chr"].astype(str) + ":" + d["strand"].astype(str) + ":"
        + d[coord_cols].astype(str).agg(":".join, axis=1)
    )
    d = d.drop_duplicates(subset=["event_id"]).copy()
    d["fdr"] = pd.to_numeric(d["FDR"], errors="coerce")
    d["PValue"] = pd.to_numeric(d.get("PValue", d["FDR"]), errors="coerce")
    d["dpsi"] = pd.to_numeric(d["IncLevelDifference"], errors="coerce")

    r1_m, r2_m, base_p = [], [], []
    for a1, b1, a2, b2, inc2 in zip(
        d["IJC_SAMPLE_1"], d["SJC_SAMPLE_1"], d["IJC_SAMPLE_2"], d["SJC_SAMPLE_2"], d["IncLevel2"]
    ):
        x1 = np.fromstring(str(a1), sep=",")
        y1 = np.fromstring(str(b1), sep=",")
        x2 = np.fromstring(str(a2), sep=",")
        y2 = np.fromstring(str(b2), sep=",")
        p2 = np.fromstring(str(inc2).replace("NA", "nan"), sep=",")
        if len(x1) >= 1 and len(y1) >= 1 and len(x2) >= 1 and len(y2) >= 1:
            r1_m.append(float(np.mean(x1 + y1)))
            r2_m.append(float(np.mean(x2 + y2)))
        else:
            r1_m.append(np.nan)
            r2_m.append(np.nan)
        fin_p = p2[np.isfinite(p2)]
        base_p.append(float(np.mean(fin_p)) if len(fin_p) >= 1 else np.nan)

    d["depth"] = np.minimum(r1_m, r2_m)
    d["baseline_PSI"] = base_p
    d["comparison_id"] = comparison_id
    return d


def load_rmats_measurable_events(
    rmats_source: Path | str,
    fasta_path: Path | str | None = None,
    *,
    min_depth: float = 10.0,
    comparison_id: str = "rmats_contrast",
) -> pd.DataFrame:
    """Load measurable splicing events from an rMATS output directory, `.tar.gz` archive, or `.tsv`/`.parquet` file."""
    src = Path(rmats_source)
    valid_contigs: set[str] | None = None
    if fasta_path is not None:
        fai = Path(str(fasta_path) + ".fai")
        if fai.is_file():
            valid_contigs = {ln.split("\t")[0] for ln in fai.read_text(encoding="utf-8").splitlines() if ln.strip()}

    frames: list[pd.DataFrame] = []
    if src.is_dir():
        for kind in EVENT_COORD_COLS:
            cand = src / f"{kind}.MATS.JC.txt"
            if not cand.is_file():
                cand = src / f"{kind}.MATS.JunctionCountOnly.txt"
            if cand.is_file():
                raw = pd.read_csv(cand, sep="\t", low_memory=False)
                frames.append(_standardize_rmats_df(raw, kind, valid_contigs, comparison_id))
    elif str(src).endswith((".tar.gz", ".tgz", ".tar")):
        with tarfile.open(src, "r:*") as tar:
            for member in tar:
                kind = next(
                    (k for k in EVENT_COORD_COLS if member.name.endswith((f"{k}.MATS.JC.txt", f"{k}.MATS.JunctionCountOnly.txt"))),
                    None,
                )
                if kind is None:
                    continue
                fobj = tar.extractfile(member)
                if fobj is not None:
                    raw = pd.read_csv(fobj, sep="\t", low_memory=False)
                    frames.append(_standardize_rmats_df(raw, kind, valid_contigs, comparison_id))
    elif src.suffix == ".parquet":
        return pd.read_parquet(src)
    else:
        return pd.read_csv(src, sep="\t")

    frames = [f for f in frames if not f.empty]
    if not frames:
        raise ValueError(f"No valid rMATS JC event tables found in {src}")
    all_ev = pd.concat(frames, ignore_index=True).drop_duplicates(subset=["event_id"])
    ok = (
        all_ev["fdr"].notna()
        & all_ev["dpsi"].notna()
        & all_ev["baseline_PSI"].notna()
        & (all_ev["depth"] >= float(min_depth))
        & (all_ev["gene_id"] != "")
    )
    return all_ev.loc[ok].reset_index(drop=True)
