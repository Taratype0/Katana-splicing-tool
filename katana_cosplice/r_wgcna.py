from __future__ import annotations

import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd


def exact_wgcna_labels(
    matrix: pd.DataFrame,
    *,
    soft_power: int,
    min_module_size: int,
    merge_cut_height: float,
    deep_split: int,
    pam_respects_dendro: bool,
    rscript: str | Path | None = None,
    event_reliability: pd.Series | None = None,
    reliability_power: float = 0.5,
) -> np.ndarray:
    """Run the production signed WGCNA/TOM backbone in official R WGCNA."""

    if rscript is None:
        from katana_cosplice.benchmark import _rscript_with_wgcna

        executable, diagnostic = _rscript_with_wgcna()
        if executable is None:
            prefix = Path(sys.prefix)
            env_roots = [prefix / "envs"]
            if prefix.parent.name.lower() == "envs":
                env_roots.append(prefix.parent)
            candidates = [prefix / "Scripts" / "Rscript.exe"]
            for root in env_roots:
                if root.is_dir():
                    candidates.extend(root.glob("*/Scripts/Rscript.exe"))
            executable = next(
                (str(path) for path in candidates if path.exists()), None
            )
        if executable is None:
            raise RuntimeError(
                "Production CoSplice requires R WGCNA. " + diagnostic
            )
        rscript = executable
    executable = Path(rscript)
    if not executable.exists():
        raise FileNotFoundError(f"Configured Rscript does not exist: {executable}")

    with tempfile.TemporaryDirectory(prefix="katana_cosplice_wgcna_") as temp_dir:
        temp = Path(temp_dir)
        matrix_path = temp / "matrix.csv"
        reliability_path = temp / "reliability.csv"
        output_path = temp / "modules.csv"
        script_path = temp / "network.R"
        matrix.T.to_csv(matrix_path, index=True, index_label="sample_id")
        if event_reliability is None:
            reliability = pd.Series(1.0, index=matrix.index, name="reliability")
        else:
            reliability = (
                pd.to_numeric(event_reliability.reindex(matrix.index), errors="coerce")
                .fillna(0.0)
                .clip(0.0, 1.0)
                .rename("reliability")
            )
        reliability.to_csv(reliability_path, index=True, index_label="event_key")
        script_path.write_text(
            """
args <- commandArgs(trailingOnly=TRUE)
suppressPackageStartupMessages(library(WGCNA))
suppressPackageStartupMessages(library(dynamicTreeCut))
dat <- read.csv(args[1], row.names=1, check.names=FALSE)
rel <- read.csv(args[2], row.names=1, check.names=FALSE)[, 1]
rel <- rel[colnames(dat)]
corMat <- bicor(dat, use='pairwise.complete.obs')
weightPower <- as.numeric(args[9])
if (weightPower <= 0) {
  net <- blockwiseModules(
    dat,
    power=as.numeric(args[4]),
    corType='bicor',
    networkType='signed',
    TOMType='signed',
    minModuleSize=as.integer(args[5]),
    mergeCutHeight=as.numeric(args[6]),
    numericLabels=TRUE,
    pamRespectsDendro=as.logical(args[8]),
    deepSplit=as.integer(args[7]),
    maxBlockSize=ncol(dat),
    verbose=0
  )
  colors <- net$colors
} else {
  similarity <- (1 + corMat) / 2
  adj <- similarity ^ as.numeric(args[4])
  relFactor <- outer(rel, rel, function(x, y) (x * y) ^ (weightPower / 2))
  adj <- adj * relFactor
  diag(adj) <- 1
  tom <- TOMsimilarity(adj, TOMType='signed', TOMDenom='min', verbose=0)
  diss <- 1 - tom
  diag(diss) <- 0
  tree <- hclust(as.dist(diss), method='average')
  raw <- cutreeDynamic(
    dendro=tree,
    distM=diss,
    deepSplit=as.integer(args[7]),
    pamRespectsDendro=as.logical(args[8]),
    minClusterSize=as.integer(args[5]),
    verbose=0
  )
  colors <- labels2colors(raw)
  if (length(unique(colors[colors != 'grey'])) > 1) {
    merged <- mergeCloseModules(
      dat,
      colors,
      cutHeight=as.numeric(args[6]),
      verbose=0
    )
    colors <- merged$colors
  }
}
write.csv(
  data.frame(event_key=colnames(dat), raw_module=colors),
  args[3],
  row.names=FALSE,
  quote=TRUE
)
""".strip(),
            encoding="utf-8",
        )
        from katana_cosplice.benchmark import _r_subprocess_env

        completed = subprocess.run(
            [
                str(executable),
                str(script_path),
                str(matrix_path),
                str(reliability_path),
                str(output_path),
                str(soft_power),
                str(min_module_size),
                str(merge_cut_height),
                str(deep_split),
                "TRUE" if pam_respects_dendro else "FALSE",
                str(reliability_power if event_reliability is not None else 0.0),
            ],
            capture_output=True,
            text=True,
            timeout=1800,
            check=False,
            env=_r_subprocess_env(executable),
        )
        if completed.returncode != 0 or not output_path.exists():
            raise RuntimeError(
                "Official R WGCNA backbone failed: "
                + (completed.stderr.strip() or completed.stdout.strip())
            )
        result = pd.read_csv(output_path).set_index("event_key").reindex(matrix.index)
        if result["raw_module"].isna().any():
            raise RuntimeError("R WGCNA output did not preserve the event universe.")

    raw = result["raw_module"].astype(str).replace({"0": "grey", "0.0": "grey"})
    counts = raw.loc[raw.ne("grey")].value_counts()
    retained = counts.loc[counts.ge(min_module_size)].index.tolist()
    retained.sort(key=lambda value: (-int(counts[value]), str(value)))
    names = {value: f"M{index:02d}" for index, value in enumerate(retained, start=1)}
    return raw.map(names).fillna("grey").to_numpy(dtype=object)
