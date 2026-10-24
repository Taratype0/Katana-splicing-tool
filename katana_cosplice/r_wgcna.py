from __future__ import annotations

import subprocess
import shutil
import sys
import tempfile
import os
import threading
from collections import deque
from pathlib import Path

import numpy as np
import pandas as pd


def allocated_wgcna_threads(threads: int | None = None) -> int:
    """Use a declared CPU budget, bounded by the scheduler and CPU affinity.

    An unconfigured local call stays single threaded.  In particular, the
    machine's total core count is not a safe allocation for an individual task.
    """
    def positive(value: object) -> int | None:
        try:
            parsed = int(str(value))
        except (TypeError, ValueError):
            return None
        return parsed if parsed > 0 else None

    if threads is not None and positive(threads) != threads:
        raise ValueError("WGCNA threads must be a positive integer.")
    scheduler = [positive(os.environ.get(key)) for key in
                 ("KATANA_CPUS", "SLURM_CPUS_PER_TASK", "PBS_NP", "NSLOTS")]
    scheduler = [value for value in scheduler if value is not None]
    requested = positive(threads) or positive(os.environ.get("KATANA_NETWORK_THREADS"))
    if requested is None:
        requested = min(scheduler) if scheduler else next(
            (value for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")
             if (value := positive(os.environ.get(key))) is not None), 1)
    caps = scheduler.copy()
    if hasattr(os, "sched_getaffinity"):
        try:
            caps.append(len(os.sched_getaffinity(0)))
        except OSError:
            pass
    return max(1, min([requested, *caps]))


def _threaded_r_env(executable: Path, threads: int) -> dict[str, str]:
    from katana_cosplice.benchmark import _r_subprocess_env

    environment = _r_subprocess_env(executable)
    for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                "BLIS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "ALLOW_WGCNA_THREADS"):
        environment[key] = str(threads)
    environment["OMP_DYNAMIC"] = "FALSE"
    return environment


def _run_wgcna(command: list[str], *, environment: dict[str, str], timeout: int) -> tuple[int, str]:
    """Stream progress to the job log while retaining a bounded failure tail."""
    tail: deque[str] = deque(maxlen=100)
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                               text=True, encoding="utf-8", errors="replace", env=environment,
                               bufsize=1)

    def read_progress() -> None:
        assert process.stdout is not None
        for line in process.stdout:
            tail.append(line)
            sys.stderr.write(line)
            sys.stderr.flush()

    reader = threading.Thread(target=read_progress, name="katana-wgcna-progress", daemon=True)
    reader.start()
    try:
        code = process.wait(timeout=timeout)
    except subprocess.TimeoutExpired:
        process.kill()
        process.wait()
        reader.join()
        raise
    finally:
        if process.poll() is not None:
            reader.join()
            if process.stdout is not None:
                process.stdout.close()
    return code, "".join(tail)


def exact_wgcna_labels(
    matrix: pd.DataFrame,
    *,
    soft_power: int,
    min_module_size: int,
    merge_cut_height: float,
    deep_split: int,
    pam_respects_dendro: bool,
    cor_type: str = "bicor",
    signed: bool = True,
    rscript: str | Path | None = None,
    event_reliability: pd.Series | None = None,
    reliability_power: float = 0.5,
    precomputed_correlation: np.ndarray | None = None,
    return_tom: bool = False,
    timeout_seconds: int = 1800,
    threads: int | None = None,
) -> np.ndarray | tuple[np.ndarray, np.ndarray]:
    """Run official WGCNA with the declared correlation and network signedness."""

    if cor_type not in {"bicor", "pearson"}:
        raise ValueError("Official WGCNA module discovery supports explicit bicor or pearson; unsupported methods must not be silently replaced.")
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
    executable = Path(shutil.which(str(rscript)) or rscript).expanduser().absolute()
    if not executable.is_file():
        raise FileNotFoundError(f"Configured Rscript does not exist: {executable}")
    if not matrix.index.is_unique or matrix.empty:
        raise ValueError("WGCNA requires a nonempty matrix with unique event IDs.")
    if reliability_power < 0:
        raise ValueError("reliability_power must be nonnegative.")
    if timeout_seconds <= 0:
        raise ValueError("timeout_seconds must be positive.")
    if return_tom and event_reliability is None:
        event_reliability = pd.Series(1.0, index=matrix.index)
    thread_count = allocated_wgcna_threads(threads)

    with tempfile.TemporaryDirectory(prefix="katana_cosplice_wgcna_") as temp_dir:
        temp = Path(temp_dir)
        matrix_path = temp / "matrix.csv"
        reliability_path = temp / "reliability.csv"
        correlation_path = temp / "correlation.bin"
        output_path = temp / "modules.csv"
        tom_path = temp / "tom.bin"
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
        if precomputed_correlation is not None:
            correlation = np.asarray(precomputed_correlation, dtype=float)
            if correlation.shape != (len(matrix), len(matrix)):
                raise ValueError("precomputed correlation must match the WGCNA event universe")
            if not np.isfinite(correlation).all() or not np.allclose(correlation, correlation.T):
                raise ValueError("precomputed correlation must be finite and symmetric")
            # A symmetric float64 matrix has identical row-major and R
            # column-major byte order after transposition. Symmetry therefore
            # lets R read this compact binary directly without a huge CSV.
            np.asarray(correlation, dtype="<f8", order="C").tofile(correlation_path)
        script_path.write_text(
            """
args <- commandArgs(trailingOnly=TRUE)
suppressPackageStartupMessages(library(WGCNA))
suppressPackageStartupMessages(library(dynamicTreeCut))
threads <- as.integer(args[16])
if (threads > 1) enableWGCNAThreads(nThreads=threads) else disableWGCNAThreads()
progress <- function(stage) {
  message(sprintf("[Katana WGCNA %s] %s", format(Sys.time(), tz="UTC", usetz=TRUE), stage))
}
progress(sprintf("runtime threads=%d; WGCNA=%d; BLAS=%s; OMP=%s", threads,
  WGCNAnThreads(), extSoftVersion()[["BLAS"]], Sys.getenv("OMP_NUM_THREADS")))
progress("read input matrix and event reliability")
dat <- read.csv(args[1], row.names=1, check.names=FALSE)
relData <- read.csv(args[2], row.names=1, check.names=FALSE)
rel <- setNames(relData[[1]], rownames(relData))
rel <- rel[colnames(dat)]
if (anyNA(rel) || any(!is.finite(rel))) stop("Reliability did not align to event IDs")
usePrecomputed <- as.logical(args[10])
weightPower <- as.numeric(args[9])
needsCorMat <- weightPower > 0 || usePrecomputed || as.logical(args[12])
if (usePrecomputed) {
  progress("read precomputed correlation")
  handle <- file(args[11], "rb")
  values <- readBin(handle, what="double", n=ncol(dat)^2, size=8, endian="little")
  close(handle)
  if (length(values) != ncol(dat)^2) stop("Incomplete precomputed correlation matrix")
  corMat <- matrix(values, nrow=ncol(dat), ncol=ncol(dat))
  rownames(corMat) <- colnames(dat)
  colnames(corMat) <- colnames(dat)
} else if (needsCorMat) {
  progress("correlation")
  corMat <- if (args[14] == 'bicor') bicor(dat, use='pairwise.complete.obs', nThreads=threads) else WGCNA::cor(dat, use='pairwise.complete.obs', method='pearson', nThreads=threads)
}
if (weightPower <= 0 && !usePrecomputed && !as.logical(args[12])) {
  progress("official blockwiseModules (single complete event block)")
  net <- blockwiseModules(
    dat,
    power=as.numeric(args[4]),
    corType=args[14],
    networkType=args[15],
    TOMType='signed',
    minModuleSize=as.integer(args[5]),
    mergeCutHeight=as.numeric(args[6]),
    numericLabels=TRUE,
    pamRespectsDendro=as.logical(args[8]),
    deepSplit=as.integer(args[7]),
    maxBlockSize=ncol(dat),
    nThreads=threads,
    verbose=0
  )
  colors <- net$colors
} else {
  progress("adjacency and reliability weights")
  similarity <- if (args[15] == "signed") (1 + corMat) / 2 else abs(corMat)
  adj <- similarity ^ as.numeric(args[4])
  relFactor <- outer(rel, rel, function(x, y) (x * y) ^ (weightPower / 2))
  adj <- adj * relFactor
  diag(adj) <- 1
  progress("TOMsimilarity")
  tom <- TOMsimilarity(adj, TOMType='signed', TOMDenom='min', verbose=0)
  diss <- 1 - tom
  diag(diss) <- 0
  progress("average linkage")
  tree <- hclust(as.dist(diss), method='average')
  progress("dynamic tree cut")
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
    progress("mergeCloseModules")
    merged <- mergeCloseModules(
      dat,
      colors,
      cutHeight=as.numeric(args[6]),
      verbose=0
    )
    colors <- merged$colors
  }
}
if (as.logical(args[12])) {
  progress("write TOM")
  handle <- file(args[13], "wb")
  writeBin(as.double(tom), handle, size=8, endian='little')
  close(handle)
}
progress("write module assignments")
write.csv(
  data.frame(event_key=colnames(dat), raw_module=colors),
  args[3],
  row.names=FALSE,
  quote=TRUE
)
progress("complete")
""".strip(),
            encoding="utf-8",
        )
        returncode, diagnostics = _run_wgcna(
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
                "TRUE" if precomputed_correlation is not None else "FALSE",
                str(correlation_path),
                "TRUE" if return_tom else "FALSE",
                str(tom_path),
                cor_type,
                "signed" if signed else "unsigned",
                str(thread_count),
            ],
            timeout=timeout_seconds,
            environment=_threaded_r_env(executable, thread_count),
        )
        if returncode != 0 or not output_path.exists():
            raise RuntimeError(
                "Official R WGCNA backbone failed: "
                + diagnostics.strip()
            )
        result = pd.read_csv(output_path).set_index("event_key").reindex(matrix.index)
        if result["raw_module"].isna().any():
            raise RuntimeError("R WGCNA output did not preserve the event universe.")
        if return_tom:
            values = np.fromfile(tom_path, dtype="<f8")
            if values.size != len(matrix) ** 2:
                raise RuntimeError("R WGCNA returned an incomplete TOM matrix.")
            exported_tom = values.reshape((len(matrix), len(matrix)), order="F")
            if not np.isfinite(exported_tom).all() or not np.allclose(exported_tom, exported_tom.T):
                raise RuntimeError("R WGCNA returned a non-finite or asymmetric TOM matrix.")

    raw = result["raw_module"].astype(str).replace({"0": "grey", "0.0": "grey"})
    counts = raw.loc[raw.ne("grey")].value_counts()
    retained = counts.loc[counts.ge(min_module_size)].index.tolist()
    retained.sort(key=lambda value: (-int(counts[value]), str(value)))
    names = {value: f"M{index:02d}" for index, value in enumerate(retained, start=1)}
    labels = raw.map(names).fillna("grey").to_numpy(dtype=object)
    return (labels, exported_tom) if return_tom else labels
