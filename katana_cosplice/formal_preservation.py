from __future__ import annotations

import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from katana_cosplice.benchmark import _r_subprocess_env
from katana_cosplice.utils import bh_adjust


NETREP_STATISTICS = (
    "avg.weight",
    "coherence",
    "cor.cor",
    "cor.degree",
    "cor.contrib",
    "avg.cor",
    "avg.contrib",
)


@dataclass
class FormalPreservationResult:
    wgcna: pd.DataFrame
    netrep: pd.DataFrame
    audit: pd.DataFrame


def prepare_formal_preservation_inputs(
    reference: pd.DataFrame,
    target: pd.DataFrame,
    assignments: pd.DataFrame,
    *,
    min_module_size: int = 5,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Align one fixed discovery partition to a common external universe."""

    if not reference.index.is_unique or not target.index.is_unique:
        raise ValueError("Formal preservation requires unique event keys.")
    assignment = assignments.set_index("event_key")["consensus_module"].astype(str)
    target_quality = target.notna().sum(axis=1).ge(3) & target.std(
        axis=1, skipna=True
    ).gt(1e-8)
    shared = pd.Index(
        [
            event
            for event in reference.index.astype(str)
            if event in target.index and bool(target_quality.get(event, False))
        ],
        dtype=object,
    )
    if len(shared) < min_module_size:
        raise ValueError(
            f"Only {len(shared)} usable events are shared with the target cohort."
        )
    reference_aligned = reference.reindex(shared).apply(
        lambda row: row.fillna(row.mean()), axis=1
    )
    target_aligned = target.reindex(shared).apply(
        lambda row: row.fillna(row.mean()), axis=1
    )
    labels = assignment.reindex(shared).fillna("grey")
    counts = labels.loc[labels.ne("grey")].value_counts()
    retained = set(counts.loc[counts.ge(min_module_size)].index.astype(str))
    labels = labels.where(labels.isin(retained), "grey")
    label_frame = pd.DataFrame(
        {
            "event_key": shared.astype(str),
            "wgcna_module": labels.to_numpy(dtype=str),
            "netrep_module": labels.replace("grey", "0").to_numpy(dtype=str),
        }
    )
    if not retained:
        raise ValueError("No module has enough shared events for formal preservation.")
    return reference_aligned, target_aligned, label_frame


def _formal_preservation_r_script() -> str:
    return r"""
args <- commandArgs(trailingOnly=TRUE)
suppressPackageStartupMessages(library(WGCNA))
if (!requireNamespace('NetRep', quietly=TRUE)) {
  stop('NetRep is required for formal preservation validation')
}
reference <- read.delim(args[1], row.names=1, check.names=FALSE)
target <- read.delim(args[2], row.names=1, check.names=FALSE)
labels <- read.delim(args[3], stringsAsFactors=FALSE, check.names=FALSE)
events <- labels$event_key
reference <- t(as.matrix(reference[events, , drop=FALSE]))
target <- t(as.matrix(target[events, , drop=FALSE]))
colnames(reference) <- events
colnames(target) <- events
storage.mode(reference) <- 'double'
storage.mode(target) <- 'double'
moduleNames <- sort(unique(labels$wgcna_module[labels$wgcna_module != 'grey']))
moduleMap <- setNames(seq_along(moduleNames), moduleNames)
wgcnaLabels <- as.integer(moduleMap[labels$wgcna_module])
wgcnaLabels[labels$wgcna_module == 'grey'] <- 0L
names(wgcnaLabels) <- events
netrepLabels <- labels$netrep_module
names(netrepLabels) <- events
nPerm <- as.integer(args[6])
power <- as.numeric(args[7])
seed <- as.integer(args[8])
nThreads <- as.integer(args[9])
set.seed(seed)

multiExpr <- list(
  discovery=list(data=reference),
  test=list(data=target)
)
colorList <- list(discovery=wgcnaLabels)
wgcna <- modulePreservation(
  multiExpr,
  colorList,
  referenceNetworks=1,
  nPermutations=nPerm,
  networkType='signed',
  corFnc='bicor',
  randomSeed=seed,
  verbose=0,
  quickCor=0
)
wgcnaZ <- as.data.frame(
  wgcna$preservation$Z$ref.discovery$inColumnsAlsoPresentIn.test
)
wgcnaZ$module_index <- rownames(wgcnaZ)
numericModuleIndex <- suppressWarnings(as.integer(wgcnaZ$module_index))
keepModules <- !is.na(numericModuleIndex) & numericModuleIndex > 0L & numericModuleIndex <= length(moduleNames)
wgcnaZ <- wgcnaZ[keepModules, , drop=FALSE]
wgcnaZ$module <- moduleNames[numericModuleIndex[keepModules]]
write.table(wgcnaZ, args[4], sep='\t', quote=FALSE, row.names=FALSE)

corDiscovery <- bicor(reference, use='pairwise.complete.obs')
corTest <- bicor(target, use='pairwise.complete.obs')
adjDiscovery <- ((1 + corDiscovery) / 2) ^ power
adjTest <- ((1 + corTest) / 2) ^ power
diag(adjDiscovery) <- 1
diag(adjTest) <- 1
netrep <- NetRep::modulePreservation(
  network=list(discovery=adjDiscovery, test=adjTest),
  data=list(discovery=reference, test=target),
  correlation=list(discovery=corDiscovery, test=corTest),
  moduleAssignments=list(discovery=netrepLabels),
  backgroundLabel='0',
  discovery='discovery',
  test='test',
  nThreads=nThreads,
  nPerm=nPerm,
  null='overlap',
  alternative='greater',
  simplify=TRUE,
  verbose=FALSE
)
observed <- as.data.frame(netrep$observed)
pvalues <- as.data.frame(netrep$p.values)
colnames(observed) <- paste0('observed_', colnames(observed))
colnames(pvalues) <- paste0('pvalue_', colnames(pvalues))
netrepTable <- cbind(
  module=rownames(observed),
  nVarsPresent=as.integer(netrep$nVarsPresent[rownames(observed)]),
  propVarsPresent=as.numeric(netrep$propVarsPresent[rownames(observed)]),
  observed,
  pvalues
)
write.table(netrepTable, args[5], sep='\t', quote=FALSE, row.names=FALSE)
write.table(
  data.frame(
    WGCNA=as.character(packageVersion('WGCNA')),
    NetRep=as.character(packageVersion('NetRep'))
  ),
  args[10], sep='\t', quote=FALSE, row.names=FALSE
)
""".strip()


def run_formal_preservation(
    reference: pd.DataFrame,
    target: pd.DataFrame,
    assignments: pd.DataFrame,
    *,
    cohort: str,
    arm: str,
    rscript: str | Path,
    permutations: int = 1_000,
    soft_power: int = 6,
    min_module_size: int = 5,
    random_seed: int = 43,
    n_threads: int = 1,
) -> FormalPreservationResult:
    """Run official WGCNA and NetRep preservation on identical inputs."""

    if permutations < 1:
        raise ValueError("Formal preservation requires at least one permutation.")
    if n_threads < 1:
        raise ValueError("Formal preservation requires at least one thread.")
    executable = Path(rscript)
    if not executable.exists():
        raise FileNotFoundError(f"Configured Rscript does not exist: {executable}")
    reference_aligned, target_aligned, labels = prepare_formal_preservation_inputs(
        reference,
        target,
        assignments,
        min_module_size=min_module_size,
    )
    with tempfile.TemporaryDirectory(prefix="katana_formal_preservation_") as raw_temp:
        temp = Path(raw_temp)
        reference_path = temp / "reference.tsv"
        target_path = temp / "target.tsv"
        labels_path = temp / "labels.tsv"
        wgcna_path = temp / "wgcna.tsv"
        netrep_path = temp / "netrep.tsv"
        script_path = temp / "formal_preservation.R"
        version_path = temp / "versions.tsv"
        reference_aligned.to_csv(reference_path, sep="\t")
        target_aligned.to_csv(target_path, sep="\t")
        labels.to_csv(labels_path, sep="\t", index=False)
        script_path.write_text(_formal_preservation_r_script(), encoding="utf-8")
        completed = subprocess.run(
            [
                str(executable),
                str(script_path),
                str(reference_path),
                str(target_path),
                str(labels_path),
                str(wgcna_path),
                str(netrep_path),
                str(permutations),
                str(soft_power),
                str(random_seed),
                str(n_threads),
                str(version_path),
            ],
            capture_output=True,
            text=True,
            timeout=7_200,
            check=False,
            env=_r_subprocess_env(str(executable)),
            cwd=temp,
        )
        if (
            completed.returncode != 0
            or not wgcna_path.exists()
            or not netrep_path.exists()
            or not version_path.exists()
        ):
            diagnostic = completed.stderr.strip() or completed.stdout.strip()
            raise RuntimeError(f"Formal WGCNA/NetRep preservation failed: {diagnostic}")
        wgcna = pd.read_csv(wgcna_path, sep="\t")
        netrep = pd.read_csv(netrep_path, sep="\t")
        versions = pd.read_csv(version_path, sep="\t").iloc[0]

    wgcna.insert(0, "arm", arm)
    wgcna.insert(0, "cohort", cohort)
    netrep.insert(0, "arm", arm)
    netrep.insert(0, "cohort", cohort)
    pvalue_columns = [column for column in netrep if column.startswith("pvalue_")]
    for column in pvalue_columns:
        netrep[column.replace("pvalue_", "FDR_")] = bh_adjust(
            pd.to_numeric(netrep[column], errors="coerce")
        )
    netrep["significant_statistics_FDR05"] = netrep[
        [column.replace("pvalue_", "FDR_") for column in pvalue_columns]
    ].le(0.05).sum(axis=1)
    netrep["all_seven_statistics_FDR05"] = netrep[
        "significant_statistics_FDR05"
    ].eq(len(NETREP_STATISTICS))
    audit = pd.DataFrame(
        [
            {
                "cohort": cohort,
                "arm": arm,
                "shared_background_events": len(reference_aligned),
                "tested_modules": len(wgcna),
                "permutations": permutations,
                "soft_power": soft_power,
                "threads": n_threads,
                "rscript": str(executable.resolve()),
                "wgcna_version": str(versions["WGCNA"]),
                "netrep_version": str(versions["NetRep"]),
            }
        ]
    )
    return FormalPreservationResult(wgcna=wgcna, netrep=netrep, audit=audit)
