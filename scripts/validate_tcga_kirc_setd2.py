from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests
import seaborn as sns
import statsmodels.api as sm
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests


STUDY = "kirc_tcga_pan_can_atlas_2018"
PROFILE = f"{STUDY}_mutations"
SAMPLE_LIST = f"{STUDY}_sequenced"
SETD2_ENTREZ = 29072


def _patient_id(value: str) -> str:
    parts = str(value).replace("_", "-").split("-")
    return "-".join(parts[:3])


def _fetch_setd2_status(output: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    base = "https://www.cbioportal.org/api"
    headers = {"Accept": "application/json"}
    sample_payload = requests.get(
        f"{base}/sample-lists/{SAMPLE_LIST}", headers=headers, timeout=120
    ).json()
    sample_ids = sample_payload["sampleIds"]
    response = requests.post(
        f"{base}/molecular-profiles/{PROFILE}/mutations/fetch",
        params={"projection": "DETAILED"},
        headers={**headers, "Content-Type": "application/json"},
        json={"entrezGeneIds": [SETD2_ENTREZ], "sampleIds": sample_ids},
        timeout=180,
    )
    response.raise_for_status()
    mutations = pd.DataFrame(response.json())
    if mutations.empty:
        raise RuntimeError("cBioPortal returned no SETD2 mutation records for TCGA KIRC.")
    mutations.to_csv(output / "tcga_kirc_SETD2_mutations.tsv", sep="\t", index=False)

    damaging = {
        "Frame_Shift_Del",
        "Frame_Shift_Ins",
        "Nonsense_Mutation",
        "Splice_Site",
        "Translation_Start_Site",
        "Nonstop_Mutation",
    }
    mutations["patient_id"] = mutations["sampleId"].map(_patient_id)
    mutations["is_truncating_or_splice"] = mutations["mutationType"].isin(damaging)
    by_patient = mutations.groupby("patient_id", observed=True).agg(
        n_SETD2_mutations=("mutationType", "size"),
        SETD2_truncating_or_splice=("is_truncating_or_splice", "max"),
        SETD2_mutation_types=("mutationType", lambda values: ";".join(sorted(set(values)))),
    )
    sequenced = pd.DataFrame({"sample_id": sample_ids})
    sequenced["patient_id"] = sequenced["sample_id"].map(_patient_id)
    status = sequenced.drop_duplicates("patient_id").set_index("patient_id")
    status = status.join(by_patient, how="left")
    status["n_SETD2_mutations"] = status["n_SETD2_mutations"].fillna(0).astype(int)
    status["SETD2_any_mutation"] = status["n_SETD2_mutations"].gt(0)
    status["SETD2_truncating_or_splice"] = (
        status["SETD2_truncating_or_splice"].eq(True).astype(bool)
    )
    status["SETD2_mutation_types"] = status["SETD2_mutation_types"].fillna("")
    status.reset_index().to_csv(output / "tcga_kirc_SETD2_patient_status.tsv", sep="\t", index=False)
    return status, mutations


def _stage_numeric(value: object) -> float:
    text = str(value).upper()
    for number, token in ((4, "IV"), (3, "III"), (2, "II"), (1, "I")):
        if token in text:
            return float(number)
    return np.nan


def _module_scores(matrix: pd.DataFrame, assignments: pd.DataFrame, module_column: str) -> pd.DataFrame:
    scores: dict[str, pd.Series] = {}
    for module, part in assignments.loc[assignments[module_column].ne("grey")].groupby(
        module_column, observed=True
    ):
        events = part["event_key"].astype(str)
        block = matrix.reindex(events).dropna(how="all")
        if len(block) < 3:
            continue
        block = block.T
        block = block.fillna(block.median(axis=0))
        scale = block.std(axis=0, ddof=1).replace(0, np.nan)
        standardized = (block - block.mean(axis=0)) / scale
        score = standardized.mean(axis=1, skipna=True)
        scores[str(module)] = score
    return pd.DataFrame(scores, index=matrix.columns)


def _association_rows(
    scores: pd.DataFrame,
    status: pd.DataFrame,
    clinical: pd.DataFrame,
    *,
    run_id: str,
    event_type: str,
    method: str,
) -> list[dict[str, object]]:
    sample = pd.DataFrame({"sample_id": scores.index.astype(str)})
    sample["patient_id"] = sample["sample_id"].map(_patient_id)
    sample = sample.set_index("patient_id").join(
        status.drop(columns=["sample_id"], errors="ignore"), how="inner"
    )
    clinical = clinical.copy()
    clinical["patient_id"] = clinical["sample_id"].map(_patient_id)
    clinical = clinical.drop_duplicates("patient_id").set_index("patient_id")
    sample = sample.join(
        clinical[["Absolute_Purity_Value", "Pathologic_Tumor_Stage"]], how="left"
    )
    sample["stage_numeric"] = sample["Pathologic_Tumor_Stage"].map(_stage_numeric)
    sample_scores = scores.copy()
    sample_scores["patient_id"] = sample_scores.index.map(_patient_id)
    sample_scores = sample_scores.set_index("patient_id")
    sample = sample.join(sample_scores, how="inner")
    rows: list[dict[str, object]] = []
    for phenotype in ("SETD2_any_mutation", "SETD2_truncating_or_splice"):
        for module in scores.columns:
            frame = sample[[module, phenotype, "Absolute_Purity_Value", "stage_numeric"]].copy()
            frame[module] = pd.to_numeric(frame[module], errors="coerce")
            frame["Absolute_Purity_Value"] = pd.to_numeric(
                frame["Absolute_Purity_Value"], errors="coerce"
            ).replace(0, np.nan)
            frame = frame.dropna(subset=[module, phenotype])
            mutant = frame.loc[frame[phenotype].astype(bool), module]
            wildtype = frame.loc[~frame[phenotype].astype(bool), module]
            if len(mutant) < 3 or len(wildtype) < 10:
                continue
            pooled = np.sqrt(
                ((len(mutant) - 1) * mutant.var(ddof=1) + (len(wildtype) - 1) * wildtype.var(ddof=1))
                / max(1, len(mutant) + len(wildtype) - 2)
            )
            cohen_d = float((mutant.mean() - wildtype.mean()) / pooled) if pooled > 0 else np.nan
            mw_p = float(mannwhitneyu(mutant, wildtype, alternative="two-sided").pvalue)

            model_frame = frame.dropna(
                subset=[module, phenotype, "Absolute_Purity_Value", "stage_numeric"]
            )
            design = pd.DataFrame(
                {
                    "SETD2": model_frame[phenotype].astype(float),
                    "purity": model_frame["Absolute_Purity_Value"].astype(float),
                    "stage": model_frame["stage_numeric"].astype(float),
                },
                index=model_frame.index,
            )
            design = sm.add_constant(design, has_constant="add")
            fit = sm.OLS(model_frame[module].astype(float), design).fit(cov_type="HC3")
            rows.append(
                {
                    "run_id": run_id,
                    "event_type": event_type,
                    "method": method,
                    "module": module,
                    "phenotype": phenotype,
                    "n_mutant": len(mutant),
                    "n_wildtype": len(wildtype),
                    "mean_score_mutant": float(mutant.mean()),
                    "mean_score_wildtype": float(wildtype.mean()),
                    "cohen_d": cohen_d,
                    "mann_whitney_pvalue": mw_p,
                    "adjusted_SETD2_beta": float(fit.params["SETD2"]),
                    "adjusted_SETD2_pvalue": float(fit.pvalues["SETD2"]),
                    "adjusted_model_n": int(fit.nobs),
                    "covariates": "tumor_purity + pathologic_stage",
                }
            )
    return rows


def _render(
    associations: pd.DataFrame,
    preservation: pd.DataFrame,
    functional: pd.DataFrame,
    status: pd.DataFrame,
    output: Path,
) -> None:
    sns.set_theme(style="whitegrid", context="talk")
    fig, axes = plt.subplots(2, 2, figsize=(17, 12))
    counts = pd.Series(
        {
            "Sequenced WT": int((~status["SETD2_any_mutation"]).sum()),
            "Any SETD2 mutation": int(status["SETD2_any_mutation"].sum()),
            "Truncating/splice": int(status["SETD2_truncating_or_splice"].sum()),
        }
    )
    axes[0, 0].bar(counts.index, counts.values, color=["#9CA3AF", "#2D6CDF", "#C44536"])
    axes[0, 0].set_title("TCGA KIRC SETD2 phenotype", fontweight="bold")
    axes[0, 0].set_ylabel("Patients")
    axes[0, 0].tick_params(axis="x", rotation=18)

    plotted = associations.loc[
        associations["phenotype"].eq("SETD2_any_mutation")
        & associations["method"].eq("three_objective_cosplice")
    ].copy()
    plotted["label"] = plotted["event_type"] + " " + plotted["module"]
    plotted = plotted.sort_values("cohen_d")
    colors = np.where(plotted["adjusted_SETD2_FDR"].le(0.05), "#C44536", "#2D6CDF")
    axes[0, 1].barh(plotted["label"], plotted["cohen_d"], color=colors)
    axes[0, 1].axvline(0, color="#17212B", linewidth=1)
    axes[0, 1].set_title("Frozen-module association with SETD2 mutation", fontweight="bold")
    axes[0, 1].set_xlabel("Cohen's d (mutant - WT)")

    part = preservation.loc[preservation["method"].eq("three_objective_cosplice")].copy()
    if not part.empty:
        sns.stripplot(
            data=part,
            x="event_type",
            y="Zsummary_wgcna_like",
            hue="preservation_class",
            palette={"preserved": "#0F766E", "partially_preserved": "#E07A00", "not_preserved": "#9CA3AF"},
            size=8,
            ax=axes[1, 0],
        )
        axes[1, 0].axhline(2, color="#E07A00", linestyle="--", linewidth=1.5)
        axes[1, 0].axhline(10, color="#0F766E", linestyle="--", linewidth=1.5)
    axes[1, 0].set_title("Independent KICH module preservation", fontweight="bold")
    axes[1, 0].set_ylabel("Zsummary-like")

    fcounts = (
        functional.loc[
            functional["method"].eq("three_objective_cosplice") & functional["FDR"].le(0.05)
        ]
        .groupby("event_type", observed=True)
        .size()
        .reindex(sorted(functional["event_type"].unique()), fill_value=0)
    )
    axes[1, 1].bar(fcounts.index, fcounts.values, color="#0F766E")
    axes[1, 1].set_title("Discovery GO/KEGG support", fontweight="bold")
    axes[1, 1].set_ylabel("FDR < 0.05 terms")

    fig.suptitle(
        "Independent human SETD2-loss validation in TCGA KIRC",
        fontsize=22,
        fontweight="bold",
        x=0.02,
        ha="left",
    )
    fig.tight_layout()
    fig.savefig(output, dpi=220, bbox_inches="tight")
    plt.close(fig)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Post-lock SETD2 phenotype validation of frozen TCGA KIRC modules."
    )
    parser.add_argument("--locked-root", type=Path, required=True)
    parser.add_argument("--panel-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)

    locked = args.locked_root.resolve()
    panel = args.panel_root.resolve()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    status, mutations = _fetch_setd2_status(output)
    clinical = pd.read_csv(
        panel / "KIRC" / "prepared" / "KIRC_sample_metadata.tsv", sep="\t", low_memory=False
    )

    run_map = {
        "ES": "primary_ES500_KICH",
        "AA": "sensitivity_AA500_KICH",
        "AD": "sensitivity_AD500_KICH",
        "AP": "sensitivity_AP500_KICH",
        "AT": "sensitivity_AT500_KICH",
        "ME": "sensitivity_ME68_KICH",
        "RI": "sensitivity_RI500_KICH",
    }
    association_rows: list[dict[str, object]] = []
    preservation_frames: list[pd.DataFrame] = []
    functional_frames: list[pd.DataFrame] = []
    for event_type, run_id in run_map.items():
        run = locked / run_id
        matrix = pd.read_csv(
            run / "01_prelock" / "frozen_discovery_matrix.tsv", sep="\t", index_col=0
        )
        method_specs = {
            "three_objective_cosplice": (
                run / "01_prelock" / "three_objective_cosplice" / "discovery_assignments.tsv",
                "consensus_module",
            ),
            "event_level_wgcna": (
                run / "01_prelock" / "event_level_wgcna" / "discovery_assignments.tsv",
                "wgcna_module",
            ),
        }
        for method, (path, module_column) in method_specs.items():
            assignment = pd.read_csv(path, sep="\t")
            scores = _module_scores(matrix, assignment, module_column)
            scores.to_csv(output / f"{run_id}_{method}_module_scores.tsv", sep="\t")
            association_rows.extend(
                _association_rows(
                    scores,
                    status,
                    clinical,
                    run_id=run_id,
                    event_type=event_type,
                    method=method,
                )
            )
        preservation = pd.read_csv(
            run / "02_holdout_evaluation" / "holdout_module_preservation_all_methods.tsv",
            sep="\t",
        )
        preservation.insert(0, "run_id", run_id)
        preservation.insert(1, "event_type", event_type)
        preservation_frames.append(preservation)
        functional = pd.read_csv(
            run / "02_holdout_evaluation" / "local_GO_KEGG_enrichment.tsv", sep="\t"
        )
        functional.insert(0, "run_id", run_id)
        functional.insert(1, "event_type", event_type)
        functional_frames.append(functional)

    associations = pd.DataFrame(association_rows)
    if associations.empty:
        raise RuntimeError("No KIRC module had enough sequenced SETD2-mutant and WT patients.")
    associations["adjusted_SETD2_FDR"] = multipletests(
        associations["adjusted_SETD2_pvalue"], method="fdr_bh"
    )[1]
    associations["mann_whitney_FDR"] = multipletests(
        associations["mann_whitney_pvalue"], method="fdr_bh"
    )[1]
    preservation = pd.concat(preservation_frames, ignore_index=True)
    functional = pd.concat(functional_frames, ignore_index=True)
    associations.to_csv(output / "tcga_kirc_SETD2_module_associations.tsv", sep="\t", index=False)
    preservation.to_csv(output / "tcga_kirc_to_kich_module_preservation.tsv", sep="\t", index=False)
    functional.to_csv(output / "tcga_kirc_module_GO_KEGG_enrichment.tsv", sep="\t", index=False)

    significant = associations.loc[associations["adjusted_SETD2_FDR"].le(0.05)]
    significant_functions = functional.loc[functional["FDR"].le(0.05)]
    _render(
        associations,
        preservation,
        functional,
        status,
        output / "tcga_kirc_SETD2_external_validation.png",
    )
    parameters = {
        "design": "post-lock phenotype validation",
        "study": STUDY,
        "mutation_profile": PROFILE,
        "setd2_entrez": SETD2_ENTREZ,
        "locked_root": str(locked),
        "panel_root": str(panel),
        "module_selection_used_SETD2_status": False,
        "event_types": list(run_map),
        "association_model": "OLS HC3: module_score ~ SETD2 + tumor_purity + pathologic_stage",
    }
    (output / "parameters.json").write_text(json.dumps(parameters, indent=2), encoding="utf-8")
    report = [
        "# Independent human SETD2-loss validation",
        "",
        "Frozen KIRC modules were defined before SETD2 mutation status was downloaded. Mutation status was used only for post-lock association testing.",
        "",
        f"- Sequenced patients: {len(status)}",
        f"- Any SETD2 mutation: {int(status['SETD2_any_mutation'].sum())}",
        f"- Truncating/splice SETD2 mutation: {int(status['SETD2_truncating_or_splice'].sum())}",
        f"- FDR-significant adjusted module associations: {len(significant)}",
        f"- FDR-significant local GO/KEGG rows: {len(significant_functions)}",
        "",
        "## Significant phenotype associations",
        "",
        significant.to_markdown(index=False, floatfmt=".4g") if not significant.empty else "None passed FDR 0.05.",
        "",
        "## Interpretation boundary",
        "",
        "This is an independent human ccRCC cohort with SETD2 mutation annotation. It tests phenotype association and KIRC-to-KICH preservation, but it is not an exact cross-species projection of the mouse Setd2 event catalogue.",
    ]
    (output / "TCGA_KIRC_SETD2_EXTERNAL_VALIDATION_REPORT.md").write_text(
        "\n".join(report), encoding="utf-8"
    )
    print(associations.sort_values("adjusted_SETD2_FDR").head(20).to_string(index=False))
    print(f"Output: {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
