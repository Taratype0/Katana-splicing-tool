from __future__ import annotations

import math

import pandas as pd


ARM_FULL = "wgcna_all_members"
ARM_MATCHED = "wgcna_top_kme_matched"
ARM_CORE = "cosplice_split_half_core"
ARM_CONFIDENCE = "cosplice_confidence_matched"
ARMS = (ARM_FULL, ARM_MATCHED, ARM_CORE, ARM_CONFIDENCE)


def add_continuous_membership_confidence(
    assignments: pd.DataFrame,
    event_reliability: pd.Series | None = None,
    *,
    evidence_floor: float = 0.05,
) -> pd.DataFrame:
    """Add an auditable confidence score without changing module membership."""

    required = {
        "event_key",
        "consensus_module",
        "consensus_membership_confidence",
        "consensus_confidence_margin",
        "bootstrap_q10_abs_kME",
        "bootstrap_kME_sign_consistency",
    }
    missing = required.difference(assignments.columns)
    if missing:
        raise ValueError(
            "Continuous confidence requires membership evidence columns: "
            f"{sorted(missing)}"
        )
    result = assignments.copy()

    def numeric(column: str) -> pd.Series:
        return pd.to_numeric(result[column], errors="coerce").clip(0.0, 1.0)

    own = numeric("consensus_membership_confidence")
    margin = numeric("consensus_confidence_margin")
    q10_kme = numeric("bootstrap_q10_abs_kME")
    sign = numeric("bootstrap_kME_sign_consistency")
    result["confidence_coassignment_component"] = (own * margin).pow(0.5)
    result["confidence_kme_component"] = (q10_kme * sign).pow(0.5)

    split_columns = {
        "split_half_membership_role_a",
        "split_half_membership_role_b",
    }
    if split_columns.issubset(result.columns):
        first = result["split_half_membership_role_a"].astype(str)
        second = result["split_half_membership_role_b"].astype(str)
        first_supported = first.isin(["core", "peripheral"])
        second_supported = second.isin(["core", "peripheral"])
        split = pd.Series(0.05, index=result.index, dtype=float)
        split.loc[first_supported | second_supported] = 0.35
        split.loc[first_supported & second_supported] = 0.75
        split.loc[first.eq("core") & second.eq("core")] = 1.0
        split.loc[
            first.eq("not_evaluated") | second.eq("not_evaluated")
        ] = float("nan")
    else:
        split = pd.Series(float("nan"), index=result.index, dtype=float)
    result["confidence_split_half_component"] = split

    if event_reliability is not None:
        reliability = pd.to_numeric(
            event_reliability.reindex(result["event_key"].astype(str)).set_axis(
                result.index
            ),
            errors="coerce",
        ).clip(0.0, 1.0)
    elif "network_reliability" in result.columns:
        reliability = numeric("network_reliability")
    else:
        reliability = pd.Series(float("nan"), index=result.index, dtype=float)
    result["confidence_read_component"] = reliability

    component_columns = [
        "confidence_coassignment_component",
        "confidence_kme_component",
        "confidence_split_half_component",
        "confidence_read_component",
    ]
    components = result[component_columns].apply(pd.to_numeric, errors="coerce")
    available = components.notna().sum(axis=1)
    clipped = components.clip(lower=evidence_floor)
    logs = clipped.apply(lambda column: column.map(math.log, na_action="ignore"))
    score = logs.sum(axis=1, skipna=True).div(
        available.where(available > 0)
    ).map(math.exp, na_action="ignore")
    assigned = result["consensus_module"].astype(str).ne("grey")
    result["continuous_membership_confidence"] = score.where(
        assigned, 0.0
    ).clip(0.0, 1.0)
    result["confidence_evidence_domains"] = available.where(
        assigned, 0
    ).astype(int)
    result["continuous_confidence_rank_within_module"] = (
        result.loc[assigned]
        .groupby("consensus_module", observed=True)[
            "continuous_membership_confidence"
        ]
        .rank(method="first", ascending=False)
        .reindex(result.index)
    )
    result["continuous_confidence_percentile_within_module"] = (
        result.loc[assigned]
        .groupby("consensus_module", observed=True)[
            "continuous_membership_confidence"
        ]
        .rank(method="average", pct=True, ascending=True)
        .reindex(result.index)
    )
    return result


def build_membership_arms(
    assignments: pd.DataFrame,
) -> tuple[dict[str, pd.DataFrame], pd.DataFrame]:
    """Create full, count-matched kME, and split-half core member sets."""

    required = {
        "event_key",
        "consensus_module",
        "membership_role",
        "full_data_kME",
    }
    missing = required.difference(assignments.columns)
    if missing:
        raise ValueError(
            "Consensus assignments lack required membership columns: "
            f"{sorted(missing)}"
        )
    scored = (
        assignments
        if "continuous_membership_confidence" in assignments.columns
        else add_continuous_membership_confidence(assignments)
    )
    visible = scored.loc[
        scored["consensus_module"].astype(str).ne("grey")
    ].copy()
    visible["abs_full_data_kME"] = (
        pd.to_numeric(visible["full_data_kME"], errors="coerce")
        .abs()
        .fillna(0.0)
    )
    core = visible.loc[
        visible["membership_role"].astype(str).eq("core")
    ].copy()
    matched_parts: list[pd.DataFrame] = []
    confidence_parts: list[pd.DataFrame] = []
    audit_rows: list[dict[str, object]] = []
    for module, group in visible.groupby(
        "consensus_module", observed=True, sort=True
    ):
        n_core = int(
            core["consensus_module"].astype(str).eq(str(module)).sum()
        )
        ranked = group.sort_values(
            ["abs_full_data_kME", "event_key"], ascending=[False, True]
        )
        selected = ranked.head(n_core).copy()
        matched_parts.append(selected)
        confidence_selected = group.sort_values(
            ["continuous_membership_confidence", "event_key"],
            ascending=[False, True],
        ).head(n_core).copy()
        confidence_parts.append(confidence_selected)
        audit_rows.append(
            {
                "consensus_module": str(module),
                ARM_FULL: len(group),
                ARM_CORE: n_core,
                ARM_MATCHED: len(selected),
                ARM_CONFIDENCE: len(confidence_selected),
                "core_fraction": n_core / len(group) if len(group) else 0.0,
            }
        )
    matched = (
        pd.concat(matched_parts, ignore_index=True)
        if matched_parts
        else visible.iloc[0:0].copy()
    )
    confidence_matched = (
        pd.concat(confidence_parts, ignore_index=True)
        if confidence_parts
        else visible.iloc[0:0].copy()
    )

    def compact(frame: pd.DataFrame) -> pd.DataFrame:
        return frame[["event_key", "consensus_module"]].copy().reset_index(
            drop=True
        )

    arms = {
        ARM_FULL: compact(visible),
        ARM_MATCHED: compact(matched),
        ARM_CORE: compact(core),
        ARM_CONFIDENCE: compact(confidence_matched),
    }
    return arms, pd.DataFrame(audit_rows)
