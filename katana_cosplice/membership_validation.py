from __future__ import annotations

import pandas as pd


ARM_FULL = "wgcna_all_members"
ARM_MATCHED = "wgcna_top_kme_matched"
ARM_CORE = "cosplice_split_half_core"
ARMS = (ARM_FULL, ARM_MATCHED, ARM_CORE)


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
    visible = assignments.loc[
        assignments["consensus_module"].astype(str).ne("grey")
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
        audit_rows.append(
            {
                "consensus_module": str(module),
                ARM_FULL: len(group),
                ARM_CORE: n_core,
                ARM_MATCHED: len(selected),
                "core_fraction": n_core / len(group) if len(group) else 0.0,
            }
        )
    matched = (
        pd.concat(matched_parts, ignore_index=True)
        if matched_parts
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
    }
    return arms, pd.DataFrame(audit_rows)
