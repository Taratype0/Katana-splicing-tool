"""Cumulative filtering with explicit zero cells and first-exclusion reasons."""
import pandas as pd


def cumulative_attrition(frame, stages, group="comparison_id"):
    groups = pd.Index(sorted(frame[group].dropna().unique()), name=group)
    previous = pd.Series(True, index=frame.index)
    first_failure = pd.Series("retained", index=frame.index, dtype="string")
    initial = frame.groupby(group).size().reindex(groups, fill_value=0)
    previous_counts = initial
    ledgers = []
    for order, (label, raw_mask) in enumerate(stages):
        mask = pd.Series(raw_mask, index=frame.index).fillna(False).astype(bool)
        if (mask & ~previous).any():
            raise ValueError(f"Attrition stage {label} reintroduces excluded rows; use a separate analysis branch")
        first_failure.loc[previous & ~mask] = label
        counts = frame.loc[mask].groupby(group).size().reindex(groups, fill_value=0)
        table = counts.rename("events").to_frame()
        table["input_events"] = initial
        table["excluded_at_stage"] = previous_counts-counts
        table["retained_fraction"] = counts/initial.replace(0,float("nan"))
        table["stage"], table["stage_order"] = label, order
        ledgers.append(table.reset_index())
        previous, previous_counts = mask, counts
    columns = [group,"events","input_events","excluded_at_stage","retained_fraction","stage","stage_order"]
    ledger = pd.concat(ledgers, ignore_index=True) if ledgers else pd.DataFrame(columns=columns)
    ledger["row_id"] = ledger[group].astype(str)+"::"+ledger.stage.astype(str)
    return ledger, first_failure
