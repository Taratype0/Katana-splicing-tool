from __future__ import annotations

import numpy as np
import pandas as pd


def observation_units(frame: pd.DataFrame) -> pd.DataFrame:
    """Keep donor-condition observations distinct while retaining donor blocks."""
    frame = frame.copy()
    keys = ["donor_id", "condition"]
    for key in keys:
        if frame[key].isna().any() or frame[key].astype(str).str.strip().eq("").any():
            raise ValueError(f"Network observations require nonempty {key}.")
    frame["observation_id"] = frame[keys].astype(str).apply(
        lambda row: "::".join(f"{len(value)}:{value}" for value in row), axis=1
    )
    if frame.duplicated(["observation_id", "event_id"]).any():
        raise ValueError("Duplicate donor-condition-event observations: aggregate technical partitions before networking.")
    return frame


def donor_covariates(frame: pd.DataFrame, *, adjust_for_batch: bool = True) -> tuple[np.ndarray | None, bool, bool]:
    """Paired donor intercepts, plus nonredundant optional batch effects."""
    paired = bool(frame.groupby("donor_id")["condition"].nunique().gt(1).any())
    columns = []
    if paired:
        columns.extend(pd.get_dummies(frame["donor_id"].astype(str), drop_first=True, dtype=float).to_numpy().T)
    batch_added = False
    if adjust_for_batch and "batch" in frame:
        for column in pd.get_dummies(frame["batch"].fillna("").astype(str), drop_first=True, dtype=float).to_numpy().T:
            old = np.column_stack([np.ones(len(frame)), *columns])
            if np.linalg.matrix_rank(np.column_stack([old, column])) > np.linalg.matrix_rank(old):
                columns.append(column)
                batch_added = True
    return (np.column_stack(columns) if columns else None), paired, batch_added


def bootstrap_blocks(metadata: pd.DataFrame, rng: np.random.Generator) -> np.ndarray:
    """Resample independent donors, keeping their entire condition profile together."""
    frame = metadata.reset_index(drop=True)
    profiles = frame.groupby("donor_id", sort=False)["condition"].agg(lambda x: tuple(sorted(set(x.astype(str)))))
    blocks: list[int] = []
    for profile in dict.fromkeys(profiles):
        donors = profiles.index[profiles.map(lambda value: value == profile)].to_numpy()
        for donor in rng.choice(donors, size=len(donors), replace=True):
            blocks.extend(np.flatnonzero(frame["donor_id"].to_numpy() == donor))
    return np.asarray(blocks, dtype=int)
