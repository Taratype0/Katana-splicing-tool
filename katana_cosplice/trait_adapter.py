"""Design-aware adapter for already frozen scores, not network discovery.

Independent samples retain the historical three-row/panel production family.
Repeated units use explicit production donor-fixed omnibus models, whose BH
family includes nuisance terms. Unsupported designs fail rather than guessing.
"""
from dataclasses import dataclass

import numpy as np
import pandas as pd

from katana_cosplice.association import (
    _aligned_metadata, associate_program_scores, factorial_program_association,
)


@dataclass
class TraitResult:
    condition_tests: pd.DataFrame
    nuisance_tests: pd.DataFrame
    coefficients: pd.DataFrame
    audit: dict


def biological_unit_design(metadata, unit_column=None):
    """Use declared IDs only; never infer pairing from replicate-name suffixes."""
    candidates = [key for key in ('donor_id', 'biological_sample_id') if key in metadata]
    if unit_column:
        if unit_column not in metadata:
            raise ValueError(f'Biological unit column is absent: {unit_column}')
        candidates = list(dict.fromkeys([unit_column, *candidates]))
    for key in candidates:
        values = metadata[key]
        if values.isna().any() or values.astype(str).str.strip().isin({'', 'nan', 'None'}).any():
            raise ValueError(f'Declared biological identifiers must be complete: {key}')
    selected = candidates[0] if candidates else None
    if selected:
        # An alias must not group distinct declared donors into one biological ID.
        for other in candidates[1:]:
            if metadata.groupby(other)[selected].nunique().gt(1).any():
                raise ValueError(f'Conflicting biological identity columns: {selected}, {other}')
    units = metadata[selected].astype(str) if selected else pd.Series(metadata.index.astype(str), index=metadata.index)
    return {'unit_column': selected or 'sample_id', 'units': units,
            'repeated_units': bool(units.duplicated().any()),
            'n_independent_units': int(units.nunique()),
            'identity_status': 'declared' if selected else 'sample_id assumed independent; donor identity unavailable'}


def require_independent_association_units(samples):
    """Fail-fast guard proposed for the legacy unblocked production function."""
    design = biological_unit_design(samples)
    if design['repeated_units']:
        raise ValueError('Unblocked ANOVA/Welch requires independent biological units; repeated '
                         f"{design['unit_column']} detected. Use the donor-aware SC workflow or an "
                         'explicit estimable blocked frozen-score model; batch_column does not '
                         'block associate_program_scores.')
    return {key: value for key, value in design.items() if key != 'units'}


def safe_trait_association(scores, samples, *, factor='condition', unit_column=None, batch_column=None):
    metadata = _aligned_metadata(scores, samples).copy()
    if factor not in metadata:
        raise ValueError(f'Trait factor is absent: {factor}')
    if metadata[factor].isna().any() or metadata[factor].astype(str).str.strip().eq('').any():
        raise ValueError('Trait factor must be complete')
    if metadata[factor].nunique() != 2:
        raise ValueError('This bounded adapter supports one predeclared binary factor only')
    numeric = scores.apply(pd.to_numeric, errors='coerce')
    if not np.isfinite(numeric.to_numpy()).all():
        raise ValueError('Freeze finite scores first; no silent per-panel donor deletion')
    if numeric.var(axis=0).le(1e-12).any():
        raise ValueError('Constant program scores are non-testable')
    design = biological_unit_design(metadata, unit_column)
    audit = {key: value for key, value in design.items() if key != 'units'}
    audit.update({'factor': factor, 'network_discovery': 'not performed',
                  'scope': 'frozen-score test layer only; no post-discovery independent validation'})
    metadata = metadata.reset_index()
    blocked = design['repeated_units']
    nuisance_column = batch_column
    if blocked:
        unit = design['unit_column']
        if metadata.duplicated([unit, factor]).any():
            raise ValueError('Aggregate technical partitions: more than one score per biological unit-factor level')
        if not metadata.groupby(unit)[factor].nunique().eq(2).all():
            raise ValueError('Incomplete pairs are not supported; no silent dropping or guessed pairing')
        if design['n_independent_units'] < 3:
            raise ValueError('Blocked production procedure needs at least three complete independent pairs')
        if batch_column and batch_column != unit:
            raise ValueError('Additional batch with donor blocking is unsupported by this single-nuisance adapter; '
                             'use the existing SC donor+batch design or an explicit model')
        nuisance_column = unit
    if not blocked and not nuisance_column:
        table = associate_program_scores(numeric, metadata, [factor])
        audit.update({'procedure': 'production independent ANOVA plus two Welch rows',
                      'BH_family': 'all 3*K rows including duplicate two-sided Welch hypotheses',
                      'BH_family_size': len(table)})
        return TraitResult(table, pd.DataFrame(), pd.DataFrame(), audit)
    if nuisance_column not in metadata:
        raise ValueError(f'Nuisance column is absent: {nuisance_column}')
    if metadata[nuisance_column].isna().any() or metadata[nuisance_column].astype(str).str.strip().eq('').any():
        raise ValueError('Nuisance column must be complete')
    if nuisance_column == factor:
        raise ValueError('Trait and nuisance factor cannot be identical')
    dummies = pd.get_dummies(metadata[[factor, nuisance_column]].astype(str), drop_first=True, dtype=float)
    full_design = np.column_stack([np.ones(len(metadata)), dummies.to_numpy()])
    if np.linalg.matrix_rank(full_design) != full_design.shape[1]:
        raise ValueError('Trait/nuisance design is rank deficient or confounded')
    residual_df = len(metadata) - full_design.shape[1]
    if residual_df <= 0:
        raise ValueError('No positive residual degrees of freedom')
    all_terms, coefficients = factorial_program_association(numeric, metadata, [factor], nuisance_column)
    condition = all_terms.loc[all_terms.model_term.eq(f'C(Q("{factor}"), Sum)')].copy()
    nuisance = all_terms.loc[~all_terms.model_term.eq(f'C(Q("{factor}"), Sum)')].copy()
    if len(condition) != numeric.shape[1] or not np.isfinite(condition.pvalue).all():
        raise ValueError('Production blocked fit returned incomplete/nonfinite condition tests')
    audit.update({'procedure': 'production donor-fixed Type-III OLS' if blocked else 'production batch-adjusted Type-III OLS',
                  'nuisance_column': nuisance_column, 'residual_df': residual_df,
                  'BH_family': 'all condition AND nuisance omnibus terms; unchanged production family',
                  'BH_family_size': len(all_terms), 'unblocked_rows_emitted': False})
    return TraitResult(condition, nuisance, coefficients, audit)
