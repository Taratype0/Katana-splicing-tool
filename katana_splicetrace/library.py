"""Curated v3.0 (170 RBPs: Splicing + Epitranscriptome m6A + ENCODE RBNS/eCLIP) & Extensible Motif Database Loader."""
from __future__ import annotations

import gzip
from importlib import resources
import io
import json
from pathlib import Path
import pandas as pd

from katana_splicetrace.extensible_db import (
    CURATED_EXPANSION_ENTRIES_V3,
    NON_SPLICING_CONTAMINANTS_V3 as NON_SPLICING_CONTAMINANTS,
    CuratedMotifLibrary,
    build_curated_library_v3,
    compute_pwm_ic_bits,
    export_katana_motif_pack,
    get_user_motif_packs_dir,
    human_to_mouse_rbp_symbol,
    iupac_to_pwm,
    load_external_motif_entries,
    merge_motif_entries_into_library,
    parse_homer_motifs,
    parse_iupac_table,
    parse_jaspar_pfm,
    parse_meme_motifs,
    pwm_consensus,
    validate_pwm_matrix,
)


def load_curated_human_rbp_library(
    models_json: Path,
    assignments_tsv: Path,
    *,
    include_v3_expansion: bool = True,
    extra_motif_files: list[Path | str] | None = None,
    include_installed_user_packs: bool = True,
    species: str = "human",
    min_width: int = 4,
    max_width: int = 15,
    min_ic_bits: float = 4.0,
    max_intra_rbp_similarity: float = 0.92,
) -> CuratedMotifLibrary:
    raw_models = json.loads(Path(models_json).read_text(encoding="utf-8"))
    assign = pd.read_csv(assignments_tsv, sep="\t", dtype=str)
    return build_curated_library_v3(
        raw_models,
        assign,
        include_v3_expansion=include_v3_expansion,
        extra_motif_files=extra_motif_files,
        include_installed_user_packs=include_installed_user_packs,
        species=species,
        min_width=min_width,
        max_width=max_width,
        min_ic_bits=min_ic_bits,
        max_intra_rbp_similarity=max_intra_rbp_similarity,
    )


def load_bundled_human_rbp_library(
    *,
    include_v3_expansion: bool = True,
    extra_motif_files: list[Path | str] | None = None,
    include_installed_user_packs: bool = True,
    species: str = "human",
    min_width: int = 4,
    max_width: int = 15,
    min_ic_bits: float = 4.0,
    max_intra_rbp_similarity: float = 0.92,
) -> CuratedMotifLibrary:
    """Load the bundled 170-RBP Splicing & Epitranscriptome Library v3.0 (+ any installed/extra motif packs)."""
    data_pkg = resources.files("katana_splicetrace.data")
    with gzip.open(io.BytesIO((data_pkg / "MODELS.json.gz").read_bytes()), "rt", encoding="utf-8") as fm:
        raw_models = json.load(fm)
    with gzip.open(io.BytesIO((data_pkg / "ASSIGNMENTS.tsv.gz").read_bytes()), "rt", encoding="utf-8") as fa:
        assign = pd.read_csv(fa, sep="\t", dtype=str)
    return build_curated_library_v3(
        raw_models,
        assign,
        include_v3_expansion=include_v3_expansion,
        extra_motif_files=extra_motif_files,
        include_installed_user_packs=include_installed_user_packs,
        species=species,
        min_width=min_width,
        max_width=max_width,
        min_ic_bits=min_ic_bits,
        max_intra_rbp_similarity=max_intra_rbp_similarity,
    )
