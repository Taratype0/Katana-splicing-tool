import copy
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from katana_next.science.motiftrace_custom import add_custom_results, fixed_pattern, load_custom_library


HEADER = "motif_id\tRBP\tmotif\tspecies\tevidence\tsource\n"


def library_file(tmp_path, word="CUCUYY", species="Mus_musculus"):
    path = tmp_path / "user.tsv"
    path.write_text(HEADER + f"test1\tPTBP1\t{word}\t{species}\tuser_defined\ttest fixture only\n")
    return path


def payload():
    rows = []
    for i in range(60):
        role = ("up", "dn", "bg")[i//20]
        word = "CUCUCU" if role == "dn" else "AAAAAA"
        forms = [{"region": region, "RNA": word * 5, "blocks": [[100, 130]], "junctions": [] if region == "cassette_exon" else [15]}
                 for region in ("cassette_exon", "inclusion_junctions", "skipping_junction")]
        rows.append(dict(event_id=f"synthetic_{i}", role=role, cassette_length=30, forms=forms))
    return dict(cases=["test"], events={"test": rows}, species="Mus_musculus", genome_build="TEST_NOT_REAL_DATA",
                family_size=100, model_count=0, source_count=0, models={}, sources={}, scalar={"test": {"public": [["keep", 1]]}},
                summary=[], backgrounds={"test": [.25]*4})


@pytest.mark.parametrize("word,pattern,width", [("CUCUYY", "CTCT[CT][CT]", 6), ("AU[CU]N", "AT[CT][ACGT]", 4), ("A", "A", 1)])
def test_safe_rna_iupac_words(word, pattern, width):
    assert fixed_pattern(word) == (pattern, width)


@pytest.mark.parametrize("word", ["A.*G", "A{1,4}", "(?=A)", "A|C", "X", "[AG", "", "A"*26])
def test_unsafe_or_unbounded_regex_rejected(word):
    with pytest.raises(ValueError):
        fixed_pattern(word)


def test_species_source_and_duplicate_ids_are_not_silently_fixed(tmp_path):
    path = library_file(tmp_path, species="Homo_sapiens")
    with pytest.raises(ValueError, match="species"):
        load_custom_library(path, species="Mus_musculus")
    path = library_file(tmp_path)
    path.write_text(path.read_text() + path.read_text().splitlines()[1] + "\n")
    with pytest.raises(ValueError, match="duplicate"):
        load_custom_library(path, species="Mus_musculus")


def test_dedup_preserves_all_user_source_mappings(tmp_path):
    path = library_file(tmp_path)
    path.write_text(path.read_text() + "test2\tOtherRBP\tCTCT[CT][CT]\tMus_musculus\tunknown\tsecond source\n")
    result = load_custom_library(path, species="Mus_musculus")
    assert result["unique_models"] == 1 and result["source_rows"] == 2
    assert "OtherRBP" in next(iter(result["models"].values()))["RBPs"]
    assert all(row["evidence"].startswith("user_declared_") for row in next(iter(result["sources"].values())))


def test_preview_does_not_create_inferential_results_or_mutate_public_family(tmp_path):
    library = load_custom_library(library_file(tmp_path), species="Mus_musculus")
    data = payload()
    original = copy.deepcopy(data["scalar"])
    result = add_custom_results(data, library, mode="preview_only", output=tmp_path / "preview")
    assert result["scalar"] == original and result["family_size"] == 100
    assert result["custom"]["family_size"] is None
    assert not (tmp_path / "preview/CUSTOM_RESULTS.json").exists()
    with pytest.raises(FileExistsError):
        add_custom_results(data, library, mode="preview_only", output=tmp_path / "preview")


def test_explicit_custom_inference_keeps_originals_and_short_NA(tmp_path):
    library = load_custom_library(library_file(tmp_path), species="Mus_musculus")
    data = payload()
    result = add_custom_results(data, library, mode="fixed_enrichment", output=tmp_path / "run")
    table = pd.read_parquet(tmp_path / "run/CUSTOM_ALL_TESTS_WITH_NA.parquet")
    assert data["scalar"]["test"]["public"] == [["keep", 1]]
    assert data["family_size"] == 100 and result["custom"]["family_size"] > 0
    assert table[table.scope.eq("short_exon_1_9")].raw_P.isna().all()
    assert "full_exon_sliding" in set(table.scope)
    assert table.raw_P.isna().equals(table.q_custom_library_family.isna())
    assert json.loads((tmp_path / "run/CUSTOM_ANALYSIS_PLAN.json").read_text())["selection_posthoc_possible"] is True


def test_minimal_meme_requires_explicit_metadata_and_native_probabilities(tmp_path):
    path = tmp_path / "user.meme"
    path.write_text("MEME version 4\nALPHABET= ACGU\nMOTIF native1 Example\nletter-probability matrix: alength= 4 w= 3\n.1 .2 .3 .4\n1 0 0 0\n0 1 0 0\n")
    sidecar = Path(str(path) + ".metadata.tsv")
    with pytest.raises(ValueError):
        load_custom_library(path, species="Mus_musculus")
    sidecar.write_text("motif_id\tRBP\tspecies\tevidence\tsource\nnative1\tExample\tMus_musculus\tunknown\ttest only\n")
    result = load_custom_library(path, species="Mus_musculus")
    assert len(result["inputs"]) == 2
    assert np.allclose(next(iter(result["models"].values()))["matrix"][0], [.1, .2, .3, .4])
    path.write_text(path.read_text().replace(".1 .2 .3 .4", "1 2 3 4"))
    with pytest.raises(ValueError, match="probabilities"):
        load_custom_library(path, species="Mus_musculus")


def test_normal_ui_names_and_legacy_ids_stay_separate():
    from katana_next.workflows.registry import WorkflowRegistry
    registry = WorkflowRegistry()
    # Retain the accepted scientific identity and the independently approved
    # researcher-facing copy; packaging does not rename either layer.
    labels = {
        "prioritize.rbp_mechanism": "RBP sequence evidence: regulatory clues / consequences",
        "prioritize.rbp_sequence_results": "Katana MotifTrace: sequence workbench",
        "prioritize.se_rbp_impact": "Katana MotifTrace: legacy SE fixed-site enrichment",
    }
    for key, label in labels.items():
        assert registry.get(key).definition.label.startswith("Katana MotifTrace")
        assert registry.get(key).definition.to_dict()["label"] == label
    main = registry.get("prioritize.rbp_mechanism").normalize_parameters({})
    assert main["sequence_scope"] == "boundary_association"
    assert main["custom_motif_mode"] == "preview_only"


def test_new_legacy_display_names_leave_statistical_files_untouched(tmp_path):
    from katana_next.workflows.se_rbp_impact import _brand_new_report
    report = tmp_path / "REPORT.html"
    report.write_text("<title>Katana SE RBP real evidence</title><h1>Katana：SE 序列影响与 RBP motif</h1><p>P=0.123</p>", encoding="utf-8")
    (tmp_path / "CONCLUSION.md").write_text("# comparison SE RBP\n\nOriginal negative.\n", encoding="utf-8")
    counts = tmp_path / "enrichment.tsv"
    counts.write_bytes(b"raw_P\tq\n0.123\t0.456\n")
    _brand_new_report(tmp_path)
    assert "Katana MotifTrace" in report.read_text(encoding="utf-8")
    assert "P=0.123" in report.read_text(encoding="utf-8")
    assert counts.read_bytes() == b"raw_P\tq\n0.123\t0.456\n"
