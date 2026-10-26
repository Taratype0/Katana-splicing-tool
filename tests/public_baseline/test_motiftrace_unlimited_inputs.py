"""Exercise removed input-size gates without allocating production-sized matrices."""
from itertools import islice, product
from pathlib import Path

import pandas as pd
import pytest

from katana_next.science import motiftrace_custom as custom
from katana_next.science.motiftrace_events import extract_event_rna, trace_instances, validate_models


def _rna_fixture(tmp_path):
    fasta = tmp_path / "generated.fa"
    sequence = "ACGT" * 100
    fasta.write_text(">chr1\n" + sequence + "\n", encoding="ascii", newline="\n")
    Path(str(fasta) + ".fai").write_text("chr1\t400\t6\t400\t401\n", newline="\n")
    events = pd.DataFrame([dict(event_id="event", gene_id="gene", comparison_id="B_vs_A",
        event_type="SE", strand="+", chr="chr1", dpsi=.3, fdr=.001,
        upstreamES=0, upstreamEE=40, exonStart_0base=100, exonEnd=135,
        downstreamES=300, downstreamEE=340, IncLevel1=".6,.7", IncLevel2=".3,.4")])
    return events, fasta


class _CounterOnlyModels(dict):
    """Large declared scheduling count, one iterable model: no giant science run."""
    def __len__(self):
        return 1_000_001


def test_legacy_base_limit_is_accepted_without_dropping_any_form(tmp_path):
    events, fasta = _rna_fixture(tmp_path)
    expected = extract_event_rna(events, fasta)
    assert expected.length.sum() > 10
    actual = extract_event_rna(events, fasta, max_total_bases=10)
    pd.testing.assert_frame_equal(actual, expected)
    assert set(actual.region) >= {"inclusion", "exclusion", "affected_inclusion"}


def test_more_than_256_valid_models_are_all_validated():
    models = {f"model_{i}": {"kind": "consensus", "pattern": "AC", "width": 2} for i in range(257)}
    validate_models(models)
    models["model_256"]["width"] = 3
    with pytest.raises(ValueError, match="width disagrees"):
        validate_models(models)


def test_large_declared_motif_sequence_work_does_not_block_locus_tracing(tmp_path):
    events, fasta = _rna_fixture(tmp_path)
    rna = extract_event_rna(events, fasta)
    ordinary = {"test": {"kind": "consensus", "pattern": "AC", "width": 2}}
    counters = _CounterOnlyModels(ordinary)
    assert rna.loc[rna.layer.eq("mature_local_form"), "length"].sum() * len(counters) > 25_000_000
    expected = trace_instances(rna, events, ordinary)
    actual = trace_instances(rna, events, counters)
    for first, second in zip(actual, expected):
        pd.testing.assert_frame_equal(first, second)


def test_large_declared_local_library_file_is_accepted(tmp_path, monkeypatch):
    path = tmp_path / "small_input_with_size_counter.tsv"
    path.write_text("fixture only\n")
    original_stat = Path.stat
    class _FileSizeCounter:
        def __init__(self, actual):
            self.actual = actual
            self.st_size = 16 * 1024**2 + 1
        def __getattr__(self, name):
            return getattr(self.actual, name)
    def counted_stat(self, *args, **kwargs):
        result = original_stat(self, *args, **kwargs)
        return _FileSizeCounter(result) if self == path else result
    monkeypatch.setattr(Path, "stat", counted_stat)
    assert custom._file(path) == path.resolve(strict=True)


def test_more_than_1024_sources_and_256_unique_models_preserve_every_mapping(tmp_path):
    path = tmp_path / "custom.tsv"
    words = ["".join(value) for value in islice(product("ACGT", repeat=6), 1025)]
    header = "motif_id\tRBP\tmotif\tspecies\tevidence\tsource\n"
    path.write_text(header + "".join(
        f"native_{i}\tRBP_{i}\t{word}\tMus_musculus\tuser_defined\tgenerated fixture\n"
        for i, word in enumerate(words)))
    library = custom.load_custom_library(path, species="Mus_musculus")
    assert library["source_rows"] == library["unique_models"] == len(words)
    assert len(library["models"]) == len(library["sources"]) == 1025
    assert {r["native_id"] for rows in library["sources"].values() for r in rows} == {
        f"native_{i}" for i in range(1025)}
    validate_models(library["models"])


def test_custom_test_count_counter_does_not_reject_complete_small_inference(tmp_path):
    path = tmp_path / "custom.tsv"
    path.write_text("motif_id\tRBP\tmotif\tspecies\tevidence\tsource\n"
                    "native\tRBP\tAC\tMus_musculus\tuser_defined\tgenerated fixture\n")
    library = custom.load_custom_library(path, species="Mus_musculus")
    library["models"] = _CounterOnlyModels(library["models"])
    rows = []
    for i in range(60):
        rows.append({"event_id": f"event_{i}", "role": ("up", "dn", "bg")[i // 20],
            "cassette_length": 30, "forms": [
                {"region": region, "RNA": "ACGU" * 7 + "AC", "blocks": [[0, 30]],
                 "junctions": [] if region == "cassette_exon" else [15]}
                for region in ("cassette_exon", "inclusion_junctions", "skipping_junction")]})
    payload = dict(cases=["B_vs_A"], events={"B_vs_A": rows}, species="Mus_musculus",
        genome_build="GENERATED_FIXTURE", family_size=123, model_count=0, source_count=0,
        models={}, sources={}, scalar={"B_vs_A": {}}, summary=[], backgrounds={"B_vs_A": [.25] * 4})
    output = tmp_path / "run"
    result = custom.add_custom_results(payload, library, mode="fixed_enrichment", output=output)
    table = pd.read_parquet(output / "CUSTOM_ALL_TESTS_WITH_NA.parquet")
    assert not table.empty and set(table.model_id) == set(library["models"])
    assert set(table.window_nt) >= {0, 10, 20}
    assert payload["family_size"] == 123
    assert result["custom"]["family_size"] == int(table.raw_P.notna().sum())
