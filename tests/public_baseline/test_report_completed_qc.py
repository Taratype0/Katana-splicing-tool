"""QC imports preserve native units and cannot silently claim absent gene-body QC."""
from hashlib import sha256
import json
from pathlib import Path
import zipfile

import pytest

from katana_next.reports.qc import build_qc_assets, parse_fastqc, parse_star_log, parse_trim_galore


FASTQC = """##FastQC\t0.12.1
>>Basic Statistics\tpass
#Measure\tValue
Filename\tSample_1.fastq.gz
Total Sequences\t1000
Sequence length\t150
%GC\t50
>>END_MODULE
>>Per base sequence quality\twarn
#Base\tMean\tMedian
1\t39.123456789\t40
10-14\t38.999999999\t40
>>END_MODULE
>>Sequence Length Distribution\tpass
#Length\tCount
150\t1000
>>END_MODULE
"""


def _fastqc(path):
    with zipfile.ZipFile(path, "w") as archive:
        archive.writestr("Sample_1_fastqc/fastqc_data.txt", FASTQC)


def test_native_fastqc_status_cells_and_units_preserved(tmp_path):
    path = tmp_path / "Sample_1_fastqc.zip"
    _fastqc(path)
    result = parse_fastqc(path)
    assert result["version"] == "0.12.1"
    assert result["basic"]["Total Sequences"] == "1000"
    assert result["basic"]["Sequence length"] == "150"
    quality = result["modules"]["Per base sequence quality"]
    assert quality["status"] == "warn"
    assert quality["rows"] == [["1", "39.123456789", "40"], ["10-14", "38.999999999", "40"]]


def test_trim_summary_does_not_claim_final_paired_retention(tmp_path):
    path = tmp_path / "Sample_1.fastq.gz_trimming_report.txt"
    path.write_text("Total reads processed: 1,000\nReads with adapters: 250 (25.0%)\n"
                    "Reads written (passing filters): 1,000 (100.0%)\n"
                    "Quality-trimmed: 7,000 bp (1%)\n", encoding="utf-8")
    result = parse_trim_galore(path)
    assert result["reads_processed"] == 1000
    assert result["reads_with_adapters"] == 250
    assert result["cutadapt_reads_written"] == 1000
    assert "not final paired-read retention" in result["count_scope"]
    assert result["quality_cutoff"] is None


def test_native_star_percentages_are_not_rescaled(tmp_path):
    path = tmp_path / "Log.final.out"
    path.write_text("Number of input reads | 1000\nUniquely mapped reads % | 79.66%\n"
                    "Average input read length | 297\n", encoding="utf-8")
    result = parse_star_log(path)
    assert result["Number of input reads"] == "1000"
    assert result["Uniquely mapped reads %"] == "79.66%"
    assert result["Average input read length"] == "297"


def test_build_qc_outputs_without_reading_bam_or_rerunning_tools(tmp_path):
    root = tmp_path / "original-qc"
    root.mkdir()
    path = root / "Sample_1_fastqc.zip"
    _fastqc(path)
    html = root / "Sample_1_fastqc.html"
    html.write_text("<html>Complete original QC</html>")
    before = {p.name: (p.stat().st_mtime_ns, p.read_bytes()) for p in root.iterdir()}
    out = tmp_path / "qc-report"
    result = build_qc_assets(root, out)
    assert result["receipt"]["passed"]
    assert result["receipt"]["upstream_jobs"] == 0 and result["receipt"]["BAM_files"] == 0
    assert result["receipt"]["fastqc_files"] == 1
    assert result["receipt"]["original_html_reports"] == 1
    assert result["receipt"]["figures"] >= 4
    assert result["unavailable"][0]["status"] == "not_computed"
    assert "RSeQC" in result["unavailable"][0]["reason"]
    assert {p.name: (p.stat().st_mtime_ns, p.read_bytes()) for p in root.iterdir()} == before
    raw = (out / "tables/fastqc_all_native_module_values.tsv").read_text()
    assert "39.123456789" in raw and "38.999999999" in raw
    assert all(row["analysis_id"] == "quality.reads" for row in result["artifacts"])
    figures = {Path(a['path']).stem: a for a in result['artifacts'] if a['path'].endswith('.png')}
    assert all({l['role'] for l in a['metadata']['data_links']} == {'plot_data', 'full_results', 'methods'} for a in figures.values())
    assert all(Path(l['path']).is_file() for a in figures.values() for l in a['metadata']['data_links'])
    import pandas as pd
    def plot(stem):
        return pd.read_csv(next(l['path'] for l in figures[stem]['metadata']['data_links'] if l['role'] == 'plot_data'), sep='\t')
    assert plot('fastqc_read_counts').reads_millions.tolist() == [.001]
    quality = plot('fastqc_base_quality')
    assert quality.x_bin_midpoint.tolist() == [1., 12.]
    assert quality.plotted_y.tolist() == pytest.approx([39.123456789, 38.999999999])
    assert quality.native_y.tolist() == quality.plotted_y.tolist()
    length = plot('fastqc_read_lengths')
    assert length.native_y.tolist() == [1000.]
    assert length.plotted_y.tolist() == [100.]


def test_transfer_hash_mismatch_is_rejected_before_rendering(tmp_path):
    root = tmp_path / "original"
    root.mkdir()
    path = root / "Sample_1_fastqc.zip"
    _fastqc(path)
    (root / "QC_FETCH_MANIFEST.json").write_text(json.dumps({"files": [{
        "relative_path": path.name, "recorded_size": path.stat().st_size,
        "sha256": "wrong", "sample_id": "Sample"}]}))
    with pytest.raises(ValueError, match="verification failed"):
        build_qc_assets(root, tmp_path / "report")
    assert not (tmp_path / "report").exists()
