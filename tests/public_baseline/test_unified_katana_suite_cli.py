"""End-to-end verification of the dual-entrypoint & converging `katana` Python suite (`run-bulk`, `run-sc`, `--bundle`)."""
from __future__ import annotations

import json
from pathlib import Path
import zipfile

import numpy as np
import pandas as pd

from katana import __version__, run_bulk, run_single_cell
from katana.cli import main as katana_cli_main


def _write_public_rmats_fixture(rdir: Path, n_events: int = 55) -> None:
    rdir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(42)
    genes = [
        "PTBP1", "RBFOX2", "SRSF1", "SRSF3", "HNRNPA1", "HNRNPK", "U2AF2", "TRA2A",
        "SNRNP70", "MAPT", "PTPRC", "NUMB", "FLNA", "PKM", "CD44", "FGFR2", "TPM1", "ACTN1",
    ]
    rows = []
    for idx in range(n_events):
        gs = genes[idx % len(genes)]
        p1 = np.clip(rng.normal(0.72 if idx < 14 else 0.45, 0.05, size=3), 0.05, 0.95)
        p2 = np.clip(rng.normal(0.28 if idx < 14 else 0.46, 0.05, size=3), 0.05, 0.95)
        rows.append({
            "ID": idx + 1,
            "GeneID": f"ENSG00000{100000 + idx}",
            "geneSymbol": gs,
            "chr": "chr1",
            "strand": "+",
            "exonStart_0base": 10000 + idx * 500,
            "exonEnd": 10120 + idx * 500,
            "upstreamES": 9500 + idx * 500,
            "upstreamEE": 9600 + idx * 500,
            "downstreamES": 10400 + idx * 500,
            "downstreamEE": 10500 + idx * 500,
            "IJC_SAMPLE_1": "24,22,26" if idx < 45 else "2,1,2",
            "SJC_SAMPLE_1": "8,9,7" if idx < 45 else "1,1,0",
            "IJC_SAMPLE_2": "8,7,9",
            "SJC_SAMPLE_2": "22,25,21",
            "IncFormLen": 148,
            "SkipFormLen": 74,
            "PValue": 1e-6 if idx < 14 else 0.65,
            "FDR": 0.002 if idx < 14 else 0.80,
            "IncLevel1": ",".join(f"{x:.3f}" for x in p1),
            "IncLevel2": ",".join(f"{x:.3f}" for x in p2),
            "IncLevelDifference": round(float(p1.mean() - p2.mean()), 4),
        })
    pd.DataFrame(rows).to_csv(rdir / "SE.MATS.JC.txt", sep="\t", index=False)


def test_standalone_bulk_cli_and_nextflow_bundle_convergence(tmp_path: Path):
    assert __version__ == "2.7.0"

    # 1. Route A: Standalone Python package from rMATS + samples + MultiQC
    rmats_dir = tmp_path / "rmats_input" / "KD_vs_Ctrl"
    _write_public_rmats_fixture(rmats_dir)
    samples_tsv = tmp_path / "samples.tsv"
    pd.DataFrame([
        {"sample_id": "KD_1", "condition": "KD", "biological_sample_id": "KD_1"},
        {"sample_id": "KD_2", "condition": "KD", "biological_sample_id": "KD_2"},
        {"sample_id": "KD_3", "condition": "KD", "biological_sample_id": "KD_3"},
        {"sample_id": "Ctrl_1", "condition": "Ctrl", "biological_sample_id": "Ctrl_1"},
        {"sample_id": "Ctrl_2", "condition": "Ctrl", "biological_sample_id": "Ctrl_2"},
        {"sample_id": "Ctrl_3", "condition": "Ctrl", "biological_sample_id": "Ctrl_3"},
    ]).to_csv(samples_tsv, sep="\t", index=False)

    mqc_dir = tmp_path / "multiqc_data"
    mqc_dir.mkdir()
    mqc_json = {
        "report_plot_data": {
            "rseqc_gene_body_coverage_plot": {
                "datasets": [[
                    {
                        "name": s,
                        "data": [[pct, float(80 + 20 * np.sin(pct / 100 * np.pi))] for pct in range(1, 101)],
                    }
                    for s in ("KD_1", "KD_2", "KD_3", "Ctrl_1", "Ctrl_2", "Ctrl_3")
                ]]
            }
        }
    }
    (mqc_dir / "multiqc_data.json").write_text(json.dumps(mqc_json), encoding="utf-8")

    out_standalone = tmp_path / "out_standalone_bulk"
    rc = katana_cli_main([
        "run-bulk",
        "--rmats", str(tmp_path / "rmats_input"),
        "--samples", str(samples_tsv),
        "--multiqc", str(mqc_dir),
        "--study-name", "Public ENCODE RBP Demo",
        "--outdir", str(out_standalone),
    ])
    assert rc == 0
    assert (out_standalone / "REPORT.html").is_file()
    assert (out_standalone / "REPORT.zip").is_file()
    assert (out_standalone / "outputs/landscape/psi-heatmap.png").is_file()
    assert (out_standalone / "outputs/landscape/psi-pca.png").is_file()
    assert (out_standalone / "figures/qc/multiqc-5to3-coverage.png").is_file()
    with zipfile.ZipFile(out_standalone / "REPORT.zip") as zf:
        assert "REPORT.html" in zf.namelist()
        assert "assets/report.css" in zf.namelist()

    # 2. Route B Convergence: Completed Nextflow upstream bundle (`--bundle`)
    bundle_dir = tmp_path / "nextflow_bundle"
    (bundle_dir / "rmats" / "KD_vs_Ctrl").mkdir(parents=True)
    _write_public_rmats_fixture(bundle_dir / "rmats" / "KD_vs_Ctrl")
    pd.read_csv(samples_tsv, sep="\t").to_csv(bundle_dir / "samplesheet.csv", index=False)
    pd.DataFrame([{"contrast": "KD_vs_Ctrl", "treatment": "KD", "control": "Ctrl"}]).to_csv(bundle_dir / "contrasts.csv", index=False)

    out_converged = tmp_path / "out_converged_bundle"
    receipt = run_bulk(bundle=bundle_dir, outdir=out_converged, study_name="Converged Nextflow Bundle Demo")
    assert receipt["status"] == "completed"
    assert Path(receipt["report_zip"]).is_file()


def test_standalone_single_cell_suite_and_event_psi_umap(tmp_path: Path):
    rng = np.random.default_rng(7)
    cells = []
    for idx in range(28):
        cond = "Activated" if idx < 14 else "Resting"
        ct = "Th1" if idx % 2 == 0 else "Treg"
        donor = f"Donor_{(idx % 4) + 1}"
        cells.append({"cell_id": f"Cell_{idx+1}", "condition": cond, "donor_id": donor, "cell_type": ct})
    cells_df = pd.DataFrame(cells)
    cells_tsv = tmp_path / "sc_cells.tsv"
    cells_df.to_csv(cells_tsv, sep="\t", index=False)

    cids = cells_df["cell_id"].tolist()
    eids = [f"PTPRC|SE|{i}" for i in range(30)]
    psi_mat = pd.DataFrame(rng.uniform(0.1, 0.9, size=(30, len(cids))), index=eids, columns=cids)
    sup_mat = pd.DataFrame(rng.integers(5, 35, size=(30, len(cids))), index=eids, columns=cids)
    expr_mat = pd.DataFrame(rng.normal(3.0, 1.0, size=(40, len(cids))), index=[f"GENE_{i}" for i in range(40)], columns=cids)
    psi_tsv = tmp_path / "sc_psi.tsv"
    sup_tsv = tmp_path / "sc_sup.tsv"
    expr_tsv = tmp_path / "sc_expr.tsv"
    psi_mat.to_csv(psi_tsv, sep="\t")
    sup_mat.to_csv(sup_tsv, sep="\t")
    expr_mat.to_csv(expr_tsv, sep="\t")

    ev_df = pd.DataFrame([
        {
            "event_id": eids[i],
            "gene_symbol": ["PTPRC", "SRSF1", "TRA2A", "HNRPLL", "CD44", "MAP4K1"][i % 6],
            "comparison_id": "Activated_vs_Resting",
            "dpsi": 0.28 if i < 10 else 0.03,
            "fdr": 0.004 if i < 10 else 0.72,
        }
        for i in range(30)
    ])
    ev_tsv = tmp_path / "sc_events.tsv"
    ev_df.to_csv(ev_tsv, sep="\t", index=False)

    out_sc = tmp_path / "out_sc"
    receipt = run_single_cell(
        outdir=out_sc,
        events=ev_tsv,
        cells=cells_tsv,
        psi_matrix=psi_tsv,
        support_matrix=sup_tsv,
        expression=expr_tsv,
        study_name="Public Smart-seq3 T-Cell Demo",
    )
    assert receipt["status"] == "completed"
    assert (out_sc / "REPORT.html").is_file()
    assert (out_sc / "REPORT.zip").is_file()
    assert (out_sc / "outputs/landscape/sc-dual-psi-umap.png").is_file()
