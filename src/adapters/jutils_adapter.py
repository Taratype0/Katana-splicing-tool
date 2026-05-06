from __future__ import annotations

from pathlib import Path

import pandas as pd
import subprocess
import sys


def _safe_name(value: str) -> str:
    return "".join(ch if ch.isalnum() or ch in ("-", "_", ".") else "_" for ch in value).strip("_") or "comparison"


class JutilsAdapter:
    def __init__(self, script_path: Path) -> None:
        self.script_path = script_path

    def build_tsv_list(self, converted_root: Path, comparison_names: list[str], output_dir: Path) -> Path:
        output_dir.mkdir(parents=True, exist_ok=True)
        rows = []
        for comparison_name in comparison_names:
            rows.append(
                {
                    "tsv_file": str(converted_root / _safe_name(comparison_name) / "rmats_JC_results.tsv"),
                    "label": comparison_name,
                }
            )
        out = output_dir / "tsv_file_list.tsv"
        pd.DataFrame(rows).to_csv(out, sep="\t", index=False, header=False)
        return out

    def build_tsv_list_frame(self, comparison_names: list[str], converted_root: Path) -> pd.DataFrame:
        rows = []
        for comparison_name in comparison_names:
            rows.append(
                {
                    "comparison_name": comparison_name,
                    "converted_tsv": str(converted_root / _safe_name(comparison_name) / "rmats_JC_results.tsv"),
                }
            )
        return pd.DataFrame(rows)

    def build_meta_template(
        self,
        rmats_post_dir: Path,
        output_file: Path,
        control_label: str,
        experiment_label: str,
    ) -> Path:
        sample_file = None
        for candidate in rmats_post_dir.glob("*.MATS.JC.txt"):
            sample_file = candidate
            break
        if sample_file is None:
            raise FileNotFoundError(f"No JC file found in {rmats_post_dir}")

        df = pd.read_csv(sample_file, sep="\t", nrows=1, low_memory=False)
        if df.empty:
            raise ValueError(f"Unable to infer sample counts from {sample_file}")

        n1 = len(str(df.iloc[0]["IJC_SAMPLE_1"]).split(","))
        n2 = len(str(df.iloc[0]["IJC_SAMPLE_2"]).split(","))
        output_file.parent.mkdir(parents=True, exist_ok=True)
        with output_file.open("w", encoding="utf-8") as handle:
            for idx in range(1, n1 + 1):
                handle.write(f"{control_label}_rep_{idx}\t{control_label}\n")
            for idx in range(1, n2 + 1):
                handle.write(f"{experiment_label}_rep_{idx}\t{experiment_label}\n")
        return output_file

    def run_convert_results(self, rmats_post_dir: Path, output_dir: Path) -> subprocess.CompletedProcess:
        output_dir.mkdir(parents=True, exist_ok=True)
        return subprocess.run(
            [sys.executable, str(self.script_path), "convert-results", "--rmats-dir", str(rmats_post_dir), "--out-dir", str(output_dir)],
            check=True,
            capture_output=True,
            text=True,
        )

    def run_heatmap(self, tsv_file: Path, meta_file: Path, output_dir: Path, prefix: str) -> subprocess.CompletedProcess:
        output_dir.mkdir(parents=True, exist_ok=True)
        return subprocess.run(
            [
                sys.executable,
                str(self.script_path),
                "heatmap",
                "--tsv-file",
                str(tsv_file),
                "--meta-file",
                str(meta_file),
                "--p-value",
                "1.0",
                "--q-value",
                "0.05",
                "--dpsi",
                "0.20",
                "--prefix",
                prefix,
                "--out-dir",
                str(output_dir),
                "--pdf",
            ],
            check=True,
            capture_output=True,
            text=True,
        )

    def run_pca(self, tsv_file: Path, meta_file: Path, output_dir: Path, prefix: str) -> subprocess.CompletedProcess:
        output_dir.mkdir(parents=True, exist_ok=True)
        return subprocess.run(
            [
                sys.executable,
                str(self.script_path),
                "pca",
                "--tsv-file",
                str(tsv_file),
                "--meta-file",
                str(meta_file),
                "--p-value",
                "1.0",
                "--q-value",
                "0.05",
                "--dpsi",
                "0.20",
                "--prefix",
                prefix,
                "--out-dir",
                str(output_dir),
                "--pdf",
            ],
            check=True,
            capture_output=True,
            text=True,
        )

    def run_venn(self, tsv_list_file: Path, output_dir: Path, prefix: str) -> subprocess.CompletedProcess:
        output_dir.mkdir(parents=True, exist_ok=True)
        return subprocess.run(
            [
                sys.executable,
                str(self.script_path),
                "venn-diagram",
                "--tsv-file-list",
                str(tsv_list_file),
                "--prefix",
                prefix,
                "--out-dir",
                str(output_dir),
                "--pdf",
            ],
            check=True,
            capture_output=True,
            text=True,
        )
