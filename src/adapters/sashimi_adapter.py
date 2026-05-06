from __future__ import annotations

from pathlib import Path

import pandas as pd
import subprocess
import sys


class SashimiAdapter:
    def __init__(self, script_path: Path) -> None:
        self.script_path = script_path

    def run_manifest(self, manifest: pd.DataFrame, output_root: Path) -> tuple[str, pd.DataFrame]:
        logs: list[str] = []
        failures: list[dict[str, object]] = []
        for _, row in manifest.iterrows():
            outdir = Path(str(row["outdir"]))
            if not outdir.is_absolute():
                outdir = output_root / "06_sashimi" / "plots" / outdir
            outdir.mkdir(parents=True, exist_ok=True)
            b1 = self._bam_arg(Path(str(row["b1_txt"])))
            b2 = self._bam_arg(Path(str(row["b2_txt"])))
            cmd = [
                sys.executable,
                str(self.script_path),
                "-o",
                str(outdir),
                "--b1",
                b1,
                "--b2",
                b2,
                "--l1",
                str(row["label1"]),
                "--l2",
                str(row["label2"]),
                "--event-type",
                str(row["event_type"]),
                "-e",
                str(row["event_file"]),
            ]
            try:
                result = subprocess.run(cmd, check=True, capture_output=True, text=True)
                logs.append(result.stdout.strip())
            except subprocess.CalledProcessError as exc:
                failures.append(
                    {
                        "gene": str(row.get("geneSymbol", "")),
                        "event_type": str(row.get("event_type", "")),
                        "event_id": str(row.get("event_uid", row.get("event_file", ""))),
                        "comparison_id": str(row.get("comparison_id", "")),
                        "command": " ".join(cmd),
                        "script_path": str(self.script_path),
                        "output_folder": str(outdir),
                        "return_code": exc.returncode,
                        "stdout": exc.stdout or "",
                        "stderr": exc.stderr or "",
                        "error_message": str(exc),
                    }
                )
            except Exception as exc:
                failures.append(
                    {
                        "gene": str(row.get("geneSymbol", "")),
                        "event_type": str(row.get("event_type", "")),
                        "event_id": str(row.get("event_uid", row.get("event_file", ""))),
                        "comparison_id": str(row.get("comparison_id", "")),
                        "command": " ".join(cmd),
                        "script_path": str(self.script_path),
                        "output_folder": str(outdir),
                        "return_code": "",
                        "stdout": "",
                        "stderr": "",
                        "error_message": str(exc),
                    }
                )
        return "\n".join(line for line in logs if line), pd.DataFrame(failures)

    def _bam_arg(self, bam_list_file: Path) -> str:
        with bam_list_file.open("r", encoding="utf-8") as handle:
            lines = [line.strip() for line in handle if line.strip()]
        return ",".join(lines)
