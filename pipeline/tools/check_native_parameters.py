#!/usr/bin/env python3
"""Exercise flat options with real Nextflow and explicit transport-only tasks."""
import argparse
import csv
import json
import os
from pathlib import Path
import subprocess

PIPELINE = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--nextflow-command-json", type=Path)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    prefix = json.loads(args.nextflow_command_json.read_text()) if args.nextflow_command_json else ["nextflow"]
    root = args.outdir.resolve()
    root.mkdir(parents=True, exist_ok=False)
    inputs = root / "inputs"
    inputs.mkdir()
    (inputs / "project.json").write_text("{}\n")
    results = {}
    def execute(label, extra, success=True):
        command = [*prefix, "run", str(PIPELINE / "main.nf"), *extra,
                   "--flow_executable", str(PIPELINE / "tests/protocol_fixture.py"),
                   "-c", str(PIPELINE / "tests/transport.config"),
                   "-ansi-log", "false", "--outdir", str(root / "output"),
                   "-work-dir", str(root / "work")]
        result = subprocess.run(command, cwd=root, env={**os.environ, "NXF_OFFLINE": "true"},
                                text=True, capture_output=True, timeout=180)
        (root / (label + ".stdout.txt")).write_text(result.stdout)
        (root / (label + ".stderr.txt")).write_text(result.stderr)
        assert (result.returncode == 0) == success, result.stdout[-6000:] + result.stderr[-2000:]
        traces = list((root / "output/pipeline_info").glob("*trace*"))
        # Use the trace path selected by the workflow's own config.
        traces = [p for p in traces if p.is_file()]
        latest = max(traces, key=lambda p: p.stat().st_mtime_ns)
        rows = list(csv.DictReader(latest.open(), delimiter="\t"))
        (root / (label + ".trace.json")).write_text(json.dumps(rows, indent=2))
        results[label] = {"exit_code": result.returncode, "tasks": len(rows),
                          "cached": sum(r["status"] == "CACHED" for r in rows)}
        return result, rows
    base = ["--mode", "downstream-only", "--project", str(inputs / "project.json"),
            "--skip_network", "--skip_mapp", "--skip_motiftrace", "--skip_motif_binding",
            "--skip_enrichment", "--skip_rbp_identity", "--deg_fdr", "0.01"]
    _, rows = execute("first", base)
    assert rows and all(r["status"] == "COMPLETED" for r in rows)
    assert any(r["name"].endswith("KATANA_FLOW_CONFIGURE") for r in rows)
    config = root / "output/pipeline_info/configuration/flow.json"
    assert json.loads(config.read_text())["parameters"]["differential.deg"]["fdr"] == .01
    params_file = root / "parameters.yaml"
    params_file.write_text(json.dumps({"mode": "downstream-only", "project": str(inputs / "project.json"),
        "skip_network": True, "skip_mapp": True, "skip_motiftrace": True, "skip_motif_binding": True,
        "skip_enrichment": True, "skip_rbp_identity": True, "deg_fdr": .01}, indent=2))
    _, resumed = execute("resume", ["-params-file", str(params_file), "-resume"])
    assert resumed and all(r["status"] == "CACHED" for r in resumed)
    changed = list(base)
    changed[changed.index("0.01")] = "0.02"
    _, updated = execute("changed-threshold", [*changed, "-resume"])
    assert any(r["name"].endswith("KATANA_FLOW_CONFIGURE") and r["status"] == "COMPLETED" for r in updated)
    assert json.loads(config.read_text())["parameters"]["differential.deg"]["fdr"] == .02
    # Actual raw main.nf wiring: missing optional-analysis resources stop the
    # generated preflight configuration before any read-processing task starts.
    for name, content in (("samples.csv", "sample,condition,strandedness\na,A,reverse\n"),
                          ("contrasts.csv", "contrast,treatment,control\nB_vs_A,B,A\n"),
                          ("genome.fa", ">test\nACGT\n"), ("genes.gtf", "# fixture\n")):
        (inputs / name).write_text(content)
    failed, raw_rows = execute("missing-resources", ["--input", str(inputs / "samples.csv"),
        "--contrasts", str(inputs / "contrasts.csv"), "--fasta", str(inputs / "genome.fa"),
        "--gtf", str(inputs / "genes.gtf"), "--species", "Mus musculus"], success=False)
    assert "--gene_sets" in failed.stdout + failed.stderr
    assert raw_rows and all(r["name"].endswith("KATANA_FLOW_CONFIGURE") for r in raw_rows)
    (root / "acceptance.json").write_text(json.dumps({"status": "passed", "checks": results,
        "scope": "Native Nextflow configuration generation, task handoff, early rejection and resume; no biological calculations or container execution."}, indent=2))
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
