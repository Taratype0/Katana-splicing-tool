"""Parse and render existing FastQC, Trim Galore and STAR QC without reruns."""
from __future__ import annotations

from datetime import datetime, timezone
from hashlib import sha256
import json
from pathlib import Path
import re
import zipfile

VERSION = "completed-qc-report-v1"
BLUE, RED = "#2863A5", "#C63F4A"


def parse_fastqc(path: str | Path) -> dict:
    """Preserve native module status, raw table cells and basic-statistic units."""
    path = Path(path)
    with zipfile.ZipFile(path) as archive:
        names = [name for name in archive.namelist() if name.endswith("/fastqc_data.txt") or name == "fastqc_data.txt"]
        if len(names) != 1:
            raise ValueError("Expected exactly one fastqc_data.txt in " + path.name)
        if archive.getinfo(names[0]).file_size > 32 * 1024 * 1024:
            raise ValueError("FastQC text exceeds the explicit report parser size limit")
        text = archive.read(names[0]).decode("utf-8-sig")
    modules, current, version = {}, None, "not_recorded"
    for line in text.splitlines():
        if line.startswith("##FastQC"):
            version = line.split("\t", 1)[1]
        elif line.startswith(">>END_MODULE"):
            current = None
        elif line.startswith(">>"):
            name, status = line[2:].split("\t", 1)
            current = {"status": status, "columns": [], "rows": [], "comments": []}
            modules[name] = current
        elif current is not None and line.startswith("#"):
            current["comments"].append(line[1:])
            current["columns"] = line[1:].split("\t")
        elif current is not None and line.strip():
            current["rows"].append(line.split("\t"))
    basic = {row[0]: row[1] for row in modules.get("Basic Statistics", {}).get("rows", []) if len(row) == 2}
    return {"version": version, "basic": basic, "modules": modules}


def _number(value):
    try:
        return float(str(value).replace(",", "").rstrip("%"))
    except (TypeError, ValueError):
        return None


def _position(value):
    match = re.fullmatch(r"(\d+(?:\.\d+)?)(?:-(\d+(?:\.\d+)?))?", str(value))
    if not match:
        return None
    return (float(match[1]) + float(match[2] or match[1])) / 2


def parse_trim_galore(path: str | Path) -> dict:
    text = Path(path).read_text(encoding="utf-8-sig")
    values = {}
    patterns = {
        "reads_processed": r"Total reads processed:\s*([\d,]+)",
        "reads_with_adapters": r"Reads with adapters:\s*([\d,]+)",
        "cutadapt_reads_written": r"Reads written \(passing filters\):\s*([\d,]+)",
        "bases_processed": r"Total basepairs processed:\s*([\d,]+)",
        "quality_trimmed_bases": r"Quality-trimmed:\s*([\d,]+)",
        "bases_written": r"Total written \(filtered\):\s*([\d,]+)",
        "quality_cutoff": r"Quality Phred score cutoff:\s*(\d+)",
        "min_paired_length": r"Minimum required sequence length for both reads.*?:\s*(\d+)",
    }
    for key, pattern in patterns.items():
        found = re.search(pattern, text)
        values[key] = int(found[1].replace(",", "")) if found else None
    for key, pattern in {"trim_galore_version": r"Trim Galore version:\s*([^\n]+)",
                         "cutadapt_version": r"Cutadapt version:\s*([^\n]+)"}.items():
        match = re.search(pattern, text)
        values[key] = match[1].strip() if match else None
    values["count_scope"] = "Per-mate cutadapt summary; reads-written is not final paired-read retention"
    return values


def parse_star_log(path: str | Path) -> dict:
    """Retain all native STAR fields; numeric percentages are not recomputed."""
    result = {}
    for line in Path(path).read_text(encoding="utf-8-sig").splitlines():
        if "|" in line:
            key, value = line.split("|", 1)
            result[key.strip()] = value.strip()
    return result


def _identity(path: Path, sample_id: str | None = None) -> tuple[str, str]:
    match = re.search(r"(?:_|\.)R?([12])(?:_fastqc|\.fastq|\.fq|\.gz)", path.name)
    mate = "R" + match[1] if match else "not_recorded"
    if sample_id:
        return str(sample_id), mate
    sample = re.sub(r"(?:_|\.)R?[12](?:_fastqc|\.fastq|\.fq|\.gz).*$", "", path.name)
    sample = re.sub(r"(?:_fastqc\.(?:zip|html)|Log\.final\.out|_trimming_report\.txt)$", "", sample).rstrip("_.")
    return sample or path.parent.name, mate


def build_qc_assets(qc_root: str | Path, output_dir: str | Path) -> dict:
    """Create report-only figures/tables from existing completed QC files.

    An optional QC_FETCH_MANIFEST.json binds relocated files to exact sample IDs,
    historical job IDs, source byte sizes and SHA256s. Without it, ordinary FastQC
    names are accepted and sample/mate naming is explicitly recorded. No raw read,
    BAM, alignment or QC executable is read or launched by this function.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap
    from matplotlib.patches import Patch
    import numpy as np
    import pandas as pd

    root, out = Path(qc_root).resolve(), Path(output_dir).resolve()
    if out == root or out.is_relative_to(root) or root.is_relative_to(out):
        raise ValueError("QC report rendering must use an independent output directory")
    if out.exists():
        raise FileExistsError("Existing QC report output is preserved")
    manifest = root / "QC_FETCH_MANIFEST.json"
    inventory = json.loads(manifest.read_text(encoding="utf-8"))["files"] if manifest.is_file() else []
    by_path = {(root / row["relative_path"]).resolve(): row for row in inventory}
    paths = [path for path in root.rglob("*") if path.is_file() and
             (path.name.endswith(("_fastqc.zip", "_fastqc.html", "_trimming_report.txt", "Log.final.out")))]
    source_audit = []
    for path in paths:
        stat = path.stat()
        digest = sha256(path.read_bytes()).hexdigest()
        row = by_path.get(path, {})
        if row and (stat.st_size != row["recorded_size"] or digest != row["sha256"]):
            raise ValueError("Completed QC source verification failed: " + path.name)
        source_audit.append({"relative_path": path.relative_to(root).as_posix(), "bytes": stat.st_size,
                             "mtime_ns": stat.st_mtime_ns, "sha256": digest,
                             "sample_id": row.get("sample_id"), "source_native_id": row.get("native_id")})
    out.mkdir(parents=True)
    figures = out / "figures"
    tables = out / "tables"
    figures.mkdir()
    tables.mkdir()
    artifacts, basic_rows, status_rows, module_rows, trim_rows, star_rows, curves = [], [], [], [], [], [], {}
    method = "Report-only rendering of original completed FastQC, Trim Galore and STAR outputs; no QC, trimming or alignment rerun."

    def record(path, title, kind="file", **meta):
        identifier = "qc-" + sha256((title + str(path)).encode()).hexdigest()[:18]
        data = {"artifact_id": identifier, "analysis_id": "quality.reads", "name": title,
                "path": str(path), "kind": kind, "status": "complete", "run_fingerprint": VERSION,
                "contrast_ids": [], "parameters": {},
                "metadata": {"workflow_version": VERSION, "method": method, "no_upstream_rerun": True,
                             "native_QC_values_preserved": True, **meta}}
        artifacts.append(data)
        return data

    def table(frame, stem, title):
        path = tables / (stem + ".tsv")
        frame.to_csv(path, sep="\t", index=False)
        record(path, title, "table")
        return path

    for path in paths:
        info = by_path.get(path, {})
        sample, mate = _identity(path, info.get("sample_id"))
        if path.name.endswith("_fastqc.html"):
            record(path, "Full original FastQC report | " + sample + " | " + mate,
                   report_asset_kind="html", display_in_main=False, sample_id=sample, mate=mate)
        elif path.name.endswith("_fastqc.zip"):
            parsed = parse_fastqc(path)
            basic_rows.append({"sample_id": sample, "mate": mate, "FastQC_version": parsed["version"], **parsed["basic"]})
            for name, module in parsed["modules"].items():
                status_rows.append({"sample_id": sample, "mate": mate, "module": name, "native_status": module["status"]})
                for row in module["rows"]:
                    module_rows.append({"sample_id": sample, "mate": mate, "module": name,
                                        "columns_json": json.dumps(module["columns"]), "values_json": json.dumps(row)})
                    if len(row) >= 2:
                        curves.setdefault(name, []).append((sample, mate, module["columns"], row))
        elif path.name.endswith("_trimming_report.txt"):
            trim_rows.append({"sample_id": sample, "mate": mate, **parse_trim_galore(path)})
        elif path.name.endswith("Log.final.out"):
            sample = str(info.get("sample_id") or path.parent.name)
            star_rows.append({"sample_id": sample, **parse_star_log(path)})

    basic, status, trim, star = map(pd.DataFrame, (basic_rows, status_rows, trim_rows, star_rows))
    native_tables = {}
    for frame, stem, title in [(basic, "fastqc_basic_statistics", "Complete FastQC basic statistics"),
                               (status, "fastqc_module_status", "Complete native FastQC module statuses"),
                               (pd.DataFrame(module_rows), "fastqc_all_native_module_values", "Complete original FastQC module table cells"),
                               (trim, "trim_galore_metrics", "Complete Trim Galore / cutadapt metrics"),
                               (star, "star_mapping_metrics", "Complete native STAR alignment metrics")]:
        if len(frame):
            native_tables[stem] = table(frame, stem, title)

    def save(fig, stem, title, caption, *, plot_data, unit, full_results):
        fig.suptitle(title, fontsize=14, color="#1E2630")
        for ext in ("png", "svg"):
            fig.savefig(figures / f"{stem}.{ext}", dpi=170, bbox_inches="tight", facecolor="white")
        plt.close(fig)
        plotted = table(plot_data, stem + '_plot_data', 'Exact plot data | ' + title)
        methods = tables / (stem + '_methods.json')
        methods.write_text(json.dumps({'schema': 'katana.report_plot_methods/v1', 'figure_id': stem,
                                      'method': method, 'display_transform': caption, 'unit': unit,
                                      'contrast_id': 'all', 'policy': None, 'native_QC_values_preserved': True,
                                      'upstream_jobs': 0}, indent=2), encoding='utf-8')
        def link(role, path, label, declared_unit):
            return dict(role=role, path=str(path), label=label, contrast_id='all', policy=None, unit=declared_unit)
        links = [link('plot_data', plotted, 'Plot data', unit),
                 link('full_results', full_results, 'Full native QC results', 'native_QC_measurement'),
                 link('methods', methods, 'Methods', 'method_record')]
        record(figures / (stem + ".png"), title, caption=caption,
               data_links=links, exports={"png": str(figures / (stem + ".png")), "svg": str(figures / (stem + ".svg"))})

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9,
                         "axes.spines.top": False, "axes.spines.right": False})
    if len(basic):
        fig, ax = plt.subplots(figsize=(13, 5), layout="constrained")
        counts = basic.pivot(index="sample_id", columns="mate", values="Total Sequences").apply(pd.to_numeric) / 1e6
        counts.plot.bar(ax=ax, color=[BLUE, RED][:len(counts.columns)], width=.8)
        ax.set(xlabel="Sample", ylabel="Original FASTQ reads (millions)")
        ax.tick_params(axis="x", labelrotation=65)
        plotted = counts.rename_axis('sample_id').reset_index().melt(id_vars='sample_id', var_name='mate', value_name='reads_millions')
        save(fig, "fastqc_read_counts", "Sequencing read counts", "Native FastQC total sequences per FASTQ file. R1/R2 are paired read files, not additional biological replicates.",
             plot_data=plotted, unit='million_reads_per_FASTQ', full_results=native_tables['fastqc_basic_statistics'])
    if len(status):
        order = status.assign(library=status.sample_id + " | " + status.mate).pivot(index="library", columns="module", values="native_status")
        coded = order.apply(lambda column: column.map({"pass": 0, "warn": 1, "fail": 2})).to_numpy(dtype=float)
        fig, ax = plt.subplots(figsize=(14, 12), layout="constrained")
        ax.imshow(coded, cmap=ListedColormap([BLUE, "#D49A25", RED]), vmin=0, vmax=2, aspect="auto", interpolation="nearest")
        ax.set_xticks(range(len(order.columns)), order.columns, rotation=65, ha="right")
        ax.set_yticks(range(len(order)), order.index, fontsize=7)
        ax.legend(handles=[Patch(color=c, label=s) for c, s in [(BLUE, "Pass"), ("#D49A25", "Warn"), (RED, "Fail")]], loc="upper left", bbox_to_anchor=(1, 1))
        plotted = status.assign(library=status.sample_id + ' | ' + status.mate).copy()
        plotted['plot_row'] = plotted.library.map({value: i for i, value in enumerate(order.index)})
        plotted['plot_column'] = plotted.module.map({value: i for i, value in enumerate(order.columns)})
        plotted['plot_status_code'] = plotted.native_status.map({'pass': 0, 'warn': 1, 'fail': 2})
        save(fig, "fastqc_checks", "Native FastQC module checks", "Original FastQC pass / warn / fail flags. All modules and libraries are shown; these technical flags do not automatically exclude a sample.",
             plot_data=plotted, unit='FastQC_native_module_status', full_results=native_tables['fastqc_module_status'])

    specifications = [
        ("Per base sequence quality", "Mean", "Base position (native bin midpoint)", "Mean Phred quality", "fastqc_base_quality", "Sequencing quality along reads", False),
        ("Per sequence quality scores", "Count", "Mean read Phred quality", "Reads (% within FASTQ)", "fastqc_read_quality", "Per-read quality distribution", True),
        ("Per sequence GC content", "Count", "GC content (%)", "Reads (% within FASTQ)", "fastqc_gc_distribution", "Read GC content distribution", True),
        ("Sequence Length Distribution", "Count", "Read length (native bin midpoint, bp)", "Reads (% within FASTQ)", "fastqc_read_lengths", "Read length distribution", True),
        ("Per base N content", "N-Count", "Base position (native bin midpoint)", "N content (%)", "fastqc_N_content", "Ambiguous-base content along reads", False),
    ]
    for name, value_column, xlabel, ylabel, stem, title, normalize in specifications:
        items = curves.get(name, [])
        if not items:
            continue
        grouped = {}
        for sample, mate, columns, row in items:
            if value_column not in columns:
                continue
            index = columns.index(value_column)
            x, y = _position(row[0]), _number(row[index]) if len(row) > index else None
            if x is not None and y is not None:
                grouped.setdefault((sample, mate), []).append((x, y))
        if not grouped:
            continue
        fig, ax = plt.subplots(figsize=(12, 5), layout="constrained")
        plotted = []
        for (sample, mate), values in grouped.items():
            values = sorted(values)
            xs, ys = np.asarray(values, dtype=float).T
            native_y = ys.copy()
            if normalize and ys.sum() > 0:
                ys = ys / ys.sum() * 100
            plotted.extend(dict(sample_id=sample, mate=mate, module=name, native_value_column=value_column,
                                x_bin_midpoint=float(x), native_y=float(ny), plotted_y=float(y),
                                display_normalized_within_FASTQ=normalize)
                           for x, ny, y in zip(xs, native_y, ys))
            ax.plot(xs, ys, color=BLUE if mate == "R1" else RED, alpha=.4,
                    marker="o" if len(xs) < 3 else None, markersize=3, lw=1)
        ax.set(xlabel=xlabel, ylabel=ylabel)
        ax.legend(handles=[Patch(color=BLUE, label="Read 1"), Patch(color=RED, label="Read 2")], frameon=False)
        save(fig, stem, title, "One curve per FASTQ file. Native FastQC bins are plotted at their midpoint; complete raw values remain in outputs. " + ("Counts are divided only by the file's total in this native distribution for display." if normalize else "Original QC values are preserved."),
             plot_data=pd.DataFrame(plotted), unit=ylabel, full_results=native_tables['fastqc_all_native_module_values'])
    if len(trim) and "reads_processed" in trim:
        good = trim.dropna(subset=["reads_processed", "reads_with_adapters", "cutadapt_reads_written"])
        if len(good):
            labels = good.sample_id + " | " + good.mate
            fig, axes = plt.subplots(2, 1, figsize=(14, 8), sharex=True, layout="constrained")
            axes[0].bar(range(len(good)), good.reads_with_adapters / good.reads_processed * 100, color=BLUE)
            axes[1].bar(range(len(good)), good.cutadapt_reads_written / good.reads_processed * 100, color=RED)
            axes[0].set(ylabel="Reads with adapter (%)")
            axes[1].set(ylabel="Cutadapt reads written (%)", xticks=range(len(good)), xticklabels=labels)
            axes[1].tick_params(axis="x", labelrotation=90, labelsize=7)
            plotted = good.copy()
            plotted['plot_index'] = range(len(good))
            plotted['adapter_reads_percent'] = good.reads_with_adapters / good.reads_processed * 100
            plotted['cutadapt_reads_written_percent'] = good.cutadapt_reads_written / good.reads_processed * 100
            save(fig, "trimming", "Adapter trimming and cutadapt output", "Per-mate native cutadapt summaries. Adapter-containing reads can be retained after trimming; reads-written here is not final Trim Galore paired-read retention.",
                 plot_data=plotted, unit='per_mate_read_percentage', full_results=native_tables['trim_galore_metrics'])
    mapping = [("Uniquely mapped reads %", "Unique", BLUE), ("% of reads mapped to multiple loci", "Multiple loci", RED),
               ("% of reads mapped to too many loci", "Too many loci", "#9A64B0"),
               ("% of reads unmapped: too many mismatches", "Too many mismatches", "#D49A25"),
               ("% of reads unmapped: too short", "Too short", "#CDD1D6"),
               ("% of reads unmapped: other", "Other unmapped", "#218D87")]
    if len(star) and all(key in star for key, _, _ in mapping):
        values = star.sort_values("sample_id")
        fig, ax = plt.subplots(figsize=(13, 6), layout="constrained")
        bottom = np.zeros(len(values))
        plotted = []
        for field, label, color in mapping:
            y = values[field].map(_number).to_numpy(dtype=float)
            ax.bar(range(len(values)), y, bottom=bottom, color=color, label=label)
            plotted.extend(dict(sample_id=sample, category=label, native_field=field,
                                native_text=text, plotted_percentage=float(height), stack_bottom_percentage=float(base))
                           for sample, text, height, base in zip(values.sample_id, values[field], y, bottom))
            bottom += y
        ax.set(ylabel="Native STAR input-read percentage", xticks=range(len(values)), xticklabels=values.sample_id)
        ax.tick_params(axis="x", labelrotation=65)
        ax.legend(frameon=False, ncol=3)
        save(fig, "mapping", "Genome alignment outcomes", "Original rounded STAR percentages; categories can sum slightly away from 100% because of native rounding. STAR input-read units follow the original paired-end alignment log.",
             plot_data=pd.DataFrame(plotted), unit='native_STAR_input_read_percentage', full_results=native_tables['star_mapping_metrics'])
    unavailable = [{"chapter": 1, "title": "Transcript 5′ / 3′ coverage bias", "status": "not_computed",
                    "reason": "No completed RSeQC gene-body coverage / transcript 5′–3′ bias result was supplied. Read-position FastQC quality is a different measurement; no substitute plot or upstream rerun is created."}]
    for row in source_audit:
        path = root / row["relative_path"]
        if path.stat().st_size != row["bytes"] or path.stat().st_mtime_ns != row["mtime_ns"] or sha256(path.read_bytes()).hexdigest() != row["sha256"]:
            raise RuntimeError("Original completed QC source changed during report rendering")
    receipt = {"schema": "katana.completed_qc_report/v1", "version": VERSION,
               "observed_utc": datetime.now(timezone.utc).isoformat(), "passed": True,
               "source_files": source_audit, "fastqc_files": len(basic), "trim_files": len(trim),
               "alignment_logs": len(star), "figures": len([r for r in artifacts if r["path"].endswith(".png")]),
               "original_html_reports": len([r for r in artifacts if r["path"].endswith(".html")]),
               "sources_preserved": True, "upstream_jobs": 0, "BAM_files": 0,
               "gene_body_coverage": "not_computed"}
    (out / "qc-report.receipt.json").write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    (out / "qc-artifacts.json").write_text(json.dumps(artifacts, indent=2) + "\n", encoding="utf-8")
    return {"artifacts": artifacts, "unavailable": unavailable, "receipt": receipt}
