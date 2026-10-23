"""Draw completed, explicitly scoped ORA tables without calculating enrichment.

Each library keeps its original correction family. GO category filters and the
minimum overlap below affect the displayed panel only, never the saved tests.
"""
from __future__ import annotations

import hashlib
from pathlib import Path
import shutil

import numpy as np
import pandas as pd

from .comparisons import BLUE, RED, _axes, _safe
from .sources import read_table, table_schema


CONTRACT = "katana.report_scoped_ora/v1"


def _scope_inputs(metadata, scope, frame):
    """Resolve the producer's exact set declarations, never infer from names."""
    declared = metadata.get("report_enrichment_scope_audit")
    if not declared:
        return []
    root = Path(metadata["report_enrichment_root"]).resolve()
    audit = pd.read_csv(declared, sep="\t", dtype=str)
    libraries = frame.library.dropna().astype(str).unique() if "library" in frame else []
    if len(libraries) != 1:
        raise ValueError("A scoped ORA panel requires one explicit library identity")
    chosen = audit
    for key in ("comparison_id", "evidence", "report_policy"):
        chosen = chosen.loc[chosen[key].eq(str(scope[key]))]
    chosen = chosen.loc[chosen.library.eq(libraries[0])]
    if len(chosen) != 1 or chosen.iloc[0].fdr_family != frame.fdr_family.iloc[0]:
        raise ValueError("ORA foreground/background declaration does not match its comparison and FDR family")
    row = chosen.iloc[0]
    links = []
    for key, role, label in (("foreground_path", "foreground", "Input genes"),
                             ("background_path", "background", "Detected / testable genes"),
                             ("actual_foreground_path", "foreground", "Input genes mapped to this library"),
                             ("actual_background_path", "background", "Actual enrichment background")):
        path = (root / row[key]).resolve()
        if not path.is_relative_to(root) or not path.is_file():
            raise ValueError("Missing or unbounded declared enrichment set: " + str(path))
        links.append({"role": role, "label": label, "path": str(path),
                      "contrast_id": str(scope["comparison_id"]),
                      "policy": None if scope["report_policy"] == "common" else scope["report_policy"],
                      "unit": "gene", "fdr_family": row.fdr_family,
                      "library": libraries[0]})
    for key, label in (("report_enrichment_mapping", "Gene mapping and library identifiers"),
                       ("report_enrichment_methods", "Methods and resource versions")):
        if metadata.get(key):
            path = Path(metadata[key]).resolve()
            if not path.is_file():
                raise FileNotFoundError(path)
            links.append({"role": "methods", "label": label, "path": str(path),
                          "contrast_id": "all", "policy": None, "scope": "shared", "unit": "resource"})
    return links


def enrichment_attachment_tables(metadata):
    """Expose the producer's explicit complete attachment manifest.

    No directory is guessed or scanned: the completed producer names every
    foreground, background, actual library intersection, mapping and provenance
    file. JSON/HTML receipts remain file links, never parsed as result tables.
    """
    rows = []
    for raw in metadata.get("report_enrichment_attachments", []):
        record = dict(raw)
        path = Path(record.get("path", "")).expanduser().resolve()
        if not path.is_file():
            raise FileNotFoundError("Declared enrichment attachment is unavailable: " + str(path))
        identity = hashlib.sha256(str(path).encode("utf-8")).hexdigest()[:16]
        record.update({"id": "ora-attachment-" + identity,
                       "path": str(path), "source": str(path),
                       "name": record.get("name", path.name), "chapter": record.get("chapter", 5)})
        record["kind"] = record.get("kind", "file" if path.suffix.lower() in {".json", ".html"} else "table")
        if record["kind"] != "file":
            record.update(table_schema(path))
        rows.append(record)
    return rows


def build_enrichment_views(project, artifacts, destination):
    """Read only certified report-scoped completed tables; do not run ORA/BH."""
    output = Path(destination).resolve()
    panels, tables, audit = [], [], []
    for raw in artifacts:
        artifact = raw.model_dump(mode="json") if hasattr(raw, "model_dump") else dict(raw)
        meta = artifact.get("metadata", {})
        if meta.get("report_enrichment_contract") != CONTRACT:
            continue
        tables.extend(enrichment_attachment_tables(meta))
        source = Path(artifact.get("path", ""))
        if not source.is_file() or artifact.get("status", "complete") not in ("complete", "completed"):
            audit.append({"name": artifact.get("name"), "status": "not_computed", "reason": "Completed enrichment table unavailable."})
            continue
        required = {"comparison_id", "evidence", "report_policy", "library_kind", "ontology_category",
                    "term", "pvalue", "fdr", "fdr_family", "overlap_genes", "selected_genes", "background_genes"}
        schema = table_schema(source)
        if not required.issubset(schema["columns"]):
            raise ValueError("Scoped ORA table lacks explicit comparison, policy, resource, or correction family.")
        stamp = (source.stat().st_size, source.stat().st_mtime_ns)
        with source.open("rb") as stream:
            digest = hashlib.file_digest(stream, "sha256").hexdigest()
        directory = output / "outputs/enrichment"
        directory.mkdir(parents=True, exist_ok=True)
        target = directory / (digest[:12] + "-" + source.name)
        if source.resolve() != target.resolve():
            shutil.copy2(source, target)
        tables.append({"id": "scoped-ora-" + digest[:12], "name": "Complete measured-background ORA tests",
                       "path": target.relative_to(output).as_posix(), "source": target.relative_to(output).as_posix(),
                       "rows": schema["rows"], "columns": schema["columns"], "unit": "library_term_test", "chapter": 5})
        keys = ["comparison_id", "evidence", "report_policy", "library_kind", "ontology_category"]
        scopes = read_table(source, columns=keys).drop_duplicates()
        for scope in scopes.to_dict("records"):
            if str(scope["comparison_id"]) not in {c["contrast_id"] for c in project.get("contrasts", [])} and project.get("contrasts"):
                continue
            filters = [(key, "=", value) for key, value in scope.items()]
            frame = read_table(source, filters=filters) if source.suffix == ".parquet" else read_table(source)
            if source.suffix != ".parquet":
                for key, value in scope.items():
                    frame = frame.loc[frame[key].eq(value)]
            # A heterogeneous correction family cannot silently become one
            # claimed family in a single scope/panel.
            families = list(frame.fdr_family.dropna().astype(str).unique())
            if len(families) != 1:
                raise ValueError("A displayed ORA scope must declare exactly one correction family.")
            visible = frame.loc[pd.to_numeric(frame.overlap_genes, errors="coerce").ge(2)].copy()
            visible = visible.sort_values(["fdr", "pvalue", "term"], kind="stable").head(12)
            category = str(scope["ontology_category"])
            evidence = str(scope["evidence"])
            comparison = str(scope["comparison_id"])
            policy = scope["report_policy"] if scope["report_policy"] != "common" else None
            chapter = 7 if scope["library_kind"] == "RBP_identity" else 5
            title = f"{category} | {evidence} | {comparison}"
            fig, ax = _axes(title, "Contributing genes / selected genes mapped to this library", "")
            fig.set_size_inches(10.5, max(4.2, min(8.5, .38 * len(visible) + 2)))
            if visible.empty:
                ax.text(.5, .5, "No terms meet the display-only overlap ≥2 rule\nAll eligible tests remain in the complete table.",
                        transform=ax.transAxes, ha="center", va="center")
            else:
                from matplotlib.colors import LinearSegmentedColormap
                visible = visible.iloc[::-1]
                ratio = pd.to_numeric(visible.overlap_genes).div(pd.to_numeric(visible.selected_genes).replace(0, np.nan))
                q = pd.to_numeric(visible.fdr, errors="coerce")
                color = -np.log10(q.clip(lower=np.nextafter(0., 1.)))
                hits = pd.to_numeric(visible.overlap_genes)
                scatter = ax.scatter(ratio, np.arange(len(visible)), s=35 + 350 * hits / max(1, hits.max()),
                                     c=color, cmap=LinearSegmentedColormap.from_list("ora_blue_red", [BLUE, RED]),
                                     edgecolors="white", linewidths=.5)
                terms = visible.term.astype(str).str.split(" / ", n=1).str[-1].str.slice(0, 105)
                ax.set_yticks(np.arange(len(visible)), terms, fontsize=7)
                fig.colorbar(scatter, ax=ax, shrink=.7, label="−log10(ORA FDR)")
            identity = _safe("ORA-" + "-".join(str(scope[key]) for key in keys))
            data_directory = output / "outputs/enrichment/plot-data"
            data_directory.mkdir(parents=True, exist_ok=True)
            plotted = visible.copy()
            plotted["plot_order"] = np.arange(len(plotted))
            plotted["plot_gene_ratio"] = pd.to_numeric(plotted.overlap_genes).div(pd.to_numeric(plotted.selected_genes).replace(0, np.nan))
            plotted["plot_neglog10_fdr"] = -np.log10(pd.to_numeric(plotted.fdr).clip(lower=np.nextafter(0., 1.)))
            plot_path = data_directory / (identity + "-plotted.tsv")
            plotted.to_csv(plot_path, sep="\t", index=False)
            full_path = data_directory / (identity + "-complete-tests.parquet")
            frame.to_parquet(full_path, index=False)
            links = [{"role": "plot_data", "label": "Plotted terms and contributing genes", "path": str(plot_path),
                      "contrast_id": comparison, "policy": policy, "unit": "library_term_test", "fdr_family": families[0]},
                     {"role": "full_results", "label": "All tested terms in this panel", "path": str(full_path),
                      "contrast_id": comparison, "policy": policy, "unit": "library_term_test", "fdr_family": families[0]},
                     *_scope_inputs(meta, scope, frame)]
            image = output / "figures/enrichment" / (identity + ".png")
            image.parent.mkdir(parents=True, exist_ok=True)
            fig.savefig(image, dpi=145, bbox_inches="tight", facecolor="white")
            fig.savefig(image.with_suffix(".svg"), bbox_inches="tight", facecolor="white")
            fig.clear()
            rbp_note = (" RBP identity membership only; this is not an RBP-target, motif-binding, or causal-regulation test."
                        if chapter == 7 else " GO categories are display subsets of the original GO-library correction family.")
            panels.append({"id": identity, "chapter": chapter, "title": title,
                           "path": image.relative_to(output).as_posix(), "contrast_id": comparison,
                           "scope": "comparison", "policy": policy,
                           "data_links": links,
                           "caption": "Completed downstream hypergeometric ORA; foreground and native detected/testable background files are in outputs. "
                                      + "BH family: " + families[0] + ". Top 12 displayed by FDR then P; overlap ≥2 is display only. "
                                      + "All eligible zero-overlap and nonsignificant tests remain in the complete table." + rbp_note})
            audit.append({**scope, "status": "rendered_completed_tests", "eligible_tests": len(frame),
                          "displayed_terms": len(visible), "fdr_family": families[0], "new_inference": False})
        if stamp != (source.stat().st_size, source.stat().st_mtime_ns):
            raise RuntimeError("Completed enrichment input changed during rendering.")
    return {"panels": panels, "tables": tables, "audit": audit, "new_scientific_fits": 0}
