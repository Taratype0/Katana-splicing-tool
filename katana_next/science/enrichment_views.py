"""Display-only ontology filters: full enrichment tables and BH families remain intact."""
from pathlib import Path

import numpy as np
import pandas as pd

from katana_next.workflows.bulk import _plot


CATEGORY_NAMES = {
    "biological_process": "GO BP", "molecular_function": "GO MF", "cellular_component": "GO CC",
    "BP": "GO BP", "MF": "GO MF", "CC": "GO CC", "KEGG": "KEGG",
    "GO BP": "GO BP", "GO MF": "GO MF", "GO CC": "GO CC",
}


def read_go_categories(path: Path) -> dict[str, str]:
    categories, current = {}, {}
    for raw in path.read_text(encoding="utf-8").splitlines() + ["[end]"]:
        line = raw.strip()
        if line.startswith("["):
            if current.get("id") and not current.get("obsolete"):
                category = CATEGORY_NAMES.get(current.get("namespace"), "Unclassified")
                categories[current["id"]] = category
                for alias in current.get("aliases", []):
                    categories[alias] = category
            current = {}
        elif line.startswith("id: GO:"):
            current["id"] = line[4:]
        elif line.startswith("namespace: "):
            current["namespace"] = line[11:]
        elif line == "is_obsolete: true":
            current["obsolete"] = True
        elif line.startswith("alt_id: GO:"):
            current.setdefault("aliases", []).append(line[8:])
    return categories


def annotate_categories(result: pd.DataFrame, manifest) -> pd.DataFrame:
    result = result.copy()
    namespaces = {}
    library_categories = {}
    from katana_next.resources import sha256
    for asset in manifest.assets:
        if asset.kind == "go_ontology_obo":
            expected = asset.metadata.get("sha256")
            if not expected or sha256(asset.path) != expected:
                raise ValueError("GO ontology category resource needs a matching declared SHA256")
            for term, category in read_go_categories(asset.path).items():
                if term in namespaces and namespaces[term] != category:
                    raise ValueError("Conflicting GO ontology namespaces")
                namespaces[term] = category
        if asset.kind == "gene_sets_gmt":
            value = asset.metadata.get("ontology", asset.metadata.get("category", ""))
            if value in CATEGORY_NAMES:
                library_categories[asset.asset_id] = CATEGORY_NAMES[value]
    def category(term):
        library, _, label = str(term).partition(" / ")
        identifier = label.split()[0] if label else ""
        return namespaces.get(identifier, library_categories.get(library, "Unclassified"))
    result["ontology_category"] = result.term.map(category)
    result["gene_ratio"] = pd.to_numeric(result.overlap_genes).div(pd.to_numeric(result.selected_genes).replace(0, np.nan))
    return result


def enrichment_bubble_plot(result, *, category="All", x_axis="gene_ratio", top_terms=12, labels=None):
    if category not in {"All", "GO BP", "GO MF", "GO CC", "KEGG", "Unclassified"}:
        raise ValueError("Unknown enrichment display category")
    if x_axis not in {"gene_ratio", "fold_enrichment"}:
        raise ValueError("Choose gene_ratio or fold_enrichment")
    visible = result.loc[result.get("passes_display_overlap", pd.Series(False, index=result.index)).fillna(False)].copy()
    if category != "All":
        visible = visible.loc[visible.ontology_category.eq(category)]
    if "comparison_id" not in visible:
        visible["comparison_id"] = "selection"
    ranking = ["fdr", *(["pvalue"] if "pvalue" in visible else []), "term"]
    visible = visible.sort_values(ranking).groupby("comparison_id", sort=False).head(int(top_terms))
    visible["term_label"] = visible.term.str.split(" / ", n=1).str[-1]
    visible["comparison_label"] = visible.comparison_id.map(labels or {}).fillna(visible.comparison_id)
    color_max = max(1.3, float(pd.to_numeric(visible.get("minus_log10_fdr", pd.Series(dtype=float)), errors="coerce").max())) if len(visible) else 1.3
    chart = _plot(visible.to_dict("records"), {"type": "circle", "opacity": .9, "stroke": "#52616b", "strokeWidth": .4}, {
        "x": {"field": x_axis, "type": "quantitative", "title": "Gene ratio (hits / input genes)" if x_axis == "gene_ratio" else "Fold enrichment"},
        "y": {"field": "term_label", "type": "nominal", "sort": {"field": "fdr"}, "title": None, "axis": {"labelLimit": 340}},
        "size": {"field": "overlap_genes", "type": "quantitative", "title": "Contributing genes", "scale": {"range": [35, 650]}},
        "color": {"field": "minus_log10_fdr", "type": "quantitative", "title": "-log10(FDR)", "scale": {"scheme": "yelloworangered", "domain": [0, color_max]}},
        "facet": {"field": "comparison_label", "type": "nominal", "columns": 2, "title": "Exact comparison"},
        "tooltip": [{"field": c} for c in ("comparison_label", "term", "ontology_category", "overlap_genes", "selected_genes", "background_genes", "gene_ratio", "fold_enrichment", "fdr", "genes")],
    }, title="Functional enrichment: size = genes, color = FDR", description="Display-only ontology filtering after the declared BH correction. Full tables retain significant and nonsignificant terms.", height=max(300, min(620, 28 * int(top_terms))))
    chart["background"] = "#ffffff"
    chart["resolve"] = {"scale": {"y": "independent"}}
    chart["config"]["font"] = "Georgia"
    for key in ("axis", "legend", "header", "title"):
        chart["config"].setdefault(key, {})
        if key == "title":
            chart["config"][key]["font"] = "Georgia"
        else:
            chart["config"][key].update(labelFont="Georgia", titleFont="Georgia")
    return chart
