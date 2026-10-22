from __future__ import annotations

import json
from urllib.request import Request, urlopen

import pandas as pd


GPROFILER_ENDPOINT = "https://biit.cs.ut.ee/gprofiler/api/gost/profile/"


def gprofiler_enrichment(
    module_genes: dict[str, list[str]],
    *,
    organism: str,
    background: list[str] | None = None,
    timeout_seconds: int = 60,
) -> tuple[pd.DataFrame, dict[str, object]]:
    rows: list[dict[str, object]] = []
    versions: dict[str, object] = {}
    custom_background = sorted({str(gene).strip() for gene in (background or []) if str(gene).strip()})
    for module, genes in sorted(module_genes.items()):
        query = sorted({str(gene).strip() for gene in genes if str(gene).strip()})
        if len(query) < 3:
            continue
        payload = {
            "organism": organism,
            "query": query,
            "sources": ["GO:BP", "GO:MF", "GO:CC", "KEGG"],
            "user_threshold": 1.0,
            "significance_threshold_method": "fdr",
            "no_evidences": True,
        }
        if custom_background:
            payload["domain_scope"] = "custom"
            payload["background"] = custom_background
        request = Request(
            GPROFILER_ENDPOINT,
            data=json.dumps(payload).encode("utf-8"),
            headers={"Content-Type": "application/json", "Accept": "application/json"},
            method="POST",
        )
        with urlopen(request, timeout=timeout_seconds) as response:  # noqa: S310 - fixed HTTPS endpoint
            result = json.loads(response.read().decode("utf-8"))
        versions = result.get("meta", {}).get("data_version", versions)
        for item in result.get("result", []):
            rows.append({
                "module": module,
                "term_id": item.get("native", ""),
                "term_name": item.get("name", ""),
                "source": item.get("source", ""),
                "pvalue": item.get("p_value"),
                "FDR": item.get("p_value"),
                "intersection_size": item.get("intersection_size"),
                "term_size": item.get("term_size"),
                "query_size": item.get("query_size"),
                "effective_domain_size": item.get("effective_domain_size"),
                "intersection_genes": ";".join(map(str, item.get("intersections", []))),
                "method": "g:Profiler g:GOSt with FDR correction",
            })
    columns = [
        "module", "term_id", "term_name", "source", "pvalue", "FDR",
        "intersection_size", "term_size", "query_size", "effective_domain_size",
        "intersection_genes", "method",
    ]
    return pd.DataFrame(rows, columns=columns), {
        "provider": "g:Profiler g:GOSt", "endpoint": GPROFILER_ENDPOINT,
        "organism": organism, "data_version": versions,
        "background": "custom" if background else "annotated_genome",
        "background_size": len(custom_background) if background else None,
    }
