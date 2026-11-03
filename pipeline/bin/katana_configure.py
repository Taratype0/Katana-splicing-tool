#!/usr/bin/env python3
"""Generate internal Katana configuration from public Nextflow parameters.

Standard-library only: runs in the control task, before scientific calculations.
Resource bytes remain at their original paths and are staged by Nextflow.
"""
from __future__ import annotations

import argparse
import base64
import hashlib
import json
import math
import os
import urllib.request
import urllib.error
import tempfile
import gzip
from collections import defaultdict
from pathlib import Path

ANALYSES = ("landscape", "cross_comparison", "enrichment", "rbp_identity",
            "motiftrace", "motif_binding", "mapp", "network", "interaction", "sequencing_qc")
RESOURCE_KINDS = {"gene_sets": "gene_sets_gmt", "rbp_gene_sets": "rbp_gene_sets",
                  "motif_database": "rbp_motif_database", "clip_peaks": "rbp_binding_peaks",
                  "fimo_runtime": "fimo_runtime"}


def fetch_resource(item, cache):
    """Versioned public downloads: atomic files, one SHA check, cached reuse."""
    from katana_next.flow.file_io import digest, remember, identity
    checksum = item['metadata']['sha256']
    if len(checksum) != 64 or any(c not in '0123456789abcdef' for c in checksum):
        raise ValueError('Invalid resource checksum')
    name = item['filename']
    if Path(name).name != name or name in ('.', '..'):
        raise ValueError('Invalid resource filename')
    folder = Path(cache) / checksum
    folder.mkdir(parents=True, exist_ok=True)
    target = folder / name
    if target.is_file():
        if digest(target) != checksum:
            raise ValueError('Cached resource changed: ' + str(target))
        return target
    url = item['url']
    if not url.startswith(('https://', 'http://')):
        raise ValueError('Resource downloads require HTTP(S)')
    temporary = None
    try:
        req = urllib.request.Request(url, headers={'User-Agent': 'Katana-resource-download/1'})
        with urllib.request.urlopen(req, timeout=120) as response, tempfile.NamedTemporaryFile(dir=folder, delete=False) as stream:
            temporary = Path(stream.name)
            h = hashlib.sha256()
            if item.get('transform') == 'mouse_kegg_symbols':
                raw = response.read()
                if hashlib.sha256(raw).hexdigest() != item['source_sha256']:
                    raise ValueError('KEGG upstream source differs from pinned version')
                mapping = fetch_resource(item['mapping_resource'], cache)
                folded = defaultdict(set)
                for symbol in mapping.read_text().splitlines():
                    folded[symbol.casefold()].add(symbol)
                sets = {}
                for line in raw.decode().splitlines():
                    fields = line.strip().split('\t')
                    if len(fields) < 3:
                        continue
                    genes = set()
                    for gene in fields[2:]:
                        matches = folded[gene.strip().casefold()]
                        if len(matches) == 1:
                            genes.update(matches)
                    sets[fields[0]] = genes
                content = ''.join(term + '\t' + item.get('content_source_url', item['url']) + '\t' + '\t'.join(sorted(genes)) + '\n'
                                  for term, genes in sorted(sets.items()) if genes).replace('\n', '\r\n').encode()
                stream.write(content)
                h.update(content)
            else:
                reader = gzip.GzipFile(fileobj=response) if item.get('compression') == 'gzip' else response
                for block in iter(lambda: reader.read(1024 * 1024), b''):
                    stream.write(block)
                    h.update(block)
        if h.hexdigest() != checksum:
            raise ValueError('Downloaded resource checksum differs from the pinned release: ' + name)
        # A concurrent process may have finished the same immutable download.
        if target.exists():
            if digest(target) != checksum:
                raise ValueError('Concurrent resource download differs: ' + name)
        else:
            os.replace(temporary, target)
            remember(target, checksum, identity(target))
        return target
    except (urllib.error.URLError, TimeoutError) as exc:
        raise RuntimeError(
            f"Unable to download resource {name} from {url}. "
            "KATANA_FLOW_CONFIGURE should run locally on the Nextflow launch host; "
            "that host needs network access or a pre-populated --resource_cache. "
            f"Existing cached resources were retained. Cause: {exc}"
        ) from exc
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def resource_entries(request, resources):
    names = request.get('resource_names', [])
    if len(names) != len(resources) or set(names) - set(RESOURCE_KINDS):
        raise ValueError('Resource names and staged inputs do not match')
    declarations = request.get('resource_entries')
    if declarations is not None and len(declarations) != len(names):
        raise ValueError('Resource declarations and staged inputs do not match')
    entries = []
    for index, (name, filename) in enumerate(zip(names, resources)):
        declared = declarations[index] if declarations is not None else {}
        entries.append(dict(parameter=name, path=str(filename), asset_id=declared.get('asset_id', name),
                            metadata=declared.get('metadata', {})))
    catalog = request.get('resource_catalog') or []
    supplied = set(names)
    for item in catalog:
        if item['parameter'] in supplied:
            continue
        entries.append(dict(parameter=item['parameter'], asset_id=item['asset_id'],
            path=str(fetch_resource(item, request['resource_cache'])), metadata=item['metadata']))
    return entries


def boolean(value, name):
    if not isinstance(value, bool):
        raise ValueError(f"--{name} must be true or false")
    return value


def number(request, name, default, minimum, maximum=None):
    value = request.get(name, default)
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ValueError(f"--{name} must be a finite number")
    if value < minimum or (maximum is not None and value > maximum):
        raise ValueError(f"--{name} must be >= {minimum}" + (f" and <= {maximum}" if maximum is not None else ""))
    return value


def generate(request, resources, outdir):
    factors = request.get('design_factors', ['condition'])
    if isinstance(factors, str):
        factors = [item.strip() for item in factors.split(',')]
    if (not isinstance(factors, list) or not factors or
            any(not isinstance(item, str) or not item.strip() for item in factors) or
            len(set(factors)) != len(factors)):
        raise ValueError('--design_factors requires unique nonempty factor names')
    analyses = {key: not boolean(request.get("skip_" + key, False), "skip_" + key) for key in ANALYSES}
    analyses["coverage"] = boolean(request.get("with_genebody_coverage", False), "with_genebody_coverage")
    event_set = request.get("event_set", "stable_significant")
    if event_set not in ("stable_significant", "significant"):
        raise ValueError("--event_set must be stable_significant (RS-DAS) or significant (ordinary DAS)")
    rmats_mode = request.get("rmats_mode", "JCEC")
    if rmats_mode not in ("JC", "JCEC"):
        raise ValueError("--rmats_mode must be JC or JCEC")
    settings = {"schema": "katana.flow.config/v1", "event_set": event_set, "analyses": analyses,
                "network": {"engine": "whole_study_4_3", "scope": "whole_study", "bootstrap": {"enabled": boolean(
                    request.get("with_network_bootstrap", False), "with_network_bootstrap")},
                    "parameters": {}}, "parameters": {}}
    size = number(request, "network_min_module_size", 20, 2)
    height = number(request, "network_merge_cut_height", .25, 0, 1)
    if int(size) != size or not 0 < height < 1:
        raise ValueError("network_min_module_size must be an integer >=2; network_merge_cut_height must be strictly between 0 and 1")
    settings["network"]["parameters"].update(min_module_size=int(size), merge_cut_height=height)
    mappings = (("deg_fdr", "differential.deg", "fdr", .05, 0, 1),
                ("deg_log2fc", "differential.deg", "min_abs_log2fc", 1., 0, None),
                ("das_fdr", "quality.stable_as", "fdr", .05, 0, 1),
                ("das_dpsi", "quality.stable_as", "min_abs_dpsi", .1, 0, 1),
                ("dtu_fdr", "differential.dtu", "fdr", .05, 0, 1))
    for name, workflow, key, default, low, high in mappings:
        settings["parameters"].setdefault(workflow, {})[key] = number(request, name, default, low, high)
    # Candidate evidence uses the same expression selection as the DEG branch.
    settings["parameters"]["prioritize.evidence"] = {
        "deg_fdr": settings["parameters"]["differential.deg"]["fdr"],
        "dtu_fdr": settings["parameters"]["differential.dtu"]["fdr"],
        "deg_min_abs_log2fc": settings["parameters"]["differential.deg"]["min_abs_log2fc"]}
    selection = request.get("contrast_ids")
    if selection:
        selection = selection.split(",") if isinstance(selection, str) else selection
        if not isinstance(selection, list) or any(not isinstance(x, str) or not x.strip() for x in selection):
            raise ValueError("--contrast_ids must be a comma-separated list or YAML list")
        selection = [x.strip() for x in selection]
        if len(set(selection)) != len(selection):
            raise ValueError("--contrast_ids contains duplicates")
        settings["contrast_ids"] = selection

    settings['design_factors'] = factors
    template = {"name": request.get("study_name") or "Katana study",
                "rmats_mode": rmats_mode, "metadata": {}, "assets": []}
    for option, key in (("species", "species"), ("reference_build", "genome_build"),
                        ("annotation_release", "annotation_release")):
        if request.get(option):
            template["metadata"][key] = str(request[option])
    metadata = request.get("resource_metadata") or {}
    if not isinstance(metadata, dict) or set(metadata) - set(RESOURCE_KINDS):
        raise ValueError("--resource_metadata accepts metadata objects keyed by resource parameter")
    entries = resource_entries(request, resources)
    names = [entry['parameter'] for entry in entries]
    ids = [entry['asset_id'] for entry in entries]
    if len(set(ids)) != len(ids):
        raise ValueError('Resource asset IDs must be unique')
    for entry in entries:
        name, filename = entry['parameter'], entry['path']
        path = Path(filename).resolve(strict=True)
        if not path.is_file():
            raise ValueError(f"--{name} must be a file")
        try:
            from katana_next.flow.file_io import digest as file_digest
            checksum = file_digest(path)
        except ImportError:
            with path.open('rb') as stream:
                checksum = hashlib.file_digest(stream, 'sha256').hexdigest()
        declared = {**(metadata.get(name) or {}), **entry['metadata']}
        if declared.get("sha256", checksum) != checksum:
            raise ValueError(f"--{name}: declared resource SHA256 does not match its bytes")
        for key in ("species", "genome_build"):
            if declared.get(key) and template["metadata"].get(key) and str(declared[key]).replace(" ", "_") != str(template["metadata"][key]).replace(" ", "_"):
                raise ValueError(f"--{name}: {key} conflicts with project declaration")
        record = {"source_url": path.as_uri(), "version": "sha256:" + checksum,
                  "version_kind": "content_identity", **template["metadata"], **declared, "sha256": checksum}
        template["assets"].append({"asset_id": entry['asset_id'], "kind": RESOURCE_KINDS[name],
                                    "path": str(path), "metadata": record})
    raw = request.get("mode", "complete") != "downstream-only"
    if raw:
        required = set()
        if analyses["enrichment"] or analyses["network"]: required.add("gene_sets")
        if analyses["rbp_identity"]: required.add("rbp_gene_sets")
        if analyses["motiftrace"] or analyses["motif_binding"]:
            required.add("motif_database")
            if not request.get("species"):
                raise ValueError("Supply --species for motif analysis, or select --skip_motiftrace --skip_motif_binding")
        if analyses["motif_binding"]: required.add("clip_peaks")
        missing = sorted(required - set(names))
        if missing:
            raise ValueError("Selected analyses need " + ", ".join("--" + x for x in missing) +
                             ". Supply the resources or explicitly skip the corresponding analyses.")
    elif names:
        raise ValueError("Existing-project mode uses the resources in --project; do not supply replacement resource paths")
    outdir = Path(outdir)
    if outdir.exists():
        raise FileExistsError(f"Configuration destination already exists: {outdir}")
    outdir.mkdir(parents=True)
    def write(name, document):
        (outdir / name).write_text(json.dumps(document, indent=2, ensure_ascii=True) + "\n", encoding="utf-8")
    if raw and (analyses["motiftrace"] or analyses["motif_binding"]) and "fimo_runtime" not in names:
        # The installed worker imports pymemesuite in the selected task image.
        # No host Python path or WSL assumption is embedded in this descriptor.
        write("fimo-runtime.json", {"python": "python3"})
        template["assets"].append({"asset_id": "fimo_runtime", "kind": "fimo_runtime", "path": "fimo-runtime.json"})
    write("project-template.json", template)
    write("flow.json", settings)
    write("parameters.json", request)
    return template, settings


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--request-base64", required=True)
    parser.add_argument("--resources", nargs="*", default=[])
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    generate(json.loads(base64.b64decode(args.request_base64)), args.resources, args.outdir)
