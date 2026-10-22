"""Explicit, versioned annotation resources. Opening a project never downloads data."""
from __future__ import annotations

from collections import defaultdict
from datetime import datetime, timezone
from functools import lru_cache
import gzip
import hashlib
import json
from pathlib import Path
import shutil
from urllib.request import Request, urlopen
from zipfile import ZipFile

import pandas as pd

from katana_next.domain.models import InputAsset, ProjectManifest
from katana_next.infrastructure.manifest_store import ManifestStore


def sha256(path: Path) -> str:
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def download(url: str, destination: Path) -> dict:
    destination.parent.mkdir(parents=True, exist_ok=True)
    if not destination.is_file():
        partial = destination.with_suffix(destination.suffix + ".part")
        with urlopen(Request(url, headers={"User-Agent": "Katana/2 resource preparation"}), timeout=120) as source, partial.open("wb") as output:
            shutil.copyfileobj(source, output)
        if partial.stat().st_size == 0:
            raise ValueError(f"Empty annotation download: {url}")
        partial.replace(destination)
    return {"url": url, "path": str(destination.resolve()), "sha256": sha256(destination),
            "bytes": destination.stat().st_size}


def go_gene_sets(obo: Path, gaf: Path, taxon: str) -> tuple[dict[str, set[str]], dict]:
    """Propagate only is_a / part_of, excluding NOT and no-data annotations."""
    names, parents, aliases = {}, defaultdict(set), {}
    ontology_version = ""
    current, obsolete = {}, False
    def finish():
        if current.get("id") and not obsolete:
            key = current["id"]
            names[key] = current.get("name", key)
            parents[key].update(current.get("parents", []))
            for alias in current.get("alt_id", []):
                aliases[alias] = key
    with obo.open(encoding="utf-8") as stream:
        for raw in stream:
            line = raw.strip()
            if line.startswith("data-version:"):
                ontology_version = line.partition(":")[2].strip()
            if line.startswith("["):
                finish()
                current, obsolete = {}, False
            elif line.startswith("id: GO:"):
                current["id"] = line[4:]
            elif line.startswith("name: "):
                current["name"] = line[6:]
            elif line == "is_obsolete: true":
                obsolete = True
            elif line.startswith("alt_id: "):
                current.setdefault("alt_id", []).append(line[8:])
            elif line.startswith("is_a: "):
                current.setdefault("parents", []).append(line.split()[1])
            elif line.startswith("relationship: part_of "):
                current.setdefault("parents", []).append(line.split()[2])
    finish()
    @lru_cache(None)
    def ancestors(term):
        return frozenset({term}).union(*(ancestors(p) for p in parents[term] if p in names))
    members, excluded, used = defaultdict(set), 0, 0
    annotation_headers = []
    opener = gzip.open if gaf.suffix == ".gz" else open
    with opener(gaf, "rt", encoding="utf-8") as stream:
        for line in stream:
            if line.startswith("!"):
                annotation_headers.append(line.strip())
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 15:
                raise ValueError("Malformed GO GAF annotation row")
            source_taxon = fields[12].split("|")[0].split(":")[-1]
            if "NOT" in fields[3].split("|") or fields[6] == "ND" or source_taxon != taxon or fields[0] == "ComplexPortal":
                excluded += 1
                continue
            term = aliases.get(fields[4], fields[4])
            if term not in names or not fields[2]:
                excluded += 1
                continue
            used += 1
            for parent in ancestors(term):
                members[parent].add(fields[2])
    result = {f"{term} {names[term]}": genes for term, genes in sorted(members.items())}
    return result, {"used_annotations": used, "excluded_annotations": excluded,
                    "terms": len(result), "genes": len(set().union(*members.values())) if members else 0,
                    "ontology_version": ontology_version, "annotation_headers": annotation_headers,
                    "propagation": "is_a and part_of only; NOT, ND and protein-complex entities excluded"}


def write_gmt(sets: dict[str, set[str]], path: Path, source: str) -> None:
    if not sets:
        raise ValueError(f"No gene sets available for {path.name}")
    path.write_text("".join(f"{term}\t{source}\t" + "\t".join(sorted(genes)) + "\n"
                            for term, genes in sorted(sets.items()) if genes), encoding="utf-8")


def prepare_mouse_resources(root: Path, attract_zip: Path | None = None) -> Path:
    """Build an immutable local mouse bundle; no human-symbol case conversion."""
    root = root.resolve()
    index = root / "resources.json"
    if index.exists():
        payload = json.loads(index.read_text(encoding="utf-8"))
        if payload.get("builder_version") != 2:
            raise ValueError("This annotation bundle predates taxon/symbol validation. Prepare a new versioned directory.")
        for item in payload["assets"]:
            if sha256(Path(item["path"])) != item["metadata"]["sha256"]:
                raise ValueError("Resource content changed; prepare a new versioned directory")
        return index
    root.mkdir(parents=True, exist_ok=True)
    sources = [download("https://purl.obolibrary.org/obo/go/go-basic.obo", root / "go-basic.obo"),
               download("https://current.geneontology.org/annotations/gaf/MOUSE-mod.gaf.gz", root / "gene_association.mgi.gz"),
               download("https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=KEGG_2019_Mouse", root / "KEGG_2019_Mouse.gmt")]
    sets, audit = go_gene_sets(root / "go-basic.obo", root / "gene_association.mgi.gz", "10090")
    write_gmt(sets, root / "GO_MGI_mouse.gmt", sources[1]["url"])
    from katana_next.science.enrichment import read_gene_sets
    symbols = set().union(*sets.values())
    folded = defaultdict(set)
    for symbol in symbols:
        folded[symbol.casefold()].add(symbol)
    kegg, unmapped, ambiguous = {}, set(), set()
    for term, genes in read_gene_sets(root / "KEGG_2019_Mouse.gmt").items():
        mapped = set()
        for gene in genes:
            matches = folded[gene.casefold()]
            if len(matches) == 1:
                mapped.update(matches)
            elif not matches:
                unmapped.add(gene)
            else:
                ambiguous.add(gene)
        kegg[term] = mapped
    write_gmt(kegg, root / "KEGG_2019_Mouse_symbols.gmt", sources[2]["url"])
    rbp = {term: genes for term, genes in sets.items() if term.startswith(("GO:0003723 ", "GO:0008380 "))}
    assets = []
    def add(kind, path, semantics):
        assets.append({"asset_id": f"mouse::{path.stem}", "kind": kind, "path": str(path),
                       "metadata": {"species": "Mus_musculus", "taxon": "10090", "sha256": sha256(path),
                                    "semantics": semantics, "resource_bundle": str(index)}})
    if attract_zip:
        target = root / "ATtRACT.zip"
        if not target.exists():
            shutil.copy2(attract_zip, target)
        with ZipFile(target) as archive:
            frame = pd.read_csv(archive.open("ATtRACT_db.txt"), sep="\t")
        mouse = frame.loc[frame.Organism.eq("Mus_musculus") & frame.Gene_id.astype(str).str.startswith("ENSMUSG") & frame.Mutated.eq("no")]
        rbp["ATtRACT mouse RBPs with experimental motifs"] = set(mouse.Gene_id.astype(str).str.replace(r"\.\d+$", "", regex=True))
        add("rbp_motif_database", target, "ATtRACT PWMs, filtered by organism during scanning; not target gene sets")
        sources.append({"url": "https://attract.cnic.es/download", "local_source": str(attract_zip.resolve()), "sha256": sha256(target)})
    write_gmt(rbp, root / "RBP_mouse.gmt", "MGI GO annotation and species-matched ATtRACT membership")
    add("gene_sets_gmt", root / "GO_MGI_mouse.gmt", "mouse MGI GO annotations, propagated is_a/part_of")
    add("go_ontology_obo", root / "go-basic.obo", "GO BP/MF/CC category definitions; not gene membership")
    add("gene_sets_gmt", root / "KEGG_2019_Mouse_symbols.gmt", "KEGG 2019 Mouse via Enrichr; unique case-insensitive mapping to mouse MGI symbols; frozen legacy library, not current KEGG")
    assets[-1]["metadata"]["ontology"] = "KEGG"
    add("rbp_gene_sets", root / "RBP_mouse.gmt", "RBP/splicing-factor membership, NOT RBP target sets or causal regulators")
    index.write_text(json.dumps({"builder_version": 2, "created_at": datetime.now(timezone.utc).isoformat(), "species": "Mus_musculus",
                                "kegg_symbol_audit": {"rule": "unique casefold match to mouse MGI symbols only; no orthology conversion", "unmapped": sorted(unmapped), "ambiguous": sorted(ambiguous)},
                                "sources": sources, "go_audit": audit, "assets": assets}, indent=2), encoding="utf-8")
    return index


def add_ontology_categories(bundle: Path, destination: Path) -> Path:
    """Explicitly expose a bundle's already-hashed OBO; no downloads or mutation."""
    if destination.exists():
        raise FileExistsError("Preserve the resource index; choose a new destination")
    payload = json.loads(bundle.read_text(encoding="utf-8"))
    for asset in payload.get("assets", []):
        if sha256(Path(asset["path"])) != asset["metadata"]["sha256"]:
            raise ValueError("Resource content changed")
    if not any(a["kind"] == "go_ontology_obo" for a in payload["assets"]):
        obos = [s for s in payload.get("sources", []) if str(s.get("path", "")).endswith(".obo")]
        if len(obos) != 1:
            raise ValueError("Declare exactly one existing hashed GO OBO in bundle sources")
        source = obos[0]
        path = Path(source["path"])
        if sha256(path) != source["sha256"]:
            raise ValueError("GO ontology fingerprint mismatch")
        payload["assets"].append({"asset_id": "go::ontology_categories", "kind": "go_ontology_obo", "path": str(path),
                                  "metadata": {"sha256": source["sha256"], "source_url": source["url"],
                                               "semantics": "GO category definitions only; original membership and statistics unchanged"}})
    for asset in payload["assets"]:
        if asset["kind"] == "gene_sets_gmt" and str(asset["metadata"].get("semantics", "")).startswith("KEGG "):
            asset["metadata"]["ontology"] = "KEGG"
    payload["category_extension"] = {"source_bundle": str(bundle.resolve()), "source_sha256": sha256(bundle),
                                     "created_at": datetime.now(timezone.utc).isoformat(), "membership_changed": False}
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    return destination


def attach_resources(source: Path, bundle: Path, destination: Path, output_root: Path | None = None) -> Path:
    _, manifest = ManifestStore().load(source)
    payload = json.loads(bundle.read_text(encoding="utf-8"))
    if not payload.get("assets"):
        raise ValueError("Resource bundle contains no assets")
    if any(a.get("asset_id", "").startswith("mouse::") for a in payload["assets"]) and payload.get("builder_version", 0) < 2:
        raise ValueError("Mouse resource bundle requires rebuild with the validated taxon/symbol parser")
    if payload["species"] == "Mus_musculus" and manifest.annotation.genome_build not in {"GRCm39", "GRCm38", "mm39", "mm10"}:
        raise ValueError("Mouse resource bundle does not match the declared genome build")
    assets = [InputAsset.model_validate(item) for item in payload["assets"]]
    for item in assets:
        if item.kind in {"rbp_binding_peaks", "uniprot_protein_json"}:
            from katana_next.science.versioned_evidence import verified_resource
            verified_resource(item, species=manifest.metadata.get("species", ""),
                              genome_build=manifest.annotation.genome_build if item.kind == "rbp_binding_peaks" else None)
        if sha256(item.path) != item.metadata["sha256"]:
            raise ValueError(f"Resource fingerprint mismatch: {item.path}")
    kinds = {item.kind for item in assets}
    # Explicit replacement of these resource kinds; unrelated analysis inputs stay intact.
    update = {"assets": tuple(a for a in manifest.assets if a.kind not in kinds) + tuple(assets)}
    if output_root:
        from uuid import NAMESPACE_URL, uuid5
        update.update(output_root=output_root.resolve(), project_id=str(uuid5(NAMESPACE_URL, str(destination.resolve()))))
    if destination.exists():
        raise FileExistsError(f"Preserve the existing manifest; use a new destination: {destination}")
    return ManifestStore().save(ProjectManifest.model_validate({**manifest.model_dump(), **update}), destination)
