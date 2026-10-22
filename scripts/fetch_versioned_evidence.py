"""Explicit opt-in downloads with frozen provider records and checksums."""
import argparse
import hashlib
import json
from pathlib import Path
import requests


def fetch(output, uniprot=(), encode=()):
    output.mkdir(parents=True, exist_ok=True)
    assets = []
    def get(url, path):
        if path.exists():
            return path.read_bytes()
        response = requests.get(url, headers={"Accept": "application/json"}, timeout=90)
        response.raise_for_status()
        if len(response.content) > 100_000_000:
            raise ValueError("Resource exceeds 100 MB budget")
        path.write_bytes(response.content)
        return response.content
    for accession in uniprot:
        if not accession.isalnum():
            raise ValueError("Use an exact canonical UniProt accession")
        url = f"https://rest.uniprot.org/uniprotkb/{accession}.json"
        path = output/f"{accession}.json"
        raw = get(url, path); entry = json.loads(raw)
        assets.append(dict(asset_id=f"uniprot:{accession}", kind="uniprot_protein_json", path=str(path.resolve()),
            metadata=dict(source_url=url, species=entry["organism"]["scientificName"], sha256=hashlib.sha256(raw).hexdigest(),
                version=f"entry:{entry['entryAudit']['entryVersion']};sequence:{entry['entryAudit']['sequenceVersion']}",
                last_annotation_update=entry["entryAudit"]["lastAnnotationUpdateDate"])))
    for accession in encode:
        if not accession.startswith("ENCFF") or not accession.isalnum():
            raise ValueError("Use an ENCODE file accession")
        url = f"https://www.encodeproject.org/files/{accession}/?format=json"
        record = json.loads(get(url, output/f"{accession}.metadata.json"))
        experiment = json.loads(get("https://www.encodeproject.org"+record["dataset"]+"?format=json", output/f"{accession}.experiment.json"))
        if record.get("status") != "released" or experiment.get("assay_title") not in {"eCLIP", "iCLIP", "PAR-CLIP", "HITS-CLIP"}:
            raise ValueError("Require released CLIP data")
        path = output/f"{accession}.bed.gz"
        raw = get("https://www.encodeproject.org"+record["href"], path)
        if hashlib.md5(raw).hexdigest() != record["md5sum"]:
            raise ValueError("Provider MD5 mismatch")
        assets.append(dict(asset_id=f"encode:{accession}", kind="rbp_binding_peaks", path=str(path.resolve()),
            metadata=dict(source_url=url, species=experiment["target"]["organism"]["scientific_name"], genome_build=record["assembly"],
                sha256=hashlib.sha256(raw).hexdigest(), version=record["date_created"], coordinate_system="0-based-half-open",
                assay=experiment["assay_title"], biosample=experiment["biosample_ontology"]["term_name"], target=experiment["target"]["label"],
                accession=accession, peak_calling=record["output_type"], biological_replicates=record.get("biological_replicates", []))))
    (output/"input_assets.json").write_text(json.dumps(assets, indent=2), encoding="utf-8")
    for species in sorted({a["metadata"]["species"] for a in assets}):
        bundle = {"builder_version": 1, "species": species, "assets": [a for a in assets if a["metadata"]["species"] == species]}
        (output/f"resources_{species.replace(' ', '_')}.json").write_text(json.dumps(bundle, indent=2), encoding="utf-8")
    print(json.dumps(assets, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--uniprot", nargs="*", default=[])
    parser.add_argument("--encode", nargs="*", default=[])
    args = parser.parse_args()
    fetch(args.output, args.uniprot, args.encode)
