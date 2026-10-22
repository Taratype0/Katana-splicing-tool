"""Explicit public-resource coverage ledger; not a claim of global completeness."""
import json
from pathlib import Path

import pandas as pd

RESOURCES = [
    ("CISBP-RNA", "https://cisbp-rna.ccbr.utoronto.ca/bulk.php", "PWM; direct / homology / JPLE model", "scanned", "2.00 mouse archive SHA verified; unavailable native matrices audited; evidence layers separate"),
    ("ATtRACT", "https://attract.cnic.es/download", "curated experimental motifs/PWM", "scanned", "mouse, non-mutated, Ensembl-mapped entries; frozen archive SHA verified"),
    ("RBPDB", "https://rbpdb.ccbr.utoronto.ca/help.php", "curated experiments, PWM and individual binding sequences", "upstream_in_aggregate", "ATtRACT R provenance and CISBP source; not claiming every native RBPDB experiment has a usable mouse PWM"),
    ("RNAcompete / EuPRI", "https://www.nature.com/articles/s41587-025-02733-6", "experimental preferences and PWMs", "upstream_in_CISBP", "included as direct or explicitly inferred mouse assignment; raw array intensities are not rescored"),
    ("HTR-SELEX", "https://pmc.ncbi.nlm.nih.gov/articles/PMC12807593/", "experimental PWM source", "upstream_in_CISBP", "species assignment and native method retained in source archive; not a separate independent replicate"),
    ("ENCODE RBNS", "https://www.encodeproject.org/rbns/", "in-vitro binding specificity / k-mer preference / PWM", "upstream_in_CISBP", "RNA-Bind-n-Seq native method present; mouse assignments may be inferred, NOT recast as direct mouse experiments"),
    ("RBPmap", "https://rbpmap.technion.ac.il/manual.html", "literature consensus/PSSM and weighted-rank prediction", "audited_not_separately_integrated", "native download timed out in web access; no fabricated matrix, no claim of running RBPmap; new method uses FIMO"),
    ("mCrossBase", "https://zhanglab.c2b2.columbia.edu/mCrossBase/", "in-vivo crosslink-registered motifs", "audited_not_scanned", "human/cell-line in-vivo specificity is not silently transferred to mouse; separate matrix/species mapping still needed"),
    ("SpliceAid / SpliceAid2 / SpliceAid-F", "https://pmc.ncbi.nlm.nih.gov/articles/PMC3531144/", "curated human splicing-factor binding sequences / activity", "upstream_partial_not_full_native", "ATtRACT integrates literature resources; human individual sequences are not invented into calibrated mouse PWMs"),
    ("ASD / AEDB", "https://attract.cnic.es/documentation", "curated splicing regulatory elements", "upstream_partial_in_ATtRACT", "AEDB provenance retained; individual functional elements are not independent PWM experiments"),
    ("oRNAment", "https://rnabiology.ircm.qc.ca/oRNAment/", "predicted transcriptome sites from RNAcompete/RBNS PWMs", "audited_not_used_as_CLIP", "official 453 PWMs / 223 RBPs; predictions exclude introns; not independent experimental binding and not exhaustive new-junction annotation"),
    ("MotifMap-RNA", "https://motifmap-rna.ics.uci.edu/", "comparative-genomic predicted RBP binding sites", "audited_not_scanned", "371 motif collection advertised; build/species/native PWM identity not yet matched to GRCm39; not evidence of occupancy"),
    ("ENCODE eCLIP", "https://www.nature.com/articles/s41586-020-2077-3", "reference-cell measured binding peaks and RBP perturbation", "separate_binding_layer_not_computed", "not a universal PWM catalog; previously attached mouse author CLIP remains available in legacy binding workflow"),
    ("POSTAR3 / CLIPdb", "https://pmc.ncbi.nlm.nih.gov/articles/PMC8728292/", "CLIP peaks, structure and post-transcriptional annotation", "audited_not_used_as_PWM", "no complete current species/build-matched local snapshot in this run; gene exports are not genome-wide matched background"),
    ("ENCORI / starBase", "https://rnasysu.com/encori/rbpMotif.php", "CLIP-derived interactions and RBP motif module", "audited_not_scanned", "public PWM/region export exists; native species/build/source identity and nonredundancy not independently bound here"),
    ("MEME RNA collections", "https://meme-suite.org/meme/db/motifs", "redistributed RNA motif libraries", "audited_mirror_not_independent_experiment", "collection overlap with CISBP/RNAcompete/RBPDB; no DNA TF database substitution such as JASPAR/HOCOMOCO"),
]


def write_coverage(directory):
    directory = Path(directory)
    table = pd.DataFrame(RESOURCES, columns=["resource", "primary_url", "evidence_type", "this_run_status", "coverage_and_gap"])
    frozen = json.loads((directory/"ANALYSIS_FREEZE.json").read_text(encoding="utf-8"))
    species = frozen.get("species", "Mus_musculus")
    table["declared_target_species"] = species
    table["review_date"] = "2026-10-15"
    actual = " ".join(str(db.get("asset_id", "")) for db in frozen["databases"]).casefold()
    for label, marker in [("CISBP-RNA", "cisbp"), ("ATtRACT", "attract")]:
        if marker not in actual:
            table.loc[table.resource.eq(label), "this_run_status"] = "not_loaded_in_this_run"
    table.to_csv(directory/"PUBLIC_DATABASE_COVERAGE.tsv", sep="\t", index=False)
    sources = []
    from zipfile import ZipFile
    for db in frozen["databases"]:
        if Path(db["path"]).suffix.lower() != ".zip":
            sources.append({"database": db["asset_id"], "native_method_or_source": "custom numeric PWM resource",
                            "native_assignment": "see declared resource metadata", "native_rows_not_necessarily_usable_PWMs": None})
            continue
        with ZipFile(db["path"]) as z:
            info = next((n for n in z.namelist() if n.endswith("RBP_Information_all_motifs.txt")), None)
            if info:
                meta = pd.read_csv(z.open(info), sep="\t", dtype=str, keep_default_na=False)
                meta = meta[meta.RBP_Species.str.replace(" ", "_").eq(species.replace(" ", "_"))]
                for (method, assignment), group in meta.groupby(["Motif_Type", "RBP_Status"]):
                    sources.append({"database": db["asset_id"], "native_method_or_source": method,
                                    "native_assignment": assignment, "native_rows_not_necessarily_usable_PWMs": len(group)})
            elif "ATtRACT_db.txt" in z.namelist():
                meta = pd.read_csv(z.open("ATtRACT_db.txt"), sep="\t")
                meta = meta[meta.Organism.str.replace(" ", "_").eq(species.replace(" ", "_"))]
                for source, group in meta.groupby("Database"):
                    sources.append({"database": db["asset_id"], "native_method_or_source": source,
                                    "native_assignment": "native metadata before Ensembl/nonmutated/PWM filters",
                                    "native_rows_not_necessarily_usable_PWMs": len(group)})
    pd.DataFrame(sources).to_csv(directory/"AGGREGATED_SOURCE_AUDIT.tsv", sep="\t", index=False)
    return table
