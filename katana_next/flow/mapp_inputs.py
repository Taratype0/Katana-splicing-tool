"""Explicit SE transcript-TPM bridge for the standalone native MAPP MAEI model.

This is an adapter for existing quantification, not native SUPPA discovery.
Its exact-triplet bridge retains the previously exercised research algorithm;
it does not substitute event PSI or motif display scores for native quantities.
"""
from collections import defaultdict
import gzip
import json
from pathlib import Path
import re

import numpy as np
import pandas as pd

from katana_next.science.isoform_consequences import normalized_id
from katana_cosplice.motif_scanner import IndexedFasta, reverse_complement


COUNTER = "modules/CREATE_SITECOUNT_MATRICES/scripts/create-kmers-sitecount-matrix.py"
BRIDGE_POLICY = "katana.exact-se-transcript-tpm/v1"
APPLICABILITY_SCHEMA = "katana.flow.mapp_applicability/v1"


def validate_inapplicability(path, *, sha, scope=None):
    """Verify a measured zero-eligible universe, never a missing-input skip."""
    proof_path = Path(path)
    proof = json.loads(proof_path.read_text(encoding="utf-8"))
    if (proof.get("schema") != APPLICABILITY_SCHEMA or proof.get("status") != "not_applicable"
            or proof.get("reason_code") != "no_exact_mapped_se"
            or proof.get("complete_measured_universe") is not True
            or proof.get("eligible_events") != 0 or proof.get("native_fitting_performed") is not False):
        raise ValueError("Native MAEI applicability requires a proved complete zero-eligible measured universe")
    audit_path = Path(proof["eligibility_table"])
    if not audit_path.is_file() or sha(audit_path) != proof.get("eligibility_table_sha256"):
        raise ValueError("Native MAEI complete eligibility audit changed")
    audit = pd.read_csv(audit_path, sep="\t")
    allowed = {"no_exact_annotated_included_or_skipped_transcript", "outside_native_exon_length_range",
               "native_min_total_TPM_in_every_sample", "native_uninformative_exon_in_all_samples"}
    if (not {"event_id", "gene_id", "status"} <= set(audit)
            or len(audit) != proof.get("measured_events") or audit.event_id.duplicated().any()
            or audit.status.isna().any() or audit.status.eq("eligible").any()
            or not set(audit.status) <= allowed):
        raise ValueError("Native MAEI eligibility audit does not prove an actual complete zero-eligible universe")
    if not proof.get("source_inputs") or not proof.get("sample_ids") or not proof.get("comparison_id"):
        raise ValueError("Native MAEI applicability needs its actual source, sample and comparison identities")
    for item in proof["source_inputs"]:
        source = Path(item["path"])
        if not source.is_file() or source.stat().st_size != item["size_bytes"] or sha(source) != item["sha256"]:
            raise ValueError("Native MAEI applicability source input changed or is missing")
    if scope:
        if (scope.get("contrast_ids") != [proof["comparison_id"]]
                or any(scope.get(key) != proof.get(key) for key in ("site", "window_start", "window_end"))):
            raise ValueError("Native MAEI applicability belongs to another comparison or native window")
    return proof, audit_path


def native_windows():
    """The existing 14-window family: 50 nt wide, step 25, both splice sites."""
    return [{"site": site, "window_start": start, "window_end": start + 50}
            for site, first, last in (("3ss", -150, 50), ("5ss", -50, 150))
            for start in range(first, last - 49, 25)]


def window_id(start, end):
    if start >= end:
        raise ValueError("MAPP window start must precede end")
    if end <= 0:
        return f"u{-start}to{-end}.d0to0"
    if start >= 0:
        return f"u0to0.d{start}to{end}"
    return f"u{-start}to0.d0to{end}"


def prepare_inputs(spec, out, *, run, sha, executable):
    """Build one comparison/window using all measurable, exactly mapped SEs."""
    inputs, parameters = spec["inputs"], spec.get("parameters") or {}
    samples = list(parameters["sample_ids"])
    if len(samples) < 4 or len(set(samples)) != len(samples):
        raise ValueError("Native MAEI needs at least four distinct declared biological samples")
    if set(samples) != set(inputs["quantification"]):
        raise ValueError("Declare exactly one existing transcript quantification file per selected sample")
    start, end = int(parameters["window_start"]), int(parameters["window_end"])
    name = window_id(start, end)
    site = parameters["site"]
    if site not in {"3ss", "5ss"}:
        raise ValueError("Native MAEI site must be 3ss or 5ss")
    native_counter = Path(inputs["native_root"]) / COUNTER
    if sha(native_counter) != inputs["counter_sha256"]:
        raise ValueError("Official MAPP kmer counter source changed")
    event_path = Path(inputs["events"])
    events = pd.read_parquet(event_path) if event_path.suffix == ".parquet" else pd.read_csv(event_path, sep="\t")
    source_event_rows = len(events)
    comparison = str(parameters.get("comparison_id", ""))
    if "comparison_id" in events:
        if not comparison:
            raise ValueError("Declare the comparison_id; MAPP inputs are never pooled")
        if not events.empty and not events.comparison_id.eq(comparison).any():
            raise ValueError("The requested MAPP comparison has no measured source universe")
        events = events.loc[events.comparison_id.eq(comparison)].copy()
    comparison_event_rows = len(events)
    if "event_type" in events:
        events = events.loc[events.event_type.eq("SE")].copy()
    if "GeneID" in events and "gene_id" not in events:
        events["gene_id"] = events.GeneID
    if "event_id" not in events:
        raise ValueError("Declare canonical exact SE event_id values before the MAEI bridge")
    required = {"event_id", "gene_id", "chr", "strand", "upstreamEE", "exonStart_0base", "exonEnd", "downstreamES"}
    if not required <= set(events) or events.event_id.duplicated().any():
        raise ValueError("MAPP preparation requires one complete exact SE event row per canonical event_id")
    if not events.strand.isin(["+", "-"]).all():
        raise ValueError("MAPP event strands must be + or -")
    events["gene_id"] = events.gene_id.astype(str).str.strip('"').map(normalized_id)
    genes = set(events.gene_id)
    tx_exons, tx_locus = defaultdict(set), {}
    reference_exons = 0
    gtf = Path(inputs["gtf"])
    opener = gzip.open if gtf.suffix == ".gz" else open
    with opener(gtf, "rt", encoding="utf-8") as stream:
        for line in stream:
            if line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            if len(fields) != 9 or fields[2] != "exon":
                continue
            reference_exons += 1
            attributes = dict(re.findall(r'(\S+)\s+"([^"]*)"', fields[8]))
            gene = normalized_id(attributes.get("gene_id", ""))
            if (not gene or not attributes.get("transcript_id") or fields[6] not in {"+", "-"}
                    or int(fields[3]) < 1 or int(fields[4]) < int(fields[3])):
                raise ValueError("MAPP requires a valid strand-oriented exon annotation")
            if gene not in genes:
                continue
            raw_tx = attributes.get("transcript_id", "")
            tx = normalized_id(raw_tx)
            locus = (gene, fields[0], fields[6], raw_tx)
            if not tx or (tx in tx_locus and tx_locus[tx] != locus):
                raise ValueError(f"Ambiguous GTF transcript identity: {raw_tx}")
            tx_locus[tx] = locus
            tx_exons[tx].add((int(fields[3])-1, int(fields[4])))
    if not reference_exons or (genes and not tx_exons):
        raise ValueError("MAPP inputs do not match a measured gene universe in the declared exon annotation")
    included, skipped = defaultdict(set), defaultdict(set)
    for tx, exons in tx_exons.items():
        gene, chrom, strand, _ = tx_locus[tx]
        blocks = sorted(exons)
        for left, right in zip(blocks, blocks[1:]):
            skipped[(gene,chrom,strand,left[1],right[0])].add(tx)
        for left, target, right in zip(blocks, blocks[1:], blocks[2:]):
            included[(gene,chrom,strand,left[1],target[0],target[1],right[0])].add(tx)
    quantities = {}
    for sample in samples:
        source = Path(inputs["quantification"][sample])
        table = pd.read_csv(source, sep="\t")
        tx_column = "Name" if "Name" in table else "transcript_id"
        if tx_column not in table or "TPM" not in table:
            raise ValueError(f"{sample}: existing Salmon/RSEM transcript-level TPM is required")
        raw_ids = table[tx_column].astype(str)
        keys = raw_ids.map(normalized_id)
        if keys.duplicated().any():
            raise ValueError(f"{sample}: ambiguous transcript IDs after version normalization")
        values = pd.to_numeric(table.TPM, errors="raise")
        if not np.isfinite(values).all() or values.lt(0).any():
            raise ValueError(f"{sample}: transcript TPM must be finite and nonnegative")
        for raw_id, key in zip(raw_ids, keys):
            if key in tx_locus:
                annotated = tx_locus[key][3]
                # Where both sources declare versions, they must agree.
                if re.search(r"\.\d+$", raw_id) and re.search(r"\.\d+$", annotated) and raw_id != annotated:
                    raise ValueError(f"{sample}: transcript quantification/GTF version mismatch: {raw_id} / {annotated}")
        quantities[sample] = pd.Series(values.to_numpy(), index=keys)
    fasta = IndexedFasta(Path(inputs["genome_fasta"]))
    if not fasta.index or any(length < 1 or offset < 0 or bases < 1 or width < bases
                              for length,offset,bases,width in fasta.index.values()):
        raise ValueError("MAPP requires a valid, nonempty matched genome FASTA index")
    audit, mapping, expression, sequences = [], [], [], []
    minimum = float(parameters.get("min_total_tpm", 1.0))
    min_length = int(parameters.get("min_exon_length", 100))
    max_length = int(parameters.get("max_exon_length", 1000000))
    with fasta.fasta.open("rb") as stream:
        for row in events.to_dict("records"):
            key = (row["gene_id"], str(row["chr"]), row["strand"], int(row["upstreamEE"]),
                   int(row["exonStart_0base"]), int(row["exonEnd"]), int(row["downstreamES"]))
            inc = included.get(key,set())
            exc = skipped.get((key[0],key[1],key[2],key[3],key[6]),set()) - inc
            total = inc | exc
            status = "eligible"
            if not inc or not exc:
                status = "no_exact_annotated_included_or_skipped_transcript"
            elif not min_length <= key[5]-key[4] <= max_length:
                status = "outside_native_exon_length_range"
            elif any(not total <= set(quantities[sample].index) for sample in samples):
                status = "missing_quantified_transcript_no_zero_filling"
            if status == "eligible":
                i = {s: float(quantities[s].loc[sorted(inc)].sum()) for s in samples}
                t = {s: float(quantities[s].loc[sorted(total)].sum()) for s in samples}
                if any(t[s] < minimum for s in samples):
                    status = "native_min_total_TPM_in_every_sample"
                elif all(i[s] < .00001 for s in samples) or all(t[s]-i[s] < .00001 for s in samples):
                    status = "native_uninformative_exon_in_all_samples"
            if status == "eligible":
                at = (key[4] if key[2] == "+" else key[5]-1) if site == "3ss" else (key[5]-1 if key[2] == "+" else key[4])
                a,b = (at+start,at+end) if key[2] == "+" else (at+1-end,at+1-start)
                sequence = fasta.fetch(key[1], a,b,stream=stream)
                if a < 0 or len(sequence) != end-start:
                    status = "reference_window_not_fully_available"
                elif key[2] == "-":
                    sequence = reverse_complement(sequence)
            audit.append({"event_id": row["event_id"], "gene_id": key[0], "status": status})
            if status != "eligible":
                continue
            expression.append({"coordinates": row["event_id"], **{s+"_included": i[s] for s in samples},
                               **{s+"_total": t[s] for s in samples}})
            mapping.append({"event_id": row["event_id"], "gene_id": key[0], "chr": key[1], "strand": key[2],
                            "window_start0": a, "window_end0": b,
                            "included_transcripts": ",".join(sorted(inc)), "skipped_transcripts": ",".join(sorted(exc))})
            sequences.append((row["event_id"],sequence))
    audit_path = out / "all_event_eligibility.tsv"
    pd.DataFrame(audit, columns=["event_id", "gene_id", "status"]).to_csv(audit_path,sep="\t",index=False)
    if not expression:
        invalid = {row["status"] for row in audit} & {
            "missing_quantified_transcript_no_zero_filling", "reference_window_not_fully_available"}
        if invalid:
            raise ValueError("Native MAEI has incomplete quantified/reference inputs, not scientific inapplicability: " + ", ".join(sorted(invalid)))
        sources = [event_path, gtf, Path(inputs["genome_fasta"]), Path(str(inputs["genome_fasta"])+".fai"),
                   native_counter, *[Path(inputs["quantification"][s]) for s in samples]]
        input_records = [{"path": str(path.resolve()), "size_bytes": path.stat().st_size, "sha256": sha(path)}
                         for path in sources]
        proof = {"schema": APPLICABILITY_SCHEMA, "status": "not_applicable", "reason_code": "no_exact_mapped_se",
                 "complete_measured_universe": True, "eligible_events": 0, "measured_events": len(events),
                 "source_event_rows": source_event_rows, "comparison_event_rows": comparison_event_rows,
                 "comparison_id": comparison, "sample_ids": samples, "site": site,
                 "window_start": start, "window_end": end, "window_id": name,
                 "native_fitting_performed": False, "source_inputs": input_records,
                 "eligibility_table": str(audit_path.resolve()), "eligibility_table_sha256": sha(audit_path),
                 "loss_counts": pd.Series([row["status"] for row in audit], dtype=str).value_counts().to_dict(),
                 "quantity_provenance": {"policy": BRIDGE_POLICY, "not_native_SUPPA_discovery": True,
                    "min_total_tpm_every_sample": minimum, "min_exon_length": min_length,
                    "max_exon_length": max_length, "selection": "all supplied measurable SEs; no P/FDR or motif-based preselection"}}
        proof_path = out / "native_maei_applicability.json"
        proof_path.write_text(json.dumps(proof, indent=2), encoding="utf-8")
        return {"status": "not_applicable", "native_fitting_performed": False,
                "prepared_events": 0, "all_events": len(events), "applicability_proof": str(proof_path),
                "maei_inputs": {"applicability_proof": str(proof_path)}}, [audit_path, proof_path]
    inclusion_path, design_path = out / "INCLUSION_TABLE.tsv", out / "DESIGN.tsv"
    pd.DataFrame(expression).set_index("coordinates").to_csv(inclusion_path,sep="\t")
    pd.DataFrame({"sample": samples}).to_csv(design_path,sep="\t",index=False)
    mapping_path = out / "event_transcript_coordinate_mapping.tsv"
    pd.DataFrame(mapping).to_csv(mapping_path,sep="\t",index=False)
    sequences_path = out / "sequences.fasta"
    with sequences_path.open("x",encoding="ascii",newline="\n") as stream:
        for identity, sequence in sequences:
            stream.write(f">{identity}\n{sequence}\n")
    counts = out / "NATIVE_SITECOUNTS.tsv"
    python = executable(parameters.get("python"), __import__("sys").executable)
    run([python,"-B",str(native_counter.resolve()),"--sequences",str(sequences_path),"--matrix",str(counts),
         "--window-id",name,"--kmer-min","3","--kmer-max","6"],out,"MAPP_native_kmer_counter")
    actual = pd.read_csv(counts,sep="\t",index_col=0)
    if actual.shape != (len(expression),5440) or list(actual.index) != [r["coordinates"] for r in expression]:
        raise ValueError("Official native counter output does not match all 5440 models and exact exon order")
    provenance = {"policy": BRIDGE_POLICY, "not_native_SUPPA_discovery": True,
                  "site": site, "window_id": name, "window_start": start, "window_end": end,
                  "min_total_tpm_every_sample": minimum, "min_exon_length": min_length,
                  "max_exon_length": max_length, "counter_sha256": sha(native_counter),
                  "gtf_sha256": sha(gtf), "events_sha256": sha(event_path),
                  "quantification_sha256": {s: sha(Path(inputs["quantification"][s])) for s in samples},
                  "sample_ids": samples, "comparison_id": comparison,
                  "selection": "all supplied measurable SEs; no P/FDR or motif-based preselection"}
    provenance_path = out / "quantity_coordinate_provenance.json"
    provenance_path.write_text(json.dumps(provenance,indent=2),encoding="utf-8")
    return {"prepared_events": len(expression), "all_events": len(events), "models": 5440,
            "quantity_provenance": provenance, "maei_inputs": {
                "inclusion_table": str(inclusion_path), "sitecount_matrix": str(counts),
                "design": str(design_path), "quantity_provenance": provenance}}, [
                    inclusion_path, counts, design_path, mapping_path, audit_path, sequences_path,
                    provenance_path, out / "MAPP_native_kmer_counter.command.json"]
