"""Resource binding around the existing, unchanged native Sashimi algorithm."""
from __future__ import annotations

import importlib.util
from pathlib import Path
import re
import shutil
import subprocess
import tempfile

import pandas as pd

from ._execution import PLOT_LOCK, sha256, write_json


GEOMETRY_COLUMNS = (
    "chr", "strand", "exonStart_0base", "exonEnd", "upstreamES", "upstreamEE",
    "downstreamES", "downstreamEE", "riExonStart_0base", "riExonEnd",
    "longExonStart_0base", "longExonEnd", "shortES", "shortEE", "flankingES",
    "flankingEE", "1stExonStart_0base", "1stExonEnd", "2ndExonStart_0base", "2ndExonEnd",
)


def native_reader():
    # bamnostic accepts a separately bound index path. Its installed native
    # adapter reader is reused verbatim, including all original read filters.
    if importlib.util.find_spec("bamnostic") is not None:
        return "bamnostic"
    if importlib.util.find_spec("pysam") is not None:
        return "pysam"
    if shutil.which("samtools"):
        return "samtools"
    return None


def _event_key(event):
    return str(event.get("event_id", event.get("event_uid", event.get("ID", ""))))


def _reference_lengths(reader, bam, index):
    if reader == "bamnostic":
        import bamnostic
        alignment = bamnostic.AlignmentFile(str(bam), "rb", filepath_index=str(index))
        try:
            return dict(zip(alignment.references, alignment.lengths))
        finally:
            alignment.close()
    if reader == "pysam":
        import pysam
        with pysam.AlignmentFile(str(bam), "rb", index_filename=str(index)) as alignment:
            return dict(zip(alignment.references, alignment.lengths))
    header = subprocess.run([shutil.which("samtools"), "view", "-H", str(bam)], check=True,
                            text=True, capture_output=True).stdout
    result = {}
    for line in header.splitlines():
        if line.startswith("@SQ\t"):
            fields = dict(x.split(":", 1) for x in line.split("\t")[1:] if ":" in x)
            if "SN" in fields and "LN" in fields:
                result[fields["SN"]] = int(fields["LN"])
    return result


def generate_sashimi(request, target, events, resources, output_dir):
    """Run exact events and original arm ordering; never download/copy BAMs.

    Missing resources produce explicit unavailable entries, not empty read
    tracks. The return value records every requested event/comparison.
    """
    from src.adapters.native_sashimi_adapter import NativeSashimiAdapter

    destination = Path(output_dir)
    destination.mkdir(parents=True, exist_ok=True)
    samples = {str(s["sample_id"]): s for s in request["project"].get("samples", [])}
    contrasts = {str(c["contrast_id"]): c for c in request["project"].get("contrasts", [])}
    reader = native_reader()
    index_map, sample_map = {}, {}
    repeated_binding = False
    for sample_id, sample in samples.items():
        bam = resources.get(sample.get("bam_resource_id"), {})
        bai = resources.get(sample.get("bai_resource_id"), {})
        if bam.get("available") and bai.get("available"):
            key = str(Path(bam["path"]).resolve())
            if key in sample_map and sample_map[key] != sample_id:
                repeated_binding = True
            index_map[key] = Path(bai["path"])
            sample_map[key] = sample_id

    class BoundNativeSashimi(NativeSashimiAdapter):
        @staticmethod
        def _bam_index_path(bam):
            return index_map.get(str(Path(bam).resolve()))

        def _sample_track(self, bam, geometry, mapper, samtools):
            if reader == "bamnostic":
                track = super()._sample_track_bamnostic(bam, geometry, mapper)
            else:
                # The original pysam/samtools routes discover adjacent indexes.
                expected = {Path(str(bam) + ".bai").resolve(), bam.with_suffix(".bai").resolve()}
                if self._bam_index_path(bam).resolve() not in expected:
                    raise RuntimeError("This bound non-adjacent BAM index requires bamnostic")
                track = super()._sample_track(bam, geometry, mapper, samtools)
            track["sample"] = sample_map[str(Path(bam).resolve())]
            return track

    adapter = BoundNativeSashimi()
    outcome = []
    requested_cids = set(target.get("comparison_ids", []))
    seen = set()
    for event in events:
        cid = str(event.get("comparison_id", event.get("contrast_id", "")))
        eid = _event_key(event)
        key = (eid, cid, event.get("native_statistical_family"), tuple(str(event.get(c, "")) for c in GEOMETRY_COLUMNS))
        if key in seen:
            continue
        seen.add(key)
        if requested_cids and cid not in requested_cids:
            continue
        record = {"event_id": eid, "comparison_id": cid, "event_type": event.get("event_type"),
                  "reader": reader, "status": "unavailable", "files": [],
                  "count_mode": event.get("count_mode", "unknown"),
                  "source_evidence_id": event.get("source_evidence_id"),
                  "native_statistical_family": event.get("native_statistical_family")}
        outcome.append(record)
        if not reader:
            record["reason"] = "No supported BAM reader is installed (bamnostic, pysam or samtools)"
            continue
        if repeated_binding:
            record["reason"] = "One BAM path was bound to different declared biological samples"
            continue
        if request["project"].get("references", {}).get("coordinate_system") != "0-based-half-open":
            record["reason"] = "Native event geometry requires a declared 0-based-half-open coordinate system"
            continue
        contrast = contrasts.get(cid)
        if not contrast:
            record["reason"] = "No exact selected comparison arms are registered"
            continue
        arms = [contrast.get("numerator", {}), contrast.get("denominator", {})]
        ordered_samples = [list(arm.get("sample_ids", [])) for arm in arms]
        unavailable = []
        assay = str(request["project"].get("assay", "bulk")).lower().replace("-", "_")
        single_cell = any(name in assay for name in ("scrna", "snrna", "sc_rna", "sn_rna", "single_cell", "single_nucleus"))
        if single_cell:
            for sid in [sid for ids in ordered_samples for sid in ids]:
                attrs = samples.get(str(sid), {}).get("attributes", {})
                cell_type = attrs.get("cell_type")
                if (attrs.get("bam_scope") != "cell_type_filtered" or not cell_type
                        or attrs.get("bam_cell_type") != cell_type):
                    unavailable.append("Cell-type Sashimi requires a proven cell-type-filtered BAM or the dedicated single-cell Sashimi workflow")
                    break
        bam_lists = []
        for ids in ordered_samples:
            paths = []
            if not ids:
                unavailable.append("Empty comparison arm")
            for sid in ids:
                sample = samples.get(str(sid), {})
                bam_id, bai_id = sample.get("bam_resource_id"), sample.get("bai_resource_id")
                for rid, kind in ((bam_id, "BAM"), (bai_id, "BAM index")):
                    if not rid or not resources.get(rid, {}).get("available"):
                        unavailable.append(f"{kind} unavailable for {sid}")
                if bam_id and resources.get(bam_id, {}).get("available"):
                    paths.append(resources[bam_id]["path"])
            bam_lists.append(paths)
        if unavailable:
            record["reason"] = "; ".join(unavailable)
            continue
        safe_name = re.sub(r"[^A-Za-z0-9_.-]", "_", f"event-{len(outcome):05d}")
        event_dir = destination / safe_name
        event_dir.mkdir()
        native_event = dict(event)
        if not native_event.get("chr"):
            native_event["chr"] = native_event.get("chrom", native_event.get("chromosome", ""))
        if any(str(native_event.get(flag, "")).lower() in {"1", "true", "yes"}
               for flag in ("splicing_direction_flipped", "reverse_direction")):
            native_event["splicing_direction_flipped"] = True
            record["direction_normalization"] = "Native raw IncLevel/IJC/SJC arms swapped to declared comparison; stored effect unchanged"
        if target.get("target_kind") == "event":
            for name, value in target.get("event_geometry", {}).items():
                native_event[name] = value
            native_event["event_type"] = target.get("event_type", native_event.get("event_type"))
        try:
            geometry = adapter._event_geometry(pd.Series(native_event), str(native_event.get("event_type", "")))
            if geometry.strand not in {"+", "-"}:
                raise ValueError("Exact event strand must be declared as + or -")
            headers = []
            for bam in [Path(p) for paths in bam_lists for p in paths]:
                lengths = _reference_lengths(reader, bam, index_map[str(bam.resolve())])
                if geometry.chrom not in lengths:
                    raise ValueError(f"Event chromosome {geometry.chrom} is absent from a bound BAM header")
                if any(exon.start < 1 or exon.end > lengths[geometry.chrom] or exon.start > exon.end for exon in geometry.exons):
                    raise ValueError("Event coordinates are outside a bound BAM reference contig")
                headers.append({"sample_id": sample_map[str(bam.resolve())], "chromosome": geometry.chrom,
                                "contig_length": int(lengths[geometry.chrom])})
            record["reference_preflight"] = headers
            with tempfile.TemporaryDirectory(prefix="katana-target-sashimi-") as scratch:
                scratch = Path(scratch)
                event_path = scratch / "event.tsv"
                pd.DataFrame([native_event]).to_csv(event_path, sep="\t", index=False)
                paths = []
                for i, values in enumerate(bam_lists):
                    path = scratch / f"arm-{i}.txt"
                    # One BAM per line permits commas in actual file names.
                    # The unchanged legacy reader splits commas, so reject those
                    # names rather than silently read a different path.
                    if any("," in p or "\n" in p or "\r" in p for p in values):
                        raise ValueError("Native Sashimi cannot read BAM path names containing commas or newlines")
                    path.write_text("\n".join(values) + "\n", encoding="utf-8")
                    paths.append(path)
                row = pd.Series({**native_event, "event_file": str(event_path),
                                 "event_uid": eid, "event_id": eid,
                                 "GeneID": target.get("gene_id", native_event.get("GeneID", "")),
                                 "geneSymbol": target.get("gene_symbol", native_event.get("geneSymbol", "")),
                                 "comparison_id": cid, "comparison_name": contrast.get("label", cid) + " | " + record["count_mode"],
                                 "label1": arms[0].get("label", "numerator"),
                                 "label2": arms[1].get("label", "denominator"),
                                 "b1_txt": str(paths[0]), "b2_txt": str(paths[1]),
                                 "outdir": str(event_dir.resolve()),
                                 "dPSI": native_event.get("dPSI", native_event.get("dpsi", native_event.get("IncLevelDifference"))),
                                 "FDR": native_event.get("FDR", native_event.get("fdr"))})
                gtf_id = request["project"].get("references", {}).get("gtf_resource_id")
                gtf_info = resources.get(gtf_id, {})
                gtf = Path(gtf_info["path"]) if gtf_info.get("available") else None
                with PLOT_LOCK:
                    adapter._run_one(row, destination, gtf)
            record.update(status="complete", sample_ids=ordered_samples,
                          files=[p.relative_to(destination).as_posix() for p in sorted(event_dir.rglob("*")) if p.is_file()],
                          native_adapter_sha256=sha256(Path(__import__("src.adapters.native_sashimi_adapter", fromlist=["__file__"]).__file__)))
        except Exception as exc:
            record.update(status="failed", reason=f"{type(exc).__name__}: {exc}")
    if not outcome:
        outcome.append({"status": "unavailable", "reason": "No exact supported event geometry and comparison was exported",
                        "reader": reader, "files": []})
    write_json(destination / "sashimi-results.json", outcome)
    return outcome
