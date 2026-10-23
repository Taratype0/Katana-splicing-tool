"""Verify an installed downstream runtime using tiny, generated inputs only.

Run outside the checkout with PYTHONPATH unset. This is installation acceptance,
not study validation. It never downloads references or accesses user datasets.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import importlib
import importlib.metadata
import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys


PYTHON_MODULES = (
    "numpy", "pandas", "scipy", "sklearn", "statsmodels", "matplotlib",
    "seaborn", "Bio.Seq", "Bio.SeqIO", "bamnostic", "yaml", "networkx",
    "igraph", "leidenalg", "dynamicTreeCut", "umap", "pydantic", "pyarrow",
    "duckdb", "vl_convert", "bs4", "requests", "jsonschema", "fastapi",
    "uvicorn", "reportlab", "matplotlib_venn", "katana_next.api.app",
    "katana_motiftrace", "katana_cosplice", "katana_next.flow.stages",
)
R_CHECK = r'''
packages <- c('WGCNA','dynamicTreeCut','Rcpp','optparse','DESeq2','DEXSeq',
              'fgsea','limma','NetRep','dplyr','magrittr','tidyr','jsonlite',
              'GO.db','GenomeInfoDbData')
stopifnot(all(vapply(packages, requireNamespace, logical(1), quietly=TRUE)))
library_root <- paste0(normalizePath(R.home('library')), .Platform$file.sep)
package_paths <- vapply(packages, function(x)normalizePath(find.package(x)), '')
stopifnot(all(startsWith(package_paths, library_root)))
set.seed(1729)
x <- matrix(rnorm(160), nrow=8)
adj <- WGCNA::adjacency(x, type='signed', power=2)
tom <- WGCNA::TOMsimilarity(adj, verbose=0)
stopifnot(all(is.finite(tom)), identical(dim(tom), c(20L,20L)))
Rcpp::cppFunction('double katana_sum(Rcpp::NumericVector x) { return Rcpp::sum(x); }')
stopifnot(katana_sum(c(1,2,3)) == 6)
counts <- matrix(rnbinom(800, mu=80, size=10),nrow=100,
                 dimnames=list(paste0('g',1:100),paste0('s',1:8)))
meta <- data.frame(condition=factor(rep(c('A','B'),each=4)),row.names=colnames(counts))
dds <- DESeq2::DESeqDataSetFromMatrix(counts,meta,~condition)
dds <- DESeq2::DESeq(dds,fitType='mean',quiet=TRUE)
res <- as.data.frame(DESeq2::results(dds))
stopifnot(nrow(res)==100,any(is.finite(res$pvalue)))
dx <- DEXSeq::DEXSeqDataSet(counts,sampleData=meta,design=~sample+exon+condition:exon,
           featureID=rep(paste0('E',1:2),50),groupID=rep(paste0('g',1:50),each=2))
dx <- DESeq2::estimateSizeFactors(dx)
stopifnot(all(is.finite(DESeq2::sizeFactors(dx))))
rank <- setNames(rnorm(100),paste0('g',1:100))
fg <- fgsea::fgseaSimple(list(set1=paste0('g',1:10)),rank,nperm=100,nproc=1)
stopifnot(nrow(fg)==1,is.finite(fg$pval))
jsonlite::write_json(list(status='passed',R=R.version.string,
     versions=as.list(setNames(vapply(packages,function(x)as.character(packageVersion(x)),''),packages)),
     package_paths=as.list(package_paths),
     checks=c('WGCNA TOM','Rcpp native compilation','DESeq2 tiny fit',
              'DEXSeq constructor and size factors','fgsea tiny permutation')),
     'r-check.json',auto_unbox=TRUE,pretty=TRUE)
'''


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--portable-only", action="store_true",
                        help="Check Python/browser/CLI only; does not certify native Linux tools")
    args = parser.parse_args()
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    os.environ.setdefault("MPLBACKEND", "Agg")
    os.environ["MPLCONFIGDIR"] = str(out / "matplotlib-config")
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
    report = {"schema": "katana.runtime_acceptance/v1", "started_at": datetime.now(timezone.utc).isoformat(),
              "python": sys.version, "prefix": sys.prefix, "executable": sys.executable,
              "scope": "portable_python" if args.portable_only else "linux_downstream_runtime",
              "study_analysis": False, "checks": [], "failures": []}

    def check(name, operation):
        try:
            detail = operation()
            report["checks"].append({"name": name, "status": "passed", "detail": detail})
        except Exception as exc:
            report["failures"].append({"name": name, "error": str(exc)})

    def run(name, argv, timeout=180):
        result = subprocess.run(argv, cwd=out, text=True, capture_output=True, timeout=timeout)
        (out / (name + ".stdout.log")).write_text(result.stdout, encoding="utf-8")
        (out / (name + ".stderr.log")).write_text(result.stderr, encoding="utf-8")
        if result.returncode:
            raise RuntimeError(f"exit {result.returncode}: {result.stderr[-1500:]}")
        return {"argv": list(map(str, argv)), "returncode": result.returncode}

    def installed_source():
        import katana_next
        path = Path(katana_next.__file__).resolve()
        if not path.is_relative_to(Path(sys.prefix).resolve()):
            raise RuntimeError(f"Katana resolves outside the installed environment: {path}")
        if importlib.util.find_spec("PySide6") is not None:
            raise RuntimeError("Headless acceptance unexpectedly has PySide6 installed")
        root = path.parent
        required = [root / "science" / name for name in (
            "expression_response.R", "dexseq_interaction.R", "hub_metrics.R", "djexpress_driver.R")]
        required += [root / "vendor/DJExpress/R/JCNA1pass.R",
                     root / "flow/_pipeline/main.nf"]
        missing = [str(p) for p in required if not p.is_file()]
        if missing:
            raise RuntimeError(f"Wheel runtime resources missing: {missing}")
        return {"package": str(path), "resources": [str(p.relative_to(root)) for p in required]}

    check("installed_headless_wheel", installed_source)
    for module in PYTHON_MODULES:
        check("import:" + module, lambda module=module: str(importlib.import_module(module).__file__))
    for command in ("katana-flow", "katana-report", "katana-targets", "katana-motiftrace",
                    "katana-splice-network", "katana-contracts", "katana-open"):
        check("cli:" + command, lambda command=command: run(command, [command, "--help"]))
    check("pip_check", lambda: run("pip-check", [sys.executable, "-m", "pip", "check"]))

    def graphics():
        import matplotlib.pyplot as plt
        import vl_convert
        fig, ax = plt.subplots()
        ax.plot([0, 1], [0, 1]); fig.savefig(out / "matplotlib.png"); plt.close(fig)
        svg = vl_convert.vegalite_to_svg({"$schema": "https://vega.github.io/schema/vega-lite/v5.json",
                 "data": {"values": [{"x": 1, "y": 2}]}, "mark": "point",
                 "encoding": {"x": {"field": "x", "type": "quantitative"},
                              "y": {"field": "y", "type": "quantitative"}}})
        (out / "vega.svg").write_text(svg, encoding="utf-8")
        return {"matplotlib_bytes": (out / "matplotlib.png").stat().st_size, "vega_bytes": len(svg)}

    check("headless_graphics", graphics)
    if not args.portable_only:
        def native_paths():
            paths = {name: shutil.which(name) for name in ("Rscript", "samtools", "geneBody_coverage.py")}
            outside = [name for name, path in paths.items() if not path or
                       not Path(path).resolve().is_relative_to(Path(sys.prefix).resolve())]
            if outside:
                raise RuntimeError(f"Native tools outside this installed environment: {outside}")
            return paths

        check("native_tools_in_environment", native_paths)
        for module in ("pysam", "pymemesuite.fimo", "pyhmmer"):
            check("import:" + module, lambda module=module: str(importlib.import_module(module).__file__))
        (out / "r-check.R").write_text(R_CHECK, encoding="utf-8")
        check("R_native", lambda: run("r-native", ["Rscript", "--vanilla", str(out / "r-check.R")], 300))

        def mapp_compile():
            from katana_next.flow.bundled_runtime import verified_mapp_sources
            source, descriptor = verified_mapp_sources()
            copied = out / "native-mapp"
            shutil.copytree(source / "modules/MAEI/scripts", copied)
            script = out / "compile-mapp.R"
            script.write_text("Rcpp::sourceCpp(commandArgs(TRUE)[1], rebuild=TRUE)\n"
                              "stopifnot(is.function(fit_model_parameters))\n"
                              "cat('MAPP_NATIVE_MODEL_COMPILED\\n')\n", encoding="utf-8")
            command = run("mapp-native-compile", ["Rscript", "--vanilla", str(script),
                          str(copied / "model.cpp")], 300)
            # Both compilation input and original packaged source must retain
            # the exact author fingerprints after Rcpp's private build.
            _, after = verified_mapp_sources()
            if descriptor != after:
                raise RuntimeError("Packaged MAPP identity changed during private compilation")
            for relative, digest in descriptor["source_sha256"].items():
                if hashlib.sha256((copied / Path(relative).name).read_bytes()).hexdigest() != digest:
                    raise RuntimeError("Private MAPP source changed during compilation")
            return {"native_commit": descriptor["native_commit"], "compilation": command,
                    "scientific_fit": False}

        check("installed_MAPP_native_compilation", mapp_compile)

        def network_bridge():
            import numpy as np
            import pandas as pd
            from katana_cosplice.config import NetworkConfig
            from katana_cosplice.network import fit_network
            rng = np.random.default_rng(1729)
            samples = [f"synthetic_{i}" for i in range(6)]
            signals = np.array([[.15, .20, .25, .70, .75, .80],
                                [.20, .65, .30, .80, .40, .60]])
            matrix = pd.DataFrame(np.vstack([signal + rng.normal(0, .015, (30, 6))
                                              for signal in signals]),
                                  index=[f"event_{i}" for i in range(60)], columns=samples)
            metrics = pd.DataFrame({"stability_score": 1.0}, index=matrix.index)
            metadata = pd.DataFrame({"sample_id": samples, "condition": ["A"] * 3 + ["B"] * 3})
            settings = NetworkConfig(minimum_network_samples=6,
                                     wgcna_rscript=shutil.which("Rscript"))
            result = fit_network(matrix, metrics, settings, samples=metadata)
            if len(result.assignments) != 60 or result.linkage_matrix.shape != (59, 4):
                raise RuntimeError("Installed native network bridge lost events or linkage rows")
            result.assignments.to_csv(out / "synthetic-network-assignments.tsv", sep="\t", index=False)
            return {"samples": 6, "events": 60, "linkage_rows": 59,
                    "interpretation": "synthetic installation check; six independent samples"}

        check("installed_WGCNA_network_bridge", network_bridge)

        def fimo():
            from scripts.fimo_worker import scan
            motif = out / "tiny.meme"
            motif.write_text("MEME version 4\n\nALPHABET= ACGT\n\nstrands: +\n\n"
                 "Background letter frequencies\nA 0.25 C 0.25 G 0.25 T 0.25\n\n"
                 "MOTIF test_ACG\nletter-probability matrix: alength= 4 w= 3 nsites= 10 E= 0\n"
                 "0.97 0.01 0.01 0.01\n0.01 0.97 0.01 0.01\n0.01 0.01 0.97 0.01\n", encoding="utf-8")
            windows = out / "tiny-windows.json"
            windows.write_text(json.dumps([{"window_id": "target", "sequence": "ACGACGACG", "role": "target"},
                 {"window_id": "control", "sequence": "ATGCTAGCT", "role": "background"}]), encoding="utf-8")
            scan(str(motif), str(windows), str(out / "fimo.json"), 0.1)
            result = json.loads((out / "fimo.json").read_text())
            if not result["sites"] or result["tested_positions_all_motifs"] != 14:
                raise RuntimeError("Synthetic FIMO scan produced unexpected site/test results")
            return {"version": result["version"], "hits": len(result["sites"]), "tests": 14}

        check("native_FIMO_worker", fimo)

        def hmmer():
            import pyhmmer
            alphabet = pyhmmer.easel.Alphabet.amino()
            sequence = pyhmmer.easel.TextSequence(
                name=b"synthetic_protein", sequence="ACDEFGHIKLMNPQRSTVWY" * 3).digitize(alphabet)
            model, _, _ = pyhmmer.plan7.Builder(alphabet).build(
                sequence, pyhmmer.plan7.Background(alphabet))
            hits = list(pyhmmer.hmmer.hmmsearch(
                model, pyhmmer.easel.DigitalSequenceBlock(alphabet, [sequence]), cpus=1))
            if len(hits) != 1 or len(hits[0]) != 1:
                raise RuntimeError("Synthetic HMM search did not recover its exact training sequence")
            return {"version": pyhmmer.__version__, "hits": len(hits[0]),
                    "scope": "native engine check; no Pfam resource or biological inference"}

        check("native_PyHMMER", hmmer)

        def bam_qc():
            import pysam
            bam = out / "tiny.bam"
            header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chrTiny", "LN": 2000}]}
            with pysam.AlignmentFile(str(bam), "wb", header=header) as handle:
                for i in range(1000):
                    read = pysam.AlignedSegment(); read.query_name = f"synthetic_{i}"
                    read.query_sequence = "A" * 75; read.flag = 0; read.reference_id = 0
                    read.reference_start = 100 + i; read.mapping_quality = 60
                    read.cigar = [(0, 75)]; read.query_qualities = pysam.qualitystring_to_array("I" * 75)
                    handle.write(read)
            pysam.index(str(bam))
            run("samtools", ["samtools", "quickcheck", "-v", str(bam)])
            bed = out / "tiny.bed"
            bed.write_text("chrTiny\t100\t1200\ttiny_gene\t0\t+\t100\t1200\t0\t1\t1100\t0\n")
            run("rseqc", ["geneBody_coverage.py", "-i", str(bam), "-r", str(bed), "-o", str(out / "coverage")])
            table = out / "coverage.geneBodyCoverage.txt"
            if not table.is_file() or table.stat().st_size == 0:
                raise RuntimeError("RSeQC did not create a coverage table")
            return {"synthetic_reads": 1000, "coverage_bytes": table.stat().st_size}

        check("native_BAM_RSeQC", bam_qc)
    report["finished_at"] = datetime.now(timezone.utc).isoformat()
    report["status"] = "passed" if not report["failures"] else "failed"
    report["files"] = [{"name": p.name, "size": p.stat().st_size,
                         "sha256": hashlib.sha256(p.read_bytes()).hexdigest()}
                        for p in sorted(out.iterdir()) if p.is_file()]
    (out / "acceptance.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps({"status": report["status"], "checks": len(report["checks"]),
                      "failures": report["failures"], "receipt": str(out / "acceptance.json")}, indent=2))
    return 1 if report["failures"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
