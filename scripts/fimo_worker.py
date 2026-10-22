"""Small MEME FIMO worker; only stdlib and pymemesuite, also runnable in WSL."""
import argparse
import json
from pathlib import Path

from pymemesuite.common import MotifFile, Sequence, Background, Alphabet, Array
from pymemesuite.fimo import FIMO
import pymemesuite


def scan(meme, windows, output, threshold):
    data = json.loads(Path(windows).read_text(encoding="utf-8"))
    # PyMEME 0.1.0a4 scans range(length-width), omitting the last valid site.
    # A terminal unscoreable sentinel restores it without changing real bases.
    sequences = [Sequence(row["sequence"] + "N", name=row["window_id"].encode()) for row in data]
    control = "".join(row["sequence"] for row in data if row["role"] == "background")
    counts = [control.count(base)+1 for base in "ACGT"]
    frequencies = [count/sum(counts) for count in counts]
    background = Background(Alphabet.dna(), Array(frequencies))
    rows, tests = [], 0
    # The exported PWM already includes the declared background pseudocount.
    # Pass background explicitly: this reader otherwise defaults to NRDB.
    with MotifFile(meme, pseudocount=0.0) as handle:
        while True:
            motif = handle.read()
            if motif is None:
                break
            # Counting ambiguous sites too only makes global BH conservative.
            possible = sum(max(0, len(row["sequence"]) - motif.width + 1) for row in data)
            tests += possible
            fimo = FIMO(both_strands=False, threshold=threshold, max_stored_scores=max(100, possible + 1))
            pattern = fimo.score_motif(motif, sequences, background, background)
            for hit in pattern.matched_elements:
                rows.append({"motif_token": motif.accession.decode(), "window_id": hit.source.name.decode(),
                             "start_1based": hit.start, "stop_1based": hit.stop,
                             "score": hit.score, "site_pvalue": hit.pvalue})
    Path(output).write_text(json.dumps({"engine": "MEME FIMO via pymemesuite", "version": pymemesuite.__version__,
                                       "tested_positions_all_motifs": tests, "threshold": threshold,
                                       "both_strands": False, "pseudocount": .1,
                                       "background_ACGT": frequencies,
                                       "terminal_site_guard": "unscoreable N sentinel; test counts and genomic coordinates use unpadded sequence",
                                       "sites": rows}), encoding="utf-8")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("meme")
    parser.add_argument("windows")
    parser.add_argument("output")
    parser.add_argument("--threshold", type=float, default=1e-4)
    args = parser.parse_args()
    scan(args.meme, args.windows, args.output, args.threshold)
