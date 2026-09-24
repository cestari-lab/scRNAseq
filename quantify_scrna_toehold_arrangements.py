#!/usr/bin/env python3
# Author: Lissa Cruz-Saavedra
# Date: 24-09-2026


"""
quantify_scrna_toehold_arrangements.py

Quantifies barcode chain completeness for the scRNA-seq toehold/linker
scaffold design: SCAFFOLD_TAIL -> BC1_TOEHOLD -> [BC1] -> BC1_LINKERA ->
BC2_TOEHOLD -> [BC2] -> BC2_LINKERB -> BC3_TOEHOLD -> [BC3] -> polyT/insert.

The scaffold is strictly sequential/ligated in order, so
reads are categorized by HOW FAR the chain was successfully found, not by
which specific round is missing:

    complete            full scaffold->BC1->BC2->BC3 chain found
    stalled_after_BC2   scaffold+BC1+BC2 found, BC3_TOEHOLD not found after
    stalled_after_BC1   scaffold+BC1 found, BC2_TOEHOLD not found after
    stalled_after_scaffold   SCAFFOLD_TAIL found, BC1_TOEHOLD not found after
    no_scaffold         SCAFFOLD_TAIL not found anywhere in the read

Anchor matching allows up to --max-anchor-mismatch mismatches (default 1),
same tolerance as extract_scrna_toehold_barcodes.py, so results are
directly comparable/consistent with that script's "chain_found" rate.

Usage:
    python3 quantify_scrna_toehold_arrangements.py \
        --r1 ScRNA1_R1.fastq.gz \
        --stats-out toehold_arrangement_stats.json \
        --examples-out toehold_arrangement_examples.json \
        [--n-reads 2000000] [--max-anchor-mismatch 1]
"""
import argparse
import gzip
import json
from collections import Counter

SCAFFOLD_TAIL = "GACGCT"
BC1_TOEHOLD   = "TGTAGC"
BC1_LINKERA   = "CTACAG"
BC2_TOEHOLD   = "GGGTAA"
BC2_LINKERB   = "GAGAAC"
BC3_TOEHOLD   = "CTGACT"

BC1_LEN_RANGE = (7, 12)
BC2_LEN_RANGE = (7, 13)
BC3_LEN_RANGE = (6, 12)
UMI_LEN = 10                # fixed UMI between BC3 and poly-T -- was missing
                            # from the polyT search offset in an earlier
                            # version, misplacing chain_end by up to 10bp
POLYT_CHECK_LEN = 8
POLYT_MIN_COUNT = 6


def open_maybe_gz(path, mode="rt"):
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)


def fastq_seqs(path, n_reads=None):
    fh = open_maybe_gz(path)
    i = 0
    for line_no, line in enumerate(fh):
        if line_no % 4 == 1:
            yield line.rstrip("\n")
            i += 1
            if n_reads and i >= n_reads:
                break
    fh.close()


def find_approx(seq, motif, start=0, max_mm=1):
    L = len(motif)
    for p in range(start, len(seq) - L + 1):
        window = seq[p:p + L]
        mm = sum(1 for a, b in zip(window, motif) if a != b)
        if mm <= max_mm:
            return p, mm
    return -1, None


def polyT_score(seq, pos, check_len=POLYT_CHECK_LEN):
    window = seq[pos:pos + check_len]
    if len(window) < check_len:
        return -1
    return window.count("T")


def classify_read(seq, max_mm=1):
    """Returns (category, annotation_dict). Categories, in order of chain
    progress:
        complete_with_polyA   full chain + real poly-A/T tail confirmed
        complete_no_polyA     full chain found, but no real poly-A/T tail
                               right after BC3 (structure present, but not
                               genuine mRNA capture -- e.g. empty bead oligo)
        stalled_after_BC2 / stalled_after_BC1 / stalled_after_scaffold /
        no_scaffold           as before
    """
    i0, _ = find_approx(seq, SCAFFOLD_TAIL, 0, max_mm)
    if i0 == -1:
        return "no_scaffold", None

    i1, _ = find_approx(seq, BC1_TOEHOLD, i0 + len(SCAFFOLD_TAIL), max_mm)
    if i1 == -1:
        return "stalled_after_scaffold", {"i0": i0}
    bc1_start = i1 + len(BC1_TOEHOLD)

    i1b, _ = find_approx(seq, BC1_LINKERA, bc1_start + BC1_LEN_RANGE[0], max_mm)
    if i1b == -1 or i1b - bc1_start > BC1_LEN_RANGE[1]:
        return "stalled_after_scaffold", {"i0": i0, "i1": i1}
    bc1 = seq[bc1_start:i1b]

    i2, _ = find_approx(seq, BC2_TOEHOLD, i1b + len(BC1_LINKERA), max_mm)
    if i2 == -1:
        return "stalled_after_BC1", {"i0": i0, "i1": i1, "i1b": i1b, "bc1": bc1}
    bc2_start = i2 + len(BC2_TOEHOLD)

    i2b, _ = find_approx(seq, BC2_LINKERB, bc2_start + BC2_LEN_RANGE[0], max_mm)
    if i2b == -1 or i2b - bc2_start > BC2_LEN_RANGE[1]:
        return "stalled_after_BC1", {"i0": i0, "i1": i1, "i1b": i1b, "bc1": bc1}
    bc2 = seq[bc2_start:i2b]

    i3, _ = find_approx(seq, BC3_TOEHOLD, i2b + len(BC2_LINKERB), max_mm)
    if i3 == -1:
        return "stalled_after_BC2", {"i0": i0, "i1": i1, "i1b": i1b, "i2": i2,
                                       "i2b": i2b, "bc1": bc1, "bc2": bc2}
    bc3_start = i3 + len(BC3_TOEHOLD)

    best_len, best_score = None, -1
    for L in range(BC3_LEN_RANGE[0], BC3_LEN_RANGE[1] + 1):
        cand_end = bc3_start + L + UMI_LEN   # skip past the UMI before checking polyT
        score = polyT_score(seq, cand_end)
        if score > best_score:
            best_score, best_len = score, L

    if best_len is None:
        return "stalled_after_BC2", {"i0": i0, "i1": i1, "i1b": i1b, "i2": i2,
                                       "i2b": i2b, "bc1": bc1, "bc2": bc2}
    bc3 = seq[bc3_start:bc3_start + best_len]
    umi = seq[bc3_start + best_len: bc3_start + best_len + UMI_LEN]
    chain_end = bc3_start + best_len + UMI_LEN
    ann = {"i0": i0, "i1": i1, "i1b": i1b, "i2": i2, "i2b": i2b, "i3": i3,
           "bc1": bc1, "bc2": bc2, "bc3": bc3, "umi": umi, "chain_end": chain_end,
           "polyT_score": best_score}

    if best_score >= POLYT_MIN_COUNT:
        return "complete_with_polyA", ann
    else:
        return "complete_no_polyA", ann


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--r1", required=True)
    ap.add_argument("--stats-out", required=True)
    ap.add_argument("--examples-out", required=True)
    ap.add_argument("--n-reads", type=int, default=None)
    ap.add_argument("--n-examples", type=int, default=5)
    ap.add_argument("--max-anchor-mismatch", type=int, default=1)
    args = ap.parse_args()

    counts = Counter()
    examples = {}
    total = 0

    for seq in fastq_seqs(args.r1, args.n_reads):
        total += 1
        cat, ann = classify_read(seq, max_mm=args.max_anchor_mismatch)
        counts[cat] += 1
        if len(examples.get(cat, [])) < args.n_examples:
            examples.setdefault(cat, []).append({"seq": seq, "annotation": ann})
        if total % 5_000_000 == 0:
            print(f"  ...{total:,} reads scanned", flush=True)

    stats = {
        "total_reads": total,
        "counts": dict(counts),
        "percent": {k: 100 * v / total for k, v in counts.items()} if total else {},
    }
    with open(args.stats_out, "w") as f:
        json.dump(stats, f, indent=2)
    with open(args.examples_out, "w") as f:
        json.dump(examples, f, indent=2)

    order = ["complete_with_polyA", "complete_no_polyA", "stalled_after_BC2", "stalled_after_BC1",
             "stalled_after_scaffold", "no_scaffold"]
    print(f"\nScanned {total:,} reads.")
    for cat in order:
        if cat in counts:
            print(f"  {cat:24s} {counts[cat]:>10,}  ({100*counts[cat]/total:.3f}%)")
    print(f"\nWrote {args.stats_out} and {args.examples_out}")


if __name__ == "__main__":
    main()
