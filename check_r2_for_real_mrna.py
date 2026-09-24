#!/usr/bin/env python3
# Author: Lissa Cruz-Saavedra
# Date: 24-09-2026

"""
check_r2_for_real_mrna.py

Tests whether R2 (the cDNA-insert mate) shows evidence of real mRNA capture
for R1 reads in each barcode-completeness category from
quantify_scrna_toehold_arrangements.py -- specifically to check whether
"complete_no_polyA" R1 reads (no poly-T found immediately after BC3) still
have real transcript signal sitting in their R2 mate, which the R1-only
poly-T heuristic would have missed entirely.

For a sample of read pairs, this:
  1. classifies R1 exactly as quantify_scrna_toehold_arrangements.py does
  2. for each category, samples up to --n-per-category R2 mates
  3. reports: R2 poly-A content (a quick, rough signal) and writes out a
     small per-category R2 fastq so you can map each one separately with
     your aligner and directly compare mapping rates -- mapping rate is
     the real, unambiguous test, poly-A content is just a fast preview.

Usage:
    python3 check_r2_for_real_mrna.py \
        --r1 ScRNA1_R1.fastq.gz --r2 ScRNA1_R2.fastq.gz \
        --out-prefix r2_check --n-per-category 20000
"""
import argparse
import gzip
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
POLYT_CHECK_LEN = 8
POLYT_MIN_COUNT = 6


def open_maybe_gz(path, mode="rt"):
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)


def fastq_reader(path):
    fh = open_maybe_gz(path)
    while True:
        name = fh.readline().rstrip("\n")
        if not name:
            break
        seq = fh.readline().rstrip("\n")
        fh.readline()
        qual = fh.readline().rstrip("\n")
        yield name, seq, qual
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


def classify_r1(seq, max_mm=1):
    i0, _ = find_approx(seq, SCAFFOLD_TAIL, 0, max_mm)
    if i0 == -1:
        return "no_scaffold"
    i1, _ = find_approx(seq, BC1_TOEHOLD, i0 + len(SCAFFOLD_TAIL), max_mm)
    if i1 == -1:
        return "stalled_after_scaffold"
    bc1_start = i1 + len(BC1_TOEHOLD)
    i1b, _ = find_approx(seq, BC1_LINKERA, bc1_start + BC1_LEN_RANGE[0], max_mm)
    if i1b == -1 or i1b - bc1_start > BC1_LEN_RANGE[1]:
        return "stalled_after_scaffold"
    i2, _ = find_approx(seq, BC2_TOEHOLD, i1b + len(BC1_LINKERA), max_mm)
    if i2 == -1:
        return "stalled_after_BC1"
    bc2_start = i2 + len(BC2_TOEHOLD)
    i2b, _ = find_approx(seq, BC2_LINKERB, bc2_start + BC2_LEN_RANGE[0], max_mm)
    if i2b == -1 or i2b - bc2_start > BC2_LEN_RANGE[1]:
        return "stalled_after_BC1"
    i3, _ = find_approx(seq, BC3_TOEHOLD, i2b + len(BC2_LINKERB), max_mm)
    if i3 == -1:
        return "stalled_after_BC2"
    bc3_start = i3 + len(BC3_TOEHOLD)
    best_len, best_score = None, -1
    for L in range(BC3_LEN_RANGE[0], BC3_LEN_RANGE[1] + 1):
        score = polyT_score(seq, bc3_start + L)
        if score > best_score:
            best_score, best_len = score, L
    if best_len is None:
        return "stalled_after_BC2"
    return "complete_with_polyA" if best_score >= POLYT_MIN_COUNT else "complete_no_polyA"


def r2_polyA_signal(seq, window=15):
    """Rough, fast check: does R2 START with or contain an early poly-A/T
    run (either strand, since cDNA orientation depends on which strand was
    sequenced)? This is just a quick preview -- mapping is the real test."""
    start_window = seq[:window]
    a_count = max(start_window.count("A"), start_window.count("T"))
    return a_count >= window * 0.7


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--r1", required=True)
    ap.add_argument("--r2", required=True)
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--n-per-category", type=int, default=20000,
                     help="Max read pairs to sample per R1 category")
    ap.add_argument("--n-reads-scan", type=int, default=None,
                     help="Limit total pairs scanned (default: all)")
    args = ap.parse_args()

    r1_reader = fastq_reader(args.r1)
    r2_reader = fastq_reader(args.r2)

    per_cat_counts = Counter()
    per_cat_r2_polyA = Counter()
    out_handles = {}

    total = 0
    for (n1, s1, q1), (n2, s2, q2) in zip(r1_reader, r2_reader):
        total += 1
        if args.n_reads_scan and total > args.n_reads_scan:
            break

        cat = classify_r1(s1)
        if per_cat_counts[cat] >= args.n_per_category:
            continue
        per_cat_counts[cat] += 1

        if r2_polyA_signal(s2):
            per_cat_r2_polyA[cat] += 1

        if cat not in out_handles:
            out_handles[cat] = open_maybe_gz(f"{args.out_prefix}.{cat}.R2.fastq.gz", "wt")
        out_handles[cat].write(f"{n2}\n{s2}\n+\n{q2}\n")

        if total % 5_000_000 == 0:
            print(f"  ...{total:,} pairs scanned", flush=True)

    for h in out_handles.values():
        h.close()

    print(f"\nScanned {total:,} read pairs.\n")
    print(f"{'category':24s} {'n_sampled':>10s} {'R2 quick polyA%':>16s}")
    for cat in ["complete_with_polyA", "complete_no_polyA", "stalled_after_BC2",
                "stalled_after_BC1", "stalled_after_scaffold", "no_scaffold"]:
        n = per_cat_counts.get(cat, 0)
        if n == 0:
            continue
        pct = 100 * per_cat_r2_polyA.get(cat, 0) / n
        print(f"{cat:24s} {n:>10,} {pct:>15.1f}%")

    print(f"\nPer-category R2 fastqs written as {args.out_prefix}.<category>.R2.fastq.gz")
    print("Map each one separately (e.g. with minimap2/bwa against the T. brucei genome)")
    print("and compare mapping rates -- that's the definitive test, not the quick polyA check above.")
    print("Example:")
    print(f"  for f in {args.out_prefix}.*.R2.fastq.gz; do")
    print(f"    minimap2 -ax sr genome.mmi \"$f\" | samtools flagstat -")
    print(f"  done")


if __name__ == "__main__":
    main()
