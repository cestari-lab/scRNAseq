#!/usr/bin/env python3
# Author: Lissa Cruz-Saavedra
# Date: 24-09-2026

"""
highlight_illumina_reads.py

Finds the scaffold->BC1->BC2->BC3 connector chain in Illumina R1 reads and
produces colored output (HTML + ANSI terminal), same color scheme as the
Nanopore version:
    green  = SCAFFOLD_TAIL + BC1_TOEHOLD
    red    = BC1_LINKERA + BC2_TOEHOLD
    blue   = BC2_LINKERB + BC3_TOEHOLD
    orange = poly-T run after BC3

Illumina data matches the matching logic in
quantify_scrna_toehold_arrangements.py / extract_scrna_toehold_barcodes.py,
already validated against the real production run:
  1. Fixed-window Hamming (substitution-only) matching, not edit-distance --
     Illumina's error profile is dominated by substitutions, not indels,
     so Hamming matching is both sufficient and much faster at this scale.
  2. Forward orientation only -- R1's orientation is fixed by the
     sequencing primer, unlike Nanopore reads which can come from either
     strand.
  3. Searches near the start of the read only (chain begins around
     position 0), not the whole read -- R1's barcode is always at a
     predictable position, unlike a long Nanopore read where it could be
     anywhere.

Usage:
    python3 highlight_illumina_reads.py \
        --fastq ScRNA1_R1.fastq.gz \
        --out-html highlighted_illumina_reads.html \
        --n-reads 100 --max-mismatch 1
"""
import argparse
import gzip
import html as htmlmod

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
POLYT_SEARCH_WINDOW = 10

ANSI = {"green": "\033[92m", "red": "\033[91m", "blue": "\033[94m",
        "orange": "\033[93m", "reset": "\033[0m"}
HTML_COLOR = {"green": "#2ca02c", "red": "#d62728", "blue": "#1f77b4", "orange": "#e69138"}


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


def find_approx(seq, motif, start=0, max_mm=1, max_gap=0):
    """Fixed-window Hamming search: first position in [start, start+max_gap]
    where motif matches within max_mm substitutions. Returns (pos, end) or
    None. max_gap=0 means the match must start EXACTLY at `start` (tight
    adjacency to the previous anchor) -- this matters because constant
    flanking sequence (e.g. the P5 adapter) can coincidentally contain a
    near-match to a short anchor motif; without enforcing adjacency, a
    false early match can be picked up and the chain stretched all the way
    to the true, correctly-spaced anchors further downstream."""
    L = len(motif)
    for p in range(start, min(start + max_gap, len(seq) - L) + 1):
        window = seq[p:p + L]
        mm = sum(1 for a, b in zip(window, motif) if a != b)
        if mm <= max_mm:
            return p, p + L
    return None


def find_chain(seq, max_mm=1):
    """Forward orientation only. Tries successive SCAFFOLD_TAIL candidates
    (in case an earlier one is a coincidental false match, e.g. inside the
    P5 adapter's own constant sequence) until one yields a FULL chain with
    every subsequent anchor tightly adjacent (no gap, beyond the allowed
    variable barcode length) -- this avoids both false matches (stretching
    across real sequence) and silently giving up on the first candidate."""
    search_from = 0
    while True:
        r0 = find_approx(seq, SCAFFOLD_TAIL, search_from, max_mm, max_gap=len(seq))
        if r0 is None:
            return None
        i0, i0_end = r0

        r1 = find_approx(seq, BC1_TOEHOLD, i0_end, max_mm, max_gap=0)
        if r1 is None:
            search_from = i0 + 1
            continue
        i1, i1_end = r1
        bc1_start = i1_end

        r1b = find_approx(seq, BC1_LINKERA, bc1_start + BC1_LEN_RANGE[0], max_mm,
                           max_gap=BC1_LEN_RANGE[1] - BC1_LEN_RANGE[0])
        if r1b is None:
            search_from = i0 + 1
            continue
        i1b, i1b_end = r1b
        bc1 = seq[bc1_start:i1b]

        r2 = find_approx(seq, BC2_TOEHOLD, i1b_end, max_mm, max_gap=0)
        if r2 is None:
            search_from = i0 + 1
            continue
        i2, i2_end = r2
        bc2_start = i2_end

        r2b = find_approx(seq, BC2_LINKERB, bc2_start + BC2_LEN_RANGE[0], max_mm,
                           max_gap=BC2_LEN_RANGE[1] - BC2_LEN_RANGE[0])
        if r2b is None:
            search_from = i0 + 1
            continue
        i2b, i2b_end = r2b
        bc2 = seq[bc2_start:i2b]

        r3 = find_approx(seq, BC3_TOEHOLD, i2b_end, max_mm, max_gap=0)
        if r3 is None:
            search_from = i0 + 1
            continue
        i3, i3_end = r3
        bc3_start = i3_end

        best_len, best_score = None, -1
        for L in range(BC3_LEN_RANGE[0], BC3_LEN_RANGE[1] + 1):
            cand_end = bc3_start + L
            for offset in range(POLYT_SEARCH_WINDOW):
                window = seq[cand_end + offset:cand_end + offset + POLYT_CHECK_LEN]
                if len(window) < POLYT_CHECK_LEN:
                    continue
                score = window.count("T")
                if score > best_score:
                    best_score, best_len = score, L

        if best_len is None:
            best_len = 9
        bc3 = seq[bc3_start:bc3_start + best_len]
        chain_end = bc3_start + best_len

        return {
            "i0": i0, "i0_end": i0_end, "i1": i1, "i1_end": i1_end,
            "i1b": i1b, "i1b_end": i1b_end, "i2": i2, "i2_end": i2_end,
            "i2b": i2b, "i2b_end": i2b_end, "i3": i3, "i3_end": i3_end,
            "bc1": bc1, "bc2": bc2, "bc3": bc3, "chain_end": chain_end,
            "polyT_score": best_score,
        }


def build_colored_html(seq, ann):
    i0, i1_end = ann["i0"], ann["i1_end"]
    i1b, i2_end = ann["i1b"], ann["i2_end"]
    i2b, i3_end = ann["i2b"], ann["i3_end"]
    ce = ann["chain_end"]

    spans = [(0, i0, None), (i0, i1_end, "green"), (i1_end, i1b, None),
             (i1b, i2_end, "red"), (i2_end, i2b, None), (i2b, i3_end, "blue"),
             (i3_end, ce, None)]

    trail = seq[ce:min(ce + 200, len(seq))]
    j = 0
    while j < len(trail) and trail[j] != "T":
        j += 1
    k = j
    while k < len(trail) and trail[k] == "T":
        k += 1
    spans.append((ce, ce + j, None))
    if k > j:
        spans.append((ce + j, ce + k, "orange"))
    spans.append((ce + k, len(seq), None))

    html_parts, ansi_parts = [], []
    for s, e, color in spans:
        if e <= s:
            continue
        text = seq[s:e]
        esc = htmlmod.escape(text)
        if color:
            html_parts.append(f'<span style="color:{HTML_COLOR[color]}">{esc}</span>')
            ansi_parts.append(f"{ANSI[color]}{text}{ANSI['reset']}")
        else:
            html_parts.append(esc)
            ansi_parts.append(text)
    return "".join(html_parts), "".join(ansi_parts)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--fastq", required=True)
    ap.add_argument("--out-html", required=True)
    ap.add_argument("--n-reads", type=int, default=100,
                     help="Stop after finding this many reads with a valid chain")
    ap.add_argument("--n-scan-max", type=int, default=None,
                     help="Cap total reads scanned looking for matches (default: unlimited)")
    ap.add_argument("--max-mismatch", type=int, default=1,
                     help="Max mismatches per anchor (substitutions only)")
    ap.add_argument("--min-polyt-score", type=int, default=0,
                     help="Only keep reads with a poly-T score >= this (out of 8). "
                          "Use 6+ to see only reads with a REAL poly-T tail, not just "
                          "complete barcode structure with no confirmed mRNA capture.")
    args = ap.parse_args()

    found = []
    scanned = 0
    for seq in fastq_seqs(args.fastq):
        scanned += 1
        ann = find_chain(seq, max_mm=args.max_mismatch)
        if ann is not None and ann["polyT_score"] >= args.min_polyt_score:
            found.append((seq, ann))
        if len(found) >= args.n_reads:
            break
        if args.n_scan_max and scanned >= args.n_scan_max:
            break

    print(f"Scanned {scanned:,} reads, found chain in {len(found)}.")

    html_blocks = []
    for seq, ann in found:
        colored_html, colored_ansi = build_colored_html(seq, ann)
        html_blocks.append(
            f'<p style="margin:0 0 1.5rem;font-family:monospace;font-size:13px;'
            f'word-break:break-all;">{colored_html}</p>'
        )
        print(colored_ansi)

    with open(args.out_html, "w") as f:
        f.write("<html><body style='background:white;color:#111;padding:20px;'>\n")
        f.write("\n".join(html_blocks))
        f.write("\n</body></html>\n")

    print(f"\nWrote {args.out_html}")


if __name__ == "__main__":
    main()
