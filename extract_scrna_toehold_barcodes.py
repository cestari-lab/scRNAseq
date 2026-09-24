#!/usr/bin/env python3
# Author: Lissa Cruz-Saavedra
# Date: 24-10-2026
"""
extract_scrna_toehold_barcodes.py

Extracts the scaffold -> BC1 -> BC2 -> BC3 -> polyT structure from scRNA-seq
reads built on the toehold/linker bead-oligo design (NOT the GATC/Sau3AI
combinatorial scheme used elsewhere in this project for Hi-C -- these are
two different barcode designs, confirmed empirically on FV74MG data).

Structure (5'->3'):
    [P5 adapter] SCAFFOLD_TAIL [BC1] BC1_LINKERA [BC2] BC2_LINKERB [BC3] [polyT/insert]

Anchors are matched with up to --max-anchor-mismatch mismatches (default 1)
at every position in a search window, not just exact substring search --
important since real reads have sequencing errors and exact-only matching
undercounts real structure (as seen when a quick exact-only check found
~14-19% of reads with structure; this tolerant version should find more).

If --plates is given (an .xlsx with BC1/BC2/BC3 well sequences, one sheet
per round), extracted barcodes are matched against the real well sequences
(Hamming <=1) and reads get a real cell_id. Without it, raw extracted
barcode strings are still tallied and reported -- recurring identical
strings across many reads are themselves evidence of real, discrete
barcode identities, just not yet mapped to a physical well.

Usage:
    python3 extract_scrna_toehold_barcodes.py \
        --r1 ScRNA1_R1.fastq.gz \
        --out-prefix scrna_toehold \
        [--plates Primers_singel_cell.xlsx] \
        [--n-reads 2000000] [--max-anchor-mismatch 1] [--min-insert-len 20]
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

BC1_LEN_RANGE = (7, 12)   # observed variable-region length range in real reads
BC2_LEN_RANGE = (7, 13)
BC3_LEN_RANGE = (6, 12)    # BC3 is ALSO variable length (was wrongly fixed at 9
                            # in an earlier version -- that misplaced chain_end
                            # inside the adapter for most reads). Correct length
                            # is chosen as whichever gives the best poly-T score
                            # immediately after it.
POLYT_CHECK_LEN = 8
POLYT_MIN_COUNT = 6         # require >=6/8 T's right after BC3 to call it real
                            # poly-A capture, not just barcode structure

COMP = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}


def revcomp(s):
    return "".join(COMP[b] for b in reversed(s))


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
    """Find the first position >= start where `motif` matches within
    max_mm mismatches. Returns (pos, mismatches) or (-1, None)."""
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


def extract_chain(seq, max_mm=1, require_polyA=True):
    """Returns a dict with extracted BC1/BC2/BC3 and anchor positions, or
    None if the full chain (and, if require_polyA, a real poly-A/T tail)
    isn't found. BC3's length is chosen from BC3_LEN_RANGE as whichever
    gives the best poly-T score immediately after it -- BC3 is variable
    length, not fixed, so this also fixes correct chain_end placement even
    when require_polyA=False."""
    i0, _ = find_approx(seq, SCAFFOLD_TAIL, 0, max_mm)
    if i0 == -1:
        return None
    i1, _ = find_approx(seq, BC1_TOEHOLD, i0 + len(SCAFFOLD_TAIL), max_mm)
    if i1 == -1:
        return None
    bc1_start = i1 + len(BC1_TOEHOLD)

    i1b, _ = find_approx(seq, BC1_LINKERA, bc1_start + BC1_LEN_RANGE[0], max_mm)
    if i1b == -1 or i1b - bc1_start > BC1_LEN_RANGE[1]:
        return None
    bc1 = seq[bc1_start:i1b]

    i2, _ = find_approx(seq, BC2_TOEHOLD, i1b + len(BC1_LINKERA), max_mm)
    if i2 == -1:
        return None
    bc2_start = i2 + len(BC2_TOEHOLD)

    i2b, _ = find_approx(seq, BC2_LINKERB, bc2_start + BC2_LEN_RANGE[0], max_mm)
    if i2b == -1 or i2b - bc2_start > BC2_LEN_RANGE[1]:
        return None
    bc2 = seq[bc2_start:i2b]

    i3, _ = find_approx(seq, BC3_TOEHOLD, i2b + len(BC2_LINKERB), max_mm)
    if i3 == -1:
        return None
    bc3_start = i3 + len(BC3_TOEHOLD)

    best_len, best_score = None, -1
    for L in range(BC3_LEN_RANGE[0], BC3_LEN_RANGE[1] + 1):
        cand_end = bc3_start + L
        score = polyT_score(seq, cand_end)
        if score > best_score:
            best_score, best_len = score, L

    if require_polyA and best_score < POLYT_MIN_COUNT:
        return None
    if best_len is None:
        return None
    bc3 = seq[bc3_start:bc3_start + best_len]

    return {
        "bc1": bc1, "bc2": bc2, "bc3": bc3,
        "polyT_score": best_score,
        "chain_end": bc3_start + best_len,   # insert/polyT starts here
    }


def load_plate_lookup(xlsx_path, sheet, prefix_len=6, suffix_len=6, umi_marker=None):
    import openpyxl
    wb = openpyxl.load_workbook(xlsx_path, data_only=True)
    ws = wb[sheet]
    entries = {}
    for row in ws.iter_rows(min_row=3, values_only=True):
        well, name, seq = row[0], row[1], row[2]
        if not seq:
            continue
        seq = seq.replace("/5Phos/", "").replace(" ", "")
        if umi_marker and umi_marker in seq:
            rest = seq[prefix_len:]
            bc = rest[:rest.find(umi_marker)]
        else:
            bc = seq[prefix_len:-suffix_len] if suffix_len else seq[prefix_len:]
        entries[bc] = name
    return entries


def hamming_leq(a, b, max_mm):
    if len(a) != len(b):
        return None
    mm = sum(1 for x, y in zip(a, b) if x != y)
    return mm if mm <= max_mm else None


def match_well(bc, lookup, max_mm=1):
    if bc in lookup:
        return lookup[bc]
    best, best_mm = None, max_mm + 1
    for k, v in lookup.items():
        mm = hamming_leq(bc, k, max_mm)
        if mm is not None and mm < best_mm:
            best, best_mm = v, mm
    return best


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--r1", required=True)
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--plates", default=None,
                     help="Optional .xlsx with BC1/BC2/BC3 sheets for real well matching")
    ap.add_argument("--bc3-sheet", default="BC3-PT")
    ap.add_argument("--n-reads", type=int, default=None)
    ap.add_argument("--max-anchor-mismatch", type=int, default=1)
    ap.add_argument("--min-insert-len", type=int, default=20)
    ap.add_argument("--require-polyA", action="store_true", default=True,
                     help="Only keep reads with a real poly-A/T tail right after BC3 "
                          "(default: on). Use --no-require-polyA to disable.")
    ap.add_argument("--no-require-polyA", dest="require_polyA", action="store_false")
    args = ap.parse_args()

    bc1_lookup = bc2_lookup = bc3_lookup = None
    if args.plates:
        bc1_lookup = load_plate_lookup(args.plates, "BC1")
        bc2_lookup = load_plate_lookup(args.plates, "BC2")
        bc3_lookup = load_plate_lookup(args.plates, args.bc3_sheet, umi_marker="N")
        print(f"Loaded plate lookups: BC1={len(bc1_lookup)}, BC2={len(bc2_lookup)}, BC3={len(bc3_lookup)} wells")

    total = 0
    n_chain_found = 0
    n_insert_kept = 0
    combo_counts = Counter()
    well_combo_counts = Counter()

    out_fastq = open_maybe_gz(f"{args.out_prefix}.tagged_insert.fastq.gz", "wt")

    for seq in fastq_seqs(args.r1, args.n_reads):
        total += 1
        chain = extract_chain(seq, max_mm=args.max_anchor_mismatch, require_polyA=args.require_polyA)
        if chain is None:
            continue
        n_chain_found += 1

        raw_combo = f"{chain['bc1']}_{chain['bc2']}_{chain['bc3']}"
        combo_counts[raw_combo] += 1

        cell_id = raw_combo
        if bc1_lookup:
            w1 = match_well(chain["bc1"], bc1_lookup)
            w2 = match_well(chain["bc2"], bc2_lookup)
            w3 = match_well(chain["bc3"], bc3_lookup)
            if w1 and w2 and w3:
                cell_id = f"{w1}_{w2}_{w3}"
                well_combo_counts[cell_id] += 1

        insert = seq[chain["chain_end"]:]
        # trim leading poly-T (mRNA capture) before keeping the rest as insert
        j = 0
        while j < len(insert) and insert[j] == "T":
            j += 1
        insert_trimmed = insert[j:]
        if len(insert_trimmed) >= args.min_insert_len:
            n_insert_kept += 1
            out_fastq.write(f"@read{total} CB:Z:{cell_id}\n{insert_trimmed}\n+\n{'I'*len(insert_trimmed)}\n")

        if total % 5_000_000 == 0:
            print(f"  ...{total:,} reads scanned", flush=True)

    out_fastq.close()

    stats = {
        "total_reads": total,
        "chain_found": n_chain_found,
        "chain_found_pct": 100 * n_chain_found / total if total else 0,
        "insert_kept": n_insert_kept,
        "unique_raw_barcode_combos": len(combo_counts),
        "top_raw_combos": combo_counts.most_common(30),
    }
    if bc1_lookup:
        stats["unique_well_combos"] = len(well_combo_counts)
        stats["top_well_combos"] = well_combo_counts.most_common(30)

    with open(f"{args.out_prefix}.stats.json", "w") as f:
        json.dump(stats, f, indent=2)

    print(f"\nTotal reads: {total:,}")
    print(f"Full scaffold->BC1->BC2->BC3 chain found: {n_chain_found:,} ({stats['chain_found_pct']:.2f}%)")
    print(f"Reads with usable insert after polyT (>= {args.min_insert_len}bp): {n_insert_kept:,}")
    print(f"Unique raw barcode combos: {len(combo_counts):,}")
    print(f"\nTop 10 raw barcode combos:")
    for combo, cnt in combo_counts.most_common(10):
        print(f"  {combo}: {cnt}")
    if bc1_lookup:
        print(f"\nUnique real-well combos: {len(well_combo_counts):,}")
        print("Top 10 real-well combos (cells):")
        for combo, cnt in well_combo_counts.most_common(10):
            print(f"  {combo}: {cnt}")

    print(f"\nWrote {args.out_prefix}.tagged_insert.fastq.gz and {args.out_prefix}.stats.json")


if __name__ == "__main__":
    main()
