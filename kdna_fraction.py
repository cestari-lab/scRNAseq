#!/usr/bin/env python3
# Author: Lissa Cruz-Saavedra
# Date: 24-09-2026

"""
Per-cell kDNA (maxicircle) fraction - a T. brucei equivalent of "% mitochondrial"
viability QC, computed directly by genomic position since no gene is annotated
on the maxicircle contig (featureCounts can never see it).

Usage:
    python3 kdna_fraction.py <mapped_bam> <output_csv>

Requires: pysam (pip install --user pysam)
"""

import sys
import re
import pysam
import csv
from collections import defaultdict

MAXICIRCLE_CONTIG = "unitig_851_maxicircle_Tb427v9"
CB_RE = re.compile(r"_CB:(.+)_UMI:([A-Za-z]+)$")


def get_cell_barcode(read):
    """Return the cell barcode from the CB BAM tag if present, otherwise
    fall back to parsing it out of the read name (pre-tagging BAMs)."""
    if read.has_tag("CB"):
        return read.get_tag("CB")
    m = CB_RE.search(read.query_name)
    if m:
        return m.group(1)
    return None


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        sys.exit(1)

    bam_path, out_csv = sys.argv[1], sys.argv[2]

    total_reads = defaultdict(int)
    kdna_reads = defaultdict(int)

    bam = pysam.AlignmentFile(bam_path, "rb")

    contigs = set(bam.references)
    if MAXICIRCLE_CONTIG not in contigs:
        print(f"WARNING: '{MAXICIRCLE_CONTIG}' not found in BAM header contigs.")
        print("Check the exact contig name with: samtools view -H <bam> | grep -i maxi")
        sys.exit(1)

    n_reads = 0
    n_no_cb = 0
    for read in bam:
        if read.is_unmapped:
            continue
        cb = get_cell_barcode(read)
        if cb is None:
            n_no_cb += 1
            continue
        n_reads += 1
        total_reads[cb] += 1
        if read.reference_name == MAXICIRCLE_CONTIG:
            kdna_reads[cb] += 1
    bam.close()

    if n_no_cb > 0:
        print(f"NOTE: {n_no_cb} mapped reads had no recoverable cell barcode "
              f"(neither a CB tag nor a parseable read name) and were skipped.")

    cells = sorted(total_reads.keys())
    with open(out_csv, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["cell_id", "total_reads", "kdna_reads", "kdna_fraction"])
        for cell in cells:
            tot = total_reads[cell]
            kdna = kdna_reads.get(cell, 0)
            frac = kdna / tot if tot > 0 else 0.0
            writer.writerow([cell, tot, kdna, round(frac, 5)])

    print(f"Total mapped, tagged reads: {n_reads}")
    print(f"Cells: {len(cells)}")
    print(f"Cells with any kDNA reads: {sum(1 for c in cells if kdna_reads.get(c, 0) > 0)}")
    print(f"Written to: {out_csv}")


if __name__ == "__main__":
    main()
