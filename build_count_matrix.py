#!/usr/bin/env python3
"""
build_count_matrix.py

Builds the cell x gene UMI count matrix from a featureCounts-tagged BAM
(the *.featureCounts.bam produced by `featureCounts -R BAM`, run on the
output of align_and_tag.sh). Deduplicates PCR copies by (cell barcode,
gene, UMI) -- multiple reads sharing the same CB+gene+UMI are collapsed
to a single count, since they most likely originated from the same
original captured molecule, amplified by PCR.

Reads BAM tags:
    CB  cell barcode (set by extract_scrna_toehold_barcodes.py)
    UR  raw UMI (set by extract_scrna_toehold_barcodes.py)
    XT  assigned gene (set by featureCounts -R BAM)
    XS  assignment status (only "Assigned" reads are counted)

Usage:
    python3 build_count_matrix.py <tagged.featureCounts.bam> <out_matrix.csv>

Requires: pysam (pip install --user pysam)
"""
import sys
import csv
from collections import defaultdict


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        sys.exit(1)

    bam_path, out_csv = sys.argv[1], sys.argv[2]

    import pysam
    bam = pysam.AlignmentFile(bam_path, "rb")

    # (cell, gene) -> set of UMIs seen (for dedup) -> final count = len(set)
    seen_umis = defaultdict(set)
    genes_seen = set()
    cells_seen = set()

    n_total = 0
    n_assigned = 0
    n_no_cb = 0
    n_no_umi = 0

    for read in bam:
        n_total += 1
        if read.is_unmapped:
            continue
        if not read.has_tag("XS") or read.get_tag("XS") != "Assigned":
            continue
        if not read.has_tag("XT"):
            continue
        if not read.has_tag("CB"):
            n_no_cb += 1
            continue
        if not read.has_tag("UR"):
            n_no_umi += 1
            continue

        cell = read.get_tag("CB")
        umi = read.get_tag("UR")
        gene = read.get_tag("XT")
        # XT can be a comma-separated list if a read overlaps >1 gene
        # (shouldn't happen with default featureCounts settings, but guard
        # against it rather than silently mis-assigning)
        if "," in gene:
            continue

        n_assigned += 1
        seen_umis[(cell, gene)].add(umi)
        genes_seen.add(gene)
        cells_seen.add(cell)

    bam.close()

    cells = sorted(cells_seen)
    genes = sorted(genes_seen)
    print(f"Total BAM records: {n_total:,}")
    print(f"Assigned to a gene with CB+UMI tags: {n_assigned:,}")
    print(f"Skipped (no CB tag): {n_no_cb:,}")
    print(f"Skipped (no UR/UMI tag): {n_no_umi:,}")
    print(f"Cells: {len(cells):,}  Genes: {len(genes):,}")

    with open(out_csv, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["cell_id"] + genes)
        for cell in cells:
            row = [cell]
            for gene in genes:
                row.append(len(seen_umis.get((cell, gene), ())))
            writer.writerow(row)

    print(f"Wrote {out_csv} ({len(cells)} cells x {len(genes)} genes)")


if __name__ == "__main__":
    main()
