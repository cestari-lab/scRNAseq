#!/bin/bash
# Author: Lissa Cruz-Saavedra
# Date: 24-09-2026
# align_and_tag.sh
#
# Aligns the CB:Z:/UMI:Z:-tagged insert fastq (output of
# extract_scrna_toehold_barcodes.py) to the reference genome, preserving
# the CB/UMI tags from the fastq read-header comment as real BAM tags
# (minimap2 -y copies read-name comments straight into the SAM record,
# and since our comments are already in "TAG:TYPE:VALUE" format, they
# land as proper BAM tags with no extra parsing needed).
#
# Usage:
#   ./align_and_tag.sh <genome.fa> <tagged_insert.fastq.gz> <output_prefix>
#
# Requires: minimap2, samtools

set -euo pipefail

GENOME="$1"
FASTQ="$2"
OUT_PREFIX="$3"

if [ -z "${GENOME}" ] || [ -z "${FASTQ}" ] || [ -z "${OUT_PREFIX}" ]; then
    echo "Usage: $0 <genome.fa> <tagged_insert.fastq.gz> <output_prefix>"
    exit 1
fi

if [ ! -f "${GENOME}.mmi" ]; then
    echo "Building minimap2 index for ${GENOME} ..."
    minimap2 -d "${GENOME}.mmi" "${GENOME}"
fi

echo "Aligning ${FASTQ} ..."
# -y: copy read-name comments (our CB:Z:.../UMI:Z:... tags) into the SAM record as real tags
# -ax sr: short-read preset (adjust to -ax map-ont / -ax splice if ever using long reads)
minimap2 -ax sr -y -t 8 "${GENOME}.mmi" "${FASTQ}" | \
    samtools sort -@ 4 -m 512M -o "${OUT_PREFIX}.sorted.bam"
samtools index "${OUT_PREFIX}.sorted.bam"

echo "Verifying CB/UMI tags survived alignment (first 3 mapped reads):"
samtools view -F 4 "${OUT_PREFIX}.sorted.bam" | head -3 | cut -f1,12-

samtools flagstat "${OUT_PREFIX}.sorted.bam" > "${OUT_PREFIX}.flagstat.txt"
echo "Wrote ${OUT_PREFIX}.sorted.bam (+.bai) and ${OUT_PREFIX}.flagstat.txt"
