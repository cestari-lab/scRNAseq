# scRNA-seq Analysis Pipeline (Preliminary, Poly-T Capture)

Analysis pipeline for split-and-pool combinatorial-barcoded scRNA-seq data
(bead-based capture, poly-dT priming, TSO/template-switching chemistry).

## Pipeline stages

```
raw fastq (R1 = barcode+UMI chain, R2 = cDNA insert)
    |
    v
[1] barcode/UMI extraction + arrangement quantification
    extract_scrna_toehold_barcodes.py   -> tagged, barcode-trimmed R1
                                            insert fastq + per-barcode stats
    quantify_scrna_toehold_arrangements.py -> full funnel: how many reads
                                            reach each stage of the
                                            scaffold->BC1->BC2->BC3->UMI->
                                            polyA chain, and how many show
                                            a confirmed real poly-A tail
                                            vs. barcode-complete-but-empty
    |
    v
[2] QC / diagnostics
    check_r2_for_real_mrna.py  -> tests whether R2 shows real mRNA/mapping
                                   evidence for reads R1 alone couldn't
                                   confirm (e.g. barcode-complete-no-polyA)
    highlight_illumina_reads.py  -> colored visual inspection of real R1
                                   reads (which connector/barcode segment
                                   is where), for spot-checking
    highlight_nanopore_reads.py  -> same, adapted for long-read/indel-
                                   tolerant data if ever applicable
    |
    v
[3] alignment + gene assignment + UMI-collapsed count matrix
    align_and_tag.sh       -> minimap2 alignment of the trimmed insert
                               fastq, preserving CB/UR tags as real BAM
                               tags (via minimap2 -y)
    (run featureCounts -R BAM yourself, see Usage below)
    build_count_matrix.py  -> UMI-deduplicated cell x gene count matrix
                               from the featureCounts-tagged BAM
    |
    v
[4] per-cell QC
    kdna_fraction.py  -> kinetoplast DNA (kDNA) read fraction per cell,
                          computed by genomic position on the maxicircle
                          contig (no gene model exists there)
    qc_flags.py        -> Scrublet doublet scores + kDNA-based viability
                          flag + (caveated) tubulin signal, merged into
                          one qc_flags.csv
    |
    v
[5] downstream analysis
    scanpy_analysis.py  -> cell/gene filtering, normalization, HVG
                          selection, PCA, Leiden clustering, UMAP,
                          marker-gene visualization, depth-confound check,
                          per-cluster marker genes (Wilcoxon)
```

## Requirements

```bash
pip install --user pandas numpy scanpy leidenalg igraph scrublet pysam matplotlib edlib
```

## Usage

```bash
# 1. Extract barcodes/UMI, trim R1 to the real insert, quantify completeness
python3 extract_scrna_toehold_barcodes.py --r1 R1.fastq.gz --out-prefix sample1
python3 quantify_scrna_toehold_arrangements.py --r1 R1.fastq.gz \
    --stats-out sample1_stats.json --examples-out sample1_examples.json

# 2. Align, preserving CB/UR tags; assign genes; build the UMI-deduplicated matrix
./align_and_tag.sh genome.fa sample1.tagged_insert.fastq.gz sample1_aligned
featureCounts -a annotation.gtf -o sample1_fc.txt -R BAM sample1_aligned.sorted.bam
python3 build_count_matrix.py sample1_aligned.sorted.bam.featureCounts.bam \
    sample1_cell_by_gene_counts.csv

# 3. Per-cell QC
python3 kdna_fraction.py sample1_aligned.sorted.bam kdna_fraction.csv
python3 qc_flags.py sample1_cell_by_gene_counts.csv qc_output/ kdna_fraction.csv

# 4. Clustering / downstream analysis
python3 scanpy_analysis.py sample1_cell_by_gene_counts.csv scanpy_output/
```

Requires (system packages, e.g. `apt install`): `minimap2`, `samtools`, `subread` (for `featureCounts`).
