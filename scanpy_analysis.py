#!/usr/bin/env python3
# Author: Lissa Cruz-Saavedra
# Date: 24-09-2026
"""
Scanpy single-cell analysis pipeline for the cell x gene UMI count matrix.

Usage:
    python3 scanpy_analysis.py ScRNA1_cell_by_gene_counts.csv output_dir/

Requires: scanpy, leidenalg, igraph
    pip install --user scanpy leidenalg igraph
"""

import sys
import os
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sc.settings.verbosity = 1

# ---- Config: adjust these as needed ----
MIN_GENES_PER_CELL = 200      # cell-calling threshold: drop barcodes below this
MIN_CELLS_PER_GENE = 3        # drop genes detected in fewer than this many cells
N_TOP_HVG = 2000               # number of highly variable genes to use
N_PCS = 20                     # principal components to use for neighbors/clustering
N_NEIGHBORS = 15
LEIDEN_RESOLUTION = 0.5

# Known marker genes to plot if present (edit as needed)
MARKER_GENES = {
    "VSG": "T.brucei.00g089810",
    "XRN_exonuclease": "T.brucei.00g103590",
}


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        sys.exit(1)

    csv_path, out_dir = sys.argv[1], sys.argv[2]
    os.makedirs(out_dir, exist_ok=True)

    # ---- 1. Load ----
    print(f"Loading {csv_path} ...")
    df = pd.read_csv(csv_path, index_col="cell_id")
    print(f"Loaded matrix: {df.shape[0]} cells x {df.shape[1]} genes")

    adata = sc.AnnData(df)
    adata.var_names_make_unique()
    sc.pp.calculate_qc_metrics(adata, inplace=True, percent_top=None)

    # ---- 2. Cell-calling / gene filtering ----
    print(f"\nBefore filtering: {adata.shape}")
    sc.pp.filter_cells(adata, min_genes=MIN_GENES_PER_CELL)
    sc.pp.filter_genes(adata, min_cells=MIN_CELLS_PER_GENE)
    print(f"After filtering (min_genes={MIN_GENES_PER_CELL}, min_cells_per_gene={MIN_CELLS_PER_GENE}): {adata.shape}")

    # ---- 3. Normalize + log transform ----
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # ---- 4. Highly variable genes ----
    sc.pp.highly_variable_genes(adata, n_top_genes=min(N_TOP_HVG, adata.shape[1]))
    print(f"Highly variable genes selected: {adata.var['highly_variable'].sum()}")

    adata.raw = adata
    adata_hvg = adata[:, adata.var.highly_variable].copy()

    # ---- 5. Scale + PCA ----
    sc.pp.scale(adata_hvg, max_value=10)
    n_pcs = min(N_PCS, adata_hvg.shape[0] - 1, adata_hvg.shape[1] - 1)
    sc.tl.pca(adata_hvg, n_comps=n_pcs)

    # ---- 6. Neighbors, Leiden clustering, UMAP ----
    sc.pp.neighbors(adata_hvg, n_neighbors=N_NEIGHBORS, n_pcs=n_pcs)
    sc.tl.leiden(adata_hvg, resolution=LEIDEN_RESOLUTION)
    sc.tl.umap(adata_hvg)

    print(f"\nClusters found: {adata_hvg.obs['leiden'].nunique()}")
    print(adata_hvg.obs["leiden"].value_counts())

    # transfer results back onto the full (unscaled) adata for gene-level plotting
    adata.obs["leiden"] = adata_hvg.obs["leiden"]
    adata.obsm["X_umap"] = adata_hvg.obsm["X_umap"]
    adata.obsm["X_pca"] = adata_hvg.obsm["X_pca"]

    adata.write(f"{out_dir}/adata_clustered.h5ad")

    # ---- 7. Depth-confound check ----
    # IMPORTANT: always check whether clusters are just splitting on sequencing
    # depth rather than real biology, especially at low cell numbers.
    depth_by_cluster = adata.obs.groupby("leiden")["total_counts"].median()
    genes_by_cluster = adata.obs.groupby("leiden")["n_genes_by_counts"].median()
    print("\nMedian total_counts per cluster (check for depth confound):")
    print(depth_by_cluster)
    print("\nMedian n_genes_by_counts per cluster:")
    print(genes_by_cluster)

    max_ratio = depth_by_cluster.max() / depth_by_cluster.min()
    if max_ratio > 1.5:
        print(f"\nWARNING: {max_ratio:.1f}x difference in median depth between clusters.")
        print("This split may be driven by sequencing depth, not biology.")
        print("Consider re-running with sc.pp.regress_out(adata_hvg, ['total_counts']) before PCA.")

    # ---- 8. Plots ----
    fig, ax = plt.subplots(figsize=(6, 5))
    sc.pl.umap(adata, color="leiden", ax=ax, show=False, title=f"Leiden clusters (n={adata.shape[0]} cells)")
    plt.tight_layout()
    plt.savefig(f"{out_dir}/umap_clusters.png", dpi=150)
    plt.close()

    fig, ax = plt.subplots(figsize=(6, 5))
    sc.pl.umap(adata, color="total_counts", ax=ax, show=False, title="UMI depth per cell")
    plt.tight_layout()
    plt.savefig(f"{out_dir}/umap_depth.png", dpi=150)
    plt.close()

    genes_present = [g for g in MARKER_GENES.values() if g in adata.var_names]
    if genes_present:
        fig, axes = plt.subplots(1, len(genes_present), figsize=(6 * len(genes_present), 5))
        if len(genes_present) == 1:
            axes = [axes]
        for ax, gene in zip(axes, genes_present):
            sc.pl.umap(adata, color=gene, ax=ax, show=False, title=gene)
        plt.tight_layout()
        plt.savefig(f"{out_dir}/umap_markers.png", dpi=150)
        plt.close()
        print(f"\nMarker gene UMAPs saved for: {genes_present}")
    else:
        print("\nNo configured marker genes found in the filtered matrix (check MARKER_GENES dict).")

    # ---- 9. Rank genes per cluster (top markers) ----
    sc.tl.rank_genes_groups(adata, groupby="leiden", method="wilcoxon", use_raw=True)
    markers_df = sc.get.rank_genes_groups_df(adata, group=None)
    markers_df.to_csv(f"{out_dir}/cluster_marker_genes.csv", index=False)
    print(f"\nTop markers per cluster written to: {out_dir}/cluster_marker_genes.csv")

    print(f"\nDone. All outputs in: {out_dir}/")


if __name__ == "__main__":
    main()
