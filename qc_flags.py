#!/usr/bin/env python3
# Author: Lissa Cruz-Saavedra
# Date: 24-09-2026

"""
QC flags: tubulin signal (unreliable - see caveat below), Scrublet doublet
detection, and merge with kDNA fraction (from kdna_fraction.py) if provided.

Usage:
    python3 qc_flags.py ScRNA1_cell_by_gene_counts.csv output_dir/ [kdna_fraction.csv]

Requires: scanpy, scrublet
    pip install --user scanpy scrublet

"""

import sys
import os
import pandas as pd
import numpy as np
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Confirmed structural alpha/beta-tubulin genes (true Tubulin IPR000217 domain,
# Pfam PF00091/PF10644 - excludes tubulin-modifying enzymes like ligases)
TUBULIN_GENES = [
    "T.brucei.00g012280",
    "T.brucei.00g024870",
    "T.brucei.00g120250",
    "T.brucei.00g131990",
]

SCRUBLET_EXPECTED_DOUBLET_RATE = 0.06


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        sys.exit(1)

    csv_path, out_dir = sys.argv[1], sys.argv[2]
    kdna_csv = sys.argv[3] if len(sys.argv) > 3 else None
    os.makedirs(out_dir, exist_ok=True)

    df = pd.read_csv(csv_path, index_col="cell_id")
    print(f"Loaded matrix: {df.shape[0]} cells x {df.shape[1]} genes")

    total_umi = df.sum(axis=1)

    # ---- Tubulin fraction (exploratory only - see caveat in docstring) ----
    present_tub = [g for g in TUBULIN_GENES if g in df.columns]
    tubulin_umi = df[present_tub].sum(axis=1) if present_tub else pd.Series(0, index=df.index)
    tubulin_frac = (tubulin_umi / total_umi.replace(0, np.nan)).fillna(0)

    n_zero_tubulin = (tubulin_umi == 0).sum()
    print(f"\nTubulin genes found in matrix: {present_tub}")
    print(f"Cells with zero tubulin UMIs: {n_zero_tubulin} / {len(df)} "
          f"({100*n_zero_tubulin/len(df):.1f}%)")
    if n_zero_tubulin / len(df) > 0.5:
        print("WARNING: majority of cells have zero tubulin signal.")
        print("Likely due to multi-mapping reads being discarded at the tandem")
        print("tubulin array, not true absence of expression. Do not use this")
        print("as a parasites-per-cell measure without addressing multi-mapping.")

    # ---- Scrublet doublet detection ----
    adata = sc.AnnData(df)
    adata.var_names_make_unique()

    import scrublet as scr
    n_pcs = min(10, adata.shape[0] - 1)
    scrub = scr.Scrublet(adata.X, expected_doublet_rate=SCRUBLET_EXPECTED_DOUBLET_RATE)
    doublet_scores, predicted_doublets = scrub.scrub_doublets(
        min_counts=1, min_cells=1, n_prin_comps=n_pcs
    )
    if predicted_doublets is None:
        # threshold detection failed - fall back to manual flag at score >0.2
        predicted_doublets = doublet_scores > 0.2
        print("\nNOTE: Scrublet's automatic threshold detection failed "
              "(common at low cell counts). Using a manual cutoff of 0.2 instead - "
              "treat this as a rough guide, not a validated threshold.")

    n_doublets = int(np.sum(predicted_doublets))
    print(f"\nScrublet flagged doublets: {n_doublets} / {len(df)} "
          f"({100*n_doublets/len(df):.1f}%)")
    print("NOTE: at <1000 cells Scrublet's simulation-based correction has low "
          "power. Trust the flagged cell list; do not trust any extrapolated "
          "'overall doublet rate' Scrublet prints to console.")

    # ---- Build combined QC table ----
    qc = pd.DataFrame({
        "total_umi": total_umi,
        "n_genes": (df > 0).sum(axis=1),
        "tubulin_umi": tubulin_umi,
        "tubulin_fraction": tubulin_frac,
        "scrublet_score": doublet_scores,
        "flag_doublet": predicted_doublets,
    })

    # ---- Merge kDNA fraction if provided ----
    if kdna_csv and os.path.exists(kdna_csv):
        kdna = pd.read_csv(kdna_csv, index_col="cell_id")
        qc = qc.join(kdna[["kdna_fraction"]], how="left")
        # flag high kDNA% as potential low-viability cell (threshold: adjust after
        # inspecting the distribution - this is a starting point, not a validated cutoff)
        kdna_thresh = qc["kdna_fraction"].quantile(0.95)
        qc["flag_high_kdna"] = qc["kdna_fraction"] > kdna_thresh
        print(f"\nMerged kDNA fraction. 95th percentile threshold: {kdna_thresh:.4f}")
        print(f"Cells flagged as high-kDNA (potential low viability): "
              f"{qc['flag_high_kdna'].sum()}")
    else:
        print("\nNo kDNA fraction file provided/found - skipping viability flag. "
              "Run kdna_fraction.py first and pass its output as the 3rd argument.")
        qc["kdna_fraction"] = np.nan
        qc["flag_high_kdna"] = False

    # combined "cell of concern" flag: doublet OR high kDNA
    qc["flag_any_concern"] = qc["flag_doublet"] | qc["flag_high_kdna"]

    qc.to_csv(f"{out_dir}/qc_flags.csv")
    print(f"\nCombined QC table written to: {out_dir}/qc_flags.csv")
    print(f"Cells flagged for any concern: {qc['flag_any_concern'].sum()} / {len(qc)}")

    # ---- Plots ----
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.scatter(qc["total_umi"], qc["scrublet_score"],
               c=qc["flag_doublet"].map({True: "#C44E52", False: "#4C72B0"}), s=15)
    ax.set_xlabel("Total UMIs per cell")
    ax.set_ylabel("Scrublet doublet score")
    ax.set_title("Doublet score vs depth (red = flagged doublet)")
    plt.tight_layout()
    plt.savefig(f"{out_dir}/doublet_score_vs_depth.png", dpi=150)
    plt.close()

    if kdna_csv and os.path.exists(kdna_csv):
        fig, ax = plt.subplots(figsize=(7, 5))
        ax.hist(qc["kdna_fraction"].dropna(), bins=50, color="#55A868")
        ax.axvline(kdna_thresh, color="#C44E52", linestyle="--", label="95th percentile")
        ax.set_xlabel("kDNA fraction of total mapped reads")
        ax.set_ylabel("Number of cells")
        ax.set_title("kDNA fraction per cell")
        ax.legend()
        plt.tight_layout()
        plt.savefig(f"{out_dir}/kdna_fraction_hist.png", dpi=150)
        plt.close()

    print(f"\nDone. Outputs in: {out_dir}/")


if __name__ == "__main__":
    main()
