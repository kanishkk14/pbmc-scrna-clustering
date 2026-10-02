"""
main.py — Single-Cell PBMC Clustering and Immune-Cell Annotation
          Dataset : PBMC3k (10x Genomics, ~2,700 peripheral blood mononuclear cells)
          Tool    : Scanpy

Pipeline
--------
1.  Load PBMC3k dataset
2.  Quality-control filtering  (cells and genes)
3.  QC metric calculation + violin plots
4.  Normalization, log-transformation, highly variable gene selection
5.  Scaling, PCA, neighbourhood graph, UMAP, Leiden clustering
6.  Wilcoxon differential expression + marker filtering
7.  Dotplot of canonical marker genes across clusters
8.  Biological cell-type annotation
9.  Export annotated UMAP, marker table, AnnData object
"""

import os
import scanpy as sc
import pandas as pd

# ── Directories ───────────────────────────────────────────────────────────────
SCRIPT_DIR  = os.path.dirname(os.path.abspath(__file__))
RESULTS_DIR = os.path.join(SCRIPT_DIR, "results")
DATA_DIR    = os.path.join(SCRIPT_DIR, "data")
os.makedirs(RESULTS_DIR, exist_ok=True)
os.makedirs(DATA_DIR,    exist_ok=True)

sc.set_figure_params(dpi=150, facecolor="white", frameon=False)
sc.settings.figdir    = RESULTS_DIR
sc.settings.datasetdir = DATA_DIR   # pbmc3k downloads here, not cwd

# ── Canonical marker genes per cell type ──────────────────────────────────────
MARKER_GENES = {
    "T cells":        ["CD3D", "CD3E", "IL7R", "CCR7", "CD8A"],
    "B cells":        ["MS4A1", "CD79A", "CD37"],
    "Monocytes":      ["LYZ", "S100A8", "S100A9", "FCGR3A", "MS4A7"],
    "NK cells":       ["NKG7", "GNLY", "PRF1"],
    "Dendritic cells":["FCER1A", "CST3"],
    "Megakaryocytes": ["PPBP", "PF4"],
}

# Annotation applied after inspecting the dotplot.
# Keys must match actual Leiden cluster IDs produced at resolution=0.5.
CLUSTER_ANNOTATION = {
    "0": "CD4+ T Cells",
    "1": "NK Cells",
    "2": "B Cells",
    "3": "CD14+ Monocytes",
    "4": "Dendritic Cells",
    "5": "Megakaryocytes",
}


# ── Step 1: Load dataset ───────────────────────────────────────────────────────

def load_data() -> sc.AnnData:
    print("=== Step 1: Loading PBMC3k dataset (~6 MB download on first run) ===")
    adata = sc.datasets.pbmc3k()
    print(f"  Cells: {adata.n_obs}  |  Genes: {adata.n_vars}")
    return adata


# ── Step 2: QC filtering ───────────────────────────────────────────────────────

def quality_control(adata: sc.AnnData) -> sc.AnnData:
    print("\n=== Step 2: QC filtering ===")

    # Cell and gene filters
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)

    # Mitochondrial gene annotation
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt"], percent_top=None, inplace=True
    )

    print(f"  After gene/cell filter — Cells: {adata.n_obs}  Genes: {adata.n_vars}")

    # Violin plot of QC metrics
    sc.pl.violin(
        adata,
        ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
        jitter=0.4,
        multi_panel=True,
        show=False,
        save="_qc_violin.png",
    )

    # Apply thresholds
    adata = adata[adata.obs.n_genes_by_counts < 2500].copy()
    adata = adata[adata.obs.pct_counts_mt < 5.0].copy()
    print(f"  After threshold filter  — Cells: {adata.n_obs}")

    return adata


# ── Step 3: Normalisation, log-transform, HVG selection ───────────────────────

def preprocess(adata: sc.AnnData) -> sc.AnnData:
    print("\n=== Step 3: Normalisation, log-transform, HVG selection ===")

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    sc.pp.highly_variable_genes(
        adata, min_mean=0.0125, max_mean=3, min_disp=0.5
    )
    n_hvg = adata.var.highly_variable.sum()
    print(f"  Highly variable genes selected: {n_hvg}")

    # Save full normalised matrix before subsetting
    adata.raw = adata
    adata = adata[:, adata.var.highly_variable].copy()

    return adata


# ── Step 4: Dimensionality reduction and clustering ───────────────────────────

def cluster(adata: sc.AnnData) -> sc.AnnData:
    print("\n=== Step 4: Scaling, PCA, UMAP, Leiden clustering ===")

    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, svd_solver="arpack")
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
    sc.tl.umap(adata)
    sc.tl.leiden(
        adata,
        resolution=0.5,
        flavor="igraph",
        directed=False,
        key_added="leiden_clusters",
    )

    n_clusters = adata.obs["leiden_clusters"].nunique()
    print(f"  Leiden clusters found: {n_clusters}")
    print(f"  Cluster IDs: {sorted(adata.obs['leiden_clusters'].unique().tolist())}")

    return adata


# ── Step 5: Differential expression ──────────────────────────────────────────

def differential_expression(adata: sc.AnnData) -> pd.DataFrame:
    print("\n=== Step 5: Wilcoxon differential expression ===")

    sc.tl.rank_genes_groups(
        adata,
        groupby="leiden_clusters",
        method="wilcoxon",
        use_raw=True,
        pts=True,
    )

    deg_df = sc.get.rank_genes_groups_df(adata, group=None)

    # Filter: adjusted p-value, log-fold change, and expression fraction
    mask = (
        (deg_df["pvals_adj"] < 0.05) &
        (deg_df["logfoldchanges"] > 1.0) &
        (deg_df["pct_nz_group"] > 0.25)   # >25% cells in group express the gene
    )
    filtered = deg_df[mask].copy()
    print(f"  Significant markers (padj<0.05, LFC>1, pct>25%): {len(filtered)}")

    return filtered


# ── Step 6: Dotplot + Heatmap ─────────────────────────────────────────────────

def plot_markers(adata: sc.AnnData) -> None:
    print("\n=== Step 6: Marker gene dotplot and heatmap ===")

    # Dotplot — expression level + fraction per cluster
    sc.pl.dotplot(
        adata,
        MARKER_GENES,
        groupby="leiden_clusters",
        use_raw=True,
        save="_pbmc_marker_dotplot.png",
        show=False,
    )

    # Heatmap — top 5 DEGs per cluster, cells ordered by cluster
    sc.pl.rank_genes_groups_heatmap(
        adata,
        n_genes=5,
        groupby="leiden_clusters",
        use_raw=True,
        swap_axes=True,
        show_gene_labels=True,
        save="_pbmc_marker_heatmap.png",
        show=False,
    )


# ── Step 7: Cell-type annotation ──────────────────────────────────────────────

def annotate(adata: sc.AnnData) -> sc.AnnData:
    print("\n=== Step 7: Cell-type annotation ===")

    actual_clusters = set(adata.obs["leiden_clusters"].unique())
    missing = actual_clusters - set(CLUSTER_ANNOTATION.keys())
    if missing:
        print(f"  WARNING: clusters {missing} have no annotation — labelled 'Unknown'.")
        print("  Update CLUSTER_ANNOTATION at the top of main.py after reviewing the dotplot.")

    mapped = adata.obs["leiden_clusters"].astype(str).map(CLUSTER_ANNOTATION)
    mapped = mapped.fillna("Unknown")
    adata.obs["cell_type"] = pd.Categorical(mapped)

    counts = adata.obs["cell_type"].value_counts()
    for cell_type, count in counts.items():
        print(f"  {cell_type:<25s}: {count} cells")

    return adata


# ── Step 8: Export ────────────────────────────────────────────────────────────

def export(adata: sc.AnnData, markers: pd.DataFrame) -> None:
    print("\n=== Step 8: Saving outputs ===")

    # Annotated UMAP
    sc.pl.umap(
        adata,
        color=["leiden_clusters", "cell_type"],
        save="_pbmc_annotated.png",
        show=False,
    )

    # Marker table
    marker_path = os.path.join(RESULTS_DIR, "pbmc_significant_markers.csv")
    markers.to_csv(marker_path, index=False)
    print(f"  Saved marker table  → {marker_path}")


# ── Entry point ───────────────────────────────────────────────────────────────

def main():
    adata   = load_data()
    adata   = quality_control(adata)
    adata   = preprocess(adata)
    adata   = cluster(adata)
    markers = differential_expression(adata)
    plot_markers(adata)
    adata   = annotate(adata)
    export(adata, markers)
    print("\nPipeline complete. All outputs saved to results/")


if __name__ == "__main__":
    main()
