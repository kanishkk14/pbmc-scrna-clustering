# Single-Cell PBMC Clustering and Immune-Cell Annotation

End-to-end scRNA-seq analysis pipeline using **Scanpy** on the real **PBMC3k dataset** (10x Genomics) — from raw counts to annotated immune cell populations.

---

## Background

Peripheral Blood Mononuclear Cells (PBMCs) are a standard benchmark dataset for single-cell RNA sequencing analysis. This project implements a complete analysis pipeline including quality control, normalisation, dimensionality reduction, unsupervised clustering, and biological cell-type annotation using canonical marker genes.

---

## Pipeline

```
1. Load PBMC3k  (2,700 cells × 32,738 genes)
        ↓
2. QC filtering  →  remove low-quality cells (n_genes < 200, MT% > 5%)
        ↓
3. Normalise (total-count) → log1p → select highly variable genes
        ↓
4. Scale → PCA (50 PCs) → k-NN graph → UMAP → Leiden clustering
        ↓
5. Wilcoxon DEG per cluster  (padj < 0.05, LFC > 1, pct_nz > 25%)
        ↓
6. Marker dotplot → manual biological annotation
        ↓
7. Export annotated UMAP, marker CSV, AnnData (.h5ad)
```

---

## Identified Cell Populations

| Leiden Cluster | Cell Type | Key Markers |
|---|---|---|
| 0 | CD4+ T Cells | CD3D, IL7R, CCR7 |
| 1 | NK Cells | NKG7, GNLY, PRF1 |
| 2 | B Cells | MS4A1, CD79A, CD37 |
| 3 | CD14+ Monocytes | LYZ, S100A8, S100A9 |
| 4 | Dendritic Cells | FCER1A, CST3 |
| 5 | Megakaryocytes | PPBP, PF4 |

---

## Key Outputs

| File | Description |
|---|---|
| `results/umap_pbmc_annotated.png` | UMAP coloured by cluster and cell type |
| `results/dotplot_pbmc_marker_dotplot.png` | Canonical marker expression dotplot |
| `results/heatmap_pbmc_marker_heatmap.png` | Top 5 DEGs per cluster heatmap |
| `results/violin_qc_violin.png` | QC metric distributions |
| `results/pbmc_significant_markers.csv` | Significant DEGs per cluster |

---

## Project Structure

```
pbmc_scrna/
├── main.py              — full analysis pipeline (single script)
├── requirements.txt     — Python dependencies
├── results/             — output plots and tables
└── README.md
```

---

## Getting Started

```bash
# 1. Clone
git clone https://github.com/kanishkkadyan/pbmc-scrna-clustering.git
cd pbmc-scrna-clustering

# 2. Install dependencies
pip install -r requirements.txt

# 3. Run — PBMC3k downloads automatically (~6 MB) on first run
python main.py
```

> **Note on cluster annotation:** Leiden clustering is stochastic. If cluster IDs differ from those in the README, update the `CLUSTER_ANNOTATION` dictionary at the top of `main.py` after inspecting the marker dotplot.

---

## Methods

### Quality Control
- Cells with < 200 detected genes removed
- Genes expressed in < 3 cells removed
- Cells with > 2,500 genes (potential doublets) removed
- Cells with > 5% mitochondrial counts (potential dead cells) removed

### Normalisation
- Total-count normalisation to 10,000 counts per cell
- Log1p transformation

### Dimensionality Reduction
- Highly variable gene selection (mean: 0.0125–3, dispersion > 0.5)
- Z-score scaling (max value = 10)
- PCA with ARPACK solver
- k-NN graph: 10 neighbours, 40 PCs
- UMAP embedding

### Clustering and Annotation
- Leiden algorithm (resolution = 0.5, igraph flavour)
- Wilcoxon rank-sum test for differential expression
- Marker filtering: adjusted p-value < 0.05, log-fold change > 1.0, expression fraction > 25%
- Manual annotation using canonical PBMC marker genes

---

## References

- Wolf, F. A. et al. (2018). SCANPY: large-scale single-cell gene expression data analysis. *Genome Biology*, 19, 15.
- Zheng, G. X. et al. (2017). Massively parallel digital transcriptional profiling of single cells. *Nature Communications*, 8, 14049.
- Traag, V. A. et al. (2019). From Louvain to Leiden: guaranteeing well-connected communities. *Scientific Reports*, 9, 5233.
