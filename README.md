# Single-Cell RNA-seq Analysis


## 1. **Loading and Preprocessing**
- **Loading Data**: Raw count data for both `hPSC organoids` and `hPSC fetal lung cells` were loaded as `AnnData` objects.
- **Removing Duplicates**: Duplicate cells were removed using `scv.pp.remove_duplicate_cells()`.
- **Filtering Genes**: Genes expressed in fewer than 30 cells were filtered out using `sc.pp.filter_genes()`.
- **Normalization and Log Transformation**: Data was normalized and log-transformed using `sc.pp.normalize_total()` and `sc.pp.log1p()`.
- **Highly Variable Genes (HVG) Selection**: Highly variable genes were identified using `sc.pp.highly_variable_genes()`.

## 2. **Dimensionality Reduction & Clustering**
- **Principal Component Analysis (PCA)**: PCA was applied using `sc.pp.pca()` to reduce the dimensionality of the data.
- **Neighbor Graph Construction**: A neighbor graph was constructed using `sc.pp.neighbors()` based on the first 30 principal components.
- **Clustering**: The Louvain algorithm (`sc.tl.louvain()`) was applied to cluster the cells based on the neighbor graph.
- **UMAP Visualization**: Uniform Manifold Approximation and Projection (UMAP) was applied using `sc.tl.umap()` to visualize the clusters.

## 3. **Differential Expression (DEG) Analysis**
- **Identifying DEGs**: Differentially expressed genes were identified between clusters using the `sc.tl.rank_genes_groups()` method with the Wilcoxon test.
- **Saving DEG Results**: DEG results for both datasets were saved for further analysis.

## 4. **Marker Gene Overlap**
- **Reference Marker Genes**: A reference list of marker genes was loaded from an external file (`quach_degs`).
- **Marker Gene Overlap**: The overlap between DEG results and reference marker genes was calculated using `sc.tl.marker_gene_overlap()`.

## 5. **Heatmap Visualization**
- **Overlap Heatmaps**: Heatmaps were generated to visualize the overlap coefficient between DEG results and reference markers for both datasets. These heatmaps were plotted using `seaborn.heatmap()`.

---

## Requirements
- `scanpy`
- `scvelo`
- `pandas`
- `numpy`
- `seaborn`
- `matplotlib`

---

