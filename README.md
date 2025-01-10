# SC-RNA-Analysis

1. Loading and Preprocessing:

Loaded the raw count data for both hPSC organoids and hPSC fetal lung cells.

Removed duplicate cells using scv.pp.remove_duplicate_cells().

Filtered genes that are expressed in fewer than 30 cells using sc.pp.filter_genes().

Normalized and log-transformed the data using sc.pp.normalize_total() and sc.pp.log1p().

Performed Highly Variable Genes (HVG) selection using sc.pp.highly_variable_genes().

2. Dimensionality Reduction & Clustering:

Applied PCA (sc.pp.pca()) and neighbors (sc.pp.neighbors()) to both datasets.

Used the Louvain method (sc.tl.louvain()) to cluster the cells.

Applied UMAP (sc.tl.umap()) for visualization.

3. Differential Expression (DEG) Analysis:

Used sc.tl.rank_genes_groups() with the Wilcoxon test to identify differentially expressed genes between clusters.

Saved the DEG results for both datasets.

4. Marker Gene Overlap:

Loaded a reference list of marker genes from an external file (quach_degs).

Used sc.tl.marker_gene_overlap() to measure the overlap between your DEG results and the reference markers.

5. Heatmap Visualization:

Plotted heatmaps showing the overlap coefficient between your DEG results and the reference markers for both datasets.
