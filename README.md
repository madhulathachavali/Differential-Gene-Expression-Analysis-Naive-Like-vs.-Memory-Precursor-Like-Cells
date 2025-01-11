# Single-Cell RNA-seq Analysis

1. Read Data:
Command: d0 <- DGEList(counts)
Purpose: Create a DGEList object to store raw counts and sample information, the required input format for differential expression analysis.
Command: anno <- read.delim("annotation.txt"); metadata <- read.csv("metadata_for_course.csv")
Purpose: Load annotation and metadata files, which include gene annotations and sample group information needed for analysis.

2. Normalize Data:
Command: d0 <- calcNormFactors(d0)
Purpose: Normalize the counts to account for differences in library sizes between samples, ensuring fair comparison of gene expression.
3. Filter Genes:
Command: keep <- filterByExpr(d0, mm); sum(keep)
Purpose: Remove genes with low expression across all samples that are unlikely to be informative, improving the statistical power of the analysis.
4. Visualize Data with MDS Plot:
Command: plotMDS(d, col = as.numeric(factor(metadata$simplified_cell_type)), cex = 1)
Purpose: Check for sample clustering and variability between groups, ensuring that the groups of interest (e.g., naive-like vs. memory precursor-like) separate well.
5. Model Design Matrix:
Command: mm <- model.matrix(~0 + simplified_cell_type, data = metadata)
Purpose: Create a design matrix that encodes group labels for comparison, specifying the experimental conditions.
6. Transform Data (Voom):
Command: y <- voom(d, mm, plot = T)
Purpose: Apply the voom transformation to estimate mean-variance relationships, converting count data to log2 counts per million (logCPM) with appropriate weights for linear modeling.
7. Fit Linear Model:
Command: fit <- lmFit(y, mm)
Purpose: Fit a linear model to the voom-transformed data, allowing estimation of group-specific expression levels for each gene.
8. Set Contrasts:
Command: contr <- makeContrasts(simplified_cell_typenaive_like - simplified_cell_typememory_precursor_like, levels = colnames(coef(fit)))
Purpose: Define the specific contrast for comparison (e.g., naive-like vs. memory precursor-like), specifying the groups to compare.
9. Estimate Contrasts:
Command: tmp <- contrasts.fit(fit, contr); tmp <- eBayes(tmp)
Purpose: Apply the contrast to estimate differential expression for each gene and use empirical Bayes shrinkage to stabilize variance estimates.
10. Multiple Testing Adjustment:
Command: top.table <- topTable(tmp, adjust.method = "BH", sort.by = "P", n = Inf)
Purpose: Adjust p-values for multiple testing using the Benjamini-Hochberg method, controlling the false discovery rate (FDR).
11. Merge Annotation:
Command: Add gene annotations (e.g., names and descriptions) to the results table.
Purpose: Link Ensembl gene IDs to meaningful gene names and descriptions for easier interpretation of results.
12. Create Volcano Plot:
Command: volcanoplot(fit2, ...)
Purpose: Visualize the differential expression results, showing the relationship between fold changes and statistical significance.
13. Create Heatmap:
Command: heatmap.2(logcpm[rownames(top.table),], ...)
Purpose: Visualize expression patterns of top genes across samples, showing clustering of genes and samples.
