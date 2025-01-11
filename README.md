# **Differential Gene Expression Analysis**

## **Aim**
The goal of this analysis is to identify differentially expressed genes between **Naive-Like** and **Memory Precursor-Like** cell types using RNA-seq data. 

---

## **Data Source**
- **Publication**: [PMC6336113](https://pmc.ncbi.nlm.nih.gov/articles/PMC6336113/#S2)  
- **Dataset**: RNA-seq count data from mouse samples, with metadata describing cell type classifications and gene annotations.
- **Training Reference**: This analysis follows the  [UC Davis Bioinformatics Training](https://ucdavis-bioinformatics-training.github.io/2022-April-GGI-DE-in-R/data_analysis/R_code_for_quizzes).

---

## **Analysis Workflow**

### **1. Data Preprocessing**
- **Input**: Raw RNA-seq count data and metadata.  
- Normalized the data to account for library size differences using `calcNormFactors`.  
- Filtered lowly expressed genes using `filterByExpr` to improve statistical power.  
- Generated MDS plot to visualize sample clustering and group separation.

#### **MDS Plot**:
<img width="461" alt="image" src="https://github.com/user-attachments/assets/7fcae73d-e64b-4209-add1-d1d7ca7bebec" />

**Interpretation**: Samples cluster according to their cell types, indicating distinct expression profiles for Naive-Like and Memory Precursor-Like cells.

---

### **2. Differential Expression Analysis**
- Designed a model matrix to encode experimental groups.  
- Applied `voom` transformation to stabilize variance and generate logCPM values.  
- Fitted a linear model using `lmFit`.  
- Defined contrasts to compare Naive-Like vs. Memory Precursor-Like cells.  
- Used `eBayes` to calculate log fold changes and adjusted p-values.
  
#### **voom transformation**:
<img width="449" alt="image" src="https://github.com/user-attachments/assets/596acc2e-13e6-4be7-82ed-edbe0c0061a8" />


#### **Volcano Plot**:
<img width="465" alt="image" src="https://github.com/user-attachments/assets/3348d8de-7c86-4401-bcab-2f06b76f32ac" />


**Interpretation**:
- Significant genes (highlighted) show substantial differential expression.  
- Examples include:
  - **Upregulated in Memory Precursor-Like**: `Tcf7`, `S1pr1`, `Irf1`.  
  - **Upregulated in Naive-Like**: `Gzmk`, `Gzma`, `Pik3ap1`.

---

### **3. Visualization of Results**
- Generated a heatmap of top differentially expressed genes to visualize expression patterns across samples.

#### **Heatmap**:
<img width="470" alt="image" src="https://github.com/user-attachments/assets/b4bf04b7-22b2-4e2b-bc72-c847eaa3930e" />


**Interpretation**:
- Genes cluster into distinct groups, highlighting differences in expression between the two cell types.  
- Sample clustering reflects group separation, validating the experimental design.

---

## **Key Results**
- Identified **top differentially expressed genes**:
  - **Naive-Like**: `Gzmk`, `Gzma`, `Pik3ap1`.  
  - **Memory Precursor-Like**: `Tcf7`, `S1pr1`, `Irf1`.
    
- Visualizations highlight distinct gene expression patterns between groups.

---

## **Conclusion**
This analysis demonstrates clear transcriptional differences between Naive-Like and Memory Precursor-Like cells. Genes such as `Tcf7` and `Gzmk` are strong candidates for further investigation as biomarkers or regulators of these states. 

---

### **Requirements**
- **R packages**:
  - `edgeR`  
  - `limma`  
  - `gplots`  
  - `RColorBrewer`

---

## **References**
- **Publication**: [PMC6336113](https://pmc.ncbi.nlm.nih.gov/articles/PMC6336113/#S2)  
- **Tools and Libraries**:
  - [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html)  
  - [limma](https://bioconductor.org/packages/release/bioc/html/limma.html)

---

## **Acknowledgments**
Special thanks to the authors of [PMC6336113](https://pmc.ncbi.nlm.nih.gov/articles/PMC6336113/#S2) for providing the data and resources and [UC Davis Bioinformatics Training](https://ucdavis-bioinformatics-training.github.io/2022-April-GGI-DE-in-R/data_analysis/R_code_for_quizzes).

---









