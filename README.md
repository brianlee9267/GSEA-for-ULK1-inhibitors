Overview
This script performs differential expression analysis and pathway enrichment for RNA-Seq data, derived from HEL cells treated with vehicle or ULK1 inhibitors. It uses the edgeR package for identifying differentially expressed genes (DEGs) and the clusterProfiler package for Gene Set Enrichment Analysis (GSEA). The results include a comprehensive CSV file of gene-level statistics and enriched pathways for downstream analysis.

Dependencies
The following R packages are required:

edgeR
tidyverse
clusterProfiler
org.Hs.eg.db
AnnotationDbi
enrichplot
DOSE
cowplot
ggplot2
dplyr
Install missing packages using:

For ease of use, we provide example files for both the gene counts table and sample classification CSV on GitHub. These can be downloaded and adapted to your own data by replacing the content with your individual gene counts and sample metadata.

Converts the counts matrix to a DGEList object.
Filters lowly expressed genes and normalizes library sizes.
Fits a negative binomial model to the data and performs a likelihood ratio test.
Outputs a table of DEGs, including:
Log2 fold change (logFC)
Raw p-value (PValue)
Adjusted p-value (adj_p_value)
False Discovery Rate (FDR)
Average expression in log-CPM (logCPM)

Generates a PCA plot to visualize sample clustering based on log-transformed counts.

Performs GSEA using pathways defined in the GMT file.
Outputs the top 1000 enriched pathways in a CSV file.

Usage
Run the script: Execute the script in an R environment. Ensure all file paths are correct, and required packages are installed.
