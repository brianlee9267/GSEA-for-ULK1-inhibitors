Overview
This script performs differential expression analysis and pathway enrichment for RNA-Seq data. It uses the edgeR package for identifying differentially expressed genes (DEGs) and the clusterProfiler package for Gene Set Enrichment Analysis (GSEA). The results include a comprehensive CSV file of gene-level statistics and enriched pathways for downstream analysis.

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

R
Copy code
BiocManager::install(c("edgeR", "clusterProfiler", "org.Hs.eg.db", "AnnotationDbi", "enrichplot", "cowplot"))
install.packages(c("tidyverse", "DOSE", "ggplot2", "dplyr"))
Input Files
Gene Counts Table: A CSV file with genes as rows and samples as columns.

Path: ~/Documents/Platanias lab/Diana ULK1 RNASeq GSEA Analysis EdgeR/htseq.all.counts.csv
Example:
Copy code
GeneID,Sample1,Sample2,Sample3
ENSG000001,10,20,30
ENSG000002,5,15,25
Sample Classification: A CSV file specifying sample metadata.

Path: ~/Documents/Platanias lab/Diana ULK1 RNASeq GSEA Analysis EdgeR/sample_classification.csv
Required Columns:
ID: Sample IDs matching column names in the counts file.
Treatment: Treatment group (e.g., CTRL, ULK1i).
Gene Set File (GMT): A file defining gene sets for enrichment analysis.

Path: ~/Documents/Platanias lab/Diana ULK1 RNASeq GSEA Analysis EdgeR/c2.all.v2024.1.Hs.symbols.gmt
Example Files on GitHub
For ease of use, we provide example files for both the gene counts table and sample classification CSV on GitHub. These can be downloaded and adapted to your own data by replacing the content with your individual gene counts and sample metadata.

Gene Counts Example CSV
Sample Classification Example CSV
Steps in the Script
Data Preparation:

Filters sample metadata (colData) to match column names in the counts file.
Subsets the counts matrix to include only samples with metadata.
Differential Expression Analysis (EdgeR):

Converts the counts matrix to a DGEList object.
Filters lowly expressed genes and normalizes library sizes.
Fits a negative binomial model to the data and performs a likelihood ratio test.
Outputs a table of DEGs, including:
Log2 fold change (logFC)
Raw p-value (PValue)
Adjusted p-value (adj_p_value)
False Discovery Rate (FDR)
Average expression in log-CPM (logCPM)
Principal Component Analysis (PCA):

Generates a PCA plot to visualize sample clustering based on log-transformed counts.
Gene Set Enrichment Analysis (GSEA):

Prepares a ranked gene list from DEGs (sorted by logFC).
Performs GSEA using pathways defined in the GMT file.
Outputs the top 1000 enriched pathways.
Visualization:

Creates a PCA plot.
Visualizes enrichment results for selected pathways.
Outputs
Comprehensive DEG Results:
A CSV file containing all genes with their statistics.

Path: ~/Documents/Platanias lab/Diana ULK1 RNASeq GSEA Analysis EdgeR/EdgeR_results_full_analysis.csv
Columns:
symbol: Gene symbol.
logFC: Log2 fold change.
PValue: Raw p-value.
adj_p_value: Adjusted p-value (Benjamini-Hochberg).
FDR: False Discovery Rate.
logCPM: Average log-transformed counts per million.
Top Enriched Pathways:
A CSV file containing the top 1000 enriched pathways.

Path: ~/Documents/Platanias lab/Diana ULK1 RNASeq GSEA Analysis EdgeR/top_1000_pathways.csv
PCA Plot:
A scatterplot of samples in PCA space, saved as a figure or displayed in the R environment.

Usage
Run the script: Execute the script in an R environment. Ensure all file paths are correct, and required packages are installed.

Adapt to Your Data:

Replace the example gene counts table and sample classification files with your own data.
Use the example CSVs provided on GitHub as templates for adapting your data to the format required by the script.
Interpret the Outputs:

Use the comprehensive DEG results for further validation or visualization.
Explore the enriched pathways to identify biological processes or molecular functions of interest.
Key Notes
Gene Symbol Conversion: ENSEMBL IDs in the counts file are converted to gene symbols using org.Hs.eg.db. Ensure this matches your species (use org.Mm.eg.db for mouse).
Custom GMT Files: Replace the GMT file path with your own if analyzing custom gene sets.
Adjust Parameters: Modify the pvalueCutoff, minGSSize, and maxGSSize in GSEA to suit your dataset.
Contact
For questions or issues, contact Diana Saleiro at diana.saleiro@northwestern.edu
