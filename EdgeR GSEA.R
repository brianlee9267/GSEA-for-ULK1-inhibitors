# Install and load necessary packages
BiocManager::install("edgeR")
BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("AnnotationDbi")
BiocManager::install("enrichplot")
BiocManager::install("cowplot")
BiocManager::install("cli")

library(edgeR)
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(enrichplot)
library(DOSE)
library(cowplot)
library(ggplot2)
library(dplyr)

# Load your data
counts <- read.csv('~/Documents/Platanias lab/Diana ULK1 RNASeq GSEA Analysis EdgeR/htseq.all.counts.csv', row.names = 1) # Gene counts table
colData <- read.csv('~/Documents/Platanias lab/Diana ULK1 RNASeq GSEA Analysis EdgeR/sample_classification.csv') #Classification of sample names with appropriate treatments

# Filter colData to only include rows for which IDs are in the column names of counts
colData_filtered <- colData[colData$ID %in% colnames(counts), ]
rownames(colData_filtered) <- colData_filtered$ID

# Ensure compatibility by subsetting counts
counts <- counts[, rownames(colData_filtered)]

# Verification step
if (ncol(counts) == nrow(colData_filtered)) {
  print("The datasets are now compatible for EdgeR analysis.")
} else {
  stop("There is still a mismatch between the datasets. Check your sample IDs.")
}

# Convert counts to a DGEList object for EdgeR
group <- colData_filtered$Treatment # Treatment column
dge <- DGEList(counts = counts, group = group)

# Filter lowly expressed genes
keep <- filterByExpr(dge)
dge <- dge[keep, , keep.lib.sizes = FALSE]

# Normalize library sizes
dge <- calcNormFactors(dge)

# Create a design matrix
design <- model.matrix(~ group)
rownames(design) <- colnames(dge)

# Estimate dispersions
dge <- estimateDisp(dge, design)

# Fit a negative binomial generalized log-linear model
fit <- glmFit(dge, design)

# Perform likelihood ratio test for the comparison
lrt <- glmLRT(fit, coef = 2) # coef = 2 compares 'ULK1i' vs 'CTRL'

# Extract results and order by log2FC
res_EDGE <- topTags(lrt, n = Inf)$table
res_EDGE <- res_EDGE[order(-res_EDGE$logFC), ] # Order by descending log2FC

# Convert ENSEMBL IDs to gene symbols
symbols <- mapIds(org.Hs.eg.db, keys = rownames(res_EDGE), 
                  column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")

res_EDGE$symbol <- symbols

# Ensure FDR and adjusted p-value are calculated
res_EDGE$FDR <- p.adjust(res_EDGE$PValue, method = "BH") # False Discovery Rate
res_EDGE$adj_p_value <- res_EDGE$FDR # For clarity, duplicate FDR as adj_p_value

# Calculate log-CPM for all genes
log_cpm <- cpm(dge, log = TRUE) # Log-transformed CPM values
res_EDGE$logCPM <- rowMeans(log_cpm[rownames(res_EDGE), ])

# Remove rows with NA symbols
res_EDGE <- na.omit(res_EDGE)

# Save a comprehensive CSV with all gene information
final_results_all <- res_EDGE %>%
  dplyr::select(symbol, logFC, PValue, adj_p_value, FDR, logCPM) %>%
  arrange(desc(logFC)) # Sort by log2FC

output_file_all <- "~/Documents/Platanias lab/Diana ULK1 RNASeq GSEA Analysis EdgeR/EdgeR_results_full_analysis.csv"
write.csv(final_results_all, file = output_file_all, row.names = TRUE)

# Print confirmation message
cat("Comprehensive results for all genes saved to:", output_file_all, "\n")

# PCA plot for sample clustering
pca_res <- prcomp(t(log_cpm)) # Perform PCA

# Create PCA plot
pca_df <- data.frame(PC1 = pca_res$x[, 1], PC2 = pca_res$x[, 2], Treatment = group)
ggplot(pca_df, aes(x = PC1, y = PC2, color = Treatment)) +
  geom_point(size = 3) +
  labs(title = "PCA of Log-Transformed Counts", x = "PC1", y = "PC2") +
  theme_minimal()

# Prepare gene list for GSEA
gene_list <- res_EDGE$logFC
names(gene_list) <- res_EDGE$symbol
gene_list <- sort(gene_list, decreasing = TRUE)

# Remove duplicates and NA symbols
gene_list <- gene_list[!is.na(names(gene_list)) & !duplicated(names(gene_list))]

# Load the GMT file for GSEA
gmt_file <- "~/Documents/Platanias lab/Diana ULK1 RNASeq GSEA Analysis EdgeR/c2.all.v2024.1.Hs.symbols.gmt"
gmt <- read.gmt(gmt_file)

# Perform GSEA
set.seed(123)
gsea_results <- GSEA(geneList = gene_list,
                     TERM2GENE = gmt,
                     minGSSize = 10,
                     maxGSSize = 500,
                     pvalueCutoff = 0.1,
                     eps = 0)

# Extract top 1000 enriched pathways
top_1000_pathways <- head(as.data.frame(gsea_results), 1000)

# Save top pathways to CSV
output_file <- "~/Documents/Platanias lab/Diana ULK1 RNASeq GSEA Analysis EdgeR/top_1000_pathways.csv"
write.csv(top_1000_pathways, file = output_file, row.names = FALSE)

# Print confirmation for GSEA
if (length(gsea_results@result$ID) == 0) {
  print("No terms enriched under the specified pvalueCutoff.")
} else {
  # Extract specific pathways for visualization
  top_pathways <- gsea_results@result[c(9, 12, 28, 32, 34),]
  top_ids <- top_pathways$ID
  
  # Plot GSEA results for the selected pathways
  gseaplot2(gsea_results, geneSetID = top_ids, title = "Myc Transcriptional Signature", base_size = 2)
}

