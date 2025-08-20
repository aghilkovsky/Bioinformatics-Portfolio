if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("TCGAbiolinks", "SummarizedExperiment", "DESeq2", "clusterProfiler", "org.Hs.eg.db"))
install.packages("tidyverse")  # for dplyr, ggplot2, etc.

library(TCGAbiolinks)
library(SummarizedExperiment)
library(DESeq2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)

query <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "HTSeq - Counts",
  sample.type = c("Primary Tumor")
)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("TCGAbiolinks")

library(TCGAbiolinks)

query <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "HTSeq - Counts",
  sample.type = c("Primary Tumor")
)

query <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = c("Primary Tumor")
)

GDCdownload(query)
brca_data <- GDCprepare(query)

# Define groups
er_pos <- c("Luminal A", "Luminal B")
tnbc <- "Basal-like"

# Filter only relevant samples
brca_sub <- brca_data[, colData(brca_data)$paper_BRCA_Subtype_PAM50 %in% c(er_pos, tnbc)]

# Double-check group counts
table(colData(brca_sub)$paper_BRCA_Subtype_PAM50)

colnames(colData(brca_data))

table(brca_data$paper_BRCA_Subtype_PAM50)

brca_sub <- brca_data[, brca_data$paper_BRCA_Subtype_PAM50 %in% c("LumA", "Basal")]
table(brca_sub$paper_BRCA_Subtype_PAM50)

# Create condition factor
colData(brca_sub)$condition <- factor(brca_sub$paper_BRCA_Subtype_PAM50, levels = c("LumA", "Basal"))

# Build DESeq2 object
dds <- DESeqDataSet(brca_sub, design = ~ condition)

# Run DESeq2
dds <- DESeq(dds)
res <- results(dds)

# View top DE genes
head(res[order(res$pvalue), ])



# Strip version numbers from ENSEMBL IDs
res$ensembl <- gsub("\\..*", "", rownames(res))

# Now map to gene symbols
library(org.Hs.eg.db)
library(AnnotationDbi)

res_annot <- as.data.frame(res)
res_annot$symbol <- mapIds(org.Hs.eg.db,
                           keys = res_annot$ensembl,
                           column = "SYMBOL",
                           keytype = "ENSEMBL",
                           multiVals = "first")
# Remove rows with NA gene symbols (optional)
res_annot_clean <- res_annot[!is.na(res_annot$symbol), ]

# Remove duplicate gene symbols, keeping the one with lowest padj
res_annot_clean <- res_annot_clean[order(res_annot_clean$padj), ]
res_annot_clean <- res_annot_clean[!duplicated(res_annot_clean$symbol), ]

# Filter for significantly up/downregulated genes
res_sig <- res_annot_clean[res_annot_clean$padj < 0.05 & abs(res_annot_clean$log2FoldChange) > 1, ]
dim(res_sig)  # how many genes passed the filter

# Save all results
write.csv(as.data.frame(res_annot_clean), "DEG_results_Basal_vs_LumA_full.csv")

# Save only significant ones
write.csv(as.data.frame(res_sig), "DEG_results_Basal_vs_LumA_significant.csv")

library(ggplot2)

res_annot_clean$threshold <- as.factor(res_annot_clean$padj < 0.05 & abs(res_annot_clean$log2FoldChange) > 1)

ggplot(res_annot_clean, aes(x = log2FoldChange, y = -log10(padj), color = threshold)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = c("grey", "red")) +
  theme_minimal() +
  labs(title = "Volcano Plot: Basal vs LumA",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-value") +
  theme(legend.position = "none")



# Get Entrez IDs for significant genes
entrez_ids <- bitr(res_sig$symbol, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

# Run GO enrichment
ego <- enrichGO(gene = entrez_ids$ENTREZID,
                OrgDb = org.Hs.eg.db,
                keyType = "ENTREZID",
                ont = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.2)

# Visualize
barplot(ego, showCategory = 15, title = "GO Biological Processes")

write.csv(as.data.frame(ego), "GO_enrichment_Basal_vs_LumA.csv", row.names = FALSE)

dotplot(ego, showCategory = 20, title = "GO Dotplot: Basal vs LumA")



# Get your DE gene IDs
deg_ids <- rownames(subset(res, padj < 0.05))

# Remove version suffix (e.g., ".14") from ENSEMBL IDs
deg_ids_clean <- gsub("\\.\\d+$", "", deg_ids)

library(org.Hs.eg.db)

gene_annot <- AnnotationDbi::select(org.Hs.eg.db,
                                    keys = deg_ids_clean,
                                    keytype = "ENSEMBL",
                                    columns = c("ENTREZID", "SYMBOL"))

# Remove NAs and duplicates
entrez_ids <- unique(na.omit(gene_annot$ENTREZID))

kegg_enrich <- enrichKEGG(
  gene         = entrez_ids,
  organism     = 'hsa',      # 'hsa' is the KEGG code for Homo sapiens
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05
)

head(kegg_enrich)

dotplot(kegg_enrich, showCategory = 20, title = "KEGG Pathways: Basal vs LumA")


barplot(kegg_enrich, showCategory = 20, title = "KEGG Enrichment", font.size = 10)


BiocManager::install("ReactomePA")
library(ReactomePA)



reactome_enrich <- enrichPathway(
  gene         = entrez_ids,
  organism     = "human",
  pvalueCutoff = 0.1,
  qvalueCutoff = 0.2,
  readable     = TRUE
)

as.data.frame(reactome_enrich)

dotplot(reactome_enrich, showCategory = 20, title = "Reactome Pathway Enrichment: Basal vs LumA")
# Save as a large PNG
png("Reactome_Enrichment_Basal_vs_LumA.png", width = 1600, height = 1000, res = 150)
dotplot(reactome_enrich, showCategory = 20, title = "Reactome Pathway Enrichment: Basal vs LumA")
dev.off()

barplot(reactome_enrich, showCategory = 20, title = "Reactome Enrichment")
ggsave("reactome_barplot.png", width = 10, height = 6, dpi = 300)

# Order by adjusted p-value and get top 50 genes
top_genes <- rownames(res[order(res$padj), ])[1:50]

vsd <- vst(dds, blind = FALSE)  # or rlog(dds, blind = FALSE)
mat <- assay(vsd)[top_genes, ]

subtype_annot <- as.data.frame(colData(vsd)[, "paper_BRCA_Subtype_PAM50", drop = FALSE])
colnames(subtype_annot) <- "PAM50"

library(pheatmap)

pheatmap(mat,
         annotation_col = subtype_annot,
         show_rownames = FALSE,
         show_colnames = FALSE,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         fontsize = 10,
         main = "Top 50 DE Genes (Basal vs LumA)")


# Identify top genes by padj
top_genes <- head(order(res$padj), 10)
annotated_gene_names <- rownames(res)[top_genes]

library(ggplot2)



# Extract top 50 genes by adjusted p-value
topgenes <- rownames(res[order(res$padj), ])[1:50]
head(topgenes)

vsd_mat <- assay(vsd)[topgenes, ]
pca <- prcomp(t(vsd_mat))

pc_df <- data.frame(
  PC1 = pca$x[,1],
  PC2 = pca$x[,2],
  Subtype = brca_data$paper_BRCA_Subtype_PAM50[match(colnames(vsd_mat), colnames(brca_data))]
)

library(ggplot2)

ggplot(pc_df, aes(PC1, PC2, color = Subtype)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal() +
  labs(title = "PCA: Top 50 DEGs (Basal vs LumA)", color = "PAM50 Subtype")


# Extract gene IDs for "Cell cycle" pathway
kegg_cell_cycle_genes <- kegg_enrich@result$geneID[kegg_enrich@result$Description == "Cell cycle"]

# Split string of ENTREZ IDs into a vector
kegg_gene_list <- unique(unlist(strsplit(kegg_cell_cycle_genes, "/")))

# Map ENTREZ to ENSEMBL IDs
library(org.Hs.eg.db)
gene_map <- AnnotationDbi::select(org.Hs.eg.db, 
                                  keys = kegg_gene_list,
                                  keytype = "ENTREZID",
                                  columns = "ENSEMBL")

# Get unique ENSEMBL IDs for Cell Cycle genes
ensembl_kegg_genes <- unique(na.omit(gene_map$ENSEMBL))


# Check that ENSEMBL IDs match rownames of VSD matrix
common_genes <- intersect(ensembl_kegg_genes, rownames(vsd))

# Subset expression matrix
kegg_expr <- assay(vsd)[common_genes, ]

library(pheatmap)

pheatmap(kegg_expr,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         annotation_col = as.data.frame(colData(vsd)[, "paper_BRCA_Subtype_PAM50", drop = FALSE]),
         show_rownames = FALSE,
         main = "KEGG: Cell Cycle Genes")


length(common_genes)

# Remove version numbers from ENSEMBL IDs
rownames(vsd) <- gsub("\\..*", "", rownames(vsd))

# Now re-check intersection
common_genes <- intersect(ensembl_kegg_genes, rownames(vsd))

# Subset expression matrix
kegg_expr <- assay(vsd)[common_genes, ]

pheatmap(kegg_expr,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         annotation_col = as.data.frame(colData(vsd)[, "paper_BRCA_Subtype_PAM50", drop = FALSE]),
         show_rownames = FALSE,
         main = "KEGG: Cell Cycle Genes")
show_rownames = FALSE

# Order by row variance and take top 50–100
top_kegg_genes <- head(order(rowVars(kegg_expr), decreasing = TRUE), 50)
kegg_expr_subset <- kegg_expr[top_kegg_genes, ]

library(matrixStats)

# Make sure kegg_expr is a matrix
kegg_expr_mat <- as.matrix(kegg_expr)

# Get the top 50 most variable genes
top_kegg_genes <- head(order(rowVars(kegg_expr_mat), decreasing = TRUE), 50)

# Subset expression matrix
kegg_expr_subset <- kegg_expr_mat[top_kegg_genes, ]

# Plot
pheatmap(kegg_expr_subset,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         annotation_col = as.data.frame(colData(vsd)[, "paper_BRCA_Subtype_PAM50", drop = FALSE]),
         show_rownames = FALSE,
         main = "Top 50 Variable KEGG Cell Cycle Genes")

library(matrixStats)

# Ensure the matrix is numeric and has rownames
kegg_expr_mat <- as.matrix(kegg_expr)

# Check structure if needed
# str(kegg_expr_mat)

# Get top 50 most variable genes (explicitly using matrixStats)
top_kegg_genes <- head(order(matrixStats::rowVars(kegg_expr_mat), decreasing = TRUE), 50)

# Subset matrix to those genes
kegg_expr_subset <- kegg_expr_mat[top_kegg_genes, ]

# Plot heatmap
pheatmap(kegg_expr_subset,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         annotation_col = as.data.frame(colData(vsd)[, "paper_BRCA_Subtype_PAM50", drop = FALSE]),
         show_rownames = FALSE,
         main = "Top 50 Variable KEGG Cell Cycle Genes")

library(matrixStats)
kegg_expr_mat <- as.matrix(kegg_expr)
str(kegg_expr_mat)

# Make sure matrixStats is loaded
library(matrixStats)

# Get the top 50 most variable genes
top_kegg_genes <- head(order(rowVars(kegg_expr_mat), decreasing = TRUE), 50)

# Subset the expression matrix
kegg_expr_subset <- kegg_expr_mat[top_kegg_genes, ]

# Create the heatmap
pheatmap(kegg_expr_subset,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         annotation_col = as.data.frame(colData(vsd)[, "paper_BRCA_Subtype_PAM50", drop = FALSE]),
         show_rownames = FALSE,
         main = "Top 50 Variable KEGG Cell Cycle Genes")

library(matrixStats)

is.matrix(kegg_expr_mat)  # should return TRUE
top_kegg_genes <- head(order(rowVars(kegg_expr_mat, useNames = FALSE), decreasing = TRUE), 50)

kegg_expr_subset <- kegg_expr_mat[top_kegg_genes, ]

pheatmap(kegg_expr_subset,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         annotation_col = as.data.frame(colData(vsd)[, "paper_BRCA_Subtype_PAM50", drop = FALSE]),
         show_rownames = FALSE,
         main = "Top 50 Variable KEGG Cell Cycle Genes")


colnames(colData(brca_data))


# Load survival packages
install.packages("survminer")

library(survival)
library(survminer)

# Clean up survival time variable
brca_data$paper_days_to_death <- as.numeric(as.character(brca_data$paper_days_to_death))

# Check the structure to confirm
str(brca_data$paper_days_to_death)

sum(!is.na(brca_data$paper_days_to_death))

# Force numeric conversion
brca_data$survival_time <- as.numeric(as.character(ifelse(
  !is.na(brca_data$paper_days_to_death),
  brca_data$paper_days_to_death,
  brca_data$paper_days_to_last_followup
)))

# Recreate survival event
brca_data$survival_event <- ifelse(
  brca_data$paper_vital_status == "Dead", 1, 0
)

# Create survival object again
surv_obj <- Surv(time = brca_data$survival_time, event = brca_data$survival_event)

# Fit KM model
fit <- survfit(surv_obj ~ paper_BRCA_Subtype_PAM50, data = colData(brca_data))

# Plot
km_plot <- ggsurvplot(
  fit,
  data = colData(brca_data),
  pval = TRUE,
  risk.table = TRUE,
  title = "Survival by PAM50 Subtype",
  legend.title = "Subtype",
  palette = "Dark2"
)

# Show the plot
print(km_plot)

km_plot <- ggsurvplot(
  fit,
  data = colData(brca_data),
  pval = TRUE,
  risk.table = TRUE,
  title = "Survival by PAM50 Subtype",
  legend.title = "Subtype",
  palette = "Dark2",
  risk.table.height = 0.2,       # Shrink the risk table
  font.main = 12,                # Smaller title
  font.x = 10,
  font.y = 10,
  font.tickslab = 8,
  ggtheme = theme_minimal() + 
    theme(legend.position = "bottom")  # Move legend to save horizontal space
)

print(km_plot)

ggsave("KM_PAM50_survival_plot_clean.png", km_plot$plot, width = 8, height = 5, dpi = 300)



library(ggplot2)

# 1. Get variance for each gene
gene_vars <- rowVars(assay(vsd), useNames = FALSE)

# 2. Select top 500 most variable genes
top500 <- head(order(gene_vars, decreasing = TRUE), 500)

# 3. Subset the matrix
vsd_top500 <- assay(vsd)[top500, ]

# 4. Run PCA
pca <- prcomp(t(vsd_top500), scale. = TRUE)

# 5. Prepare metadata for plotting
pca_data <- as.data.frame(pca$x)
pca_data$Subtype <- colData(vsd)$paper_BRCA_Subtype_PAM50

# 6. Plot PCA
ggplot(pca_data, aes(x = PC1, y = PC2, color = Subtype)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal() +
  labs(title = "PCA of Top 500 Variable Genes", x = "PC1", y = "PC2")


library(pheatmap)

# Subset expression matrix to top 500 genes
heatmap_mat <- assay(vsd)[top500, ]

# Sample annotations
annotation_col <- as.data.frame(colData(vsd)[, "paper_BRCA_Subtype_PAM50", drop = FALSE])
colnames(annotation_col) <- "PAM50 Subtype"

# Plot heatmap
pheatmap(heatmap_mat,
         show_rownames = FALSE,
         show_colnames = FALSE,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         annotation_col = annotation_col,
         main = "Heatmap: Top 500 Variable Genes")



#Rank Genes
# Step 1: Rank genes by log2 fold change
gene_list <- res$log2FoldChange
names(gene_list) <- rownames(res)
gene_list <- sort(gene_list, decreasing = TRUE)

# Optional: remove NAs if any
gene_list <- gene_list[!is.na(gene_list)]


# Install if needed
if (!require("msigdbr")) install.packages("msigdbr")
library(msigdbr)

# Remove version numbers from ENSEMBL IDs
cleaned_ids <- sub("\\..*", "", names(gene_list))

# Convert to ENTREZ IDs
gene_df <- bitr(cleaned_ids,
                fromType = "ENSEMBL",
                toType = "ENTREZID",
                OrgDb = org.Hs.eg.db)

# Match with gene_list values
gene_list_entrez <- gene_list
names(gene_list_entrez) <- cleaned_ids

# Keep only successfully mapped ENTREZ IDs
gene_list_entrez <- gene_list_entrez[gene_df$ENSEMBL]
names(gene_list_entrez) <- gene_df$ENTREZID

gsea_hallmark <- GSEA(
  geneList = gene_list_entrez,
  TERM2GENE = hallmark_df,
  pvalueCutoff = 0.05,
  verbose = TRUE
)


# Dotplot of enriched Hallmark pathways
dotplot(gsea_hallmark, showCategory = 20, title = "GSEA: Hallmark Pathways")

# View top results
head(as.data.frame(gsea_hallmark))

# Export as CSV
write.csv(as.data.frame(gsea_hallmark), "GSEA_Hallmark_Results.csv", row.names = FALSE)

# Save plot
ggsave("GSEA_Hallmark_dotplot.png", width = 8, height = 6, dpi = 300)

library(enrichplot)

gseaplot2(gsea_hallmark, geneSetID = "HALLMARK_ESTROGEN_RESPONSE_EARLY")

gseaplot2(gsea_hallmark, geneSetID = "HALLMARK_E2F_TARGETS")



# Filter DEGs
top_deg <- res[which(res$padj < 0.05 & abs(res$log2FoldChange) > 1), ]
top_deg <- top_deg[order(top_deg$padj), ]
head(top_deg)

library(ggplot2)

# How many genes are differentially expressed?
summary(res)

# Filter significant DEGs with adjusted p-value < 0.05 and |log2FC| > 1
top_deg <- res[which(res$padj < 0.05 & abs(res$log2FoldChange) > 1), ]
top_deg <- top_deg[order(top_deg$padj), ]
dim(top_deg)  # Should return number of DEGs

# Pick top DEG
top_gene <- rownames(top_deg)[1]

# Extract expression for that gene
gene_expr <- data.frame(
  expression = assay(vsd)[top_gene, ],
  subtype = colData(vsd)$paper_BRCA_Subtype_PAM50
)

# Strip version numbers from res rownames
rownames(res) <- sub("\\..*", "", rownames(res))
rownames(top_deg) <- sub("\\..*", "", rownames(top_deg))

rownames(vsd) <- sub("\\..*", "", rownames(vsd))

top_gene <- rownames(top_deg)[1]

gene_expr <- data.frame(
  expression = assay(vsd)[top_gene, ],
  subtype = colData(vsd)$paper_BRCA_Subtype_PAM50
)

head(gene_expr)

ggplot(gene_expr, aes(x = subtype, y = expression, fill = subtype)) +
  geom_boxplot() +
  labs(title = paste("Expression of", top_gene),
       y = "VST Normalized Expression", x = "Subtype") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Create the plot and assign it to a variable
p <- ggplot(gene_expr, aes(x = subtype, y = expression, fill = subtype)) +
  geom_boxplot() +
  labs(title = paste("Expression of", top_gene),
       y = "VST Normalized Expression", x = "Subtype") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save the plot
ggsave("boxplot_top_gene_expression.png", plot = p, width = 8, height = 6, dpi = 300)


# Select top 50 DEGs by adjusted p-value
top50_genes <- rownames(top_deg)[1:50]

# Extract VST expression data for top 50 DEGs
top50_expr <- assay(vsd)[top50_genes, ]

# Get subtype info for annotation
subtype_anno <- data.frame(Subtype = colData(vsd)$paper_BRCA_Subtype_PAM50)
rownames(subtype_anno) <- colnames(top50_expr)

library(pheatmap)

# Plot
pheatmap(top50_expr,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         annotation_col = subtype_anno,
         show_rownames = FALSE,
         show_colnames = FALSE,
         main = "Top 50 DEGs Across BRCA Subtypes")


png("top50_DEG_heatmap.png", width = 1000, height = 800)
pheatmap(top50_expr,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         annotation_col = subtype_anno,
         show_rownames = FALSE,
         show_colnames = FALSE,
         main = "Top 50 DEGs Across BRCA Subtypes")


library(pheatmap)

pheatmap(top50_expr,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         annotation_col = subtype_anno,
         show_rownames = FALSE,
         show_colnames = FALSE,
         main = "Top 50 DEGs Across BRCA Subtypes")

library(biomaRt)

# Extract top 50 DEGs again
top50_genes <- rownames(top_deg)[1:50]  

# Subset the VST matrix
top50_expr_mat <- assay(vsd)[top50_genes, ]

             
# Remove version numbers from ENSEMBL IDs
gene_ids_clean <- gsub("\\..*", "", rownames(top50_expr_mat))











