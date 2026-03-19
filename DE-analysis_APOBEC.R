#!/usr/bin/env Rscript 
#
#------------------------------ Data preparations ----------------------------#

# Clear data 
rm(list = ls())

# Load libraries
library("logging")
library("DESeq2")
library("ggplot2")
library("ggsignif")
library("cowplot")
library("dplyr")
library("readxl")
library("pdftools")
library("PCAtools")
library("pheatmap")
library("UpSetR")
library("tidyr")
library("clusterProfiler")
library("MSnID")
library("Rtsne")
library("MutationalPatterns")

# Set global variables and path to data directories
PROJECT_DIR <- "/Users/m.m.kleisman/Projects/HypermutatedALL_project/"
RESULTS_DIR <- paste0(PROJECT_DIR, "RESULTS/DE_analysis/APOBEC/")
RNASEQ_SUBTYPE_INFO <- "rnaSeq_diseaseSubtypes.xlsx"
WGS_SUBTYPE_INFO <- "ALL_ETV6RUNX1_Coded_Num.xlsx"
if (!dir.exists(RESULTS_DIR)){dir.create(RESULTS_DIR)}
setwd(RESULTS_DIR)

# Define name of count matrix (will be created if CREATE_COUNT_MATRIX == TRUE)
CREATE_COUNT_MATRIX <- FALSE
RNASEQ_COUNT_MAT <- "count_data_ALLsamples_PMC_duplicatesAndIsoformsMerged.Rda"

# Define extension of RNAseq count files (only required if CREATE_COUNT_MATRIX == TRUE)
COUNT_EXT <- "_RNA-Seq.gene_id.exon.counts.txt"

# Load ETV6::RUNX1 sample info
full_subtype_info <- as.data.frame(read_excel(tail(list.files(path = PROJECT_DIR,
                                                              pattern = WGS_SUBTYPE_INFO,
                                                              full.names = TRUE,
                                                              recursive = TRUE), 1)))
subtype_info <- full_subtype_info[full_subtype_info$Subtype == "ETV6::RUNX1", ]
rownames(subtype_info) <- gsub("P_", "", subtype_info$Sample)
wgs_sample_ids <- rownames(subtype_info) 
random_ids <- subtype_info$Anonym

# Parse subtype info
subtype_info <- subtype_info[,c("APOBEC_group", "Gender")]
subtype_info$APOBEC_group <- gsub("N/A", "Bootstrap <100", subtype_info$APOBEC_group)
subtype_info$APOBEC_group <- factor(subtype_info$APOBEC_group, 
                                    levels = c("Clonal", "Subclonal", "None", "Bootstrap <100"))
colnames(subtype_info)[1] <- "SBS2/13 level"

# Load RNAseq subtype info
rnaseq_subtypes <- as.data.frame(read_excel(tail(list.files(path = PROJECT_DIR,
                                                    pattern = RNASEQ_SUBTYPE_INFO,
                                                    full.names = TRUE,
                                                    recursive = TRUE), n = 1)))



#--------------------------- Prepare count matrix -----------------------------#

# Create count matrix of all samples
if (CREATE_COUNT_MATRIX){
  
  # Define paths to the txt files with count data
  all_count_data <- list.files(path = RESULTS_DIR,
                               pattern = COUNT_EXT,
                               full.names = TRUE,
                               recursive = TRUE)
  
  # Define sample IDs
  sample_ids <- gsub(COUNT_EXT, "", all_count_data)
  sample_ids <- basename(sample_ids)
  loginfo(paste("found:", length(sample_ids), "count files"))
  
  # Extract column names of the count data
  header <- readLines(all_count_data[1], n = 4)[4]
  header <- unlist(strsplit(gsub("# ", "", header), "\t"))
  
  # Generate a data frame containing all count data
  for (i in 1:length(sample_ids)){
    sample_count_data <- read.table(all_count_data[i],
                                    col.names = header)
    
    # Parse gene info to GENE SYMBOL (GENE ID)
    gene_row_names <- paste0(sample_count_data$GeneName, 
                             " (", sample_count_data$ID, ")")
    
    if (i == 1){
      # Add row names of gene IDs to the data frame 
      count_data_df <- data.frame(row.names = gene_row_names)
      
      # Create a data frame of gene lengths
      gene_lengths <- data.frame(Length = sample_count_data$Length,
                                 row.names = gene_row_names)
    } 
    
    # Create sample data frame with raw counts and FPKM values 
    sample_counts <- data.frame(sample_count_data$Counts, 
                                row.names = gene_row_names)
    colnames(sample_counts) <- sample_ids[i]
    
    # Add sample counts to the data frame
    count_data_df <- cbind(count_data_df, sample_counts[rownames(count_data_df),,drop=FALSE])
  }
  
  # Sum all duplicate entries in the count data frame
  unique_count_data_df <- count_data_df
  unique_count_data_df["gene"] <- sapply(strsplit(rownames(unique_count_data_df), "\\.| "), `[`, 1)
  unique_count_data_df <- unique_count_data_df %>% 
    group_by(gene) %>% 
    summarise_all(funs(sum))
  unique_count_data_df <- as.data.frame(unique_count_data_df)
  rownames(unique_count_data_df) <- unique_count_data_df$gene
  unique_count_data_df <- unique_count_data_df[, -1]
  
  # Get the mean of all duplicate entries in the gene lengths data frame
  unique_gene_lengths <- gene_lengths
  unique_gene_lengths["gene"] <- sapply(strsplit(rownames(unique_gene_lengths), "\\.| "), `[`, 1)
  unique_gene_lengths <- unique_gene_lengths %>% 
    group_by(gene) %>% 
    summarise_all(funs(mean))
  unique_gene_lengths <- as.data.frame(unique_gene_lengths)
  rownames(unique_gene_lengths) <- unique_gene_lengths$gene
  unique_gene_lengths <- unique_gene_lengths[, -1, drop = FALSE]
  
  # Calculate TPM values
  tpm_data_df <- convertCounts(t(t(unique_count_data_df)), 
                               unit = "TPM", 
                               geneLength = unique_gene_lengths$Length)
  
  # Combine raw counts and TPM values and write to out file
  temp_counts <- unique_count_data_df
  colnames(temp_counts) <- paste0(temp_counts, "_rawCounts")
  temp_tpm <- tpm_data_df
  colnames(temp_tpm) <- paste0(temp_tpm, "_TPM")
  count_data_df <- cbind(temp_counts, temp_tpm)
  save(count_data_df, file = RNASEQ_COUNT_MAT)
  
  # Clean up temporary files
  rm(temp_counts)
  rm(temp_tpm)
  count_data_df <- unique_count_data_df

# Load and parse RNAseq count matrix
} else {
  load(RNASEQ_COUNT_MAT)
  
  # Extract raw counts and TPM values
  tpm_data_df <- count_data_df[,grepl("TPM", colnames(count_data_df))]
  count_data_df <- count_data_df[,grepl("rawCounts", colnames(count_data_df))]
  colnames(tpm_data_df) <- as.character(gsub("_TPM", "", colnames(tpm_data_df)))
  colnames(count_data_df) <- gsub("_rawCounts", "", colnames(count_data_df))
}

# Calculate log TPM values
log_tpm_data_df <- log2(tpm_data_df + 1)

# Obtain count data of ETV6::RUNX1 samples
target_count_df <- count_data_df[,rownames(subtype_info)]
target_tpm_df <- tpm_data_df[,rownames(subtype_info)]
target_log_tpm_df <- log_tpm_data_df[,rownames(subtype_info)]

# Anonymize IDs
colnames(target_count_df) <- random_ids
colnames(target_tpm_df) <- random_ids
colnames(target_log_tpm_df) <- random_ids
rownames(subtype_info) <- random_ids



#----------------------------- Variance analysis -----------------------------#

## TPM count distribution
data_df <- target_count_df
mean_df <- data.frame(Mean = colMeans(data_df),
                      Sample = names(colMeans(data_df)),
                      Tumor = "B-ALL")

# Calculate means and count outliers
ordered_means <- mean_df[order(mean_df$Mean), 2]
n_row <- length(ordered_means)
outliers <- c(ordered_means[1:2], ordered_means[(n_row - 1):n_row])

# Generate plot
ggplot(mean_df, aes(x = Tumor, y = Mean)) +
  geom_boxplot(width = 0.5, position = "identity") +
  ggtitle("Mean counts of all genes") +
  ylab("Mean counts") +
  geom_text(aes(label = ifelse(Sample %in% outliers, Sample, "")), nudge_x = 0.25) +
  geom_point(stat = "identity", size = 2) +
  theme(legend.position = "none",
        axis.text = element_text(size = 11.5),
        axis.title = element_text(size = 12),
        plot.title = element_text(size = 14, face = "bold"))


## Hierarchical cluster dendrogram
data_df <- target_log_tpm_df
colnames(data_df) <- paste0(colnames(data_df), ", SBS2/13 level: ", 
                            subtype_info$`SBS2/13 level`)
hc <- hclust(dist(t(data_df)))
clusDendro <- as.dendrogram(hc)
new.par <- par(mar = par("mar") + c(5,10,5,10), cex = 0.85)
plot(clusDendro, horiz = TRUE, main = "Hierarchical cluster dendogram", 
     xlab = "Sample distance", ylab = "Height", 
     cex.axis = 1.2, cex.lab = 1.2, cex.main = 1.4)
default.par <- par(mar = c(0, 0, 0, 0))


## PCA plot
data_df <- target_log_tpm_df
pca_plot <- pca(data_df, metadata = subtype_info[,1, drop=F], removeVar = 0.1)
plot(biplot(pca_plot, x = "PC1", y = "PC2", 
            lab = NULL, 
            colby = "SBS2/13 level",
            legendPosition = 'right', 
            pointSize = 2, axisLabSize = 12, 
            titleLabSize = 15, subtitleLabSize = 13,
            legendLabSize = 10, legendIconSize = 3, 
            legendTitleSize = 12, title = 'PCA bi-plot',
            subtitle = 'PC1 versus PC2, all genes'))



#--------------------------- Gene expression profiling ------------------------#

# General function to generate heatmaps
generate_heatmap <- function(count_df, 
                             meta_data_df, 
                             meta_data_cols,
                             row_meta_data_df = NA,
                             breaks = seq(0, 10, 0.1),
                             show_rownames = TRUE, 
                             show_colnames = TRUE,
                             cluster_rows = TRUE,
                             border_col = "grey60",
                             gaps_row = NULL){
  n_cols <- (breaks[length(breaks)] * 10) - 3
  heatmap <- pheatmap(count_df, 
                      color = c(rep("black", 3), 
                                colorRampPalette(colors = c("#000081", "lightyellow", "lightyellow", "#b40000", 
                                                            "#9b0000", "#9a0000", "#690000", "#4f0000"))(n_cols)), 
                      breaks = breaks, scale = "none",
                      fontsize = 11, treeheight_row = 50, treeheight_col = 50, 
                      fontsize_row = 10, fontsize_col = 10,
                      show_rownames = show_rownames, show_colnames = show_colnames,
                      gaps_row = gaps_row, border_color = border_col,
                      cluster_rows = cluster_rows,
                      annotation_col = meta_data_df,
                      annotation_row = row_meta_data_df,
                      annotation_colors = meta_data_cols)
  return(heatmap)
}

# General annotation colors
group_cols <- c(Clonal = "#000067",
                Subclonal = "#6bb6ff",
                `Bootstrap <100` = "plum2",
                None = "grey85")


## Heatmap of the most variable genes in the data set
# Calculate variance and stdevs for log(TPM) values and sort
log_tpm_stdevs <- apply(target_log_tpm_df, 1, sd)
log_tpm_stdevs <- sort(log_tpm_stdevs, decreasing = TRUE)

# Visualize results for top X most variable genes
pdf("variableGenesHeatmap_logTPM_top100-200-250-500-1000-2000-2500-5000_ETV6RUNX1samples_121125.pdf", 
    width = 10, height = 8)
for (n_genes in c(100, 200, 250, 500, 1000, 2000, 2500, 5000)){
  
  # Select target genes
  target_genes <- names(head(log_tpm_stdevs, n_genes))
  
  # Subset counts for target genes and samples
  subset_log_tpm_df <- target_log_tpm_df[target_genes, ]
  
  # Generate heatmap
  var_plot <- generate_heatmap(count_df = subset_log_tpm_df, 
                               meta_data_df = subtype_info[,1, drop=F],
                               meta_data_cols = list(`SBS2/13 level` = group_cols),
                               show_rownames = FALSE)
}
dev.off()

# Write stats of top 200 most variable genes to Excel
var_genes <- names(head(log_tpm_stdevs, 200))
stdev_genes <- log_tpm_stdevs[names(log_tpm_stdevs) %in% var_genes]
mean_genes <- rowMeans(log_tpm_data_df[var_genes,])
stats_genes <- data.frame(Gene = var_genes, Mean = mean_genes, Stdev = stdev_genes)
openxlsx::write.xlsx(stats_genes[sort(rownames(stats_genes)),], file = "variableGenes.xlsx")


## Generate PCA
# Create data frame with genes to visualise (e.g. 200 most variable genes)
subset_log_tpm_df <- target_log_tpm_df[names(head(log_tpm_stdevs, 200)), ]

# Conduct PCA and visualise
pca_plot <- pca(subset_log_tpm_df, 
                metadata = subtype_info[,1, drop=F], 
                removeVar = NULL)
var_pca <- plot(biplot(pca_plot, x = "PC1", y = "PC2", 
                       lab = NULL, 
                       colby = "SBS2/13 level", 
                       colkey = group_cols,
                       legendPosition = 'right', 
                       pointSize = 2, axisLabSize = 12, 
                       labSize = 2.5,
                       titleLabSize = 15, subtitleLabSize = 13,
                       legendLabSize = 10, legendIconSize = 3, 
                       legendTitleSize = 12, title = 'PCA bi-plot'))

# Save to PDF
pdf("variableGenesPCA_top200Genes_ETV6RUNX1samples.pdf", width = 7.5, height = 7)
var_pca
dev.off()



#-------------------------------- DE-analysis ---------------------------------#

## DE-analysis preparations
#
# Rename column names for correct formula definition
names(subtype_info) <- make.names(names(subtype_info))

#  Define comparison groups
pos_control <- c("Clonal", "Subclonal")               # multiple groups possible
neg_control <- "None"

# Subset counts and meta data for comparison groups
target_samples <- rownames(subtype_info)[subtype_info[,"SBS2.13.level"] %in% 
                                         c(pos_control, neg_control)]
subset_meta_df <- subtype_info[target_samples, ]
subset_count_df <- target_count_df[,target_samples]
identical(rownames(subset_meta_df), colnames(subset_count_df))  # should be TRUE

# Rename pos_control variable in meta data if multiple positive controls are specified
subset_meta_df$SBS2.13.level <- as.character(subset_meta_df$SBS2.13.level)
if (length(pos_control) > 1){
  pos_control <- "Yes"
  subset_meta_df[!subset_meta_df[,"SBS2.13.level"] %in% neg_control, "SBS2.13.level"] <- pos_control
}

# Filter data set for interesting genes (count of 10 in at least 5 samples)
#keep <- rowSums(subset_count_df >= 10) >= 10
#subset_count_df <- subset_count_df[keep,]


## Conduct DE-analysis
#
#  Set factor levels to specify the reference level (= negative control)
subset_meta_df$SBS2.13.level <- factor(subset_meta_df$SBS2.13.level, 
                                       levels = c("None", pos_control))
subset_meta_df$Gender <- factor(subset_meta_df$Gender, levels = c("M", "F"))

# Run function 
dds <- DESeqDataSetFromMatrix(countData = subset_count_df, 
                              colData = subset_meta_df, 
                              design = ~Gender + SBS2.13.level)  
diff_expr <- DESeq(dds)

# Extract DE results
DE_result <- results(diff_expr, contrast = c("SBS2.13.level", pos_control, "None"))



#----------------------- Visualizing DE-analysis results ----------------------#

#  Print general DEG statistics 
up_genes <- subset(DE_result, log2FoldChange >= 1 & padj <= 0.05)
down_genes <- subset(DE_result, log2FoldChange <= -1 & padj <= 0.05)
DEGs <- rbind(up_genes[order(rownames(up_genes)), ], 
              down_genes[order(rownames(down_genes)), ])
cat(paste0("Summary statistics for ", pos_control, " vs ", neg_control, "\n",
           "Total DEGs: ", nrow(DEGs), "\n",
           "Upregulated DEGs (padj <= 0.05, LFC >= 1): ", nrow(up_genes), "\n",
           "Downregulated DEGs (padj <= 0.05, LFC <= -1): ", nrow(down_genes), "\n\n",
           "Top 20 most significantly upregulated DEGs:\n", 
           paste(rownames(head(up_genes[order(up_genes$padj), ], 20)), collapse = " "), "\n\n",
           "Top 20 most significantly downregulated DEGs:\n", 
           paste(rownames(head(down_genes[order(down_genes$padj), ], 20)), collapse = " ")))

# Add info to data frame
deg_df <- data.frame()
deg_df <- rbind(deg_df, data.frame(Gene = rownames(DEGs),
                                   LFC = DEGs$log2FoldChange,
                                   padj = DEGs$padj,
                                   Expression = c(rep("Upregulated", nrow(up_genes)), 
                                                  rep("Downregulated", nrow(down_genes))),
                                   Comparison = paste(pos_control, "vs None")))
openxlsx::write.xlsx(deg_df[,1:4], file = "DEGs-ETV6RUNX1samples_121125.xlsx")


## Volcano plot
volcano_df <- as.data.frame(DE_result) %>%
  mutate(gene = rownames(DE_result),
         log_pval = -log10(padj), 
         significance = case_when(padj < 0.05 & log2FoldChange >= 0.5 ~ "Upregulated",
                                  padj < 0.05 & log2FoldChange <= -0.5 ~ "Downregulated",
                                  TRUE ~ "Not significant"))
top_genes <- volcano_df %>% filter((padj < 0.001 & abs(log2FoldChange) > 2) | 
                                     abs(log2FoldChange) > 5)

volcano_plot <- ggplot(volcano_df, aes(x = log2FoldChange, y = log_pval, color = significance)) +
  geom_point(alpha = 0.8, size = 2) + 
  geom_text_repel(data = top_genes, aes(label = gene), size = 2.2, show.legend = FALSE) + 
  scale_color_manual(values = c("Upregulated" = "darkgreen", "Downregulated" = "blue", 
                                "Not significant" = "gray"),
                     breaks = c("Upregulated", "Downregulated", "Not significant")) +
  labs(title = "Volcano plot, bulk RNAseq data", x = "Log2 fold change",
       y = "-Log10 adjusted p-value", color = "Significance") +
  theme_classic() +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  theme(plot.title = element_text(face = "bold"),
        axis.text = element_text(size = 11),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 10.2))
pdf("DEGvolcanoPlot_ClonalAndSubconalVsNone_121125.pdf", 
    width = 10, height = 9)
volcano_plot
dev.off()


## Heatmap of DEGs
sorted_DEGs <- rbind(up_genes[order(up_genes$padj), ], 
                     down_genes[order(down_genes$padj), ])
  
# Subset log TPM counts for DEGs
#col_meta_df <- subtype_info[subtype_info$SBS2.13.level %in% c("Clonal", "Subclonal", "None"),"SBS2.13.level", drop = FALSE]
subset_log_tpm_data_df <- target_log_tpm_df[rownames(sorted_DEGs), rownames(subset_meta_df)]

# Generate heatmap
col_meta_df <- subtype_info[rownames(subset_meta_df), "SBS2.13.level", drop = FALSE]
row_meta_df <- data.frame(Expression = c(rep("Upregulated", nrow(up_genes)), 
                                         rep("Downregulated", nrow(down_genes))),
                          row.names = rownames(sorted_DEGs))
deg_heatmap <- generate_heatmap(count_df = subset_log_tpm_data_df, 
                                meta_data_df = col_meta_df,
                                meta_data_cols = list(`SBS2.13.level` = group_cols[names(group_cols) %in% col_meta_df$SBS2.13.level],
                                                      Expression = c(Upregulated = "darkgreen",
                                                                     Downregulated = "darkred")),
                                row_meta_data_df = row_meta_df,
                                gaps_row = nrow(up_genes),
                                show_rownames = FALSE,
                                cluster_rows = FALSE)
pdf("DEGheatmap_sortedByPadj_ClonalAndSubconalVsNone_121125.pdf", 
    width = 10.5, height = 10.5)
deg_heatmap
dev.off()


# Switch P0761 to left side right cluster
#dend2 <- as.dendrogram(deg_heatmap$tree_col)
#new.order <- unlist(dend2)[c(1:14, 27:15)]
#wts <- order(new.order)
#dend.reorder2 <- reorder(dend2, wts, mean)
#reordered_deg_heatmap <- pheatmap(subset_log_tpm_data_df, 
#                                  color = c(rep("black", 3), 
#                                            colorRampPalette(colors = c("#000081", "lightyellow", "lightyellow", "#b40000", 
#                                                                        "#9b0000", "#9a0000", "#690000", "#4f0000"))(97)),
#                                  fontsize = 11, treeheight_row = 50, treeheight_col = 50, 
#                                  fontsize_row = 10, fontsize_col = 10, border_color = NA,
#                                  breaks = seq(0, 10, 0.1), scale = "none", cluster_rows = FALSE,
#                                  cluster_cols = as.hclust(dend.reorder2),
#                                  show_rownames = TRUE, show_colnames = TRUE,
#                                  annotation_row = row_meta_df, gaps_row = nrow(up_genes),
#                                  annotation_col = col_meta_df, 
#                                  annotation_colors = list(APOBEC_signal = signal_cols[names(signal_cols) %in% col_meta_df$APOBEC_signal],
#                                                           Expression = c(Upregulated = "darkgreen",
#                                                                          Downregulated = "darkred")))
#pdf("DEGheatmap_HighMidVsNo_sortedByPadj.pdf", 
#    width = 10, height = 16)
#reordered_deg_heatmap
#dev.off()



#----------------------------- Enrichment analysis ----------------------------#

# Link entrez IDs to gene symbols
conv_tbl <- fetch_conversion_table(organism_name = "Homo sapiens", 
                                   from = "SYMBOL", to = "ENTREZID")

# Add ENTREZ IDs column to result
DE_result <- as.data.frame(DE_result) %>% 
  mutate(SYMBOL = rownames(.)) %>% 
  left_join(conv_tbl, by = "SYMBOL")

# Create a named ENTREZ ID vector of -log10(p-value) * sign(LFC) for GSEA
gene_ranks_entrez <- DE_result %>%
  filter(!is.na(ENTREZID), !is.na(log2FoldChange)) %>% 
  #mutate(ranking_metric = -log10(pvalue)*sign(log2FoldChange)) %>% 
  mutate(ranking_metric = log2FoldChange) %>% 
  group_by(ENTREZID) %>% 
  summarise(ranking_metric = mean(ranking_metric, na.rm = TRUE)) %>% 
  tibble::deframe() # convert to named vector

# Sort by rank
gene_ranks_entrez <- sort(gene_ranks_entrez, decreasing = TRUE)


## Over-representation analysis
# Subset gene_ranks_entrez for DEGs
DEGs <- subset(DE_result, abs(log2FoldChange) >= 1 & padj <= 0.05 & !is.na(ENTREZID))
DEG_ranks <- gene_ranks_entrez[names(gene_ranks_entrez) %in% DEGs$ENTREZID]
up_genes <- names(DEG_ranks)[DEG_ranks > 0]
down_genes <- names(DEG_ranks)[DEG_ranks < 0]

# Conduct over-representation analysis
overrepr_test <- enrichKEGG(gene = down_genes,    # Change to up-genes or down_genes
                            keyType = "kegg",
                            organism = "hsa",
                            pvalueCutoff = 0.05,
                            pAdjustMethod = "fdr")
dotplot(overrepr_test)

pdf("kegg_overrepresentTest_downGenes_ClonalAndSubconalVsNone_121125.pdf",
    height = 6, width = 9)
dotplot(overrepr_test)
dev.off()


## GSEA
# Conduct KEGG pathway enrichment analysis
kegg_gsea <- gseKEGG(geneList = gene_ranks_entrez,
                     organism = "hsa",
                     pvalueCutoff = 0.05,
                     pAdjustMethod = "fdr",
                     keyType = "ncbi-geneid")

# Visualise enriched KEGG pathways
kegg_plot <- dotplot(kegg_gsea, 
                     x = "GeneRatio",
                     showCategory = 100, 
                     font.size = 12, 
                     label_format = 100,
                     split = ".sign",
                     title = "GSEA, KEGG pathways") + 
  facet_grid(.sign ~ ., scales = "free", space = "free") +
  theme(strip.text.y = element_text(angle = 0, size = 10),
        plot.title = element_text(face = "bold", size = 12),
        axis.title = element_text(size = 11),
        axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11))

pdf("kegg_plot_logFC_ClonalAndSubconalVsNone_121125.pdf",
    height = 10, width = 9)
kegg_plot
dev.off()



#----------------------- Complete RNAseq cohort analysis ----------------------#

# Parse RNAseq subtype info
rnaseq_subtypes <- rnaseq_subtypes[rnaseq_subtypes$`Disease sub-class` == "B-ALL" &
                                   !grepl("_", rnaseq_subtypes$SKION_ID), c(1,5)]
rnaseq_subtypes <- rnaseq_subtypes[!rnaseq_subtypes$Subtype == "Unknown",]
rnaseq_subtypes <- rnaseq_subtypes %>% distinct(SKION_ID, .keep_all = TRUE)
rnaseq_subtypes$Subtype <- gsub(" positive", "", rnaseq_subtypes$Subtype)
rownames(rnaseq_subtypes) <- rnaseq_subtypes$SKION_ID
rnaseq_subtypes <- rnaseq_subtypes[,2,drop = FALSE]

# Subset log TPM data for B-ALL RNAseq samples
subset_log_tpm_df <- log_tpm_data_df[,colnames(log_tpm_data_df) %in% 
                                      rownames(rnaseq_subtypes)]
subset_subtypes <- rnaseq_subtypes[colnames(subset_log_tpm_df),,drop=F]


## tSNE of the most variable genes in the data set
# Calculate variance and stdevs for log(TPM) values and sort
log_tpm_stdevs <- apply(subset_log_tpm_df, 1, sd)
log_tpm_stdevs <- sort(log_tpm_stdevs, decreasing = TRUE)

# Create data frame with the top 1000 most variable genes 
subset_log_tpm_df_sd <- subset_log_tpm_df[names(head(log_tpm_stdevs, 1000)), ]

# Perform t-SNE
data <- t(subset_log_tpm_df_sd)
tsne_results <- Rtsne(data, dims = 2, perplexity = 30, verbose = TRUE, max_iter = 500)

# Convert results to a data frame for plotting
tsne_df <- as.data.frame(tsne_results$Y)
tsne_df$Sample <- rownames(data)  
tsne_df$Subtype <- subset_subtypes$Subtype[match(tsne_df$Sample, rownames(subset_subtypes))]

# Define plot colors
tsne_df$Color <- "Other"
tsne_df$Color[grep("ETV6", tsne_df$Subtype)] <- tsne_df$Subtype[grep("ETV6", tsne_df$Subtype)]

# Reorder data frame so that the ETV6::RUNX1 subtypes are at the end
tsne_df <- tsne_df %>%
  clusterProfiler::arrange((Color %in% c("ETV6::RUNX1", "ETV6::RUNX1-like")))

# Color samples by SBS2/SBS13 status
# tsne_df <- as.data.frame(readxl::read_excel("metaData_tSNEplot.xlsx"))
other_pos_samples <- full_subtype_info[!full_subtype_info$Subtype %in% c("ETV6::RUNX1", "Relapsed ETV6::RUNX1") &
                                        full_subtype_info$`SBS2/13` == "Yes", "Sample"]
tsne_df[tsne_df$Sample %in% gsub("P_", "", other_pos_samples), "Color"] <- "Other, SBS2/SBS13 positive"
tsne_df$Color <- ifelse(grepl("ETV6::RUNX1$|positive", tsne_df$Color), tsne_df$Color, "Other, SBS2/SBS13 negative\nor SBS2/SBS13 unknown")

# Reorder data frame so that the SBS2/13+ samples are at the end
tsne_df <- tsne_df %>%
  clusterProfiler::arrange((Color %in% c("ETV6::RUNX1", "Other, SBS2/SBS13 positive")))

# Plot t-SNE results
pdf("tSNE-completeBALLwithoutUnknowns-top1000logTPM_v3.pdf", width = 9, height = 6.5)
ggplot(tsne_df, aes(x = V1, y = V2, color = Color)) +
  geom_point(size = 0.6) +
  labs(x = "t-SNE 1", y = "t-SNE 2", color = "Subtype") +
  scale_color_manual(values = c(`Other, SBS2/SBS13 negative\nor SBS2/SBS13 unknown` = "grey80", 
                                `ETV6::RUNX1` = "darkred", `Other, SBS2/SBS13 positive` = "blue"),
                     breaks = c("ETV6::RUNX1", "Other, SBS2/SBS13 positive", "Other, SBS2/SBS13 negative\nor SBS2/SBS13 unknown")) +
  theme_bw() + theme(panel.grid.major = element_line(color = "white"),
                     panel.grid.minor = element_line(color = "white"),
                     legend.key.height = unit(0.9, "cm"),
                     legend.text = element_text(size = 10))
dev.off()
#openxlsx::write.xlsx(tsne_df, file = "metaData_tSNEplot_v2.xlsx")

# Generate plot with SBS2/SBS13 annotations
cluster_samples <- tsne_df[!tsne_df$Color == "Other",]
other_samples <- tsne_df[tsne_df$Color == "Other",]
cluster_samples$Signature <- full_subtype_info$`SBS2/13`[match(paste0("P_", cluster_samples$Sample), full_subtype_info$Sample)]
cluster_samples$Signature[is.na(cluster_samples$Signature)] <- "Unknown"
cluster_samples$Color <- paste0(ifelse(cluster_samples$Subtype == "ETV6::RUNX1", "ETV6::RUNX1", "ETV6::RUNX1-like"), 
                                       str_replace_all(cluster_samples$Signature, c("Unknown" = "\nUnknown APOBEC mutagenesis", 
                                                                                    "Yes" = "\nWith APOBEC mutagenesis", 
                                                                                    "No" = "\nWithout APOBEC mutagenesis")))
parsed_tsne_df <- rbind(other_samples, cluster_samples[,c(1:5)])
parsed_tsne_df <- parsed_tsne_df %>% clusterProfiler::arrange(grepl("With", Color))

# Plot t-SNE results
pdf("tSNE-completeBALLwithoutUnknowns-top1000logTPM_signatureAnnotations.pdf", width = 10.5, height = 7)
ggplot(parsed_tsne_df, aes(x = V1, y = V2)) +
  geom_point(size = 1.2, shape = 21,  stroke = 0.3, aes(color = Color, fill = Color)) +
  labs(x = "t-SNE 1", y = "t-SNE 2", fill = "Subtype and presence\nof SBS2/SBS13") +
  scale_fill_manual(values = c(Other = "grey80", `ETV6::RUNX1\nWith APOBEC mutagenesis` = "darkred", 
                               `ETV6::RUNX1\nWithout APOBEC mutagenesis` = "#e95277", `ETV6::RUNX1\nUnknown APOBEC mutagenesis` = "white", 
                               `ETV6::RUNX1-like\nWith APOBEC mutagenesis` = "#7c00a7", `ETV6::RUNX1-like\nUnknown APOBEC mutagenesis` = "white"),
                    breaks = c("ETV6::RUNX1\nWith APOBEC mutagenesis", "ETV6::RUNX1\nWithout APOBEC mutagenesis",
                               "ETV6::RUNX1\nUnknown APOBEC mutagenesis", "ETV6::RUNX1-like\nWith APOBEC mutagenesis", 
                                "ETV6::RUNX1-like\nUnknown APOBEC mutagenesis", "Other")) +
  scale_color_manual(values = c(Other = "grey80",  `ETV6::RUNX1\nWith APOBEC mutagenesis` = "darkred", 
                                `ETV6::RUNX1\nWithout APOBEC mutagenesis` = "#e95277", `ETV6::RUNX1\nUnknown APOBEC mutagenesis` = "#e95277", 
                                `ETV6::RUNX1-like\nWith APOBEC mutagenesis` = "#7c00a7", `ETV6::RUNX1-like\nUnknown APOBEC mutagenesis` = "#db74ff"),
                     guide = "none") +
  guides(fill = guide_legend(override.aes = list(shape = 21, size = 1.5, stroke = 0.5, 
                                                 color = c("darkred", "#e95277", "#e95277", "#7c00a7", "#db74ff", "grey80")))) +
  theme_bw() + theme(panel.grid.major = element_line(color = "white"),
                     panel.grid.minor = element_line(color = "white"),
                     legend.key.height = unit(2, "lines"), 
                     legend.text = element_text(size = 10),
                     legend.title = element_text(face = "bold", size = 10))
dev.off()

# Create mutational profile of 2 outlier samples of which we have WGS data
outliers_mut_mat <- read.table("outlierTSNEsubtypeSamples-mut-mat.tsv")
colnames(outliers_mut_mat) <- paste0(colnames(outliers_mut_mat), "\n(n=", colSums(outliers_mut_mat), ")")
pdf("mutProfile-outlierSubtypeSamplesTSNE.pdf", width = 8.5, height = 4.5)
plot_96_profile(outliers_mut_mat, ymax = 0.25) +
  xlab("Mutational context") + theme(axis.title = element_text(size = 10),
                                     strip.text.x = element_text(size = 10),
                                     strip.text.y = element_text(size = 10))
dev.off()

## ETV6::RUNX1 unsupervised clustering
colnames(subtype_info)[1] <- "SBS2/13 level"
subtype_info[["SBS2/13 level"]] <- as.character(subtype_info[["SBS2/13 level"]])

# Define ETV6::RUNX1 samples
subtype_samples <- rownames(rnaseq_subtypes)[rnaseq_subtypes == "ETV6::RUNX1"]
apobec_info <- data.frame(APOBEC = subtype_info$`SBS2/13 level`[match(subtype_samples, wgs_sample_ids)])
apobec_info$APOBEC[is.na(apobec_info$APOBEC)] <- "Unknown"
rownames(apobec_info) <- subtype_samples

# Subset log TPM data for B-ALL ETV6::RUNX1 RNAseq samples
subset_log_tpm_df <- log_tpm_data_df[,colnames(log_tpm_data_df) %in% 
                                       subtype_samples]
apobec_info <- apobec_info[colnames(subset_log_tpm_df),,drop=F]

# Define most variable genes
log_tpm_stdevs <- apply(subset_log_tpm_df, 1, sd)
log_tpm_stdevs <- sort(log_tpm_stdevs, decreasing = TRUE)

# Visualize results for top X most variable genes
pdf("variableGenesHeatmap_logTPM_top100-200-250-500-1000-2000-2500-5000_allETV6RUNX1samples.pdf", 
    width = 10, height = 8)
for (n_genes in c(100, 200, 250, 500, 1000, 2000, 2500, 5000)){
  target_genes <- names(head(log_tpm_stdevs, n_genes))
  
  # Subset counts for target genes and samples
  target_log_tpm_df <- subset_log_tpm_df[target_genes, ]
  
  # Generate heatmap
  var_plot <- generate_heatmap(count_df = target_log_tpm_df, 
                               meta_data_df = apobec_info,
                               meta_data_cols = list(APOBEC = c(group_cols, Unknown = "white")),
                               border_col = NA, show_rownames = FALSE, show_colnames = FALSE)
}
dev.off()



#----------------------- ETV6::RUNX1 specific analyses ------------------------#

# Parse gene list information of Jun Yang et al. 2025
target_genes <- as.data.frame(read_excel(list.files(path = PROJECT_DIR,
                                                    pattern = "JunYang2025_supplementaryTable.xlsx",
                                                    recursive = TRUE,
                                                    full.names = TRUE)))[6:257,2]

## Perform VST
target_data_df <- count_data_df[,colnames(count_data_df) %in% subtype_samples]
dds <- DESeqDataSetFromMatrix(countData = target_data_df, 
                              colData = apobec_info, 
                              design = ~ APOBEC) 
vst_data <- varianceStabilizingTransformation(dds, blind = TRUE)
vst_matrix <- assay(vst_data)
target_data_df <- vst_matrix[rownames(vst_matrix) %in% target_genes,]


## Perform UMAP
umap_results <- uwot::umap(t(target_data_df), 
                           n_neighbors = 20, 
                           min_dist = 0.01, 
                           metric = "euclidean",
                           seed = 12345)

# Convert results to a data frame for plotting
umap_df <- as.data.frame(umap_results)
colnames(umap_df) <- c("UMAP1", "UMAP2")
umap_df["APOBEC"] <- apobec_info$APOBEC
umap_df$APOBEC <- factor(umap_df$APOBEC, levels = c("Clonal", "Subclonal", "None", "Bootstrap <100", "Unknown"))

# Plot UMAP results
pdf("JY1.pdf", width = 9.5, height = 7)
umap_plot <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, fill = APOBEC)) +
  geom_point(shape = 21, size = 1.8, stroke = 0.7, aes(color = APOBEC)) +
  labs(fill = "APOBEC-associated\ndamage") +  
  scale_fill_manual(values = c(group_cols, Unknown = "white")) +
  scale_color_manual(values = group_cols, guide = "none") +
  scale_x_continuous(limits = c(-2.2, 2.2)) +
  guides(fill = guide_legend(override.aes = list(shape = 21, size = 2, stroke = 0.7, 
                                                 color = c("#000067", "#6bb6ff",  "grey80", "plum2", "grey80")))) +
  theme_bw() +
  theme(panel.grid.major = element_line(color = "white"),
        panel.grid.minor = element_line(color = "white"),
        legend.text = element_text(size = 10))
umap_plot

# Plot UMAP results for samples with WGS info only
subset_wgs_samples <- rownames(apobec_info)[apobec_info$APOBEC %in% c("Clonal", "Subclonal", "None")]
subset_group_cols <- group_cols[names(group_cols) %in% c("Clonal", "Subclonal", "None")]
umap_plot <- ggplot(umap_df[subset_wgs_samples,], aes(x = UMAP1, y = UMAP2, color = APOBEC)) +
  geom_point(size = 1.8) +
  labs(color = "APOBEC-associated\ndamage") +
  scale_color_manual(values = subset_group_cols) +
  scale_x_continuous(limits = c(-2.2, 2.2)) +
  theme_bw() + theme(panel.grid.major = element_line(color = "white"),
                     panel.grid.minor = element_line(color = "white"),
                     legend.text = element_text(size = 10))
umap_plot
dev.off()


## Clustering
umap_df <- as.data.frame(umap_results)
colnames(umap_df) <- c("UMAP1", "UMAP2")
kmeans_res <- kmeans(umap_df, centers = 2, nstart = 25)
umap_df$Cluster <- ifelse(kmeans_res$cluster == "1", "C1", "C2")

# Visualize for all RNAseq samples
pdf("JY2.pdf", width = 8.5, height = 7)
ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Cluster)) +
  geom_point(size = 1.8) +
  scale_color_manual(values = c(C1 = "#1eb6ff", C2 = "#f96a59")) +
  scale_x_continuous(limits = c(-2.2, 2.2)) +
  theme_bw() + theme(panel.grid.major = element_line(color = "white"),
                     panel.grid.minor = element_line(color = "white"),
                     legend.text = element_text(size = 10))

# Visualize for samples with WGS info only
ggplot(umap_df[subset_wgs_samples,], 
       aes(x = UMAP1, y = UMAP2, color = Cluster)) +
  geom_point(size = 1.8) +
  scale_color_manual(values = c(C1 = "#1eb6ff", C2 = "#f96a59")) +
  scale_x_continuous(limits = c(-2.2, 2.2)) +
  theme_bw() + theme(panel.grid.major = element_line(color = "white"),
                     panel.grid.minor = element_line(color = "white"),
                     legend.text = element_text(size = 10))
dev.off()


## Barplot of APOBEC signal per cluster
wgs_cluster_info <- apobec_info[subset_wgs_samples,, drop = FALSE]
wgs_cluster_info["Cluster"] <- umap_df$Cluster[match(subset_wgs_samples, rownames(umap_df))]
wgs_cluster_info$APOBEC <- factor(wgs_cluster_info$APOBEC, levels = c("Clonal", "Subclonal", "None"))

# Calculate significance
contingency_table <- wgs_cluster_info %>%
  dplyr::count(Cluster, APOBEC) %>%
  pivot_wider(names_from = APOBEC, values_from = n, values_fill = 0) %>%
  mutate(Yes = Clonal + Subclonal) %>%
  dplyr::select(Cluster, Yes, None)
contingency_matrix <- as.matrix(contingency_table[, -1])
rownames(contingency_matrix) <- contingency_table$Cluster
fisher_test <- fisher.test(contingency_matrix)
p_val <- round(fisher_test$p.value, digits = 4)

# Calculate percentages
wgs_cluster_info_summary <- wgs_cluster_info %>%
  group_by(Cluster, APOBEC) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(Cluster) %>%
  mutate(percentage = count / sum(count) * 100)
wgs_cluster_info_summary$APOBEC <- factor(wgs_cluster_info_summary$APOBEC,
                                          levels = rev(unique(wgs_cluster_info_summary$APOBEC)))

# Create stacked bar plot (relative)
freq_rel <- ggplot(wgs_cluster_info_summary, aes(x = Cluster, y = percentage, fill = APOBEC)) +
  geom_bar(stat = "identity", width = 0.6) +
  labs(x = "Cluster", y = "Percentage", fill = "SBS2/13 level", title = "SBS2/13 level") +
  theme_minimal() +
  scale_y_continuous(limits = c(0,103.4)) +
  geom_signif(comparisons = list(c("C1", "C2")), 
              annotations = paste0("P=", p_val), 
              y_position = 101.05,  
              tip_length = 0.03) +
  scale_fill_manual(values = group_cols, breaks = c("Clonal", "Subclonal", "None")) +
  theme(panel.grid.minor = element_line(color = "white"),
        legend.text = element_text(size = 10),
        plot.title = element_text(face = "bold"))

# Create stacked bar plot (absolute)
freq_abs <- ggplot(wgs_cluster_info_summary, aes(x = Cluster, y = count, fill = APOBEC)) +
  geom_bar(stat = "identity", width = 0.6) +
  labs(x = "Cluster", y = "Count", fill = "SBS2/13 level", title = "SBS2/13 level") +
  theme_minimal() +
  geom_signif(comparisons = list(c("C1", "C2")), 
              annotations = paste0("P=", p_val), 
              y_position = 30,  
              tip_length = 0.02) +
  scale_y_continuous(breaks = seq(0, 30, 5)) +
  scale_fill_manual(values = group_cols, breaks = c("Clonal", "Subclonal", "None")) +
  theme(panel.grid.minor = element_line(color = "white"),
        legend.text = element_text(size = 10),
        plot.title = element_text(face = "bold"))


## PAX5
pax5_df <- umap_df[, 3, drop = FALSE]
mut_info <- as.data.frame(read_excel(list.files(path = PROJECT_DIR,
                                                pattern = "ALL_SKION-IDs_okt2023.xlsx",
                                                recursive = TRUE,
                                                full.names = TRUE)))

# Parse PAX5 information
pax5_status <- mut_info$Pax5[match(rownames(pax5_df), mut_info$`SKION nr`)]
pax5_status[is.na(pax5_status)] <- "Unknown"
pax5_status[grep("inconclusive", pax5_status)] <- "Inconclusive"
pax5_status[grep("normal", pax5_status)] <- "Normal"
pax5_status[grep("del", pax5_status)] <- "Focal deletion"
pax5_df$PAX5 <- pax5_status

# Calculate significance
pax5_contingency_table <- pax5_df %>%
  dplyr::count(Cluster, PAX5) %>%
  pivot_wider(names_from = PAX5, values_from = n, values_fill = 0) %>%
  mutate(Other = Unknown + Inconclusive + Normal) %>%
  dplyr::select(Cluster, `Focal deletion`, Other)
pax5_contingency_matrix <- as.matrix(pax5_contingency_table[, -1])
rownames(pax5_contingency_matrix) <- pax5_contingency_table$Cluster
pax5_fisher_test <- fisher.test(pax5_contingency_matrix)
pax5_p_val <- signif(pax5_fisher_test$p.value, 3)

# Calculate percentages
pax5_df_summary <- pax5_df %>%
  group_by(Cluster, PAX5) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(Cluster) %>%
  mutate(percentage = count / sum(count) * 100)
pax5_df_summary$PAX5 <- factor(pax5_df_summary$PAX5,
                                levels = rev(c("Focal deletion", "Inconclusive", "Unknown", "Normal")))

# Create stacked bar plot (relative)
pax5_rel <- ggplot(pax5_df_summary, aes(x = Cluster, y = percentage, fill = PAX5)) +
  geom_bar(stat = "identity", width = 0.56) +
  labs(x = "Cluster", y = "Percentage", 
       fill = "Alteration", title = "PAX5 mutagenesis") +
  theme_minimal() +
  geom_signif(comparisons = list(c("C1", "C2")), 
              annotations = paste0("P=", pax5_p_val), 
              y_position = 99.5,  
              tip_length = 0.02) +
  scale_fill_manual(values = c(`Focal deletion` = "#1c9c95",
                              Inconclusive = "#bfa0dd",
                              Unknown = "orange2",
                              Normal = "grey90"),
                    breaks = c("Focal deletion", "Inconclusive", "Unknown", "Normal")) +
  theme(panel.grid.minor = element_line(color = "white"),
        legend.text = element_text(size = 10),
        plot.title = element_text(face = "bold"))

# Create stacked bar plot (absolute)
pax5_abs <- ggplot(pax5_df_summary, aes(x = Cluster, y = count, fill = PAX5)) +
  geom_bar(stat = "identity", width = 0.56) +
  labs(x = "Cluster", y = "Count", 
       fill = "Alteration", title = "PAX5 mutagenesis") +
  theme_minimal() +
  geom_signif(comparisons = list(c("C1", "C2")), 
              annotations = paste0("P=", pax5_p_val), 
              y_position = 48.5,  
              tip_length = 0.02) +
  scale_y_continuous(breaks = seq(0, 50, 10),
                     limits = c(0,51)) +
  scale_fill_manual(values = c(`Focal deletion` = "#1c9c95",
                               Inconclusive = "#bfa0dd",
                               Unknown = "orange2",
                               Normal = "grey90"),
                    breaks = c("Focal deletion", "Inconclusive", "Unknown", "Normal")) +
  theme(panel.grid.minor = element_line(color = "white"),
        legend.text = element_text(size = 10),
        plot.title = element_text(face = "bold"))



## Age
age_df <- umap_df[, 3, drop = FALSE]
age_info <- as.data.frame(read_excel(list.files(path = PROJECT_DIR,
                                                pattern = "ETV6-RUNX1-rnaSeq_diseaseSubtypes.xlsx",
                                                recursive = TRUE,
                                                full.names = TRUE)))

# Parse age information
age_df$Age <- age_info$Age_at_diagnosis[match(rownames(age_df), age_info$SKION_ID)]
age_df$Age <- as.numeric(age_df$Age)
age_df$Cluster <- factor(age_df$Cluster)
age_df["ParsedAge"] <- ifelse(age_df$Age <= 5, "fiveorless", "abovefive")

# Calculate significance
#age_wilcox_rank_sum <- wilcox.test(Age ~ Cluster, 
#                                   data = age_df, 
#                                   alternative = "less")
#age_p_val <- round(age_wilcox_rank_sum$p.value, digits = 4)
age_contingency_matrix <- table(age_df$Cluster, age_df$ParsedAge)
age_p_val <- fisher.test(age_contingency_matrix)$p.value

# Create dotplot
age_plot <- ggplot(age_df, aes(x = Cluster, y = Age, color = Cluster)) +
  geom_jitter(width = 0.2, height = 0) +
  geom_signif(comparisons = list(c("C1", "C2")), 
              annotations = paste0("P=", round(age_p_val, digits = 4)), 
              y_position = 19.82,  color = "black",
              tip_length = 0.02) +
  scale_y_continuous(breaks = seq(0, 20, 5),
                     limits = c(0,20.68)) +
  labs(x = "Cluster", y = "Age (in years)", title = "Age at diagnosis") +
  scale_color_manual(values = c(C1 = "#1eb6ff", C2 = "#f96a59")) +
  theme_minimal() +
  theme(panel.grid.minor = element_line(color = "white"),
        legend.text = element_text(size = 10),
        legend.position = "none",
        plot.title = element_text(face = "bold"))
age_plot

# Combine figures
pdf("JY3.pdf", width = 16, height = 7)
plot_grid(freq_abs, pax5_abs, age_plot, rel_widths = c(2,2,1.5), ncol = 3)
plot_grid(freq_rel, pax5_rel, age_plot, rel_widths = c(2,2,1.5), ncol = 3)
dev.off()

# Combine all output files into a single PDF
pdf_combine(Sys.glob("JY*.pdf"), output = "JunYangPaper-validations.pdf")
file.remove(Sys.glob("JY*.pdf"))


