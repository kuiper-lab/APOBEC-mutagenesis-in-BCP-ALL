#------------------------------ Data preparations -----------------------------# 

rm(list = ls())

# Load libraries
library(Seurat) 
library(SCutils)
library(qs)
library(readr)
library(readxl)
library(dplyr)
library(tidyr)
library(patchwork)
library(viridis)
library(ggplot2)
library(ggsignif)
library(ggrepel)
library(ggh4x)
library(reshape2)
library(cowplot)
library(VennDiagram)
library(clusterProfiler)
library(org.Hs.eg.db)
library(MSnID)
library(AnnotationHub)

# List global variables
RESULTS_DIR <- "/hpc/pmc_kuiper/HypermutatedALL_project/RESULTS/scRNAseq/"
CFG <- list() 
CFG$random_seed <- 2033
SEURAT_OBJECT <- "APOBEC-final-patientCellsOnly-121224.qs"
MERGED_CLUSTERS <- "merged_lib_clusters-patientCellsOnly-121224.tsv"
CLUSTER_META_DATA <- "scRNAseq-APOBEC-clusterAssignment-121224.xlsx"
TEST_UMAP_PARAMS <- FALSE
setwd(RESULTS_DIR)

# Load the Seurat object and clusters
srat <- qread(SEURAT_OBJECT)
snp_clusters <- as.data.frame(read_table(MERGED_CLUSTERS))
rownames(snp_clusters) <- snp_clusters$barcode

# Load cluster meta data
cluster_info <- as.data.frame(read_excel(CLUSTER_META_DATA))

# Create objects to store sample information
library_assignment <- c()
sample_assignment <- c()
barcode_status <- c()
apobec_status <- c()

# Add sample information to each cell
for (n_row in 1:nrow(srat@meta.data)) {
  
  # Load cell barcode
  current_row <- srat@meta.data[n_row, ]
  barcode <- rownames(current_row)
  snp_cluster <- snp_clusters[barcode,]
  
  # Check if the cell is assigned to multiple clusters
  if (grepl("/", snp_cluster$assignment)){
    lib_cluster <- cluster <- "Mixed clusters"
    apobec_status <- append(apobec_status, NA)
    
    # Add sample ID to the cluster  
  } else {
    lib_cluster <- paste(unlist(strsplit(snp_cluster$barcode, "_"))[1], 
                         snp_cluster$assignment, 
                         sep = "_")
    cluster_meta <- cluster_info[match(lib_cluster, cluster_info$Cluster), ]
    cluster <- cluster_meta$Random_ID
    apobec_status <- append(apobec_status, cluster_meta$APOBEC)
  }
  
  # Add info to object
  library_assignment <- append(library_assignment, lib_cluster)
  sample_assignment <- append(sample_assignment, cluster)
  barcode_status <- append(barcode_status, snp_cluster$status)
}

# Add meta data to seurat object
srat@meta.data$library_assignment <- library_assignment
srat@meta.data$sample_assignment <- sample_assignment
srat@meta.data$barcode_status <- barcode_status
srat@meta.data$apobec_status <- apobec_status

# Show all clusters and cell- and droplet types in the data 
table(srat@meta.data$celltype)
table(srat@meta.data$barcode_status)
table(srat@meta.data$library_assignment)

# Load cellranger data set and define custom colors
data(refdata_cellranger_GRCh38_3.0.0)
CFG$gradient.colors <- viridis(101, direction = -1)



#-------------------------------- Data filtering ------------------------------#

# Report median and stdev transcript counts
cell_counts_raw <- colSums(srat[["RNA"]]@counts)
median(cell_counts_raw)
sd(cell_counts_raw)

# Remove non-singlets from the data set
CFG$barcode_status <- "singlet"
dim(subset(srat, barcode_status != CFG$barcode_status))   ## [1] 31157  8713
srat <- subset(srat, barcode_status == CFG$barcode_status)

# Not all droplets contain high quality cells, filter before normalizing. Judge
# by 1) nr of transcripts (too low: no cell or dead; too high: may contain a 
# multiplet) and 2) percentage of mitochondria (too high: cell may be degraded)
srat@meta.data$sample_assignment <- factor(srat@meta.data$sample_assignment,
                                           levels = unique(cluster_info$Random_ID))

# Violin plots in linear (higher bound) and log scale (lower bound limit)
v <- VlnPlot(srat, "nCount_RNA", group.by = "sample_assignment") & 
             labs(x = NULL, x = "sample_assignment")
hlines <- geom_hline(yintercept = seq(from = 0, to = 1e+05, by = 250),
                     col = "black", 
                     linetype = 1,
                     lwd = 0.1)
v + hlines | v + scale_y_continuous(trans = "log2") + hlines

# Plot percentage mitochondrial reads. Any cell with more than 25% is suspect
# If a cell’s integrity is compromised, the cytoplasmic RNA is degraded earlier 
# than the mitochondrial (protected by mitochondrial membranes)
v <- VlnPlot(srat, "pct_mito")
v

# Inspect logNormalized expression of hemoglobin genes (>1-2 is suspect)
FeaturePlot(srat, 
            pt.size = 1, 
            feature = "hemo",
            order = TRUE, 
            cols = CFG$gradient.colors)
VlnPlot(srat, features = "hemo", pt.size = 0.1)


## Set filter cut-offs
CFG$min_txpts <- 1000
CFG$max_txpts <- 25000
CFG$max_pctmito <- 25
CFG$max_hemo <- 1

# Plot final selection
f <- FeatureScatter(srat, 
                    feature1 = "nCount_RNA", 
                    feature2 = "pct_mito", 
                    pt.size = 2) + 
  geom_hline(yintercept = CFG$max_pctmito, linetype = 2) +
  geom_vline(xintercept = c(CFG$min_txpts, CFG$max_txpts), linetype = 2)
f | f + scale_x_continuous(trans = "log2")

# Subset Seurat object
dim(srat)
## [1] 31157 82231
srat <- subset(srat, nCount_RNA >= CFG$min_txpts & 
                     nCount_RNA <= CFG$max_txpts &
                     pct_mito <= CFG$max_pctmito &
                     hemo <= CFG$max_hemo &
                     celltype == "Pro-B_cell_CD34+") 
table(srat@meta.data$library_assignment)
dim(srat)
## [1] 30089 70816

# Report median and stdev transcript counts
cell_counts <- colSums(srat[["RNA"]]@counts)
median(cell_counts)
sd(cell_counts)



#--------------------------- Dimensionality reduction -------------------------#

set.seed(CFG$random_seed)

## Normalization
# Perform log-normalization, stored in 'RNA' assay
srat <- NormalizeData(srat, 
                      normalization.method = "LogNormalize")
srat <- ScaleData(srat, 
                  features = rownames(srat), 
                  verbose = FALSE)
srat <- FindVariableFeatures(srat)

# Perform normalization with SCTransform to account for the number of transcripts
# per cell in a non-linear way. Stored in 'SCT' assay. 
#
# Less appropriate than LogNorm when comparing individual genes but works better 
# to see overal data structure in dimensional reduction and clustering
#
# DefaultAssay is automatically set to SCT. To retrieve logNormalized data, use 
# DefaultAssay(srat) <- 'RNA' or specify the assay as a function argument.
#! Note: assay argument sometimes still retrieves SCT values, better to use 1st method
srat <- SCTransform(srat, 
                    vars.to.regres=NULL, verbose=FALSE,
                    vst.flavor='v2', # version that solved some problems
                    seed.use=CFG$random_seed)
DefaultAssay(srat)
#srat <- qread(file = "filtered-APOBEC-final-patientCellsOnly-txptsBetween1000and25000-max25pctmito-max1hemo-singletCD34only-normalizedAndScaled-140725.qs")


## PCA
# Input needs to be centered and usually scaled. First dimensions (~30) contain
# most of the relevant info 
srat <- RunPCA(srat, npcs = 30)

# Generate plots
DimPlot(srat, 
        reduction = "pca", 
        pt.size = 1,
        group.by = "sample_assignment")

# Deciding the number of dimensions, 'knee' id often a good first guess
ElbowPlot(srat, ndims = 100, reduction = "pca")

# Plot PCA loadings, i.e. the weights with which a gene contributes to a PC
# Plot the top 10 genes for the 100 cells with the strongest scores in these PCs
DimHeatmap(srat, dims = 1:12, cells = 100)


## UMAP
# Compress the PC dimensions to 2 dimensions for visualization purposes
#
# Specify nr of neighbors to determine how globally the UMAP will preserve overall 
# data structure. A low nr (e.g. 10) preserves only very local structure. A high 
# nr (e.g. 50) preserves global structure better. Default of 30 often works well
#
# Check UMAP output for different settings
if (TEST_UMAP_PARAMS){
  CFG$ndims <- 30

  # Play with the n.neighbors parameter to see how things change
  for (n_ngbs in c(10, 15, 20, 25, 30, 50)){
    srat_custom <- RunUMAP(srat, 
                           dims = 1:CFG$ndims, 
                           n.neighbors = n_ngbs)
    
    # Create UMAP overviews
    p_sample_custom <- DimPlot(srat_custom, 
                               pt.size = 1, 
                               group.by = "sample_assignment") + labs(title = paste0("n.ngbs=", n_ngbs))
    p_apobec_custom <- DimPlot(srat_custom, 
                               pt.size = 1, 
                               group.by = "apobec_status") + labs(title = paste0("n.ngbs=", n_ngbs))
    
    # Assign plots to variables
    assign(paste0("sample_ngbs", n_ngbs), p_sample_custom)
    assign(paste0("apobec_ngbs", n_ngbs), p_apobec_custom)
    rm(srat_custom)
  }
  l <- theme(legend.position = "none")
  print(sample_ngbs10 + l | sample_ngbs15 + l | sample_ngbs20 + l | 
          sample_ngbs25 + l | sample_ngbs30 + l | sample_ngbs50) 
  
  # Ask which ngbs the user wants to use 
  CFG$n_ngbs <- as.numeric(readline(prompt = "Enter the value you want to use for n.ngbs: "))
  
  # Play with the n.components parameter to see how things change
  for (n_comp in c(2,3,5)){
    srat_custom <- RunUMAP(srat, 
                           dims = 1:CFG$ndims, 
                           n.neighbors = CFG$n_ngbs,
                           n.components = n_comp)
    
    # Create UMAP overview
    p_sample_custom <- DimPlot(srat_custom, 
                               pt.size = 1, 
                               group.by = "sample_assignment") + 
      labs(title = paste0("n.ngbs=", CFG$n_ngbs, ", n.comp=", n_comp))
    
    # Assign plots to variables
    assign(paste0("sample_ncomp", n_comp), p_sample_custom)
    rm(srat_custom)
  }
  print(sample_ncomp2 + l | sample_ncomp3 + l | sample_ncomp5)
  
  # Ask which ncomp the user wants to use 
  CFG$n_comp <- as.numeric(readline(prompt = "Enter the value you want to use for n.comp: "))
  
  # Play with the dims parameter to see how things change
  for (n_dims in c(20, 25, 30)){
    srat_custom <- RunUMAP(srat, 
                           dims = 1:n_dims, 
                           n.neighbors = CFG$n_ngbs,
                           n.components = CFG$n_comp)
    
    # Create UMAP overview
    p_sample_custom <- DimPlot(srat_custom, 
                               pt.size = 1, 
                               group.by = "sample_assignment") + 
      labs(title = paste0("n.ngbs=", CFG$n_ngbs, ", dims=", n_dims))
    
    # Assign plots to variables
    assign(paste0("sample_dims", n_dims), p_sample_custom)
    rm(srat_custom)
  }
  print(sample_dims20 + l | sample_dims25 + l | sample_dims30)
  
  # Ask which ncomp the user wants to use 
  CFG$ndims <- as.numeric(readline(prompt = "Enter the value you want to use for ndims: "))
  
  # Remove all tmp generated plots
  rm(list = ls(pattern = "^sample_"))
  rm(list = ls(pattern = "^apobec_"))
  
} else {
  print("Optimizing UMAP parameters not requested, previously optimized values will be used")
  CFG$n_ngbs <- 30
  CFG$n_comp <- 2
  CFG$ndims <- 30
}

## Conduct UMAP
srat <- RunUMAP(srat, 
                dims = 1:CFG$ndims, 
                seed.use = CFG$random_seed,
                n.neighbors = CFG$n_ngbs,
                n.components = CFG$n_comp)

# Generate plot
p_cells <- DimPlot(srat, 
                   reduction = "umap", # default
                   pt.size = 1, 
                   group.by = "celltype")
p_sample <- DimPlot(srat, 
                    pt.size = 1, 
                    group.by = "sample_assignment") 
p_apobec <- DimPlot(srat, 
                    pt.size = 1, 
                    group.by = "apobec_status")
p_sample | p_apobec | p_cells


## Create an overview figure with SBS2/13 status and APOBEC3A expression
# Plot manually with ggplot to manipulate plotting order
#
# Fetch UMAP coordinates and APOBEC3A expression values
DefaultAssay(srat) <- "RNA"
plot_df <- cbind(Embeddings(srat, "umap"),
                 FetchData(srat, assay = "RNA", vars = c("APOBEC3A", "apobec_status")))

# Assign colors based on SBS2/13 status and APOBEC3A expression
plot_df$cols <- ifelse(plot_df$apobec_status == "Positive", 
                       "SBS2/SBS13 positive", "SBS2/SBS13 negative")
plot_df$cols[plot_df$APOBEC3A > 0] <- "APOBEC3A+ cells"

# Reorder so that APOBEC3A-expressing cells are displayed on top
plot_df <- plot_df %>% arrange(cols == "APOBEC3A+ cells")
plot_df$cols <- factor(plot_df$cols, levels = c("SBS2/SBS13 positive", 
                                                "SBS2/SBS13 negative", 
                                                "APOBEC3A+ cells"))

# Generate plot (10x8 in)
ggplot(plot_df, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = cols), size = 0.2) +
  scale_color_manual(name = "Cell status",
                     values = c(`SBS2/SBS13 positive` = "plum",
                                `SBS2/SBS13 negative` = "lightsalmon",
                                `APOBEC3A+ cells` = "darkblue")) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  theme_classic() +
  ggtitle("UMAP colored by APOBEC3A expression and SBS2/SBS13 status") +
  theme(plot.title = element_text(face = "bold", size = 13),
        axis.text = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 11))



#------------------------- Highlight specific clusters -----------------------#

## Visualise all characteristics for an individual cluster in the original UMAP space
samples <- cluster_info$Random_ID[cluster_info$Random_ID %in% unique(srat$sample_assignment)]

# Get Seurat result of sample
for (i in 1:length(samples)){
  sample_srat <- subset(srat, sample_assignment == samples[i])
  
  # Visualise cluster characteristics
  p_sample <- DimPlot(sample_srat, 
                      pt.size = 0.7, 
                      group.by = "sample_assignment")
  p_stress <- FeaturePlot(sample_srat, 
                          pt.size = 1, 
                          feature = "stress1",
                          order = TRUE, 
                          cols = CFG$gradient.colors) + labs(title = "stress")
  p_male <- FeaturePlot(sample_srat, 
                        reduction = "umap", 
                        features='log2M', 
                        order=TRUE) + ggtitle("log2Male") + 
    scale_color_gradient(low="lightgrey", high="blue", limits=c(0,3), oob = scales::squish)
  p_female <- FeaturePlot(sample_srat, 
                          reduction = "umap", 
                          features='log2F', 
                          order=TRUE) + ggtitle("log2Female") + 
    scale_color_gradient(low="lightgrey", high="blue", limits=c(0,3), oob = scales::squish)
  print((p_sample | p_stress) / (p_male | p_female))
  rm(sample_srat)
}



#-------------------------- APOBEC-specific analyses --------------------------#

DefaultAssay(srat) <- "SCT"

## Rerun UMAP and visualise all characteristics for an individual cluster
# Define APOBEC markers
apobec_markers <- c("APOBEC3A", "APOBEC3B")

# Get Seurat result of sample
apobec_df <- data.frame()
for (i in 1:length(samples)){
  sample_srat <- subset(srat, sample_assignment == samples[i])

  # Run PCA and UMAP for sample-specific cells
  sample_srat <- RunPCA(sample_srat, npcs = 30)
  sample_srat <- RunUMAP(sample_srat, 
                         dims = 1:CFG$ndims, 
                         seed.use = CFG$random_seed,
                         n.neighbors = CFG$n_ngbs)
  
  # Visualise cluster characteristics
  p_sample <- DimPlot(sample_srat, 
                      pt.size = 0.6, 
                      group.by = "sample_assignment")
  print(p_sample)
  
  # Visualise relative expression of APOBEC genes
  p_apobec_rel <- FeaturePlot(sample_srat,
                              ncol = 3,
                              pt.size = 0.7, 
                              features = apobec_markers, 
                              order = TRUE,
                              cols = CFG$gradient.colors)
  
  # Use following code to use absolute values instead of 0-100 relative scaling 
  # Obtain data for ggplot
  DefaultAssay(sample_srat) <- "RNA"
  plot_data <- cbind(Embeddings(sample_srat, "umap"),
                     FetchData(sample_srat, assay = "RNA", vars = apobec_markers))
  
  # Reshape the data into long format for facetting and order by expression
  plot_data_long <- plot_data %>%
    pivot_longer(cols = apobec_markers, names_to = "gene", values_to = "expression") %>%
    arrange(expression)
  
  # Create the plot with facetting
  p_apobec_abs <- ggplot(plot_data_long, aes(x = UMAP_1, y = UMAP_2, color = expression)) + 
    geom_point(size = 0.6) + 
    scale_colour_gradientn(colors = CFG$gradient.colors) +
    facet_grid(~gene) + 
    theme_bw(base_size = 12) +
    labs(color = "Expression")
  print(p_apobec_rel | p_apobec_abs)
  
  # Determine the percentage of cells with APOBEC3A/3B expression
  n_cells <- dim(sample_srat)[2]
  apobec3a_barc <- rownames(plot_data)[which(plot_data$APOBEC3A > 0)]
  apobec3b_barc <- rownames(plot_data)[which(plot_data$APOBEC3B > 0)]
  n_intersect <- length(intersect(apobec3a_barc, apobec3b_barc))
  apobec_df <- rbind(apobec_df, data.frame(Sample = samples[i],
                                           APOBEC3A = length(apobec3a_barc) * 100 / n_cells,
                                           APOBEC3B = length(apobec3b_barc) * 100 / n_cells,
                                           n_total = n_cells, 
                                           n_intersect = n_intersect,
                                           n_apobec3a = length(apobec3a_barc),
                                           n_apobec3b = length(apobec3b_barc),
                                           signal = cluster_info$APOBEC[cluster_info$Random_ID == samples[i]]))
}

# Melt data frame for visualization purposes
#load("R_cache/apobec_info_df.Rda")
apobec_df_melt <- melt(apobec_df[,c("APOBEC3A", "APOBEC3B", "signal")])
apobec_df_melt$signal <- factor(apobec_df_melt$signal, levels = c("Positive", "Negative"))

# Perform statistical test for each gene
sig_results <- apobec_df_melt %>%
  group_by(variable) %>%
  summarise(p_value = wilcox.test(value ~ signal)$p.value) %>%
  mutate(significance = ifelse(p_value < 0.05, "*", ""))  # Add asterisk if p < 0.05

# Create a boxplot of percentage of cells that express APOBEC3A (7x6 in)
ggplot(apobec_df_melt, aes(x = signal, y = value)) + 
  geom_boxplot(width = 0.4, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6, color = "black") +
  facet_grid(. ~ variable) +
  theme_bw() +
  geom_signif(data = apobec_df_melt[apobec_df_melt$variable == "APOBEC3A",],
                                              comparisons = list(c("Positive", "Negative")), 
              annotations = paste0("P=", as.numeric(round(sig_results[1,"p_value"], digits = 4))), 
              y_position = 9.5,  color = "grey30",
              tip_length = 0.01, textsize = 3.7) +
  scale_y_continuous(limits = c(-0.005, 10.2)) +
  labs(x = "SBS2/SBS13 status", y = "Percentage of gene-expressing cells",    
       title = "APOBEC gene expression") +
  theme(axis.text = element_text(size = 11),
        panel.border = element_rect(colour = "grey80"),
        strip.background = element_rect(colour = "grey80"),
        strip.text = element_text(size = 11))


## VENN diagram (8x7.5 in)
venn.plot <- draw.pairwise.venn(area1 = sum(apobec_df$n_apobec3a),       
                                area2 = sum(apobec_df$n_apobec3b),     
                                cat.dist = c(0.05, 0.04),
                                cross.area = sum(apobec_df$n_intersect),
                                ext.text = F,
                                category = c("APOBEC3A", "APOBEC3B"),
                                fill = c("dodgerblue4", "salmon"),
                                alpha = c(0.4, 0.4), 
                                cex = 1.5,  
                                offset = 0.1,
                                cat.cex = 1.5,
                                margin = 0.05,
                                cat.col = c("dodgerblue4", "salmon"))
grid.newpage()
grid.draw(venn.plot)


## APOBEC expression to cell phase (all cells)
a3a_expr <- ifelse(srat[["RNA"]]@data["APOBEC3A", ] > 0, 1, 0)
a3b_expr <- ifelse(srat[["RNA"]]@data["APOBEC3B", ] > 0, 1, 0)
phase_expr_df <- data.frame(APOBEC3A = a3a_expr,
                            APOBEC3B = a3b_expr,  
                            Phase = srat$Phase)
phase_expr_df$Phase <- factor(phase_expr_df$Phase, levels = c("G1", "S", "G2M"))

# Calculate percentage of cells with APOBEC expression in each phase
a3a_df <- phase_expr_df %>% group_by(Phase) %>%
  summarise(pct_expr = mean(APOBEC3A) * 100) %>% as.data.frame()
a3b_df <- phase_expr_df %>% group_by(Phase) %>%
  summarise(pct_expr = mean(APOBEC3B) * 100) %>% as.data.frame()

# Generate plots (9x5.5 in)
plots <- list()
for (gene in c("APOBEC3A", "APOBEC3B")){
  ifelse(gene == "APOBEC3A", target_df <- a3a_df, target_df <- a3b_df)
  
  # Create contingency table for p-value pairwise comparisons
  table_df <- table(phase_expr_df$Phase, phase_expr_df[[gene]])
  colnames(table_df) <- c("Negative", "Positive")
  pval_sg1 <- fisher.test(table_df[c("S", "G1"),])
  pval_sg2m <- fisher.test(table_df[c("S", "G2M"),])
  pval_g1g2m <- fisher.test(table_df[c("G1", "G2M"),])
  
  # Generate plot
  plots <- append(plots, list(ggplot(target_df, aes(x = Phase, y = pct_expr)) +
                                geom_col(fill = ifelse(gene == "APOBEC3A", "steelblue4", "steelblue"), width = 0.6) +
                                labs(y = paste0("Percentage of ", gene, "-positive cells"), 
                                     x = "Cell cycle phase", title = gene) + ylim(0, 3.5) +
                                geom_signif(comparisons = list(c("S", "G1"), c("S", "G2M"), c("G1", "G2M")),
                                            annotations = c(ifelse(pval_sg1$p.value < 0.001, "***", "N.S."), 
                                                            ifelse(pval_sg2m$p.value < 0.001, "***", "N.S."), 
                                                            ifelse(pval_g1g2m$p.value < 0.001, "***", "N.S.")),
                                            y_position = c(3, 3.15, 3.3), size = 0.2, textsize = 3.2,
                                            tip_length = ifelse(gene == "APOBEC3A", 0.03, 0.02)) +
                                theme_minimal() + theme(plot.title = element_text(face = "bold", size = 12),
                                                        axis.title = element_text(size = 10),
                                                        panel.grid.minor = element_blank())))
}
do.call(plot_grid, plots)


## APOBEC expression to cell phase (A3A+ cells)
# Subset srat object for APOBEC3A positive cells
a3a_positive_cells <- names(which(srat[["RNA"]]@data["APOBEC3A", ] > 0))
srat_subset <- subset(srat, cells = a3a_positive_cells)

a3a_df <- (prop.table(table(phase_expr_df[phase_expr_df$APOBEC3A > 0,"Phase", drop = F])) * 100) %>% as.data.frame()
a3b_df <- (prop.table(table(phase_expr_df[phase_expr_df$APOBEC3B > 0,"Phase", drop = F])) * 100) %>% as.data.frame()

plots <- list()
for (gene in c("APOBEC3A", "APOBEC3B")){
  ifelse(gene == "APOBEC3A", target_df <- a3a_df, target_df <- a3b_df)
  
  # Generate plot
  plots <- append(plots, list(ggplot(target_df, aes(x = Phase, y = Freq)) +
                                geom_col(fill = ifelse(gene == "APOBEC3A", "steelblue4", "steelblue"), width = 0.6) +
                                labs(y = paste0("Percentage of ", gene, "-positive cells"), 
                                     x = "Cell cycle phase", title = gene) + 
                                theme_minimal() + theme(plot.title = element_text(face = "bold", size = 12),
                                                        axis.title = element_text(size = 10),
                                                        panel.grid.minor = element_blank())))
}
do.call(plot_grid, plots)


ggplot(a3a_df, aes(x = "", fill = Phase)) +
  geom_bar(position = "fill", color = NA) +
  scale_y_continuous(labels = scales::percent) 


#-------------------------------- Marker genes --------------------------------#

DefaultAssay(srat) <- "RNA"

# Marker genes can be used to make sense of the data
# Selection from literature and/or (cluster-specific) DE-analysis
genes <- c("APOBEC3A", "APOBEC3B")
FeaturePlot(srat, 
            pt.size = 0.7, 
            features = genes, 
            cols = CFG$gradient.colors,
            order = TRUE) + labs(color = "Relative\nexpression")

# Absolute values are more informative but cannot be visualized with FeaturePlot
plot_data <- FetchData(srat, assay = "RNA",
                       vars = c("UMAP_1", "UMAP_2", genes))

# Reshape the data into long format for facetting and order by expression
plot_data_long <- plot_data %>%
  pivot_longer(cols = genes, names_to = "gene", values_to = "expression") %>%
  arrange(expression)

# Create the plot with facetting
p_gene_absolute <- ggplot(plot_data_long, aes(x = UMAP_1, y = UMAP_2, color = expression)) + 
  geom_point(size = 0.6) + 
  scale_colour_gradientn(colors = CFG$gradient.colors) +
  facet_grid(~gene) + 
  theme_bw(base_size = 12) +
  labs(color= "Absolute\nexpression")
p_apobec | p_gene_absolute


## Genes that are upregulated in bulk RNAseq data
bulk_markers <- c("ABCA4", "ACY3", "AHRR", "AP001057.1", "CKM", "COL4A1", "KLC3",  "NDRG4", 
                  "NGFR", "OPLAH", "PDCD6IP", "POTEE", "RASD2", "SMC1B", "AC004034.1", 
                  "BMS1", "IFI44", "INPP4B", "LINC01262", "NXN",  "RGMB", "SLC35D2")
FeaturePlot(srat, features = bulk_markers, cols = CFG$gradient.colors)

# Determine nr of positive cells per bulk marker
marker_df <- data.frame()
for (marker in bulk_markers){
  marker_pos_bc <- names(which(srat[["RNA"]]@data[marker, ] > 0))
  
  for (i in 1:length(samples)){
    
    # Determine the percentage of sample cells with the marker gene expression
    sample_bc <- names(which(srat$sample_assignment == samples[i]))
    n_cells <- length(sample_bc)
    n_pos_cells <- sum(sample_bc %in% marker_pos_bc)
    
    # Store in data frame
    marker_df <- rbind(marker_df, data.frame(Sample = samples[i],
                                             Gene = marker,
                                             Percentage = n_pos_cells * 100 / n_cells,
                                             Signal = cluster_info$APOBEC[cluster_info$Random_ID == samples[i]],
                                             Expression = ""))
  }
}

# Generate barplot (16.5x6 in)
marker_df$Signal <- factor(marker_df$Signal, levels = c("Positive", "Negative"))
marker_df$Gene <- factor(marker_df$Gene, levels = bulk_markers)
marker_df$Expression <- c(rep("Upregulated", 14 * length(samples)),
                          rep("Downregulated", 8 * length(samples)))
ggplot(marker_df, aes(x = Signal, y = Percentage, color = Signal)) + 
    geom_boxplot(width = 0.4, outlier.shape = NA, show.legend = FALSE) +
    geom_jitter(width = 0.2, alpha = 0.6) +
    facet_nested_wrap(. ~ Gene, scales = "free", ncol = 14) +
    theme_bw() +
    scale_color_manual(values = c(Positive = "black", Negative = "grey30")) +
    labs(x = "SBS2/SBS13", y = "Percentage of gene-expressing cells",    
         color = "SBS2/SBS13") +
    theme(axis.text = element_text(size = 10),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.border = element_rect(colour = "grey80"),
          strip.background = element_rect(colour = "grey80"))

# Perform statistical test for each gene
sig_results <- marker_df %>%
  group_by(Gene) %>%
  summarise(p_value = wilcox.test(Percentage ~ Signal)$p.value) %>%
  mutate(significance = ifelse(p_value < 0.05, "*", ""))
print(sig_results[sig_results$significance == "*",])



#--------------------------------- DE-analysis --------------------------------#

DefaultAssay(srat) <- "RNA"

## DE analysis APOBEC3A pos cells vs neg cells
apobec3a_positive_cells <- names(which(srat[["RNA"]]@data["APOBEC3A", ] > 0))
diff_expr <- FindMarkers(object = srat,
                         assay = "RNA",
                         ident.1 = apobec3a_positive_cells,
                         ident.2 = setdiff(Cells(srat), apobec3a_positive_cells))
#load(file = "R_cache/global_diff_expr.Rda")
 
# Define significance criteria
diff_expr$gene <- rownames(diff_expr)
diff_expr <- diff_expr %>%
  mutate(log_pval = -log10(p_val_adj), 
         significance = case_when(p_val_adj <= 0.05 & avg_log2FC >= 0.5 ~ "Upregulated",
                                  p_val_adj <= 0.05 & avg_log2FC <= -0.5 ~ "Downregulated",
                                  TRUE ~ "Not Significant"))
table(diff_expr$significance)

# Determine significant genes 
sign_genes <- diff_expr[!diff_expr$significance == "Not Significant", c("gene", "avg_log2FC", "p_val_adj", "significance")]
sign_genes %>% arrange(desc(significance), gene)

# Feature plot
up_sig_genes <- rownames(sign_genes[sign_genes$significance == "Upregulated",])
down_sig_genes <- rownames(sign_genes[sign_genes$significance == "Downregulated",])
FeaturePlot(srat, features = up_sig_genes, ncol = 4)
FeaturePlot(srat, features = down_sig_genes, ncol = 4)

# Determine nr of positive cells per DEG
deg_df <- data.frame()
for (gene in c(up_sig_genes, down_sig_genes)){
  
  # Determine all cells with gene expression
  expression <- ifelse(gene %in% up_sig_genes, "Upregulated", "Downregulated")
  gene_pos_bc <- names(which(srat[["RNA"]]@data[gene, ] > 0))
  for (i in 1:length(samples)){
    
    # Determine the percentage of sample cells with the DEG expression
    sample_bc <- names(which(srat$sample_assignment == samples[i]))
    n_cells <- length(sample_bc)
    n_pos_cells <- sum(sample_bc %in% gene_pos_bc)

    # Store in data frame
    deg_df <- rbind(deg_df, data.frame(Sample = samples[i],
                                       Gene = gene,
                                       Pos_perc = n_pos_cells * 100 / n_cells,
                                       signal = cluster_info$APOBEC[cluster_info$Random_ID == samples[i]],
                                       Expression = expression))
  }
}

# Generate barplot (19x5/8 in)
deg_df$signal <- factor(deg_df$signal, levels = c("Positive", "Negative"))
deg_df$Gene <- factor(deg_df$Gene, levels = c(sort(up_sig_genes), sort(down_sig_genes)))
for (expression in c("Upregulated", "Downregulated")){
  print(ggplot(deg_df[deg_df$Expression == expression,], aes(x = signal, y = Pos_perc, color = signal)) + 
    geom_boxplot(width = 0.4, outlier.shape = NA, show.legend = FALSE) +
    geom_jitter(width = 0.2, alpha = 0.6) +
    facet_nested_wrap(. ~ Gene, scales = "free", ncol = 17) +
    theme_bw() +
    scale_color_manual(values = c(Positive = "black", Negative = "grey30")) +
    labs(x = "SBS2/SBS13", y = "Percentage of gene-expressing cells",   
         title = paste(expression, "genes"), color = "SBS2/SBS13") +
    theme(axis.text = element_text(size = 10),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.border = element_rect(colour = "grey80"),
          strip.background = element_rect(colour = "grey80")))
}



#----------------------- DE-analysis with balanced groups ---------------------#

DefaultAssay(srat) <- "RNA"

# Select APOBEC3A+ cells among APOBEC-positive samples
apobec3a_pos_cells <- WhichCells(srat, expression = APOBEC3A > 0 & 
                                       apobec_status == "Positive")

# Select APOBEC3A- cells among APOBEC-negative samples with a min. transcr. count of 2000
apobec3a_neg_cells <- WhichCells(srat, expression = APOBEC3A == 0 & 
                                       apobec_status == "Negative" & 
                                       nCount_RNA >= 2000)

# Create balanced groups
apobec3a_neg_cells_balanced <- sample(apobec3a_neg_cells, 
                                      length(apobec3a_pos_cells), 
                                      replace = FALSE)

# Subset srat object for selected cells
srat_subset <- subset(srat, cells = c(apobec3a_pos_cells, 
                                      apobec3a_neg_cells_balanced))

# Perform normalization
srat_subset <- NormalizeData(srat_subset, 
                             normalization.method = "LogNormalize")
srat_subset <- ScaleData(srat_subset, 
                         features = rownames(srat_subset))
srat_subset <- FindVariableFeatures(srat_subset)
srat_subset <- SCTransform(srat_subset, 
                           vars.to.regres=NULL,
                           vst.flavor='v2', verbose = FALSE,
                           seed.use=CFG$random_seed)

# Perform DE-analysis
DefaultAssay(srat_subset) <- "RNA"
diff_expr <- FindMarkers(object = srat_subset,
                         assay = "RNA",
                         logfc.threshold = 0, min.pct = 0, # to obtain DE-results for all genes
                         ident.1 = apobec3a_pos_cells,
                         ident.2 = apobec3a_neg_cells_balanced)
#load("R_cache/diff_expr-balanced.Rda")
#diff_expr <- DE_list[[1]]; srat_subset <- DE_list[[2]]; apobec3a_pos_cells <- DE_list[[3]]; apobec3a_neg_cells_balanced <- DE_list[[4]]; rm(DE_list)
diff_expr <- diff_expr[rownames(diff_expr) %in% rownames(srat_subset),]

# Define significance criteria
diff_expr$gene <- rownames(diff_expr)
diff_expr <- diff_expr %>%
  mutate(log_pval = -log10(p_val_adj), 
    significance = case_when(p_val_adj < 0.05 & avg_log2FC >= 0.5 ~ "Upregulated",
                             p_val_adj < 0.05 & avg_log2FC <= -0.5 ~ "Downregulated",
                             TRUE ~ "Not significant"))
table(diff_expr$significance)

# Cap infinite values
diff_expr$log_pval[!is.finite(diff_expr$log_pval)] <- 145

# Determine the top DEGs (p_val_adj < 0.01 & abs(log2FC > 1))
top_genes <- diff_expr %>% filter((p_val_adj <= 0.01 & abs(avg_log2FC) >= 1))
label_genes <- rbind(top_genes, diff_expr %>% filter(log_pval > 65)) %>% unique()

# Create volcano plot (8.5x7 in)
ggplot(diff_expr, aes(x = avg_log2FC, y = log_pval, color = significance)) +
  geom_point(alpha = 0.8, size = 2) + 
  geom_text_repel(data = label_genes, aes(label = gene), size = 2.2, show.legend = FALSE) + 
  scale_x_continuous(limits = c(-2, 2)) + scale_y_continuous(limits = c(0, 150)) +
  scale_color_manual(values = c("Upregulated" = "darkgreen", "Downregulated" = "blue", 
                                "Not significant" = "gray"),
                     breaks = c("Upregulated", "Downregulated", "Not significant")) +
  labs(title = "Volcano plot, scRNAseq data", x = "Log2 fold change",
       y = "-Log10 adjusted p-value", color = "Significance") +
  theme_classic() +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "black") +
  theme(plot.title = element_text(face = "bold"),
        axis.text = element_text(size = 11),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 10.2))

# Write significant genes to out file
sign_genes <- diff_expr[!diff_expr$significance == "Not significant", 
                        c("gene", "p_val", "p_val_adj", "avg_log2FC", "significance")]
sign_genes <- sign_genes %>% arrange(desc(significance), gene)
write.table(sign_genes, file = "DEanalysis-balanced-APOBECposVSneg.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)



#----------------------------- Correlation analysis ---------------------------#

DefaultAssay(srat_subset) <- "RNA"

## Correlation analysis
# Extract APOBEC3A expression
apobec3a_expr <- srat_subset@assays$RNA@data["APOBEC3A", ]

# Calculate correlations with all other genes
gene_correlations <- apply(srat_subset@assays$RNA@data, 1, 
    function(gene_expr) {cor(gene_expr, apobec3a_expr, method = "pearson")})

# Sort genes by correlation values
sorted_genes <- sort(gene_correlations, decreasing = TRUE)

# View top positively and negatively correlated genes
head(sorted_genes)
tail(sorted_genes)



#--------------------------- Module Scores for DEGs ---------------------------#

DefaultAssay(srat) <- "RNA"

# Parse bulk RNAseq DE-analysis results
bulk_DEGs <- as.data.frame(read_excel("DEGs-ETV6RUNX1samples_complete_121125.xlsx"))
bulk_DEGs$gene <- as.character(sapply(bulk_DEGs$Gene, function(x) 
  {ifelse(x %in% rownames(srat), x, rownames(srat)[grepl(x, rownames(srat))])}))
colnames(bulk_DEGs)[4] <- "significance"
bulk_DEGs <- bulk_DEGs[!is.na(bulk_DEGs$gene),]

# Subset for genes with abs(LFC) >= 2 and p-adj < 0.01
subset_DEGs <- subset(bulk_DEGs, padj <= 0.01 & abs(LFC) >= 2)

# Define barcodes of APOBEC3A positive cells
apobec3a_pos_cells_complete <- WhichCells(srat, expression = APOBEC3A > 0)
  
# Calculate a ModuleScore for DEGs for visualization purposes
#! Note: change diff_expr to top_genes for top DEGs module score analysis
#        or to subset_DEGs for bulk module score analysis
up_deg_list <- diff_expr[diff_expr$significance == "Upregulated", "gene"]
down_deg_list <- diff_expr[diff_expr$significance == "Downregulated", "gene"]
srat <- AddModuleScore(srat, 
                       assay = "RNA",
                       features = list(up_deg_list, down_deg_list), 
                       name = c("UpDEGs", "DownDEGs"))
colnames(srat@meta.data)[grepl("DEG", colnames(srat@meta.data))] <- c("UpregulatedGenes", "DownregulatedGenes")

# Define the top 10% quantile of module scores
up_cutoff <- quantile(srat$UpregulatedGenes, probs = 0.9)
down_cutoff <- quantile(srat$DownregulatedGenes, probs = 0.1)

# Add columns whether barcode meets cutoff scores
srat@meta.data <- srat@meta.data %>%
  mutate(`ParsedModuleScore\nUpregulatedGenes` = ifelse(UpregulatedGenes >= up_cutoff, 1, 0),
         `ParsedModuleScore\nDownregulatedGenes` = ifelse(DownregulatedGenes <= down_cutoff, 1, 0))


## Create violin plots (7.5x6 in)
p <- FetchData(srat, vars = c("UpregulatedGenes", "DownregulatedGenes", "ident")) 
p$color_group <- ifelse(rownames(p) %in% apobec3a_pos_cells_complete, "black", "grey80")
p <- p %>% arrange(color_group == "black")
p1 <- ggplot(p, aes(x = ident, y = UpregulatedGenes)) +
  geom_jitter(aes(color = color_group), width = 0.2, size = 0.3) +
  geom_violin(fill = NA, color = "grey50", width  = 0.55) +
  scale_color_manual(values = c("grey80" = "grey80", "black" = "black")) +
  labs(x = "All cells", y = "ModuleScore", title = "ModuleScore\nUpregulated DEGs") +
  geom_hline(yintercept = up_cutoff, linetype = "dashed", color = "blue", linewidth = 0.9) +
  theme_minimal() + theme(legend.position = "none", axis.text.x = element_blank(),
                          plot.title = element_text(hjust = 0.5, face = "bold", size = 11.5))
p2 <- ggplot(p, aes(x = ident, y = DownregulatedGenes)) +
  geom_jitter(aes(color = color_group), width = 0.2, size = 0.3) +
  geom_violin(fill = NA, color = "grey50", width = 0.55) +
  scale_color_manual(values = c("grey80" = "grey80", "black" = "black")) +
  labs(x = "All cells", y = "ModuleScore", title = "ModuleScore\nDownregulated DEGs") +
  geom_hline(yintercept = down_cutoff, linetype = "dashed", color = "blue", linewidth = 0.9) +
  theme_minimal() + theme(legend.position = "none", axis.text.x = element_blank(),
                          plot.title = element_text(hjust = 0.5, face = "bold", size = 11.5))
plot_grid(p1, p2)


## Relative expression overview
# Create a FeaturePlot separated by APOBEC status
srat$apobec_status_plot <- srat$apobec_status
srat$apobec_status_plot[srat$apobec_status_plot == "Positive"] <- "SBS2/SBS13 positive"
srat$apobec_status_plot[srat$apobec_status_plot == "Negative"] <- "SBS2/SBS13 negative"
srat$apobec_status_plot <- factor(srat$apobec_status_plot, 
                                  levels = c("SBS2/SBS13 positive", "SBS2/SBS13 negative"))
rel_overview <- FeaturePlot(srat, 
                            reduction = "umap", 
                            features = c("APOBEC3A", "ParsedModuleScore\nUpregulatedGenes", 
                                         "ParsedModuleScore\nDownregulatedGenes"), 
                            split.by = "apobec_status_plot",
                            order = TRUE, 
                            cols = CFG$gradient.colors)

# Extract legend and combine figures (10x9 in)
rel_legend <- get_legend(FeaturePlot(srat, 
                                     reduction = "umap", 
                                     features = "APOBEC3A",
                                     cols = CFG$gradient.colors) + 
                                    labs(color = "Relative\nExpression") +
                                    theme(legend.title = element_text(size = 12)))
plot_grid(rel_overview, rel_legend, rel_widths = c(1,0.1))


## Correlation plot (7x6 in)
df <- FetchData(srat, vars = c("UpregulatedGenes", "DownregulatedGenes"))
df$Color <- ifelse(df$UpregulatedGenes >= up_cutoff & df$DownregulatedGenes <= down_cutoff,
                   "yes", "no")
nrow(df[rownames(df) %in% apobec3a_pos_cells_complete & df$Color == "yes", ])
df$Color[rownames(df) %in% apobec3a_pos_cells_complete] <- "apobec"
df <- df %>% arrange(Color == "apobec")
p3 <- ggplot(df, aes(x = DownregulatedGenes, y = UpregulatedGenes, color = Color)) +
  geom_point(alpha = 0.5, size = 0.6) +  
  annotate("segment", x = down_cutoff, xend = down_cutoff, y = up_cutoff, yend = 1.1, # use yend=2.1 for top_genes and yend=0.29 for bulk genes
           color = "blue", linetype = "dashed", linewidth = 1) +
  annotate("segment", x = -0.4, xend = down_cutoff, y = up_cutoff, yend = up_cutoff,  # use x=-1.6 for top_genes and x=-0.12 for bulk genes
           color = "blue", linetype = "dashed", linewidth = 1) +
  scale_y_continuous(limits = c(-0.6, 1.1), expand = c(0.001, 0.001)) +           # use limits c(-1.6, 2.1) for top_genes and c(-0.15, 0.29) for bulk genes
  scale_x_continuous(limits = c(-0.4, 1.3), expand = c(0.001, 0.001)) +           # use limits c(-1.6, 2.2) for top_genes and c(-0.12, 0.33) for bulk genes
  theme_minimal() +
  scale_color_manual(values = c(yes = "lightgrey", no = "lightgrey", apobec = "black")) +
  labs(title = "ModuleScore Comparison", y = "Upregulated Module Score", x = "Downregulated Module Score") +
  theme(legend.position = "none")
p3

# Remove ModuleScores from meta data
srat@meta.data <- srat@meta.data[,-c(which(colnames(srat@meta.data) == "UpregulatedGenes"):ncol(srat@meta.data))]

# Create combined figure for manuscript (11x5.5 in)
plot_grid(p1 + labs(title = "Upregulated DEGs") + theme(plot.title = element_text(size = 11.5)),
          p2 + labs(title = "Downregulated DEGs") + theme(plot.title = element_text(size = 11.5)),
          p3 + theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 11.5)), 
          nrow = 1, rel_widths = c(1, 1, 2.3))



#---------------------------- Enrichment analysis -----------------------------#

# Fetch conversion table from gene symbol to ENTREZ IDs
Sys.setenv(ANNOTATIONHUB_CACHE = "/hpc/pmc_kuiper/Rstudio_Server_Libs/Rstudio_4.3.1_libs/AnnotationHub")
conversion_table <- fetch_conversion_table(organism_name = "Homo sapiens", 
                                           from = "SYMBOL", 
                                           to = "ENTREZID")
#conversion_table <- qread("/hpc/pmc_kuiper/HypermutatedALL_project/CODE/SHELL/scRNAseq/MSnID-conversionTable.qs")

# Add entrez IDs to DE-analysis result
diff_expr$gene <- sapply(strsplit(rownames(diff_expr), "\\."), `[`, 1)
diff_expr <- as.data.frame(diff_expr) %>% 
  mutate(SYMBOL = gene) %>% 
  left_join(conversion_table, by = "SYMBOL")

# Remove genes without an ENTREZ ID
diff_expr <- diff_expr[!is.na(diff_expr$ENTREZID),]

# Conduct overrepresentation test among DEGs
up_genes_entrez <- diff_expr[diff_expr$significance == "Upregulated", "ENTREZID"]
kegg_up <- enrichKEGG(gene = up_genes_entrez,
                      universe = diff_expr$ENTREZID, 
                      organism = "hsa", 
                      pvalueCutoff = 0.05)
down_genes_entrez <- diff_expr[diff_expr$significance =="Downregulated", "ENTREZID"]
kegg_down <- enrichKEGG(gene = down_genes_entrez,
                        universe = diff_expr$ENTREZID, 
                        organism = "hsa", 
                        pvalueCutoff = 0.05)
dotplot(kegg_up) + ggtitle("KEGG overrepresentation test, upregulated")
dotplot(kegg_down) + ggtitle("KEGG overrepresentation test, downregulated")

# Get gene info of genes present in key pathways
gene_ids <- c()
for (path in c("Human T-cell leukemia virus 1 infection", 
               "Antigen processing and presentation")){
  gene_ids <- c(gene_ids, kegg_up@result$geneID[kegg_up@result$Description == path])
}
gene_freq <- as.data.frame(sort(table(unlist(strsplit(gene_ids, "/"))), decreasing = TRUE))
names(gene_freq) <- c("ENTREZ", "Freq")
gene_freq["Gene"] <- diff_expr$SYMBOL[match(gene_freq$ENTREZ, diff_expr$ENTREZID)]


## GSEA
# Create a named ENTREZ ID vector of log2FC for GSEA
gene_ranks_entrez <- diff_expr %>%
  filter(!is.na(ENTREZID), !is.na(avg_log2FC)) %>% 
  #mutate(ranking_metric = -log10(p_val)*sign(avg_log2FC)) %>% 
  mutate(ranking_metric = avg_log2FC) %>%
  group_by(ENTREZID) %>% 
  summarise(ranking_metric = mean(ranking_metric, na.rm = TRUE)) %>% 
  tibble::deframe() 

# Sort by rank
gene_ranks_entrez <- sort(gene_ranks_entrez, decreasing = TRUE)

# Conduct KEGG pathway enrichment analysis
kegg_gsea <- gseKEGG(geneList = gene_ranks_entrez,
                     organism = "hsa",
                     pvalueCutoff = 0.05,
                     pAdjustMethod = "fdr",
                     keyType = "ncbi-geneid")

# Visualise enriched KEGG pathways (10x8 in)
dotplot(kegg_gsea, 
        x = "GeneRatio",
        showCategory = 100, 
        font.size = 12, 
        label_format = 40,
        split = ".sign") + 
  facet_grid(.sign ~ ., scales = "free", space = "free") +
  theme(strip.text.y = element_text(size = 10.5),
        axis.title = element_text(size = 11.5),
        axis.text.x = element_text(size = 11.5),
        axis.text.y = element_text(size = 11.5))



#--------------------------- Single-sample DE-analysis ------------------------#

select_samples <- c("P0773", "P0763", "P0757", "P0679", "P0745", "P0677")
diff_expr_df <- data.frame()
for (sample in select_samples){
  
  # Obtain barcodes with- and without APOBEC3A expression
  sample_bc <- names(which(srat$sample_assignment == sample))
  pos_cells <- sample_bc[sample_bc %in% apobec3a_pos_cells_complete]
  neg_cells <- sample_bc[!sample_bc %in% apobec3a_pos_cells_complete]
  
  # Create balanced groups
  neg_cells_balanced <- sample(neg_cells, length(pos_cells), replace = FALSE)
  sample_srat_subset <- subset(srat, cells = c(pos_cells, neg_cells))
  
  # Perform normalization
  sample_srat_subset <- NormalizeData(sample_srat_subset, 
                                      normalization.method = "LogNormalize")
  sample_srat_subset <- ScaleData(sample_srat_subset, 
                                  features = rownames(sample_srat_subset))
  sample_srat_subset <- FindVariableFeatures(sample_srat_subset)
  
  # Perform DE-analysis
  DefaultAssay(sample_srat_subset) <- "RNA"
  sample_diff_expr <- FindMarkers(object = sample_srat_subset,
                                  assay = "RNA",
                                  ident.1 = pos_cells,
                                  ident.2 = neg_cells)
  
  # Define significance criteria
  sample_diff_expr$gene <- rownames(sample_diff_expr)
  sample_diff_expr <- sample_diff_expr %>%
    mutate(log_pval = -log10(p_val_adj), 
           significance = case_when(p_val_adj < 0.05 & avg_log2FC >= 0.5 ~ "Upregulated",
                                    p_val_adj < 0.05 & avg_log2FC <= -0.5 ~ "Downregulated",
                                    TRUE ~ "Not significant"))
  
  # Determine significant genes
  sign_genes <- sample_diff_expr[!sample_diff_expr$significance == "Not significant", c("gene", "p_val_adj", "avg_log2FC", "significance")]
  sign_genes$sample <- sample
  
  # Add to data frame
  diff_expr_df <- rbind(diff_expr_df, sign_genes)
}

table(diff_expr_df[diff_expr_df$significance == "Upregulated", "gene"])
table(diff_expr_df[diff_expr_df$significance == "Downregulated", "gene"])



#-------------------------------- Miscellaneous -------------------------------#

### Identify genes that explain differences between 2 clusters
# Obtain srat info from sample
sample_srat <- subset(srat, sample_assignment == samples[i])
sample_srat <- RunPCA(sample_srat, npcs = 30)
sample_srat <- RunUMAP(sample_srat, dims = 1:CFG$ndims, 
                       seed.use = CFG$random_seed, n.neighbors = 20)
sample_srat$seurat_clusters <- as.numeric(as.character(sample_srat$seurat_clusters))
DimPlot(sample_srat, group.by = "seurat_clusters", label = TRUE) 

# Define seurat cluster numbers for both sides
left_cells <- names(lib_srat$seurat_clusters)[lib_srat$seurat_clusters %in% c(0, 1, 3, 4, 5, 6, 7, 10, 11, 13, 14, 16, 19)]
right_cells <- names(lib_srat$seurat_clusters)[lib_srat$seurat_clusters %in% c(12, 17, 20)]
lib_srat$Side <- ifelse(names(lib_srat$seurat_clusters) %in% left_cells, "Left", "Right")
table(lib_srat$Side)
DimPlot(lib_srat, group.by = "Side", label = TRUE) + NoLegend()

# Determine DEGs
markers_left_vs_right <- FindMarkers(object = lib_srat,
                                     ident.1 = left_cells,
                                     ident.2 = right_cells)
cluster_markers <- markers_left_vs_right
top_genes <- head(rownames(cluster_markers), 20)
up_sig_genes <- rownames(cluster_markers[cluster_markers$p_val_adj < 0.05 & cluster_markers$avg_log2FC > 1, ])
down_sig_genes <- rownames(cluster_markers[cluster_markers$p_val_adj < 0.05 & cluster_markers$avg_log2FC < -1, ])


### Determine the nr of transcripts per cell
# Extract the transcript count per gene per cell
raw_counts <- GetAssayData(srat, slot = "counts", assay = "RNA")
norm_counts <- GetAssayData(srat, slot = "data", assay = "RNA")

# Sum the counts for each barcode/cell
norm_transcript_counts_per_cell <- colSums(norm_counts)

# Add the transcript counts as metadata
srat <- AddMetaData(srat, 
                    metadata = colSums(raw_counts), 
                    col.name = "TranscriptCount") # this column is identical to nCount_RNA
srat <- AddMetaData(srat, 
                    metadata = norm_transcript_counts_per_cell, 
                    col.name = "NormTranscriptCount")

# Visualize
VlnPlot(srat, "TranscriptCount", group.by = "library_assignment") & 
  labs(x = NULL, x = "library_assignment")
VlnPlot(srat, "NormTranscriptCount", group.by = "library_assignment") & 
  labs(x = NULL, x = "library_assignment")
VlnPlot(srat, features = c("nFeature_RNA",  "TranscriptCount", "NormTranscriptCount"))



#--------------------- APOBEC3A transcript count inspection -------------------#

## Determine the nr of transcripts vs the nr of APOBEC3A per cell
# Get transcript counts for APOBEC3A
transcript_count <- srat@meta.data$nCount_RNA
apobec3a_raw_counts <- GetAssayData(srat, slot = "counts", assay = "RNA")["APOBEC3A", ]
apobec3a_norm_counts <- GetAssayData(srat, slot = "data", assay = "RNA")["APOBEC3A", ]

# Combine with total transcript counts in a data frame
plot_data <- data.frame(TranscriptCount = transcript_count,
                        NormTranscriptCount = norm_transcript_counts_per_cell,
                        APOBEC3A_Counts = apobec3a_raw_counts,
                        APOBEC3A_Norm_Counts = apobec3a_norm_counts)

# Visualize
ggplot(plot_data, aes(x = TranscriptCount, y = APOBEC3A_Counts)) +
  geom_point(alpha = 0.6, color = "black") + 
  theme_minimal() +                    
  scale_y_continuous(breaks = seq(0, 17, 1)) +
  labs(title = "Raw transcript count vs raw APOBEC3A counts per cell",
       x = "Total transcript count",
       y = "APOBEC3A transcript count") + theme(panel.grid.minor = element_blank())
ggplot(plot_data, aes(x = NormTranscriptCount, y = APOBEC3A_Norm_Counts)) +
  geom_point(alpha = 0.6, color = "black") +
  theme_minimal() +                        
  labs(title = "Normalized transcript count vs normalized APOBEC3A counts per cell",
       x = "Total transcript count",
       y = "APOBEC3A transcript count") + theme(panel.grid.minor = element_blank())



## Determine genes with a corresponding transcript count per cell as APOBEC3A
corresp_genes <- data.frame()
for (i in 1:length(samples)){
  sample_srat <- subset(srat, sample_assignment == samples[i])
  sample_norm_counts_per_gene <- rowSums(GetAssayData(sample_srat, slot = "data", assay = "RNA"))
  
  # Obtain APOBEC3A norm counts
  sample_apobec_3a <- sample_norm_counts_per_gene[["APOBEC3A"]]
  lower_limit <- sample_apobec_3a * 0.9
  upper_limit <- sample_apobec_3a * 1.1
  
  # Check which genes have a corresponding normalized count as APOBEC3A
  select_genes <- sample_norm_counts_per_gene[sample_norm_counts_per_gene <= upper_limit &
                                              sample_norm_counts_per_gene >= lower_limit]
  # Add to data frame
  corresp_genes <- rbind(corresp_genes, data.frame(Gene = names(select_genes),
                                                   Norm_count = as.numeric(select_genes),
                                                   Sample = samples[i]))
}

# Select genes that are found in at least 4 samples
select_genes <- sort(table(corresp_genes$Gene), decreasing = T)
select_genes <- names(select_genes[select_genes > 3])

# Determine the percentage of cells of each sample with target_gene expression
select_gene_df <- data.frame()
for (gene in select_genes){
  gene_pos_bc <- names(which(srat[["RNA"]]@data[gene, ] > 0))
  
  # Determine the nr of positive cells per sample
  for (i in 1:length(samples)){
    sample_bc <- names(which(srat$sample_assignment == samples[i]))
    n_cells <- length(sample_bc)
    n_pos_cells <- sum(sample_bc %in% gene_pos_bc)
    select_gene_df <- rbind(select_gene_df, data.frame(Sample = samples[i],
                                                       Gene = gene,
                                                       Gene_perc = n_pos_cells * 100 / n_cells,
                                                       n_total = n_cells, 
                                                       n_gene = n_pos_cells,
                                                       signal = cluster_info$APOBEC[cluster_info$Random_ID == samples[i]]))
  }
}

# Create a boxplot of percentage of cells that express the target genes
select_gene_df_melt <- melt(select_gene_df[,c("Gene", "Gene_perc", "signal")])
select_gene_df_melt$signal <- factor(select_gene_df_melt$signal, levels = c("Positive", "Negative"))
ggplot(select_gene_df_melt, aes(x = signal, y = value)) + 
  geom_boxplot(width = 0.4, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6, color = "black") +
  facet_wrap(. ~ Gene, ncol = 12, scales = "free") +
  theme_bw() +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  labs(x = "SBS2/SBS13", y = "Percentage of cells",    
       title = "Gene expression") +
  theme(axis.text = element_text(size = 10),
        panel.border = element_rect(colour = "grey80"),
        strip.background = element_rect(colour = "grey80"))

# Perform statistical test for each gene
sig_results <- select_gene_df_melt %>%
  group_by(Gene) %>%
  summarise(p_value = wilcox.test(value ~ signal)$p.value) %>%
  mutate(significance = ifelse(p_value < 0.05, "*", ""))
print(sig_results[sig_results$significance == "*",])

