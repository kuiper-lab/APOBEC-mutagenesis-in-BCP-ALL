## ----prepare, include=TRUE----------------------------------------------
rm(list = ls())

# Load libraries
library(Seurat) 
library(qs)
library(readr)
library(readxl)
library(ggplot2)
library(gridExtra)
library(patchwork)

# List global variables of library LX343
RESULTS_DIR <- "/hpc/pmc_kuiper/HypermutatedALL_project/RESULTS/scRNAseq/"
SEURAT_OBJECT_LX343 <- "LX343-final.qs"

# List global variables of all combined libraries
MERGED_SEURAT_OBJECT <- "APOBEC-final.qs"
MERGED_CLUSTERS <- "merged_lib_clusters-111024.tsv"
CLUSTER_META_DATA <- "scRNAseq-APOBEC-clusterAssignment-111024.xlsx"
setwd(RESULTS_DIR)

# Load the Seurat objects and cluster assignments
lib_srat <- qread(SEURAT_OBJECT_LX343)
srat <- qread(MERGED_SEURAT_OBJECT)
snp_clusters <- as.data.frame(read_table(MERGED_CLUSTERS))
rownames(snp_clusters) <- snp_clusters$barcode

# Load cluster meta data
cluster_info <- as.data.frame(read_excel(CLUSTER_META_DATA))

# Create objects to store sample information
library_assignment <- c()
sample_assignment <- c()
barcode_status <- c()

# Add sample information to each cell
for (n_row in 1:nrow(lib_srat@meta.data)) {
  
  # Load cell barcode
  current_row <- lib_srat@meta.data[n_row, ]
  barcode <- rownames(current_row)
  snp_cluster <- snp_clusters[barcode,]
  
  # Check if the cell is assigned to multiple clusters
  if (grepl("/", snp_cluster$assignment)){
    lib_cluster <- cluster <- snp_cluster$assignment
    
    # Add sample ID to the cluster  
  } else {
    lib_cluster <- paste(unlist(strsplit(snp_cluster$barcode, "_"))[1], 
                         snp_cluster$assignment, 
                         sep = "_")
    cluster_meta <- cluster_info[match(lib_cluster, cluster_info$Cluster), ]
    cluster <- cluster_meta$Sample_ID
  }
  
  # Add info to object
  library_assignment <- append(library_assignment, lib_cluster)
  sample_assignment <- append(sample_assignment, cluster)
  barcode_status <- append(barcode_status, snp_cluster$status)
}

# Add meta data to Seurat object
lib_srat@meta.data$library_assignment <- library_assignment
lib_srat@meta.data$sample_assignment_compl <- sample_assignment
lib_srat@meta.data$barcode_status <- barcode_status

# Parse cell annotations
lib_srat$Sample_Descr <- "NA"
for (i in 1:nrow(lib_srat@meta.data)){
  lib_row <- lib_srat@meta.data[i,]
  
  if (lib_row$sample_assignment == "NALM-6"){
    sample_descr <- "NALM-6 cell"
  } else if (lib_row$barcode_status == "singlet"){
    if (grepl("^P", lib_row$sample_assignment)){
      sample_descr <- "Patient cell"
    } else {
      sample_descr <- "NALM-6 cell"
    }
  } else {
    if (lib_row$library_assignment %in% c("0/3", "3/0")){
      sample_descr <- "Doublet among NALM-6"
    } else if (grepl("^0/|^3/", lib_row$library_assignment)) {
      sample_descr <- "Doublet among NALM-6/patient cells"
    } else {
      sample_descr <- "Doublet among patient cells"
    }
  }
  lib_srat@meta.data[i,"Sample_Descr"] <- sample_descr
}
sample_assignment[grepl("/", sample_assignment)] <- "Mixed clusters"
lib_srat$sample_assignment <- sample_assignment



#---------------------- Define barcodes of NALM-6 cells -----------------------#

lib_samples <- unique(cluster_info$Sample_ID[grep("LX343", cluster_info$Cluster)])


## Run PCA and UMAP
lib_srat <- RunPCA(lib_srat, npcs = 50)
lib_srat <- RunUMAP(lib_srat, 
                    dims = 1:25, 
                    seed.use = 2033)

# Extract UMAP coordinates
umap_coords <- Embeddings(lib_srat, reduction = "umap")
umap_df <- as.data.frame(umap_coords)


## Generate library overview plots
p_sample <- DimPlot(lib_srat, 
                    pt.size = 1, 
                    group.by = "sample_assignment") 
p_phase <- DimPlot(lib_srat, 
                   pt.size = 1, 
                   group.by = "Phase")
p_cells <- DimPlot(lib_srat, 
                   pt.size = 1, 
                   group.by = "celltype") 
p_sample | p_phase | p_cells

# Generate plots for each sample in the target library
plots <- rep(list(NA), length(lib_samples))
for (i in 1:length(lib_samples)){
  sample <- lib_samples[i]
  
  # Get Seurat object for sample
  sample_srat <- subset(lib_srat, sample_assignment == sample)
  
  # Generate plots
  p_phase_sample <- DimPlot(sample_srat, 
                            pt.size = 1, 
                            group.by = "Phase") + labs(title = sample)
  p_cells_sample <- DimPlot(sample_srat, 
                            pt.size = 1, 
                            group.by = "celltype") + labs(title = paste("Celltype", sample))
  plots[[i]] <- arrangeGrob(p_phase_sample, p_cells_sample, nrow = 1)
}
do.call("grid.arrange", plots)


## Generate plots displaying cell annotations
gradient_cols <- viridis::viridis(101, direction = -1)

# Complete overview
p_descr <- DimPlot(lib_srat, 
                   pt.size = 0.5, 
                   group.by = "Sample_Descr") 
p_descr | p_cells

# NALM-6 specific overview
nalm6_srat <- subset(lib_srat, Sample_Descr %in% c("Doublet among NALM-6/patient cells", "Doublet among NALM-6", "NALM-6 cell"))
p_nalm6 <- DimPlot(nalm6_srat, 
                   pt.size = 0.5, 
                   group.by = "Sample_Descr") 
p_cells_nalm6 <- DimPlot(nalm6_srat, 
                         reduction = "umap", # you can leave this out, default
                         pt.size = 0.5, 
                         group.by = "celltype")
p_apobec_nalm6 <- FeaturePlot(nalm6_srat, 
                              pt.size = 1, 
                              features = "APOBEC3A", 
                              cols = gradient_cols,
                              order = TRUE) + labs(title = "APOBEC3A")
p_nalm6 | p_cells_nalm6 | p_apobec_nalm6

# Patient specific overview
patient_srat <- subset(lib_srat, Sample_Descr %in% c("Doublet among patient cells", "Patient cell"))
p_patient <- DimPlot(patient_srat, 
                   pt.size = 0.5, 
                   group.by = "Sample_Descr") 
p_cells_patient <- DimPlot(patient_srat, 
                         reduction = "umap", # you can leave this out, default
                         pt.size = 0.5, 
                         group.by = "celltype")
p_apobec_patient <- FeaturePlot(patient_srat, 
                              pt.size = 1, 
                              features = "APOBEC3A", 
                              cols = gradient_cols,
                              order = TRUE) + labs(title = "APOBEC3A")
p_patient | p_cells_patient | p_apobec_patient


## Obtain barcodes from NALM-6 cells
nalm6_annotated_cells <- rownames(nalm6_srat@meta.data)
nalm6_cluster_cells <- rownames(umap_df[umap_df$UMAP_1 < -6,])

# Check that the used UMAP coords are correct for determining the NALM-6 cluster
lib_srat$Cluster <- ifelse(rownames(lib_srat@meta.data) %in% nalm6_cluster_cells, "NALM-6", "NA")
DimPlot(lib_srat, group.by = "Cluster", label = TRUE) + NoLegend()

# Create vector of final selection of cell line-specific barcodes
nalm6_barcodes <- unique(c(nalm6_cluster_cells, nalm6_annotated_cells))

# Visualize location of NALM-6 cells in the UMAP overview
lib_srat$Cluster <- ifelse(rownames(lib_srat@meta.data) %in% nalm6_barcodes, "NALM-6", "Patients")
DimPlot(lib_srat, group.by = "Cluster", label = TRUE) + NoLegend()



#--------------------- Subset data sets for patient cells -------------------#

# Write patient-specific barcodes of LX343 to out file
LX343_patient_barcodes <- rownames(lib_srat@meta.data[!rownames(lib_srat@meta.data) %in% nalm6_barcodes, ])
LX343_patient_barcodes <- paste0(gsub("LX343_", "", LX343_patient_barcodes), "-1")
write.table(LX343_patient_barcodes, file = "LX343_nonCellLineBarcodes.tsv", 
            quote = FALSE, row.names = FALSE, col.names = FALSE)

# Remove NALM-6 cells from Seurat object of all combined libraries
cells_to_keep <- rownames(srat@meta.data[!rownames(srat@meta.data) %in% nalm6_barcodes, ])
srat <- subset(srat, cells = cells_to_keep)
qsave(srat, file = "APOBEC-final-cellineFiltered-121224.qs")



## Subsequent steps
# 1. Gzip file LX343_nonCellLineBarcodes.tsv and use this as the barcode file for 
#    cluster assignment with Souporcell (k=4). Script: generateSouporcellJobs.sh
# 2. Use output file cluster_genotypes.vcf to determine the corresponding genotype 
#    of each cluster using VCF-eval. Adjust the meta data Excel file accordingly
# 3. Merge output file clusters.tsv file with the cluster files of the other 
#    libraries using script scRNAanalysis_combineClusterFiles.R
# 4. Continue with scRNAseq analysis using the meta data Excel file, the merged 
#    cluster files and file APOBEC-final-cellineFiltered-121224.qs as input 
#    parameters. Script: APOBEC-scRNAseq-analysis.R

