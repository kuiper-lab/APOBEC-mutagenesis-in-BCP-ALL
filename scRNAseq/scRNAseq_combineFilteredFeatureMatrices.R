#!/usr/bin/env Rscript

# Load libraries
library(Seurat) 
library(SCutils) 

# List global variables
PROJECT_DIR <- "/hpc/pmc_kuiper/HypermutatedALL_project/RESULTS/scRNAseq"
COUNT_MATRIX_EXT <- "_filtered_feature_bc_matrix.h5"
setwd(PROJECT_DIR)

# List the input count matrices
count_matrices <- list.files(path = PROJECT_DIR,
                             pattern = COUNT_MATRIX_EXT,
                             full.names = TRUE,
                             recursive = TRUE)

# Combine all count matrices (run on HPC due to memory overload)
COMBINE <- TRUE
if (COMBINE){
  for (i in 1:length(count_matrices)){
    library_id <- gsub(COUNT_MATRIX_EXT, "", basename(count_matrices[i]))
    
    # Read the library count matrix and add library ID to cell barcodes
    lib_content <- SCutils::read10xlibs(count_matrices[i])
    colnames(lib_content) <- gsub("scRNAseq", library_id,
                                  colnames(lib_content))
    print(colnames(lib_content))

    # Merge count matrix
    if (i == 1){
      merged_lib_contents <- lib_content
    } else {
      merged_lib_contents <- merge(merged_lib_contents, lib_content)
    }
    rm(lib_content)
  }
  save(merged_lib_contents, file = "merged_filtered_feature_matrices_220724.Rda")
  
} else {
  load("merged_filtered_feature_matrices_220724.Rda")
}

