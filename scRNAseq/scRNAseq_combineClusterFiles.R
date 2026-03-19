#------------------------------ Data preparations ----------------------------#

rm(list = ls())

# Load libraries
library("readr")

# List global variables
PROJECT_DIR <- "/Users/m.m.kleisman/Projects/HypermutatedALL_project/"
RESULTS_DIR <- paste0(PROJECT_DIR, "/RESULTS/scRNAseq/clusters/")
OUT_FILE <- "merged_lib_clusters-cellineFiltered-121224.tsv"
CLUSTERS_EXT <- "-clusters.tsv"
setwd(RESULTS_DIR)



#------------------------------ Combine clusters ------------------------------#

# List the input clusters
cluster_files <- list.files(pattern = CLUSTERS_EXT,
                            full.names = TRUE,
                            recursive = FALSE)

combined_lib_clusters <- data.frame()
for (i in 1:length(cluster_files)){
  library_id <- unlist(strsplit(basename(cluster_files[i]), "-"))[1]
  
  # Read the library clusters 
  lib_clusters <- as.data.frame(read_tsv(cluster_files[i]))[,c("barcode", "status", "assignment")]
  
  # Add library ID to cell barcodes
  lib_clusters$barcode <- paste0(library_id, "_", gsub("-1", "", lib_clusters$barcode))
  
  # Define the log loss value of each assigned cluster
  #lib_clusters["log_loss"] <- NA
  #for (n_row in 1:nrow(lib_clusters)){
  #  assignment <- lib_clusters[n_row, "assignment"]
  #  log_loss <- ifelse(grepl("/", assignment), 
  #                     lib_clusters[n_row, "log_prob_doublet"],
  #                     lib_clusters[n_row, paste0("cluster", assignment)])
  #  lib_clusters[n_row, "log_loss"] <- log_loss
  #}
  
  # Note: log_prob_singleton is identical to the log loss value of each assigned cluster
  # (doublets will be ignored anyways in subsequent analyses)
  #lib_clusters["log_loss"] <- lib_clusters$log_prob_singleton
  #lib_clusters["probability"] <- exp(lib_clusters$log_loss)
  
  # Store in object
  combined_lib_clusters <- rbind(combined_lib_clusters, lib_clusters)
}
write.table(combined_lib_clusters, file = OUT_FILE,
            quote = FALSE, col.names = TRUE, row.names = FALSE)


