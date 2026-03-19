#---------------------------- Data preparations -----------------------------#

# Clear data 
rm(list = ls())

# Load libraries
library("VariantAnnotation")
library("logging")

# Set parameters and the path to data directories
PROJECT_DIR <- "/Users/m.m.kleisman/Projects/HypermutatedALL_project/"
VCF_FILE_EXT <- "_WGS.HyperExomeRegions.vcf.gz"
TAB_FILE_EXT <- "_WGS.HyperExomeRegions.tab.txt"
OUT_DIR <- paste0(PROJECT_DIR, "RESULTS/scRNAseq/combinedAPOBECresults/filteredVCFs/")
ref_genome <- "BSgenome.Hsapiens.UCSC.hg38"

setwd(OUT_DIR)

# Load custom functions
FUNCTIONS_SCRIPT <- list.files(path = "/Users/m.m.kleisman/Projects/git/pmc_kuiper_projects/", 
                               pattern = "createdFunctions_HM-ALL.R",
                               full.names = TRUE, recursive = TRUE)
source(FUNCTIONS_SCRIPT)



#---------------------------- Variant filtering -----------------------------#

# List input files
unfiltered_vcf_path <- list.files(path = PROJECT_DIR,
                                  pattern = VCF_FILE_EXT,
                                  recursive = TRUE,
                                  full.names = TRUE)
unfiltered_tab_path <- list.files(path = PROJECT_DIR,
                                  pattern = TAB_FILE_EXT,
                                  recursive = TRUE,
                                  full.names = TRUE)

# Define sample and library IDs (file name syntax: LIBRARY_SKION_BMID_VCFEXT)
file_ids <- gsub(VCF_FILE_EXT, "", basename(unfiltered_vcf_path))
sample_ids <- sapply(strsplit(file_ids, "_"), `[`, 2)
library_ids <- sapply(strsplit(file_ids, "_"), `[`, 1)

# Create objects to store info
filtered_vcf_entry_list <- rep(list(NA), length(unfiltered_vcf_path))
filtered_vcf_object_list <- rep(list(NA), length(unfiltered_vcf_path))
names(filtered_vcf_entry_list) <- sample_ids
names(filtered_vcf_object_list) <- sample_ids

# Subset VCF for filtered entries
for (i in (1:length(unfiltered_tab_path))){
  DEPTH <- 20
  
  # Read the VCF and tab converted VCF
  sample_vcf <- loadVcfFile(unfiltered_vcf_path[i],
                            ref_genome,
                            "chromosomes")
  sample_tab <- read.delim(unfiltered_tab_path[i])
  GenomeInfoDb::genome(sample_vcf)  <- "hg38"
  
  # Select SNVs with depth >= 20 and ALT >= 3
  alt_idx <- grep(".ALTCOUNT$", colnames(sample_tab))[1]
  depth_idx <- grep("\\.DP$", colnames(sample_tab))[1]
  filtered_tab <- sample_tab[sample_tab$VARIANT_CLASS == "SNV" &
                             sample_tab[,alt_idx] >= 3 &
                             sample_tab[,depth_idx] >= DEPTH , ]
  
  # Decrease depth if filtered file contains less than 10000 variants
  if (nrow(filtered_tab) < 10000){
    DEPTH <- 10
    filtered_tab <- sample_tab[sample_tab$VARIANT_CLASS == "SNV" &
                                 sample_tab[,alt_idx] >= 3 &
                                 sample_tab[,depth_idx] >= DEPTH , ]
  }
  # Parse filtered info
  filtered_entries <- paste0(filtered_tab$CHROM, ":", 
                             filtered_tab$POS, "_", 
                             filtered_tab$REF, "/", 
                             filtered_tab$ALT)
  
  # Subset VCF for filtered entries
  filtered_vcf <- sample_vcf[rownames(sample_vcf) %in% filtered_entries, ]
  writeVcf(obj = filtered_vcf, filename = gsub(VCF_FILE_EXT, 
                                               paste0("_WGS.HyperExomeFiltered.DPgte", DEPTH, ".ALTgte3.filtered.vcf"), 
                                               basename(unfiltered_vcf_path[i])))
  
  # Add entries to object
  filtered_vcf_entry_list[[i]] <- rownames(filtered_vcf)
  filtered_vcf_object_list[[i]] <- filtered_vcf
}

# Define paths to filtered VCF
filtered_vcf_path <- list.files(pattern = "filtered.vcf",
                                recursive = TRUE,
                                full.names = TRUE)

# Subset VCF for unique entries
for (lib_id in unique(library_ids)){
  
  # Define which samples are present in the library and select those samples' entries
  samples_in_lib <- sample_ids[grep(lib_id, library_ids)]
  samples_entries_list <- filtered_vcf_entry_list[samples_in_lib]
  n_samples <- length(samples_in_lib)
  print(samples_in_lib)

  # Loop over each sample within library to determine unique entries
  for (i in (1:n_samples)){
    sample_id <- samples_in_lib[i]
    
    # Read sample info and info of all other samples
    sample_vcf <- filtered_vcf_object_list[[sample_id]]
    sample_entries <- filtered_vcf_entry_list[[sample_id]]
    other_sample_entries <- unlist(sapply((1:n_samples)[c(!1:n_samples %in% i)], function(x) {samples_entries_list[[x]]}))
    
    # Subset VCF for entries unique to the sample
    unique_sample_entries <- setdiff(sample_entries, other_sample_entries)
    filtered_vcf <- sample_vcf[rownames(sample_vcf) %in% unique_sample_entries, ]
    
    # Write to out file
    sample_vcf_ext <- filtered_vcf_path[grep(sample_id, filtered_vcf_path)]
    writeVcf(obj = filtered_vcf, filename = gsub("vcf", 
                                                 "uniqueEntries.vcf", 
                                                 sample_vcf_ext))
    print(paste0("Unique entries for sample ", sample_id, ": ", length(unique_sample_entries)))
  }
}


# Print number of depth-filtered variants
for (sample in sample_ids){
  print(sample)
  print(length(filtered_vcf_entry_list[[sample]]))
}




