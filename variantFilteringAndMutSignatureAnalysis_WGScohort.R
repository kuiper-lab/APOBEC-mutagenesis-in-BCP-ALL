#------------------------------ Data preparations ----------------------------#

set.seed(1234)

# Clear data 
rm(list = ls())

# Load libraries
ref_genome <- "BSgenome.Hsapiens.UCSC.hg38"
library(ref_genome, character.only = TRUE)
library("MutationalPatterns")
library("ggplot2")
library("logging")
library("dplyr")
library("readxl")
library("openxlsx")
library("stringr")
library("gtools")
library("gridExtra")

# Load custom functions
FUNCTIONS_SCRIPT <- list.files(path = "/Users/m.m.kleisman/Projects/git/pmc_kuiper_projects/", 
                               pattern = "createdFunctions_HM-ALL.R",
                               full.names = TRUE, recursive = TRUE)
source(FUNCTIONS_SCRIPT)

# Set global parameters
PROJECT_DIR <- "/Users/m.m.kleisman/Projects/HypermutatedALL_project/"
RESULTS_DIR <- paste0(PROJECT_DIR, "/RESULTS/mutect2/B-ALL-DxCohort/")
SUBTYPE_INFO <- "rnaSeq_diseaseSubtypes.xlsx"
SAMPLE_ID_INFO <- "ALL_Relapse_Coded_Num.xlsx"
SNV_MUT_MAT_RDATA_FILE <- "biobank_snv_mutation_matrix_B-ALLDx_121125.rdata"
FIT_RES_SIGNATURES <- "fit_res_strict_final_17032024.rdata"

setwd(RESULTS_DIR)

# Define extension of unfiltered VCFs of new samples
# Only required if CONDUCT_VARIANT_FILTERING == TRUE
CONDUCT_SNV_FILTERING <- FALSE           
REMOVE_GERMLINE_VARIANTS <- FALSE
WRITE_FILTERED_VCF_TO_OUTFILE <- FALSE  
SNV_EXT <- "-merged-mutect2Calls.passFiltered.snp.vepAnnotated.vcf.gz" 

# Load input files
load(list.files(path = PROJECT_DIR,
                pattern = FIT_RES_SIGNATURES,
                recursive = TRUE,
                full.names = TRUE))
random_id_info <- as.data.frame(read_excel(list.files(path = PROJECT_DIR,
                                                       pattern = SAMPLE_ID_INFO,
                                                       recursive = TRUE,
                                                       full.names = TRUE)))
centromeric_variants <- list.files(path = PROJECT_DIR,
                                   pattern = "cytoBand.hg38.centromeresOnly.txt",
                                   full.names = TRUE,
                                   recursive = TRUE)
subtype_file <- as.data.frame(read_excel(tail(list.files(path = PROJECT_DIR,
                                                    pattern = SUBTYPE_INFO,
                                                    recursive = TRUE,
                                                    full.names = TRUE),n=1)))


#------------------------------ Variant filtering ----------------------------#

# Conduct SNV filtering  
if (CONDUCT_SNV_FILTERING){
  
  # List input SNV files
  snv_files <- list.files(path = RESULTS_DIR,
                          pattern = SNV_EXT,
                          recursive = TRUE,
                          full.names = TRUE)
  
  # Get sample IDs
  sample_ids <- gsub(SNV_EXT, "", basename(snv_files))
  sample_ids <- gsub("P_", "", sample_ids)
  loginfo(paste("found:", length(sample_ids), "SNV VCF files"))
  
  # Define output directory in case filtered VCFs need to be written
  out_dir <- paste0(paste0(strsplit(snv_files[1], "/")[[1]][1:(length(strsplit(snv_files[1], "/")[[1]]) - 2)], 
                           collapse = "/"), "/filteredData/")
  if (!dir.exists(out_dir)){dir.create(out_dir)}
  
  # Create object to store filtered VCFs
  filtered_vcf_object_list <- rep(list(NA), length(snv_files))
  names(filtered_vcf_object_list) <- sample_ids

  # Filter each VCF file
  for (i in 1:length(snv_files)){
    loginfo(paste("Processing VCF of sample:", sample_ids[i]))
    
    # Read the VCF file
    vcf_file <- loadVcfFile(vcfFilePath = snv_files[i],
                            bsgGenomeName = ref_genome,
                            chromosomeSet = "chromosomes")
    GenomeInfoDb::genome(vcf_file) <- "hg38"
    
    # Get tumor ID from VCF file
    tumor_sample_id_vcf_line <- grep("##tumor_sample=", 
                                     readLines(snv_files[i]), 
                                     value = TRUE)
    tumor_sample_id <- gsub("##tumor_sample=", "", tumor_sample_id_vcf_line)
    loginfo(paste("Found (tumor) sample ID:", tumor_sample_id))
    
    # Exclude variants in the centromeric region
    vcf_no_centromeric_variants <- excludeVariantsInCentromericRegions(vcf_file,
                                                                       centromeric_variants)
    
    # Filter variants on population frequency
    vcf_gnomad_filtered <- filterOnGnomadAlleleFrequency(vcf_no_centromeric_variants,
                                                         0.01)
    
    vcf_gonl_filtered <- filterOnGonlAlleleFrequency(vcf_gnomad_filtered,
                                                     0.01)
    
    # Convert the population frequency filtered VCF to data frame
    df_gonl_filtered <- convertVcfObjectToDataFrame(vcf_gonl_filtered)
    
    # Remove germline variants if requested
    normal_sample_id_vcf_line <- grep("##normal_sample=",
                                      readLines(snv_files[i]),
                                      value = TRUE)
    normal_sample_id <- gsub("##normal_sample=", "", normal_sample_id_vcf_line)
    if (REMOVE_GERMLINE_VARIANTS){
      df_somatic_filtered <- df_gonl_filtered[df_gonl_filtered[,paste0(normal_sample_id, "_read_count_alt")] == 0, ]
    } else {
      df_somatic_filtered <- df_gonl_filtered
    }
    
    # Get column idx for tumor AF, ref read count and alt read count
    AF_idx <- which(names(df_somatic_filtered) %in% paste0(tumor_sample_id, "_read_based_AF"))
    read_count_ref_idx <- which(names(df_somatic_filtered) %in% 
                                paste0(tumor_sample_id, "_read_count_ref"))
    read_count_alt_idx <- which(names(df_somatic_filtered) %in% 
                                paste0(tumor_sample_id, "_read_count_alt"))
    
    # Compute the total depth 
    df_somatic_filtered[paste0(tumor_sample_id, "_total_read_count")] <- 
      df_somatic_filtered[read_count_ref_idx] + df_somatic_filtered[read_count_alt_idx]
    
    # Filter variants on read depth
    total_read_count_idx <- ncol(df_somatic_filtered) 
    df_read_based_filtered <- df_somatic_filtered %>% filter(.data[[c(names(df_somatic_filtered)[read_count_alt_idx])]] >= 3 &
                                       .data[[c(names(df_somatic_filtered)[total_read_count_idx])]] >= 20 &
                                       .data[[c(names(df_somatic_filtered)[AF_idx])]] >= 0.25)
    
    # Subset the VCF-file for popfreq and read depth filtered variants
    filtered_vcf_object <- vcf_file[rownames(vcf_file) %in% 
                                    rownames(df_read_based_filtered), ]
    
    # Store VCF object and write to out file
    filtered_vcf_object_list[[i]] <- filtered_vcf_object
    if (WRITE_FILTERED_VCF_TO_OUTFILE){
      out_ext <- gsub(".vcf.gz", ".popMaxandReadDepthFiltered.vcf", SNV_EXT)
      writeVcf(filtered_vcf_object, file = paste0(out_dir, "/", sample_ids[i], "-", tumor_sample_id, "-", 
                                                  normal_sample_id, out_ext))
    }
  } 
  
  # Convert VCF objects to a named GRangesList
  data_grl <- sapply(filtered_vcf_object_list, granges)
  names(data_grl) <- names(filtered_vcf_object_list)
  
  # Create a mutation matrix
  sbs_mutation_matrix_combined <- mut_matrix(data_grl, ref_genome = ref_genome)
  
  # Write mutation matrix to RData file
  save(sbs_mutation_matrix_combined, file = SNV_MUT_MAT_RDATA_FILE)
  
  # Remove filtered VCF files
  rm(filtered_vcf_object_list)
} else {
  print("Variant filtering not requested, mutation matrix from RData file will be used for further analyses.")
  load(list.files(path = PROJECT_DIR,
                  pattern = SNV_MUT_MAT_RDATA_FILE,
                  recursive = TRUE,
                  full.names = TRUE))
}

# Rename mutation matrix from RData file for convenience
snv_mut_mat <- sbs_mutation_matrix_combined
colnames(snv_mut_mat) <- random_id_info$Coded_num[match(colnames(snv_mut_mat), random_id_info$SKION)]
rm(sbs_mutation_matrix_combined)



#-------------------------------- Strict refit --------------------------------#

# Function definitions
# Function to calculate the absolute contribution
calculate_absolute_contribution <- function(contribution_col){
  contributions <- contribution_col[1:(length(contribution_col)-2)]
  total <- contribution_col[length(contribution_col)-1]
  originaltotal <- contribution_col[length(contribution_col)]
  absolute_contributions <- ((contributions / total) * originaltotal)
  return(absolute_contributions)
}

# Function to calculate the relative contribution
## Calculate relative contribution
calculate_relative_contribution <- function(contribution_col){
  contributions <- contribution_col[1:(length(contribution_col) - 1)]
  total <- contribution_col[length(contribution_col)]
  relative_contributions <- contributions / total
  return(relative_contributions)
}

# Extract refitted signatures
refitted_signatures <- fit_res_strict_final$signatures

# Only use diagnosis-related signatures
refitted_signatures <- refitted_signatures[,colnames(refitted_signatures) %in% 
                                             c("SBS1", "SBS2", "SBS13", "SBS7a", "SBS18", "SBSA")]

# Perform a strict refit on the mutation matrix
strict_refit <- fit_to_signatures_strict(snv_mut_mat, 
                                         refitted_signatures, 
                                         max_delta = 0.033)
fit_res_strict <- strict_refit$fit_res

# Calculate contribution totals based on refit
contribution_with_totals <- rbind(fit_res_strict$contribution, 
                                  apply(fit_res_strict$contribution, 2, sum))
rownames(contribution_with_totals)[nrow(contribution_with_totals)] <- "Total"

# Calculate the relative contribution
relative_contributions <- apply(contribution_with_totals, 2, calculate_relative_contribution)

# Calculate contribution totals based on original data
contribution_with_totals <- rbind(contribution_with_totals, 
                                  apply(snv_mut_mat, 2, sum))
rownames(contribution_with_totals)[nrow(contribution_with_totals)] <- "OriginalTotal"

# Calculate the absolute contribution
absolute_contributions <- apply(contribution_with_totals, 2, calculate_absolute_contribution)

# Round values
absolute_contributions <- round(absolute_contributions)

# Create the reconstruction
reconstruction <- diag(cos_sim_matrix(snv_mut_mat, fit_res_strict$reconstructed[]))
absolute_contributions <- rbind(absolute_contributions, reconstruction)
rownames(absolute_contributions)[nrow(absolute_contributions)] <- "Reconstruction"

# Perform bootstrapping
contri_boots <- fit_to_signatures_bootstrapped(snv_mut_mat,
                                               refitted_signatures,
                                               n_boots = 100,
                                               method = "strict",
                                               max_delta = 0.033)

# Calculate presence percentage for bootstrap
absolute_contributions_with_bootstrap <- rbind(absolute_contributions, 
                                               matrix(ncol = ncol(absolute_contributions), 
                                                      nrow = ncol(contri_boots), 
                                                      dimnames = list(paste0(colnames(contri_boots), "_bootstrappercentage"), 
                                                                      colnames(absolute_contributions))))
for (patient in colnames(fit_res_strict$contribution)){
  patient_bootstrap <- contri_boots[str_remove(rownames(contri_boots), "_\\d+$") == patient, ]
  signature_percentage <- colSums(patient_bootstrap > 0)
  absolute_contributions_with_bootstrap[grep("bootstrap", rownames(absolute_contributions_with_bootstrap)), patient] <- signature_percentage
}

# Define samples order based on abs contribution of SBS2/SBS13
sample_order <- sort(colSums(absolute_contributions[c("SBS2", "SBS13"),]), decreasing = TRUE)
fit_res_strict$contribution <- fit_res_strict$contribution[, names(sample_order)]
pos_samples_ordered <- names(sample_order)[sample_order > 0]

# Define plot colors
CUSTOM_COLOUR_PALETTE <- c(SBS1 = "#9BC0CD",
                           SBS2 = "#CE2627",
                           SBS5 = "#A9A9A9",
                           SBS7a = "#F6EB13",
                           SBS13 = "#B12325",
                           SBS18 = "#AA6AAC",
                           SBSA = "#DADADA",
                           SBS86 = "#0A0F25",
                           SBS87 = "#5C6BC0",
                           SBS17a = "#FFC2D1",
                           SBS17b = "#FF8FAB")
plot_cols <- CUSTOM_COLOUR_PALETTE[names(CUSTOM_COLOUR_PALETTE) %in% rownames(fit_res_strict$contribution)]

# Define signature order
sign_order <- c(c("SBS2", "SBS13"), setdiff(names(plot_cols), c("SBS2", "SBS13")))
sign_order <- rev(sign_order)

# Generate absolute contribution plot of SBS2/13 pos samples
abs_contr_plot <- plot_contribution(absolute_contributions[sign_order, pos_samples_ordered],
                                    coord_flip = FALSE,
                                    mode = "absolute",
                                    palette = plot_cols) +
  ylab("Absolute contribution") +
  scale_fill_manual(values = plot_cols, breaks = names(plot_cols)) +
  theme(axis.text.x = element_text(angle = -90, size = 10,
                                   vjust = 0.5, hjust = 0))

# Generate relative contribution plot of SBS2/13 pos samples
rel_contr_plot <- plot_contribution(relative_contributions[sign_order, pos_samples_ordered],
                                    coord_flip = FALSE,
                                    mode = "relative",
                                    palette = plot_cols) +
  scale_fill_manual(values = plot_cols, breaks = names(plot_cols)) +
  theme(axis.text.x = element_text(angle = -90, size = 10,
                                   vjust = 0.5, hjust = 0))

# Plot original vs reconstructed cosine similarity
reconstr_cos_sim <- plot_original_vs_reconstructed(snv_mut_mat[,pos_samples_ordered], 
                                                   fit_res_strict$reconstructed[,pos_samples_ordered], 
                                                   y_intercept = 0.9,
                                                   ylims = c(0, 1)) +
  ylab("Cos. sim. (orig. vs reconstr.)") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  theme(axis.text.x = element_text(angle = -90, size = 10,
                                   vjust = 0.5, hjust = 0))

# Save to PDF
layout_mat <- matrix(c(rep(1, 33), 
                       rep(2, 33), 
                       rep(3, 10), NA, 
                       rep(3, 10), NA), 
                     ncol = 11, byrow = TRUE)
pdf("ContributionPlots_completeWGScohort-withSBS213mutations.pdf", width = 12, height = 9)
grid.arrange(abs_contr_plot, 
             rel_contr_plot,
             reconstr_cos_sim,
             layout_matrix = layout_mat)
dev.off()

# Write absolute contributions and bootstrap values to out file
overview_df <- cbind(data.frame(Sample = names(sample_order)),
                     Skion = random_id_info$SKION[match(names(sample_order), random_id_info$Coded_num)],
                     t(absolute_contributions_with_bootstrap)[names(sample_order), ])
overview_df["Subtype"] <- subtype_file$Subtype[match(overview_df$Skion, subtype_file$SKION_ID)]
write.xlsx(overview_df, file = "abscontr-bootstrap-reconsruction-completeWGScohort.xlsx")



