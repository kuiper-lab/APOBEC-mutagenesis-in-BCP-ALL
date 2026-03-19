#!/usr/bin/env Rscript 
#
#------------------------------ Data preparations ----------------------------#

# Clear data 
rm(list = ls())
set.seed(1234)

# Load libraries
ref_genome <- "BSgenome.Hsapiens.UCSC.hg38"
library(ref_genome, character.only = TRUE)
library("MutationalPatterns")
library("ggplot2")
library("logging")
library("dplyr")
library("plyr")
library("reshape2")
library("readxl")
library("openxlsx")
library("stringr")
library("gtools")
library("gridExtra")
library("ggpubr")

# Set global parameters
PROJECT_DIR <- "/Users/m.m.kleisman/Projects/HypermutatedALL_project/"
RESULTS_DIR <- paste0(PROJECT_DIR, "/RESULTS/mutect2/APOBEC/")
RANDOM_ID_INFO <- "ALL_ETV6RUNX1_Coded_Num.xlsx"
FIT_RES_SIGNATURES <- "fit_res_strict_final_17032024.rdata"  # for mut mat of refitted signatures
RNASEQ_MUT_MAT <- paste0("/Users/m.m.kleisman/Projects/leukemiaRnaSeq/RESULTS/DE_analysis/",
                         "count_data_ALLsamples_PMC_duplicatesAndIsoformsMerged.Rda")
setwd(RESULTS_DIR)

# Load custom functions
FUNCTIONS_SCRIPT <- list.files(path = "/Users/m.m.kleisman/Projects/git/pmc_kuiper_projects/", 
                               pattern = "createdFunctions_HM-ALL.R",
                               full.names = TRUE, recursive = TRUE)
source(FUNCTIONS_SCRIPT)

# Define extension of unfiltered VCFs 
#! Only required if CONDUCT_SNV_FILTERING == TRUE
CONDUCT_SNV_FILTERING <- FALSE 
REMOVE_GERMLINE_VARIANTS <- FALSE
SNV_EXT <- "-merged-mutect2Calls.passFiltered.snp.vepAnnotated.vcf.gz" 

# Define name of the mutation matrix and GRL
#! Files will be created if CONDUCT_SNV_FILTERING == TRUE, otherwise they are loaded
#ETV6RUNX1_MUT_MAT <- "ETV6-RUNX1_mut_mat_170624.Rda"
#ETV6RUNX1_GRANGES_LIST <- "ETV6-RUNX1_granges_list_170624.Rda"
ETV6RUNX1_MUT_MAT <- "ETV6-RUNX1_mut_mat_121125.Rda"
ETV6RUNX1_GRANGES_LIST <- "ETV6-RUNX1_granges_list_121125.Rda"

# Define mut mat of variants between VAF 0.25 and 0.05
#! Obtained by variant filtering at VAF == 0.05 with CONDUCT_SNV_FILTERING == TRUE
ETV6RUNX1_MUT_MAT_CUSTOM_VAF <- "ETV6-RUNX1_mut_mat_0.05-0.25VAFfilter_121125.Rda"
#ETV6RUNX1_MUT_MAT_CUSTOM_VAF <- "ETV6-RUNX1_mut_mat_0.05-0.25VAFfilter_290924.Rda"
#ETV6RUNX1_MUT_MAT_CUSTOM_VAF <- "ETV6-RUNX1_mut_mat_0.05VAFfilter_280624.Rda"

# Define VAF threshold
VAF <- 0.05

# Load input files
load(RNASEQ_MUT_MAT)
load(list.files(path = PROJECT_DIR,
                pattern = FIT_RES_SIGNATURES,
                recursive = TRUE,
                full.names = TRUE))
centromeric_variants <- list.files(path = PROJECT_DIR,
                                   pattern = "cytoBand.hg38.centromeresOnly.txt",
                                   full.names = TRUE,
                                   recursive = TRUE)
random_id_file <- as.data.frame(read_excel(tail(list.files(path = PROJECT_DIR,
                                                           pattern = RANDOM_ID_INFO,
                                                           recursive = TRUE,
                                                           full.names = TRUE), n = 1)))

# Determine ETV6::RUNX1 positive samples
etv6_runx1_samples <- random_id_file[random_id_file$Subtype == "ETV6::RUNX1", "Sample"]



#------------------------------ Variant filtering ----------------------------#

## Visualize frequency of VAFs and depths of unfiltered variants
if (CONDUCT_SNV_FILTERING){
  
  # List input SNV files
  snv_files <- list.files(path = RESULTS_DIR,
                          pattern = SNV_EXT,
                          recursive = TRUE,
                          full.names = TRUE)
  
  # Get sample IDs
  skion_ids <- gsub(SNV_EXT, "", basename(snv_files))
  
  # Only keep files of ETV6::RUNX1 samples
  idx_keep <- which(skion_ids %in% etv6_runx1_samples)
  skion_ids <- skion_ids[idx_keep]
  snv_files <- snv_files[idx_keep]
  loginfo(paste("found:", length(skion_ids), "SNV VCF files"))
  
  # Initiate empty objects
  vaf_freq_df <- data.frame()
  depth_freq_df <- data.frame()
  
  # Read each VCF file
  for (i in 1:length(snv_files)){
    vcf_file <- loadVcfFile(vcfFilePath = snv_files[i],
                            bsgGenomeName = ref_genome,
                            chromosomeSet = "chromosomes")
    GenomeInfoDb::genome(vcf_file) <- "hg38"
    
    # Get tumor ID from VCF file
    tumor_sample_id_vcf_line <- grep("##tumor_sample=", 
                                     readLines(snv_files[i]), 
                                     value = TRUE)
    tumor_sample_id <- gsub("##tumor_sample=", "", tumor_sample_id_vcf_line)
    
    # Convert VCF to data frame
    vcf_file <- convertVcfObjectToDataFrame(vcf_file)
    
    # Get variant VAFs and round values to 2 decimals
    sample_VAFs <- vcf_file[, colnames(vcf_file) == paste0(tumor_sample_id, "_read_based_AF")]
    sample_VAFs <- round_any(sample_VAFs, 
                             accuracy = 0.01, 
                             f = ceiling)
    
    # Add to object
    sample_VAFs_freq <- data.frame(VAF = as.numeric(names(table(sample_VAFs))),
                                   Frequency = as.numeric(table(sample_VAFs)),
                                   Sample = skion_ids[i])
    vaf_freq_df <- rbind(vaf_freq_df, sample_VAFs_freq[sample_VAFs_freq$VAF != 0, ])
    
    # Get variant depths over an interval of 10 reads
    sample_depth <- as.numeric(rowSums(vcf_file[, grep(paste0(tumor_sample_id, '_read_count'), 
                                                       colnames(vcf_file))]))
    sample_depth <- round_any(sample_depth, 
                              accuracy = 10, 
                              f = ceiling)
    
    # Add to object
    sample_depth_freq <- data.frame(Depth = as.numeric(names(table(sample_depth))),
                                    Frequency = as.numeric(table(sample_depth)),
                                    Sample = skion_ids[i])
    depth_freq_df <- rbind(depth_freq_df, sample_depth_freq)
  }
  
  # Generate depth plot
  pdf("ETV6-RUNX1_VAFandDepthOverviews.pdf", width = 30, height = 25)
  print(ggplot(depth_freq_df, aes(x = Depth, y = Frequency)) +
    facet_wrap(Sample ~ ., scales = "free", ncol = 3) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_x_continuous(limits = c(0,250), breaks = seq(0,250,10), expand = c(0.001, 0.001)) +
    theme_bw() + labs(title = "Depth frequency") +
    theme(panel.grid.minor.x = element_blank(),
          panel.grid.major.x = element_line(color = "grey95"),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          strip.text = element_text(angle = 0, size = 13),
          axis.text.x = element_text(angle = -90)))
  
  # Generate VAF plots
  for (title in c("VAF frequency", "VAF frequency (log-transformed)")){
    ifelse(title == "VAF frequency", y_value <- vaf_freq_df$Frequency, 
           y_value <- log10(vaf_freq_df$Frequency))
    print(ggplot(vaf_freq_df, aes(x = VAF, y = y_value)) +
            facet_wrap(Sample ~ ., scales = "free", ncol = 3) +
            geom_bar(stat = "identity", position = "dodge") +
            scale_x_continuous(breaks = seq(0.01,1,0.01), expand = c(0.001, 0.001)) +
            theme_bw() + labs(title = title, y = "Frequency") +
            theme(panel.grid.minor.x = element_blank(),
                  panel.grid.major.x = element_line(color = "grey95"),
                  panel.grid.major.y = element_blank(),
                  panel.grid.minor.y = element_blank(),
                  plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
                  strip.text = element_text(angle = 0, size = 13),
                  axis.text.x = element_text(angle = -90)))
  }
  dev.off()
}


## Conduct variant filtering
if (CONDUCT_SNV_FILTERING){
  
  # List input SNV files
  snv_files <- list.files(path = RESULTS_DIR,
                          pattern = SNV_EXT,
                          recursive = TRUE,
                          full.names = TRUE)
  
  # Get sample IDs
  skion_ids <- gsub(SNV_EXT, "", basename(snv_files))
  
  # Only keep files of ETV6::RUNX1 samples
  idx_keep <- which(skion_ids %in% etv6_runx1_samples)
  skion_ids <- skion_ids[idx_keep]
  snv_files <- snv_files[idx_keep]
  loginfo(paste("found:", length(skion_ids), "SNV VCF files"))
  
  # Create object to store filtered VCFs
  filtered_vcf_object_list <- rep(list(NA), length(snv_files))
  names(filtered_vcf_object_list) <- skion_ids

  # Filter each VCF file
  for (i in 1:length(snv_files)){
    loginfo(paste("Processing VCF of sample:", skion_ids[i]))
    
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
    vcf_no_centromeric_variants <- excludeVariantsInCentromericRegionsUpdated(vcf_file,
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
    
    # Get column idx for tumor AF, ref count and alt count 
    AF_idx <- which(names(df_somatic_filtered) %in% 
                      paste0(tumor_sample_id, "_read_based_AF"))
    read_count_ref_idx <- which(names(df_somatic_filtered) %in% 
                                  paste0(tumor_sample_id, "_read_count_ref"))
    read_count_alt_idx <- which(names(df_somatic_filtered) %in% 
                                  paste0(tumor_sample_id, "_read_count_alt"))
    
    # Compute the total depth 
    df_somatic_filtered[paste0(tumor_sample_id, "_total_read_count")] <- 
      df_somatic_filtered[read_count_ref_idx] + 
      df_somatic_filtered[read_count_alt_idx]
    
    # Get column idx for total read count
    total_read_count_idx <- ncol(df_somatic_filtered) 
    
    # Filter variants on read depth
    df_read_based_filtered <- 
      df_somatic_filtered %>% filter(.data[[c(names(df_somatic_filtered)[read_count_alt_idx])]] >= 3 &
                                       .data[[c(names(df_somatic_filtered)[total_read_count_idx])]] >= 20 &
                                       .data[[c(names(df_somatic_filtered)[AF_idx])]] >= VAF)
    if (VAF < 0.25){
      df_read_based_filtered <- df_read_based_filtered %>% filter(.data[[c(names(df_read_based_filtered)[AF_idx])]] < 0.25)
    }

    # Store filtered VCF object 
    filtered_vcf_object_list[[i]] <- vcf_file[rownames(vcf_file) %in% 
                                              rownames(df_read_based_filtered), ]
  } 

  # Convert VCF objects to a named GRangesList
  sbs_grl_etv6runx1 <- sapply(filtered_vcf_object_list, granges)
  names(sbs_grl_etv6runx1) <- names(filtered_vcf_object_list)
  
  # Create a mutation matrix
  sbs_mutation_matrix_etv6runx1 <- mut_matrix(sbs_grl_etv6runx1, 
                                              ref_genome = ref_genome)
  # Remove tmp files
  rm(filtered_vcf_object_list)
  
  # Write grl and mutation matrix to RData file
  if (VAF == 0.25){
    save(sbs_grl_etv6runx1, file = ETV6RUNX1_GRANGES_LIST)
    save(sbs_mutation_matrix_etv6runx1, file = ETV6RUNX1_MUT_MAT)
  } else {
    save(sbs_mutation_matrix_etv6runx1, file = ETV6RUNX1_MUT_MAT_CUSTOM_VAF)
  }
  
} else {
  print("Variant filtering not requested, mutation matrix from RData file will be used for further analyses.")
  load(list.files(path = PROJECT_DIR,
                  pattern = ETV6RUNX1_GRANGES_LIST,
                  full.names = TRUE,
                  recursive = TRUE))
  load(list.files(path = PROJECT_DIR,
                  pattern = ETV6RUNX1_MUT_MAT,
                  full.names = TRUE,
                  recursive = TRUE))
  skion_ids <- colnames(sbs_mutation_matrix_etv6runx1)
}

# Rename mutation matrix from RData file for convenience
snv_mut_mat <- sbs_mutation_matrix_etv6runx1[,etv6_runx1_samples]
rm(sbs_mutation_matrix_etv6runx1)
subtype_pos_ids <- colnames(snv_mut_mat)

# Load mutation matrix at VAF 0.05-0.25
load(ETV6RUNX1_MUT_MAT_CUSTOM_VAF)
snv_mut_mat_lowVAF <- sbs_mutation_matrix_etv6runx1
snv_mut_mat_lowVAF <- snv_mut_mat_lowVAF[, subtype_pos_ids]
rm(sbs_mutation_matrix_etv6runx1)

# Randomise IDs
colnames(snv_mut_mat) <- random_id_file$Anonym[match(colnames(snv_mut_mat), random_id_file$Sample)]
colnames(snv_mut_mat_lowVAF) <- random_id_file$Anonym[match(colnames(snv_mut_mat_lowVAF), random_id_file$Sample)]
names(subtype_pos_ids) <- colnames(snv_mut_mat)



#----------------------- Signature contribution analysis ----------------------#

## Function definitions
# Function to calculate the absolute contribution
calculate_absolute_contribution <- function(contribution_col){
  contributions <- contribution_col[1:(length(contribution_col)-2)]
  total <- contribution_col[length(contribution_col)-1]
  originaltotal <- contribution_col[length(contribution_col)]
  absolute_contributions <- ((contributions / total) * originaltotal)
  return(absolute_contributions)
}

# Function to calculate the relative contribution
calculate_relative_contribution <- function(contribution_col){
  contributions <- contribution_col[1:(length(contribution_col) - 1)]
  total <- contribution_col[length(contribution_col)]
  relative_contributions <- contributions / total
  return(relative_contributions)
}

# Wrapper function to calculate signature contributions and perform bootstrapping
wrapper_calculate_contributions <- function(target_mut_mat,
                                            target_strict_refit){
  # Calculate burden based on refit
  contribution_with_totals <- rbind(target_strict_refit$contribution, 
                                    apply(target_strict_refit$contribution, 2, sum))
  rownames(contribution_with_totals)[nrow(contribution_with_totals)] <- "Total"
  
  # Calculate the relative contribution
  relative_contributions <- apply(contribution_with_totals, 2, calculate_relative_contribution)
  
  # Calculate burden based on original data
  contribution_with_totals <- rbind(contribution_with_totals, 
                                    apply(target_mut_mat, 2, sum))
  rownames(contribution_with_totals)[nrow(contribution_with_totals)] <- "OriginalTotal"
  
  # Calculate the absolute contribution and round values
  absolute_contributions <- apply(contribution_with_totals, 2, calculate_absolute_contribution)
  absolute_contributions <- round(absolute_contributions)
  
  # Create the reconstruction
  reconstruction <- diag(cos_sim_matrix(target_mut_mat, target_strict_refit$reconstructed[]))
  absolute_contributions <- rbind(absolute_contributions, reconstruction)
  rownames(absolute_contributions)[nrow(absolute_contributions)] <- "Reconstruction"
  
  return(list(absolute_contributions, relative_contributions))
}

# Wrapper function to write contributions and bootstrapping values to Excel
wrapper_complete_overview <- function(target_mut_mat,
                                      target_strict_refit,
                                      abs_contribution,
                                      rel_contribution,
                                      refitted_signatures,
                                      out_file = "",
                                      write_to_out_file = TRUE){
  # Perform bootstrapping
  contri_boots <- fit_to_signatures_bootstrapped(target_mut_mat,
                                                 refitted_signatures,
                                                 n_boots = 100,
                                                 method = "strict",
                                                 max_delta = 0.033)
  # Calculate presence percentage for bootstrap
  absolute_contributions_with_bootstrap <- rbind(abs_contribution, 
          matrix(ncol = ncol(abs_contribution), 
                 nrow = ncol(contri_boots), 
                 dimnames = list(paste0(colnames(contri_boots), "_bootstrapPercentage"), colnames(abs_contribution))))
  for (patient in colnames(target_strict_refit$contribution)){
    patient_bootstrap <- contri_boots[str_remove(rownames(contri_boots), "_\\d+$") == patient, ]
    signature_percentage <- colSums(patient_bootstrap > 0)
    absolute_contributions_with_bootstrap[grep("bootstrap", rownames(absolute_contributions_with_bootstrap)), patient] <- signature_percentage
  }
  
  # Save info to Excel
  contr_info_df <- rbind(rel_contribution, absolute_contributions_with_bootstrap)
  rownames(contr_info_df)[1:(ncol(refitted_signatures)*2)] <- paste0(rownames(contr_info_df)[1:(ncol(refitted_signatures)*2)],
                                                                     rep(c("_relativeContr", "_absoluteContr"), each = ncol(refitted_signatures)))
  contr_info_df <- cbind(data.frame(SKION = subtype_pos_ids[match(colnames(contr_info_df), names(subtype_pos_ids))],
                                    Sample = colnames(contr_info_df)), 
                         t(contr_info_df))
  if (write_to_out_file){
    write.xlsx(contr_info_df, out_file)
  }
  return(as.data.frame(t(absolute_contributions_with_bootstrap)))
}

# Function to create an overview figure of contributions and reconstruction values
create_contribution_overview_figure <- function(target_mut_mat,
                                                reconstructed_strict_refit,
                                                abs_contribution,
                                                rel_contribution,
                                                signature_colors){
  
  # Define plot colors and signature order
  plot_cols <- signature_colors[names(signature_colors) %in% 
                                rownames(fit_res_strict$contribution)]
  sign_order <- c(c("SBS2", "SBS13"), setdiff(names(plot_cols), c("SBS2", "SBS13")))
  sign_order <- rev(sign_order)
  
  # Generate absolute contribution plot of ETV6::RUNX1 samples
  abs_contr_plot <- plot_contribution(abs_contribution[sign_order,],
                                      coord_flip = FALSE,
                                      mode = "absolute",
                                      palette = plot_cols) +
    scale_fill_manual(values = plot_cols, breaks = names(plot_cols)) +
    ylab("Absolute contribution") +
    theme(axis.text.x = element_text(angle = 90, size = 10,
                                     vjust = 0.5, hjust = 0))
  
  # Generate relative contribution plot of ETV6::RUNX1 samples
  rel_contr_plot <- plot_contribution(rel_contribution[sign_order,],
                                      coord_flip = FALSE,
                                      mode = "relative",
                                      palette = plot_cols) +
    scale_fill_manual(values = plot_cols, breaks = names(plot_cols)) +
    theme(axis.text.x = element_text(angle = 90, size = 10,
                                     vjust = 0.5, hjust = 0))
  
  # Plot original vs reconstructed cosine similarity
  reconstr_cos_sim <- plot_original_vs_reconstructed(target_mut_mat, 
                                                     reconstructed_strict_refit, 
                                                     y_intercept = 0.9,
                                                     ylims = c(0, 1)) +
    ylab("Cos. sim. (orig. vs reconstr.)") +
    scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
    theme(axis.text.x = element_text(angle = 90, size = 10,
                                     vjust = 0.5, hjust = 0))
  
  # Create overview figure
  overview_figures <- list(abs_contr_plot, rel_contr_plot, reconstr_cos_sim)
  return(overview_figures)
}


## Perform strict refit
# Extract refitted signatures
refitted_signatures <- fit_res_strict_final$signatures

# Only use diagnosis-related signatures
refitted_signatures <- refitted_signatures[,colnames(refitted_signatures) %in% 
                                             c("SBS1", "SBS2", "SBS13", "SBS7a", "SBS18", "SBSA")]

# Perform a strict refit on the mutation matrices
strict_refit <- fit_to_signatures_strict(snv_mut_mat, 
                                         refitted_signatures, 
                                         max_delta = 0.033)
fit_res_strict <- strict_refit$fit_res
strict_refit <- fit_to_signatures_strict(snv_mut_mat_lowVAF, 
                                         refitted_signatures, 
                                         max_delta = 0.033)
fit_res_strict_lowVAF <- strict_refit$fit_res


## Calculate signature contributions
snv_contr <- wrapper_calculate_contributions(target_mut_mat = snv_mut_mat,
                                             target_strict_refit = fit_res_strict)
abs_contr <- snv_contr[[1]]; rel_contr <- snv_contr[[2]]
snv_contr_lowVAF <- wrapper_calculate_contributions(target_mut_mat = snv_mut_mat_lowVAF,
                                                    target_strict_refit = fit_res_strict_lowVAF)
abs_contr_lowVAF <- snv_contr_lowVAF[[1]]; rel_contr_lowVAF <- snv_contr_lowVAF[[2]]

# Create boxplot of mutational load
mut_load_df <- data.frame(Load = colSums(abs_contr[-nrow(abs_contr),]))
mut_load_df["Type"] <- ifelse(rownames(mut_load_df) %in% random_id_file$Anonym[random_id_file$APOBEC_group == "Clonal"], "Yes", "No")
mut_load_df$Type <- gsub("No", paste("No\nn =", as.numeric(table(mut_load_df$Type)[1])), mut_load_df$Type)
mut_load_df$Type <- gsub("Yes", paste("Yes\nn =", as.numeric(table(mut_load_df$Type)[2])), mut_load_df$Type)
mut_load_df$Type <- factor(mut_load_df$Type, levels = rev(unique(sort(mut_load_df$Type))))
mut_load_plot <- ggplot(mut_load_df, aes(x = Type, y = Load, fill = Type)) +
  geom_boxplot(width = 0.8) +
  geom_signif(comparisons = list(levels(mut_load_df$Type)), map_signif_level = FALSE) +
  scale_y_continuous(limits = c(-50, 6800), breaks  = seq(0,6000, 1500), expand = c(0, 0.1)) +
  scale_fill_manual(values = c("#CE2627", "grey70")) +
  labs(x = "SBS2/13 in main clone", y = "Total mutational load") +
  theme_bw(base_size = 14) +
  theme(legend.position = "none", axis.text = element_text(size = 12.5),
        panel.grid = element_blank())

# Save to PDF
pdf("MutLoadOverview_121125.pdf", width = 4.5, height = 4.5)
mut_load_plot
dev.off()

# Perform bootstrapping and write complete contribution overview to Excel
bootstr_values <- wrapper_complete_overview(target_mut_mat = snv_mut_mat,
                                            target_strict_refit = fit_res_strict,
                                            abs_contribution = abs_contr,
                                            rel_contribution = rel_contr,
                                            refitted_signatures = refitted_signatures,
                                            out_file = "signatureContributions-Bootstrap-Reconstruction-ETV6RUNX1samples_121125.xlsx",
                                            write_to_out_file = FALSE)
bootstr_values_lowVAF <- wrapper_complete_overview(target_mut_mat = snv_mut_mat_lowVAF,
                                                   target_strict_refit = fit_res_strict_lowVAF,
                                                   abs_contribution = abs_contr_lowVAF,
                                                   rel_contribution = rel_contr_lowVAF,
                                                   refitted_signatures = refitted_signatures,
                                                   out_file = "signatureContributions-Bootstrap-Reconstruction-ETV6RUNX1samples-VAF0.05-0.25_121125.xlsx",
                                                   write_to_out_file = FALSE)

# Determine which samples have clonal and subclonal mutations
clonal_samples <- rowSums(bootstr_values[,colnames(refitted_signatures)]) >= 200 & 
  (bootstr_values$SBS2_bootstrapPercentage == 100 | bootstr_values$SBS13_bootstrapPercentage == 100)
clonal_samples <- rownames(bootstr_values)[which(clonal_samples)]
subclonal_samples <- rowSums(bootstr_values_lowVAF[,colnames(refitted_signatures)]) >= 200 & 
  (bootstr_values_lowVAF$SBS2_bootstrapPercentage == 100 | bootstr_values_lowVAF$SBS13_bootstrapPercentage == 100)
subclonal_samples <- rownames(bootstr_values_lowVAF)[which(subclonal_samples)]


## Generate overview figures
CUSTOM_COLOUR_PALETTE <- c(SBS1 = "#9BC0CD",
                           SBS5 = "#A9A9A9",
                           SBS7a = "#F6EB13",
                           SBS2 = "#CE2627",
                           SBS13 = "#B12325",
                           SBS18 = "#AA6AAC",
                           SBSA = "#DADADA",
                           SBS86 = "#0A0F25",
                           SBS87 = "#5C6BC0",
                           SBS17a = "#FFC2D1",
                           SBS17b = "#FF8FAB")
hide_axis_labels <- theme(plot.title = element_text(face = "bold", size = 12),
                          axis.text.x = element_blank(), 
                          axis.ticks.x = element_blank(),
                          legend.text = element_text(size = 10))

# Define samples order based on abs contribution of clonal and subclonal SBS2/SBS13 mutations
all_ordered <- names(sort(colSums(abs_contr), decreasing = TRUE))
clonal_ordered <- sort(colSums(abs_contr[c("SBS2", "SBS13"),]), decreasing = TRUE)
subclonal_ordered <- sort(colSums(abs_contr_lowVAF[c("SBS2", "SBS13"),]), decreasing = TRUE)
clonal_ordered <- c(names(clonal_ordered)[names(clonal_ordered) %in% clonal_samples],
                    setdiff(names(clonal_ordered)[clonal_ordered > 0], clonal_samples)) 
#full_subclonal_ordered <- c(names(subclonal_ordered)[names(subclonal_ordered) %in% subclonal_samples],
#                       setdiff(names(subclonal_ordered)[subclonal_ordered > 0], subclonal_samples)) 
full_subclonal_ordered <- names(subclonal_ordered)[subclonal_ordered > 0] 
subclonal_ordered <- setdiff(full_subclonal_ordered, clonal_ordered)
none_ordered <- setdiff(all_ordered, c(clonal_ordered, subclonal_ordered)) 
sample_order <- c(clonal_ordered, subclonal_ordered, none_ordered)
#cat(sample_order, sep = "\n")

# Create overview figures
snv_overview <- create_contribution_overview_figure(target_mut_mat = snv_mut_mat[, sample_order],
                                                    reconstructed_strict_refit = fit_res_strict$reconstructed[, sample_order],
                                                    abs_contribution = abs_contr[, sample_order], 
                                                    rel_contribution = rel_contr[, sample_order],
                                                    signature_colors = CUSTOM_COLOUR_PALETTE)
snv_overview_lowVAF <- create_contribution_overview_figure(target_mut_mat = snv_mut_mat_lowVAF[, sample_order],
                                                           reconstructed_strict_refit = fit_res_strict_lowVAF$reconstructed[, sample_order],
                                                           abs_contribution = abs_contr_lowVAF[c("SBS2", "SBS13"), sample_order], 
                                                           rel_contribution = rel_contr_lowVAF[c("SBS2", "SBS13"), sample_order],
                                                           signature_colors = CUSTOM_COLOUR_PALETTE[c("SBS2", "SBS13")])

# Plot overview figures
layout_mat <- matrix(c(rep(1, 22), rep(2, 22), rep(3, 10), NA, rep(3, 10), NA), 
                     ncol = 11, byrow = TRUE)
snv_plots_arranged <- grid.arrange(snv_overview[[1]] + hide_axis_labels, 
                                   snv_overview[[2]] + hide_axis_labels,
                                   snv_overview[[3]],
                                   layout_matrix = layout_mat)
snv_plots_arranged_lowVAF <- grid.arrange(snv_overview_lowVAF[[1]] + hide_axis_labels, 
                                          snv_overview_lowVAF[[2]] + hide_axis_labels,
                                          snv_overview_lowVAF[[3]],
                                          layout_matrix = layout_mat)

# Save to PDF
pdf("ContributionPlots_allETV6RUNX1samples_ordered_121125.pdf", width = 12, height = 9)
plot(snv_plots_arranged)
plot(snv_plots_arranged_lowVAF)
dev.off()

# Create annotation data frame
annot_df <- data.frame(Sample = sample_order,
                       Clonal = ifelse(sample_order %in% clonal_samples, "Passed", "N/A"),
                       Subclonal = ifelse(sample_order %in% subclonal_samples, "Passed", "N/A"))
annot_df$Clonal <- ifelse(sample_order %in% setdiff(clonal_ordered, clonal_samples), "Bootstrap <100", annot_df$Clonal)
annot_df$Subclonal <- ifelse(sample_order %in% setdiff(full_subclonal_ordered, subclonal_samples), "Bootstrap <100", annot_df$Subclonal)


# Parse figures to that clonal and subclonal samples are highlighted
plot_cols <- CUSTOM_COLOUR_PALETTE[names(CUSTOM_COLOUR_PALETTE) %in% c("SBS1", "SBS2", "SBS13", "SBS18", "SBSA")]
clonal_plot <- snv_overview[[1]] + labs(title = "Clonal (VAF 0.25)") + 
  scale_fill_manual(paste0("Signature", strrep(" ", 14)), 
                    values = plot_cols, breaks = names(plot_cols)) +
  ggnewscale::new_scale_fill() +  # reset fill mapping for new layer
  geom_tile(data = annot_df, aes(x = Sample, y = -100, fill = Clonal), 
            height = 250, inherit.aes = FALSE, color = "black", linewidth = 0.4) +
  scale_fill_manual(values = c("Passed" = "grey30", "Bootstrap <100" = "grey70", "N/A" = "white"), guide = "none") +
  scale_y_continuous(limits = c(-260, 6200), expand = c(0.001, -0.001)) +
  theme(plot.title = element_text(face = "bold", size = 12),
        legend.title = element_text(size = 12), 
        legend.text = element_text(size = 10.5),
        axis.text.x = element_text(color = "white"), 
        axis.ticks.x = element_line(color = "white"))
subclonal_plot <- snv_overview_lowVAF[[1]] + labs(title = "Subclonal (VAF 0.05 - 0.25)") + 
  scale_fill_manual(values = plot_cols, breaks = names(plot_cols), guide = "none") +
  ggnewscale::new_scale_fill() +  # reset fill mapping for new layer
  geom_tile(data = annot_df,
            aes(x = Sample, y = -220, fill = Subclonal), 
            height = 330, inherit.aes = FALSE, color = "black", linewidth = 0.4) +
  scale_fill_manual(name = "SBS2/13 QC",  breaks = c("Passed", "Bootstrap <100", "N/A"),
                    values = c("Passed" = "grey30", "Bootstrap <100" = "grey70", "N/A" = "white")) +
  scale_y_continuous(limits = c(-440, 8500), expand = c(0.001, -0.001)) +
  theme(plot.title = element_text(face = "bold", size = 12),
        legend.title = element_text(size = 12), 
        legend.text = element_text(size = 10.5))

# Combine absolute contribution plots
pdf("AbsoluteContributionPlots_allETV6RUNX1samples_ordered_SBS213-SubclonalAndClonalVsClonal_121125.pdf", width = 13, height = 10.5)
grid.arrange(clonal_plot, subclonal_plot)
dev.off()

pdf("AbsoluteContributionPlots_allETV6RUNX1samples_ordered_SubclonalAndClonal_121125.pdf", width = 13, height = 8)
snv_overview_lowVAF[[1]] + labs(title = "All variants (VAF >=0.05)") + 
  scale_fill_manual(values = plot_cols, breaks = names(plot_cols)) +
  theme(plot.title = element_text(face = "bold", size = 12),
        legend.title = element_text(size = 11.5), 
        legend.text = element_text(size = 10.5))
dev.off()


pdf("AbsoluteContributionPlots_allETV6RUNX1samples_ordered_ClonalAndSubclonalVSsubclonal.pdf", width = 11, height = 12)
grid.arrange(snv_overview_lowVAF_2[[1]] + labs(title =  "All variants (VAF >=0.05)" ) + 
  scale_fill_manual(values = plot_cols, breaks = names(plot_cols)) +
  theme(plot.title = element_text(face = "bold", size = 12),
        legend.title = element_text(size = 12), 
        legend.text = element_text(size = 10.5),
        axis.text.x = element_text(color = "white"), 
        axis.ticks.x = element_line(color = "white")),
snv_overview_lowVAF[[1]] + labs(title = "Subclonal (VAF 0.05 - 0.25)") + 
  scale_fill_manual(values = plot_cols, breaks = names(plot_cols)) +
  theme(plot.title = element_text(face = "bold", size = 12),
        legend.title = element_text(size = 12), 
        legend.text = element_text(size = 10.5)))
dev.off()



#library(ggbreak)
#p2 <- snv_overview_lowVAF[[1]] + labs(title = "VAF 0.05 - 0.25") + 
#  scale_y_cut(breaks = c(6000, 7500), which = c(1,3), scales = c(180,500,0.000001)) + 
#  scale_y_continuous(breaks = c(seq(0,6000,1000), 8000, seq(10000, 20000, 2500 )),
#                     limits = c(0,20000))



#------------------------ Mutational profile analysis -------------------------#

# Define samples with an abs contr > 50 to SBS2/SBS13
selected_samples <- names(which(colSums(absolute_contributions[c("SBS2", "SBS13"), sample_order]) > 50))

# Generate mutational matrix of all samples
target_mut_mat <- snv_mut_mat[,sample_order]
colnames(target_mut_mat) <- paste0(colnames(target_mut_mat), " (n = ", colSums(target_mut_mat), ")")
compl_mut_profile_2 <- plot_96_profile(target_mut_mat, ymax = 0.1) + 
  xlab("Sequence context") +
  theme(strip.text.y = element_text(angle = 0, size = 10))

# Generate mutational matrix of selected samples
target_mut_mat <- snv_mut_mat[,selected_samples]
colnames(target_mut_mat) <- paste0(colnames(target_mut_mat), " (n = ", colSums(target_mut_mat), ")")
mut_profile <- plot_96_profile(target_mut_mat, ymax = 0.2) + 
  xlab("Sequence context") +
  theme(strip.text.y = element_text(angle = 0, size = 10))

# Save to PDF
pdf("mutProfiles_allETV6RUNX1samples_noVAFfilter.pdf", height = 25, width = 8)
compl_mut_profile
compl_mut_profile_2
dev.off()
pdf("mutProfiles_samplesWithAbsContrHigherThan50toSBS2_13.pdf", height = 11, width = 10)
mut_profile
dev.off()



#--------------------- Larger mutational context analysis ---------------------#

# Get larger mutational context
mut_mat_ext_context <- mut_matrix(sbs_grl_etv6runx1, ref_genome, extension = 2)

# Obtain mutation matrix of ETV6::RUNX1 samples
mut_mat_ext_context <- mut_mat_ext_context[, subtype_pos_ids]

# Randomise IDs
colnames(mut_mat_ext_context) <- random_id_file$Anonym[match(colnames(mut_mat_ext_context), 
                                                             random_id_file$Sample)]

# Generate profile
plots <- rep(list(NA), 6)
for (i in 1:6){
  selection <- c(1, 11, 21, 31, 41, 51)[i]
  to <- ifelse(i == 6, length(sample_order) - selection, 6)
  sample_subset <- sample_order[selection:(selection + to)]
  plots[[i]] <- plot_profile_heatmap(mut_mat_ext_context[,sample_subset], 
                                     by = sample_subset) +
                                     theme(strip.text.y = element_text(angle = 0),
                                           legend.position = "none",
                                           axis.text.x = element_text(size = 6),
                                           axis.text.y = element_text(size = 6))
}

# Save to PNG
png("extMutProfiles_samplesOrderedByAbsAPOBECcontr_121125.png", 
    width = 60, height = 13, units = 'in', res = 300)
grid.arrange(grobs = plots[1:6], ncol = 6)
dev.off()

# Save to Excel
mut_mat_ext_df <- cbind(data.frame(Sample = colnames(mut_mat_ext_context),
                                   SKION = subtype_pos_ids),
                        t(mut_mat_ext_context))
write.xlsx(mut_mat_ext_df, "extMutMat-ETV6RUNX1samples_121125.xlsx")


## Create correlation plot RTCA (AT or GT [C>N]NN) vs YTCA (CT or TT [C>N]NN) context
# Calclate the nr of YTCA and RTCA mutations per sample
RT_rows <- grep("^AT\\[C>.\\]..|^GT\\[C>.\\]..", rownames(mut_mat_ext_context), value = T)
YT_rows <- grep("^CT\\[C>.\\]..|^TT\\[C>.\\]..", rownames(mut_mat_ext_context), value = T)
pentanucl_df <- data.frame(Sample = colnames(mut_mat_ext_context),
                           YTCA = as.numeric(colSums(mut_mat_ext_context[YT_rows,])),
                           RTCA = as.numeric(colSums(mut_mat_ext_context[RT_rows,])))
wilcox.test(pentanucl_df$YTCA, pentanucl_df$RTCA, paired = T)
write.xlsx(pentanucl_df, 
           file = "PentanucleotideContext_APOBEC3AandAPOBEC3B_NumberOfMutatations.xlsx")

# Generate correlation plot
pentanucl_plot <- ggplot(data = pentanucl_df,
       aes(x = RTCA, y = YTCA)) +
  geom_point(size = 0.9) +
  scale_x_continuous(limits = c(0, 3550), breaks = seq(0,3550, 500), expand = c(0.01, 0)) +
  scale_y_continuous(limits = c(0, 3550), breaks = seq(0,3550, 500), expand = c(0.01, 0)) +
  geom_abline(size = 0.8, linetype = 2, color = "grey") +
  theme_classic() +
  labs(x="RT [C>N] N", y="YT [C>N] N")

# Save to PDF
pdf("PentanucleotideContext_YTCAvsRTCA_121125.pdf", width = 5.5, height = 5.5)
pentanucl_plot
dev.off()



#-------------------------- Correlation analysis DDOST ------------------------#

# Read DDOST info
DDOST_FILE <- "/Users/m.m.kleisman/Projects/HypermutatedALL_project/RESULTS/bam-readcount/completeOverview-ETV6RUNX1-DDOSThotspot-readcount.variants.xlsx"
ddost_info <- as.data.frame(read_excel(DDOST_FILE))

# Parse information
ddost_df <- data.frame(Sample = paste0("P_", ddost_info$sample),  
                       Depth = ddost_info$depth,
                       REF_count = as.numeric(sapply(strsplit(ddost_info[,8], ":"), `[`, 2)))
ddost_df["Perc_edited"] <- 100 - ((ddost_df$REF_count * 100) / ddost_df$Depth)

# Subset for ETV6::RUNX1 samples and match order to absolute contribution df
ddost_df <- ddost_df[which(ddost_df$Sample %in% etv6_runx1_samples),]
ddost_df["ID"] <- random_id_file$Anonym[match(ddost_df$Sample, random_id_file$Sample)]
ddost_df <- ddost_df[match(colnames(abs_contr), ddost_df$ID), ]

# Add clonal and subclonal SBS2/SBS13 contributions
ddost_df["Clonal"] <- colSums(abs_contr[c("SBS2", "SBS13"),])
ddost_df["Subclonal"] <- colSums(abs_contr_lowVAF[c("SBS2", "SBS13"),])
#load("ddost_df.Rda")

# Create correlation plots
clonal_plot <- ggscatter(ddost_df, 
          x = "Clonal", 
          y = "Perc_edited", 
          add = "reg.line",
          title = "DDOST editing vs absolute contribution SBS2/13 (VAF >0.25)",
          xlab = "Absolute contribution SBS2/13",
          ylab = "Percent edited RNA") +
  stat_cor(label.x = 3800, label.y = 22) +
  theme(text = element_text(size = 10.5),
        plot.title = element_text(size = 11))
subclonal_plot <- ggscatter(ddost_df, 
          x = "Subclonal", 
          y = "Perc_edited", 
          add = "reg.line",
          title = "DDOST editing vs absolute contribution SBS2/13 (VAF 0.05 - 0.25)",
          xlab = "Absolute contribution SBS2/13",
          ylab = "Percent edited RNA") +
  stat_cor(label.x = 5500, label.y = 22) +
  theme(text = element_text(size = 10.5),
        plot.title = element_text(size = 11))

pdf("DDOSTediting-vs-absContr_121125.pdf", width = 11, height = 6)
grid.arrange(clonal_plot, subclonal_plot, nrow = 1)
dev.off()

#all_plot <- ggscatter(ddost_df, 
#                      x = "All", 
#                      y = "Perc_edited", 
#                      add = "reg.line",
#                      title = "DDOST editing vs absolute contribution SBS2/13 (VAF >0.05)",
#                      xlab = "Absolute contribution SBS2/13",
#                      ylab = "Percent edited RNA") +
#  stat_cor(label.x = 7500, label.y = 22) +
#  theme(text = element_text(size = 10.5),
#        plot.title = element_text(size = 11))
#pdf("DDOSTediting-vs-absContr_inclAllMuts.pdf", width = 15.5, height = 6)
#grid.arrange(clonal_plot, subclonal_plot, all_plot, nrow = 1)
#dev.off()

#load("/Users/m.m.kleisman/Projects/HypermutatedALL_project//RESULTS/mutect2/APOBEC/mut_mat_P0917_P0918.Rda")
#pdf("~/Downloads/mut_profile_29036-29131.pdf", width = 9, height = 5)
#plot_96_profile(extra_samples_mut_mat) + xlab("Sequence context") + 
# theme(strip.text.x = element_text(size = 10), strip.text.y = element_text(size = 10))
#dev.off()



#------------------------- Comparison to RNAseq data --------------------------#

# Compute log2(TPM) values of RNAseq count data
tpm_data_df <- count_data_df[, grep("TPM", colnames(count_data_df))]
log_tpm_data_df <- log2(tpm_data_df + 1)

# Parse column names
colnames(log_tpm_data_df) <- gsub("_TPM", "", colnames(log_tpm_data_df))
colnames(log_tpm_data_df) <- paste0("P_", colnames(log_tpm_data_df))

# Subset for ETV6::RUNX1 samples with available WGS data 
subset_log_tpm_data_df <- log_tpm_data_df[,subtype_pos_ids]

# Randomize IDs and match order of WGS data
colnames(subset_log_tpm_data_df) <- names(subtype_pos_ids)[match(colnames(subset_log_tpm_data_df),
                                                                 subtype_pos_ids)]
subset_log_tpm_data_df <- subset_log_tpm_data_df[,colnames(snv_mut_mat)]

# Calculate mean APOBEC3A/3B expression for each sample
avg_gene_expr <- colMeans(subset_log_tpm_data_df[c("APOBEC3A", "APOBEC3B"),])

# Calculate contributions total for SBS2/13
abs_contr_sum <- colSums(abs_contr[c("SBS2", "SBS13"),])
rel_contr_sum <- colSums(rel_contr[c("SBS2", "SBS13"),])

# Check if all data has an identical sample order 
# Should result in 1 printed list, if not: manually adjust sample order of wrong object
unique(list(colnames(abs_contr), colnames(rel_contr), 
            colnames(snv_mut_mat), names(avg_gene_expr), colnames(subset_log_tpm_data_df)))

# Generate data frame with expression and contribution values
target_df <- data.frame(Sample = colnames(snv_mut_mat),
                        SKION = subtype_pos_ids,
                        expr_compl = avg_gene_expr,
                        expr_3a = as.numeric(subset_log_tpm_data_df["APOBEC3A",]),
                        expr_3b = as.numeric(subset_log_tpm_data_df["APOBEC3B",]),
                        abs_contr_compl = abs_contr_sum,
                        abs_contr_2 = abs_contr["SBS2", ],
                        abs_contr_13 = abs_contr["SBS13", ],
                        rel_contr_compl = rel_contr_sum,
                        rel_contr_2 = rel_contr["SBS2", ],
                        rel_contr_13 = rel_contr["SBS13", ])

# Save to Excel
excel_file <- target_df
colnames(excel_file) <- c("Random ID", "SKION nr",
                          "Mean expression APOBEC3A/B",	"Expression APOBEC3A", "Expression APOBEC3B",	
                          "Sum absolute contr. SBS2/13", "Absolute contr. SBS2",	"Absolute contr. SBS13", 
                          "Sum relative contr. SBS2/13", "Relative contr. SBS2",	"Relative contr. SBS13"	)
write.xlsx(excel_file, "correlations_APOBEC3AB-expression_vs_signatureContributions_121125.xlsx")


## Calculate correlations of expression and contribution values
# Define contribution- and expression values and descriptions
contr_types <- c("abs_contr_compl", "abs_contr_2", "abs_contr_13",
                 "rel_contr_compl", "rel_contr_2", "rel_contr_13")
expr_types <- c("expr_compl", "expr_3a", "expr_3b")
contr_descr <- c("Absolute contribution to SBS2/SBS13", "Absolute contribution to SBS2", 
                 "Absolute contribution to SBS13", "Relative contribution to SBS2/SBS13",
                 "Relative contribution to SBS2", "Relative contribution to SBS13") 
expr_descr <- c("mean APOBEC3A/B expression", "mean APOBEC3A expression", "mean APOBEC3B expression")
names(contr_descr) <- contr_types
names(expr_descr) <- expr_types

# Generate scatter plot for each contribution-expression combination 
pdf("correlations_APOBEC3AB-expression_vs_signatureContributions_121125.pdf", width = 8, height = 6)
for (contr_type in contr_types){
  for (expr_type in expr_types){
    
    # Define plot labels
    contr_type_descr <- contr_descr[names(contr_descr) == contr_type]
    expr_type_descr <- expr_descr[names(expr_descr) == expr_type]
    x_lab <- unlist(strsplit(contr_type_descr, " to"))[1]
    
    # Define text position
    if (x_lab == "Absolute contribution"){
      txt_pos <- quantile(unlist(target_df[contr_type]), probs = 0.95)
    } else {
      txt_pos <- c(0.7, 0.3, 0.4)[grep(contr_type, contr_types[4:6])]
    }
    
    # Generate scatter plot
    cor_plot <- ggscatter(target_df, 
                          x = contr_type, 
                          y = expr_type, 
                          size = 0.9,
                          add = "reg.line",
                          add.params = list(size = 0.8),
                          title = paste(contr_type_descr, "vs", expr_type_descr),
                          xlab = x_lab,
                          ylab = "Expression (in log2(TPM))") +
      stat_cor(label.x = txt_pos) +
      theme(text = element_text(size = 11))
    print(cor_plot)
  }
}
dev.off()

# Code to generate plot of single comparison (change parameters manually)
title <- "Absolute contribution to SBS2/SBS13 vs mean APOBEC3A/B expression"
x_param <- "abs_contr_compl"
y_param <- "expr_compl"
x_lab <- "Absolute contribution"
y_lab <- "Expression (in log2(TPM))"

# Define outliers
outliers <- c(target_df$Sample[order(target_df[,y_param])][1:2],
              target_df$Sample[order(target_df[,x_param], decreasing = TRUE)][1:2])

# Generate plot and save to PDF 
ggscatter(target_df, 
          x = x_param, 
          y = y_param, 
          add = "reg.line",
          title = title,
          xlab = x_lab,
          ylab = y_lab) +
  stat_cor(label.x = 2000) +
  ggrepel::geom_text_repel(aes(label = ifelse(Sample %in% outliers, Sample, "")), size = 3) +
  theme(text = element_text(size = 11))



