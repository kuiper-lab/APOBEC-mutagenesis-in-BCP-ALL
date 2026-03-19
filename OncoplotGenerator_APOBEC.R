#!/usr/bin/env Rscript 
#
#------------------------------ Data preparations ----------------------------#

# Clear data 
rm(list = ls())

# Load libraries
ref_genome <- "BSgenome.Hsapiens.UCSC.hg38"
library(ref_genome, character.only = TRUE)
library("logging")
library("data.table")
library("dplyr")
library("ComplexHeatmap")
library("RColorBrewer")
library("stringr")
library("readxl")
library("gtools")

# Set global parameters
PROJECT_DIR <- "/Users/m.m.kleisman/Projects/HypermutatedALL_project/"
RESULTS_DIR <- paste0(PROJECT_DIR, "/RESULTS/Oncoplot/APOBEC/")
MUTECT_DIR <- paste0(PROJECT_DIR, "/RESULTS/mutect2/APOBEC/")
RANDOM_ID_INFO <- "ALL_ETV6RUNX1_Coded_Num.xlsx"
DRIVER_GENE_INFO <- "NIHMS1848401-supplement-1848401_Sup_tables6.xlsx"
setwd(RESULTS_DIR)

# Load custom functions
FUNCTIONS_SCRIPT <- list.files(path = "/Users/m.m.kleisman/Projects/git/pmc_kuiper_projects/", 
                               pattern = "createdFunctions_HM-ALL.R",
                               full.names = TRUE, recursive = TRUE)
source(FUNCTIONS_SCRIPT)

# Define extension of unfiltered variant files (Only required if CONDUCT_X_FILTERING == TRUE)
CONDUCT_SNV_FILTERING <- FALSE
SNV_EXT <- "-merged-mutect2Calls.passFiltered.snp.vepAnnotated.vcf.gz" 
CONDUCT_INDEL_FILTERING <- FALSE
INDEL_EXT <- "-merged-mutect2Calls.passFiltered.indel.vepAnnotated.vcf.gz" 
CONDUCT_CNV_FILTERING <- FALSE
CNV_EXT <- "_WGS.tumor.modelFinal.annotated.txt"

# Define output files with filtered variant info
# These files will be created if CONDUCT_X_FILTERING == TRUE
# Otherwise they are loaded
COMBINED_SNV_INFO_FILE <- "Oncoplot-SNVoverview_DP20ALT5VAF025PopMaxFiltered_ETV6RUNX1_121125.xlsx"
COMBINED_INDEL_INFO_FILE <- "Oncoplot-InDelOverview_DP20ALT5VAF025PopMaxFiltered_ETV6RUNX1_121125.xlsx"
COMBINED_CNV_INFO_FILE <- "Oncoplot-CNVOverview_LOSSinDriverGenes_ETV6RUNX1_121125.xlsx"

# Load input files
centromeric_variants <- list.files(path = PROJECT_DIR,
                                   pattern = "cytoBand.hg38.centromeresOnly.txt",
                                   full.names = TRUE,
                                   recursive = TRUE)
subtype_info <- as.data.frame(read_excel(tail(list.files(path = PROJECT_DIR,
                                                         pattern = RANDOM_ID_INFO,
                                                         full.names = TRUE,
                                                         recursive = TRUE), 1)))
subtype_info <- subtype_info[subtype_info$Subtype == "ETV6::RUNX1", ]
driver_genes <- as.data.frame(read_excel(list.files(path = RESULTS_DIR,
                                                    pattern = DRIVER_GENE_INFO,
                                                    full.names = TRUE,
                                                    recursive = TRUE)))

# Parse subtype info
apobec_info <- subtype_info[, "APOBEC_group", drop = FALSE]
rownames(apobec_info) <- subtype_info$Anonym
apobec_info[apobec_info$APOBEC_group == "N/A",] <- NA

# Parse driver genes 
target_genes <- subset(driver_genes, (`Oncogene` == "yes" | 
                         `Tumor suppressor` == "yes") & `B-ALL` == "yes")[,1]
#target_genes <- c("ETV6", "TBL1XR1", "PAX5", "ATF7IP", "BTG1", "RAG2", "BTLA", "NR3C2", 
#                  "KRAS", "CDKN2A", "CDKN2B", "RB1", "EBF1", "CSF2RA", "IL3RA", "CRLF2")

# Define colours for APOBEC signal and each mutation type shown in oncoplot
apobec_col <- c(Clonal = "#000067", Subclonal = "#6bb6ff", None = "grey85")
oncoplot_mut_types <- c("focal_deletion", "frameshift_variant", "inframe_indel", "missense_variant", "stop_gained")             
#mut_type_col <- c("#419fde", "#B15928", "#5decc7",  "#FF7F00")
mut_type_col <- c("#711616", "#708fe1", "#B15928",  "#FF7F00", "gold")
names(mut_type_col) <- oncoplot_mut_types

# Define all mutation types of which consequences are accepted
all_mut_types <- c(oncoplot_mut_types, 
                   c("UTR_variant", "stop_lost", "stop_retained_variant", "start_lost", 
                     "non_coding_transcript_exon_variant", "splice_variant"))



#----------------------------- Function definitions ---------------------------#

filterVCFfile <- function(input_vcf,
                          ref_genome = "BSgenome.Hsapiens.UCSC.hg38",
                          DP = 20,
                          ALT = 3,
                          VAF = 0.25){
  # Read the VCF file
  vcf_file <- loadVcfFile(vcfFilePath = input_vcf,
                          bsgGenomeName = ref_genome,
                          chromosomeSet = "chromosomes")
  GenomeInfoDb::genome(vcf_file) <- "hg38"
  
  # Get tumor ID from VCF file
  tumor_sample_id_vcf_line <- grep("##tumor_sample=", 
                                   readLines(input_vcf), 
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
  df_somatic_filtered <- convertVcfObjectToDataFrame(vcf_gonl_filtered)
  
  # Get column idx for tumor AF, ref and alt read counts
  AF_idx <- which(names(df_somatic_filtered) %in% paste0(tumor_sample_id, "_read_based_AF"))
  read_count_ref_idx <- which(names(df_somatic_filtered) %in% 
                                paste0(tumor_sample_id, "_read_count_ref"))
  read_count_alt_idx <- which(names(df_somatic_filtered) %in% 
                                paste0(tumor_sample_id, "_read_count_alt"))
  
  # Compute the total depth 
  df_somatic_filtered[paste0(tumor_sample_id, "_total_read_count")] <- 
    df_somatic_filtered[read_count_ref_idx] + df_somatic_filtered[read_count_alt_idx]
  total_read_count_idx <- ncol(df_somatic_filtered) 
  
  # Filter variants on read depth
  df_read_based_filtered <- df_somatic_filtered %>% 
    filter(.data[[c(names(df_somatic_filtered)[read_count_alt_idx])]] >= ALT &
             .data[[c(names(df_somatic_filtered)[total_read_count_idx])]] >= DP &
             .data[[c(names(df_somatic_filtered)[AF_idx])]] >= VAF)
  return(df_read_based_filtered)
}


## Wrapper function to determine if a string starts with any pattern
# Returns TRUE if the string starts with any of the patterns 
stringStartsWith <- function(patterns, string){
  grep_pattern <- paste0("^", paste0(patterns, collapse = "|^"))
  return(grepl(grep_pattern, string))
}


## Function to determine the mutation types in all genes for each sample
FillMatrix <- function(gene.list, sample.IDs, combine.file){
  
  ## Create input matrix 
  mat <- matrix("",nrow = length(gene.list), ncol = length(sample.IDs),
                dimnames = list(gene.list, sample.IDs))

  ## Fill matrix with mutation types per gene for each sample
  for (i in 1:nrow(combine.file)){
    
    # Get the generic info 
    ID <- combine.file[i,"SampleID"]
    Gene <- combine.file[i,"SYMBOL"]
    mut.type <- combine.file[i, "Consequence"]
    
    # Determine if the variant needs to be skipped (3'-UTR, 5'-UTR, upstream, downstream, synonymous, intron)
    if (stringStartsWith(c("5_prime_UTR_variant", 
                           "3_prime_UTR_variant", 
                           "upstream_gene_variant", 
                           "downstream_gene_variant",
                           "synonymous_variant",
                           "intron_variant"), mut.type)){
      next
    }
    
    # Determine of which type the variant is
    # Missense variant
    else if (stringStartsWith("missense_variant", mut.type)){
      mut_type <- "missense_variant"
    }
    # Stop gained variant
    else if (stringStartsWith("stop_gained", mut.type)){
      mut_type <- "stop_gained"
    }
    # Splice variant
    else if (stringStartsWith(c("splice_donor_variant", 
                                "splice_region_variant", 
                                "splice_acceptor_variant"), 
                              mut.type)){
      
      # Check if the splice variant meets splice AI thresholds
      if (all(is.na(combine.file[i, grep("Splice", colnames(combine.file))]))) {
        next
      } else if ( (combine.file[i, ]$SpliceAI_pred_DS_AG > 0.5) | 
           (combine.file[i, ]$SpliceAI_pred_DS_AL > 0.5) | 
           (combine.file[i, ]$SpliceAI_pred_DS_DG > 0.5) | 
           (combine.file[i, ]$SpliceAI_pred_DS_DL > 0.5) ){
        mut_type <- "splice_variant"
      } else {
        next
      }
    }
    # Stop retained variant
    else if (stringStartsWith("stop_retained_variant", mut.type)){
      mut_type <- "stop_retained_variant"
    }
    # Non-coding transcript exon variant
    else if (stringStartsWith("non_coding_transcript_exon_variant", mut.type)){
      mut_type <- "non_coding_transcript_exon_variant"
    }
    # Start lost variant
    else if (stringStartsWith("start_lost", mut.type)){
      mut_type <- "start_lost"
    }
    # Stop lost variant
    else if (stringStartsWith("stop_lost", mut.type)){
      mut_type <- "stop_lost"
    }
    # Inframe indel
    else if (stringStartsWith(c("inframe_insertion", 
                                "inframe_deletion",
                                "inframe_indel"), 
                              mut.type)){
      mut_type <- "inframe_indel"
    } 
    # Frameshift variant
    else if (stringStartsWith("frameshift_variant", mut.type)){
      mut_type <- "frameshift_variant"
    }
    # Focal deletion (CNV loss)
    else if (stringStartsWith("LOSS", mut.type)){
      mut_type <- "focal_deletion"
    } 
    # No variant found, print warning message
    else {
      print(paste0("WARNING: mutation type not found: ", mut.type))
      next
    }
    
    # Add the defined mutation type to the mutation type matrix
    if (mat[Gene, ID] == ""){
      mat[Gene, ID] <- mut_type
    } else {
      mat[Gene, ID] <- paste0(mat[Gene, ID], ";", mut_type)
    }
  }
  return(mat)
}


## Function to count the percentage of mutation types per gene
CountMutPerGenePerc <- function (gene.list, mat){
  
  # Create empty matrix to store percentages
  n_samples = length(colnames(mat))
  n_mutated_patients <- matrix(0, nrow = length(c(gene.list)), ncol = length(all_mut_types))
  dimnames(n_mutated_patients) <- list(c(gene.list), all_mut_types)
  
  # Iterate over each gene to calculate frequencies of mutation types
  for (gene in gene.list){
    gene_mut_types <- mat[gene,][as.character(mat[gene,]) != ""]
    gene_mut_freq <- cbind(table(as.character(gene_mut_types)))
    
    # Check if a gene has multiple mutation types in a single sample
    if (any(grepl(";", gene_mut_types))){
      
      # Determine the (number of) mutation types that are found in a single sample
      multi_mut_types <- gene_mut_freq[grepl(";", rownames(gene_mut_freq)), , drop = FALSE]
      #multi_mut_types_n <- apply(multi_mut_types, 2, function(x)  str_count(names(x), ";"))
      multi_mut_types_n <- sapply(rownames(multi_mut_types), function(x)  length(unique(unlist(strsplit(x, ";")))) - 1)
      
      # Subset frequencies for mutation types that are only found once in a sample
      gene_mut_freq <- as.data.frame(gene_mut_freq[!grepl(";", rownames(gene_mut_freq)), , drop = FALSE])
      
      # Determine freq. of multiple mut. types that are found in a single sample (1 / n_mut_types)
      for (i in 1:nrow(multi_mut_types)){
        sample_mut_types <- unique(unlist(strsplit(names(multi_mut_types[i,]), ";")))
        n_multi_samples <- as.numeric(multi_mut_types[i,])
        
        for (mut_type in sample_mut_types){
          if (mut_type %in% rownames(gene_mut_freq)){
            gene_mut_freq[mut_type, ] <- gene_mut_freq[mut_type, ] + 
              ((1 / (multi_mut_types_n[i] + 1)) * n_multi_samples)
          } else {
            gene_mut_freq[nrow(gene_mut_freq) + 1, ] <- (1 / (multi_mut_types_n[i] + 1)) * n_multi_samples
            rownames(gene_mut_freq)[nrow(gene_mut_freq)] <- mut_type
          }
        }
      }
    } else {
      gene_mut_freq <- as.data.frame(gene_mut_freq)
    }
    
    # Calculate the percentage of mutation types found in the gene
    for (mut_type in rownames(gene_mut_freq)){
      mut_type_freq <- gene_mut_freq[mut_type,]
      n_mutated_patients[gene, mut_type] <- (mut_type_freq / n_samples) * 100
    }
  }
  return(n_mutated_patients)
}


# Code for a grey background in oncoplot
custom_alter_fun <- function(fill,
                             background_color = "#e5e5e5",
                             col = NA) {
  function(x, y, w, h, v) {
    n = sum(v)  # how many alterations for current gene in current sample
    h = h*0.9
    
    # use `names(which(v))` to correctly map between `v` and `col`
    if(n) {grid::grid.rect(x, 
                           y - h*0.5 + 1:n/n*h,
                           w*0.9,
                           1/n*h,
                           gp = gpar(fill = fill[names(which(v))], col = col),
                           just = "top")}
    if(!n) {grid::grid.rect(x,
                            y,
                            w*0.9,
                            h,
                            gp = gpar(fill = background_color, col = background_color),
                            just = "center")}
  }
}



#------------------------------ Variant filtering ----------------------------#

# Define target columns
target_cols <- c("chr", "start", "end", "ref", "alt", "Consequence", "IMPACT", 
                 "SYMBOL", "Gene", "Feature_type", "Feature", "BIOTYPE", "CADD_PHRED", 
                 "CADD_RAW", "SpliceAI_pred_DS_AG", "SpliceAI_pred_DS_AL", "SpliceAI_pred_DS_DG", 
                 "SpliceAI_pred_DS_DL", "gnomADg_AF", "SampleID")


## Conduct SNV filtering
if (CONDUCT_SNV_FILTERING){
  combined_snv_df <- data.frame()
  
  # List input SNV files
  snv_files <- list.files(path = MUTECT_DIR,
                          pattern = SNV_EXT,
                          recursive = TRUE,
                          full.names = TRUE)
  
  # Obtain sample IDs
  skion_ids <- gsub(SNV_EXT, "", basename(snv_files))
  sample_idx <- skion_ids %in% subtype_info$Sample
  
  # Subset files for target samples
  snv_files <- snv_files[sample_idx]
  skion_ids <- skion_ids[sample_idx]
  loginfo(paste("found:", length(skion_ids), "SNV VCF files"))

  # Filter each VCF file
  for (i in 1:length(snv_files)){
    df_read_based_filtered <- filterVCFfile(input_vcf = snv_files[i])
    
    # Subset the data frame for target columns
    df_read_based_filtered["SampleID"] <- skion_ids[i]
    filtered_df <- df_read_based_filtered[,target_cols]
    
    # Add to data frame
    combined_snv_df <- rbind(combined_snv_df, filtered_df)
  } 
  # Write filtered data frame to Excel
  openxlsx::write.xlsx(combined_snv_df, file = COMBINED_SNV_INFO_FILE)
  
} else {
  print("SNV filtering not requested, filtered Excel file will be used for further analyses.")
  combined_snv_df <- as.data.frame(read_excel(list.files(pattern = COMBINED_SNV_INFO_FILE,
                                                         full.names = TRUE,
                                                         recursive = TRUE)))
}


## Conduct InDel filtering
if (CONDUCT_INDEL_FILTERING){
  combined_indel_df <- data.frame()
  
  # List input SNV files
  indel_files <- list.files(path = MUTECT_DIR,
                            pattern = INDEL_EXT,
                            recursive = TRUE,
                            full.names = TRUE)
  
  # Obtain sample IDs
  skion_ids <- gsub(INDEL_EXT, "", basename(indel_files))
  sample_idx <- skion_ids %in% subtype_info$Sample
  
  # Subset files for target samples
  indel_files <- indel_files[sample_idx]
  skion_ids <- skion_ids[sample_idx]
  loginfo(paste("found:", length(skion_ids), "InDel VCF files"))
  
  # Filter each VCF file
  for (i in 1:length(indel_files)){
    df_read_based_filtered <- filterVCFfile(input_vcf = indel_files[i])
    
    # Subset the data frame for target columns
    df_read_based_filtered["SampleID"] <- skion_ids[i]
    filtered_df <- df_read_based_filtered[,target_cols]
    
    # Add to data frame
    combined_indel_df <- rbind(combined_indel_df, filtered_df)
  } 
  # Write filtered data frame to Excel
  openxlsx::write.xlsx(combined_indel_df, file = COMBINED_INDEL_INFO_FILE)
  
} else {
  print("InDel filtering not requested, filtered Excel file will be used for further analyses.")
  combined_indel_df <- as.data.frame(read_excel(list.files(pattern = COMBINED_INDEL_INFO_FILE,
                                                           full.names = TRUE,
                                                           recursive = TRUE)))
}


## Conduct CNV filtering
if (CONDUCT_CNV_FILTERING) {
  combined_cnv_df <- data.frame()
  
  # List input CNV files
  cnv_files <- list.files(path = RESULTS_DIR,
                            pattern = CNV_EXT,
                            recursive = TRUE,
                            full.names = TRUE)
  
  # Obtain sample IDs
  sample_ids <- gsub(CNV_EXT, "", basename(cnv_files))
  
  # Convert BM IDs to SKION IDs
  skion_ids <- sapply(sample_ids, function(id){ifelse(grepl("^PM", id), 
                                                      subtype_info$Sample[match(id, subtype_info$BM_ID)],
                                                      paste0("P_", gsub("Dx", "", id)))})
  sample_idx <- skion_ids %in% subtype_info$Sample
  
  # Subset files for target samples
  cnv_files <- cnv_files[sample_idx]
  skion_ids <- skion_ids[sample_idx]
  loginfo(paste("found:", length(skion_ids), "CNV files"))

  # Read each CNV file and check for losses in target genes
  for (i in (1:length(cnv_files))){
    cnv_file <- read.delim(file = cnv_files[i], 
                           sep = "\t", 
                           header = TRUE, 
                           row.names = NULL)
    colnames(cnv_file) <- c("GENE_NAME", "CHR", "GENE_START", "GENE_STOP", 
                            "CNV_START", "CNV_STOP", "MEAN_COPY_RATIO", 
                            "MINOR_ALLELE_FRACTION_POSTERIOR_90", "OVERLAP",
                            "CNV_CALL", "BAF_CALL", "CNV_LENGTH", "OVERLAPPING_GENES_#", 
                            "OVERLAPPING_GENES")
    cnv_file$SampleID <- skion_ids[i]
    
    # Only grep the genes of interest with
    cnv_file_sub <- cnv_file[cnv_file$GENE_NAME %in% target_genes &
                             cnv_file$CNV_CALL == "LOSS" & 
                             cnv_file$MEAN_COPY_RATIO <= 0.75, ]
    
    
    if (nrow(cnv_file_sub) > 0){
      combined_cnv_df <- rbind(combined_cnv_df, cnv_file_sub)
    } else {
      print(paste("No CNV losses found in sample", skion_ids[i]))
    }
  }
  # Write filtered data frame to Excel
  openxlsx::write.xlsx(combined_cnv_df, file = COMBINED_CNV_INFO_FILE)
  
} else {
  print("CNV filtering not requested, filtered Excel file will be used for further analyses.")
  combined_cnv_df <- as.data.frame(read_excel(list.files(pattern = COMBINED_CNV_INFO_FILE,
                                                         full.names = TRUE,
                                                         recursive = TRUE)))
}

# Parse CNV data in format of SNV and Indel data
combined_cnv_df_parsed <- as.data.frame(matrix(nrow = nrow(combined_cnv_df), ncol = length(target_cols)))
colnames(combined_cnv_df_parsed) <- target_cols
combined_cnv_df_parsed[,c("SampleID", "Consequence", "SYMBOL")] <- combined_cnv_df[,c("SampleID", "CNV_CALL", "GENE_NAME")]


## Combine filtered variants
combined_filtered_df <- rbind(combined_snv_df, combined_indel_df, combined_cnv_df_parsed) 

# Randomise IDs
combined_filtered_df$SampleID <- subtype_info$Anonym[match(combined_filtered_df$SampleID, subtype_info$Sample)]

# Subset data frame for target genes
combined_filtered_df <- combined_filtered_df[combined_filtered_df$SYMBOL %in% target_genes,]

# Rename CDKN2A and CDKN2B TO CDKN2A/B
combined_filtered_df[combined_filtered_df$SYMBOL %in% c("CDKN2A", "CDKN2B"), "SYMBOL"] <- "CDKN2A/B"
target_genes <- c(target_genes[!grepl("CDKN2", target_genes)], "CDKN2A/B")

# Rename CSF2RA, IL3RA and CRLF2 to PAR1
#combined_filtered_df[combined_filtered_df$SYMBOL %in% c("CSF2RA", "IL3RA", "CRLF2"), "SYMBOL"] <- "PAR1"
#target_genes <- c(target_genes[!target_genes %in% c("CSF2RA", "IL3RA", "CRLF2")], "PAR1")

# Remove duplicate entries from data frame due to gene renaming/concatenation
combined_filtered_df <- combined_filtered_df[!duplicated(combined_filtered_df), ]



#------------------------------ Oncoplot generation ---------------------------#

ORDER_BY_APOBEC_LEVEL <- TRUE

# Fill matrix with mutation frequencies and consequences
mat <- FillMatrix(gene.list = target_genes, 
                  sample.IDs = subtype_info$Anonym, 
                  combine.file = combined_filtered_df)


## Subset the data for genes with mutations in at least 2 samples 
n_filter <- length(subtype_info$Anonym) / 10
n_mut_occurrences <- data.frame(length(subtype_info$Anonym) - rowSums(mat == ""))
colnames(n_mut_occurrences) <- c("n_mutations")
n_mut_occurrences$gene <- rownames(n_mut_occurrences)
gene_list <- n_mut_occurrences[n_mut_occurrences$n_mutations >= n_filter ,]$gene

# Recreate matrix for genes with mutations in at least 2 samples
subset_combined_filtered_df <- combined_filtered_df[combined_filtered_df$SYMBOL %in% gene_list, ]
mat <- FillMatrix(gene.list = gene_list, 
                  sample.IDs = subtype_info$Anonym, 
                  combine.file = subset_combined_filtered_df)

# Calculate the frequency of mutations per gene for each class of mutation
n_mutated_patients_perc <- CountMutPerGenePerc(gene.list = gene_list, 
                                               mat = mat)
gene_order <- data.frame(sort(rowSums(n_mutated_patients_perc), 
                              decreasing = TRUE))

# Create two bars for the right annotation to calculate % per gene per sample group
pos_samples <- rownames(apobec_info)[which(apobec_info$APOBEC_group == "Clonal" | 
                                           apobec_info$APOBEC_group == "Subclonal")]
neg_samples <- rownames(apobec_info)[which(apobec_info$APOBEC_group == "None")]
n_mut_high <- rowSums(mat[, colnames(mat) %in% pos_samples] != "") 
n_mut_none <- rowSums(mat[, colnames(mat) %in% neg_samples] != "")
perc_mut_high <- n_mut_high / length(pos_samples) * 100
perc_mut_none <- n_mut_none / length(neg_samples) * 100
bar_data <- cbind(Pos_APOBEC = perc_mut_high, Neg_APOBEC = perc_mut_none)

# Calculate significance for each gene 
cont_table <- rbind(data.frame(Signal = "Clonal", Gene = rownames(mat),
                               Yes = n_mut_high, No = length(pos_samples) - n_mut_high),
                    data.frame(Signal = "None", Gene = rownames(mat),
                               Yes = n_mut_none, No = length(neg_samples) - n_mut_none))
by_gene <- split(cont_table, cont_table$Gene)
res_list <- lapply(by_gene, function(gene_df) {
  p_val <- data.frame(Gene = rownames(gene_df)[1],
                      p_val = fisher.test(gene_df[,c("Yes", "No")])$p.value)})
p_vals <- do.call(rbind, res_list)
print(paste("Significant genes:", paste(p_vals$Gene[which(p_vals$p_val <= 0.05)], collapse = ", ")))

# Define sample order
if (ORDER_BY_APOBEC_LEVEL) {
  pdf_ext <- "_orderedByAPOBEClevel_newComparisonGroups_121125.pdf"
  sample_order <- subtype_info$Anonym[order(as.numeric(subtype_info$APOBEC_order))]
  sample_order <- sample_order[sample_order %in% colnames(mat)]
} else {
  pdf_ext <- "_orderedByMutualExclusivity.pdf"
  sample_order <- colnames(mat)
}

# Reorder matrix and subtype info based on sample order
mat <- mat[,sample_order]
apobec_info_plot <- apobec_info[sample_order,, drop = FALSE]
colnames(apobec_info_plot) <- "APOBEC level"
apobec_info_plot$`APOBEC level` <- factor(apobec_info_plot$`APOBEC level`, 
                                          levels = c("Clonal", "Subclonal", "None"))
if (ORDER_BY_APOBEC_LEVEL){sample_order <- sample_order} else {sample_order <- NULL}

# Generate the oncoplot
onco_plot <- oncoPrint(mat,
                       alter_fun = custom_alter_fun(fill = mut_type_col), col = mut_type_col, 
                       row_names_gp = gpar(fontsize = 8),
                       column_names_gp = gpar(fontsize = 8),
                       top_annotation = NULL,   # comment this line if you want to plot top_annotation 
                       bottom_annotation = HeatmapAnnotation(df = apobec_info_plot, na_col = "white",
                                                             col = list(`APOBEC level` = apobec_col), 
                                                             annotation_name_gp = gpar(fontsize = 7, fontface = "bold"), 
                                                             simple_anno_size_adjust = TRUE, 
                                                             height = unit(0.5, "cm")), 
                       right_annotation = rowAnnotation(" " = anno_barplot(bar_data, ylim = c(0, 65),
                           gp = gpar(fill = c("#305a53", "#b2d6d0"), col = NA),  # set colors per group
                           border = FALSE, bar_width = 0.5,
                           axis = TRUE, beside = TRUE,
                           axis_param = list(side = "top")), 
                           width = unit(2, "cm")),
                       show_pct = FALSE,
                       row_gap = unit(2, "cm"), 
                       show_column_names = TRUE,
                       row_order = rownames(gene_order),
                       column_order = sample_order)

# Generate mut burden legend manually
burden_legend = Legend(labels = c("(Sub)clonal APOBEC\nmutagenesis", "No APOBEC\nmutagenesis"),
                       title = "Gene mutation burden",
                       legend_gp = gpar(fill = c("#305a53", "#b2d6d0"), col = NA))

# Save combined plots to PDF
pdf(paste0("Oncoplot_DepthFiltered-popMax0.01-CADD15-newDriverGenes_ETV6RUNX1top10", pdf_ext), width = 9.5, height = 5) 
draw(onco_plot, merge_legends = TRUE, 
     annotation_legend_list = list(burden_legend))
dev.off()

# Save gene mutation frequencies to out file
mut_freq_info <- cont_table[sort(rownames(cont_table)), c("Gene", "Signal", "Yes", "No")]
mut_freq_info$Signal <- gsub("Clonal", "(Sub)clonal", mut_freq_info$Signal)
rownames(mut_freq_info) <- NULL
mut_freq_info["Percentage"] <- (mut_freq_info$Yes * 100) / c(mut_freq_info$Yes + mut_freq_info$No)
write.xlsx(mut_freq_info, file = "GeneMutationFrequencies_SubClonalVsNone.xlsx")



