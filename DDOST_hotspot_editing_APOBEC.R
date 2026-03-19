# Load libraries
library("tidyverse")
library("openxlsx")
library("ggplot2")
library("ggpubr")
library("ggbreak")
library("patchwork")
library("pdftools")

# Define global parameters
PROJECT_DIR <- "/Users/m.m.kleisman/Projects/HypermutatedALL_project/"
RESULTS_DIR <- paste0(PROJECT_DIR, "RESULTS/bam-readcount/")
SUBTYPE_INFO <- "ALL_ETV6RUNX1_Coded_Num.xlsx"
BAM_READCOUNT_RESULT <- "completeOverview-AllSubtypes-DDOSThotspot-readcount.variants.combinedAndParsed.xlsx"
setwd(RESULTS_DIR)

# Load APOBEC status info
apobec_info <- as.data.frame(read_excel(tail(list.files(path = PROJECT_DIR,
                                                        pattern = SUBTYPE_INFO,
                                                        full.names = TRUE,
                                                        recursive = TRUE), 1)))
apobec_info <- apobec_info[apobec_info$Subtype == "ETV6::RUNX1",]
pos_samples <- apobec_info$Sample[apobec_info$APOBEC_group %in% c("Clonal", "Subclonal")]
neg_samples <- apobec_info$Sample[apobec_info$APOBEC_group == "None"]
lowqual_samples <- setdiff(apobec_info$Sample, c(pos_samples, neg_samples))

# Load combined bam-readcount data
DDOST_workfile <- as.data.frame(read_excel(BAM_READCOUNT_RESULT))

# Parse subtype annotations
DDOST_workfile$subtype <- gsub(" positive", "", DDOST_workfile$subtype)

# Each nucleotide comes with long set of data in the bam-readcount output,
# e.g. A:2:60.00:30.00:0.00:1:1:0.93:0.01:30.00:1:0.37:95.00:0.32
# Extract the number of reads per nucleotide (first integer) 
DDOST_workfile$A<- as.integer(str_match(DDOST_workfile$`base:stats...6`,"A:\\s*(.*?)\\s*:")[,2])
DDOST_workfile$C<-as.integer(str_match(DDOST_workfile$`base:stats...7`,"C:\\s*(.*?)\\s*:")[,2])
DDOST_workfile$G<-as.integer(str_match(DDOST_workfile$`base:stats...8`,"G:\\s*(.*?)\\s*:")[,2])
DDOST_workfile$T<-as.integer(str_match(DDOST_workfile$`base:stats...9`,"T:\\s*(.*?)\\s*:")[,2])


## Determine the percentage of edited RNA
# Calculate percentage of mutated bases G (ref) : A (alt)
DDOST_workfile$Mutated<-round(DDOST_workfile$A/DDOST_workfile$G*100, digits = 2)

# Add column for fill color
DDOST_workfile$isER<-"No"
DDOST_workfile$isER[DDOST_workfile$subtype=="ETV6::RUNX1"]<-"Yes"

# Add info of the number of samples per subtype
counts <- table(DDOST_workfile$subtype)
x_labeling <- scale_x_discrete(labels = setNames(paste0(names(counts), " n = ", counts),
                               names(counts)))

# Sort by number of patients per subtype
plot_order <- names(counts[order(counts, decreasing = T)])
DDOST_workfile$subtype <- factor(DDOST_workfile$subtype, levels = plot_order)

# Generate plot
p1 <- ggplot(data = DDOST_workfile, aes(x = subtype, y = Mutated, fill = isER)) +
  geom_boxplot() +
  labs(x = NULL, y = "Percent edited DDOST RNA") +
  x_labeling +
  scale_fill_manual(values = c("lightgrey", "#CE2627")) +
  scale_y_continuous(breaks = c(0, 5, 10, 20, 30), limits = c(0, 35)) +
  scale_y_break(c(10, 20), space = 0.02, scales = 0.5) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title.y = element_text(size = 10, angle = 90),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.95),
        legend.position = "none") 
p1


## Subset data for ETV6::RUNX1 samples 
ER_DDOST_workfile <- DDOST_workfile[DDOST_workfile$subtype == "ETV6::RUNX1" &
                                    !DDOST_workfile$sample %in% gsub("P_", "", lowqual_samples),]

# Add APOBEC status
ER_DDOST_workfile["APOBEC"] <- "RNAseq only"
ER_DDOST_workfile$APOBEC[ER_DDOST_workfile$sample %in% gsub("P_", "", pos_samples)] <- "APOBEC-positive"
ER_DDOST_workfile$APOBEC[ER_DDOST_workfile$sample %in% gsub("P_", "", neg_samples)] <- "APOBEC-negative"
ER_DDOST_workfile$APOBEC <- factor(ER_DDOST_workfile$APOBEC, 
                                   levels = c("APOBEC-positive", "APOBEC-negative", "RNAseq only"))

# Define the number of samples per group
x_labeling <- scale_x_discrete(labels = function(x) {
                               counts <- table(ER_DDOST_workfile$APOBEC)
                               paste0("ETV6::RUNX1\n", x, " n = ", counts[x])})

# Generate plot
p2 <- ggplot(data = ER_DDOST_workfile, aes(x = APOBEC, y = Mutated, fill = APOBEC)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#CE2627", "lightblue", "grey")) +
  stat_compare_means(comparisons = list(c("APOBEC-positive", "APOBEC-negative")),
                     method = "t.test", label = "p.format",
                     size = 3, tip.length = 0.02) +
  labs(x = NULL, y = "Percent edited DDOST RNA") +
  x_labeling +
  scale_y_continuous(breaks = c(0, 5, 10, 20, 30), limits = c(0, 35)) +
  scale_y_break(c(10, 20),space = 0.02, scales = 0.5) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title.y = element_text(size = 10, angle = 90),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.95),
        legend.position = "none") 
p2


## Combine the two plots (use patchwork to respect custom y-axis)
pdf("DDOSTediting-tmp.pdf", width = 9, height = 6.2)
p1 + p2 + plot_layout(widths = c(1, 0.35))
dev.off()

# Remove first empty page (bug of using patchwork)
pdf_subset(input = "DDOSTediting-tmp.pdf", pages = 2, 
           output = "DDOSTediting-overview_new.pdf")
file.remove("DDOSTediting-tmp.pdf")


## Significance tests for ETV6::RUNX1 samples versus other subtypes
# Only include subtypes with at least 5 samples
other_subtypes <- names(which(counts > 4))
other_subtypes <- other_subtypes[!grepl("ETV6::RUNX1", other_subtypes)]

# Determine p-values
for (subtype in other_subtypes){
  cat("\n\nETV6::RUNX1 vs", subtype, "samples\n")
  p_val <- t.test(DDOST_workfile[DDOST_workfile$subtype == "ETV6::RUNX1", "Mutated"],
                  DDOST_workfile[DDOST_workfile$subtype == subtype, "Mutated"])$p.value 
  cat("P-value t-test:", p_val)
}

