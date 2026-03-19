library(GenomicRanges)

run <- function(input_file, bed_file, aux_file, protein_coding_genes_bed, out_file) {
home_dir = "/hpc/pmc_kuiper/iTHER_germline/REFERENCE/" #project main directory
# input_file = "/Users/marcelwillemsen/workspace_R/CNV_annotation/DATA/CNV_PMC_WESKidTs/PMRBM000AHP_WXS.modelFinal.seg"
# input_file = "/Users/marcelwillemsen/workspace_Kuiper/SCNV_Janna/RESULTS/batch4/3e663fdd-330c-4060-8190-4571602b706f/call-GatherCNVresults/PMRBM000AGF_WXS.modelFinal.seg"

  #..define dirs..
  #gene_filter_dir = paste0(home_dir, "gene_panel_filter.bed")
  gene_filter_dir = paste0(bed_file)
  
  #..load libraries..
  source(paste(aux_file))
  
  #-----------------------------------------------------------------------------
  # CNV annotation
  #-----------------------------------------------------------------------------
  #..load gene filter..
  gene_filter = read.table(gene_filter_dir, header=FALSE, sep="\t")
  
  tot_g = nrow(gene_filter)
  
  #..check genes affected by CNV..
  # all_f = list.files(input_data_dir, ".modelFinal.seg$", full.names=TRUE)
  # tot_f = length(all_f)
  # for(f in 1:tot_f){
  #   print("-------")
  #   print(f)
    #..load segment data..
    # fname = all_f[f]
  
  con = file(input_file, "r")
  c_data = readLines(con)
  close(con)
  
  #..parse segment data..
  tot_at = length(grep("@", substr(c_data,1,1)))#..determine line with the header of the segment file..
  c_h = unlist(strsplit(c_data[tot_at+1], "\t"))#..get data header..
  
  #..get segment data..
  tot_h = length(c_h)
  seg_data = array(0, dim=c(0, tot_h))
  for(l in (tot_at+2):length(c_data)){
  	new_seg = unlist(strsplit(c_data[l], "\t"))
  	seg_data = rbind(seg_data, new_seg)
  }
  colnames(seg_data) = c_h
  
  
  all_res = array(0, dim=c(0,12))
  for(g in 1:tot_g){
  	# g_chr = as.character(gene_filter$Chromosome[g])
    g_chr = as.character(gene_filter[g,1]) # CHROM
  	# g_sta = gene_filter$Start_bp[g]
    g_sta = as.numeric(gene_filter[g,2]) + 1 # START, bed file zero based
    # g_sto = gene_filter$Stop_bp[g]
    g_sto = as.numeric(gene_filter[g,3]) + 1 # STOP, bed file zero based
  
  	chr_flag = seg_data[,1]%in%g_chr
  	full_gene_flag = (as.numeric(seg_data[,2])<g_sta) & (as.numeric(seg_data[,3])>g_sto)
  	out_in_gene_flag = (as.numeric(seg_data[,2])<g_sta) & (as.numeric(seg_data[,3])>=g_sta & as.numeric(seg_data[,3])<=g_sto)
  	in_out_gene_flag = (as.numeric(seg_data[,2])>=g_sta & as.numeric(seg_data[,2])<=g_sto) & (as.numeric(seg_data[,3])>g_sto)
  	in_gene_flag = (as.numeric(seg_data[,2])>=g_sta) & (as.numeric(seg_data[,3])<=g_sto)
  	g_flag = chr_flag & (full_gene_flag | out_in_gene_flag | in_out_gene_flag | in_gene_flag)
  
  	tot_in_g = sum(g_flag)
  	if(tot_in_g>0){
  		# c_g_data = c(as.character(gene_filter[g, 1]), as.character(gene_filter[g, 2]), as.character(gene_filter[g, 3:4]), as.character(gene_filter[g, 5]))
  	  c_g_data = c(as.character(gene_filter[g, 4]), as.character(gene_filter[g, 1]), as.character(gene_filter[g, 2:3]),"?")
  		c_in_g_data = seg_data[g_flag,]
  
  		if(tot_in_g==1){
  			#..get category..
  			if(full_gene_flag[g_flag]){
  				c_class = "FullGene"
  			}
  			if(out_in_gene_flag[g_flag]){
  				c_class = "Terminal5"
  			}
  			if(in_out_gene_flag[g_flag]){
  				c_class = "Terminal3"
  			}
  			if(in_gene_flag[g_flag]){
  				c_class = "InsideGene"
  			}
  
  			#..check allele frequencies..
  			MINOR_ALLELE_FRACTION_POSTERIOR_90 = as.numeric(c_in_g_data[length(c_in_g_data)])
  			c_baf_tag = get_baf_tag(MINOR_ALLELE_FRACTION_POSTERIOR_90)
  
  			#..get call tag..
  			MEAN_COPY_RATIO = 2^(as.numeric(c_in_g_data[8]))
  			c_call_tag = get_call_tag(MEAN_COPY_RATIO)
  
  			new_res = c(c_g_data, c_in_g_data[c(2, 3, 8, 10)], c_class, c_call_tag, c_baf_tag)
  			all_res = rbind(all_res, new_res)
  		}else{
  			for(s in 1:tot_in_g){
  				#..get category..
  				if((full_gene_flag[g_flag])[s]){
  					c_class = "FullGene"
  				}
  				if((out_in_gene_flag[g_flag])[s]){
  					c_class = "Terminal5"
  				}
  				if((in_out_gene_flag[g_flag])[s]){
  					c_class = "Terminal3"
  				}
  				if((in_gene_flag[g_flag])[s]){
  					c_class = "InsideGene"
  				}
  
  				#..check allele frequencies..
  				if((full_gene_flag[g_flag])[s]){
  					c_class = "FullGene"
  				}
  
  				#..check allele frequencies..
  				MINOR_ALLELE_FRACTION_POSTERIOR_90 = as.numeric(c_in_g_data[s,length(c_in_g_data[s,])])
  				c_baf_tag = get_baf_tag(MINOR_ALLELE_FRACTION_POSTERIOR_90)
  
  				#..get call tag..
  				MEAN_COPY_RATIO = 2^(as.numeric(c_in_g_data[s,8]))
  				c_call_tag = get_call_tag(MEAN_COPY_RATIO)
  
  				new_res = c(c_g_data, c_in_g_data[s, c(2, 3, 8, length(c_in_g_data[s,]))], c_class, c_call_tag, c_baf_tag)
  				all_res = rbind(all_res, new_res)
  
  			}
  		}
  	}
  }
  #..convert log2 copy ratios to absolute values..
  all_res[, 8] = signif(2^(as.numeric(all_res[, 8])), digits=3)
  
  #..round minor allele fractions..
  all_res[, 9] = signif(as.numeric(all_res[, 9]), digits=3)
  
  #..add column names..	
  colnames(all_res) = c("GENE_NAME", "CHR", "GENE_START", "GENE_STOP", "GENE_STRAND", "CNV_START", "CNV_STOP", "MEAN_COPY_RATIO", "MINOR_ALLELE_FRACTION_POSTERIOR_90", "OVERLAP", "CNV_CALL", "BAF_CALL")
  
  ##..analyse copy neutral segments..
  #cn_flag = !((all_res[,11]%in%"NEUTRAL") & (all_res[,12]%in%"HET" | all_res[,12]%in%"NA"))
  #cn_res = all_res[cn_flag,]
  cn_res = all_res
  
  # remove neutral
  index = which(cn_res[,"CNV_CALL"]=="NEUTRAL")
  if (length(index)!=0)	cn_res = cn_res[-which(cn_res[,"CNV_CALL"]=="NEUTRAL"),,drop=FALSE]

  # check if any left
  if (dim(cn_res)[1]==0){
    message("Only neutral CNVs called. There will be no annotation file!")
    system(paste0("touch ", out_file))
  } else {
    
    # remove strand
    cn_res = cn_res[,-which(colnames(cn_res)=="GENE_STRAND"),drop=FALSE]
    
    # CNV length
    cn_res = cbind(cn_res, "CNV_LENGTH"=as.numeric(cn_res[,"CNV_STOP"])-as.numeric(cn_res[,"CNV_START"]))
    
    #..replace dots by commas in decimal numbers for excel compatibility..
    # cn_res[, c(8,9)] = gsub("\\.", ",", cn_res[, c(8,9)])
    
    # Overlapping genes
    #all_genes <- read.delim(paste0(home_dir, "Homo_sapiens.GRCh38.98_protein_coding_genes.bed"), header = FALSE)
    all_genes <- read.delim(paste0(protein_coding_genes_bed), header = FALSE)
    
    gene_ranges <- GRanges(
      seqnames = Rle(all_genes[,1]),
      ranges = IRanges(as.numeric(all_genes[,2])+1, end = as.numeric(as.numeric(all_genes[,3])+1)),
      strand = Rle(strand(rep("+",dim(all_genes)[1]))))
    
    values(gene_ranges) <- DataFrame("gene_name" = as.vector(all_genes[,4]))
    
    genes = NULL
    number = NULL
    for (i in 1:dim(cn_res)[1]){
      cnv <- GRanges(seqnames=cn_res[i,"CHR"],
                  ranges=IRanges(start = as.numeric(cn_res[i,"CNV_START"]),
                                 end = as.numeric(cn_res[i,"CNV_STOP"])),
                  strand="+")
      
      overlap <- intersect(gene_ranges, cnv)
      gene_ranges.overlap <- subsetByOverlaps(gene_ranges,overlap)
  
      if (length(gene_ranges.overlap)!=0){
          genes = c(genes, paste(mcols(gene_ranges.overlap)$gene_name, collapse = ";"))
          number = c(number, length(mcols(gene_ranges.overlap)$gene_name))
      } else {
          genes = c(genes, NA)
          number = c(number, 0)
      }
    }
    
    names = c(colnames(cn_res),"OVERLAPPING_GENES_#","OVERLAPPING_GENES")
    cn_res = cbind(cn_res, number, genes)
    colnames(cn_res) = names
    
    #..save..
    #sname = gsub(".modelFinal.seg", ".modelFinal.filteredOnPanel.txt", input_file)
    sname = out_file
    write.table(cn_res, file=sname, col.names=TRUE, row.names=FALSE, quote=FALSE, sep = "\t")
  }  
}

message("Start")
message("Reading arguments")

cmd_args = commandArgs(trailingOnly = TRUE)
for (arg in cmd_args) cat("  ", arg, "\n", sep="")
#run(cmd_args[1], cmd_args[2])
run(cmd_args[1], cmd_args[2], cmd_args[3], cmd_args[4], cmd_args[5])

message("Ready")
