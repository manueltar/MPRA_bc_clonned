

suppressMessages(library("plyr", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("data.table", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("crayon", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("withr", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("ggplot2", lib.loc = "/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("optparse", lib.loc = "/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("dplyr", lib.loc = "/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("withr", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("backports", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("broom", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("rstudioapi", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("tzdb", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("cli", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("tidyverse", lib.loc="/nfs/team151/software/manuel_R_libs_4_1//"))

suppressMessages(library("R.methodsS3", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("R.oo", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("R.utils", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
# suppressMessages(library("splitstackshape", lib.loc = "/nfs/users/nfs_m/mt19/sOFTWARE/R_libs"))



opt = NULL

READ_FILES = function(option_list)
{
  
  
  opt_in = option_list
  opt <<- option_list
  
  cat("All options:\n")
  printList(opt)
  
  #### READ and transform type ----
  
  type = opt$type
  
  cat("TYPE_\n")
  cat(sprintf(as.character(type)))
  cat("\n")
  
  #### READ and transform out ----
  
  out = opt$out
  
  cat("OUT_\n")
  cat(sprintf(as.character(out)))
  cat("\n")
  
  #### READ input_sequencing_seqname_collapsed ----
  
  
  input_sequencing_seqname_collapsed<-as.data.frame(fread(opt$input_sequencing_seqname_collapsed, sep=" ", header = T, fill=TRUE), stringsAsFactors = F)
  
  colnames(input_sequencing_seqname_collapsed)<-c("number_of_barcodes","seq_name")
  cat("input_sequencing_seqname_collapsed_0\n")
  cat(str(input_sequencing_seqname_collapsed))
  cat("\n")
  
  
  input_sequencing_seqname_collapsed$REAL_TILE<-gsub("\\|.+","",input_sequencing_seqname_collapsed$seq_name, perl=T)
  
  input_sequencing_seqname_collapsed$carried_variants<-gsub("^[^\\|]+\\|","",input_sequencing_seqname_collapsed$seq_name, perl=T)
  input_sequencing_seqname_collapsed$carried_variants<-gsub("\\|.+","",input_sequencing_seqname_collapsed$carried_variants, perl=T)
  
  input_sequencing_seqname_collapsed$Flag_intra_candidate_pattern<-gsub("^[^\\|]+\\|[^\\|]+\\|","",input_sequencing_seqname_collapsed$seq_name, perl=T)
  
  #colnames(input_sequencing_seqname_collapsed)[which(colnames(input_sequencing_seqname_collapsed) == "V2")]<-"Barcode"
  
  
  cat("input_sequencing_seqname_collapsed_\n")
  cat(str(input_sequencing_seqname_collapsed))
  cat("\n")
  
  #### input_sequencing_Barcode_collapsed----
  
  input_sequencing_Barcode_collapsed<-as.data.frame(fread(opt$input_sequencing_Barcode_collapsed, 
                                                               sep="\t", header = T, fill=TRUE), stringsAsFactors = F)
  
  
  colnames(input_sequencing_Barcode_collapsed)<-c("string","Barcode")
  
  input_sequencing_Barcode_collapsed$number_of_reads<-as.integer(gsub(" .+","",input_sequencing_Barcode_collapsed$string))
  input_sequencing_Barcode_collapsed$seq_name<-gsub("^[^ ] ","",input_sequencing_Barcode_collapsed$string)
  
  indx.dep<-which(colnames(input_sequencing_Barcode_collapsed) == "string")
  
  input_sequencing_Barcode_collapsed<-input_sequencing_Barcode_collapsed[,-indx.dep]
  
  input_sequencing_Barcode_collapsed$REAL_TILE<-gsub("\\|.+","",input_sequencing_Barcode_collapsed$seq_name, perl=T)
  
  input_sequencing_Barcode_collapsed$carried_variants<-gsub("^[^\\|]+\\|","",input_sequencing_Barcode_collapsed$seq_name, perl=T)
  input_sequencing_Barcode_collapsed$carried_variants<-gsub("\\|.+","",input_sequencing_Barcode_collapsed$carried_variants, perl=T)
  
  input_sequencing_Barcode_collapsed$Flag_intra_candidate_pattern<-gsub("^[^\\|]+\\|[^\\|]+\\|","",input_sequencing_Barcode_collapsed$seq_name, perl=T)
  
  input_sequencing_Barcode_collapsed$nchar_Barcode<-nchar(input_sequencing_Barcode_collapsed$Barcode)
  
  
  cat("input_sequencing_Barcode_collapsed_0\n")
  cat(str(input_sequencing_Barcode_collapsed))
  cat("\n")
  
  #### READ MPRA_Rosetta ----
  
  
  MPRA_Rosetta<-as.data.frame(fread(opt$MPRA_Rosetta, sep="\t", header = T), stringsAsFactors = F)
  
  cat("MPRA_Rosetta_\n")
  cat(str(MPRA_Rosetta))
  cat("\n")
  
  #### READ Blacklist ----
  
  
  Blacklist<-as.data.frame(fread(opt$Blacklist, sep="\t", header = T), stringsAsFactors = F)
  
  colnames(Blacklist)<-"Barcode"
  
  cat("Blacklist_\n")
  cat(str(Blacklist))
  cat("\n")
  
  #quit(status=1)
  #### READ and transform Threshold_min_length_barcode ----
  
  Threshold_min_length_barcode = opt$Threshold_min_length_barcode
  
  cat("Threshold_min_length_barcode_\n")
  cat(sprintf(as.character(Threshold_min_length_barcode)))
  cat("\n")
  
 
  #### READ and transform Threshold_reads_per_bc ----
  
  Threshold_reads_per_bc = opt$Threshold_reads_per_bc
  
  cat("Threshold_reads_per_bc_\n")
  cat(sprintf(as.character(Threshold_reads_per_bc)))
  cat("\n")
  
  #### READ and transform Threshold_min_bc_per_tile ----
  
  Threshold_min_bc_per_tile = opt$Threshold_min_bc_per_tile
  
  cat("Threshold_min_bc_per_tile_\n")
  cat(sprintf(as.character(Threshold_min_bc_per_tile)))
  cat("\n")
  
  #### sequences in Rosetta ----
  
  input_sequencing_Barcode_collapsed_IN<-input_sequencing_Barcode_collapsed[which(input_sequencing_Barcode_collapsed$REAL_TILE%in%MPRA_Rosetta$REAL_TILE &
                                                input_sequencing_Barcode_collapsed$carried_variants%in%MPRA_Rosetta$carried_variants),]
  
  cat("input_sequencing_Barcode_collapsed_IN_\n")
  cat(str(input_sequencing_Barcode_collapsed_IN))
  cat("\n")
  
  
  
  input_sequencing_Barcode_collapsed_OUT<-input_sequencing_Barcode_collapsed[-which(input_sequencing_Barcode_collapsed$sequence%in%MPRA_Rosetta$sequence),]
  
  cat("input_sequencing_Barcode_collapsed_OUT_\n")
  cat(str(input_sequencing_Barcode_collapsed_OUT))
  cat("\n")
  
  # quit(status=1)
  
  #### Filter 0 Threshold_reads_per_bc ----
  
  input_sequencing_Barcode_collapsed_IN_Threshold_0<-input_sequencing_Barcode_collapsed_IN[which(input_sequencing_Barcode_collapsed_IN$nchar_Barcode >= Threshold_min_length_barcode),]
  
  cat("input_sequencing_Barcode_collapsed_IN_Threshold_0_\n")
  cat(str(input_sequencing_Barcode_collapsed_IN_Threshold_0))
  cat("\n")
  cat(sprintf(as.character(names(summary(input_sequencing_Barcode_collapsed_IN_Threshold_0$number_of_reads)))))
  cat("\n")
  cat(sprintf(as.character(summary(input_sequencing_Barcode_collapsed_IN_Threshold_0$number_of_reads))))
  cat("\n")
  
  
  #### Filter 1 Threshold_reads_per_bc ----
  
  input_sequencing_Barcode_collapsed_IN_Threshold_1<-input_sequencing_Barcode_collapsed_IN_Threshold_0[which(input_sequencing_Barcode_collapsed_IN_Threshold_0$number_of_reads >= Threshold_reads_per_bc),]
  
  cat("input_sequencing_Barcode_collapsed_IN_Threshold_1_\n")
  cat(str(input_sequencing_Barcode_collapsed_IN_Threshold_1))
  cat("\n")
  
  #### Filter 2 NO Barcode collisions in the FOUND sequences----
  
  
  Freq.table_barcodes<-as.data.frame(table(input_sequencing_Barcode_collapsed_IN_Threshold_1$Barcode), stringsAsFactors=F)
  
  colnames(Freq.table_barcodes)[which(colnames(Freq.table_barcodes) == "Var1")]<-"Barcode"
  
  
  cat("Freq.table_barcodes\n")
  cat(str(Freq.table_barcodes))
  cat("\n")
  
  Freq.table_barcodes_unique<-Freq.table_barcodes[which(Freq.table_barcodes$Freq == 1),]
  
  cat("Freq.table_barcodes_unique\n")
  cat(str(Freq.table_barcodes_unique))
  cat("\n")
  
  input_sequencing_Barcode_collapsed_IN_Threshold_2<-input_sequencing_Barcode_collapsed_IN_Threshold_1[which(input_sequencing_Barcode_collapsed_IN_Threshold_1$Barcode%in%Freq.table_barcodes_unique$Barcode),]
  
  cat("input_sequencing_Barcode_collapsed_IN_Threshold_2_\n")
  cat(str(input_sequencing_Barcode_collapsed_IN_Threshold_2))
  cat("\n")
  
  #### Filter 3 NO Barcode collisions with Filtered out sequences----
  
  
  input_sequencing_Barcode_collapsed_IN_Threshold_3<-input_sequencing_Barcode_collapsed_IN_Threshold_2[-which(input_sequencing_Barcode_collapsed_IN_Threshold_2$Barcode%in%Blacklist$Barcode),]
  

  cat("input_sequencing_Barcode_collapsed_IN_Threshold_3_\n")
  cat(str(input_sequencing_Barcode_collapsed_IN_Threshold_3))
  cat("\n")
  
  #### Filter 4 Range of barcodes per REAL TILE ----
  
  
  Freq.table_range<-as.data.frame(table(input_sequencing_Barcode_collapsed_IN_Threshold_3$seq_name), stringsAsFactors=F)
  
  colnames(Freq.table_range)[which(colnames(Freq.table_range) == "Var1")]<-"seq_name"
  colnames(Freq.table_range)[which(colnames(Freq.table_range) == "Freq")]<-"number_of_barcodes"
  
  
  cat("Freq.table_range\n")
  cat(str(Freq.table_range))
  cat("\n")
  cat(sprintf(as.character(names(summary(Freq.table_range$number_of_barcodes)))))
  cat("\n")
  cat(sprintf(as.character(summary(Freq.table_range$number_of_barcodes))))
  cat("\n")
  
  
  Freq.table_range_PASS<-Freq.table_range[which(Freq.table_range$number_of_barcodes >= Threshold_min_bc_per_tile),]
  
  cat("Freq.table_range_PASS\n")
  cat(str(Freq.table_range_PASS))
  cat("\n")
  
  input_sequencing_Barcode_collapsed_IN_Threshold_3<-input_sequencing_Barcode_collapsed_IN_Threshold_3[which(input_sequencing_Barcode_collapsed_IN_Threshold_3$seq_name%in%Freq.table_range_PASS$seq_name),]
  cat("input_sequencing_Barcode_collapsed_IN_Threshold_3_\n")
  cat(str(input_sequencing_Barcode_collapsed_IN_Threshold_3))
  cat("\n")
  
  # quit(status=1)
  
  #### merge with Rosetta ----
  
  DEF<-merge(input_sequencing_Barcode_collapsed_IN_Threshold_3,
                                         MPRA_Rosetta,
                                         by=c("REAL_TILE","carried_variants","Flag_intra_candidate_pattern","seq_name"),
                                         all.x=T)
  
  
  cat("DEF\n")
  cat(str(DEF))
  cat("\n")
  
  DEF_subset<-unique(DEF[,c(which(colnames(DEF) == "seq_name"),which(colnames(DEF) == "factor4"))])
  
  cat("DEF_subset\n")
  cat(str(DEF_subset))
  cat("\n")
  
  cat(sprintf(as.character(names(summary(as.factor(DEF_subset$factor4))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(DEF_subset$factor4)))))
  cat("\n")
  
 # quit(status=1)
  
  # #### SAVE TILE DEF ----
  # 
  setwd(out)

  filename18<-paste(type,".rds",sep='')

  cat("\n")
  cat(sprintf(as.character(filename18)))
  cat("\n")

  saveRDS(DEF,file=filename18)

  #quit(status = 1)

}



printList = function(l, prefix = "    ") {
  list.df = data.frame(val_name = names(l), value = as.character(l))
  list_strs = apply(list.df, MARGIN = 1, FUN = function(x) { paste(x, collapse = " = ")})
  cat(paste(paste(paste0(prefix, list_strs), collapse = "\n"), "\n"))
}
#### main script ----

main = function() {
  cmd_line = commandArgs()
  cat("Command line:\n")
  cat(paste(gsub("--file=", "", cmd_line[4], fixed=T),
            paste(cmd_line[6:length(cmd_line)], collapse = " "),
            "\n\n"))
  option_list <- list(
    make_option(c("--input_sequencing_seqname_collapsed"), type="character", default=NULL,
                metavar="FILE.txt",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--input_sequencing_Barcode_collapsed"), type="character", default=NULL,
                metavar="FILE.txt",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Blacklist"), type="character", default=NULL,
                metavar="FILE.txt",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--MPRA_Rosetta"), type="character", default=NULL,
                metavar="FILE.txt",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Threshold_min_length_barcode"), type="numeric", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Threshold_reads_per_bc"), type="numeric", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Threshold_min_bc_per_tile"), type="numeric", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
       make_option(c("--type"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--out"), type="character", default=NULL, 
                metavar="filename", 
                help="Path to tab-separated input file listing regions to analyze. Required.")
    
  )
  parser = OptionParser(usage = "140__Rscript_v106.R
                        --subset type
                        --PCHIC_original FILE.txt
                        --GeneEXP FILE.txt
                        --type type
                        --out filename",
                        option_list = option_list)
  opt <<- parse_args(parser)
  
  
  READ_FILES(opt)

   
}

###########################################################################

system.time( main() )

