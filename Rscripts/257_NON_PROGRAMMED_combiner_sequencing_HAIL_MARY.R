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


suppressMessages(library("BiocGenerics", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("S4Vectors", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("IRanges", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("GenomeInfoDb", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("GenomicRanges", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("Biobase", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("AnnotationDbi", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("GenomicFeatures", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("OrganismDbi", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
library("XVector", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")
suppressMessages(library("Biostrings", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("seqinr", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
#suppressMessages(library("seqRFLP", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("splitstackshape", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))

opt = NULL

ALL_TOGETHER = function(option_list)
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
  
  #### READ MPRA_Rosetta ----
  
  
  MPRA_Rosetta<-as.data.frame(fread(opt$MPRA_Rosetta, sep="\t", header = T), stringsAsFactors = F)
  
  cat("MPRA_Rosetta_\n")
  cat(str(MPRA_Rosetta))
  cat("\n")
  
  
  # quit(status = 1)
  
  #### Fasta_Reference ----
  
  fastaFile<-readDNAStringSet(file=opt$Fasta_Reference_1)
  
  # cat("fastaFile\n")
  # cat(str(fastaFile))
  # cat("\n")
  
  seq_name = names(fastaFile)
  sequence = paste(fastaFile)
  Fasta_Reference_1 <- data.frame(seq_name, sequence, stringsAsFactors = F)
  
  Fasta_Reference_1$REAL_TILE<-gsub("\\|.+$","",Fasta_Reference_1$seq_name)
  
  cat("Fasta_Reference_1\n")
  cat(str(Fasta_Reference_1))
  cat("\n")
  
  MPRA_Rosetta_Fasta_Reference_1<-MPRA_Rosetta[which(MPRA_Rosetta$REAL_TILE%in%Fasta_Reference_1$REAL_TILE),]
  
  cat("MPRA_Rosetta_Fasta_Reference_1\n")
  cat(str(MPRA_Rosetta_Fasta_Reference_1))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(MPRA_Rosetta_Fasta_Reference_1$factor4))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(MPRA_Rosetta_Fasta_Reference_1$factor4)))))
  cat("\n")
  
  ####
  

  
  fastaFile<-readDNAStringSet(file=opt$Fasta_Reference_2)
  
  # cat("fastaFile\n")
  # cat(str(fastaFile))
  # cat("\n")
  
  seq_name = names(fastaFile)
  sequence = paste(fastaFile)
  Fasta_Reference_2 <- data.frame(seq_name, sequence, stringsAsFactors = F)
  
  Fasta_Reference_2$REAL_TILE<-gsub("\\|.+$","",Fasta_Reference_2$seq_name)
  
  cat("Fasta_Reference_2\n")
  cat(str(Fasta_Reference_2))
  cat("\n")
  
  MPRA_Rosetta_Fasta_Reference_2<-MPRA_Rosetta[which(MPRA_Rosetta$REAL_TILE%in%Fasta_Reference_2$REAL_TILE),]
  
  cat("MPRA_Rosetta_Fasta_Reference_2\n")
  cat(str(MPRA_Rosetta_Fasta_Reference_2))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(MPRA_Rosetta_Fasta_Reference_2$factor4))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(MPRA_Rosetta_Fasta_Reference_2$factor4)))))
  cat("\n")
  
  ####
  
  
  fastaFile<-readDNAStringSet(file=opt$Fasta_Reference_3)
  
  # cat("fastaFile\n")
  # cat(str(fastaFile))
  # cat("\n")
  
  seq_name = names(fastaFile)
  sequence = paste(fastaFile)
  Fasta_Reference_3 <- data.frame(seq_name, sequence, stringsAsFactors = F)
  
  Fasta_Reference_3$REAL_TILE<-gsub("\\|.+$","",Fasta_Reference_3$seq_name)
  
  cat("Fasta_Reference_3\n")
  cat(str(Fasta_Reference_3))
  cat("\n")
  
  MPRA_Rosetta_Fasta_Reference_3<-MPRA_Rosetta[which(MPRA_Rosetta$REAL_TILE%in%Fasta_Reference_3$REAL_TILE),]
  
  cat("MPRA_Rosetta_Fasta_Reference_3\n")
  cat(str(MPRA_Rosetta_Fasta_Reference_3))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(MPRA_Rosetta_Fasta_Reference_3$factor4))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(MPRA_Rosetta_Fasta_Reference_3$factor4)))))
  cat("\n")
  
  ####
  
  ##### OVERLAPS ----
  
  
  OV_seq_name_1_2<-Fasta_Reference_1$seq_name[which(Fasta_Reference_1$seq_name%in%Fasta_Reference_2$seq_name)]
  
  cat("OV_seq_name_1_2\n")
  cat(str(OV_seq_name_1_2))
  cat("\n")
  
  
  OV_seq_name_1_3<-Fasta_Reference_1$seq_name[which(Fasta_Reference_1$seq_name%in%Fasta_Reference_3$seq_name)]
  
  cat("OV_seq_name_1_3\n")
  cat(str(OV_seq_name_1_3))
  cat("\n")
  
  OV_seq_name_2_3<-Fasta_Reference_2$seq_name[which(Fasta_Reference_2$seq_name%in%Fasta_Reference_3$seq_name)]
  
  cat("OV_seq_name_2_3\n")
  cat(str(OV_seq_name_2_3))
  cat("\n")
  
  
  ###
  
  
  OV_sequence_1_2<-Fasta_Reference_1$sequence[which(Fasta_Reference_1$sequence%in%Fasta_Reference_2$sequence)]
  
  cat("OV_sequence_1_2\n")
  cat(str(OV_sequence_1_2))
  cat("\n")
  
  
  OV_sequence_1_3<-Fasta_Reference_1$sequence[which(Fasta_Reference_1$sequence%in%Fasta_Reference_3$sequence)]
  
  cat("OV_sequence_1_3\n")
  cat(str(OV_sequence_1_3))
  cat("\n")
  
  OV_sequence_2_3<-Fasta_Reference_2$sequence[which(Fasta_Reference_2$sequence%in%Fasta_Reference_3$sequence)]
  
  cat("OV_sequence_2_3\n")
  cat(str(OV_sequence_2_3))
  cat("\n")
  
  
  
  
  
  ####REF DEF ----
  
  
  DEF_reference<-unique(rbind(Fasta_Reference_1,Fasta_Reference_2,Fasta_Reference_3))
  
  cat("DEF_reference\n")
  cat(str(DEF_reference))
  cat("\n")
  cat(str(unique(DEF_reference$seq_name)))
  cat("\n")
  cat(str(unique(DEF_reference$sequence)))
  cat("\n")
  
  
  
  DEF_reference.dt<-data.table(DEF_reference,
                                         key="sequence")
  
  
  cat("DEF_reference.dt\n")
  cat(str(DEF_reference.dt))
  cat("\n")
  
  Freq_table_Quan<-as.data.frame(DEF_reference.dt[, .N,
                                                            by=key(DEF_reference.dt)]
                                 ,stringsAsFactors=F)
  
  
  
  
  
  colnames(Freq_table_Quan)[which(colnames(Freq_table_Quan) == "N")]<-"Repetitions"
  
  
  cat("Freq_table_Quan\n")
  cat(str(Freq_table_Quan))
  cat("\n")
  
  
  Freq_table_Quan_unique_sequences<-Freq_table_Quan[which(Freq_table_Quan$Repetitions == 1),]
  
  cat("Freq_table_Quan_unique_sequences\n")
  cat(str(Freq_table_Quan_unique_sequences))
  cat("\n")
  
  
  DEF_reference_DEF<-DEF_reference[which(DEF_reference$sequence%in%Freq_table_Quan_unique_sequences$sequence),]
  
  cat("DEF_reference_DEF\n")
  cat(str(DEF_reference_DEF))
  cat("\n")
  cat(str(unique(DEF_reference_DEF$seq_name)))
  cat("\n")
  cat(str(unique(DEF_reference_DEF$sequence)))
  cat("\n")
  
  #### merge with MPRA_Rosetta ----  
  
  
  # DEF_reference_DEF<-merge(DEF_reference_DEF,
  #                          MPRA_Rosetta,
  #                          by="REAL_TILE",
  #                          all=T)
  # 
  # cat("DEF_reference_DEF\n")
  # cat(str(DEF_reference_DEF))
  # cat("\n")
  # cat(str(unique(!is.na(DEF_reference_DEF$seq_name))))
  # cat("\n")
  # cat(str(unique(!is.na(DEF_reference_DEF$sequence))))
  # cat("\n")
  # 
  
  REPRESENTED_sequences<-MPRA_Rosetta[which(MPRA_Rosetta$REAL_TILE%in%DEF_reference_DEF$REAL_TILE),]
  
  
  cat("REPRESENTED_sequences\n")
  cat(str(REPRESENTED_sequences))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(REPRESENTED_sequences$factor4))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(REPRESENTED_sequences$factor4)))))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(REPRESENTED_sequences$Label))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(REPRESENTED_sequences$Label)))))
  cat("\n")
  
  UNREPRESENTED_sequences<-MPRA_Rosetta[-which(MPRA_Rosetta$REAL_TILE%in%DEF_reference_DEF$REAL_TILE),]
  
  
  cat("UNREPRESENTED_sequences\n")
  cat(str(UNREPRESENTED_sequences))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(UNREPRESENTED_sequences$factor4))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(UNREPRESENTED_sequences$factor4)))))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(UNREPRESENTED_sequences$Label))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(UNREPRESENTED_sequences$Label)))))
  cat("\n")
  
  #### PRINT FASTA ----
  
  
  indx.fasta<-c(which(colnames(DEF_reference_DEF) == "seq_name"),which(colnames(DEF_reference_DEF) == "sequence"))
  
  
  
  Merge2_TO_FASTA<-unique(DEF_reference_DEF[,indx.fasta])
  
  cat("Merge2_TO_FASTA\n")
  str(Merge2_TO_FASTA)
  cat("\n")
  
  
  
  sequences<-list()
  
  for(i in 1:length(Merge2_TO_FASTA$sequence))
  {
    sequences[[i]]<-Merge2_TO_FASTA$sequence[i]
    
  }#i
  
  
  cat("sequences\n")
  str(sequences)
  cat("\n")
  
  names<-list()
  
  for(i in 1:length(Merge2_TO_FASTA$seq_name))
  {
    names[[i]]<-Merge2_TO_FASTA$seq_name[i]
    
  }#i
  
  
  
  cat("names\n")
  str(names)
  cat("\n")
  
  setwd(out)
  
  filename<-paste("combined_reference",'.fasta',sep='')
  
  
  write.fasta(sequences, names, filename, open = "w", nbchar = 60, as.string = FALSE)
  
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
    make_option(c("--Fasta_Reference_1"), type="character", default=NULL,
                metavar="FILE.txt",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Fasta_Reference_2"), type="character", default=NULL,
                metavar="FILE.txt",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Fasta_Reference_3"), type="character", default=NULL,
                metavar="FILE.txt",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--MPRA_Rosetta"), type="character", default=NULL,
                metavar="FILE.txt",
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
                        --TranscriptEXP FILE.txt
                        --cadd FILE.txt
                        --ncboost FILE.txt
                        --type type
                        --out filename",
                        option_list = option_list)
  opt <<- parse_args(parser)
  
  ALL_TOGETHER(opt)

  
}


###########################################################################

system.time( main() )
