

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
  
  #### READ TOTAL_ALIGNED_readIDs ----
  
  
  TOTAL_ALIGNED_readIDs<-as.data.frame(fread(opt$TOTAL_ALIGNED_readIDs, sep=" ", header = F, fill=TRUE), stringsAsFactors = F)
  
  TOTAL_ALIGNED_readIDs$readID<-gsub("_.+$","",TOTAL_ALIGNED_readIDs$V1)
  TOTAL_ALIGNED_readIDs$Barcode<-gsub("^[^_]+_","",TOTAL_ALIGNED_readIDs$V1)
  
 
  cat("TOTAL_ALIGNED_readIDs_0\n")
  cat(str(TOTAL_ALIGNED_readIDs))
  cat("\n")
  
  
  #### FILTERED_ALIGNED_readIDs----
  
  FILTERED_ALIGNED_readIDs<-as.data.frame(fread(opt$FILTERED_ALIGNED_readIDs, 
                                                               sep="\t", header = F, fill=TRUE), stringsAsFactors = F)
  
  
  FILTERED_ALIGNED_readIDs$readID<-gsub("_.+$","",FILTERED_ALIGNED_readIDs$V1)
  FILTERED_ALIGNED_readIDs$Barcode<-gsub("^[^_]+_","",FILTERED_ALIGNED_readIDs$V1)
  
  cat("FILTERED_ALIGNED_readIDs_0\n")
  cat(str(FILTERED_ALIGNED_readIDs))
  cat("\n")
  
  #### FIND FILTERED OUT ----
  
  FILTERED_OUT<-TOTAL_ALIGNED_readIDs[-which(TOTAL_ALIGNED_readIDs$readID%in%FILTERED_ALIGNED_readIDs$readID &
                                               TOTAL_ALIGNED_readIDs$Barcode%in%FILTERED_ALIGNED_readIDs$Barcode),]
  
  cat("FILTERED_OUT\n")
  cat(str(FILTERED_OUT))
  cat("\n")
  
  if(dim(FILTERED_OUT)[1]>0)
  {
    Blacklist<-unique(FILTERED_OUT$Barcode)
    
    cat("Blacklist\n")
    cat(str(Blacklist))
    cat("\n")
    
   
    
  }else{
    
    Blacklist<-"NA"
    
    
  }
  
  # #### SAVE TILE DEF ----
  # 
  setwd(out)
  
  filename18<-paste(type,".txt",sep='')
  
  cat("\n")
  cat(sprintf(as.character(filename18)))
  cat("\n")
  
  write.table(Blacklist, sep="\t", quote =F, row.names=F, file=filename18)

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
    make_option(c("--TOTAL_ALIGNED_readIDs"), type="character", default=NULL,
                metavar="FILE.txt",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--FILTERED_ALIGNED_readIDs"), type="character", default=NULL,
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

