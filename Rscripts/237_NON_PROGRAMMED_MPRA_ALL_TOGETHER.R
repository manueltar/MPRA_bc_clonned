
suppressMessages(library("data.table", lib.loc="/nfs/users/nfs_m/mt19/sOFTWARE/R_libs/"))
suppressMessages(library("ggplot2", lib.loc = "/nfs/users/nfs_m/mt19/sOFTWARE/R_libs"))
suppressMessages(library("optparse", lib.loc = "/nfs/users/nfs_m/mt19/sOFTWARE/R_libs"))
suppressMessages(library("plyr", lib.loc = "/nfs/users/nfs_m/mt19/sOFTWARE/R_libs/"))
suppressMessages(library("dplyr", lib.loc = "/nfs/users/nfs_m/mt19/sOFTWARE/R_libs/"))
suppressMessages(library("tidyverse", lib.loc = "/nfs/users/nfs_m/mt19/sOFTWARE/R_libs/"))
suppressMessages(library("svglite", lib.loc="/nfs/users/nfs_m/mt19/sOFTWARE/R_libs/"))
suppressMessages(library("ggforce", lib.loc = "/nfs/users/nfs_m/mt19/sOFTWARE/R_libs/"))
suppressMessages(library("ggeasy", lib.loc = "/nfs/users/nfs_m/mt19/sOFTWARE/R_libs/"))
suppressMessages(library("sandwich", lib.loc = "/nfs/users/nfs_m/mt19/sOFTWARE/R_libs/"))


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
  
  
 # quit(status=1)
  
  setwd(out)
  
  file_list <- list.files(path=out, include.dirs = FALSE)
  
  cat("file_list\n")
  cat(str(file_list))
  cat("\n")
  
  indexes_sel <- grep("\\.rds$",file_list)
  
  cat("indexes_sel\n")
  cat(sprintf(as.character(indexes_sel)))
  cat("\n")
  
  file_list_sel <- as.data.frame(file_list[indexes_sel], stringsAsFactors=F)
  
  colnames(file_list_sel)<-"file"
  file_list_sel$lane<-gsub('.rds$','',file_list_sel$file)
  file_list_sel$subsampling<-gsub("^lane_[0-9]+_","",file_list_sel$lane)
  file_list_sel$lane<-gsub("_$","",gsub(paste(file_list_sel$subsampling, collapse="|"),"",file_list_sel$lane))
  
  
  
  cat("file_list_sel\n")
  cat(str(file_list_sel))
  cat("\n")
  
  
  
  lanes<-c("lane_1","lane_2")
  subsampling<-rev(c("100_percent","0.9","0.75","0.5","0.25","0.1","0.01"))

  file_list_sel$lane<-factor(file_list_sel$lane,
                             levels=lanes,
                             ordered=T)
  
  file_list_sel$subsampling<-factor(file_list_sel$subsampling,
                             levels=subsampling,
                             ordered=T)
  
  
  
  
  cat("file_list_sel\n")
  cat(str(file_list_sel))
  cat("\n")
  
  # quit(status=1)
  

  
  #### LOOP READ RDS files ----
  
  List_lane<-list()
  
  
  for(i in 1:length(levels(file_list_sel$lane)))
  {
    
    lane_sel<-levels(file_list_sel$lane)[i]
    
    cat("lane_sel\n")
    cat(sprintf(as.character(lane_sel)))
    cat("\n")
    
    # #### path2 ----
    # 
    # path2<-paste(out,lane_sel, sep='')
    # 
    # cat("path2\n")
    # cat(sprintf(as.character(path2)))
    # cat("\n")
    # 
    # if (file.exists(path2)){
    #   
    #   unlink(path2, recursive = TRUE)
    #   dir.create(file.path(path2))
    #   
    #   
    # } else {
    #   dir.create(file.path(path2))
    #   
    # }
    
    
   
    file_list_subset<-file_list_sel[which(file_list_sel$lane%in%lane_sel),]
    
    cat("file_list_subset\n")
    cat(str(file_list_subset))
    cat("\n")
    
    
    
    if(dim(file_list_subset)[1]>0)
    {
      
     # quit(status = 1)
      
      List_RESULTS<-list()
      
      for(k in 1:dim(file_list_subset)[1])
      {
        subsampling_sel<-file_list_subset$subsampling[k]
        
        cat("subsampling_sel\n")
        cat(sprintf(as.character(subsampling_sel)))
        cat("\n")
        
        file_sel<-file_list_subset$file[file_list_subset$subsampling == subsampling_sel]
        
        cat("file_sel\n")
        cat(sprintf(as.character(file_sel)))
        cat("\n")
        
        if(length(file_sel) > 0)
        {
          setwd(out)
          
          df<-readRDS(file = file_sel)
          
          cat("df\n")
          cat(str(df))
          cat("\n")
          
          if(dim(df)[1] >0)
          {
            df$lane<-lane_sel
            df$subsampling<-subsampling_sel
            
            
            
            List_RESULTS[[k]]<-df
            
          }
        
          
          
          # quit(status=1)
          
        }# length(file_sel) > 0
      
        
      } #k subsampling_sel
      
      
      Lane_TABLE = as.data.frame(data.table::rbindlist(List_RESULTS, fill = T))

      # cat("Lane_TABLE\n")
      # str(Lane_TABLE)
      # cat("\n")
      # 
      # quit(status = 1)
      List_lane[[i]]<-Lane_TABLE
      
  }# dim(file_list_subset)[1]>0
   
    
  }  #i
  
  
  ####
  
  TABLE_DEF = as.data.frame(data.table::rbindlist(List_lane, fill = T))
  
  cat("TABLE_DEF\n")
  str(TABLE_DEF)
  cat("\n")
  
  
  TABLE_DEF_Total<-as.data.frame(setDT(TABLE_DEF)[, .(.N,paste(Barcode, collapse = "|"),paste(number_of_reads, collapse = "|")), by=.(lane,subsampling,seq_name)], stringsAsFactors=F)
  
  colnames(TABLE_DEF_Total)[which(colnames(TABLE_DEF_Total) == "N")]<-"Tagging_barcodes"
  colnames(TABLE_DEF_Total)[which(colnames(TABLE_DEF_Total) == "V2")]<-"Barcode_string"
  colnames(TABLE_DEF_Total)[which(colnames(TABLE_DEF_Total) == "V3")]<-"number_of_reads_string"
  
  
  
  cat("TABLE_DEF_Total\n")
  str(TABLE_DEF_Total)
  cat("\n")
  
  
  # quit(status=1)
  
  list_DEF<-list()  
  
  #### POOL A Freq table ----
  
  
  Rosetta_POOL_A<-MPRA_Rosetta[which(MPRA_Rosetta$factor4 == "High_AT"),]
  
  # cat("Rosetta_POOL_A\n")
  # str(Rosetta_POOL_A)
  # cat("\n")
  # 
  TABLE_DEF_Total_POOL_A<-TABLE_DEF_Total[which(TABLE_DEF_Total$seq_name%in%Rosetta_POOL_A$seq_name),]
  
  if(dim(TABLE_DEF_Total_POOL_A)[1]>0)
  {
    TABLE_DEF_Total_POOL_A$pool<-"POOL_A"
    
    
    cat("TABLE_DEF_Total_POOL_A\n")
    str(TABLE_DEF_Total_POOL_A)
    cat("\n")
    
    TABLE_DEF_Total_POOL_A_CAPTURED<-as.data.frame(setDT(TABLE_DEF_Total_POOL_A)[, .N, by=.(lane,subsampling,pool)], stringsAsFactors=F)
    
    colnames(TABLE_DEF_Total_POOL_A_CAPTURED)[which(colnames(TABLE_DEF_Total_POOL_A_CAPTURED) == "N")]<-"CAPTURED"
    
    TABLE_DEF_Total_POOL_A_CAPTURED$Perc_captured<-100*(TABLE_DEF_Total_POOL_A_CAPTURED$CAPTURED/dim(Rosetta_POOL_A)[1])
    
    
    
    cat("TABLE_DEF_Total_POOL_A_CAPTURED\n")
    str(TABLE_DEF_Total_POOL_A_CAPTURED)
    cat("\n")
    
    
    
    list_DEF$tagging_barcodes_df[[1]]<-TABLE_DEF_Total_POOL_A
    list_DEF$summary_df[[1]]<-TABLE_DEF_Total_POOL_A_CAPTURED
    
  }
  
  
  #quit(status=1)
  #### POOL G Freq table ----
  
  
  Rosetta_POOL_G<-MPRA_Rosetta[which(MPRA_Rosetta$factor4 == "High_GC"),]
  
  # cat("Rosetta_POOL_G\n")
  # str(Rosetta_POOL_G)
  # cat("\n")
  
  TABLE_DEF_Total_POOL_G<-TABLE_DEF_Total[which(TABLE_DEF_Total$seq_name%in%Rosetta_POOL_G$seq_name),]
  
  if(dim(TABLE_DEF_Total_POOL_G)[1]>0)
  {
    TABLE_DEF_Total_POOL_G$pool<-"POOL_G"
    
    
    cat("TABLE_DEF_Total_POOL_G\n")
    str(TABLE_DEF_Total_POOL_G)
    cat("\n")
    
    TABLE_DEF_Total_POOL_G_CAPTURED<-as.data.frame(setDT(TABLE_DEF_Total_POOL_G)[, .N, by=.(lane,subsampling,pool)], stringsAsFactors=F)
    
    colnames(TABLE_DEF_Total_POOL_G_CAPTURED)[which(colnames(TABLE_DEF_Total_POOL_G_CAPTURED) == "N")]<-"CAPTURED"
    
    TABLE_DEF_Total_POOL_G_CAPTURED$Perc_captured<-100*(TABLE_DEF_Total_POOL_G_CAPTURED$CAPTURED/dim(Rosetta_POOL_G)[1])
    
    
    
    cat("TABLE_DEF_Total_POOL_G_CAPTURED\n")
    str(TABLE_DEF_Total_POOL_G_CAPTURED)
    cat("\n")
    
    
    
    list_DEF$tagging_barcodes_df[[2]]<-TABLE_DEF_Total_POOL_G
    list_DEF$summary_df[[2]]<-TABLE_DEF_Total_POOL_G_CAPTURED
    
  }
  
  
  #### POOL M Freq table ----
  
  
  Rosetta_POOL_M<-MPRA_Rosetta[which(MPRA_Rosetta$factor4 == "Medium"),]
  
  cat("Rosetta_POOL_M\n")
  str(Rosetta_POOL_M)
  cat("\n")
  
  TABLE_DEF_Total_POOL_M<-TABLE_DEF_Total[which(TABLE_DEF_Total$seq_name%in%Rosetta_POOL_M$seq_name),]
  
  if(dim(TABLE_DEF_Total_POOL_M)[1]>0)
  {
    
    TABLE_DEF_Total_POOL_M$pool<-"POOL_M"
    
    
    cat("TABLE_DEF_Total_POOL_M\n")
    str(TABLE_DEF_Total_POOL_M)
    cat("\n")
    
    TABLE_DEF_Total_POOL_M_CAPTURED<-as.data.frame(setDT(TABLE_DEF_Total_POOL_M)[, .N, by=.(lane,subsampling,pool)], stringsAsFactors=F)
    
    colnames(TABLE_DEF_Total_POOL_M_CAPTURED)[which(colnames(TABLE_DEF_Total_POOL_M_CAPTURED) == "N")]<-"CAPTURED"
    
    TABLE_DEF_Total_POOL_M_CAPTURED$Perc_captured<-100*(TABLE_DEF_Total_POOL_M_CAPTURED$CAPTURED/dim(Rosetta_POOL_M)[1])
    
    
    
    cat("TABLE_DEF_Total_POOL_M_CAPTURED\n")
    str(TABLE_DEF_Total_POOL_M_CAPTURED)
    cat("\n")
    
    
    
    list_DEF$tagging_barcodes_df[[3]]<-TABLE_DEF_Total_POOL_M
    list_DEF$summary_df[[3]]<-TABLE_DEF_Total_POOL_M_CAPTURED
  }
 
      
  
  #### FINAL PRINT AND SAVE ----
  
  TABLE_FINAL_SUMMARY = as.data.frame(data.table::rbindlist(list_DEF$summary_df, fill = T))
  
  cat("TABLE_FINAL_SUMMARY\n")
  str(TABLE_FINAL_SUMMARY)
  cat("\n")
  
  TABLE_FINAL_TAGGING_BARCODES = as.data.frame(data.table::rbindlist(list_DEF$tagging_barcodes_df, fill = T))
  
  cat("TABLE_FINAL_TAGGING_BARCODES\n")
  str(TABLE_FINAL_TAGGING_BARCODES)
  cat("\n")
  
  
  #### SAVE ----
  
  setwd(out)
  
  filename_18<-paste("SUMMARY_TABLE_",type,'.tsv',sep='')
  
  write.table(TABLE_FINAL_SUMMARY,file=filename_18,quote=F, sep="\t", row.names = F)
  
  setwd(out)
  
  filename_18<-paste("TAGGING_BARCODES_TABLE_",type,'.tsv',sep='')
  
  write.table(TABLE_FINAL_TAGGING_BARCODES,file=filename_18,quote=F, sep="\t", row.names = F)
  
  setwd(out)
  
  # filename_18<-paste("READS_BARCODES_TABLE_",type,'.tsv',sep='')
  # 
  # write.table(TABLE_DEF,file=filename_18,quote=F, sep="\t", row.names = F)
      
      # quit(status=1)
      
  
  
}

graphs = function(option_list)
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
  
  
  # quit(status=1)
  
  setwd(out)
  
  filename_18<-paste("SUMMARY_TABLE_",type,'.tsv',sep='')
  
  SUMMARY_TABLE<-as.data.frame(fread(file=filename_18, sep="\t", header = T), stringsAsFactors = F)
  
  cat("SUMMARY_TABLE_\n")
  cat(str(SUMMARY_TABLE))
  cat("\n")
  
  filename_18<-paste("TAGGING_BARCODES_TABLE_",type,'.tsv',sep='')
  
  TAGGING_BARCODES_TABLE<-as.data.frame(fread(file=filename_18, sep="\t", header = T), stringsAsFactors = F)
  
  cat("TAGGING_BARCODES_TABLE_\n")
  cat(str(TAGGING_BARCODES_TABLE))
  cat("\n")
  
  # filename_18<-paste("READS_BARCODES_TABLE_",type,'.tsv',sep='')
  # 
  # READS_BARCODES_TABLE<-as.data.frame(fread(file=filename_18, sep="\t", header = T), stringsAsFactors = F)
  # 
  # cat("READS_BARCODES_TABLE_\n")
  # cat(str(READS_BARCODES_TABLE))
  # cat("\n")
  
  
  
  
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
    make_option(c("--LD_string"), type="character", default=NULL, 
                metavar="type", 
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
  # graphs(opt)
  
  
}


###########################################################################

system.time( main() )
