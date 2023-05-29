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

library("svglite", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")
library("cowplot",lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")
library("digest",lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")
library("farver",lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")
library("labeling",lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")
library("ggeasy",lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")
library("seqinr", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")

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
  
  cat("file_list_sel_0\n")
  cat(str(file_list_sel))
  cat("\n")
  
 
  
 
  
  file_list_sel$lane<-gsub("^N*o*Y*e*s*NM_MAPQ_[0-9]+_[0-9]+_[0-9]+_",'',file_list_sel$file)
  file_list_sel$lane<-gsub(".rds$",'',file_list_sel$lane)
  file_list_sel$lane<-gsub("_0\\.[0-9]*$",'',file_list_sel$lane)
  
  
  indexes_sel_2 <- grep("MAPQ",file_list_sel$file)
  
  cat("indexes_sel_2\n")
  cat(sprintf(as.character(indexes_sel_2)))
  cat("\n")
  
  file_list_sel<-file_list_sel[indexes_sel_2,]
  
  cat("file_list_sel_1\n")
  cat(str(file_list_sel))
  cat("\n")
  
  
  
  
  lane_substitution<-paste(as.character(names(summary(as.factor(file_list_sel$lane)))), collapse="|")
  
  cat("lane_substitution\n")
  cat(sprintf(as.character(names(summary(as.factor(lane_substitution))))))
  cat("\n")
  
  file_list_sel$lane<-gsub("lane_1_lane_1",'lane_1',file_list_sel$lane)
  file_list_sel$lane<-gsub("lane_2_lane_2",'lane_2',file_list_sel$lane)
  
  
  
  cat("lane\n")
  cat(sprintf(as.character(names(summary(as.factor(file_list_sel$lane))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(file_list_sel$lane)))))
  cat("\n")
  
  
 
  
 
  
  file_list_sel$subsampling<-gsub(lane_substitution,'',file_list_sel$file)
  file_list_sel$subsampling<-gsub("^N*o*Y*e*s*NM_MAPQ_[0-9]+_[0-9]+_[0-9]+_",'',file_list_sel$subsampling)
  file_list_sel$subsampling<-gsub(".rds$",'',file_list_sel$subsampling)
  file_list_sel$subsampling[file_list_sel$subsampling == ""]<-"100_percent"
  file_list_sel$subsampling<-gsub("^_*",'',file_list_sel$subsampling)
  
  cat("subsampling\n")
  cat(sprintf(as.character(names(summary(as.factor(file_list_sel$subsampling))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(file_list_sel$subsampling)))))
  cat("\n")
  
  subsampling_substitution<-paste(as.character(names(summary(as.factor(file_list_sel$subsampling)))), collapse="|")
  
  cat("subsampling_substitution\n")
  cat(sprintf(as.character(names(summary(as.factor(subsampling_substitution))))))
  cat("\n")
  
  file_list_sel$bc_per_tile<-gsub(lane_substitution,'',file_list_sel$file)
  
  cat("file_list_sel$bc_per_tile_1\n")
  cat(sprintf(as.character(file_list_sel$bc_per_tile)))
  cat("\n")
  
  file_list_sel$bc_per_tile<-gsub("__.+$",'',file_list_sel$bc_per_tile)
  
  cat("file_list_sel$bc_per_tile_1.5\n")
  cat(sprintf(as.character(file_list_sel$bc_per_tile)))
  cat("\n")
  
  
  # file_list_sel$bc_per_tile<-gsub(subsampling_substitution,'',file_list_sel$bc_per_tile)
  # 
  # cat("file_list_sel$bc_per_tile_2\n")
  # cat(sprintf(as.character(file_list_sel$bc_per_tile)))
  # cat("\n")
  
  file_list_sel$bc_per_tile<-gsub("^N*o*Y*e*s*NM_MAPQ_[0-9]+_",'',file_list_sel$bc_per_tile)
  
  cat("file_list_sel$bc_per_tile_3\n")
  cat(sprintf(as.character(file_list_sel$bc_per_tile)))
  cat("\n")
  
  file_list_sel$bc_per_tile<-gsub("_.+",'',file_list_sel$bc_per_tile)
  
  cat("file_list_sel$bc_per_tile_4\n")
  cat(sprintf(as.character(file_list_sel$bc_per_tile)))
  cat("\n")
  
  file_list_sel$bc_per_tile<-as.numeric(file_list_sel$bc_per_tile)
  
  
  cat("bc_per_tile\n")
  cat(sprintf(as.character(names(summary(as.factor(file_list_sel$bc_per_tile))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(file_list_sel$bc_per_tile)))))
  cat("\n")
  
  file_list_sel$Threshold_reads_per_bc<-gsub(lane_substitution,'',file_list_sel$file)
  
  cat("file_list_sel$Threshold_reads_per_bc_1\n")
  cat(sprintf(as.character(file_list_sel$Threshold_reads_per_bc)))
  cat("\n")
  
  file_list_sel$Threshold_reads_per_bc<-gsub("__.+$",'',file_list_sel$Threshold_reads_per_bc)
  
  cat("file_list_sel$Threshold_reads_per_bc_1.5\n")
  cat(sprintf(as.character(file_list_sel$Threshold_reads_per_bc)))
  cat("\n")
  
  
  
  file_list_sel$Threshold_reads_per_bc<-gsub("^N*o*Y*e*s*NM_MAPQ_[0-9]+_[0-9]+_",'',file_list_sel$Threshold_reads_per_bc)
  file_list_sel$Threshold_reads_per_bc<-gsub("_.+",'',file_list_sel$Threshold_reads_per_bc)
  
  file_list_sel$Threshold_reads_per_bc<-as.numeric(file_list_sel$Threshold_reads_per_bc)
  
  
  cat("file_list_sel\n")
  cat(str(file_list_sel))
  cat("\n")
  
  #### HERE HERE HERE
  
  
  
  cat("Threshold_reads_per_bc\n")
  cat(sprintf(as.character(names(summary(as.factor(file_list_sel$Threshold_reads_per_bc))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(file_list_sel$Threshold_reads_per_bc)))))
  cat("\n")
  
  file_list_sel$NM<-gsub(lane_substitution,'',file_list_sel$file)
  file_list_sel$NM<-gsub(subsampling_substitution,'',file_list_sel$NM)
  
  # file_list_sel$NM<-gsub("^N*o*Y*e*s*NM__",'',file_list_sel$NM)
  file_list_sel$NM<-gsub("_.+",'',file_list_sel$NM)
  # file_list_sel$NM<-paste("N*o*Y*e*s*NM__",file_list_sel$NM, sep='')
  
  cat("NM\n")
  cat(sprintf(as.character(names(summary(as.factor(file_list_sel$NM))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(file_list_sel$NM)))))
  cat("\n")
  
  file_list_sel$MAPQ<-gsub(lane_substitution,'',file_list_sel$file)
  file_list_sel$MAPQ<-gsub(subsampling_substitution,'',file_list_sel$MAPQ)
  
  file_list_sel$MAPQ<-gsub("^N*o*Y*e*s*NM_MAPQ_",'',file_list_sel$MAPQ)
  file_list_sel$MAPQ<-gsub("_.+",'',file_list_sel$MAPQ)
  file_list_sel$MAPQ<-paste("MAPQ_",file_list_sel$MAPQ, sep='')
  
  
  cat("MAPQ\n")
  cat(sprintf(as.character(names(summary(as.factor(file_list_sel$MAPQ))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(file_list_sel$MAPQ)))))
  cat("\n")
  
  
  
  file_list_sel<-file_list_sel[order(file_list_sel$bc_per_tile, decreasing = F),]
  
  file_list_sel$parameters<-paste(file_list_sel$NM,file_list_sel$MAPQ,file_list_sel$bc_per_tile,
                                  file_list_sel$Threshold_reads_per_bc, sep="_")
  
  cat("file_list_sel\n")
  cat(str(file_list_sel))
  cat("\n")
  
  # setwd(out)
  # write.table(file_list_sel,file="test.tsv", sep="\t", quote = F, row.names = F)
  # quit(status=1)
  
  
  parameters<-as.character(names(summary(as.factor(file_list_sel$parameters))))
  subsampling<-rev(c("100_percent","0.9","0.75","0.5","0.25","0.1","0.01"))
  lane<-c("lane_1","lane_2")
  
  
  file_list_sel$parameters<-factor(file_list_sel$parameters,
                             levels=parameters,
                             ordered=T)
  
  file_list_sel$subsampling<-factor(file_list_sel$subsampling,
                                    levels=subsampling,
                                    ordered=T)
  
  file_list_sel$lane<-factor(file_list_sel$lane,
                                    levels=lane,
                                    ordered=T)
  
  
  cat("file_list_sel\n")
  cat(str(file_list_sel))
  cat("\n")
  
  
  write.table(file_list_sel,file="file_list_sel.tsv",sep="\t", quote=F, row.names = F)
  
  # #########################################################################################################################################################
  # quit(status=1)
  # 
  
  
  #### LOOP READ RDS files ----
  
  List_parameters<-list()
  
  List_DEF<-list()
  
  for(i in 1:length(levels(file_list_sel$parameters)))
  {
    
    parameters_sel<-levels(file_list_sel$parameters)[i]
    
    cat("parameters_sel\n")
    cat(sprintf(as.character(parameters_sel)))
    cat("\n")
    
    # #### path2 ----
    # 
    # path2<-paste(out,parameters_sel, sep='')
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
    
    
    
    file_list_subset<-file_list_sel[which(file_list_sel$parameters%in%parameters_sel),]

    # cat("file_list_subset\n")
    # cat(str(file_list_subset))
    # cat("\n")
    
    
    if(dim(file_list_subset)[1]>0)
    {
      
      # file_list_subset$lane<-droplevels(file_list_subset$lane)
      # 
      # cat("file_list_subset\n")
      # cat(str(file_list_subset))
      # cat("\n")
      
      List_lane<-list()
      
      for(l in 1:length(levels(file_list_subset$lane)))
      {
        lane_sel<-levels(file_list_subset$lane)[l]
        
        # cat("lane_sel\n")
        # cat(sprintf(as.character(lane_sel)))
        # cat("\n")
        
        file_list_subset_lane<-file_list_subset[which(file_list_subset$lane%in%lane_sel),]
        
        # cat("file_list_subset_lane\n")
        # cat(str(file_list_subset_lane))
        # cat("\n")
        
        
        
        if(dim(file_list_subset_lane)[1]>0)
        {
          List_RESULTS<-list()
          
          for(k in 1:dim(file_list_subset_lane)[1])
          {
            subsampling_sel<-file_list_subset_lane$subsampling[k]
            
            # cat("subsampling_sel\n")
            # cat(sprintf(as.character(subsampling_sel)))
            # cat("\n")
            
            file_sel<-file_list_subset_lane$file[file_list_subset_lane$subsampling == subsampling_sel]
            
            # cat("file_sel\n")
            # cat(sprintf(as.character(file_sel)))
            # cat("\n")
            
            if(length(file_sel) > 0)
            {
              setwd(out)
              
              df<-readRDS(file = file_sel)
              
             if(dim(df)[1] >0)
              {
                df$parameters<-parameters_sel
                df$lane<-lane_sel
                df$subsampling<-subsampling_sel
                
                # cat("df\n")
                # cat(str(df))
                # cat("\n")
                
                List_RESULTS[[k]]<-df
              }
              
              
              
              # quit(status=1)
              
            }# length(file_sel) > 0
            
            
          } #k subsampling_sel
          
          
          parameters_TABLE = as.data.frame(data.table::rbindlist(List_RESULTS, fill = T))
          
          # cat("parameters_TABLE\n")
          # str(parameters_TABLE)
          # cat("\n")
          
          List_lane[[l]]<-parameters_TABLE
          
          
          # quit(status = 1)
           
        }# dim(file_list_subset_lane)[1]>0 
        
      }#l lanes
      
      parameters_TABLE_lane = as.data.frame(data.table::rbindlist(List_lane, fill = T))
      
      # cat("parameters_TABLE_lane\n")
      # str(parameters_TABLE_lane)
      # cat("\n")
     
      # quit(status = 1)
      List_DEF[[i]]<-parameters_TABLE_lane
      
    }# dim(file_list_subset)[1]>0
    
    
  }  #i
  
  
  ####
  
  TABLE_DEF = as.data.frame(data.table::rbindlist(List_DEF, fill = T))
  
  cat("TABLE_DEF\n")
  str(TABLE_DEF)
  cat("\n")
  
  TABLE_DEF$parameters<-factor(TABLE_DEF$parameters,
                                   levels=parameters,
                                   ordered=T)
  
  TABLE_DEF$subsampling<-factor(TABLE_DEF$subsampling,
                                    levels=subsampling,
                                    ordered=T)
  
  TABLE_DEF$lane<-factor(TABLE_DEF$lane,
                             levels=lane,
                             ordered=T)
  
  cat("TABLE_DEF\n")
  str(TABLE_DEF)
  cat("\n")
  
  
  TABLE_DEF_Total<-as.data.frame(setDT(TABLE_DEF)[, .(.N,paste(Barcode, collapse = "|"),
                                                      paste(number_of_reads, collapse = "|")), 
                                                  by=.(parameters,subsampling,lane,seq_name)], stringsAsFactors=F)
  
  cat("TABLE_DEF_Total\n")
  str(TABLE_DEF_Total)
  cat("\n")
  
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
    
    TABLE_DEF_Total_POOL_A_CAPTURED<-as.data.frame(setDT(TABLE_DEF_Total_POOL_A)[, .N, by=.(lane,parameters,subsampling,pool)], stringsAsFactors=F)
    
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
    
    TABLE_DEF_Total_POOL_G_CAPTURED<-as.data.frame(setDT(TABLE_DEF_Total_POOL_G)[, .N, by=.(lane,parameters,subsampling,pool)], stringsAsFactors=F)
    
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
    
    TABLE_DEF_Total_POOL_M_CAPTURED<-as.data.frame(setDT(TABLE_DEF_Total_POOL_M)[, .N, by=.(lane,parameters,subsampling,pool)], stringsAsFactors=F)
    
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
  
  filename_18<-paste("SUMMARY_TABLE_",type,'.rds',sep='')
  
  saveRDS(TABLE_FINAL_SUMMARY,file=filename_18)
  
  setwd(out)
  
  filename_18<-paste("TAGGING_BARCODES_TABLE_",type,'.rds',sep='')
  
  saveRDS(TABLE_FINAL_TAGGING_BARCODES,file=filename_18)
  
  # filename_18<-paste("READS_BARCODES_TABLE_",type,'.tsv',sep='')
  # 
  # write.table(TABLE_DEF,file=filename_18,quote=F, sep="\t", row.names = F)
  
  # quit(status=1)
  
  
  
}


Reference_fasta_printer = function(option_list)
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
  
  #### READ and transform QC_dnaALT_REF_threshold ----
  
  QC_dnaALT_REF_threshold = as.numeric(unlist(strsplit(opt$QC_dnaALT_REF_threshold, split=",")))
  
  
  
  cat("QC_dnaALT_REF_threshold\n")
  cat(sprintf(as.character(QC_dnaALT_REF_threshold)))
  cat("\n")
  
  #### READ MPRA_Rosetta ----
  
  
  MPRA_Rosetta<-as.data.frame(fread(opt$MPRA_Rosetta, sep="\t", header = T), stringsAsFactors = F)
  
  cat("MPRA_Rosetta_\n")
  cat(str(MPRA_Rosetta))
  cat("\n")
  
  indx.int<-c(which(colnames(MPRA_Rosetta) == "REAL_TILE"),which(colnames(MPRA_Rosetta) == "seq_name"),which(colnames(MPRA_Rosetta) == "carried_variants"),
              which(colnames(MPRA_Rosetta) == "chr"),which(colnames(MPRA_Rosetta) == "start"))
  
  
  
  MPRA_Rosetta_subset<-unique(MPRA_Rosetta[,indx.int])
  
  
 
  cat("MPRA_Rosetta_subset_\n")
  cat(str(MPRA_Rosetta_subset))
  cat("\n")
  
  indx.int<-c(which(colnames(MPRA_Rosetta) == "REAL_TILE"),which(colnames(MPRA_Rosetta) == "carried_variants"),
             which(colnames(MPRA_Rosetta) == "factor4"))
  
  
  
  MPRA_Rosetta_subset_2<-unique(MPRA_Rosetta[,indx.int])
  
  MPRA_Rosetta_subset_2$pool<-"NA"
  
  cat("MPRA_Rosetta_subset_2_\n")
  cat(str(MPRA_Rosetta_subset_2))
  cat("\n")
  
  
  MPRA_Rosetta_subset_2$pool[which(MPRA_Rosetta_subset_2$factor4 == "High_AT")]<-"POOL_A"
  MPRA_Rosetta_subset_2$pool[which(MPRA_Rosetta_subset_2$factor4 == "High_GC")]<-"POOL_G"
  MPRA_Rosetta_subset_2$pool[which(MPRA_Rosetta_subset_2$factor4 == "Medium")]<-"POOL_M"
  
  cat("MPRA_Rosetta_subset_2_\n")
  cat(str(MPRA_Rosetta_subset_2))
  cat("\n")
  
  
  MPRA_Rosetta_subset_2_REF<-MPRA_Rosetta_subset_2[which(MPRA_Rosetta_subset_2$carried_variants == "REF"),]
  
  cat("MPRA_Rosetta_subset_2_REF\n")
  cat(str(MPRA_Rosetta_subset_2_REF))
  cat("\n")
  
  MPRA_Rosetta_subset_2_ALT<-MPRA_Rosetta_subset_2[which(MPRA_Rosetta_subset_2$carried_variants != "REF"),]
  
  cat("MPRA_Rosetta_subset_2_ALT\n")
  cat(str(MPRA_Rosetta_subset_2_ALT))
  cat("\n")
  
  MPRA_Rosetta_RT_format<-merge(MPRA_Rosetta_subset_2_REF,MPRA_Rosetta_subset_2_ALT,
                                by=c("REAL_TILE","factor4","pool"))
  
  colnames(MPRA_Rosetta_RT_format)<-gsub(".x$",".REF",colnames(MPRA_Rosetta_RT_format))
  colnames(MPRA_Rosetta_RT_format)<-gsub(".y$",".ALT",colnames(MPRA_Rosetta_RT_format))
  
  cat("MPRA_Rosetta_RT_format\n")
  cat(str(MPRA_Rosetta_RT_format))
  cat("\n")
  
  #### READ NovaSeq_1 ----
  
  NovaSeq_1<-readRDS(file=opt$NovaSeq_1)
  
  NovaSeq_1$batch<-"NovaSeq_1"
  
  cat("NovaSeq_1_\n")
  cat(str(NovaSeq_1))
  cat("\n")
  
  #### READ NovaSeq_1_TILES ----
  
  NovaSeq_1_TILES<-readRDS(file=opt$NovaSeq_1_TILES)
  
  NovaSeq_1_TILES$batch<-"NovaSeq_1"
  
  cat("NovaSeq_1_TILES_\n")
  cat(str(NovaSeq_1_TILES))
  cat("\n")
  
 # quit(status=1)
  
  #### SUMMARY_PLOT_LEVEL ----
  
  REP<-rbind(NovaSeq_1)
  
  REP$pool<-factor(REP$pool,
                   levels=c("POOL_A","POOL_M","POOL_G"),
                   ordered=T)
  
  breaks.x<-levels(REP$subsampling)
  
  labels.x<-c('1','10','25','50','75','90','100')
  
  cat("labels.x\n")
  cat(sprintf(as.character(labels.x)))
  cat("\n")
  cat(str(breaks.x))
  cat("\n")
  cat(str(labels.x))
  cat("\n")
  
  breaks.y<-seq(0,100,by=10)
  
  labels.y<-as.character(breaks.y)
  
  cat("labels.y\n")
  cat(sprintf(as.character(labels.y)))
  cat("\n")
  cat(str(breaks.y))
  cat("\n")
  cat(str(labels.y))
  cat("\n")
  
  
list_graph_DEF<-list()
  
for(i in 1:length(levels(REP$parameters)))
{
    parameters_sel<-levels(REP$parameters)[i]
    
    
    cat(sprintf(as.character(parameters_sel)))
    cat("\n")
    
    
    REP_sel<-REP[which(REP$parameters%in%parameters_sel),]
    
    cat("REP_sel_\n")
    cat(str(REP_sel))
    cat("\n")
    
    list_graphs_per_lane<-list()
    
    for(k in 1:length(levels(REP_sel$lane)))
    {
      lane_sel<-levels(REP_sel$lane)[k]
      
      
      cat(sprintf(as.character(lane_sel)))
      cat("\n")
      
      REP_sel_lane<-REP_sel[which(REP_sel$lane == lane_sel),]
      
      cat("REP_sel_lane\n")
      cat(str(REP_sel_lane))
      cat("\n")
      
      graph<-ggplot(data=REP_sel_lane,aes(x=subsampling, y=Perc_captured, fill=pool)) +
        geom_bar(stat="identity",colour='black')+
        theme(plot.title=element_text(size=11))+
        theme(axis.text.x=element_text(angle=45,size=14,colour="black",vjust=1,hjust=1,face="bold"), 
              axis.title.y=element_text(size=16,face="bold"),
              legend.title=element_text(size=16,face="bold"),
              legend.text=element_text(size=12,face="bold",colour="black"),
              axis.text.y=element_text(colour="black",size=12,face="bold"))+
        theme(plot.title = element_text(size=11)) +
        theme(plot.title=element_text(size=11))+
        scale_x_discrete(name="% of NovaSeq analysed", breaks=breaks.x,labels=labels.x, drop=F)+
        scale_y_continuous(name="% of POOL captured",breaks=breaks.y,labels=labels.y, limits=c(breaks.y[1],breaks.y[length(breaks.y)]))+
        labs(x="",y="Percentage captured", fill=paste("Pool", sep="\n"))+
        scale_fill_manual(values=c('#32A852','#6DB2EE','#553B68','#62D07F','#1877C9','#C9244B','#D45E85','#87447B','#D6B8E6'),drop=F)+
        theme_bw() + theme(legend.position = "none")+
        ggeasy::easy_center_title()
      
      cat("HELLO_WORLD_1\n")
      
      graph <- graph + facet_grid(. ~ pool, drop=F)
      
      cat("HELLO_WORLD_2\n")
      
      
      list_graphs_per_lane[[lane_sel]]<-graph
          
    }#k lanes
    
    names_array<-names(list_graphs_per_lane)
    
    cat("names_array_\n")
    cat(str(names_array))
    cat("\n")
    
    # quit(status=1)
    
    
    lane_1_graph<-list_graphs_per_lane[['lane_1']]
    lane_2_graph<-list_graphs_per_lane[['lane_2']]
    
    
    title <- ggdraw() + draw_label(paste("NovaSeq results",parameters_sel,sep="  "))
    
    graph_overview_DEF<-plot_grid(title,lane_1_graph,lane_2_graph,
                                  nrow = 3,
                                  labels = c("","A) lane_1","B) lane_2"),
                                  label_size = 8,
                                  vjust=0.5,
                                  align = "v",
                                  rel_heights = c(0.05, 1,1))
    
    list_graph_DEF[[i]]<-graph_overview_DEF
    
    # pdf(file="test.pdf")
    # print(graph_overview_DEF)
    # dev.off()
    # 
    # quit(status=1)
    
    
}# i parameters
  
  
  # quit(status=1)
  
  

  #### BC and reads level ----
  
  BC_table<-rbind(NovaSeq_1_TILES)
  
  cat("BC_table_\n")
  cat(str(BC_table))
  cat("\n")
  
  
  list_BC_deconvolution<-list()
  
  BC_table$subsampling<-factor(BC_table$subsampling,
                          levels=c("0.01","0.1","0.25","0.5","0.75","0.9","100_percent"),
                          ordered=T)
  
  
  BC_table_subset<-BC_table[which(BC_table$subsampling == "100_percent"),]
  
  
  levels_MAPQ<-unique(as.character(BC_table_subset$parameters))
  
  BC_table_subset$parameters<-factor(BC_table_subset$parameters,
                   levels=levels_MAPQ,
                   ordered=T)
  
  BC_table_subset$pool<-factor(BC_table_subset$pool,
                   levels=c("POOL_A","POOL_M","POOL_G"),
                   ordered=T)
  
  
  BC_table_subset$batch<-factor(BC_table_subset$batch,
                    levels=c("NovaSeq_1"),
                    ordered=T)
  
  cat("BC_table_subset_\n")
  cat(str(BC_table_subset))
  cat("\n")
  
  # ##HERE HERE
  # quit(status = 1)
  
  list_Q_DEF<-list()
  list_FASTA_DEF<-list()
  
  for(i in 1:length(levels(BC_table_subset$batch)))
  {
    
    batch_sel<-levels(BC_table_subset$batch)[i]
    
    
    cat("----->\t")
    cat(sprintf(as.character(batch_sel)))
    cat("\t")
    
    
    list_Q_parameters<-list()
    list_FASTA_parameters<-list()
    
    
    for(k in 1:length(levels(BC_table_subset$parameters)))
    {
      parameters_sel<-levels(BC_table_subset$parameters)[k]
      
      
      cat(sprintf(as.character(parameters_sel)))
      cat("\n")
      
      
      BC_table_subset_sel<-BC_table_subset[which(BC_table_subset$batch%in%batch_sel &
                           BC_table_subset$parameters%in%parameters_sel),]
      
      cat("BC_table_subset_sel_\n")
      cat(str(BC_table_subset_sel))
      cat("\n")
      
      
      list_Q_bc<-list()
      list_FASTA_bc<-list()
      
        for(h in 1:dim(BC_table_subset_sel)[1])
        {
          
          bc_sel<-BC_table_subset_sel[h,-c(which(colnames(BC_table_subset_sel) == "parameters"),
                                           which(colnames(BC_table_subset_sel) == "batch"),
                                                 which(colnames(BC_table_subset_sel) == "subsampling"))]
          
          bc_sel$count = sum(as.numeric(unlist(strsplit(bc_sel$number_of_reads_string, split="\\|"))))
          
          
          
          MPRA_Rosetta_subset_sel<-MPRA_Rosetta_subset[which(MPRA_Rosetta_subset$seq_name%in%bc_sel$seq_name),]
          
          bc_sel<-merge(bc_sel,
                        MPRA_Rosetta_subset_sel,
                        by=c("seq_name"),
                        all.x=T)
          # cat("bc_sel_\n")
          # cat(str(bc_sel))
          # cat("\n")
          
          barcode_string<-unlist(strsplit(bc_sel$Barcode_string, split="\\|"))
          length_barcode_string<-length(barcode_string)
          
          # cat("barcode_string_\n")
          # cat(str(barcode_string))
          # cat("\n")
          
          a.df<-as.data.frame(cbind(rep(MPRA_Rosetta_subset_sel$chr, length_barcode_string),
                                    rep(MPRA_Rosetta_subset_sel$start, length_barcode_string),
                              rep(MPRA_Rosetta_subset_sel$REAL_TILE, length_barcode_string),
                              rep(MPRA_Rosetta_subset_sel$carried_variants, length_barcode_string),
                              rep(MPRA_Rosetta_subset_sel$seq_name, length_barcode_string),
                              barcode_string), stringsAsFactors=F)
          
          colnames(a.df)<-c("chr","start","REAL_TILE","carried_variants","seq_name","Barcode")
          
          a.df$ALLELE<-"NA"
          
          if(a.df$carried_variant=="REF")
          {
            
            a.df$ALLELE<-"REF"
          }else{
            
            a.df$ALLELE<-"ALT"
          }
          
          
          # cat("a.df_\n")
          # cat(str(a.df))
          # cat("\n")
          
          list_FASTA_bc[[h]]<-a.df
          list_Q_bc[[h]]<-bc_sel
        
          
          # cat("MPRA_Rosetta_subset_sel_\n")
          # cat(str(MPRA_Rosetta_subset_sel))
          # cat("\n")
          
          # quit(status=1)
          
        }#h
      
      TABLE_FASTA = unique(as.data.frame(data.table::rbindlist(list_FASTA_bc, fill = T)))
      
      # cat("TABLE_FASTA\n")
      # str(TABLE_FASTA)
      # cat("\n")
      
      list_FASTA_parameters[[k]]<-TABLE_FASTA
      
      TABLE_Q = unique(as.data.frame(data.table::rbindlist(list_Q_bc, fill = T)))
      
      # cat("TABLE_Q\n")
      # str(TABLE_Q)
      # cat("\n")
      
      list_Q_parameters[[k]]<-TABLE_Q
      
    #  quit(status=1)
      
     
      
    }#k
    
    TABLE_FASTA = unique(as.data.frame(data.table::rbindlist(list_FASTA_parameters, fill = T)))
    
    # cat("TABLE_FASTA\n")
    # str(TABLE_FASTA)
    # cat("\n")
    
    list_FASTA_DEF[[i]]<-TABLE_FASTA
    
    TABLE_Q = unique(as.data.frame(data.table::rbindlist(list_Q_parameters, fill = T)))
    
    # cat("TABLE_Q\n")
    # str(TABLE_Q)
    # cat("\n")
    
    list_Q_DEF[[i]]<-TABLE_Q
    
  }#i
   
  
  TABLE_FASTA = unique(as.data.frame(data.table::rbindlist(list_FASTA_DEF, fill = T)))
  
  
  TABLE_FASTA$ALLELE<-factor(TABLE_FASTA$ALLELE,
                             levels=c("REF","ALT"),
                             ordered=TRUE)
  
  cat("TABLE_FASTA\n")
  str(TABLE_FASTA)
  cat("\n")
  
  #### FASTA PRINTER----
  
  TABLE_FASTA_ordered<-TABLE_FASTA[order(TABLE_FASTA$chr, TABLE_FASTA$start, TABLE_FASTA$REAL_TILE, TABLE_FASTA$ALLELE, TABLE_FASTA$carried_variants),]
  
  TABLE_FASTA_ordered$seq_name_bc<-paste(TABLE_FASTA_ordered$seq_name,TABLE_FASTA_ordered$Barcode, sep=";")
  
  cat("TABLE_FASTA_ordered\n")
  str(TABLE_FASTA_ordered)
  cat("\n")
  
  indx.fasta<-c(which(colnames(TABLE_FASTA_ordered) == "seq_name_bc"),which(colnames(TABLE_FASTA_ordered) == "Barcode"))
  
  
  
  Merge2_TO_FASTA<-unique(TABLE_FASTA_ordered[,indx.fasta])
  
  cat("Merge2_TO_FASTA\n")
  str(Merge2_TO_FASTA)
  cat("\n")
  
 
  
  sequences<-list()
  
  for(i in 1:length(Merge2_TO_FASTA$Barcode))
  {
    sequences[[i]]<-Merge2_TO_FASTA$Barcode[i]
    
  }#i
  
  
  cat("sequences\n")
  str(sequences)
  cat("\n")
  
  names<-list()
  
  for(i in 1:length(Merge2_TO_FASTA$seq_name_bc))
  {
    names[[i]]<-Merge2_TO_FASTA$seq_name_bc[i]
    
  }#i
  
  
  
  cat("names\n")
  str(names)
  cat("\n")
  
  setwd(out)
  
  filename<-paste("reference",'.fasta',sep='')
  
  
  write.fasta(sequences, names, filename, open = "w", nbchar = 60, as.string = FALSE)
  
  
 
  #### QC processing ----
  TABLE_Q = unique(as.data.frame(data.table::rbindlist(list_Q_DEF, fill = T)))
  
  cat("TABLE_Q\n")
  str(TABLE_Q)
  cat("\n")
  
  TABLE_Q_REF<-TABLE_Q[which(TABLE_Q$carried_variants == "REF"),]
  
  cat("TABLE_Q_REF\n")
  str(TABLE_Q_REF)
  cat("\n")
  
  TABLE_Q_ALT<-TABLE_Q[which(TABLE_Q$carried_variants != "REF"),]
  
  cat("TABLE_Q_ALT\n")
  str(TABLE_Q_ALT)
  cat("\n")
  
  TABLE_Q_DEF<-merge(TABLE_Q_REF,
                     TABLE_Q_ALT,
                     by=c("pool",
                          "REAL_TILE",
                          "chr",
                          "start"),
                     all=T)
  
  cat("TABLE_Q_DEF_0\n")
  str(TABLE_Q_DEF)
  cat("\n")
  
  colnames(TABLE_Q_DEF)<-gsub(".x$",".REF",colnames(TABLE_Q_DEF))
  colnames(TABLE_Q_DEF)<-gsub(".y$",".ALT",colnames(TABLE_Q_DEF))
  
  TABLE_Q_DEF$count.REF[is.na(TABLE_Q_DEF$count.REF)]<-0.01
  
  
  TABLE_Q_DEF$count.ALT[is.na(TABLE_Q_DEF$count.ALT)]<-0
  
  
  TABLE_Q_DEF$dna_ALT.REF<-TABLE_Q_DEF$count.ALT/TABLE_Q_DEF$count.REF
  TABLE_Q_DEF$dna_pool<-TABLE_Q_DEF$count.ALT+TABLE_Q_DEF$count.REF
  
  
  TABLE_Q_DEF$Dropout_dna_ALT.REF<-"PASS"
  TABLE_Q_DEF$Dropout_dna_ALT.REF[TABLE_Q_DEF$dna_ALT.REF <= QC_dnaALT_REF_threshold[1]]<-"Dropout_ASE_TOO_REF"
  TABLE_Q_DEF$Dropout_dna_ALT.REF[TABLE_Q_DEF$dna_ALT.REF >= QC_dnaALT_REF_threshold[2]]<-"Dropout_ASE_TOO_ALT"
  
  cat("TABLE_Q_DEF\n")
  str(TABLE_Q_DEF)
  cat("\n")
  
  cat("Dropout_dna_ALT.REF\n")
  cat(sprintf(as.character(levels(as.factor(TABLE_Q_DEF$Dropout_dna_ALT.REF)))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(TABLE_Q_DEF$Dropout_dna_ALT.REF)))))
  cat("\n")
  
  
  points_mean_dna_ALT.REF<-round(as.numeric(summary(TABLE_Q_DEF$dna_ALT.REF)),3)[-4]
  
  dna_ALT.REF_min<-points_mean_dna_ALT.REF[1]
  dna_ALT.REF_max<-points_mean_dna_ALT.REF[length(points_mean_dna_ALT.REF)]
  
  points_mean_dna_ALT.REF<-c(dna_ALT.REF_min,round(as.numeric(summary(TABLE_Q_DEF$dna_ALT.REF)[-4]),3),dna_ALT.REF_max)
  
  breaks.dna_ALT.REF<-sort(unique(c(1,points_mean_dna_ALT.REF)))
  breaks.dna_ALT.REF.log<-log(breaks.dna_ALT.REF+0.0001)
  labels.dna_ALT.REF<-as.character(round(breaks.dna_ALT.REF,1))
  
  
  
  cat("labels.dna_ALT.REF\n")
  cat(sprintf(as.character(labels.dna_ALT.REF)))
  cat("\n")
  cat(str(breaks.dna_ALT.REF))
  cat("\n")
  cat(str(labels.dna_ALT.REF))
  cat("\n")
  
  points_mean_dna_pool<-round(as.numeric(summary(TABLE_Q_DEF$dna_pool)),3)[-4]
  dna_pool_min<-points_mean_dna_pool[1]
  dna_pool_max<-points_mean_dna_pool[length(points_mean_dna_pool)]
  
  points_mean_dna_pool<-c(dna_pool_min,round(as.numeric(summary(TABLE_Q_DEF$dna_pool)[-4]),3),dna_pool_max)
  
  breaks.dna_pool<-sort(unique(c(1,points_mean_dna_pool)))
  breaks.dna_pool.log<-log(breaks.dna_pool+0.0001)
  labels.dna_pool<-as.character(round(breaks.dna_pool,0))
  
  
  cat("labels.dna_pool\n")
  cat(sprintf(as.character(labels.dna_pool)))
  cat("\n")
  cat(str(breaks.dna_pool))
  cat("\n")
  cat(str(labels.dna_pool))
  cat("\n")
  
  dynamic_range_minimum<-log(breaks.dna_pool[length(breaks.dna_pool)]/170+0.0001)
  
  g4 <- ggplot(TABLE_Q_DEF, aes(x=log(TABLE_Q_DEF$dna_ALT.REF +0.0001),
                        y=log(TABLE_Q_DEF$dna_pool +0.0001),
                        color=TABLE_Q_DEF$pool)) + 
    geom_point(size=2)+
    scale_color_manual(values =c('#32A852','#6DB2EE','#553B68','#62D07F','#1877C9','#C9244B','#D45E85','#87447B','#D6B8E6'),drop=F)+
    theme(legend.position = "bottom", legend.title = element_blank())+
    scale_y_continuous(name ="DNA counts",
                       breaks =breaks.dna_pool.log,
                       labels=labels.dna_pool,
                       limits = c(breaks.dna_pool.log[1],breaks.dna_pool.log[length(breaks.dna_pool.log)]))+
    scale_x_continuous(name ="ratio DNA ALT/REF",
                       breaks =breaks.dna_ALT.REF.log,
                       labels= labels.dna_ALT.REF,
                       limits= c(breaks.dna_ALT.REF.log[1],breaks.dna_ALT.REF.log[length(breaks.dna_ALT.REF.log)]))+
    theme(axis.text.x=element_text(angle=45,size=14,colour="black",vjust=1,hjust=1,face="bold"), 
          axis.title.y=element_text(size=16,face="bold"),
          axis.title.x=element_text(size=16,face="bold"),
          legend.title=element_blank(),
          legend.text=element_text(size=12,face="bold",colour="black"),
          axis.text.y=element_text(angle=45,colour="black",size=12,face="bold"))+
    geom_vline(
      xintercept = log(QC_dnaALT_REF_threshold[1]+0.0001),
      col = "red",
      linetype = "dotted",
      size = 1
    )+
    geom_vline(
      xintercept = log(QC_dnaALT_REF_threshold[2]+0.0001),
      col = "red",
      linetype = "dotted",
      size = 1
    )+
    geom_hline(
      yintercept = dynamic_range_minimum,
      col = "green",
      linetype = "dotted",
      size = 1
    )+
    theme_bw()+
    ggtitle(paste("REAL TILE QC thresholds",sep=" "))
  
  
  #### Frequency table for REAL TILES per pool that are: not cloned, dropout too REF, drop out too ALT and well cloned----
  
  cat("TABLE_Q_DEF\n")
  str(TABLE_Q_DEF)
  cat("\n")
  
  FINAL_table<-merge(MPRA_Rosetta_RT_format,
                     TABLE_Q_DEF,
                     by=c("REAL_TILE","pool","carried_variants.REF","carried_variants.ALT"),
                     all.x=T)
  
  FINAL_table$status="NA"
  FINAL_table$status[is.na(FINAL_table$count.REF)]<-"not_detected"
  

  FINAL_table$status[which(FINAL_table$Dropout_dna_ALT.REF == "Dropout_ASE_TOO_REF")]<-"Dropout_ASE_TOO_REF"
  FINAL_table$status[which(FINAL_table$Dropout_dna_ALT.REF == "Dropout_ASE_TOO_ALT")]<-"Dropout_ASE_TOO_ALT"
  FINAL_table$status[which(FINAL_table$Dropout_dna_ALT.REF == "PASS")]<-"PASS"
  
  

  cat("FINAL_table\n")
  str(FINAL_table)
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(FINAL_table$status))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(FINAL_table$status)))))
  cat("\n")
  
  Freq_table<-as.data.frame(table(FINAL_table$pool,FINAL_table$status), stringsAsFactors=F)
  
  colnames(Freq_table)<-c("pool","status","Freq")
  
  cat("Freq_table\n")
  str(Freq_table)
  cat("\n")
  
  Freq_table_total<-as.data.frame(table(FINAL_table$pool), stringsAsFactors=F)
  
  colnames(Freq_table_total)<-c("pool","Total")
  
  cat("Freq_table_total\n")
  str(Freq_table_total)
  cat("\n")
  
  Freq_table<-merge(Freq_table,
                    Freq_table_total,
                    by="pool")
  
  Freq_table$Perc<-round(100*(Freq_table$Freq/Freq_table$Total),2)
  
  Freq_table$pool<-factor(Freq_table$pool,
                               levels=c("POOL_A","POOL_M","POOL_G"),
                               ordered=T)
  Freq_table$status<-factor(Freq_table$status,
                                 levels=rev(c("PASS","Dropout_ASE_TOO_REF","Dropout_ASE_TOO_ALT","not_detected")),
                                 ordered=T)
  
  
  cat(sprintf(as.character(names(summary(as.factor(FINAL_table$pool))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(FINAL_table$pool)))))
  cat("\n")
 
  #HERE HERE
  quit(status=1)
  
  
  breaks.Rank<-seq(0,100, by=10)
  labels.Rank<-as.character(breaks.Rank)
  
  
  
  
  graph_SC<-Freq_table %>%
    mutate(myaxis = paste0(pool, "\n", "n=", Total)) %>%
    mutate(myaxis=fct_reorder(myaxis,as.numeric(pool))) %>%
    ggplot(aes(x=myaxis, y=Perc, fill=status)) +
    geom_bar(stat="identity",colour='black')+
    theme(plot.title=element_text(size=11))+
    scale_x_discrete(name=NULL, drop=F)+
    scale_y_continuous(name=NULL,breaks=breaks.Rank,labels=labels.Rank, limits=c(breaks.Rank[1],breaks.Rank[length(breaks.Rank)]))+
    theme(axis.text.x=element_text(angle=45,size=14,colour="black",vjust=1,hjust=1,face="bold"), 
          axis.title.y=element_text(size=16,face="bold"),
          legend.title=element_text(size=16,face="bold"),
          legend.text=element_text(size=12,face="bold",colour="black"),
          axis.text.y=element_text(colour="black",size=12,face="bold"))+
    labs(x="",y="Percentage", fill=paste("Cloning","status", sep="\n"))+
    scale_fill_manual(values=c("red",'#6DB2EE','#553B68','#32A852','#62D07F','#1877C9','#C9244B','#D45E85','#87447B','#D6B8E6'),drop=F)+
    theme_bw()+
    ggeasy::easy_center_title()
  
  
 # quit(status=1)
  
 #### PRINT graphs ----

setwd(out)
pdf("SUMMARY_plots.pdf")
for(i in 1:length(list_graph_DEF))
{
  graph_sel<-list_graph_DEF[[i]]
  
  print(graph_sel)
  
  svgname<-paste("Graphs_",i,"_","NovaSeq_subsampling",".svg", sep='')
  makesvg = TRUE
  
  if (makesvg == TRUE)
  {
    ggsave(svgname, plot= graph_sel,
           device="svg",
           height=10, width=12)
  }
  
}

print(g4)

svgname<-paste("Graphs_",type,"_QC_thresholds",".svg", sep='')
makesvg = TRUE

if (makesvg == TRUE)
{
  ggsave(svgname, plot= g4,
         device="svg",
         height=10, width=12)
}

print(graph_SC)

svgname<-paste("Graphs_",type,"_FINAL_stacked_barplot",".svg", sep='')
makesvg = TRUE

if (makesvg == TRUE)
{
  ggsave(svgname, plot= graph_SC,
         device="svg",
         height=10, width=12)
}


dev.off()
  
#### SAVE ----


setwd(out)

filename=paste(type,"_summary_REAL_TILE.tsv",sep='')

write.table(file=filename,FINAL_table, sep="\t", quote=F, row.names=F)

  
 # quit(status = 1)
  
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
    make_option(c("--MPRA_Rosetta"), type="character", default=NULL,
                metavar="FILE.txt",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--NovaSeq_1"), type="character", default=NULL,
                metavar="FILE.txt",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--NovaSeq_1_TILES"), type="character", default=NULL,
                metavar="FILE.txt",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--QC_dnaALT_REF_threshold"), type="character", default=NULL, 
                metavar="type1", 
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
  Reference_fasta_printer(opt)
  
 
  
  
}


###########################################################################

system.time( main() )
