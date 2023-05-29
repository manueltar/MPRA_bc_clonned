
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
library("svglite", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")
library("cowplot",lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")
library("digest",lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")
library("farver",lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")
library("labeling",lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")
library("reshape2",lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")

library("ggeasy",lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")

library("viridisLite",lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")

library("RColorBrewer",lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")

library("viridisLite",lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")
suppressMessages(library("cowplot", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("ggforce", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("ggnewscale", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("scales", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))

opt = NULL

options(warn=1)

cumulative_CT_Label2_and_Log_Rank_test_enhancer = function(option_list)
{
  suppressMessages(library("nph", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
  
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
  
  #### READ and transform out2 ----
  
  out2 = opt$out2
  
  cat("out2_\n")
  cat(sprintf(as.character(out2)))
  cat("\n")
  
  setwd(out)
  
  
  filename<-paste('df_Cell_colors','.rds',sep='')
  
  df_Cell_colors<-readRDS(file=filename)
  
  cat("df_Cell_colors_0\n")
  cat(str(df_Cell_colors))
  cat("\n")
  
  Condition_DEBUG <- 0
  
  ### CUMMULATIVE_CLASSES #----
  
  CUMMULATIVE_CLASSES<-readRDS(file=opt$CUMMULATIVE_CLASSES)
  
  CUMMULATIVE_CLASSES<-unique(CUMMULATIVE_CLASSES[,-c(which(colnames(CUMMULATIVE_CLASSES) == "VEP_DEF_LABELS"),
                                                      which(colnames(CUMMULATIVE_CLASSES) == "factor4_CLASS"))])
  
  
  cat("CUMMULATIVE_CLASSES_0\n")
  cat(str(CUMMULATIVE_CLASSES))
  cat("\n")
  cat(str(unique(CUMMULATIVE_CLASSES$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(CUMMULATIVE_CLASSES$bin_PP)))))
  cat("\n")
  cat(sprintf(as.character(summary(CUMMULATIVE_CLASSES$bin_PP))))
  cat("\n")
  
  
  
  # CUMMULATIVE_CLASSES<-droplevels(CUMMULATIVE_CLASSES[which(CUMMULATIVE_CLASSES$Cell_Type != "ALL_CT"),]) # Keep "ALL_CT"
  
  CUMMULATIVE_CLASSES$enhancer_classif<-"NA"
  CUMMULATIVE_CLASSES$enhancer_classif<-"NA"
  
  if(Condition_DEBUG == 1)
  {
    cat("CUMMULATIVE_CLASSES_1\n")
    cat(str(CUMMULATIVE_CLASSES))
    cat("\n")
    cat(str(unique(CUMMULATIVE_CLASSES$VAR)))
    cat("\n")
    cat(str(unique(CUMMULATIVE_CLASSES$KEY_Plus_carried_variants)))
    cat("\n")
  }
  
  
  
  #### Read ALL_dB file ----
  
  ALL_dB<-as.data.frame(fread(file=opt$ALL_dB,sep="\t") , stringsAsFactors=F)
  
  cat("ALL_dB_0\n")
  cat(str(ALL_dB))
  cat("\n")
  cat(str(unique(ALL_dB$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(ALL_dB$finemap_beta)))))
  cat("\n")
  cat(sprintf(as.character(summary(ALL_dB$finemap_beta))))
  cat("\n")
  
  #### Obtain  Absolute_effect_size and change finemap_prob to PP ----
  
  ALL_dB$Absolute_effect_size<-abs(ALL_dB$finemap_beta)
  colnames(ALL_dB)[which(colnames(ALL_dB) == "finemap_prob")]<-"PP"
  colnames(ALL_dB)[which(colnames(ALL_dB) == "maf_origin")]<-"MAF"
  
  ALL_dB$Allelic_Series_ID<-paste(ALL_dB$phenotype,ALL_dB$block_no,sep='__')
  
  cat("ALL_dB_1\n")
  cat(str(ALL_dB))
  cat("\n")
  cat(str(unique(ALL_dB$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(ALL_dB$Absolute_effect_size)))))
  cat("\n")
  cat(sprintf(as.character(summary(ALL_dB$Absolute_effect_size))))
  cat("\n")
  
  
  
  
  
  ### Read Supp4_Table_CURATED_PLUS_PHENOTYPES----
  
  
  Supp4_Table_CURATED_PLUS_PHENOTYPES<-as.data.frame(readRDS(file=opt$Supp4_Table_CURATED_PLUS_PHENOTYPES) , stringsAsFactors=F)
  
  cat("Supp4_Table_CURATED_PLUS_PHENOTYPES_0\n")
  cat(str(Supp4_Table_CURATED_PLUS_PHENOTYPES))
  cat("\n")
  cat(str(unique(Supp4_Table_CURATED_PLUS_PHENOTYPES$VAR)))
  cat("\n")
  
  
  screened_variants<-unique(Supp4_Table_CURATED_PLUS_PHENOTYPES$VAR)
  
  cat("screened_variants_0\n")
  cat(str(screened_variants))
  cat("\n")
  
  ### Read VAR_Prioritization_dB----
  
  
  VAR_Prioritization_dB<-as.data.frame(readRDS(file=opt$VAR_Prioritization_dB) , stringsAsFactors=F)
  
  cat("VAR_Prioritization_dB_0\n")
  cat(str(VAR_Prioritization_dB))
  cat("\n")
  cat(str(unique(VAR_Prioritization_dB$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(VAR_Prioritization_dB$bin_PP)))))
  cat("\n")
  cat(sprintf(as.character(summary(VAR_Prioritization_dB$bin_PP))))
  cat("\n")
  cat(sprintf(as.character(names(summary(VAR_Prioritization_dB$VEP_DEF_LABELS_wCSQ)))))
  cat("\n")
  cat(sprintf(as.character(summary(VAR_Prioritization_dB$VEP_DEF_LABELS_wCSQ))))
  cat("\n")
  
  VAR_Prioritization_dB_subset_1<-unique(VAR_Prioritization_dB[,c(which(colnames(VAR_Prioritization_dB) == 'VAR'),
                                                                  which(colnames(VAR_Prioritization_dB) == 'Fig1_Annot_Category'))])
  
  cat("VAR_Prioritization_dB_subset_1_0\n")
  cat(str(VAR_Prioritization_dB_subset_1))
  cat("\n")
  cat(str(unique(VAR_Prioritization_dB_subset_1$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(VAR_Prioritization_dB_subset_1$Fig1_Annot_Category)))))
  cat("\n")
  cat(sprintf(as.character(summary(VAR_Prioritization_dB_subset_1$Fig1_Annot_Category))))
  cat("\n")
  
  
  coding_wCSQ<-c("LOF","MISS","SYN","UTR5","UTR3")
  
  
  VAR_Prioritization_dB_subset_2<-VAR_Prioritization_dB[-which(VAR_Prioritization_dB$VEP_DEF_LABELS_wCSQ%in%coding_wCSQ),]
  
  
  cat("VAR_Prioritization_dB_subset_2_0\n")
  cat(str(VAR_Prioritization_dB_subset_2))
  cat("\n")
  cat(str(unique(VAR_Prioritization_dB_subset_2$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(VAR_Prioritization_dB_subset_2$bin_PP)))))
  cat("\n")
  cat(sprintf(as.character(summary(VAR_Prioritization_dB_subset_2$bin_PP))))
  cat("\n")
  cat(sprintf(as.character(names(summary(VAR_Prioritization_dB_subset_2$VEP_DEF_LABELS_wCSQ)))))
  cat("\n")
  cat(sprintf(as.character(summary(VAR_Prioritization_dB_subset_2$VEP_DEF_LABELS_wCSQ))))
  cat("\n")
  
  
  ALL_dB_subset<-ALL_dB[which(ALL_dB$VAR%in%VAR_Prioritization_dB_subset_2$VAR),]
  
  
  cat("ALL_dB_subset_0\n")
  cat(str(ALL_dB_subset))
  cat("\n")
  cat(str(unique(ALL_dB_subset$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(ALL_dB_subset$PP)))))
  cat("\n")
  cat(sprintf(as.character(summary(ALL_dB_subset$PP))))
  cat("\n")
  
  ALL_dB_subset.dt<-data.table(ALL_dB_subset, key="VAR")
  
  PP_MAXED<-as.data.frame(ALL_dB_subset.dt[,.SD[which.max(PP)],by=key(ALL_dB_subset.dt)], stringsAsFactors=F)
  
  
  cat("PP_MAXED_0\n")
  cat(str(PP_MAXED))
  cat("\n")
  cat(str(unique(PP_MAXED$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(PP_MAXED$PP)))))
  cat("\n")
  cat(sprintf(as.character(summary(PP_MAXED$PP))))
  cat("\n")
  
  PP_cuts<-c(0,0.25,0.5,0.9,1.1)
  
  PP_MAXED$bin_PP<-cut(PP_MAXED$PP,breaks = PP_cuts,right = FALSE)
  
  cat("PP_MAXED_1\n")
  cat(str(PP_MAXED))
  cat("\n")
  cat(str(unique(PP_MAXED$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(PP_MAXED$bin_PP)))))
  cat("\n")
  cat(sprintf(as.character(summary(PP_MAXED$bin_PP))))
  cat("\n")
  
  
  
  PP_MAXED$bin_PP_DEF<-revalue(PP_MAXED$bin_PP, c("[0,0.25)"='PP < 0.25',
                                                  "[0.25,0.5)"='PP < 0.5',
                                                  "[0.5,0.9)"='PP < 0.9',
                                                  "[0.9,1.1)"='PP >= 0.9'))
  
  PP_MAXED$bin_PP_DEF<-factor(PP_MAXED$bin_PP_DEF,
                              levels=c('PP < 0.25','PP < 0.5','PP < 0.9','PP >= 0.9'),
                              ordered=T)
  
  
  cat("PP_MAXED_2\n")
  cat(str(PP_MAXED))
  cat("\n")
  cat(str(unique(PP_MAXED$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(PP_MAXED$bin_PP)))))
  cat("\n")
  cat(sprintf(as.character(summary(PP_MAXED$bin_PP))))
  cat("\n")
  cat(sprintf(as.character(names(summary(PP_MAXED$bin_PP_DEF)))))
  cat("\n")
  cat(sprintf(as.character(summary(PP_MAXED$bin_PP_DEF))))
  cat("\n")
  
  
  PP_MAXED_subset<-unique(PP_MAXED[,c(which(colnames(PP_MAXED) == 'VAR'),
                                      which(colnames(PP_MAXED) == 'bin_PP_DEF'))])
  
  cat("PP_MAXED_subset_0\n")
  cat(str(PP_MAXED_subset))
  cat("\n")
  cat(str(unique(PP_MAXED_subset$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(PP_MAXED_subset$bin_PP_DEF)))))
  cat("\n")
  cat(sprintf(as.character(summary(PP_MAXED_subset$bin_PP_DEF))))
  cat("\n")
  
  
  
  list_DEF<-list()
  
  for(i in 1:length(screened_variants))
  {
    Condition_DEBUG <- 0
    
    screened_variants_sel<-screened_variants[i]
    
    cat("------------------------------------------->\t")
    cat(sprintf(as.character(i)))
    cat("\t")
    cat(sprintf(as.character(screened_variants_sel)))
    cat("\n")
    
    if(screened_variants_sel == 'chr13_28604007_T_C')
    {
      
      Condition_DEBUG <- 1
    }
    
    ALL_dB_screened_variants_sel<-ALL_dB[which(ALL_dB$VAR%in%screened_variants_sel),]
    
    if(Condition_DEBUG == 1)
    {
      cat("ALL_dB_screened_variants_sel_0\n")
      cat(str(ALL_dB_screened_variants_sel))
      cat("\n")
      cat(str(unique(ALL_dB_screened_variants_sel$VAR)))
      cat("\n")
      cat(str(unique(ALL_dB_screened_variants_sel$Allelic_Series_ID)))
      cat("\n")
    }
    
    
    AS_ID_array<-unique(ALL_dB_screened_variants_sel$Allelic_Series_ID)
    
    if(Condition_DEBUG == 1)
    {
      cat("AS_ID_array\n")
      cat(str(AS_ID_array))
      cat("\n")
    }
    
    list_AS<-list()
    
    
    for(k in 1:length(AS_ID_array))
    {
      AS_ID_array_sel<-AS_ID_array[k]
      
      cat("--->\t")
      cat(sprintf(as.character(AS_ID_array_sel)))
      cat("\n")
      
      ALL_dB_AS_sel<-ALL_dB[which(ALL_dB$Allelic_Series_ID%in%AS_ID_array_sel),]
      
      if(Condition_DEBUG == 1)
      {
        cat("ALL_dB_AS_sel_0\n")
        cat(str(ALL_dB_AS_sel))
        cat("\n")
        cat(str(unique(ALL_dB_AS_sel$VAR)))
        cat("\n")
        cat(str(unique(ALL_dB_AS_sel$Allelic_Series_ID)))
        cat("\n")
      }
      
      CUMMULATIVE_CLASSES_sel<-CUMMULATIVE_CLASSES[which(CUMMULATIVE_CLASSES$VAR%in%ALL_dB_AS_sel$VAR),]
      
      if(dim(CUMMULATIVE_CLASSES_sel)[1])
      {
        
        ##### The file includes the screened check -----
        
        check_screened_RV_included<-CUMMULATIVE_CLASSES_sel[which(CUMMULATIVE_CLASSES_sel$VAR%in%screened_variants_sel),]
        
        
        
        if(dim(check_screened_RV_included)[1] >0)
        {
          
          
          if(Condition_DEBUG == 1)
          {
            cat("CUMMULATIVE_CLASSES_sel_0\n")
            cat(str(CUMMULATIVE_CLASSES_sel))
            cat("\n")
            cat(str(unique(CUMMULATIVE_CLASSES_sel$VAR)))
            cat("\n")
          }
          
          if(Condition_DEBUG == 1)
          {
            
            
            cat("check_screened_RV_included_0\n")
            cat(str(check_screened_RV_included))
            cat("\n")
            cat(str(unique(check_screened_RV_included$VAR)))
            cat("\n")
          }
          
          #### get carried variants ----
          
          CUMMULATIVE_CLASSES_sel$carried_variants<-gsub("[^;]+;","",CUMMULATIVE_CLASSES_sel$carried_variants)
          
          if(Condition_DEBUG == 1)
          {
            cat("CUMMULATIVE_CLASSES_sel_1\n")
            cat(str(CUMMULATIVE_CLASSES_sel))
            cat("\n")
            cat(str(unique(CUMMULATIVE_CLASSES_sel$VAR)))
            cat("\n")
            cat(str(unique(CUMMULATIVE_CLASSES_sel$carried_variants)))
            cat("\n")
          }
          
          ###### MErge with CUMMULATIVE CLASSES sel -----
          
          CUMMULATIVE_CLASSES_sel<-merge(CUMMULATIVE_CLASSES_sel,
                                         VAR_Prioritization_dB_subset_1,
                                         by="VAR",
                                         all.x=T)
          if(Condition_DEBUG == 1)
          {
            cat("CUMMULATIVE_CLASSES_sel_2\n")
            cat(str(CUMMULATIVE_CLASSES_sel))
            cat("\n")
            cat(str(unique(CUMMULATIVE_CLASSES_sel$VAR)))
            cat("\n")
          }
          
          CUMMULATIVE_CLASSES_sel<-merge(CUMMULATIVE_CLASSES_sel,
                                         PP_MAXED_subset,
                                         by="VAR",
                                         all.x=T)
          if(Condition_DEBUG == 1)
          {
            cat("CUMMULATIVE_CLASSES_sel_3\n")
            cat(str(CUMMULATIVE_CLASSES_sel))
            cat("\n")
            cat(str(unique(CUMMULATIVE_CLASSES_sel$VAR)))
            cat("\n")
          }
          
          
          
          ########### Build df ---------------
          
          Cell_Type_levels<-levels(CUMMULATIVE_CLASSES_sel$Cell_Type)
          n_Cell_Type_levels<-length(levels(CUMMULATIVE_CLASSES_sel$Cell_Type))
          n_MPRA_TILES<-5
          VAR_vector<-as.character(unique(CUMMULATIVE_CLASSES_sel$VAR))
          KEY_Plus_carried_variants_vector<-as.character(unique(CUMMULATIVE_CLASSES_sel$KEY_Plus_carried_variants))
          
          
          
          n_VAR<-length(VAR_vector)
          
          if(Condition_DEBUG == 1)
          {
            cat("Parameters\n")
            cat(str(Cell_Type_levels))
            cat("\n")
            cat(str(n_Cell_Type_levels))
            cat("\n")
            cat(str(n_MPRA_TILES))
            cat("\n")
            cat(str(VAR_vector))
            cat("\n")
            cat(str(KEY_Plus_carried_variants_vector))
            cat("\n")
            
            
          }
          
          
          #### Build LogRank matrix ----
          
          
          
          Log_rank_df<-as.data.frame(cbind(rep(VAR_vector,n_MPRA_TILES),
                                           rep(KEY_Plus_carried_variants_vector,n_MPRA_TILES),
                                           c(rep(1,n_VAR),rep(2,n_VAR),rep(3,n_VAR),rep(4,n_VAR),rep(5,n_VAR)),
                                           rep("NA",n_MPRA_TILES*n_VAR)),
                                     stringsAsFactors=F)
          
          colnames(Log_rank_df)<-c("VAR","KEY_Plus_carried_variants","TILE","ACTIVE")
          
          
          if(Condition_DEBUG == 1)
          {
            cat("Log_rank_df_0\n")
            cat(str(Log_rank_df))
            cat("\n")
            cat(str(unique(Log_rank_df$VAR)))
            cat("\n")
            cat(str(unique(Log_rank_df$KEY_Plus_carried_variants)))
            cat("\n")
          }
          
          CUMMULATIVE_CLASSES_sel_subset<-unique(CUMMULATIVE_CLASSES_sel[,c(which(colnames(CUMMULATIVE_CLASSES_sel) == "VAR"),
                                                                            which(colnames(CUMMULATIVE_CLASSES_sel) == "KEY_Plus_carried_variants"),
                                                                            which(colnames(CUMMULATIVE_CLASSES_sel) == "carried_variants"),
                                                                            which(colnames(CUMMULATIVE_CLASSES_sel) == "Fig1_Annot_Category"),
                                                                            which(colnames(CUMMULATIVE_CLASSES_sel) == "bin_PP_DEF"),
                                                                            which(colnames(CUMMULATIVE_CLASSES_sel) == "Cell_Type"),
                                                                            which(colnames(CUMMULATIVE_CLASSES_sel) == "enhancer_CLASS_TILES"))])
          
          if(Condition_DEBUG == 1)
          {
            cat("CUMMULATIVE_CLASSES_sel_subset_0\n")
            cat(str(CUMMULATIVE_CLASSES_sel_subset))
            cat("\n")
            cat(str(unique(CUMMULATIVE_CLASSES_sel_subset$VAR)))
            cat("\n")
            cat(str(unique(CUMMULATIVE_CLASSES_sel_subset$KEY_Plus_carried_variants)))
            cat("\n")
            cat(sprintf(as.character(names(summary(CUMMULATIVE_CLASSES_sel_subset$Cell_Type)))))
            cat("\n")
            cat(sprintf(as.character(summary(CUMMULATIVE_CLASSES_sel_subset$Cell_Type))))
            cat("\n")
            cat(sprintf(as.character(names(summary(CUMMULATIVE_CLASSES_sel_subset$bin_PP_DEF)))))
            cat("\n")
            cat(sprintf(as.character(summary(CUMMULATIVE_CLASSES_sel_subset$bin_PP_DEF))))
            cat("\n")
          }
          
          Log_rank_df<-merge(Log_rank_df,
                             CUMMULATIVE_CLASSES_sel_subset,
                             by=c("VAR","KEY_Plus_carried_variants"),
                             all=T)
          
          
          Log_rank_df$Fig1_Annot_Category<-factor(Log_rank_df$Fig1_Annot_Category,
                                                  levels=levels(VAR_Prioritization_dB$Fig1_Annot_Category),
                                                  ordered=T)
          
          Log_rank_df$bin_PP_DEF<-factor(Log_rank_df$bin_PP_DEF,
                                         levels=levels(PP_MAXED_subset$bin_PP_DEF),
                                         ordered=T)
          
          Log_rank_df$Cell_Type<-factor(Log_rank_df$Cell_Type,
                                        levels=Cell_Type_levels,
                                        ordered=T)
          
          Log_rank_df$TILE<-as.integer(Log_rank_df$TILE)
          
          if(Condition_DEBUG == 1)
          {
            cat("Log_rank_df_1\n")
            cat(str(Log_rank_df))
            cat("\n")
            cat(str(unique(Log_rank_df$VAR)))
            cat("\n")
            cat(str(unique(Log_rank_df$KEY_Plus_carried_variants)))
            cat("\n")
          }
          
          
          
          category_df<-Log_rank_df[,c(which(colnames(Log_rank_df) == "carried_variants"),
                                      which(colnames(Log_rank_df) == "bin_PP_DEF"),
                                      which(colnames(Log_rank_df) == "Fig1_Annot_Category"))]
          
          ################# Classification of TILES ----------------
          
          Log_rank_df$ACTIVE[which(Log_rank_df$TILE <= Log_rank_df$enhancer_CLASS_TILES)]<-"1"
          Log_rank_df$ACTIVE[which(Log_rank_df$TILE > Log_rank_df$enhancer_CLASS_TILES)]<-"0"
          
          Log_rank_df$ACTIVE<-as.numeric(Log_rank_df$ACTIVE)
          
          
          
          if(Condition_DEBUG == 1)
          {
            cat("Log_rank_df_2\n")
            cat(str(Log_rank_df))
            cat("\n")
            cat(str(unique(Log_rank_df$VAR)))
            cat("\n")
            cat(sprintf(as.character(names(summary(as.factor(Log_rank_df$ACTIVE))))))
            cat("\n")
            cat(sprintf(as.character(summary(as.factor(Log_rank_df$ACTIVE)))))
            cat("\n")
          }
          
          Log_rank_df<-Log_rank_df[order(Log_rank_df$KEY_Plus_carried_variants,Log_rank_df$Cell_Type, Log_rank_df$TILE),]
          
          
          
          
          
          
          
          ############## Freq Table ----------------------------
          
          ### TOTAL ---
          
          CUMMULATIVE_CLASSES_sel.dt<-data.table(CUMMULATIVE_CLASSES_sel, key=c("Cell_Type","carried_variants"))
          
          Freq_TOTAL<-as.data.frame(CUMMULATIVE_CLASSES_sel.dt[,.(TOTAL=.N),by=key(CUMMULATIVE_CLASSES_sel.dt)], stringsAsFactors=F)
          
          if(Condition_DEBUG == 1)
          {
            cat("Freq_TOTAL_0\n")
            cat(str(Freq_TOTAL))
            cat("\n")
            #quit(status = 1)
          }
          
          Freq_TOTAL.dt<-data.table(Freq_TOTAL, key=c("carried_variants"))
          
          
          MAXED_df<-as.data.frame(Freq_TOTAL.dt[,.(MAX=max(TOTAL)),by=key(Freq_TOTAL.dt)], stringsAsFactors=F)
          
          
          
          if(Condition_DEBUG == 1)
          {
            cat("MAXED_df_0\n")
            cat(str(MAXED_df))
            cat("\n")
            #quit(status = 1)
          }
          
          ########## Active _tiles --------------------
          
          Log_rank_df.dt<-data.table(Log_rank_df, key=c("Cell_Type","carried_variants","TILE"))
          
          Freq_table<-as.data.frame(Log_rank_df.dt[,.(Freq=sum(ACTIVE)),by=key(Log_rank_df.dt)], stringsAsFactors=F)
          
          if(Condition_DEBUG == 1)
          {
            cat("Freq_table_0\n")
            cat(str(Freq_table))
            cat("\n")
            #quit(status = 1)
          }
          
          Freq_table<-merge(MAXED_df,
                            Freq_table,
                            by=c("carried_variants"),
                            all=T)
          
          if(Condition_DEBUG == 1)
          {
            cat("Freq_table_1\n")
            cat(str(Freq_table))
            cat("\n")
            #quit(status = 1)
          }
          
          Freq_table$Perc<-round(100*(Freq_table$Freq/Freq_table$MAX),2)
          
          
          if(Condition_DEBUG == 1)
          {
            cat("Freq_table_2\n")
            cat(str(Freq_table))
            cat("\n")
            #quit(status = 1)
          }
          
          #####Add categories to color the curves-----
          
          Freq_table<-merge(Freq_table,
                            category_df,
                            by="carried_variants",
                            all=T)
          
          if(Condition_DEBUG == 1)
          {
            cat("Freq_table_3\n")
            cat(str(Freq_table))
            cat("\n")
            #quit(status = 1)
          }
          
          
          Freq_table$Fig1_Annot_Category<-factor(Freq_table$Fig1_Annot_Category,
                                                 levels=levels(VAR_Prioritization_dB$Fig1_Annot_Category),
                                                 ordered=T)
          
          Freq_table$bin_PP_DEF<-factor(Freq_table$bin_PP_DEF,
                                        levels=levels(PP_MAXED_subset$bin_PP_DEF),
                                        ordered=T)
          
          Freq_table$Cell_Type<-factor(Freq_table$Cell_Type,
                                       levels=Cell_Type_levels,
                                       ordered=T)
          
          Freq_table$TILE<-as.integer(Freq_table$TILE)
          
          Freq_table<-Freq_table[order(Freq_table$bin_PP_DEF,Freq_table$Cell_Type, Freq_table$TILE),]
          
          if(Condition_DEBUG == 1)
          {
            cat("Freq_table_4\n")
            cat(str(Freq_table))
            cat("\n")
            #quit(status = 1)
          }
          
          ############### Per CT LOOP ----
          
          
          CT_levels<-Cell_Type_levels
          
          list_CT<-list()
          
          for(h in 1:length(CT_levels))
          {
            CT_sel<-CT_levels[h]
            
            cat("CT:------->\t")
            cat(sprintf(as.character(CT_sel)))
            cat("\n")
            
            Log_rank_df_sel<-Log_rank_df[which(Log_rank_df$Cell_Type == CT_sel),]
            
            if(Condition_DEBUG == 1)
            {
              cat("Log_rank_df_sel_0\n")
              cat(str(Log_rank_df_sel))
              cat("\n")
              cat(str(unique(Log_rank_df_sel$VAR)))
              cat("\n")
            }
            
            if(dim(Log_rank_df_sel)[1] >0)
            {
              #### GRAPH -----
              
              
              vector_colors<-c("greenyellow","blue","dark cyan","black","red")
              
              Freq_table_sel<-Freq_table[which(Freq_table$Cell_Type == CT_sel),]
              
              lim_Perc<-max(Freq_table_sel$Perc)
              
              if(Condition_DEBUG == 1)
              {
                cat("Freq_table_sel_0\n")
                cat(str(Freq_table_sel))
                cat("\n")
                cat("lim_Perc\n")
                cat(str(lim_Perc))
                cat("\n")
                
              }
              
              
              
              
              path5<-paste(out2,'CT_AS_printing_CUMULATIVE_CURVES','/', sep='')
              
              if (file.exists(path5)){
                
                
                
                
              } else {
                dir.create(file.path(path5))
                
              }
              
              path6<-paste(path5,screened_variants_sel,'/', sep='')
              
              if (file.exists(path6)){
                
                
                
                
              } else {
                dir.create(file.path(path6))
                
              }
              
              path7<-paste(path6,AS_ID_array_sel,'/', sep='')
              
              if (file.exists(path7)){
                
                
                
                
              } else {
                dir.create(file.path(path7))
                
              }
              
              setwd(path7)
              
              
              #### G1
              
              v_parameter<-"NA" 
              
              indx_carried_variants<-which(colnames(Freq_table_sel) == "carried_variants")
              indx_category<-which(colnames(Freq_table_sel) == "Fig1_Annot_Category")
              
              indx.yaxis<-which(colnames(Freq_table_sel) == "Perc")
              mycols <- vector_colors
              
              pdfname<-paste(CT_sel,"_Cumulative_frequency_by_Fig1_cat_","enhancer_activity",".pdf",sep='')
              pdf(file=pdfname, width=5, height=4, pointsize=12)
              
              par(mai=c(0.9,0.9,0.3,0.2))
              lab <- as.character(unique(Freq_table_sel[,indx_carried_variants]))
              
              if(Condition_DEBUG == 1)
              {
                cat("lab\n")
                cat(sprintf(as.character(lab)))
                cat("\n")
              }
              
              
              plot(Freq_table_sel$TILE, Freq_table_sel[,indx.yaxis],
                   ty="n", xlab="TILES with enhancer activity",
                   ylab="Cumulative % of variants",
                   axes=F, cex.lab=1.2, cex.lab=1.3, ylim=c(0,lim_Perc))
              
              
              
              # cat("Hello_world1\n")
              
              # points(Freq_table_sel$TILE, Freq_table_sel[,indx.yaxis], col="darkgrey", pch=19)
              
              # cat("Hello_world2\n")
              # 
              # 
              for (iteration_graph in 1:length(lab))  {
                
                if(Condition_DEBUG == 1)
                {
                  cat("Hello_world3\t")
                  cat(sprintf(as.character(lab[iteration_graph])))
                  cat("\n")
                }
                
                ind <- which(Freq_table_sel[,indx_carried_variants]==lab[iteration_graph])
                
                if(Condition_DEBUG == 1)
                {
                  cat("------------------->ind\n")
                  cat(str(ind))
                  cat("\n")
                }
                
                category_sel<-unique(Freq_table_sel$Fig1_Annot_Category[ind])
                
                position_vector_sel<-as.numeric(category_sel)
                
                color_sel<-mycols[position_vector_sel]#######################################
                
                if(Condition_DEBUG == 1)
                {
                  cat(sprintf(as.character(color_sel)))
                  cat("\n")
                }
                
                points(Freq_table_sel$TILE[ind], Freq_table_sel[,indx.yaxis][ind],
                       pch=21,lwd=5, col=color_sel)
                lines(Freq_table_sel$TILE[ind], Freq_table_sel[,indx.yaxis][ind],
                      lty=1,lwd=5, col=color_sel)
                
                if(Condition_DEBUG == 1)
                {
                  cat("END_3\n")
                }
                
                # points(Freq_table_sel$TILE[ind], Freq_table_sel[,indx.yaxis][ind], pch=19, col=color_sel)
              }#iteration_graph in 1:length(lab)
              
              if(Condition_DEBUG == 1)
              {
                cat("END_4\n")
              }
              
              
              legend("topright", legend=lab, fill=mycols, border=mycols, bty="n")
              axis(1, at=seq(0,100))
              axis(2, las=1)
              
              dev.off()
              
              #### G2
              
              v_parameter<-"NA" 
              
              indx_carried_variants<-which(colnames(Freq_table_sel) == "carried_variants")
              indx_category<-which(colnames(Freq_table_sel) == "bin_PP_DEF")
              
              indx.yaxis<-which(colnames(Freq_table_sel) == "Perc")
              mycols <- vector_colors
              
              pdfname<-paste(CT_sel,"_Cumulative_frequency_by_bin_PP_cat_","enhancer_activity",".pdf",sep='')
              pdf(file=pdfname, width=5, height=4, pointsize=12)
              
              par(mai=c(0.9,0.9,0.3,0.2))
              lab <- as.character(unique(Freq_table_sel[,indx_carried_variants]))
              
              if(Condition_DEBUG == 1)
              {
                cat("lab\n")
                cat(sprintf(as.character(lab)))
                cat("\n")
              }
              
              
              plot(Freq_table_sel$TILE, Freq_table_sel[,indx.yaxis],
                   ty="n", xlab="TILES with enhancer activity",
                   ylab="Cumulative % of variants",
                   axes=F, cex.lab=1.2, cex.lab=1.3, ylim=c(0,lim_Perc))
              
              
              
              # cat("Hello_world1\n")
              
              # points(Freq_table_sel$TILE, Freq_table_sel[,indx.yaxis], col="darkgrey", pch=19)
              
              # cat("Hello_world2\n")
              # 
              # 
              for (iteration_graph in 1:length(lab))  {
                
                if(Condition_DEBUG == 1)
                {
                  cat("Hello_world3\t")
                  cat(sprintf(as.character(lab[iteration_graph])))
                  cat("\n")
                }
                
                ind <- which(Freq_table_sel[,indx_carried_variants]==lab[iteration_graph])
                
                if(Condition_DEBUG == 1)
                {
                  cat("------------------->ind\n")
                  cat(str(ind))
                  cat("\n")
                }
                
                category_sel<-unique(Freq_table_sel$bin_PP_DEF[ind])
                
                position_vector_sel<-as.numeric(category_sel)
                
                color_sel<-mycols[position_vector_sel]#######################################
                
                if(Condition_DEBUG == 1)
                {
                  cat(sprintf(as.character(color_sel)))
                  cat("\n")
                }
                
                points(Freq_table_sel$TILE[ind], Freq_table_sel[,indx.yaxis][ind],
                       pch=21,lwd=5, col=color_sel)
                lines(Freq_table_sel$TILE[ind], Freq_table_sel[,indx.yaxis][ind],
                      lty=1,lwd=5, col=color_sel)
                
                if(Condition_DEBUG == 1)
                {
                  cat("END_3\n")
                }
                
                # points(Freq_table_sel$TILE[ind], Freq_table_sel[,indx.yaxis][ind], pch=19, col=color_sel)
              }#iteration_graph in 1:length(lab)
              
              if(Condition_DEBUG == 1)
              {
                cat("END_4\n")
              }
              
              
              legend("topright", legend=lab, fill=mycols, border=mycols, bty="n")
              axis(1, at=seq(0,100))
              axis(2, las=1)
              
              dev.off()
              
              
              
              list_CT[[h]]<-Freq_table_sel
              
              
              
              
              
              
            }#dim(Log_rank_df_sel)[1] >0
          }#h in 1:length(CT_levels)
          
          
          if(length(list_CT) >0)
          {
            CT_df_recovered = unique(as.data.frame(data.table::rbindlist(list_CT, fill=T), stringsAsFactors=F))
            
            
            if(Condition_DEBUG == 1)
            {
              cat("CT_df_recovered_0\n")
              cat(str(CT_df_recovered))
              cat("\n")
              #quit(status = 1)
            }
            
            CT_df_recovered$Allelic_Series_ID<-AS_ID_array_sel
            
            if(Condition_DEBUG == 1)
            {
              cat("CT_df_recovered_1\n")
              cat(str(CT_df_recovered))
              cat("\n")
              #quit(status = 1)
            }
            
            list_AS[[k]]<-CT_df_recovered
          }#length(list_CT) >0
          
          
        }#dim(check_screened_RV_included)[1] >0
      }#dim(CUMMULATIVE_CLASSES_sel)[1]
    }#k in 1:length(AS_ID_array)
    
    
    
    if(length(list_AS) >0)
    {
      AS_df_recovered = unique(as.data.frame(data.table::rbindlist(list_AS, fill=T), stringsAsFactors=F))
      
      
      if(Condition_DEBUG == 1)
      {
        cat("AS_df_recovered_0\n")
        cat(str(AS_df_recovered))
        cat("\n")
        #quit(status = 1)
      }
      
      AS_df_recovered$screened_variant<-screened_variants_sel
      
      if(Condition_DEBUG == 1)
      {
        cat("AS_df_recovered_1\n")
        cat(str(AS_df_recovered))
        cat("\n")
        #quit(status = 1)
      }
      
      list_DEF[[i]]<-AS_df_recovered
    }#length(list_AS) >0
    # #####################################################################################
    # quit(status = 1)
    
  }# i in 1:length(screened_variants)
  
  if(length(list_DEF) >0)
  {
    DEF_df = unique(as.data.frame(data.table::rbindlist(list_DEF, fill=T), stringsAsFactors=F))
    
    
    if(Condition_DEBUG == 1)
    {
      cat("DEF_df_0\n")
      cat(str(DEF_df))
      cat("\n")
      cat(str(unique(DEF_df$screened_variant)))
      cat("\n")
      cat(str(unique(DEF_df$carried_variants)))
      cat("\n")
      #quit(status = 1)
    }
    
    path5<-paste(out2,'CT_AS_printing_CUMULATIVE_CURVES','/', sep='')
    
    setwd(path5)
    
    write.table(DEF_df, file="enhancer_AS_Printing_Cummulative_Freq_table.tsv", sep="\t", quote=F,row.names = F)
    
    
  }#length(list_DEF) >0
  
  
}

cumulative_CT_Label2_and_Log_Rank_test_E_Plus_ASE = function(option_list)
{
  suppressMessages(library("nph", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
  
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
  
  #### READ and transform out2 ----
  
  out2 = opt$out2
  
  cat("out2_\n")
  cat(sprintf(as.character(out2)))
  cat("\n")
  
  setwd(out)
  
  
  filename<-paste('df_Cell_colors','.rds',sep='')
  
  df_Cell_colors<-readRDS(file=filename)
  
  cat("df_Cell_colors_0\n")
  cat(str(df_Cell_colors))
  cat("\n")
  
  Condition_DEBUG <- 0
  
  ### CUMMULATIVE_CLASSES #----
  
  CUMMULATIVE_CLASSES<-readRDS(file=opt$CUMMULATIVE_CLASSES)
  
  CUMMULATIVE_CLASSES<-unique(CUMMULATIVE_CLASSES[,-c(which(colnames(CUMMULATIVE_CLASSES) == "VEP_DEF_LABELS"),
                                                      which(colnames(CUMMULATIVE_CLASSES) == "factor4_CLASS"))])
  
  
  cat("CUMMULATIVE_CLASSES_0\n")
  cat(str(CUMMULATIVE_CLASSES))
  cat("\n")
  cat(str(unique(CUMMULATIVE_CLASSES$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(CUMMULATIVE_CLASSES$bin_PP)))))
  cat("\n")
  cat(sprintf(as.character(summary(CUMMULATIVE_CLASSES$bin_PP))))
  cat("\n")
  
 
  
  # CUMMULATIVE_CLASSES<-droplevels(CUMMULATIVE_CLASSES[which(CUMMULATIVE_CLASSES$Cell_Type != "ALL_CT"),]) # Keep "ALL_CT"
  
  CUMMULATIVE_CLASSES$enhancer_classif<-"NA"
  CUMMULATIVE_CLASSES$E_Plus_ASE_classif<-"NA"
  
  if(Condition_DEBUG == 1)
  {
    cat("CUMMULATIVE_CLASSES_1\n")
    cat(str(CUMMULATIVE_CLASSES))
    cat("\n")
    cat(str(unique(CUMMULATIVE_CLASSES$VAR)))
    cat("\n")
    cat(str(unique(CUMMULATIVE_CLASSES$KEY_Plus_carried_variants)))
    cat("\n")
  }
  

  
  #### Read ALL_dB file ----
  
  ALL_dB<-as.data.frame(fread(file=opt$ALL_dB,sep="\t") , stringsAsFactors=F)
  
  cat("ALL_dB_0\n")
  cat(str(ALL_dB))
  cat("\n")
  cat(str(unique(ALL_dB$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(ALL_dB$finemap_beta)))))
  cat("\n")
  cat(sprintf(as.character(summary(ALL_dB$finemap_beta))))
  cat("\n")
  
  #### Obtain  Absolute_effect_size and change finemap_prob to PP ----
  
  ALL_dB$Absolute_effect_size<-abs(ALL_dB$finemap_beta)
  colnames(ALL_dB)[which(colnames(ALL_dB) == "finemap_prob")]<-"PP"
  colnames(ALL_dB)[which(colnames(ALL_dB) == "maf_origin")]<-"MAF"
  
  ALL_dB$Allelic_Series_ID<-paste(ALL_dB$phenotype,ALL_dB$block_no,sep='__')
                                         
  cat("ALL_dB_1\n")
  cat(str(ALL_dB))
  cat("\n")
  cat(str(unique(ALL_dB$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(ALL_dB$Absolute_effect_size)))))
  cat("\n")
  cat(sprintf(as.character(summary(ALL_dB$Absolute_effect_size))))
  cat("\n")
  
  
  
  
  
  ### Read Supp4_Table_CURATED_PLUS_PHENOTYPES----
  
  
  Supp4_Table_CURATED_PLUS_PHENOTYPES<-as.data.frame(readRDS(file=opt$Supp4_Table_CURATED_PLUS_PHENOTYPES) , stringsAsFactors=F)
  
  cat("Supp4_Table_CURATED_PLUS_PHENOTYPES_0\n")
  cat(str(Supp4_Table_CURATED_PLUS_PHENOTYPES))
  cat("\n")
  cat(str(unique(Supp4_Table_CURATED_PLUS_PHENOTYPES$VAR)))
  cat("\n")
  
  
  screened_variants<-unique(Supp4_Table_CURATED_PLUS_PHENOTYPES$VAR)
  
  cat("screened_variants_0\n")
  cat(str(screened_variants))
  cat("\n")
  
  ### Read VAR_Prioritization_dB----
  
  
  VAR_Prioritization_dB<-as.data.frame(readRDS(file=opt$VAR_Prioritization_dB) , stringsAsFactors=F)
  
  cat("VAR_Prioritization_dB_0\n")
  cat(str(VAR_Prioritization_dB))
  cat("\n")
  cat(str(unique(VAR_Prioritization_dB$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(VAR_Prioritization_dB$bin_PP)))))
  cat("\n")
  cat(sprintf(as.character(summary(VAR_Prioritization_dB$bin_PP))))
  cat("\n")
  cat(sprintf(as.character(names(summary(VAR_Prioritization_dB$VEP_DEF_LABELS_wCSQ)))))
  cat("\n")
  cat(sprintf(as.character(summary(VAR_Prioritization_dB$VEP_DEF_LABELS_wCSQ))))
  cat("\n")
  
  VAR_Prioritization_dB_subset_1<-unique(VAR_Prioritization_dB[,c(which(colnames(VAR_Prioritization_dB) == 'VAR'),
                                                                which(colnames(VAR_Prioritization_dB) == 'Fig1_Annot_Category'))])
  
  cat("VAR_Prioritization_dB_subset_1_0\n")
  cat(str(VAR_Prioritization_dB_subset_1))
  cat("\n")
  cat(str(unique(VAR_Prioritization_dB_subset_1$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(VAR_Prioritization_dB_subset_1$Fig1_Annot_Category)))))
  cat("\n")
  cat(sprintf(as.character(summary(VAR_Prioritization_dB_subset_1$Fig1_Annot_Category))))
  cat("\n")
  
  
  coding_wCSQ<-c("LOF","MISS","SYN","UTR5","UTR3")
  
  
  VAR_Prioritization_dB_subset_2<-VAR_Prioritization_dB[-which(VAR_Prioritization_dB$VEP_DEF_LABELS_wCSQ%in%coding_wCSQ),]
  
  
  cat("VAR_Prioritization_dB_subset_2_0\n")
  cat(str(VAR_Prioritization_dB_subset_2))
  cat("\n")
  cat(str(unique(VAR_Prioritization_dB_subset_2$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(VAR_Prioritization_dB_subset_2$bin_PP)))))
  cat("\n")
  cat(sprintf(as.character(summary(VAR_Prioritization_dB_subset_2$bin_PP))))
  cat("\n")
  cat(sprintf(as.character(names(summary(VAR_Prioritization_dB_subset_2$VEP_DEF_LABELS_wCSQ)))))
  cat("\n")
  cat(sprintf(as.character(summary(VAR_Prioritization_dB_subset_2$VEP_DEF_LABELS_wCSQ))))
  cat("\n")
  
  
  ALL_dB_subset<-ALL_dB[which(ALL_dB$VAR%in%VAR_Prioritization_dB_subset_2$VAR),]
  
  
  cat("ALL_dB_subset_0\n")
  cat(str(ALL_dB_subset))
  cat("\n")
  cat(str(unique(ALL_dB_subset$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(ALL_dB_subset$PP)))))
  cat("\n")
  cat(sprintf(as.character(summary(ALL_dB_subset$PP))))
  cat("\n")
  
  ALL_dB_subset.dt<-data.table(ALL_dB_subset, key="VAR")
  
  PP_MAXED<-as.data.frame(ALL_dB_subset.dt[,.SD[which.max(PP)],by=key(ALL_dB_subset.dt)], stringsAsFactors=F)
  
  
  cat("PP_MAXED_0\n")
  cat(str(PP_MAXED))
  cat("\n")
  cat(str(unique(PP_MAXED$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(PP_MAXED$PP)))))
  cat("\n")
  cat(sprintf(as.character(summary(PP_MAXED$PP))))
  cat("\n")
  
  PP_cuts<-c(0,0.25,0.5,0.9,1.1)
  
  PP_MAXED$bin_PP<-cut(PP_MAXED$PP,breaks = PP_cuts,right = FALSE)
  
  cat("PP_MAXED_1\n")
  cat(str(PP_MAXED))
  cat("\n")
  cat(str(unique(PP_MAXED$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(PP_MAXED$bin_PP)))))
  cat("\n")
  cat(sprintf(as.character(summary(PP_MAXED$bin_PP))))
  cat("\n")
  
  
  
  PP_MAXED$bin_PP_DEF<-revalue(PP_MAXED$bin_PP, c("[0,0.25)"='PP < 0.25',
                                                  "[0.25,0.5)"='PP < 0.5',
                                                  "[0.5,0.9)"='PP < 0.9',
                                                  "[0.9,1.1)"='PP >= 0.9'))
  
  PP_MAXED$bin_PP_DEF<-factor(PP_MAXED$bin_PP_DEF,
                              levels=c('PP < 0.25','PP < 0.5','PP < 0.9','PP >= 0.9'),
                              ordered=T)
  
  
  cat("PP_MAXED_2\n")
  cat(str(PP_MAXED))
  cat("\n")
  cat(str(unique(PP_MAXED$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(PP_MAXED$bin_PP)))))
  cat("\n")
  cat(sprintf(as.character(summary(PP_MAXED$bin_PP))))
  cat("\n")
  cat(sprintf(as.character(names(summary(PP_MAXED$bin_PP_DEF)))))
  cat("\n")
  cat(sprintf(as.character(summary(PP_MAXED$bin_PP_DEF))))
  cat("\n")
  
  
  PP_MAXED_subset<-unique(PP_MAXED[,c(which(colnames(PP_MAXED) == 'VAR'),
                                      which(colnames(PP_MAXED) == 'bin_PP_DEF'))])
  
  cat("PP_MAXED_subset_0\n")
  cat(str(PP_MAXED_subset))
  cat("\n")
  cat(str(unique(PP_MAXED_subset$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(PP_MAXED_subset$bin_PP_DEF)))))
  cat("\n")
  cat(sprintf(as.character(summary(PP_MAXED_subset$bin_PP_DEF))))
  cat("\n")
  
 
  
  list_DEF<-list()
  
  for(i in 1:length(screened_variants))
  {
    Condition_DEBUG <- 0
    
    screened_variants_sel<-screened_variants[i]
    
    cat("------------------------------------------->\t")
    cat(sprintf(as.character(i)))
    cat("\t")
    cat(sprintf(as.character(screened_variants_sel)))
    cat("\n")
    
    if(screened_variants_sel == 'chr13_28604007_T_C')
    {
      
      Condition_DEBUG <- 1
    }
    
    ALL_dB_screened_variants_sel<-ALL_dB[which(ALL_dB$VAR%in%screened_variants_sel),]
    
    if(Condition_DEBUG == 1)
    {
      cat("ALL_dB_screened_variants_sel_0\n")
      cat(str(ALL_dB_screened_variants_sel))
      cat("\n")
      cat(str(unique(ALL_dB_screened_variants_sel$VAR)))
      cat("\n")
      cat(str(unique(ALL_dB_screened_variants_sel$Allelic_Series_ID)))
      cat("\n")
    }
    
    
    AS_ID_array<-unique(ALL_dB_screened_variants_sel$Allelic_Series_ID)
    
    if(Condition_DEBUG == 1)
    {
      cat("AS_ID_array\n")
      cat(str(AS_ID_array))
      cat("\n")
    }
    
    list_AS<-list()
    
    
    for(k in 1:length(AS_ID_array))
    {
      AS_ID_array_sel<-AS_ID_array[k]
      
      cat("--->\t")
      cat(sprintf(as.character(AS_ID_array_sel)))
      cat("\n")
      
      ALL_dB_AS_sel<-ALL_dB[which(ALL_dB$Allelic_Series_ID%in%AS_ID_array_sel),]
      
      if(Condition_DEBUG == 1)
      {
        cat("ALL_dB_AS_sel_0\n")
        cat(str(ALL_dB_AS_sel))
        cat("\n")
        cat(str(unique(ALL_dB_AS_sel$VAR)))
        cat("\n")
        cat(str(unique(ALL_dB_AS_sel$Allelic_Series_ID)))
        cat("\n")
      }
      
      CUMMULATIVE_CLASSES_sel<-CUMMULATIVE_CLASSES[which(CUMMULATIVE_CLASSES$VAR%in%ALL_dB_AS_sel$VAR),]
      
      if(dim(CUMMULATIVE_CLASSES_sel)[1])
      {
        
        ##### The file includes the screened check -----
        
        check_screened_RV_included<-CUMMULATIVE_CLASSES_sel[which(CUMMULATIVE_CLASSES_sel$VAR%in%screened_variants_sel),]
        
     
        
        if(dim(check_screened_RV_included)[1] >0)
        {
          
          
          if(Condition_DEBUG == 1)
          {
            cat("CUMMULATIVE_CLASSES_sel_0\n")
            cat(str(CUMMULATIVE_CLASSES_sel))
            cat("\n")
            cat(str(unique(CUMMULATIVE_CLASSES_sel$VAR)))
            cat("\n")
          }
          
          if(Condition_DEBUG == 1)
          {
            
            
            cat("check_screened_RV_included_0\n")
            cat(str(check_screened_RV_included))
            cat("\n")
            cat(str(unique(check_screened_RV_included$VAR)))
            cat("\n")
          }
          
          #### get carried variants ----
          
          CUMMULATIVE_CLASSES_sel$carried_variants<-gsub("[^;]+;","",CUMMULATIVE_CLASSES_sel$carried_variants)
          
          if(Condition_DEBUG == 1)
          {
            cat("CUMMULATIVE_CLASSES_sel_1\n")
            cat(str(CUMMULATIVE_CLASSES_sel))
            cat("\n")
            cat(str(unique(CUMMULATIVE_CLASSES_sel$VAR)))
            cat("\n")
            cat(str(unique(CUMMULATIVE_CLASSES_sel$carried_variants)))
            cat("\n")
          }
          
          ###### MErge with CUMMULATIVE CLASSES sel -----
          
          CUMMULATIVE_CLASSES_sel<-merge(CUMMULATIVE_CLASSES_sel,
                                         VAR_Prioritization_dB_subset_1,
                                         by="VAR",
                                         all.x=T)
          if(Condition_DEBUG == 1)
          {
            cat("CUMMULATIVE_CLASSES_sel_2\n")
            cat(str(CUMMULATIVE_CLASSES_sel))
            cat("\n")
            cat(str(unique(CUMMULATIVE_CLASSES_sel$VAR)))
            cat("\n")
          }
          
          CUMMULATIVE_CLASSES_sel<-merge(CUMMULATIVE_CLASSES_sel,
                                         PP_MAXED_subset,
                                         by="VAR",
                                         all.x=T)
          if(Condition_DEBUG == 1)
          {
            cat("CUMMULATIVE_CLASSES_sel_3\n")
            cat(str(CUMMULATIVE_CLASSES_sel))
            cat("\n")
            cat(str(unique(CUMMULATIVE_CLASSES_sel$VAR)))
            cat("\n")
          }
          
          
        
          ########### Build df ---------------
          
          Cell_Type_levels<-levels(CUMMULATIVE_CLASSES_sel$Cell_Type)
          n_Cell_Type_levels<-length(levels(CUMMULATIVE_CLASSES_sel$Cell_Type))
          n_MPRA_TILES<-5
          VAR_vector<-as.character(unique(CUMMULATIVE_CLASSES_sel$VAR))
          KEY_Plus_carried_variants_vector<-as.character(unique(CUMMULATIVE_CLASSES_sel$KEY_Plus_carried_variants))
          
          
          
          n_VAR<-length(VAR_vector)
          
          if(Condition_DEBUG == 1)
          {
            cat("Parameters\n")
            cat(str(Cell_Type_levels))
            cat("\n")
            cat(str(n_Cell_Type_levels))
            cat("\n")
            cat(str(n_MPRA_TILES))
            cat("\n")
            cat(str(VAR_vector))
            cat("\n")
            cat(str(KEY_Plus_carried_variants_vector))
            cat("\n")
            
            
          }
          
          
          #### Build LogRank matrix ----
          
          
          
          Log_rank_df<-as.data.frame(cbind(rep(VAR_vector,n_MPRA_TILES),
                                           rep(KEY_Plus_carried_variants_vector,n_MPRA_TILES),
                                           c(rep(1,n_VAR),rep(2,n_VAR),rep(3,n_VAR),rep(4,n_VAR),rep(5,n_VAR)),
                                           rep("NA",n_MPRA_TILES*n_VAR)),
                                     stringsAsFactors=F)
          
          colnames(Log_rank_df)<-c("VAR","KEY_Plus_carried_variants","TILE","ACTIVE")
          
          
          if(Condition_DEBUG == 1)
          {
            cat("Log_rank_df_0\n")
            cat(str(Log_rank_df))
            cat("\n")
            cat(str(unique(Log_rank_df$VAR)))
            cat("\n")
            cat(str(unique(Log_rank_df$KEY_Plus_carried_variants)))
            cat("\n")
          }
          
          CUMMULATIVE_CLASSES_sel_subset<-unique(CUMMULATIVE_CLASSES_sel[,c(which(colnames(CUMMULATIVE_CLASSES_sel) == "VAR"),
                                                                            which(colnames(CUMMULATIVE_CLASSES_sel) == "KEY_Plus_carried_variants"),
                                                                            which(colnames(CUMMULATIVE_CLASSES_sel) == "carried_variants"),
                                                                            which(colnames(CUMMULATIVE_CLASSES_sel) == "Fig1_Annot_Category"),
                                                                            which(colnames(CUMMULATIVE_CLASSES_sel) == "bin_PP_DEF"),
                                                                            which(colnames(CUMMULATIVE_CLASSES_sel) == "Cell_Type"),
                                                                            which(colnames(CUMMULATIVE_CLASSES_sel) == "E_Plus_ASE_CLASS_TILES"))])
          
          if(Condition_DEBUG == 1)
          {
            cat("CUMMULATIVE_CLASSES_sel_subset_0\n")
            cat(str(CUMMULATIVE_CLASSES_sel_subset))
            cat("\n")
            cat(str(unique(CUMMULATIVE_CLASSES_sel_subset$VAR)))
            cat("\n")
            cat(str(unique(CUMMULATIVE_CLASSES_sel_subset$KEY_Plus_carried_variants)))
            cat("\n")
            cat(sprintf(as.character(names(summary(CUMMULATIVE_CLASSES_sel_subset$Cell_Type)))))
            cat("\n")
            cat(sprintf(as.character(summary(CUMMULATIVE_CLASSES_sel_subset$Cell_Type))))
            cat("\n")
            cat(sprintf(as.character(names(summary(CUMMULATIVE_CLASSES_sel_subset$bin_PP_DEF)))))
            cat("\n")
            cat(sprintf(as.character(summary(CUMMULATIVE_CLASSES_sel_subset$bin_PP_DEF))))
            cat("\n")
          }
          
          Log_rank_df<-merge(Log_rank_df,
                             CUMMULATIVE_CLASSES_sel_subset,
                             by=c("VAR","KEY_Plus_carried_variants"),
                             all=T)
          

          Log_rank_df$Fig1_Annot_Category<-factor(Log_rank_df$Fig1_Annot_Category,
                                                  levels=levels(VAR_Prioritization_dB$Fig1_Annot_Category),
                                                  ordered=T)
          
          Log_rank_df$bin_PP_DEF<-factor(Log_rank_df$bin_PP_DEF,
                                         levels=levels(PP_MAXED_subset$bin_PP_DEF),
                                         ordered=T)
          
          Log_rank_df$Cell_Type<-factor(Log_rank_df$Cell_Type,
                                        levels=Cell_Type_levels,
                                        ordered=T)
          
          Log_rank_df$TILE<-as.integer(Log_rank_df$TILE)
          
          if(Condition_DEBUG == 1)
          {
            cat("Log_rank_df_1\n")
            cat(str(Log_rank_df))
            cat("\n")
            cat(str(unique(Log_rank_df$VAR)))
            cat("\n")
            cat(str(unique(Log_rank_df$KEY_Plus_carried_variants)))
            cat("\n")
          }
          
          
          
          category_df<-Log_rank_df[,c(which(colnames(Log_rank_df) == "carried_variants"),
                                      which(colnames(Log_rank_df) == "bin_PP_DEF"),
                                      which(colnames(Log_rank_df) == "Fig1_Annot_Category"))]
          
          ################# Classification of TILES ----------------
          
          Log_rank_df$ACTIVE[which(Log_rank_df$TILE <= Log_rank_df$E_Plus_ASE_CLASS_TILES)]<-"1"
          Log_rank_df$ACTIVE[which(Log_rank_df$TILE > Log_rank_df$E_Plus_ASE_CLASS_TILES)]<-"0"
          
          Log_rank_df$ACTIVE<-as.numeric(Log_rank_df$ACTIVE)
          
          
          
          if(Condition_DEBUG == 1)
          {
            cat("Log_rank_df_2\n")
            cat(str(Log_rank_df))
            cat("\n")
            cat(str(unique(Log_rank_df$VAR)))
            cat("\n")
            cat(sprintf(as.character(names(summary(as.factor(Log_rank_df$ACTIVE))))))
            cat("\n")
            cat(sprintf(as.character(summary(as.factor(Log_rank_df$ACTIVE)))))
            cat("\n")
          }
          
          Log_rank_df<-Log_rank_df[order(Log_rank_df$KEY_Plus_carried_variants,Log_rank_df$Cell_Type, Log_rank_df$TILE),]
          
          
         
          
         
          
          
          ############## Freq Table ----------------------------
          
          ### TOTAL ---
          
          CUMMULATIVE_CLASSES_sel.dt<-data.table(CUMMULATIVE_CLASSES_sel, key=c("Cell_Type","carried_variants"))
          
          Freq_TOTAL<-as.data.frame(CUMMULATIVE_CLASSES_sel.dt[,.(TOTAL=.N),by=key(CUMMULATIVE_CLASSES_sel.dt)], stringsAsFactors=F)
          
          if(Condition_DEBUG == 1)
          {
            cat("Freq_TOTAL_0\n")
            cat(str(Freq_TOTAL))
            cat("\n")
            #quit(status = 1)
          }
          
          Freq_TOTAL.dt<-data.table(Freq_TOTAL, key=c("carried_variants"))
          
          
          MAXED_df<-as.data.frame(Freq_TOTAL.dt[,.(MAX=max(TOTAL)),by=key(Freq_TOTAL.dt)], stringsAsFactors=F)
          
          
          
          if(Condition_DEBUG == 1)
          {
            cat("MAXED_df_0\n")
            cat(str(MAXED_df))
            cat("\n")
            #quit(status = 1)
          }
          
          ########## Active _tiles --------------------
          
          Log_rank_df.dt<-data.table(Log_rank_df, key=c("Cell_Type","carried_variants","TILE"))
          
          Freq_table<-as.data.frame(Log_rank_df.dt[,.(Freq=sum(ACTIVE)),by=key(Log_rank_df.dt)], stringsAsFactors=F)
          
          if(Condition_DEBUG == 1)
          {
            cat("Freq_table_0\n")
            cat(str(Freq_table))
            cat("\n")
            #quit(status = 1)
          }
          
          Freq_table<-merge(MAXED_df,
                            Freq_table,
                            by=c("carried_variants"),
                            all=T)
          
          if(Condition_DEBUG == 1)
          {
            cat("Freq_table_1\n")
            cat(str(Freq_table))
            cat("\n")
            #quit(status = 1)
          }
          
          Freq_table$Perc<-round(100*(Freq_table$Freq/Freq_table$MAX),2)
          
          
          if(Condition_DEBUG == 1)
          {
            cat("Freq_table_2\n")
            cat(str(Freq_table))
            cat("\n")
            #quit(status = 1)
          }
          
          #####Add categories to color the curves-----
          
          Freq_table<-merge(Freq_table,
                            category_df,
                            by="carried_variants",
                            all=T)
          
          if(Condition_DEBUG == 1)
          {
            cat("Freq_table_3\n")
            cat(str(Freq_table))
            cat("\n")
            #quit(status = 1)
          }
          
         
          Freq_table$Fig1_Annot_Category<-factor(Freq_table$Fig1_Annot_Category,
                                                  levels=levels(VAR_Prioritization_dB$Fig1_Annot_Category),
                                                  ordered=T)
          
          Freq_table$bin_PP_DEF<-factor(Freq_table$bin_PP_DEF,
                                         levels=levels(PP_MAXED_subset$bin_PP_DEF),
                                         ordered=T)
          
          Freq_table$Cell_Type<-factor(Freq_table$Cell_Type,
                                        levels=Cell_Type_levels,
                                        ordered=T)
          
          Freq_table$TILE<-as.integer(Freq_table$TILE)
          
          Freq_table<-Freq_table[order(Freq_table$bin_PP_DEF,Freq_table$Cell_Type, Freq_table$TILE),]
          
          if(Condition_DEBUG == 1)
          {
            cat("Freq_table_4\n")
            cat(str(Freq_table))
            cat("\n")
            #quit(status = 1)
          }
          
          ############### Per CT LOOP ----
          
          
          CT_levels<-Cell_Type_levels
          
          list_CT<-list()
          
          for(h in 1:length(CT_levels))
          {
            CT_sel<-CT_levels[h]
            
            cat("CT:------->\t")
            cat(sprintf(as.character(CT_sel)))
            cat("\n")
            
            Log_rank_df_sel<-Log_rank_df[which(Log_rank_df$Cell_Type == CT_sel),]
            
            if(Condition_DEBUG == 1)
            {
              cat("Log_rank_df_sel_0\n")
              cat(str(Log_rank_df_sel))
              cat("\n")
              cat(str(unique(Log_rank_df_sel$VAR)))
              cat("\n")
            }
            
            if(dim(Log_rank_df_sel)[1] >0)
            {
              #### GRAPH -----
              
              
              vector_colors<-c("greenyellow","blue","dark cyan","black","red")
              
              Freq_table_sel<-Freq_table[which(Freq_table$Cell_Type == CT_sel),]
              
              lim_Perc<-max(Freq_table_sel$Perc)
              
              if(Condition_DEBUG == 1)
              {
                cat("Freq_table_sel_0\n")
                cat(str(Freq_table_sel))
                cat("\n")
                cat("lim_Perc\n")
                cat(str(lim_Perc))
                cat("\n")
                
              }
              
              
              
              
              path5<-paste(out2,'CT_AS_printing_CUMULATIVE_CURVES','/', sep='')
              
              if (file.exists(path5)){
                
                
                
                
              } else {
                dir.create(file.path(path5))
                
              }
              
              path6<-paste(path5,screened_variants_sel,'/', sep='')
              
              if (file.exists(path6)){
                
                
                
                
              } else {
                dir.create(file.path(path6))
                
              }
              
              path7<-paste(path6,AS_ID_array_sel,'/', sep='')
              
              if (file.exists(path7)){
                
                
                
                
              } else {
                dir.create(file.path(path7))
                
              }
              
              setwd(path7)
              
              
              #### G1
              
              v_parameter<-"NA" 
              
              indx_carried_variants<-which(colnames(Freq_table_sel) == "carried_variants")
              indx_category<-which(colnames(Freq_table_sel) == "Fig1_Annot_Category")
              
              indx.yaxis<-which(colnames(Freq_table_sel) == "Perc")
              mycols <- vector_colors
              
              pdfname<-paste(CT_sel,"_Cumulative_frequency_by_Fig1_cat_","E_Plus_ASE_activity",".pdf",sep='')
              pdf(file=pdfname, width=5, height=4, pointsize=12)
              
              par(mai=c(0.9,0.9,0.3,0.2))
              lab <- as.character(unique(Freq_table_sel[,indx_carried_variants]))
              
              if(Condition_DEBUG == 1)
              {
                cat("lab\n")
                cat(sprintf(as.character(lab)))
                cat("\n")
              }
              
              
              plot(Freq_table_sel$TILE, Freq_table_sel[,indx.yaxis],
                   ty="n", xlab="TILES with E_Plus_ASE activity",
                   ylab="Cumulative % of variants",
                   axes=F, cex.lab=1.2, cex.lab=1.3, ylim=c(0,lim_Perc))
              
              
              
              # cat("Hello_world1\n")
              
              # points(Freq_table_sel$TILE, Freq_table_sel[,indx.yaxis], col="darkgrey", pch=19)
              
              # cat("Hello_world2\n")
              # 
              # 
              for (iteration_graph in 1:length(lab))  {
                
                if(Condition_DEBUG == 1)
                {
                  cat("Hello_world3\t")
                  cat(sprintf(as.character(lab[iteration_graph])))
                  cat("\n")
                }
                
                ind <- which(Freq_table_sel[,indx_carried_variants]==lab[iteration_graph])
                
                if(Condition_DEBUG == 1)
                {
                  cat("------------------->ind\n")
                  cat(str(ind))
                  cat("\n")
                }
                
                category_sel<-unique(Freq_table_sel$Fig1_Annot_Category[ind])
                
                position_vector_sel<-as.numeric(category_sel)
                
                color_sel<-mycols[position_vector_sel]#######################################
                
                if(Condition_DEBUG == 1)
                {
                  cat(sprintf(as.character(color_sel)))
                  cat("\n")
                }
                
                points(Freq_table_sel$TILE[ind], Freq_table_sel[,indx.yaxis][ind],
                       pch=21,lwd=5, col=color_sel)
                lines(Freq_table_sel$TILE[ind], Freq_table_sel[,indx.yaxis][ind],
                      lty=1,lwd=5, col=color_sel)
                
                if(Condition_DEBUG == 1)
                {
                  cat("END_3\n")
                }
                
                # points(Freq_table_sel$TILE[ind], Freq_table_sel[,indx.yaxis][ind], pch=19, col=color_sel)
              }#iteration_graph in 1:length(lab)
              
              if(Condition_DEBUG == 1)
              {
                cat("END_4\n")
              }
              
              
              legend("topright", legend=lab, fill=mycols, border=mycols, bty="n")
              axis(1, at=seq(0,100))
              axis(2, las=1)
              
              dev.off()
              
              #### G2
              
              v_parameter<-"NA" 
              
              indx_carried_variants<-which(colnames(Freq_table_sel) == "carried_variants")
              indx_category<-which(colnames(Freq_table_sel) == "bin_PP_DEF")
              
              indx.yaxis<-which(colnames(Freq_table_sel) == "Perc")
              mycols <- vector_colors
              
              pdfname<-paste(CT_sel,"_Cumulative_frequency_by_bin_PP_cat_","E_Plus_ASE_activity",".pdf",sep='')
              pdf(file=pdfname, width=5, height=4, pointsize=12)
              
              par(mai=c(0.9,0.9,0.3,0.2))
              lab <- as.character(unique(Freq_table_sel[,indx_carried_variants]))
              
              if(Condition_DEBUG == 1)
              {
                cat("lab\n")
                cat(sprintf(as.character(lab)))
                cat("\n")
              }
              
              
              plot(Freq_table_sel$TILE, Freq_table_sel[,indx.yaxis],
                   ty="n", xlab="TILES with E_Plus_ASE activity",
                   ylab="Cumulative % of variants",
                   axes=F, cex.lab=1.2, cex.lab=1.3, ylim=c(0,lim_Perc))
              
              
              
              # cat("Hello_world1\n")
              
              # points(Freq_table_sel$TILE, Freq_table_sel[,indx.yaxis], col="darkgrey", pch=19)
              
              # cat("Hello_world2\n")
              # 
              # 
              for (iteration_graph in 1:length(lab))  {
                
                if(Condition_DEBUG == 1)
                {
                  cat("Hello_world3\t")
                  cat(sprintf(as.character(lab[iteration_graph])))
                  cat("\n")
                }
                
                ind <- which(Freq_table_sel[,indx_carried_variants]==lab[iteration_graph])
                
                if(Condition_DEBUG == 1)
                {
                  cat("------------------->ind\n")
                  cat(str(ind))
                  cat("\n")
                }
                
                category_sel<-unique(Freq_table_sel$bin_PP_DEF[ind])
                
                position_vector_sel<-as.numeric(category_sel)
                
                color_sel<-mycols[position_vector_sel]#######################################
                
                if(Condition_DEBUG == 1)
                {
                  cat(sprintf(as.character(color_sel)))
                  cat("\n")
                }
                
                points(Freq_table_sel$TILE[ind], Freq_table_sel[,indx.yaxis][ind],
                       pch=21,lwd=5, col=color_sel)
                lines(Freq_table_sel$TILE[ind], Freq_table_sel[,indx.yaxis][ind],
                      lty=1,lwd=5, col=color_sel)
                
                if(Condition_DEBUG == 1)
                {
                  cat("END_3\n")
                }
                
                # points(Freq_table_sel$TILE[ind], Freq_table_sel[,indx.yaxis][ind], pch=19, col=color_sel)
              }#iteration_graph in 1:length(lab)
              
              if(Condition_DEBUG == 1)
              {
                cat("END_4\n")
              }
              
              
              legend("topright", legend=lab, fill=mycols, border=mycols, bty="n")
              axis(1, at=seq(0,100))
              axis(2, las=1)
              
              dev.off()
              
              
              
              list_CT[[h]]<-Freq_table_sel
              
              
              
              
              
              
            }#dim(Log_rank_df_sel)[1] >0
          }#h in 1:length(CT_levels)
          
          
          if(length(list_CT) >0)
          {
            CT_df_recovered = unique(as.data.frame(data.table::rbindlist(list_CT, fill=T), stringsAsFactors=F))
            
           
            if(Condition_DEBUG == 1)
            {
              cat("CT_df_recovered_0\n")
              cat(str(CT_df_recovered))
              cat("\n")
              #quit(status = 1)
            }
            
            CT_df_recovered$Allelic_Series_ID<-AS_ID_array_sel
            
            if(Condition_DEBUG == 1)
            {
              cat("CT_df_recovered_1\n")
              cat(str(CT_df_recovered))
              cat("\n")
              #quit(status = 1)
            }
            
            list_AS[[k]]<-CT_df_recovered
          }#length(list_CT) >0
          
          
        }#dim(check_screened_RV_included)[1] >0
      }#dim(CUMMULATIVE_CLASSES_sel)[1]
    }#k in 1:length(AS_ID_array)
    
    
    
    if(length(list_AS) >0)
    {
      AS_df_recovered = unique(as.data.frame(data.table::rbindlist(list_AS, fill=T), stringsAsFactors=F))
      
      
      if(Condition_DEBUG == 1)
      {
        cat("AS_df_recovered_0\n")
        cat(str(AS_df_recovered))
        cat("\n")
        #quit(status = 1)
      }
      
      AS_df_recovered$screened_variant<-screened_variants_sel
      
      if(Condition_DEBUG == 1)
      {
        cat("AS_df_recovered_1\n")
        cat(str(AS_df_recovered))
        cat("\n")
        #quit(status = 1)
      }
      
      list_DEF[[i]]<-AS_df_recovered
    }#length(list_AS) >0
    # #####################################################################################
    # quit(status = 1)
    
  }# i in 1:length(screened_variants)
  
  if(length(list_DEF) >0)
  {
    DEF_df = unique(as.data.frame(data.table::rbindlist(list_DEF, fill=T), stringsAsFactors=F))
    
    
    if(Condition_DEBUG == 1)
    {
      cat("DEF_df_0\n")
      cat(str(DEF_df))
      cat("\n")
      cat(str(unique(DEF_df$screened_variant)))
      cat("\n")
      cat(str(unique(DEF_df$carried_variants)))
      cat("\n")
      #quit(status = 1)
    }
    
    path5<-paste(out2,'CT_AS_printing_CUMULATIVE_CURVES','/', sep='')
    
    setwd(path5)
    
    write.table(DEF_df, file="E_Plus_ASE_AS_Printing_Cummulative_Freq_table.tsv", sep="\t", quote=F,row.names = F)
    
   
  }#length(list_DEF) >0
  
 
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
    make_option(c("--type"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--out"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--out2"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--CUMMULATIVE_CLASSES"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--VAR_Prioritization_dB"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--ALL_dB"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Supp4_Table_CURATED_PLUS_PHENOTYPES"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required.")
  )
  parser = OptionParser(usage = "137_MPRA_normalization_and_filtering_Rscript_v2.R
                        --regular_table FILE.txt
                        --replicas charac
                        --type1 type1
                        --type2 type2
                        --pvalThreshold integer
                        --FDRThreshold integer
                        --EquivalenceTable FILE.txt
                        --sharpr2Threshold charac",
                        option_list = option_list)
  opt <<- parse_args(parser)
  
  cumulative_CT_Label2_and_Log_Rank_test_enhancer(opt)
  cumulative_CT_Label2_and_Log_Rank_test_E_Plus_ASE(opt)
 
  
}
  
  
  
 

###########################################################################

system.time( main() )
