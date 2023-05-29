
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
  
  
  ### CUMULATIVE_CLASSES_PREPARED_leukemias #----
  
  CUMULATIVE_CLASSES_PREPARED_leukemias<-readRDS(file=opt$CUMULATIVE_CLASSES_PREPARED_leukemias)
  
  CUMULATIVE_CLASSES_PREPARED_leukemias$Disease<-'leukemias'
  
  cat("CUMULATIVE_CLASSES_PREPARED_leukemias_0\n")
  cat(str(CUMULATIVE_CLASSES_PREPARED_leukemias))
  cat("\n")
  cat(str(unique(CUMULATIVE_CLASSES_PREPARED_leukemias$VAR)))
  cat("\n")
  
  ### CUMULATIVE_CLASSES_PREPARED_inflammatory_disease #----
  
  CUMULATIVE_CLASSES_PREPARED_inflammatory_disease<-readRDS(file=opt$CUMULATIVE_CLASSES_PREPARED_inflammatory_disease)
  
  CUMULATIVE_CLASSES_PREPARED_inflammatory_disease$Disease<-'inflammatory_disease'
  
  cat("CUMULATIVE_CLASSES_PREPARED_inflammatory_disease_0\n")
  cat(str(CUMULATIVE_CLASSES_PREPARED_inflammatory_disease))
  cat("\n")
  cat(str(unique(CUMULATIVE_CLASSES_PREPARED_inflammatory_disease$VAR)))
  cat("\n")
  
  ### CUMULATIVE_CLASSES_PREPARED_immunodeficiencies #----
  
  CUMULATIVE_CLASSES_PREPARED_immunodeficiencies<-readRDS(file=opt$CUMULATIVE_CLASSES_PREPARED_immunodeficiencies)
  
  CUMULATIVE_CLASSES_PREPARED_immunodeficiencies$Disease<-'immunodeficiencies'
  
  
  cat("CUMULATIVE_CLASSES_PREPARED_immunodeficiencies_0\n")
  cat(str(CUMULATIVE_CLASSES_PREPARED_immunodeficiencies))
  cat("\n")
  cat(str(unique(CUMULATIVE_CLASSES_PREPARED_immunodeficiencies$VAR)))
  cat("\n")
  
  ### CUMULATIVE_CLASSES_PREPARED_blood_disease #----
  
  CUMULATIVE_CLASSES_PREPARED_blood_disease<-readRDS(file=opt$CUMULATIVE_CLASSES_PREPARED_blood_disease)
  
  CUMULATIVE_CLASSES_PREPARED_blood_disease$Disease<-'blood_disease'
  
  
  cat("CUMULATIVE_CLASSES_PREPARED_blood_disease_0\n")
  cat(str(CUMULATIVE_CLASSES_PREPARED_blood_disease))
  cat("\n")
  cat(str(unique(CUMULATIVE_CLASSES_PREPARED_blood_disease$VAR)))
  cat("\n")
  
  ### CUMULATIVE_CLASSES_PREPARED_miscellaneous #----
  
  CUMULATIVE_CLASSES_PREPARED_miscellaneous<-readRDS(file=opt$CUMULATIVE_CLASSES_PREPARED_miscellaneous)
  
  CUMULATIVE_CLASSES_PREPARED_miscellaneous$Disease<-'miscellaneous'
  
  
  cat("CUMULATIVE_CLASSES_PREPARED_miscellaneous_0\n")
  cat(str(CUMULATIVE_CLASSES_PREPARED_miscellaneous))
  cat("\n")
  cat(str(unique(CUMULATIVE_CLASSES_PREPARED_miscellaneous$VAR)))
  cat("\n")
  
  #### REP ----
  
  REP<-rbind(CUMULATIVE_CLASSES_PREPARED_leukemias,CUMULATIVE_CLASSES_PREPARED_inflammatory_disease,CUMULATIVE_CLASSES_PREPARED_immunodeficiencies,CUMULATIVE_CLASSES_PREPARED_blood_disease,CUMULATIVE_CLASSES_PREPARED_miscellaneous)
  
  cat("REP_0\n")
  cat(str(REP))
  cat("\n")
  cat(str(unique(REP$VAR)))
  cat("\n")
  
  path5<-paste(out2,'CT_CLASS_CUMULATIVE_Disease','/', sep='')
  
  if (file.exists(path5)){
    
    
    
    
  } else {
    dir.create(file.path(path5))
    
  }
  
  setwd(path5)
  write.table(REP, file="Annotation_variants_to_diseases.tsv", sep="\t", quote=F, row.names=F)
  
  #### LOOP ----
  
  disease_vector<-unique(REP$Disease)
  
  cat("disease_vector_0\n")
  cat(str(disease_vector))
  cat("\n")
  
  Condition_DEBUG <- 1
  
  list_UBER_DEF<-list()
  list_UBER_DEF_Freq<-list()
  
  for(i in 1:length(disease_vector))
  {
    disease_vector_sel<-disease_vector[i]
    
    cat("------------------------------>\t")
    cat(sprintf(as.character(disease_vector_sel)))
    cat("\n")
    
    REP_sel<-droplevels(REP[which(REP$Disease == disease_vector_sel),])
    
    cat("REP_sel_0\n")
    cat(str(REP_sel))
    cat("\n")
    cat(str(unique(REP_sel$VAR)))
    cat("\n")
    
    ########### Build df ---------------
    
    Cell_Type_levels<-levels(REP_sel$Cell_Type)
    n_Cell_Type_levels<-length(levels(REP_sel$Cell_Type))
    n_MPRA_TILES<-5
    VAR_vector<-as.character(unique(REP_sel$VAR))
    KEY_Plus_carried_variants_vector<-as.character(unique(REP_sel$KEY_Plus_carried_variants))
    
    
    
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
    
    REP_sel_subset<-unique(REP_sel[,c(which(colnames(REP_sel) == "VAR"),
                                      which(colnames(REP_sel) == "KEY_Plus_carried_variants"),
                                      which(colnames(REP_sel) == "CLASS"),
                                      which(colnames(REP_sel) == "Cell_Type"),
                                      which(colnames(REP_sel) == "enhancer_CLASS_TILES"))])
    
    if(Condition_DEBUG == 1)
    {
      cat("REP_sel_subset_0\n")
      cat(str(REP_sel_subset))
      cat("\n")
      cat(str(unique(REP_sel_subset$VAR)))
      cat("\n")
      cat(str(unique(REP_sel_subset$KEY_Plus_carried_variants)))
      cat("\n")
      cat(sprintf(as.character(names(summary(REP_sel_subset$Cell_Type)))))
      cat("\n")
      cat(sprintf(as.character(summary(REP_sel_subset$Cell_Type))))
      cat("\n")
      cat(sprintf(as.character(names(summary(REP_sel_subset$CLASS)))))
      cat("\n")
      cat(sprintf(as.character(summary(REP_sel_subset$CLASS))))
      cat("\n")
    }
    
    Log_rank_df<-merge(Log_rank_df,
                       REP_sel_subset,
                       by=c("VAR","KEY_Plus_carried_variants"),
                       all=T)
    
    Log_rank_df$CLASS<-Log_rank_df$CLASS
    
    Log_rank_df$CLASS<-factor(Log_rank_df$CLASS,
                              levels=levels(REP_sel_subset$CLASS),
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
    
    
    setwd(out2)
    
    write.table(Log_rank_df, file="check2.tsv",sep="\t",quote=F,row.names = F)
    
    # ################################################################
    # quit(status = 1)
    
    
    ############## Freq Table ----------------------------
    
    ### TOTAL ---
    
    REP_sel.dt<-data.table(REP_sel, key=c("Cell_Type","CLASS"))
    
    Freq_TOTAL<-as.data.frame(REP_sel.dt[,.(TOTAL=.N),by=key(REP_sel.dt)], stringsAsFactors=F)
    
    if(Condition_DEBUG == 1)
    {
      cat("Freq_TOTAL_0\n")
      cat(str(Freq_TOTAL))
      cat("\n")
      #quit(status = 1)
    }
    
    Freq_TOTAL.dt<-data.table(Freq_TOTAL, key=c("CLASS"))
    
    
    MAXED_df<-as.data.frame(Freq_TOTAL.dt[,.(MAX=max(TOTAL)),by=key(Freq_TOTAL.dt)], stringsAsFactors=F)
    
    
    
    if(Condition_DEBUG == 1)
    {
      cat("MAXED_df_0\n")
      cat(str(MAXED_df))
      cat("\n")
      #quit(status = 1)
    }
    
    ########## Active _tiles --------------------
    
    Log_rank_df.dt<-data.table(Log_rank_df, key=c("Cell_Type","CLASS","TILE"))
    
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
                      by=c("CLASS"),
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
    
    Freq_table$CLASS<-factor(Freq_table$CLASS,
                             levels=levels(REP_sel_subset$CLASS),
                             ordered=T)
    
    Freq_table$Cell_Type<-factor(Freq_table$Cell_Type,
                                 levels=Cell_Type_levels,
                                 ordered=T)
    
    Freq_table$TILE<-as.integer(Freq_table$TILE)
    
    Freq_table<-Freq_table[order(Freq_table$CLASS,Freq_table$Cell_Type, Freq_table$TILE),]
    
    
    setwd(out2)
    
    write.table(Freq_table, file="check3.tsv", row.names = F,quote = F,sep="\t")
    
    
    
    ############### Per CT LOOP ----
    
    
    CT_levels<-Cell_Type_levels
    
    list_ABSOLUTE_DEF<-list()
    list_ABSOLUTE_DEF_Freq<-list()
    
    Condition_DEBUG <- 1
    
    for(k in 1:length(CT_levels))
    {
      CT_sel<-CT_levels[k]
      
      cat("CT:---------------------------------------->\t")
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
      
      
      #### Log Rank test ----
      
      Condition_DEBUG <- 0
      
      list_DEF<-list()
      
      CLASS_labels<-levels(Log_rank_df_sel$CLASS)
      
      for(h in 1 :length(CLASS_labels))
      {
        BAIT<-CLASS_labels[h]
        
        cat("-----BAIT------>\t")
        cat(sprintf(as.character(BAIT)))
        cat("\t")
        
        list_PREY<-list()
        
        for(j in 1 :length(CLASS_labels))
        {
          
          PREY<-CLASS_labels[j]
          
          cat("------PREY----->\t")
          cat(sprintf(as.character(PREY)))
          cat("\n")
          
          Mini_comparison<-unique(c(BAIT,PREY))
          
          if(Condition_DEBUG == 1)
          {
            cat("Mini_comparison_0\n")
            cat(str(Mini_comparison))
            cat("\n")
          }
          
          Log_rank_df_sel_subset<-Log_rank_df_sel[which(Log_rank_df_sel$CLASS%in%Mini_comparison),]
          
          if(Condition_DEBUG == 1)
          {
            cat("Log_rank_df_sel_subset_0\n")
            cat(str(Log_rank_df_sel_subset))
            cat("\n")
          }
          
          FLAG_diversity<-as.integer(length(unique(Log_rank_df_sel_subset$CLASS)))
          
          if(Condition_DEBUG == 1)
          {
            cat("FLAG_diversity_0\n")
            cat(sprintf(as.character(FLAG_diversity)))
            cat("\n")
            
            cat(str(FLAG_diversity))
            cat("\n")
          }
          
          if(FLAG_diversity >1)
          {
            cat("Hello_world\n")
            
            Mini_comparison_c1<-as.character(Mini_comparison[1])
            Mini_comparison_c2<-as.character(Mini_comparison[2])
            
            FLAG_diversity_2<-length(unique(Log_rank_df_sel_subset$ACTIVE))
            
            if(Condition_DEBUG == 1)
            {
              cat("FLAG_diversity_2_0\n")
              cat(str(FLAG_diversity_2))
              cat("\n")
            }
            
            if(FLAG_diversity_2 >1)
            {
              if(Condition_DEBUG == 1)
              {
                cat("----------->\t")
                cat(sprintf(as.character(Mini_comparison_c1)))
                cat("\t")
                cat(sprintf(as.character(Mini_comparison_c2)))
                cat("\n")
                
                cat("Log_rank_df_sel_subset_0\n")
                cat(str(Log_rank_df_sel_subset))
                cat("\n")
                cat(str(unique(Log_rank_df_sel_subset$VAR)))
                cat("\n")
                cat(sprintf(as.character(names(summary(as.factor(Log_rank_df_sel_subset$ACTIVE))))))
                cat("\n")
                cat(sprintf(as.character(summary(as.factor(Log_rank_df_sel_subset$ACTIVE)))))
                cat("\n")
              }
              
              #, "less", "greater"
              
              STATS_test<-logrank.test(
                Log_rank_df_sel_subset$TILE,
                Log_rank_df_sel_subset$ACTIVE,
                Log_rank_df_sel_subset$CLASS,
                alternative = c("two.sided"),
                rho = 0,
                gamma = 0,
                event_time_weights = NULL
              )
              
              # if(Condition_DEBUG == 1)
              # {
              #   cat("STATS_test_0\n")
              #   cat(str(STATS_test))
              #   cat("\n")
              # }
              
              
              p_value<-STATS_test$test$p
              minus_log_val<-round(-1*log10(p_value),2)
              
              if(Condition_DEBUG == 1)
              {
                cat("minus_log_val_0\n")
                cat(str(minus_log_val))
                cat("\n")
              }
              
              
              a.df<-as.data.frame(cbind(Mini_comparison_c1,Mini_comparison_c2,minus_log_val), stringsAsFactors=F)
              
              
              colnames(a.df)<-c("CLASS_c1","CLASS_c2","minus_logpval")
              a.df$minus_logpval<-as.numeric(a.df$minus_logpval)
              
              if(Condition_DEBUG == 1)
              {
                cat("a.df\n")
                cat(str(a.df))
                cat("\n")
              }
              
              list_PREY[[j]]<-a.df
              
              # ###############################################
              # quit(status = 1)
              
            }#FLAG_diversity_2 >1
          }#FLAG_diversity >1
        }# j in 1 :length(unique(Log_rank_df_sel$CLASS)
        
        Condition_DEBUG <- 0
        
        if(length(list_PREY) >0)
        {
          
          PREY_df = unique(as.data.frame(data.table::rbindlist(list_PREY, fill=T), stringsAsFactors=F))
          
          if(Condition_DEBUG == 1)
          {
            cat("PREY_df_0\n")
            cat(str(PREY_df))
            cat("\n")
            #quit(status = 1)
          }
          
          list_DEF[[h]]<-PREY_df
          
        }# length(list_PREY) >0
        
        Condition_DEBUG <- 0
        
        
      }#h in 1 :length(unique(Log_rank_df_sel$CLASS))
      
      
      Condition_DEBUG <- 0
      
      if(length(list_DEF) >0)
      {
        FINAL_df_per_CT = unique(as.data.frame(data.table::rbindlist(list_DEF, fill=T), stringsAsFactors=F))
        
        FINAL_df_per_CT$Cell_Type<-CT_sel
        
        if(Condition_DEBUG == 1)
        {
          cat("FINAL_df_per_CT_0\n")
          cat(str(FINAL_df_per_CT))
          cat("\n")
          #quit(status = 1)
        }
        
        list_ABSOLUTE_DEF[[k]]<-FINAL_df_per_CT
        
      }#length(list_DEF) >0
      
      Condition_DEBUG <- 0
      
      #### GRAPH -----
      
      Condition_DEBUG <- 1
      
      vector_colors_CLASS<-c("greenyellow","blue","dark cyan","black","red")
      
      Freq_table_sel<-Freq_table[which(Freq_table$Cell_Type == CT_sel),]
      
      lim_Perc<-max(Freq_table_sel$Perc) + 10
      
      if(Condition_DEBUG == 1)
      {
        cat("Freq_table_sel_0\n")
        cat(str(Freq_table_sel))
        cat("\n")
        cat("lim_Perc\n")
        cat(str(lim_Perc))
        cat("\n")
        #quit(status = 1)
      }
      
      
      
      path6<-paste(path5,disease_vector_sel,'/', sep='')
      
      if (file.exists(path6)){
        
        
        
        
      } else {
        dir.create(file.path(path6))
        
      }
      
      setwd(path6)
      
      
      v_parameter<-"NA" 
      
      indx_CLASS<-which(colnames(Freq_table_sel) == "CLASS")
      indx.yaxis<-which(colnames(Freq_table_sel) == "Perc")
      mycols <- vector_colors_CLASS
      
      pdfname<-paste(CT_sel,"_Cumulative_frequency_","enhancer_activity",".pdf",sep='')
      pdf(file=pdfname, width=5, height=4, pointsize=12)
      
      par(mai=c(0.9,0.9,0.3,0.2))
      lab <- as.character(unique(Freq_table_sel[,indx_CLASS]))
      
      cat("lab\n")
      cat(sprintf(as.character(lab)))
      cat("\n")
      
      
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
        
        cat("Hello_world3\t")
        cat(sprintf(as.character(lab[iteration_graph])))
        cat("\n")
        
        ind <- which(Freq_table_sel[,indx_CLASS]==lab[iteration_graph])
        
        cat("------------------->ind\n")
        cat(str(ind))
        cat("\n")
        color_sel<-mycols[iteration_graph]
        
        cat(sprintf(as.character(color_sel)))
        cat("\n")
        
        points(Freq_table_sel$TILE[ind], Freq_table_sel[,indx.yaxis][ind],
               pch=21,lwd=5, col=color_sel)
        lines(Freq_table_sel$TILE[ind], Freq_table_sel[,indx.yaxis][ind],
              lty=1,lwd=5, col=color_sel)
        
        cat("END_3\n")
        
        # points(Freq_table_sel$TILE[ind], Freq_table_sel[,indx.yaxis][ind], pch=19, col=color_sel)
      }#iteration_graph in 1:length(lab)
      
      cat("END_4\n")
      
      
      legend("topright", legend=lab, fill=mycols, border=mycols, bty="n")
      axis(1, at=seq(0,100))
      axis(2, las=1)
      
      dev.off()
      
      
      
      list_ABSOLUTE_DEF_Freq[[k]]<-Freq_table_sel
      
      
      Condition_DEBUG <- 0
      
      # ############################################################### HERE HERE
      # quit(status = 1)
      
    }################################ k 1:length(CT_levels)
    
    Condition_DEBUG <- 1
    
    if(length(list_ABSOLUTE_DEF) >0)
    {
      FINAL_df_ABSOLUTE = unique(as.data.frame(data.table::rbindlist(list_ABSOLUTE_DEF, fill=T), stringsAsFactors=F))
      
      FINAL_df_ABSOLUTE$Disease<-disease_vector_sel
      
      
      if(Condition_DEBUG == 1)
      {
        cat("FINAL_df_ABSOLUTE_0\n")
        cat(str(FINAL_df_ABSOLUTE))
        cat("\n")
        #quit(status = 1)
      }
      
      list_UBER_DEF[[i]]<-FINAL_df_ABSOLUTE
      
    }#length(list_ABSOLUTE_DEF) >0
    
    if(length(list_ABSOLUTE_DEF_Freq) >0)
    {
      FINAL_ABSOLUTE_DEF_Freq = unique(as.data.frame(data.table::rbindlist(list_ABSOLUTE_DEF_Freq, fill=T), stringsAsFactors=F))
      
      
      if(Condition_DEBUG == 1)
      {
        cat("FINAL_ABSOLUTE_DEF_Freq_0\n")
        cat(str(FINAL_ABSOLUTE_DEF_Freq))
        cat("\n")
        #quit(status = 1)
      }
      
      list_UBER_DEF_Freq[[i]]<-FINAL_ABSOLUTE_DEF_Freq
      
      
    }#length(list_ABSOLUTE_DEF_Freq) >0
    
    
    
  }# i in 1:length(disease_vector)
  
  
  if(length(list_UBER_DEF) >0)
  {
    UBER_df_ABSOLUTE = unique(as.data.frame(data.table::rbindlist(list_UBER_DEF, fill=T), stringsAsFactors=F))
    
    
    if(Condition_DEBUG == 1)
    {
      cat("UBER_df_ABSOLUTE_0\n")
      cat(str(UBER_df_ABSOLUTE))
      cat("\n")
      #quit(status = 1)
    }
    
    setwd(path5)
    
    
    
    write.table(UBER_df_ABSOLUTE, file="enhancer_CLASS_log_Rank_test_STAT.tsv", sep="\t", row.names = F, quote = F)
    
    
  }#length(list_UBER_DEF) >0
  
  if(length(list_ABSOLUTE_DEF_Freq) >0)
  {
    UBER_ABSOLUTE_DEF_Freq = unique(as.data.frame(data.table::rbindlist(list_UBER_DEF_Freq, fill=T), stringsAsFactors=F))
    
    
    if(Condition_DEBUG == 1)
    {
      cat("UBER_ABSOLUTE_DEF_Freq_0\n")
      cat(str(UBER_ABSOLUTE_DEF_Freq))
      cat("\n")
      #quit(status = 1)
    }
    
    setwd(path5)
    
    write.table(UBER_ABSOLUTE_DEF_Freq, file="enhancer_CLASS_Cumulative_Freq_table.tsv", sep="\t", quote=F,row.names = F)
    
  }#length(list_ABSOLUTE_DEF_Freq) >0
  
  
  
  
  
  
  
  
  
  
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
  
  
  ### CUMULATIVE_CLASSES_PREPARED_leukemias #----
  
  CUMULATIVE_CLASSES_PREPARED_leukemias<-readRDS(file=opt$CUMULATIVE_CLASSES_PREPARED_leukemias)
  
  CUMULATIVE_CLASSES_PREPARED_leukemias$Disease<-'leukemias'
  
  cat("CUMULATIVE_CLASSES_PREPARED_leukemias_0\n")
  cat(str(CUMULATIVE_CLASSES_PREPARED_leukemias))
  cat("\n")
  cat(str(unique(CUMULATIVE_CLASSES_PREPARED_leukemias$VAR)))
  cat("\n")
  
  ### CUMULATIVE_CLASSES_PREPARED_inflammatory_disease #----
  
  CUMULATIVE_CLASSES_PREPARED_inflammatory_disease<-readRDS(file=opt$CUMULATIVE_CLASSES_PREPARED_inflammatory_disease)
  
  CUMULATIVE_CLASSES_PREPARED_inflammatory_disease$Disease<-'inflammatory_disease'
  
  cat("CUMULATIVE_CLASSES_PREPARED_inflammatory_disease_0\n")
  cat(str(CUMULATIVE_CLASSES_PREPARED_inflammatory_disease))
  cat("\n")
  cat(str(unique(CUMULATIVE_CLASSES_PREPARED_inflammatory_disease$VAR)))
  cat("\n")
  
  ### CUMULATIVE_CLASSES_PREPARED_immunodeficiencies #----
  
  CUMULATIVE_CLASSES_PREPARED_immunodeficiencies<-readRDS(file=opt$CUMULATIVE_CLASSES_PREPARED_immunodeficiencies)
  
  CUMULATIVE_CLASSES_PREPARED_immunodeficiencies$Disease<-'immunodeficiencies'
  
  
  cat("CUMULATIVE_CLASSES_PREPARED_immunodeficiencies_0\n")
  cat(str(CUMULATIVE_CLASSES_PREPARED_immunodeficiencies))
  cat("\n")
  cat(str(unique(CUMULATIVE_CLASSES_PREPARED_immunodeficiencies$VAR)))
  cat("\n")
  
  ### CUMULATIVE_CLASSES_PREPARED_blood_disease #----
  
  CUMULATIVE_CLASSES_PREPARED_blood_disease<-readRDS(file=opt$CUMULATIVE_CLASSES_PREPARED_blood_disease)
  
  CUMULATIVE_CLASSES_PREPARED_blood_disease$Disease<-'blood_disease'
  
  
  cat("CUMULATIVE_CLASSES_PREPARED_blood_disease_0\n")
  cat(str(CUMULATIVE_CLASSES_PREPARED_blood_disease))
  cat("\n")
  cat(str(unique(CUMULATIVE_CLASSES_PREPARED_blood_disease$VAR)))
  cat("\n")
  
  ### CUMULATIVE_CLASSES_PREPARED_miscellaneous #----
  
  CUMULATIVE_CLASSES_PREPARED_miscellaneous<-readRDS(file=opt$CUMULATIVE_CLASSES_PREPARED_miscellaneous)
  
  CUMULATIVE_CLASSES_PREPARED_miscellaneous$Disease<-'miscellaneous'
  
  
  cat("CUMULATIVE_CLASSES_PREPARED_miscellaneous_0\n")
  cat(str(CUMULATIVE_CLASSES_PREPARED_miscellaneous))
  cat("\n")
  cat(str(unique(CUMULATIVE_CLASSES_PREPARED_miscellaneous$VAR)))
  cat("\n")
  
  #### REP ----
  
  REP<-rbind(CUMULATIVE_CLASSES_PREPARED_leukemias,CUMULATIVE_CLASSES_PREPARED_inflammatory_disease,CUMULATIVE_CLASSES_PREPARED_immunodeficiencies,CUMULATIVE_CLASSES_PREPARED_blood_disease,CUMULATIVE_CLASSES_PREPARED_miscellaneous)
  
  cat("REP_0\n")
  cat(str(REP))
  cat("\n")
  cat(str(unique(REP$VAR)))
  cat("\n")
  
  path5<-paste(out2,'CT_CLASS_CUMULATIVE_Disease','/', sep='')
  
  if (file.exists(path5)){
    
    
    
    
  } else {
    dir.create(file.path(path5))
    
  }
  
  setwd(path5)
  write.table(REP, file="Annotation_variants_to_diseases.tsv", sep="\t", quote=F, row.names=F)
  
  #### LOOP ----
  
  disease_vector<-unique(REP$Disease)
  
  cat("disease_vector_0\n")
  cat(str(disease_vector))
  cat("\n")
  
  Condition_DEBUG <- 1
  
  list_UBER_DEF<-list()
  list_UBER_DEF_Freq<-list()
  
  for(i in 1:length(disease_vector))
  {
    disease_vector_sel<-disease_vector[i]
    
    cat("------------------------------>\t")
    cat(sprintf(as.character(disease_vector_sel)))
    cat("\n")
    
    REP_sel<-droplevels(REP[which(REP$Disease == disease_vector_sel),])
    
    cat("REP_sel_0\n")
    cat(str(REP_sel))
    cat("\n")
    cat(str(unique(REP_sel$VAR)))
    cat("\n")
    
    ########### Build df ---------------
    
    Cell_Type_levels<-levels(REP_sel$Cell_Type)
    n_Cell_Type_levels<-length(levels(REP_sel$Cell_Type))
    n_MPRA_TILES<-5
    VAR_vector<-as.character(unique(REP_sel$VAR))
    KEY_Plus_carried_variants_vector<-as.character(unique(REP_sel$KEY_Plus_carried_variants))
    
    
    
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
    
    REP_sel_subset<-unique(REP_sel[,c(which(colnames(REP_sel) == "VAR"),
                                                                                    which(colnames(REP_sel) == "KEY_Plus_carried_variants"),
                                                                                    which(colnames(REP_sel) == "CLASS"),
                                                                                    which(colnames(REP_sel) == "Cell_Type"),
                                                                                    which(colnames(REP_sel) == "E_Plus_ASE_CLASS_TILES"))])
    
    if(Condition_DEBUG == 1)
    {
      cat("REP_sel_subset_0\n")
      cat(str(REP_sel_subset))
      cat("\n")
      cat(str(unique(REP_sel_subset$VAR)))
      cat("\n")
      cat(str(unique(REP_sel_subset$KEY_Plus_carried_variants)))
      cat("\n")
      cat(sprintf(as.character(names(summary(REP_sel_subset$Cell_Type)))))
      cat("\n")
      cat(sprintf(as.character(summary(REP_sel_subset$Cell_Type))))
      cat("\n")
      cat(sprintf(as.character(names(summary(REP_sel_subset$CLASS)))))
      cat("\n")
      cat(sprintf(as.character(summary(REP_sel_subset$CLASS))))
      cat("\n")
    }
    
    Log_rank_df<-merge(Log_rank_df,
                       REP_sel_subset,
                       by=c("VAR","KEY_Plus_carried_variants"),
                       all=T)
    
    Log_rank_df$CLASS<-Log_rank_df$CLASS
    
    Log_rank_df$CLASS<-factor(Log_rank_df$CLASS,
                                   levels=levels(REP_sel_subset$CLASS),
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
    
    
    setwd(out2)
    
    write.table(Log_rank_df, file="check2.tsv",sep="\t",quote=F,row.names = F)
    
    # ################################################################
    # quit(status = 1)
    
    
    ############## Freq Table ----------------------------
    
    ### TOTAL ---
    
    REP_sel.dt<-data.table(REP_sel, key=c("Cell_Type","CLASS"))
    
    Freq_TOTAL<-as.data.frame(REP_sel.dt[,.(TOTAL=.N),by=key(REP_sel.dt)], stringsAsFactors=F)
    
    if(Condition_DEBUG == 1)
    {
      cat("Freq_TOTAL_0\n")
      cat(str(Freq_TOTAL))
      cat("\n")
      #quit(status = 1)
    }
    
    Freq_TOTAL.dt<-data.table(Freq_TOTAL, key=c("CLASS"))
    
    
    MAXED_df<-as.data.frame(Freq_TOTAL.dt[,.(MAX=max(TOTAL)),by=key(Freq_TOTAL.dt)], stringsAsFactors=F)
    
    
    
    if(Condition_DEBUG == 1)
    {
      cat("MAXED_df_0\n")
      cat(str(MAXED_df))
      cat("\n")
      #quit(status = 1)
    }
    
    ########## Active _tiles --------------------
    
    Log_rank_df.dt<-data.table(Log_rank_df, key=c("Cell_Type","CLASS","TILE"))
    
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
                      by=c("CLASS"),
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
    
    Freq_table$CLASS<-factor(Freq_table$CLASS,
                                  levels=levels(REP_sel_subset$CLASS),
                                  ordered=T)
    
    Freq_table$Cell_Type<-factor(Freq_table$Cell_Type,
                                 levels=Cell_Type_levels,
                                 ordered=T)
    
    Freq_table$TILE<-as.integer(Freq_table$TILE)
    
    Freq_table<-Freq_table[order(Freq_table$CLASS,Freq_table$Cell_Type, Freq_table$TILE),]
    
    
    setwd(out2)
    
    write.table(Freq_table, file="check3.tsv", row.names = F,quote = F,sep="\t")
    
    
    
    ############### Per CT LOOP ----
    
    
    CT_levels<-Cell_Type_levels
    
    list_ABSOLUTE_DEF<-list()
    list_ABSOLUTE_DEF_Freq<-list()
    
    Condition_DEBUG <- 1
    
    for(k in 1:length(CT_levels))
    {
      CT_sel<-CT_levels[k]
      
      cat("CT:---------------------------------------->\t")
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
      
      
      #### Log Rank test ----
      
      Condition_DEBUG <- 0
      
      list_DEF<-list()
      
      CLASS_labels<-levels(Log_rank_df_sel$CLASS)
      
      for(h in 1 :length(CLASS_labels))
      {
        BAIT<-CLASS_labels[h]
        
        cat("-----BAIT------>\t")
        cat(sprintf(as.character(BAIT)))
        cat("\t")
        
        list_PREY<-list()
        
        for(j in 1 :length(CLASS_labels))
        {
          
          PREY<-CLASS_labels[j]
          
          cat("------PREY----->\t")
          cat(sprintf(as.character(PREY)))
          cat("\n")
          
          Mini_comparison<-unique(c(BAIT,PREY))
          
          if(Condition_DEBUG == 1)
          {
            cat("Mini_comparison_0\n")
            cat(str(Mini_comparison))
            cat("\n")
          }
          
          Log_rank_df_sel_subset<-Log_rank_df_sel[which(Log_rank_df_sel$CLASS%in%Mini_comparison),]
          
          if(Condition_DEBUG == 1)
          {
            cat("Log_rank_df_sel_subset_0\n")
            cat(str(Log_rank_df_sel_subset))
            cat("\n")
          }
          
          FLAG_diversity<-as.integer(length(unique(Log_rank_df_sel_subset$CLASS)))
          
          if(Condition_DEBUG == 1)
          {
            cat("FLAG_diversity_0\n")
            cat(sprintf(as.character(FLAG_diversity)))
            cat("\n")
            
            cat(str(FLAG_diversity))
            cat("\n")
          }
          
          if(FLAG_diversity >1)
          {
            cat("Hello_world\n")
            
            Mini_comparison_c1<-as.character(Mini_comparison[1])
            Mini_comparison_c2<-as.character(Mini_comparison[2])
            
            FLAG_diversity_2<-length(unique(Log_rank_df_sel_subset$ACTIVE))
            
            if(Condition_DEBUG == 1)
            {
              cat("FLAG_diversity_2_0\n")
              cat(str(FLAG_diversity_2))
              cat("\n")
            }
            
            if(FLAG_diversity_2 >1)
            {
              if(Condition_DEBUG == 1)
              {
                cat("----------->\t")
                cat(sprintf(as.character(Mini_comparison_c1)))
                cat("\t")
                cat(sprintf(as.character(Mini_comparison_c2)))
                cat("\n")
                
                cat("Log_rank_df_sel_subset_0\n")
                cat(str(Log_rank_df_sel_subset))
                cat("\n")
                cat(str(unique(Log_rank_df_sel_subset$VAR)))
                cat("\n")
                cat(sprintf(as.character(names(summary(as.factor(Log_rank_df_sel_subset$ACTIVE))))))
                cat("\n")
                cat(sprintf(as.character(summary(as.factor(Log_rank_df_sel_subset$ACTIVE)))))
                cat("\n")
              }
              
              #, "less", "greater"
              
              STATS_test<-logrank.test(
                Log_rank_df_sel_subset$TILE,
                Log_rank_df_sel_subset$ACTIVE,
                Log_rank_df_sel_subset$CLASS,
                alternative = c("two.sided"),
                rho = 0,
                gamma = 0,
                event_time_weights = NULL
              )
              
              # if(Condition_DEBUG == 1)
              # {
              #   cat("STATS_test_0\n")
              #   cat(str(STATS_test))
              #   cat("\n")
              # }
              
              
              p_value<-STATS_test$test$p
              minus_log_val<-round(-1*log10(p_value),2)
              
              if(Condition_DEBUG == 1)
              {
                cat("minus_log_val_0\n")
                cat(str(minus_log_val))
                cat("\n")
              }
              
              
              a.df<-as.data.frame(cbind(Mini_comparison_c1,Mini_comparison_c2,minus_log_val), stringsAsFactors=F)
              
              
              colnames(a.df)<-c("CLASS_c1","CLASS_c2","minus_logpval")
              a.df$minus_logpval<-as.numeric(a.df$minus_logpval)
              
              if(Condition_DEBUG == 1)
              {
                cat("a.df\n")
                cat(str(a.df))
                cat("\n")
              }
              
              list_PREY[[j]]<-a.df
              
              # ###############################################
              # quit(status = 1)
              
            }#FLAG_diversity_2 >1
          }#FLAG_diversity >1
        }# j in 1 :length(unique(Log_rank_df_sel$CLASS)
        
        Condition_DEBUG <- 0
        
        if(length(list_PREY) >0)
        {
          
          PREY_df = unique(as.data.frame(data.table::rbindlist(list_PREY, fill=T), stringsAsFactors=F))
          
          if(Condition_DEBUG == 1)
          {
            cat("PREY_df_0\n")
            cat(str(PREY_df))
            cat("\n")
            #quit(status = 1)
          }
          
          list_DEF[[h]]<-PREY_df
          
        }# length(list_PREY) >0
        
        Condition_DEBUG <- 0
        
        
      }#h in 1 :length(unique(Log_rank_df_sel$CLASS))
      
      
      Condition_DEBUG <- 0
      
      if(length(list_DEF) >0)
      {
        FINAL_df_per_CT = unique(as.data.frame(data.table::rbindlist(list_DEF, fill=T), stringsAsFactors=F))
        
        FINAL_df_per_CT$Cell_Type<-CT_sel
        
        if(Condition_DEBUG == 1)
        {
          cat("FINAL_df_per_CT_0\n")
          cat(str(FINAL_df_per_CT))
          cat("\n")
          #quit(status = 1)
        }
        
        list_ABSOLUTE_DEF[[k]]<-FINAL_df_per_CT
        
      }#length(list_DEF) >0
      
      Condition_DEBUG <- 0
      
      #### GRAPH -----
      
      Condition_DEBUG <- 1
      
      vector_colors_CLASS<-c("greenyellow","blue","dark cyan","black","red")
      
      Freq_table_sel<-Freq_table[which(Freq_table$Cell_Type == CT_sel),]
      
      lim_Perc<-max(Freq_table_sel$Perc) + 10
      
      if(Condition_DEBUG == 1)
      {
        cat("Freq_table_sel_0\n")
        cat(str(Freq_table_sel))
        cat("\n")
        cat("lim_Perc\n")
        cat(str(lim_Perc))
        cat("\n")
        #quit(status = 1)
      }
      
      
     
      path6<-paste(path5,disease_vector_sel,'/', sep='')
      
      if (file.exists(path6)){
        
        
        
        
      } else {
        dir.create(file.path(path6))
        
      }
      
      setwd(path6)
      
      
      v_parameter<-"NA" 
      
      indx_CLASS<-which(colnames(Freq_table_sel) == "CLASS")
      indx.yaxis<-which(colnames(Freq_table_sel) == "Perc")
      mycols <- vector_colors_CLASS
      
      pdfname<-paste(CT_sel,"_Cumulative_frequency_","E_Plus_ASE_activity",".pdf",sep='')
      pdf(file=pdfname, width=5, height=4, pointsize=12)
      
      par(mai=c(0.9,0.9,0.3,0.2))
      lab <- as.character(unique(Freq_table_sel[,indx_CLASS]))
      
      cat("lab\n")
      cat(sprintf(as.character(lab)))
      cat("\n")
      
      
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
        
        cat("Hello_world3\t")
        cat(sprintf(as.character(lab[iteration_graph])))
        cat("\n")
        
        ind <- which(Freq_table_sel[,indx_CLASS]==lab[iteration_graph])
        
        cat("------------------->ind\n")
        cat(str(ind))
        cat("\n")
        color_sel<-mycols[iteration_graph]
        
        cat(sprintf(as.character(color_sel)))
        cat("\n")
        
        points(Freq_table_sel$TILE[ind], Freq_table_sel[,indx.yaxis][ind],
               pch=21,lwd=5, col=color_sel)
        lines(Freq_table_sel$TILE[ind], Freq_table_sel[,indx.yaxis][ind],
              lty=1,lwd=5, col=color_sel)
        
        cat("END_3\n")
        
        # points(Freq_table_sel$TILE[ind], Freq_table_sel[,indx.yaxis][ind], pch=19, col=color_sel)
      }#iteration_graph in 1:length(lab)
      
      cat("END_4\n")
      
      
      legend("topright", legend=lab, fill=mycols, border=mycols, bty="n")
      axis(1, at=seq(0,100))
      axis(2, las=1)
      
      dev.off()
      
      
      
      list_ABSOLUTE_DEF_Freq[[k]]<-Freq_table_sel
      
      
      Condition_DEBUG <- 0
      
      # ############################################################### HERE HERE
      # quit(status = 1)
      
    }################################ k 1:length(CT_levels)
    
    Condition_DEBUG <- 1
    
    if(length(list_ABSOLUTE_DEF) >0)
    {
      FINAL_df_ABSOLUTE = unique(as.data.frame(data.table::rbindlist(list_ABSOLUTE_DEF, fill=T), stringsAsFactors=F))
      
      FINAL_df_ABSOLUTE$Disease<-disease_vector_sel
      
      
      if(Condition_DEBUG == 1)
      {
        cat("FINAL_df_ABSOLUTE_0\n")
        cat(str(FINAL_df_ABSOLUTE))
        cat("\n")
        #quit(status = 1)
      }
      
      list_UBER_DEF[[i]]<-FINAL_df_ABSOLUTE
      
    }#length(list_ABSOLUTE_DEF) >0
    
    if(length(list_ABSOLUTE_DEF_Freq) >0)
    {
      FINAL_ABSOLUTE_DEF_Freq = unique(as.data.frame(data.table::rbindlist(list_ABSOLUTE_DEF_Freq, fill=T), stringsAsFactors=F))
      
      
      if(Condition_DEBUG == 1)
      {
        cat("FINAL_ABSOLUTE_DEF_Freq_0\n")
        cat(str(FINAL_ABSOLUTE_DEF_Freq))
        cat("\n")
        #quit(status = 1)
      }
      
      list_UBER_DEF_Freq[[i]]<-FINAL_ABSOLUTE_DEF_Freq
      
      
    }#length(list_ABSOLUTE_DEF_Freq) >0
    
   
    
  }# i in 1:length(disease_vector)
  
  
  if(length(list_UBER_DEF) >0)
  {
    UBER_df_ABSOLUTE = unique(as.data.frame(data.table::rbindlist(list_UBER_DEF, fill=T), stringsAsFactors=F))
    
    
    if(Condition_DEBUG == 1)
    {
      cat("UBER_df_ABSOLUTE_0\n")
      cat(str(UBER_df_ABSOLUTE))
      cat("\n")
      #quit(status = 1)
    }
    
    setwd(path5)
    
    
    
    write.table(UBER_df_ABSOLUTE, file="E_Plus_ASE_CLASS_log_Rank_test_STAT.tsv", sep="\t", row.names = F, quote = F)
    

  }#length(list_UBER_DEF) >0
  
  if(length(list_ABSOLUTE_DEF_Freq) >0)
  {
    UBER_ABSOLUTE_DEF_Freq = unique(as.data.frame(data.table::rbindlist(list_UBER_DEF_Freq, fill=T), stringsAsFactors=F))
    
    
    if(Condition_DEBUG == 1)
    {
      cat("UBER_ABSOLUTE_DEF_Freq_0\n")
      cat(str(UBER_ABSOLUTE_DEF_Freq))
      cat("\n")
      #quit(status = 1)
    }
    
    setwd(path5)
    
    write.table(UBER_ABSOLUTE_DEF_Freq, file="E_Plus_ASE_CLASS_Cumulative_Freq_table.tsv", sep="\t", quote=F,row.names = F)
    
  }#length(list_ABSOLUTE_DEF_Freq) >0
  
  
 
  
  
  
  
  
  
 
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
    make_option(c("--CUMULATIVE_CLASSES_PREPARED_leukemias"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--CUMULATIVE_CLASSES_PREPARED_inflammatory_disease"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--CUMULATIVE_CLASSES_PREPARED_immunodeficiencies"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--CUMULATIVE_CLASSES_PREPARED_blood_disease"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--CUMULATIVE_CLASSES_PREPARED_miscellaneous"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--ALL_dB"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--VAR_Prioritization_dB"), type="character", default=NULL, 
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
