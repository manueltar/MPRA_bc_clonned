
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



opt = NULL


opt = NULL


UpSetR_like_E_Plus_ASE = function(option_list)
{
  
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
  
  #### Categories colors ----
  
  
  df_CSQ_colors<-readRDS(file=opt$CSQ_colors)
  
  df_CSQ_colors$color[df_CSQ_colors$VEP_DEF_LABELS == "TFBS"]<-"red"
  
  
  cat("df_CSQ_colors_0\n")
  cat(str(df_CSQ_colors))
  cat("\n")
  
  #### READ KEY_collapse_P----
  
  KEY_collapse_P<-readRDS(file=opt$KEY_collapse_P)
  
  KEY_collapse_P<-KEY_collapse_P[order(KEY_collapse_P$KEY_Plus_carried_variants,KEY_collapse_P$Cell_Type),]
  
  cat("KEY_collapse_P_0\n")
  cat(str(KEY_collapse_P))
  cat("\n")
  cat(str(unique(KEY_collapse_P$carried_variants)))
  cat("\n")
  
  #### READ KEY_collapse_NP_1_3----
  
  KEY_collapse_NP_1_3<-readRDS(file=opt$KEY_collapse_NP_1_3)
  
  KEY_collapse_NP_1_3<-KEY_collapse_NP_1_3[order(KEY_collapse_NP_1_3$KEY_Plus_carried_variants,KEY_collapse_NP_1_3$Cell_Type),]
  
  cat("KEY_collapse_NP_1_3_0\n")
  cat(str(KEY_collapse_NP_1_3))
  cat("\n")
  cat(str(unique(KEY_collapse_NP_1_3$carried_variants)))
  cat("\n")
  
  #### READ KEY_collapse_NP_3----
  
  KEY_collapse_NP_3<-readRDS(file=opt$KEY_collapse_NP_3)
  
  KEY_collapse_NP_3<-KEY_collapse_NP_3[order(KEY_collapse_NP_3$KEY_Plus_carried_variants,KEY_collapse_NP_3$Cell_Type),]
  
  cat("KEY_collapse_NP_3_0\n")
  cat(str(KEY_collapse_NP_3))
  cat("\n")
  cat(str(unique(KEY_collapse_NP_3$carried_variants)))
  cat("\n")
  
  #### READ df_Cell_colors----
  
  df_Cell_colors<-readRDS(file=opt$df_Cell_colors)
  

  cat("df_Cell_colors_0\n")
  cat(str(df_Cell_colors))
  cat("\n")
  
  
  #### Restrict to variants explored in P----
  
  
  KEY_collapse_NP_1_3_subset<-unique(KEY_collapse_NP_1_3[which(KEY_collapse_NP_1_3$carried_variants%in%KEY_collapse_P$carried_variants),-which(colnames(KEY_collapse_NP_1_3) == "VEP_DEF_LABELS")])
  
  
  KEY_collapse_NP_1_3_subset$EXP_GROUP<-"MPRA_NP_1_3"
  
  cat("KEY_collapse_NP_1_3_subset_0\n")
  cat(str(KEY_collapse_NP_1_3_subset))
  cat("\n")
  cat(str(unique(KEY_collapse_NP_1_3_subset$carried_variants)))
  cat("\n")
  
  KEY_collapse_NP_3_subset<-unique(KEY_collapse_NP_3[which(KEY_collapse_NP_3$carried_variants%in%KEY_collapse_P$carried_variants),-which(colnames(KEY_collapse_NP_3) == "VEP_DEF_LABELS")])
  
  KEY_collapse_NP_3_subset$EXP_GROUP<-"MPRA_NP_3"
  
  
  cat("KEY_collapse_NP_3_subset_0\n")
  cat(str(KEY_collapse_NP_3_subset))
  cat("\n")
  cat(str(unique(KEY_collapse_NP_3_subset$carried_variants)))
  cat("\n")
  
  KEY_collapse_P_subset<-unique(KEY_collapse_P[which(KEY_collapse_P$carried_variants%in%KEY_collapse_NP_3_subset$carried_variants),-which(colnames(KEY_collapse_P) == "VEP_DEF_LABELS")])
  
  KEY_collapse_P_subset$EXP_GROUP<-"MPRA_P"
  
  
  cat("KEY_collapse_P_subset_0\n")
  cat(str(KEY_collapse_P_subset))
  cat("\n")
  cat(str(unique(KEY_collapse_P_subset$carried_variants)))
  cat("\n")
  
  
  #### Merge ----
  
  ALL<-rbind(KEY_collapse_P_subset,KEY_collapse_NP_1_3_subset,KEY_collapse_NP_3_subset)
  
  ALL$EXP_GROUP<-factor(ALL$EXP_GROUP,
                        levels=c("MPRA_P","MPRA_NP_1_3","MPRA_NP_3"),
                        ordered=T)
  
  cat("ALL_0\n")
  cat(str(ALL))
  cat("\n")
  cat(str(unique(ALL$carried_variants)))
  cat("\n")
  cat(sprintf(as.character(names(summary(ALL$EXP_GROUP)))))
  cat("\n")
  cat(sprintf(as.character(summary(ALL$EXP_GROUP))))
  cat("\n")
  
  ALL_subset<-ALL[which(ALL$Label_2 == "ASSAYED_VARIANT"),]
 
  cat("ALL_subset_0\n")
  cat(str(ALL_subset))
  cat("\n")
  cat(str(unique(ALL_subset$carried_variants)))
  cat("\n")
  cat(sprintf(as.character(names(summary(ALL_subset$EXP_GROUP)))))
  cat("\n")
  cat(sprintf(as.character(summary(ALL_subset$EXP_GROUP))))
  cat("\n")
  
  ##### enhancer UpSetR_like_plots ----
  
  
  path5<-paste(out,'UpSetR_like_plots','/', sep='')
  
  cat("path5\n")
  cat(sprintf(as.character(path5)))
  cat("\n")
  
  
  if (file.exists(path5)){
    
    
    
    
  } else {
    dir.create(file.path(path5))
    
  }
  
  
  #### E_Plus_ASE_CLASS ----
  
  
  path6<-paste(out,'UpSetR_like_plots','/','E_Plus_ASE','/', sep='')
  
  cat("path6\n")
  cat(sprintf(as.character(path6)))
  cat("\n")
  
  
  if (file.exists(path6)){
    
    
    
    
  } else {
    dir.create(file.path(path6))
    
  }
  
  ##### LOOP Cell_Type ----
  
  array_CT<-levels(ALL_subset$Cell_Type)
  
  
  cat("array_CT_0\n")
  cat(str(array_CT))
  cat("\n")
  
  
  for(i in 1:length(array_CT))
  {
    
    CT_sel<-array_CT[i]
    
    cat("--------------------------------------------------------------------------------------------------->\t")
    cat(sprintf(as.character(CT_sel)))
    cat("\n")
    
    path7<-paste(out,'UpSetR_like_plots','/','E_Plus_ASE','/',CT_sel,'/', sep='')
    
    cat("path7\n")
    cat(sprintf(as.character(path7)))
    cat("\n")
    
    
    if (file.exists(path7)){
      
      
      
      
    } else {
      dir.create(file.path(path7))
      
    }
    
    
    ALL_subset_CT_sel<-ALL_subset[which(ALL_subset$Cell_Type == CT_sel),]
    
    cat("ALL_subset_CT_sel_0\n")
    cat(str(ALL_subset_CT_sel))
    cat("\n")
    cat(str(unique(ALL_subset_CT_sel$carried_variants)))
    cat("\n")
    cat(sprintf(as.character(names(summary(ALL_subset_CT_sel$EXP_GROUP)))))
    cat("\n")
    cat(sprintf(as.character(summary(ALL_subset_CT_sel$EXP_GROUP))))
    cat("\n")
    
    ALL_subset_CT_sel_melt<-melt(ALL_subset_CT_sel,id.vars=c("carried_variants","VAR","chr","Cell_Type","EXP_GROUP"))
    
    cat("ALL_subset_CT_sel_melt\n")
    cat(str(ALL_subset_CT_sel_melt))
    cat("\n")
    cat(sprintf(as.character(names(summary(as.factor(ALL_subset_CT_sel_melt$variable))))))
    cat("\n")
    cat(sprintf(as.character(summary(as.factor(ALL_subset_CT_sel_melt$variable)))))
    cat("\n")
    cat(sprintf(as.character(names(summary(ALL_subset_CT_sel_melt$EXP_GROUP)))))
    cat("\n")
    cat(sprintf(as.character(summary(ALL_subset_CT_sel_melt$EXP_GROUP))))
    cat("\n")
    
    
    
    ALL_subset_CT_sel_melt_sel_E_Plus_ASE_CLASS<-ALL_subset_CT_sel_melt[which(ALL_subset_CT_sel_melt$variable == "E_Plus_ASE_CLASS"),]
    
    
    cat("ALL_subset_CT_sel_melt_sel_E_Plus_ASE_CLASS\n")
    cat(str(ALL_subset_CT_sel_melt_sel_E_Plus_ASE_CLASS))
    cat("\n")
    
    
    ALL_subset_CT_sel_melt_sel_E_Plus_ASE_CLASS$value<-factor(as.character(ALL_subset_CT_sel_melt_sel_E_Plus_ASE_CLASS$value),
                                                           levels = c("0","1","2","3 or >"),
                                                           ordered=T)
    
    
    cat("ALL_subset_CT_sel_melt_sel_E_Plus_ASE_CLASS_2\n")
    cat(str(ALL_subset_CT_sel_melt_sel_E_Plus_ASE_CLASS))
    cat("\n")
    
    
    value_arrays<-levels(ALL_subset_CT_sel_melt_sel_E_Plus_ASE_CLASS$value)
    
    cat("value_arrays\n")
    cat(str(value_arrays))
    cat("\n")
    
    for(iteration_value_arrays in 1:length(value_arrays))
    {
      
      value_arrays_sel<-value_arrays[iteration_value_arrays]
      
      cat("-------->\t")
      cat(sprintf(as.character(value_arrays_sel)))
      cat("\n")
      
      
      if(value_arrays_sel == "0")
      {
        
        accepted_values<-"0"
      }
      if(value_arrays_sel == "1")
      {
        
        accepted_values<-c("1","2","3 or >")
      }
      if(value_arrays_sel == "2")
      {
        
        accepted_values<-c("2","3 or >")
      }
      if(value_arrays_sel == "3 or >")
      {
        
        accepted_values<-c("3 or >")
      }
      
      
      
      
      df_sel<-ALL_subset_CT_sel_melt_sel_E_Plus_ASE_CLASS[which(ALL_subset_CT_sel_melt_sel_E_Plus_ASE_CLASS$value%in%accepted_values),]
      
      cat("df_sel\n")
      cat(str(df_sel))
      cat("\n")
      cat(sprintf(as.character(names(summary(df_sel$EXP_GROUP)))))
      cat("\n")
      cat(sprintf(as.character(summary(df_sel$EXP_GROUP))))
      cat("\n")
      
      if(dim(df_sel)[1] >0)
      {
        
        df_sel_to_print<-df_sel
        
        df_sel_to_print$Active_tiles<-value_arrays_sel
        
        setwd(path7)
        
        write.table(df_sel_to_print,file=paste('ACTIVE_TILES_AT_LEAST_',value_arrays_sel,'.tsv'), sep="\t", quote=F, row.names = F)
        
       
        
        
        df_sel.dt<-data.table(df_sel,key=c("carried_variants","VAR","chr","Cell_Type"))
        
        cat("df_sel.dt\n")
        cat(str(df_sel.dt))
        cat("\n")
        
        
        
        
        df_sel_EXP_GROUP_summarised<-as.data.frame(df_sel.dt[,.(EXP_GROUP_string=paste(EXP_GROUP,collapse="|")), by=key(df_sel.dt)],stringsAsFactors=F)
        
        cat("df_sel_EXP_GROUP_summarised\n")
        cat(str(df_sel_EXP_GROUP_summarised))
        cat("\n")
        
        
        df_sel_EXP_GROUP_summarised$Classif_DEF<-"NA"
        
        df_sel_EXP_GROUP_summarised$Classif_DEF<-factor(df_sel_EXP_GROUP_summarised$EXP_GROUP_string,
                                                        levels = c("MPRA_P|MPRA_NP_1_3|MPRA_NP_3",
                                                                   "MPRA_P|MPRA_NP_1_3","MPRA_P|MPRA_NP_3","MPRA_NP_1_3|MPRA_NP_3",
                                                                   "MPRA_P","MPRA_NP_1_3","MPRA_NP_3"),
                                                        ordered=T)
        
        cat(sprintf(as.character(names(summary(df_sel_EXP_GROUP_summarised$Classif_DEF)))))
        cat("\n")
        cat(sprintf(as.character(summary(df_sel_EXP_GROUP_summarised$Classif_DEF))))
        cat("\n")
        
        df_sel_EXP_GROUP_summarised<-droplevels(df_sel_EXP_GROUP_summarised)
        
        
        cat("df_sel_EXP_GROUP_summarised_droplevels\n")
        cat(str(df_sel_EXP_GROUP_summarised))
        cat("\n")
        
        cat(sprintf(as.character(names(summary(df_sel_EXP_GROUP_summarised$Classif_DEF)))))
        cat("\n")
        cat(sprintf(as.character(summary(df_sel_EXP_GROUP_summarised$Classif_DEF))))
        cat("\n")
        
        
        
        
        
        #### Freq table ----
        
        
        
        df_sel_EXP_GROUP_summarised_Classif_DEF.dt<-data.table(df_sel_EXP_GROUP_summarised, key=c("Classif_DEF"))
        
        
        df_sel_EXP_GROUP_summarised_Classif_DEF_Fq<-as.data.frame(df_sel_EXP_GROUP_summarised_Classif_DEF.dt[,.N,by=key(df_sel_EXP_GROUP_summarised_Classif_DEF.dt)], stringsAsFactors=F)
        
        colnames(df_sel_EXP_GROUP_summarised_Classif_DEF_Fq)[which(colnames(df_sel_EXP_GROUP_summarised_Classif_DEF_Fq) == "N")]<-"instances"
        
        cat("df_sel_EXP_GROUP_summarised_Classif_DEF_Fq_0\n")
        cat(str(df_sel_EXP_GROUP_summarised_Classif_DEF_Fq))
        cat("\n")
        
        
        
        df_sel_EXP_GROUP_summarised_Classif_DEF_TOTAL.dt<-data.table(df_sel_EXP_GROUP_summarised, key=c("Classif_DEF"))
        
        
        df_sel_EXP_GROUP_summarised_Classif_DEF_TOTAL<-as.data.frame(df_sel_EXP_GROUP_summarised_Classif_DEF_TOTAL.dt[,.N,by=key(df_sel_EXP_GROUP_summarised_Classif_DEF_TOTAL.dt)], stringsAsFactors=F)
        
        colnames(df_sel_EXP_GROUP_summarised_Classif_DEF_TOTAL)[which(colnames(df_sel_EXP_GROUP_summarised_Classif_DEF_TOTAL) == "N")]<-"TOTAL"
        
        cat("df_sel_EXP_GROUP_summarised_Classif_DEF_TOTAL_\n")
        cat(str(df_sel_EXP_GROUP_summarised_Classif_DEF_TOTAL))
        cat("\n")
        
        df_sel_EXP_GROUP_summarised_Classif_DEF_Fq<-merge(df_sel_EXP_GROUP_summarised_Classif_DEF_Fq,
                                                          df_sel_EXP_GROUP_summarised_Classif_DEF_TOTAL,
                                                          by="Classif_DEF")
        
        df_sel_EXP_GROUP_summarised_Classif_DEF_Fq$Perc<-round(100*(df_sel_EXP_GROUP_summarised_Classif_DEF_Fq$instances/df_sel_EXP_GROUP_summarised_Classif_DEF_Fq$TOTAL),1)
        
        
        cat("df_sel_EXP_GROUP_summarised_Classif_DEF_Fq_1\n")
        cat(str(df_sel_EXP_GROUP_summarised_Classif_DEF_Fq))
        cat("\n")
        
        
        
        step<-round(max(df_sel_EXP_GROUP_summarised_Classif_DEF_Fq$TOTAL)/10,0)
        
        cat("--------->\t")
        cat(sprintf(as.character(step)))
        cat("\n")
        
        if(step == 0){
          
          step=1
        }
        
        
        breaks.y<-seq(0,max(df_sel_EXP_GROUP_summarised_Classif_DEF_Fq$TOTAL)+step, by=step)
        labels.y<-as.character(breaks.y)
        
        
        cat(sprintf(as.character(breaks.y)))
        cat("\n")
        
        
        
        cat("---------------------------->df_sel_EXP_GROUP_summarised_Classif_DEF_TOTAL\n")
        cat(str(df_sel_EXP_GROUP_summarised_Classif_DEF_TOTAL))
        cat("\n")
        
        graph<-ggplot(data=df_sel_EXP_GROUP_summarised_Classif_DEF_TOTAL,
                      aes(x=Classif_DEF,
                          y=TOTAL)) +
          geom_bar(stat="identity",colour='black')+
          theme_bw()+
          scale_y_continuous(name=paste("Variants with at least",value_arrays_sel,"ACTIVE TILES",sep=" "),breaks=breaks.y,labels=labels.y,
                             limits=c(breaks.y[1],breaks.y[length(breaks.y)]))+
          scale_x_discrete(name=NULL, drop=F)+
          theme(legend.position="none")+
          theme_classic()+
          theme(axis.title.y=element_text(size=18, family="sans"),
                axis.text.y=element_text(angle=0,size=18, color="black", family="sans"),
                axis.text.x=element_text(angle=45,vjust=1,hjust=1,size=6, color="black", family="sans"))+
          ggeasy::easy_center_title()
        
        ### subgraph
        
        EXP_GROUP.dt<-data.table(df_sel,key=c("EXP_GROUP","Cell_Type"))
        
        cat("EXP_GROUP.dt\n")
        cat(str(EXP_GROUP.dt))
        cat("\n")
        
        
        
        EXP_GROUP_df<-as.data.frame(EXP_GROUP.dt[,.N, by=key(EXP_GROUP.dt)],stringsAsFactors=F)
        colnames(EXP_GROUP_df)[which(colnames(EXP_GROUP_df) == "N")]<-"Total"
        
        cat("EXP_GROUP_df_0\n")
        cat(str(EXP_GROUP_df))
        cat("\n")
        
        
       
        
        EXP_GROUP_df$Total[is.na(EXP_GROUP_df$Total)]<-0
        
        
        cat("EXP_GROUP_df_WITHOUT_CTRLS0\n")
        cat(str(EXP_GROUP_df))
        cat("\n")
        
        
        
        step<-round(max(EXP_GROUP_df$Total)/5,0)
        
        cat("--------->\t")
        cat(sprintf(as.character(step)))
        cat("\n")
        
        if(step == 0){
          
          step=1
        }
        breaks.x<-rev(c(seq(0,max(EXP_GROUP_df$Total)+step, by=step)))
        labels.x<-as.character(breaks.x)
        
        
        cat(sprintf(as.character(breaks.x)))
        cat("\n")
        
        levels_CT<-rev(levels(EXP_GROUP_df$EXP_GROUP))
        
        cat("levels_CT_2\n")
        cat(str(levels_CT))
        cat("\n")
        
        EXP_GROUP_df$EXP_GROUP<-factor(EXP_GROUP_df$EXP_GROUP,
                                       levels=levels_CT,
                                       ordered=T)
        
        # EXP_GROUP_df<-EXP_GROUP_df[order(EXP_GROUP_df$EXP_GROUP),]
        
        cat("EXP_GROUP_df_2\n")
        cat(str(EXP_GROUP_df))
        cat("\n")
        
        
        subgraph<-ggplot(data=EXP_GROUP_df,
                         aes(y=EXP_GROUP,
                             x=Total)) +
          geom_bar(stat="identity",fill='steelblue', colour='steelblue')+
          theme_bw()+
          theme(axis.title.y=element_text(size=16, family="sans"),
                axis.text.y=element_text(angle=0,size=16, color="black", family="sans"),
                axis.text.x=element_text(angle=0,size=10, color="black", family="sans"))+
          scale_x_reverse(name=paste("Variants with at least",value_arrays_sel,"ACTIVE TILES",sep=" "),breaks=breaks.x,labels=labels.x,
                          limits=c(breaks.x[1],breaks.x[length(breaks.x)]))+
          scale_y_discrete(position="right",name=NULL, drop=F)+
          theme(legend.position="none")+
          theme_classic()+
          ggeasy::easy_center_title()
        
        cat("subgraph DONE\n")
        
        
        
        
        
        graph_FINAL<-plot_grid(NULL,graph,subgraph,NULL,
                               nrow = 2,
                               ncol=2,
                               rel_heights = c(1, 0.25),
                               rel_widths=c(0.5,1))
        
        
        
        
        
        setwd(path7)
        
        svglite(paste('ACTIVE_TILES_AT_LEAST_',value_arrays_sel,'.svg',sep=''), width = 8, height = 8)
        print(graph_FINAL)
        dev.off()
        
        
        
        # ################################################################################################## HERE HERE HERE -----
        # quit(status = 1)
      }#dim(df_sel)[1] >0
      
      
      
      
     
    }#iteration_value_arrays in 1:length(value_arrays)
    
    
    
    
    
  }#i i in 1:length(array_CT)
  
  
  
  # #######################################
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
    make_option(c("--type"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--out"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--KEY_collapse_P"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--KEY_collapse_NP_1_3"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--KEY_collapse_NP_3"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--CSQ_colors"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--df_Cell_colors"), type="character", default=NULL, 
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
  
  
 
  UpSetR_like_E_Plus_ASE(opt)
 
  
}
  
  
  
 

###########################################################################

system.time( main() )
