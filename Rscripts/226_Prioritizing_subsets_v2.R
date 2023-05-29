suppressMessages(library("optparse", lib.loc = "/nfs/users/nfs_m/mt19/sOFTWARE/R_libs"))

suppressMessages(library("data.table", lib.loc="/nfs/users/nfs_m/mt19/sOFTWARE/R_libs/"))
suppressMessages(library("plyr", lib.loc="/nfs/users/nfs_m/mt19/sOFTWARE/R_libs/"))
suppressMessages(library("dplyr", lib.loc="/nfs/users/nfs_m/mt19/sOFTWARE/R_libs/"))





opt = NULL

Prioritizing_subsets = function(option_list)
{
  opt_in = option_list
  opt <<- option_list
  
  cat("All options:\n")
  printList(opt)
  
  #### READ and transform type ----
  
  type_1 = opt$type_1
  
  cat("TYPE_\n")
  cat(sprintf(as.character(type_1)))
  cat("\n")
  
  #### READ and transform out ----
  
  out = opt$out
  
  cat("OUT_\n")
  cat(sprintf(as.character(out)))
  cat("\n")
  
  #### READ and transform PP_threshold ----
  
  PP_threshold = opt$PP_threshold
  
  cat("PP_threshold_\n")
  cat(sprintf(as.character(PP_threshold)))
  cat("\n")
  
  #### dB_ALL ALWAYS exists ----
  
  
  dB_ALL<-as.data.frame(fread(file = opt$dB_ALL,
                               sep="\t",
                               header=T), stringsAsFactors=F)
  
  dB_ALL$Allelic_Series_ID<-paste(dB_ALL$phenotype,dB_ALL$block_no, sep="__")
  
  
  cat("dB_ALL\n")
  cat(str(dB_ALL))
  cat("\n")
  
  #### Rank_file ALWAYS exists ----
  
  Rank_file<-as.data.frame(fread(file = opt$Rank_file,
                                 sep="\t",
                                 header=T), stringsAsFactors=F)
  
  Rank_file$Rank<-as.numeric(row.names(Rank_file))
  
  cat("Rank_file_0\n")
  cat(str(Rank_file))
  cat("\n")
  
  Rank_file<-Rank_file[which(Rank_file$phenotype != "NA"),]
  
  
  cat("Rank_file_NO_NA\n")
  cat(str(Rank_file))
  cat("\n")
  
  # Rank_file<-merge(Rank_file,
  #                  GENE.TABLE_equivalence,
  #                  by="ensembl_gene_id",
  #                  all.x=T)
  # 
  # cat("Rank_file\n")
  # cat(str(Rank_file))
  # cat("\n")
  
  Rank_file<-Rank_file[order(Rank_file$Rank, decreasing = F),]
  
  cat("Rank_file_ordered\n")
  cat(str(Rank_file))
  cat("\n")
  
  #quit(status = 1)
  
  
  #### Select variants ----
  
  
  Rank_file_MAX<-as.data.frame(setDT(Rank_file)[, .SD[which.min(Rank)], by=c("phenotype","block_no")])
  
  Rank_file_MAX$Allelic_Series_ID<-paste(Rank_file_MAX$phenotype,Rank_file_MAX$block_no, sep="__")
  
  
  cat("Rank_file_MAX_ordered\n")
  cat(str(Rank_file_MAX))
  cat("\n")
  
  
 
  ##### VARS
  
  AS<-unique(Rank_file_MAX$Allelic_Series_ID[order(Rank_file_MAX$Rank, decreasing = F)])
  
  cat("AS\n")
  cat(str(AS))
  cat("\n")
  
  dB_ALL_sel_check<-dB_ALL[which(dB_ALL$Allelic_Series_ID%in%AS),]
  
  
  cat("dB_ALL_sel_check\n")
  cat(str(dB_ALL_sel_check))
  cat("\n")
 
  
  
  ##### Apply the probability threshold----
  
  dB_ALL_sel_check_PP_threshold<-unique(dB_ALL_sel_check[which(dB_ALL_sel_check$finemap_prob >= PP_threshold),])
  
  cat("dB_ALL_sel_check_PP_threshold\n")
  cat(str(dB_ALL_sel_check_PP_threshold))
  cat("\n")
  
  
  VARS_check<-unique(dB_ALL_sel_check_PP_threshold$VAR)
  
  cat("VARS_check\n")
  cat(str(VARS_check))
  cat("\n")
  
  #### SAVE ----
  
  setwd(out)
  
  filename_1<-paste("Selection_",type_1,".tsv", sep='')
  
  write.table(dB_ALL_sel_check_PP_threshold,
              file=filename_1, sep="\t", quote=F, row.names = F)
  
  
  

  
}

Prioritizing_subsets_2 = function(option_list)
{
  opt_in = option_list
  opt <<- option_list
  
  cat("All options:\n")
  printList(opt)
  
  #### READ and transform type_2 ----
  
  type_2 = opt$type_2
  
  cat("type_2_\n")
  cat(sprintf(as.character(type_2)))
  cat("\n")
  
  #### READ and transform out ----
  
  out = opt$out
  
  cat("OUT_\n")
  cat(sprintf(as.character(out)))
  cat("\n")
  
  #### READ and transform PP_threshold ----
  
  PP_threshold = opt$PP_threshold
  
  cat("PP_threshold_\n")
  cat(sprintf(as.character(PP_threshold)))
  cat("\n")
  
  #### dB_ALL ALWAYS exists ----
  
  
  dB_ALL<-as.data.frame(fread(file = opt$dB_ALL,
                              sep="\t",
                              header=T), stringsAsFactors=F)
  
  dB_ALL$Allelic_Series_ID<-paste(dB_ALL$phenotype_2,dB_ALL$block_no, sep="__")
  
  
  cat("dB_ALL\n")
  cat(str(dB_ALL))
  cat("\n")
  
  #### Original_selection_Grouped_variants ALWAYS exists ----
  
  Original_selection_Grouped_variants<-read.table(file = opt$Original_selection_Grouped_variants,
                                 sep="\t",
                                 header=T,
                                 row.names=1, stringsAsFactors=F)
  
   
  cat("Original_selection_Grouped_variants_0\n")
  cat(str(Original_selection_Grouped_variants))
  cat("\n")
  
  
  Original_selection_Grouped_variants_Manuel_category<-unique(Original_selection_Grouped_variants[!is.na(Original_selection_Grouped_variants$Manuel_Category),])
  
  cat("Original_selection_Grouped_variants_Manuel_category_0\n")
  cat(str(Original_selection_Grouped_variants_Manuel_category))
  cat("\n")
  
  cat("Manuel_Category\n")
  cat(sprintf(as.character(names(summary(as.factor(Original_selection_Grouped_variants_Manuel_category$Manuel_Category))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(Original_selection_Grouped_variants_Manuel_category$Manuel_Category)))))
  cat("\n")
  
  #### categories_of_RRVV ----
  
  categories_of_RRVV = opt$categories_of_RRVV
  
  cat("categories_of_RRVV\n")
  cat(sprintf(as.character(categories_of_RRVV)))
  cat("\n")
  
  categories_of_RRVV_unstring<-unlist(strsplit(categories_of_RRVV, split=","))
  
  cat("categories_of_RRVV_unstring\n")
  cat(sprintf(as.character(categories_of_RRVV_unstring)))
  cat("\n")
  
  
  ##### Select the categories ----
  
  subset <-unique(Original_selection_Grouped_variants_Manuel_category[which(Original_selection_Grouped_variants_Manuel_category$Manuel_Category%in%categories_of_RRVV_unstring),])
  
  cat("subset_0\n")
  cat(str(subset))
  cat("\n")
  
  
  #quit(status = 1)
  
  
  ##### VARS
  
  dB_ALL_sel_check<-dB_ALL[which(dB_ALL$VAR%in%subset$VAR),]
  
  
  cat("dB_ALL_sel_check\n")
  cat(str(dB_ALL_sel_check))
  cat("\n")
  
  
  
  ##### Apply the probability threshold----
  
  dB_ALL_sel_check_PP_threshold<-unique(dB_ALL_sel_check[which(dB_ALL_sel_check$finemap_prob >= PP_threshold),])
  
  cat("dB_ALL_sel_check_PP_threshold\n")
  cat(str(dB_ALL_sel_check_PP_threshold))
  cat("\n")
  
  
  VARS_check<-unique(dB_ALL_sel_check_PP_threshold$VAR)
  
  cat("VARS_check\n")
  cat(str(VARS_check))
  cat("\n")
  
  #### SAVE ----
  
  setwd(out)
  
  filename_1<-paste("Selection_",type_2,".tsv", sep='')
  
  write.table(dB_ALL_sel_check_PP_threshold,
              file=filename_1, sep="\t", quote=F, row.names = F)
}

Prioritizing_subsets_3 = function(option_list)
{
  opt_in = option_list
  opt <<- option_list
  
  cat("All options:\n")
  printList(opt)
  
  #### READ and transform type ----
  
  type_3 = opt$type_3
  
  cat("TYPE_\n")
  cat(sprintf(as.character(type_3)))
  cat("\n")
  
  #### READ and transform out ----
  
  out = opt$out
  
  cat("OUT_\n")
  cat(sprintf(as.character(out)))
  cat("\n")
  
  #### READ and transform PP_threshold ----
  
  PP_threshold = opt$PP_threshold
  
  cat("PP_threshold_\n")
  cat(sprintf(as.character(PP_threshold)))
  cat("\n")
  
  #### dB_ALL ALWAYS exists ----
  
  
  dB_ALL<-as.data.frame(fread(file = opt$dB_ALL,
                              sep="\t",
                              header=T), stringsAsFactors=F)
  
  dB_ALL$Allelic_Series_ID<-paste(dB_ALL$phenotype,dB_ALL$block_no, sep="__")
  
  
  cat("dB_ALL\n")
  cat(str(dB_ALL))
  cat("\n")
  
  #### Rank_file_3 ALWAYS exists ----
  
  Rank_file_3<-as.data.frame(fread(file = opt$Rank_file_3,
                                 sep="\t",
                                 header=T), stringsAsFactors=F)
  
  Rank_file_3$Rank<-as.numeric(row.names(Rank_file_3))
  
  cat("Rank_file_3_0\n")
  cat(str(Rank_file_3))
  cat("\n")
  
  Rank_file_3<-Rank_file_3[which(Rank_file_3$phenotype != "NA"),]
  
  
  cat("Rank_file_3_NO_NA\n")
  cat(str(Rank_file_3))
  cat("\n")
  
  # Rank_file_3<-merge(Rank_file_3,
  #                  GENE.TABLE_equivalence,
  #                  by="ensembl_gene_id",
  #                  all.x=T)
  # 
  # cat("Rank_file_3\n")
  # cat(str(Rank_file_3))
  # cat("\n")
  
  Rank_file_3<-Rank_file_3[order(Rank_file_3$Rank, decreasing = F),]
  
  cat("Rank_file_3_ordered\n")
  cat(str(Rank_file_3))
  cat("\n")
  
  #quit(status = 1)
  
  
  #### Select variants ----
  
  
  Rank_file_3_MAX<-as.data.frame(setDT(Rank_file_3)[, .SD[which.min(Rank)], by=c("phenotype","block_no")])
  
  Rank_file_3_MAX$Allelic_Series_ID<-paste(Rank_file_3_MAX$phenotype,Rank_file_3_MAX$block_no, sep="__")
  
  
  cat("Rank_file_3_MAX_ordered\n")
  cat(str(Rank_file_3_MAX))
  cat("\n")
  
  
  
  ##### VARS
  
  AS<-unique(Rank_file_3_MAX$Allelic_Series_ID[order(Rank_file_3_MAX$Rank, decreasing = F)])
  
  cat("AS\n")
  cat(str(AS))
  cat("\n")
  
  indx.mono.178<-grep("mono__178",AS)
  
  cat("indx.mono.178\n")
  cat(sprintf(as.character(indx.mono.178)))
  cat("\n")
  
  #### READ and transform PP_threshold ----
  
  cut_3 = opt$cut_3
  
  cat("cut_3_\n")
  cat(sprintf(as.character(cut_3)))
  cat("\n")
  
  AS_adapted<-AS[1:cut_3]
  
  
  cat("AS_adapted\n")
  cat(str(AS_adapted))
  cat("\n")
  
  dB_ALL_sel_check<-dB_ALL[which(dB_ALL$Allelic_Series_ID%in%AS_adapted),]
  
  
  cat("dB_ALL_sel_check\n")
  cat(str(dB_ALL_sel_check))
  cat("\n")
  
  #quit(status = 1)
  
  ##### Apply the probability threshold----
  
  dB_ALL_sel_check_PP_threshold<-unique(dB_ALL_sel_check[which(dB_ALL_sel_check$finemap_prob >= PP_threshold),])
  
  cat("dB_ALL_sel_check_PP_threshold\n")
  cat(str(dB_ALL_sel_check_PP_threshold))
  cat("\n")
  
  
  VARS_check<-unique(dB_ALL_sel_check_PP_threshold$VAR)
  
  cat("VARS_check\n")
  cat(str(VARS_check))
  cat("\n")
  
  #### SAVE ----
  
  setwd(out)
  
  filename_1<-paste("Selection_",type_3,"_",cut_3,".tsv", sep='')
  
  write.table(dB_ALL_sel_check_PP_threshold,
              file=filename_1, sep="\t", quote=F, row.names = F)
  
}

Prioritizing_subsets_4 = function(option_list)
{
  opt_in = option_list
  opt <<- option_list
  
  cat("All options:\n")
  printList(opt)
  
  #### READ and transform type ----
  
  type_4 = opt$type_4
  
  cat("TYPE_\n")
  cat(sprintf(as.character(type_4)))
  cat("\n")
  
  #### READ and transform out ----
  
  out = opt$out
  
  cat("OUT_\n")
  cat(sprintf(as.character(out)))
  cat("\n")
  
  #### READ and transform PP_threshold ----
  
  PP_threshold = opt$PP_threshold
  
  cat("PP_threshold_\n")
  cat(sprintf(as.character(PP_threshold)))
  cat("\n")
  
  #### dB_ALL ALWAYS exists ----
  
  
  dB_ALL<-as.data.frame(fread(file = opt$dB_ALL,
                              sep="\t",
                              header=T), stringsAsFactors=F)
  
  dB_ALL$Allelic_Series_ID<-paste(dB_ALL$phenotype,dB_ALL$block_no, sep="__")
  
  
  cat("dB_ALL\n")
  cat(str(dB_ALL))
  cat("\n")
  
  #### Rank_file_4 ALWAYS exists ----
  
  Rank_file_4<-as.data.frame(fread(file = opt$Rank_file_4,
                                   sep="\t",
                                   header=T), stringsAsFactors=F)
  
  Rank_file_4$Rank<-as.numeric(row.names(Rank_file_4))
  
  cat("Rank_file_4_0\n")
  cat(str(Rank_file_4))
  cat("\n")
  
  Rank_file_4<-Rank_file_4[which(Rank_file_4$phenotype != "NA"),]
  
  
  cat("Rank_file_4_NO_NA\n")
  cat(str(Rank_file_4))
  cat("\n")
  
  # Rank_file_3<-merge(Rank_file_3,
  #                  GENE.TABLE_equivalence,
  #                  by="ensembl_gene_id",
  #                  all.x=T)
  # 
  # cat("Rank_file_3\n")
  # cat(str(Rank_file_3))
  # cat("\n")
  
  Rank_file_4<-Rank_file_4[order(Rank_file_4$Rank, decreasing = F),]
  
  cat("Rank_file_4_ordered\n")
  cat(str(Rank_file_4))
  cat("\n")
  
  #quit(status = 1)
  
  
  
  Rank_file_4_MAX<-as.data.frame(setDT(Rank_file_4)[, .SD[which.min(Rank)], by=c("phenotype","block_no")])
  
  Rank_file_4_MAX$Allelic_Series_ID<-paste(Rank_file_4_MAX$phenotype,Rank_file_4_MAX$block_no, sep="__")
  
  
  cat("Rank_file_4_MAX_ordered\n")
  cat(str(Rank_file_4_MAX))
  cat("\n")
  
  ##### Open_targets_GLOBAL_file -----
  
  Open_targets_GLOBAL_file<-as.data.frame(fread(file = opt$Open_targets_GLOBAL_file,
                              sep="\t",
                              header=T), stringsAsFactors=F)
  
  
  cat("Open_targets_GLOBAL_file\n")
  cat(str(Open_targets_GLOBAL_file))
  cat("\n")
  
  Open_targets_GLOBAL_file_anchored<-Open_targets_GLOBAL_file[which(Open_targets_GLOBAL_file$VAR_IN_AS == "ANCHORED_IN_AS"),]
  
  cat("Open_targets_GLOBAL_file_anchored\n")
  cat(str(Open_targets_GLOBAL_file_anchored))
  cat("\n")
  
  ###### Conditions: ANCHORED IN AN INTERESTING DISEASE-----
  
  blood_diseases<-c("anemia","anemia (disease)","hemolytic anemia","megaloblastic anemia (disease)","leukocyte disease","Leukocytosis","Leukopenia","Rare genetic coagulation disorder",
                    "Rare hemorrhagic disorder due to a  constitutional platelet anomaly","response to anticoagulant","response to platelet aggregation inhibitor",
                    "thrombocytopenia","Thrombotic disease","Thrombophlebitis","venous thromboembolism","deep vein thrombosis","Abnormality of blood and blood-forming tissues",
                    "Elevated erythrocyte sedimentation rate","blood coagulation disease","Eosinophilia","hematologic disease","Hemoglobinopathy","Pancytopenia","pernicious anemia","primary thrombocytopenia")
  
  inflammatory_diseases<-c("allergy","asthma","Bronchospasm","psoriasis vulgaris","psoriatic arthritis","purpura (disease)","type II diabetes mellitus","type I diabetes mellitus",
                           "Crohn's disease","lupus erythematosus","inflammatory bowel disease","rheumatoid arthritis","systemic lupus erythematosus","systemic inflammatory response syndrome","multiple sclerosis",
                           "Inflammation","inflammation","Sarcoidosis","systemic inflammatory response syndrome")
  
  immunodeficiencies<-c("immunodeficiency disease","abnormality of the immune system","Immunodeficiency",'immune system disease')
  
  miscellaneous<-c("atherosclerosis","coronary atherosclerosis","gastrointestinal hemorrhage","Hematemesis","Hematuria","hemorrhage","Hemorrhage","hemorrhagic disease","Hematochezia","Hemoptysis",
                   "Sepsis","Abnormal glucose homeostasis","arthritis","lymphatic system disease","plasma protein metabolism disease","bacteriemia",
                   "hyperlipidemia","Hypercalcemia",'Hyperkalemia','Hypokalemia',"Hypocalcemia",'Hypernatremia','Hyponatremia',"Hypovolemia","Hypervolemia","Hypoglycemia","Polymenorrhea","polyuria")
  
  leukemias<-c("lymphoid leukemia","non-Hodgkins lymphoma","leukemia")
  
  INTERESTING_DISEASE<-unique(c(blood_diseases,inflammatory_diseases,immunodeficiencies,miscellaneous,leukemias))
  
  cat("INTERESTING_DISEASE\n")
  cat(str(INTERESTING_DISEASE))
  cat("\n")
  
  Open_targets_GLOBAL_file_anchored_filtered<-Open_targets_GLOBAL_file_anchored[which(Open_targets_GLOBAL_file_anchored$Disease_or_phenotype%in%INTERESTING_DISEASE),]
  
  cat("Open_targets_GLOBAL_file_anchored_filtered\n")
  cat(str(Open_targets_GLOBAL_file_anchored_filtered))
  cat("\n")
  
  
  ##### filter Rank_file_4_MAX ----
  
  Rank_file_4_MAX_filtered<-Rank_file_4_MAX[which(Rank_file_4_MAX$Allelic_Series_ID%in%Open_targets_GLOBAL_file_anchored_filtered$Allelic_Series_ID),]
 
   cat("Rank_file_4_MAX_filtered\n")
  cat(str(Rank_file_4_MAX_filtered))
  cat("\n")
  
  
  ###### Retrieve variants ----
  
  AS<-unique(Rank_file_4_MAX_filtered$Allelic_Series_ID[order(Rank_file_4_MAX_filtered$Rank, decreasing = F)])
  
  cat("AS\n")
  cat(str(AS))
  cat("\n")
  
  dB_ALL_sel_check<-dB_ALL[which(dB_ALL$Allelic_Series_ID%in%AS),]
  
  
  cat("dB_ALL_sel_check\n")
  cat(str(dB_ALL_sel_check))
  cat("\n")
  
  
  
  ##### Apply the probability threshold----
  
  dB_ALL_sel_check_PP_threshold_type_4<-unique(dB_ALL_sel_check[which(dB_ALL_sel_check$finemap_prob >= PP_threshold),])
  
  cat("dB_ALL_sel_check_PP_threshold_type_4\n")
  cat(str(dB_ALL_sel_check_PP_threshold_type_4))
  cat("\n")
  
  
  VARS_check<-unique(dB_ALL_sel_check_PP_threshold_type_4$VAR)
  
  cat("VARS_check\t")
  cat(sprintf(as.character(type_4)))
  cat("\t")
  cat(sprintf(as.character(length(VARS_check))))
  cat("\n")
  
  setwd(out)
  
  filename_1<-paste("Selection_",type_4,".tsv", sep='')
  
  write.table(dB_ALL_sel_check_PP_threshold_type_4,
              file=filename_1, sep="\t", quote=F, row.names = F)
  
  
  
  #### Rank_file_5 ALWAYS exists ----
  
  type_5 = opt$type_5
  
  cat("TYPE_\n")
  cat(sprintf(as.character(type_5)))
  cat("\n")
  
  Rank_file_5<-as.data.frame(fread(file = opt$Rank_file_5,
                                   sep="\t",
                                   header=T), stringsAsFactors=F)
  
  Rank_file_5$Rank<-as.numeric(row.names(Rank_file_5))
  
  cat("Rank_file_5_0\n")
  cat(str(Rank_file_5))
  cat("\n")
  
  Rank_file_5<-Rank_file_5[which(Rank_file_5$phenotype != "NA"),]
  
  
  cat("Rank_file_5_NO_NA\n")
  cat(str(Rank_file_5))
  cat("\n")
  
  # Rank_file_3<-merge(Rank_file_3,
  #                  GENE.TABLE_equivalence,
  #                  by="ensembl_gene_id",
  #                  all.x=T)
  # 
  # cat("Rank_file_3\n")
  # cat(str(Rank_file_3))
  # cat("\n")
  
  Rank_file_5<-Rank_file_5[order(Rank_file_5$Rank, decreasing = F),]
  
  cat("Rank_file_5_ordered\n")
  cat(str(Rank_file_5))
  cat("\n")
  
  #quit(status = 1)
  
  
 
  
  
  Rank_file_5_MAX<-as.data.frame(setDT(Rank_file_5)[, .SD[which.min(Rank)], by=c("phenotype","block_no")])
  
  Rank_file_5_MAX$Allelic_Series_ID<-paste(Rank_file_5_MAX$phenotype,Rank_file_5_MAX$block_no, sep="__")
  
  
  cat("Rank_file_5_MAX_ordered\n")
  cat(str(Rank_file_5_MAX))
  cat("\n")
  
  ##### filter Rank_file_5_MAX ----
  
  Rank_file_5_MAX_filtered<-Rank_file_5_MAX[which(Rank_file_5_MAX$Allelic_Series_ID%in%Open_targets_GLOBAL_file_anchored_filtered$Allelic_Series_ID),]
  
  cat("Rank_file_5_MAX_filtered_1\n")
  cat(str(Rank_file_5_MAX_filtered))
  cat("\n")
  
  Rank_file_5_MAX_filtered<-Rank_file_5_MAX_filtered[-which(Rank_file_5_MAX_filtered$Allelic_Series_ID%in%dB_ALL_sel_check_PP_threshold_type_4$Allelic_Series_ID),]
  
  cat("Rank_file_5_MAX_filtered_2\n")
  cat(str(Rank_file_5_MAX_filtered))
  cat("\n")
  
  
  ###### Retrieve variants ----
  
  AS<-unique(Rank_file_5_MAX_filtered$Allelic_Series_ID[order(Rank_file_5_MAX_filtered$Rank, decreasing = F)])
  
  cat("AS\n")
  cat(str(AS))
  cat("\n")
  
  dB_ALL_sel_check<-dB_ALL[which(dB_ALL$Allelic_Series_ID%in%AS),]
  
  
  cat("dB_ALL_sel_check\n")
  cat(str(dB_ALL_sel_check))
  cat("\n")
  
  
  
  ##### Apply the probability threshold----
  
  dB_ALL_sel_check_PP_threshold_type_5<-unique(dB_ALL_sel_check[which(dB_ALL_sel_check$finemap_prob >= PP_threshold),])
  
  cat("dB_ALL_sel_check_PP_threshold_type_5\n")
  cat(str(dB_ALL_sel_check_PP_threshold_type_5))
  cat("\n")
  
  
  VARS_check<-unique(dB_ALL_sel_check_PP_threshold_type_5$VAR)
  
  cat("VARS_check\t")
  cat(sprintf(as.character(type_5)))
  cat("\t")
  cat(sprintf(as.character(length(VARS_check))))
  cat("\n")
  
  setwd(out)
  
  filename_1<-paste("Selection_",type_5,".tsv", sep='')
  
  write.table(dB_ALL_sel_check_PP_threshold_type_5,
              file=filename_1, sep="\t", quote=F, row.names = F)
  
  
  
  
  
  
  #### Rank_file_6 ALWAYS exists ----
  
  type_6 = opt$type_6
  
  cat("TYPE_\n")
  cat(sprintf(as.character(type_6)))
  cat("\n")
  
  Rank_file_6<-as.data.frame(fread(file = opt$Rank_file_6,
                                   sep="\t",
                                   header=T), stringsAsFactors=F)
  
  Rank_file_6$Rank<-as.numeric(row.names(Rank_file_6))
  
  cat("Rank_file_6_0\n")
  cat(str(Rank_file_6))
  cat("\n")
  
  Rank_file_6<-Rank_file_6[which(Rank_file_6$phenotype != "NA"),]
  
  
  cat("Rank_file_6_NO_NA\n")
  cat(str(Rank_file_6))
  cat("\n")
  
  # Rank_file_3<-merge(Rank_file_3,
  #                  GENE.TABLE_equivalence,
  #                  by="ensembl_gene_id",
  #                  all.x=T)
  # 
  # cat("Rank_file_3\n")
  # cat(str(Rank_file_3))
  # cat("\n")
  
  Rank_file_6<-Rank_file_6[order(Rank_file_6$Rank, decreasing = F),]
  
  cat("Rank_file_6_ordered\n")
  cat(str(Rank_file_6))
  cat("\n")
  
  #quit(status = 1)
  
  
  
  
  
  Rank_file_6_MAX<-as.data.frame(setDT(Rank_file_6)[, .SD[which.min(Rank)], by=c("phenotype","block_no")])
  
  Rank_file_6_MAX$Allelic_Series_ID<-paste(Rank_file_6_MAX$phenotype,Rank_file_6_MAX$block_no, sep="__")
  
  
  cat("Rank_file_6_MAX_ordered\n")
  cat(str(Rank_file_6_MAX))
  cat("\n")
  
  ##### filter Rank_file_6_MAX ----
  
  Rank_file_6_MAX_filtered<-Rank_file_6_MAX[which(Rank_file_6_MAX$Allelic_Series_ID%in%Open_targets_GLOBAL_file_anchored_filtered$Allelic_Series_ID),]
  
  cat("Rank_file_6_MAX_filtered_1\n")
  cat(str(Rank_file_6_MAX_filtered))
  cat("\n")
  
  Rank_file_6_MAX_filtered<-Rank_file_6_MAX_filtered[-which(Rank_file_6_MAX_filtered$Allelic_Series_ID%in%dB_ALL_sel_check_PP_threshold_type_4$Allelic_Series_ID),]
  
  cat("Rank_file_6_MAX_filtered_2\n")
  cat(str(Rank_file_6_MAX_filtered))
  cat("\n")
  
  Rank_file_6_MAX_filtered<-Rank_file_6_MAX_filtered[-which(Rank_file_6_MAX_filtered$Allelic_Series_ID%in%dB_ALL_sel_check_PP_threshold_type_5$Allelic_Series_ID),]
  
  cat("Rank_file_6_MAX_filtered_3\n")
  cat(str(Rank_file_6_MAX_filtered))
  cat("\n")
  
  
  ###### Retrieve variants ----
  
  AS<-unique(Rank_file_6_MAX_filtered$Allelic_Series_ID[order(Rank_file_6_MAX_filtered$Rank, decreasing = F)])
  
  cat("AS\n")
  cat(str(AS))
  cat("\n")
  
  dB_ALL_sel_check<-dB_ALL[which(dB_ALL$Allelic_Series_ID%in%AS),]
  
  
  cat("dB_ALL_sel_check\n")
  cat(str(dB_ALL_sel_check))
  cat("\n")
  
  
  
  ##### Apply the probability threshold----
  
  dB_ALL_sel_check_PP_threshold_type_6<-unique(dB_ALL_sel_check[which(dB_ALL_sel_check$finemap_prob >= PP_threshold),])
  
  cat("dB_ALL_sel_check_PP_threshold_type_6\n")
  cat(str(dB_ALL_sel_check_PP_threshold_type_6))
  cat("\n")
  
  
  VARS_check<-unique(dB_ALL_sel_check_PP_threshold_type_6$VAR)
  
  cat("VARS_check\t")
  cat(sprintf(as.character(type_6)))
  cat("\t")
  cat(sprintf(as.character(length(VARS_check))))
  cat("\n")
  
  setwd(out)
  
  filename_1<-paste("Selection_",type_6,".tsv", sep='')
  
  write.table(dB_ALL_sel_check_PP_threshold_type_6,
              file=filename_1, sep="\t", quote=F, row.names = F)
  
  
  
  
  
  
  #### Rank_file_7 ALWAYS exists ----
  
  type_7 = opt$type_7
  
  cat("TYPE_\n")
  cat(sprintf(as.character(type_7)))
  cat("\n")
  
  Rank_file_7<-as.data.frame(fread(file = opt$Rank_file_7,
                                   sep="\t",
                                   header=T), stringsAsFactors=F)
  
  Rank_file_7$Rank<-as.numeric(row.names(Rank_file_7))
  
  cat("Rank_file_7_0\n")
  cat(str(Rank_file_7))
  cat("\n")
  
  Rank_file_7<-Rank_file_7[which(Rank_file_7$phenotype != "NA"),]
  
  
  cat("Rank_file_7_NO_NA\n")
  cat(str(Rank_file_7))
  cat("\n")
  
  # Rank_file_3<-merge(Rank_file_3,
  #                  GENE.TABLE_equivalence,
  #                  by="ensembl_gene_id",
  #                  all.x=T)
  # 
  # cat("Rank_file_3\n")
  # cat(str(Rank_file_3))
  # cat("\n")
  
  Rank_file_7<-Rank_file_7[order(Rank_file_7$Rank, decreasing = F),]
  
  cat("Rank_file_7_ordered\n")
  cat(str(Rank_file_7))
  cat("\n")
  
  #quit(status = 1)
  
  
  
  
  
  Rank_file_7_MAX<-as.data.frame(setDT(Rank_file_7)[, .SD[which.min(Rank)], by=c("phenotype","block_no")])
  
  Rank_file_7_MAX$Allelic_Series_ID<-paste(Rank_file_7_MAX$phenotype,Rank_file_7_MAX$block_no, sep="__")
  
  
  cat("Rank_file_7_MAX_ordered\n")
  cat(str(Rank_file_7_MAX))
  cat("\n")
  
  ##### filter Rank_file_7_MAX ----
  
  Rank_file_7_MAX_filtered<-Rank_file_7_MAX[which(Rank_file_7_MAX$Allelic_Series_ID%in%Open_targets_GLOBAL_file_anchored_filtered$Allelic_Series_ID),]
  
  cat("Rank_file_7_MAX_filtered_1\n")
  cat(str(Rank_file_7_MAX_filtered))
  cat("\n")
  
  Rank_file_7_MAX_filtered<-Rank_file_7_MAX_filtered[-which(Rank_file_7_MAX_filtered$Allelic_Series_ID%in%dB_ALL_sel_check_PP_threshold_type_4$Allelic_Series_ID),]
  
  cat("Rank_file_7_MAX_filtered_2\n")
  cat(str(Rank_file_7_MAX_filtered))
  cat("\n")
  
  Rank_file_7_MAX_filtered<-Rank_file_7_MAX_filtered[-which(Rank_file_7_MAX_filtered$Allelic_Series_ID%in%dB_ALL_sel_check_PP_threshold_type_5$Allelic_Series_ID),]
  
  cat("Rank_file_7_MAX_filtered_3\n")
  cat(str(Rank_file_7_MAX_filtered))
  cat("\n")
  
  Rank_file_7_MAX_filtered<-Rank_file_7_MAX_filtered[-which(Rank_file_7_MAX_filtered$Allelic_Series_ID%in%dB_ALL_sel_check_PP_threshold_type_6$Allelic_Series_ID),]
  
  cat("Rank_file_7_MAX_filtered_4\n")
  cat(str(Rank_file_7_MAX_filtered))
  cat("\n")
  
  
  ###### Retrieve variants ----
  
  AS<-unique(Rank_file_7_MAX_filtered$Allelic_Series_ID[order(Rank_file_7_MAX_filtered$Rank, decreasing = F)])
  
  cat("AS\n")
  cat(str(AS))
  cat("\n")
  
  dB_ALL_sel_check<-dB_ALL[which(dB_ALL$Allelic_Series_ID%in%AS),]
  
  
  cat("dB_ALL_sel_check\n")
  cat(str(dB_ALL_sel_check))
  cat("\n")
  
  
  
  ##### Apply the probability threshold----
  
  dB_ALL_sel_check_PP_threshold_type_7<-unique(dB_ALL_sel_check[which(dB_ALL_sel_check$finemap_prob >= PP_threshold),])
  
  cat("dB_ALL_sel_check_PP_threshold_type_7\n")
  cat(str(dB_ALL_sel_check_PP_threshold_type_7))
  cat("\n")
  
  
  VARS_check<-unique(dB_ALL_sel_check_PP_threshold_type_7$VAR)
  
  cat("VARS_check\t")
  cat(sprintf(as.character(type_7)))
  cat("\t")
  cat(sprintf(as.character(length(VARS_check))))
  cat("\n")
  
  setwd(out)
  
  filename_1<-paste("Selection_",type_7,".tsv", sep='')
  
  write.table(dB_ALL_sel_check_PP_threshold_type_7,
              file=filename_1, sep="\t", quote=F, row.names = F)
  
  
  
  
  
  
}


collapsing_not_programmed = function(option_list)
{
  opt_in = option_list
  opt <<- option_list
  
  cat("All options:\n")
  printList(opt)
  
  #### READ and transform type_1 ----
  
  type_1 = opt$type_1
  
  cat("type_1\n")
  cat(sprintf(as.character(type_1)))
  cat("\n")
  
  filename_1<-paste("Selection_",type_1,".tsv", sep='')
  
  df_1<-read.table(file=filename_1, sep="\t", header=T, stringsAsFactors = F)
  
  cat("type_1\n")
  cat(str(type_1))
  cat("\n")
  
  #### READ and transform type_2 ----
  
  type_2 = opt$type_2
  
  cat("type_2\n")
  cat(sprintf(as.character(type_2)))
  cat("\n")
  
  filename_1<-paste("Selection_",type_2,".tsv", sep='')
  
  df_2<-read.table(file=filename_1, sep="\t", header=T, stringsAsFactors = F)
  
  cat("type_2\n")
  cat(str(type_2))
  cat("\n")
  
  #### READ and transform type_4 ----
  
  type_4 = opt$type_4
  
  cat("type_4\n")
  cat(sprintf(as.character(type_4)))
  cat("\n")
  
  filename_1<-paste("Selection_",type_4,".tsv", sep='')
  
  df_4<-read.table(file=filename_1, sep="\t", header=T, stringsAsFactors = F)
  
  cat("type_4\n")
  cat(str(type_4))
  cat("\n")
  
  #### READ and transform type_5 ----
  
  type_5 = opt$type_5
  
  cat("type_5\n")
  cat(sprintf(as.character(type_5)))
  cat("\n")
  
  filename_1<-paste("Selection_",type_5,".tsv", sep='')
  
  df_5<-read.table(file=filename_1, sep="\t", header=T, stringsAsFactors = F)
  
  cat("type_5\n")
  cat(str(type_5))
  cat("\n")
  
  #### READ and transform type_6 ----
  
  type_6 = opt$type_6
  
  cat("type_6\n")
  cat(sprintf(as.character(type_6)))
  cat("\n")
  
  filename_1<-paste("Selection_",type_6,".tsv", sep='')
  
  df_6<-read.table(file=filename_1, sep="\t", header=T, stringsAsFactors = F)
  
  cat("type_6\n")
  cat(str(type_6))
  cat("\n")
  
  #### READ and transform type_7 ----
  
  type_7 = opt$type_7
  
  cat("type_7\n")
  cat(sprintf(as.character(type_7)))
  cat("\n")
  
  filename_1<-paste("Selection_",type_7,".tsv", sep='')
  
  df_7<-read.table(file=filename_1, sep="\t", header=T, stringsAsFactors = F)
  
  cat("type_7\n")
  cat(str(type_7))
  cat("\n")
  
  #### READ and transform out ----
  
  out = opt$out
  
  cat("OUT_\n")
  cat(sprintf(as.character(out)))
  cat("\n")
  
  #### Rbind all dfs----
  
  NP<-rbind(df_1,df_2)
 
  cat("NP_1\n")
  cat(str(NP))
  cat("\n")
 
  
  NP<-rbind(NP,df_4)
  
  cat("NP_2\n")
  cat(str(NP))
  cat("\n")
  
  NP<-rbind(NP,df_5)
  
  cat("NP_4\n")
  cat(str(NP))
  cat("\n")
  
  NP<-rbind(NP,df_6)
  
  cat("NP_5\n")
  cat(str(NP))
  cat("\n")
  
  NP<-rbind(NP,df_7)
  
  cat("NP_6\n")
  cat(str(NP))
  cat("\n")
  
  NP_unique<-unique(NP)
  
  cat("NP_unique\n")
  cat(str(NP_unique))
  cat("\n")
  
  VARS_NP_unique<-unique(NP_unique$VAR)
  
  cat("VARS_NP_unique\n")
  cat(str(VARS_NP_unique))
  cat("\n")
  
  #### SAVE ----
  
  setwd(out)
  
  filename_1<-paste("Selection_","NON_PROGRAMMED",".tsv", sep='')
  
  write.table(NP_unique,
              file=filename_1, sep="\t", quote=F, row.names = F)
  
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
  traceback()
  options(show.error.locations = TRUE)
  
  option_list <- list(
    make_option(c("--dB_ALL"), type="character", default=NULL,
                metavar="FILE.txt",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Rank_file"), type="character", default=NULL,
                metavar="FILE.txt",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Rank_file_3"), type="character", default=NULL,
                metavar="FILE.txt",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Rank_file_4"), type="character", default=NULL,
                metavar="FILE.txt",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Rank_file_5"), type="character", default=NULL,
                metavar="FILE.txt",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Rank_file_6"), type="character", default=NULL,
                metavar="FILE.txt",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Rank_file_7"), type="character", default=NULL,
                metavar="FILE.txt",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--PP_threshold"), type="numeric", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Original_selection_Grouped_variants"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--categories_of_RRVV"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--type_1"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--type_2"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--type_3"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--cut_3"), type="numeric", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--type_4"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--type_5"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--type_6"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--type_7"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Open_targets_GLOBAL_file"), type="character", default=NULL, 
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
  
  
  Prioritizing_subsets(opt)
  Prioritizing_subsets_2(opt)
  Prioritizing_subsets_3(opt)
  Prioritizing_subsets_4(opt)
  collapsing_not_programmed(opt)
  
}




###########################################################################

system.time( main() )


