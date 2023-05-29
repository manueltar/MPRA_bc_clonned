
suppressMessages(library("Sushi"))


#suppressMessages(library("dplyr", lib.loc="/nfs/users/nfs_m/mt19/sOFTWARE/R_libs/"))
suppressMessages(library("data.table", lib.loc="/nfs/users/nfs_m/mt19/sOFTWARE/R_libs/"))
suppressMessages(library("ggplot2", lib.loc = "/nfs/users/nfs_m/mt19/sOFTWARE/R_libs"))
suppressMessages(library("optparse", lib.loc = "/nfs/users/nfs_m/mt19/sOFTWARE/R_libs"))
#suppressMessages(library("plyr", lib.loc="/nfs/users/nfs_m/mt19/sOFTWARE/R_libs/"))
suppressMessages(library("splitstackshape", lib.loc="/nfs/users/nfs_m/mt19/sOFTWARE/R_libs/"))
suppressMessages(library("tidyverse"))
suppressMessages(library("ggnewscale", lib.loc="/nfs/users/nfs_m/mt19/sOFTWARE/R_libs/"))
suppressMessages(library("scales", lib.loc="/nfs/users/nfs_m/mt19/sOFTWARE/R_libs/"))
suppressMessages(library("udunits2", lib.loc = "/nfs/users/nfs_m/mt19/sOFTWARE/R_libs/"))
suppressMessages(library("ggforce", lib.loc = "/nfs/users/nfs_m/mt19/sOFTWARE/R_libs/"))
suppressMessages(library("viridis", lib.loc="/nfs/users/nfs_m/mt19/sOFTWARE/R_libs/"))
suppressMessages(library("gridBase", lib.loc="/nfs/users/nfs_m/mt19/sOFTWARE/R_libs/"))
suppressMessages(library("gridExtra", lib.loc="/nfs/users/nfs_m/mt19/sOFTWARE/R_libs/"))
suppressMessages(library("grid", lib.loc="/nfs/users/nfs_m/mt19/sOFTWARE/R_libs/"))


opt = NULL

Analysis_fmt3 = function(option_list)
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
 
  #### Read master tiles ----
  
  setwd(out)
  
  filename<-paste(type,"_MPRA_Rosetta.tsv", sep='')
  
  Rosetta<-read.table(file=filename, sep="\t", header=T, stringsAsFactors = F)
  
  cat("Rosetta\n")
  cat(str(Rosetta))
  cat("\n")
  
  indx.int<-c(which(colnames(Rosetta) == "seq_name"),which(colnames(Rosetta) == "carried_variants"))
  
  Rosetta_subset<-unique(Rosetta[,indx.int])
  
  cat("Rosetta_subset\n")
  cat(str(Rosetta_subset))
  cat("\n")
  
  Rosetta_subset<-as.data.frame(cSplit(Rosetta_subset, splitCols = "carried_variants", 
                                                     sep = "|", direction = "long", drop = TRUE),stringsAsFactors=F)
  
  Rosetta_subset$carried_variants<-as.character(Rosetta_subset$carried_variants)
  
  cat("Rosetta_subset_2\n")
  str(Rosetta_subset)
  cat("\n")
  
  #### Read fmt3 result ----
  
  setwd(out)
  
  filename<-'global_file_fmt3_parsed.tsv'
  
  fmt3<-read.table(file=filename, sep="\t", header=F, stringsAsFactors = F)
 
  colnames(fmt3) <-c("chr","pos","ID","ref","alt","QUAL","FILTER","seq_name")
  
  fmt3$chr<-gsub("^Chr","chr",fmt3$chr)
  
  cat("fmt3\n")
  cat(str(fmt3))
  cat("\n")
  
  #### Subset which fmt3 in Rosetta (should be only the ALT)----
 
  fmt3_check<-fmt3[which(fmt3$seq_name%in%Rosetta$seq_name),]
  
  cat("fmt3_check\n")
  cat(str(fmt3_check))
  cat("\n")
  
  fmt3_lost<-fmt3[-which(fmt3$seq_name%in%Rosetta$seq_name),]
  
  cat("fmt3_lost\n")
  cat(str(fmt3_lost))
  cat("\n")
  
  ####----
  
  indx.int<-c(which(colnames(fmt3_check) == "seq_name"),which(colnames(fmt3_check) == "chr"),
              which(colnames(fmt3_check) == "pos"),which(colnames(fmt3_check) == "ref"),which(colnames(fmt3_check) == "alt"))
  
  
  fmt3_check_subset<-unique(fmt3_check[,indx.int])
  
  fmt3_check_subset$carried_variants<-paste(fmt3_check_subset$chr,fmt3_check_subset$pos,fmt3_check_subset$ref,fmt3_check_subset$alt,sep="_")
  
  cat("fmt3_check_subset\n")
  cat(str(fmt3_check_subset))
  cat("\n")
  
  #### MERGE ---
  
  check_MERGE<-merge(fmt3_check_subset,
                     Rosetta_subset,
                     by=c("seq_name","carried_variants"))
  
  cat("check_MERGE\n")
  cat(str(check_MERGE))
  cat("\n")
  
  #### fmt3 Classif ----
  
  Rosetta$fmt3_classif<-"NA"
  
  Rosetta$fmt3_classif[which(Rosetta$seq_name%in%fmt3_check_subset$seq_name)]<-"carried_variants_recovered_totally_or_partially"
  Rosetta$fmt3_classif[-which(Rosetta$seq_name%in%fmt3_check_subset$seq_name)]<-"carried_variants_NOT_recovered"
  Rosetta$fmt3_classif[which(Rosetta$carried_variants == "REF")]<-"REF"
  
  
  #### SAVE MASTER TABLE ----
  
 #  quit(status = 1)
  
  setwd(out)

  filename<-paste(type,"_fmt3_lost",".tsv", sep='')
  write.table(fmt3_lost, file=filename, sep="\t", row.names = F, quote=F)
  # 
  # filename<-paste(type,"_LONG_matrix_seq_names",".tsv", sep='')
  # write.table(Long_df_ordered, file=filename, sep="\t", row.names = F, quote=F)
  
  filename<-paste(type,"_MPRA_Rosetta_Plus_fmt3.tsv", sep='')
  
  write.table(Rosetta, file=filename, sep="\t", row.names = F, quote=F)
  
}

Analysis_fmt7 = function(option_list)
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
  
  #### Read master tiles ----
  
  setwd(out)
  
  filename<-paste(type,"_MPRA_Rosetta_Plus_fmt3.tsv", sep='')
  
  Rosetta<-read.table(file=filename, sep="\t", header=T, stringsAsFactors = F)
  
  cat("Rosetta\n")
  cat(str(Rosetta))
  cat("\n")
  
  # indx.int<-c(which(colnames(Rosetta) == "seq_name"),which(colnames(Rosetta) == "carried_variants"))
  # 
  # Rosetta_subset<-unique(Rosetta[,indx.int])
  # 
  # cat("Rosetta_subset\n")
  # cat(str(Rosetta_subset))
  # cat("\n")
  # 
  # Rosetta_subset<-as.data.frame(cSplit(Rosetta_subset, splitCols = "carried_variants", 
  #                                      sep = "|", direction = "long", drop = TRUE),stringsAsFactors=F)
  # 
  # Rosetta_subset$carried_variants<-as.character(Rosetta_subset$carried_variants)
  # 
  # cat("Rosetta_subset_2\n")
  # str(Rosetta_subset)
  # cat("\n")
  
  #### Read fmt7 result ----
  
  setwd(out)
  
  filename<-'global_file_fmt7.txt'
  
  fmt7<-as.data.frame(fread(file=filename, sep='\t', header=F, fill=TRUE)[!grepl( "^#", V1 )], stringsAsFactors = F)
  
  
  
  cat("fmt7\n")
  cat(str(fmt7))
  cat("\n")
  
  # indx.alignments<-grep ("^Element_[0-9]+__TILE",fmt7$V1)
  # 
  # fmt7_subset<-as.data.frame(unique(fmt7[indx.alignments,]), stringsAsFactors = F)
  # colnames(fmt7_subset)<-"V1"
  # rm(fmt7)
  # 
  # cat("fmt7_subset_0\n")
  # cat(str(fmt7_subset))
  # cat("\n")
  
  fmt7_subset<-as.data.frame(cSplit(fmt7, splitCols = "V1",
                                       sep = '\t', direction = "wide", drop = TRUE),stringsAsFactors=F)

  # cat("fmt7_subset_2\n")
  # cat(str(fmt7_subset))
  # cat("\n")
  
  colnames(fmt7_subset)<-c('seq_name','subject acc.ver','% identity','alignment length','mismatches','gap opens','q. start','q. end','s. start','s. end','evalue','bit score')
  fmt7_subset[,1]<-as.character(fmt7_subset[,1])
  fmt7_subset[,2]<-as.character(fmt7_subset[,2])
  
  
  # 
  cat("fmt7_subset_2\n")
  cat(str(fmt7_subset))
  cat("\n")
  
  freq<-as.data.frame(table(fmt7_subset$seq_name),stringsAsFactors=F)
  
  colnames(freq)<-c('seq_name','Genome_Alignments')
  
  cat("freq_2\n")
  cat(str(freq))
  cat("\n")
  
  freq$REAL_TILE<-gsub("\\|.+$","",freq$seq_name)
  
  cat("freq_3\n")
  cat(str(freq))
  cat("\n")
  
  freq_subset<-unique(freq[,-which(colnames(freq) == "seq_name")])
  
  cat("freq_subset_3\n")
  cat(str(freq_subset))
  cat("\n")
  
  
  Rosetta<-merge(Rosetta,
                 freq_subset,
                 by='REAL_TILE',
                 all.x=T)
  
  Rosetta$Genome_Alignments[is.na(Rosetta$Genome_Alignments)]<-"NOT_ALIGNED"
  
  cat("Rosetta\n")
  cat(str(Rosetta))
  cat("\n")
  
  #### SAVE MASTER TABLE ----
  
  #  quit(status = 1)
  
  setwd(out)
  
 
  filename<-paste(type,"_MPRA_Rosetta_Plus_fmt3_Plus_fmt7.tsv", sep='')
  
  write.table(Rosetta, file=filename, sep="\t", row.names = F, quote=F)
  
  filename<-paste(type,"_Multialignment_fmt7.tsv", sep='')
  
  
  write.table(freq_subset, file=filename, sep="\t", row.names = F, quote=F)
  
  
 
}

Per_AS_report = function(option_list)
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
  
  #### Selection_file ----
  
  Selection_file<-as.data.frame(fread(file = opt$Selection_file,
                                      sep="\t",
                                      header=T), stringsAsFactors=F)
  
  Selection_file$Allelic_Series_ID<-paste(Selection_file$phenotype,Selection_file$block_no, sep="__")
  
  
  cat("Selection_file\n")
  cat(str(Selection_file))
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

  #### Read master tiles ----
  
  setwd(out)
  
  filename<-paste(type,"_MPRA_Rosetta_Plus_fmt3_Plus_fmt7.tsv", sep='')
  
  Rosetta<-read.table(file=filename, sep="\t", header=T, stringsAsFactors = F)
  
  cat("Rosetta\n")
  cat(str(Rosetta))
  cat("\n")
  
  indx.dep<-c(which(colnames(Rosetta) == "barcode"),
              which(colnames(Rosetta) == "proposed"),
              which(colnames(Rosetta) == "seq_name_bc"))
  
  if(length(indx.dep) >0){
    
    Rosetta_subset<-unique(Rosetta[,-indx.dep])
    
  }#length(indx.dep) >0
  else{
    
    Rosetta_subset<-Rosetta
    
    
  }#length(indx.dep) <0
  
  cat("Rosetta_subset_0\n")
  cat(str(Rosetta_subset))
  cat("\n")
  
  Rosetta_subset<-Rosetta_subset[which(Rosetta_subset$carried_variants != "REF"),]
  
  cat("Rosetta_subset_1\n")
  cat(str(Rosetta_subset))
  cat("\n")
  
  #### Important points x and y ----
  
  ### GC_content
  
  pointsGC<-round(as.numeric(summary(Rosetta_subset$GC_content)),1)
  GC_min<-pointsGC[1]-1
  GC_max<-pointsGC[length(pointsGC)]+1
  
  pointsGC<-c(GC_min,round(as.numeric(summary(Rosetta_subset$GC_content)),1),GC_max)
  
  breaks.GC<-sort(unique(c(pointsGC)))
  labels.GC<-as.character(breaks.GC)
  
  cat("labels.GC\n")
  cat(sprintf(as.character(labels.GC)))
  cat("\n")
  
  ### Genome_Alignments
  
  levels_Genome_Alignments<-levels(as.factor(Rosetta_subset$Genome_Alignments))
  
  cat("levels_Genome_Alignments\n")
  cat(sprintf(as.character(names(summary(as.factor(Rosetta_subset$Genome_Alignments))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(Rosetta_subset$Genome_Alignments)))))
  cat("\n")
  
  
  Rosetta_subset$mA<-Rosetta_subset$Genome_Alignments
  Rosetta_subset$mA[Rosetta_subset$Genome_Alignments == "NOT_ALIGNED"]<-0
  
  Rosetta_subset$mA<-as.numeric(Rosetta_subset$mA)
  
  pointsmA<-round(as.numeric(summary(Rosetta_subset$mA)),0)
  mA_min<-0
  mA_max<-pointsmA[length(pointsmA)]+1

  pointsmA<-c(mA_min,round(as.numeric(summary(Rosetta_subset$mA)),0),mA_max)

  cat("pointsmA\n")
  cat(sprintf(as.character(pointsmA)))
  cat("\n")

  breaks.mA<-sort(unique(pointsmA))
  log_breaks.mA<-log(breaks.mA+1)
  labels.mA<-as.character(breaks.mA)
  
  cat("log_breaks.mA\n")
  cat(sprintf(as.character(log_breaks.mA)))
  cat("\n")

  cat("labels.mA\n")
  cat(sprintf(as.character(labels.mA)))
  cat("\n")
  
  ### levels_FLAG
  
  Rosetta_subset$Flag_factor<-gsub("NO_FLAG","",Rosetta_subset$Flag_intra_candidate_pattern)
  Rosetta_subset$Flag_factor<-gsub("\\|","",Rosetta_subset$Flag_factor)
  #Rosetta_subset$Flag_factor<-gsub("/\/","",Rosetta_subset$Flag_factor)
  
  Rosetta_subset$Flag_factor[Rosetta_subset$Flag_factor == ""]<-"NO_PROBLEM"
  
  levels_FLAG<-levels(as.factor(Rosetta_subset$Flag_factor))
  
  
  Rosetta_subset$Flag_factor<-factor(Rosetta_subset$Flag_factor,
                                     levels=levels_FLAG,
                                     ordered=T)
  
  cat("levels_FLAG\n")
  cat(sprintf(as.character(names(summary(Rosetta_subset$Flag_factor)))))
  cat("\n")
  cat(sprintf(as.character(summary(Rosetta_subset$Flag_factor))))
  cat("\n")
    
  ### levels_fmt3
  
  levels_fmt3<-levels(as.factor(Rosetta_subset$fmt3_classif))
  
  Rosetta_subset$fmt3_factor<-factor(Rosetta_subset$fmt3_classif,
                                     levels=levels_fmt3,
                                     ordered=T)
  
  cat("levels_fmt3\n")
  cat(sprintf(as.character(names(summary(Rosetta_subset$fmt3_factor)))))
  cat("\n")
  cat(sprintf(as.character(summary(Rosetta_subset$fmt3_factor))))
  cat("\n")
  

  ### levels_factor4
  
  levels_factor4<-levels(as.factor(Rosetta_subset$factor4))
  
  Rosetta_subset$factor4<-factor(Rosetta_subset$factor4,
                                     levels=c("High_AT","Medium","High_GC"),
                                     ordered=T)
  
  cat("levels_factor4\n")
  cat(sprintf(as.character(names(summary(Rosetta_subset$factor4)))))
  cat("\n")
  cat(sprintf(as.character(summary(Rosetta_subset$factor4))))
  cat("\n")
  
  ### levels_Label
  
  levels_Label<-levels(as.factor(Rosetta_subset$Label))
  
  Rosetta_subset$Label<-factor(Rosetta_subset$Label,
                               levels=c("Negative_Control_Genomic_Regions","UNDETERMINED_CTRL","ASE_CTRL","ASSAYED_VARIANT"),
                                 ordered=T)
  
  cat("levels_Label\n")
  cat(sprintf(as.character(names(summary(Rosetta_subset$Label)))))
  cat("\n")
  cat(sprintf(as.character(summary(Rosetta_subset$Label))))
  cat("\n")
  
  #### Rank
  
  breaks.Rank<-seq(0,1,by=0.2)
  labels.Rank<-as.character(breaks.Rank)
  
  #### finemap_beta
  max_value<-max(Selection_file$finemap_beta)
  min_value<-min(Selection_file$finemap_beta)
  
  breaks.finemap_beta<-round(sort(c(min_value-0.01,seq(min_value,max_value,by=0.5),0,max_value+0.01)),2)
  labels.finemap_beta<-as.character(breaks.finemap_beta)
  
  max_finemap_z<-110
  min_finemap_z<--125
  
  breaks.finemap_z<-round(sort(c(min_finemap_z,seq(min_finemap_z,max_finemap_z,by=50),0,max_finemap_z)),2)
  labels.finemap_z<-as.character(breaks.finemap_z)
  
  breaks.maf=seq(0,0.5,by=0.05)
  labels.maf=as.character(breaks.maf)
  
  #quit(status =1)
  ####
  
  # filename<-paste(type,"_Multialignment_fmt7.tsv", sep='')
  # 
  # multi_alignment<-read.table(file=filename, sep="\t", header=T, stringsAsFactors = F)
  # 
  # cat("multi_alignment\n")
  # cat(str(multi_alignment))
  # cat("\n")
  
  filename<-paste(type,"_fmt3_lost",".tsv", sep='')
  
  fmt3_lost<-read.table(file=filename, sep="\t", header=T, stringsAsFactors = F)
  
  cat("fmt3_lost\n")
  cat(str(fmt3_lost))
  cat("\n")
  
 
  
  ###### Conditions: ANCHORED IN AN INTERESTING DISEASE-----
  
  blood_diseases<-c("anemia","anemia (disease)","hemolytic anemia","megaloblastic anemia (disease)","pernicious anemia","Elevated erythrocyte sedimentation rate","Hemoglobinopathy",
                    "Rare genetic coagulation disorder","blood coagulation disease","Rare hemorrhagic disorder due to a  constitutional platelet anomaly","response to anticoagulant","response to platelet aggregation inhibitor",
                    "thrombocytopenia","Thrombotic disease","Thrombophlebitis","venous thromboembolism","deep vein thrombosis","primary thrombocytopenia",
                    "leukocyte disease","Leukocytosis","Leukopenia","Eosinophilia",
                    "hematologic disease","Abnormality of blood and blood-forming tissues","Pancytopenia")
  
  inflammatory_diseases<-c("allergy","asthma","Bronchospasm","purpura (disease)","type II diabetes mellitus","type I diabetes mellitus","inflammatory bowel disease",
                           "Crohn's disease","rheumatoid arthritis","systemic lupus erythematosus","lupus erythematosus","multiple sclerosis","Sarcoidosis",
                           "psoriasis vulgaris","psoriatic arthritis",
                           "Inflammation","inflammation","systemic inflammatory response syndrome")
  
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
  
  
  ##### Attribute AS ----
  
  VARS<-unique(Rosetta_subset$VAR)
  
  cat("VARS\n")
  cat(str(VARS))
  cat("\n")
  
  ATTRIBUTED_VARS<-VARS[which(VARS%in%Selection_file$VAR)]
 
  cat("ATTRIBUTED_VARS\n")
  cat(str(ATTRIBUTED_VARS))
  cat("\n")
  
  NON_ATTRIBUTED_VARS<-VARS[-which(VARS%in%Selection_file$VAR)]
  
  cat("NON_ATTRIBUTED_VARS\n")
  cat(str(NON_ATTRIBUTED_VARS))
  cat("\n")
  
  #Selection_file<-VARS[-which(VARS%in%Selection_file$VAR)]
  
  cat("Selection_file\n")
  cat(str(Selection_file))
  cat("\n")
  
  header_line<-c("VAR","variant_id","rs","chr","pos37","ref","alt","maf_origin","maf_dbsnp","vep_worst_cons","block_no","phenotype","ncondind",
                  "credset_size","credset_prob","finemap_beta","finemap_se","finemap_z","finemap_prob","finemap_log10bf","info","is_cond_ind","chromstates","Allelic_Series_ID")
  
  
  NON_ATTRIBUTED_VARS_df<- data.frame(matrix(vector(), length(NON_ATTRIBUTED_VARS), length(header_line),
                             dimnames=list(c(),
                                           header_line)),
                      stringsAsFactors=F)
  
  cat("NON_ATTRIBUTED_VARS_df_0\n")
  cat(str(NON_ATTRIBUTED_VARS_df))
  cat("\n")
  
  NON_ATTRIBUTED_VARS_df$VAR<-NON_ATTRIBUTED_VARS
  
  cat("NON_ATTRIBUTED_VARS_df_1\n")
  cat(str(NON_ATTRIBUTED_VARS_df))
  cat("\n")
  
  NON_ATTRIBUTED_VARS_df$chr<-gsub("_.+$","",NON_ATTRIBUTED_VARS_df$VAR)
  NON_ATTRIBUTED_VARS_df$pos37<-gsub("^[^_]+_.","",NON_ATTRIBUTED_VARS_df$VAR)
  NON_ATTRIBUTED_VARS_df$pos37<-gsub("_.$","",NON_ATTRIBUTED_VARS_df$pos37)
  NON_ATTRIBUTED_VARS_df$ref<-gsub("^[^_]+_[^_]+_","",NON_ATTRIBUTED_VARS_df$VAR)
  NON_ATTRIBUTED_VARS_df$ref<-gsub("_.$","",NON_ATTRIBUTED_VARS_df$ref)
  NON_ATTRIBUTED_VARS_df$alt<-gsub("^.+_","",NON_ATTRIBUTED_VARS_df$VAR)
  
  cat("NON_ATTRIBUTED_VARS_df_2\n")
  cat(str(NON_ATTRIBUTED_VARS_df))
  cat("\n")
  
  NON_ATTRIBUTED_VARS_df$Allelic_Series_ID<-"NON_ATTRIBUTED"
  
  cat("NON_ATTRIBUTED_VARS_df_3\n")
  cat(str(NON_ATTRIBUTED_VARS_df))
  cat("\n")
  
  
  composite_df<-rbind(Selection_file,NON_ATTRIBUTED_VARS_df)
  
  cat("composite_df\n")
  cat(str(composite_df))
  cat("\n")
  
  #quit(status=1)
  ##### LOOP AS -----
  
  AS<-unique(composite_df$Allelic_Series_ID)

  cat("AS\n")
  cat(str(AS))
  cat("\n")
  
  path_reports<-paste(out,'AS_pdf_reports', sep='')
  
 
  
  list_charac_tiles<-list()
  list_hmap_Sel<-list()
  list_blood_disease<-list()
  list_inflammatory_diseases<-list()
  list_immunodeficiencies<-list()
  list_leukemias<-list()
  list_miscellaneous<-list()
  
  list_hmap_Rosetta<-list()
  
  
  if (file.exists(path_reports)){
    setwd(path_reports)
  } else {
    dir.create(file.path(path_reports))
    setwd(path_reports)
    
    cat("paths_definitive\n")
    cat(sprintf(as.character(path_reports)))
    cat("\n")
    
  }
  
  for(i in 1:length(AS)){
    
    AS_sel<-AS[i]
    
    cat("--------------------------------------------------------------------->AS_sel:\t")
    cat(sprintf(as.character(i)))
    cat("\t")
    cat(sprintf(as.character(AS_sel)))
    cat("\n")
    
    composite_df_sel<-composite_df[which(composite_df$Allelic_Series_ID == AS_sel),]
    
    # cat("composite_df_sel\n")
    # cat(str(composite_df_sel))
    # cat("\n")
    
    #### MPRA VARIANTS -----
    
    Rosetta_subset_sel<-Rosetta_subset[which(Rosetta_subset$VAR%in%composite_df_sel$VAR),]
    
    # cat("Rosetta_subset_sel\n")
    # cat(str(Rosetta_subset_sel))
    # cat("\n")
   
    
    ##### Per element plots ----
    
    KEYS<-unique(Rosetta_subset_sel$KEY)
    
    # cat("KEYS\n")
    # cat(str(KEYS))
    # cat("\n")
    
    pdfname<-paste("Graphs_",AS_sel,".pdf", sep='')
    makepdf = TRUE
    
    if (makepdf == TRUE)
    {
      pdf ( pdfname , height=10, width=12)
    }
    
    #### like heatmap with discrete info ----
    
    list_hmap_Rosetta[[i]]<-unique(Rosetta_subset_sel[,-c(which(colnames(Rosetta_subset_sel) == "Flag_intra_candidate_pattern"),
                                                          which(colnames(Rosetta_subset_sel) == "fmt3_classif"),
                                                          which(colnames(Rosetta_subset_sel) == "sequence"),
                                                          which(colnames(Rosetta_subset_sel) == "UP_FWD_Primer"),
                                                          which(colnames(Rosetta_subset_sel) == "UP_REV_Primer"),
                                                          which(colnames(Rosetta_subset_sel) == "BamHI"),
                                                          which(colnames(Rosetta_subset_sel) == "KpnI"),
                                                          which(colnames(Rosetta_subset_sel) == "GC_content_intervals"),
                                                          which(colnames(Rosetta_subset_sel) == "Genome_Alignments"),
                                                          which(colnames(Rosetta_subset_sel) == "Flag_factor"),
                                                          which(colnames(Rosetta_subset_sel) == "fmt3_factor"))])
    
    # cat("list_hmap_Rosetta\n")
    # cat(str(list_hmap_Rosetta))
    # cat("\n")
    # 
    # quit(status = 1)
    
    Rosetta_subset_sel_REP<-unique(Rosetta_subset_sel[,-c(which(colnames(Rosetta_subset_sel) == "REAL_TILE"),
                                                          which(colnames(Rosetta_subset_sel) == "VAR"),
                                                          which(colnames(Rosetta_subset_sel) == "KEY"),
                                                          which(colnames(Rosetta_subset_sel) == "sequence"),
                                                          which(colnames(Rosetta_subset_sel) == "Flag_intra_candidate_pattern"),
                                                          which(colnames(Rosetta_subset_sel) == "GC_content"),
                                                          which(colnames(Rosetta_subset_sel) == "GC_content_intervals"),
                                                          which(colnames(Rosetta_subset_sel) == "UP_FWD_Primer"),
                                                          which(colnames(Rosetta_subset_sel) == "UP_REV_Primer"),
                                                          which(colnames(Rosetta_subset_sel) == "BamHI"),
                                                          which(colnames(Rosetta_subset_sel) == "KpnI"),
                                                          which(colnames(Rosetta_subset_sel) == "fmt3_classif"),
                                                          which(colnames(Rosetta_subset_sel) == "Genome_Alignments"),
                                                          which(colnames(Rosetta_subset_sel) == "mA"))])
                                                          
                                        
    
    Rosetta_subset_sel_REP$Tile<-factor(Rosetta_subset_sel_REP$Tile,
                                               levels=c("TILE_1","TILE_2","TILE_3","TILE_4","TILE_5"),
                                    ordered=T)
    Rosetta_subset_sel_REP$carried_variants<-factor(Rosetta_subset_sel_REP$carried_variants,
                                    levels=levels(as.factor(Rosetta_subset_sel_REP$carried_variants)),
                                    ordered=T)
    
    # cat("Rosetta_subset_sel_REP\n")
    # cat(str(Rosetta_subset_sel_REP))
    # cat("\n")
    
    Rosetta_subset_sel_REP<-Rosetta_subset_sel_REP[order(Rosetta_subset_sel_REP$start, Rosetta_subset_sel_REP$Tile,Rosetta_subset_sel_REP$carried_variants, decreasing = F),]
    
    # cat("Rosetta_subset_sel_REP_ordered\n")
    # cat(str(Rosetta_subset_sel_REP))
    # cat("\n")
    
    
    #Rosetta_subset_sel_REP.df<-melt(setDT(Rosetta_subset_sel_REP), id.vars =c("chr","start","end","Tile","carried_variants","seq_name","Label","factor4","Flag_factor","fmt3_factor"))
    
    Rosetta_subset_sel_REP.df<-melt(setDT(Rosetta_subset_sel_REP), id.vars =c("chr","start","end","Tile","carried_variants","seq_name"))
    
    
    # cat("Rosetta_subset_sel_REP.df\n")
    # cat(str(Rosetta_subset_sel_REP.df))
    # cat("\n")
    
    # setwd(out)
    # 
    # write.table(Rosetta_subset_sel_REP.df, file="test.tsv", quote=F,row.names = F, sep="\t")
    
    
    g<-ggplot() +
      scale_y_discrete(name=NULL)+
      scale_x_discrete(name=NULL)+
      geom_point(data=Rosetta_subset_sel_REP.df[which(Rosetta_subset_sel_REP.df$variable == "Label"),],
                 aes(seq_name, variable,color=value), size=3)+
      geom_point(data=Rosetta_subset_sel_REP.df[which(Rosetta_subset_sel_REP.df$variable == "factor4"),],
                 aes(seq_name, variable,color=value), size=3)+
      geom_point(data=Rosetta_subset_sel_REP.df[which(Rosetta_subset_sel_REP.df$variable== "Flag_factor"),],
                 aes(seq_name, variable,color=value), size=3)+
      geom_point(data=Rosetta_subset_sel_REP.df[which(Rosetta_subset_sel_REP.df$variable== "fmt3_factor"),],
                 aes(seq_name, variable,color=value), size=3)+
      scale_shape(name = "Flags", solid = TRUE)+
      theme(axis.text.x  = element_text(angle=45, vjust = 1, hjust = 1))+
      theme(axis.text.y  = element_text(angle=45, vjust = 1, hjust = 1))+
      theme(legend.position="bottom")+
      guides(fill=guide_legend(nrow=2,byrow=TRUE))+
      geom_vline(xintercept = c(5,10,15,20,25,30,40,50,60,70,80,90,100,110,120), linetype="dotted", color="black")+
      ggtitle("Tiles characteristics")
    
    print(g)
    
    list_charac_tiles[[i]]<-Rosetta_subset_sel_REP.df
    
    #### like heatmap with finemapo info ----
    
    composite_df_sel_REP<-unique(composite_df_sel[,-c(which(colnames(composite_df_sel) == "variant_id"),
                                                          which(colnames(composite_df_sel) == "maf_dbsnp"),
                                                          which(colnames(composite_df_sel) == "vep_worst_cons"),
                                                          which(colnames(composite_df_sel) == "block_no"),
                                                          which(colnames(composite_df_sel) == "phenotype"),
                                                          which(colnames(composite_df_sel) == "ncondind"),
                                                          which(colnames(composite_df_sel) == "credset_size"),
                                                          which(colnames(composite_df_sel) == "credset_prob"),
                                                          which(colnames(composite_df_sel) == "finemap_log10bf"),
                                                          which(colnames(composite_df_sel) == "chromstates"),
                                                          which(colnames(composite_df_sel) == "finemap_se"),
                                                          which(colnames(composite_df_sel) == "KpnI"))])
    
    
   
    composite_df_sel_REP.df<-melt(setDT(composite_df_sel_REP), id.vars =c("VAR","rs","chr","pos37","ref","alt","Allelic_Series_ID"))
    
    
    # cat("composite_df_sel_REP.df\n")
    # cat(str(composite_df_sel_REP.df))
    # cat("\n")
    
    # setwd(out)
    # 
    # write.table(composite_df_sel_REP.df, file="test.tsv", quote=F,row.names = F, sep="\t")
    
    Ranks_0_to_1<-c("info","is_cond_ind")
    Ranks_0_to_0.5<-"maf_origin"
    Ranks_0_to_1_PP<-"finemap_prob"
    Ranks_finemap_beta<-"finemap_beta"
    Ranks_finemap_z<-"finemap_z"
    
    
    
    g<-ggplot() +
      geom_tile(data=composite_df_sel_REP.df[which(composite_df_sel_REP.df$variable%in%Ranks_0_to_1),],
                aes(reorder(VAR,pos37), variable,fill = value))+
      scale_fill_gradient(name = paste("0_to_1", sep="\n"),
                          low = "white", high = "red",
                          na.value = NA,
                          breaks=breaks.Rank,
                          labels=labels.Rank,
                          limits=c(breaks.Rank[1],
                                   breaks.Rank[length(breaks.Rank)]))+
      new_scale("fill")+
      geom_tile(data=composite_df_sel_REP.df[which(composite_df_sel_REP.df$variable%in%Ranks_0_to_0.5),],
                aes(reorder(VAR,pos37), variable,fill = value))+
      scale_fill_gradient(name = paste("maf", sep="\n"),
                          low = "blue", high = "white",
                          na.value = NA,
                          breaks=seq(0,0.5,by=0.05),
                          labels=as.character(seq(0,0.5,by=0.05)),
                          limits=c(0,0.5))+
      new_scale("fill")+
      geom_tile(data=composite_df_sel_REP.df[which(composite_df_sel_REP.df$variable%in%Ranks_0_to_1_PP),],
                aes(reorder(VAR,pos37), variable,fill = value))+
      scale_fill_gradient(name = paste("PP", sep="\n"),
                          low = "white", high = "green",
                          na.value = NA,
                          breaks=breaks.Rank,
                          labels=labels.Rank,
                          limits=c(breaks.Rank[1],
                                   breaks.Rank[length(breaks.Rank)]))+
      new_scale("fill")+
      geom_tile(data=composite_df_sel_REP.df[which(composite_df_sel_REP.df$variable%in%Ranks_finemap_beta),],
                aes(reorder(VAR,pos37), variable,fill = value))+
      scale_fill_gradient(name = paste("Beta", sep="\n"),
                          low = "red", high = "green",
                          na.value = NA,
                          breaks=breaks.finemap_beta,
                          labels=labels.finemap_beta,
                          limits=c(breaks.finemap_beta[1],
                                   breaks.finemap_beta[length(breaks.finemap_beta)]))+
      new_scale("fill")+
      geom_tile(data=composite_df_sel_REP.df[which(composite_df_sel_REP.df$variable%in%Ranks_finemap_z),],
                aes(reorder(VAR,pos37), variable,fill = value))+
      scale_fill_gradient(name = paste("Z-score", sep="\n"),
                          low = "red", high = "green",
                          na.value = NA,
                          breaks=breaks.finemap_z,
                          labels=labels.finemap_z,
                          limits=c(breaks.finemap_z[1],
                                   breaks.finemap_z[length(breaks.finemap_z)]))+
      scale_y_discrete(name=NULL, drop=F)+
      scale_x_discrete(name=NULL)+
      theme(axis.text.x  = element_text(angle=45, vjust = 1, hjust = 1))+
      theme(axis.text.y  = element_text(angle=45, vjust = 1, hjust = 1))+
      theme(legend.position="right")
    
    print(g)
    
    list_hmap_Sel[[i]]<-composite_df_sel_REP.df
    
    # if(AS_sel == "NON_ATTRIBUTED")
    # {
    #   
    #   quit(status = 1)
    # }
    
    
    
    #### GWAS disease anchored ----
    
    GWAS_anchored<-Open_targets_GLOBAL_file_anchored_filtered[which(Open_targets_GLOBAL_file_anchored_filtered$Allelic_Series_ID == AS_sel),]
    
    
    if(dim(GWAS_anchored)[1] > 0)
    {
      
      # cat("GWAS_anchored\n")
      # cat(str(GWAS_anchored))
      # cat("\n")
      
      blood_disease.df<-GWAS_anchored[which(GWAS_anchored$Disease_or_phenotype%in%blood_diseases),]
      
      blood_disease.df$Disease_or_phenotype<-factor(blood_disease.df$Disease_or_phenotype,
                                                    levels=blood_diseases,
                                                    ordered=T)
      
      # cat("blood_disease.df\n")
      # cat(str(blood_disease.df))
      # cat("\n")
      
      if(dim(blood_disease.df)[1] >0)
      {
        g <- ggplot(data=blood_disease.df, 
                    aes(x=reorder(VAR,pos37),
                        y=Disease_or_phenotype,
                        color=hgnc))+ 
          geom_jitter(size = 2,height = 0.15, width=0.15)+
          scale_color_manual(values = c("#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#3F4921", "#C0717C", "#CBD588", "#5F7FC7",
                                        "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD",
                                        "#D14285", "#6DDE88", "#652926", "#D1A33D", "#C84248",
                                        "azure2","chartreuse","cyan","blue4"),
                             drop=T)+
          scale_y_discrete(name=NULL, drop =F) +
          theme(axis.text.y  = element_text(angle=0, vjust = 0, hjust = 1,
                                            colour = c(rep("red",7),rep("steelblue",11),rep("#89C5DA",4),rep("#673770",3))))+
          scale_x_discrete(name=element_blank(), drop =FALSE) + 
          theme(axis.text.x  = element_text(angle=45, 
                                            vjust = 1, hjust = 1))+
          theme(legend.key=element_blank(), legend.key.size=unit(1,"point"), legend.title = element_blank(),
                legend.position="bottom")+
          geom_hline(yintercept = c(3,6,7,8,9,11,16), linetype="dotted", color="black")+
          geom_vline(xintercept = c(5,10,15,20,25,30,40,50,60,70,80,90,100,110,120), linetype="dotted", color="black")+
          ggtitle("Anchored GWAS disease variants Blood_diseases")
        
        g <-g + aes(reorder(VAR,pos37),stringr::str_wrap(Disease_or_phenotype, 50)) + ylab(NULL)
        
        print(g)
        
        list_blood_disease[[i]]<-blood_disease.df
        
        
      }
      
      inflammatory_diseases.df<-GWAS_anchored[which(GWAS_anchored$Disease_or_phenotype%in%inflammatory_diseases),]
      
      inflammatory_diseases.df$Disease_or_phenotype<-factor(inflammatory_diseases.df$Disease_or_phenotype,
                                                    levels=inflammatory_diseases,
                                                    ordered=T)
      
      # cat("inflammatory_diseases.df\n")
      # cat(str(inflammatory_diseases.df))
      # cat("\n")
      
      if(dim(inflammatory_diseases.df)[1] >0)
      {
        g <- ggplot(data=inflammatory_diseases.df, 
                    aes(x=reorder(VAR,pos37),
                        y=Disease_or_phenotype,
                        color=hgnc))+ 
          geom_jitter(size = 2,height = 0.15, width=0.15)+
          scale_color_manual(values = c("#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#3F4921", "#C0717C", "#CBD588", "#5F7FC7",
                                        "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD",
                                        "#D14285", "#6DDE88", "#652926", "#D1A33D", "#C84248",
                                        "azure2","chartreuse","cyan","blue4"),
                             drop=T)+
          scale_y_discrete(name=NULL, drop =F) +
          theme(axis.text.y  = element_text(angle=0, vjust = 0, hjust = 1,
                                            colour = c(rep("red",7),rep("steelblue",11),rep("#89C5DA",4),rep("#673770",3))))+
          scale_x_discrete(name=element_blank(), drop =FALSE) + 
          theme(axis.text.x  = element_text(angle=45, 
                                            vjust = 1, hjust = 1))+
          theme(legend.key=element_blank(), legend.key.size=unit(1,"point"), legend.title = element_blank(),
                legend.position="bottom")+
          geom_hline(yintercept = c(3,6,7,8,9,11,16), linetype="dotted", color="black")+
          geom_vline(xintercept = c(5,10,15,20,25,30,40,50,60,70,80,90,100,110,120), linetype="dotted", color="black")+
          ggtitle("Anchored GWAS disease variants inflammatory_diseases")
        
        g <-g + aes(reorder(VAR,pos37),stringr::str_wrap(Disease_or_phenotype, 50)) + ylab(NULL)
        
        print(g)
        
        list_inflammatory_diseases[[i]]<-inflammatory_diseases.df
        
        
      }
      
      immunodeficiencies.df<-GWAS_anchored[which(GWAS_anchored$Disease_or_phenotype%in%immunodeficiencies),]
      
      immunodeficiencies.df$Disease_or_phenotype<-factor(immunodeficiencies.df$Disease_or_phenotype,
                                                            levels=immunodeficiencies,
                                                            ordered=T)
      
      # cat("immunodeficiencies.df\n")
      # cat(str(immunodeficiencies.df))
      # cat("\n")
      
      if(dim(immunodeficiencies.df)[1] >0)
      {
        g <- ggplot(data=immunodeficiencies.df, 
                    aes(x=reorder(VAR,pos37),
                        y=Disease_or_phenotype,
                        color=hgnc))+ 
          geom_jitter(size = 2,height = 0.15, width=0.15)+
          scale_color_manual(values = c("#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#3F4921", "#C0717C", "#CBD588", "#5F7FC7",
                                        "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD",
                                        "#D14285", "#6DDE88", "#652926", "#D1A33D", "#C84248",
                                        "azure2","chartreuse","cyan","blue4"),
                             drop=T)+
          scale_y_discrete(name=NULL, drop =F) +
          theme(axis.text.y  = element_text(angle=0, vjust = 0, hjust = 1,
                                            colour = c(rep("red",7),rep("steelblue",11),rep("#89C5DA",4),rep("#673770",3))))+
          scale_x_discrete(name=element_blank(), drop =FALSE) + 
          theme(axis.text.x  = element_text(angle=45, 
                                            vjust = 1, hjust = 1))+
          theme(legend.key=element_blank(), legend.key.size=unit(1,"point"), legend.title = element_blank(),
                legend.position="bottom")+
          geom_hline(yintercept = c(3,6,7,8,9,11,16), linetype="dotted", color="black")+
          geom_vline(xintercept = c(5,10,15,20,25,30,40,50,60,70,80,90,100,110,120), linetype="dotted", color="black")+
          ggtitle("Anchored GWAS disease variants immunodeficiencies")
        
        g <-g + aes(reorder(VAR,pos37),stringr::str_wrap(Disease_or_phenotype, 50)) + ylab(NULL)
        
        print(g)
        
        list_immunodeficiencies[[i]]<-immunodeficiencies.df
        
      }
      
      leukemias.df<-GWAS_anchored[which(GWAS_anchored$Disease_or_phenotype%in%leukemias),]
      
      leukemias.df$Disease_or_phenotype<-factor(leukemias.df$Disease_or_phenotype,
                                                         levels=leukemias,
                                                         ordered=T)
      
      # cat("leukemias.df\n")
      # cat(str(leukemias.df))
      # cat("\n")
      
      if(dim(leukemias.df)[1] >0)
      {
        g <- ggplot(data=leukemias.df, 
                    aes(x=reorder(VAR,pos37),
                        y=Disease_or_phenotype,
                        color=hgnc))+ 
          geom_jitter(size = 2,height = 0.15, width=0.15)+
          scale_color_manual(values = c("#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#3F4921", "#C0717C", "#CBD588", "#5F7FC7",
                                        "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD",
                                        "#D14285", "#6DDE88", "#652926", "#D1A33D", "#C84248",
                                        "azure2","chartreuse","cyan","blue4"),
                             drop=T)+
          scale_y_discrete(name=NULL, drop =F) +
          theme(axis.text.y  = element_text(angle=0, vjust = 0, hjust = 1,
                                            colour = c(rep("red",7),rep("steelblue",11),rep("#89C5DA",4),rep("#673770",3))))+
          scale_x_discrete(name=element_blank(), drop =FALSE) + 
          theme(axis.text.x  = element_text(angle=45, 
                                            vjust = 1, hjust = 1))+
          theme(legend.key=element_blank(), legend.key.size=unit(1,"point"), legend.title = element_blank(),
                legend.position="bottom")+
          geom_hline(yintercept = c(3,6,7,8,9,11,16), linetype="dotted", color="black")+
          geom_vline(xintercept = c(5,10,15,20,25,30,40,50,60,70,80,90,100,110,120), linetype="dotted", color="black")+
          ggtitle("Anchored GWAS disease variants leukemias")
        
        g <-g + aes(reorder(VAR,pos37),stringr::str_wrap(Disease_or_phenotype, 50)) + ylab(NULL)
        
        print(g)
        
        list_leukemias[[i]]<-leukemias.df
        
        
      }
    
      miscellaneous.df<-GWAS_anchored[which(GWAS_anchored$Disease_or_phenotype%in%miscellaneous),]
      
      miscellaneous.df$Disease_or_phenotype<-factor(miscellaneous.df$Disease_or_phenotype,
                                                         levels=miscellaneous,
                                                         ordered=T)
      
      # cat("miscellaneous.df\n")
      # cat(str(miscellaneous.df))
      # cat("\n")
      
      if(dim(miscellaneous.df)[1] >0)
      {
        g <- ggplot(data=miscellaneous.df, 
                    aes(x=reorder(VAR,pos37),
                        y=Disease_or_phenotype,
                        color=hgnc))+ 
          geom_jitter(size = 2,height = 0.15, width=0.15)+
          scale_color_manual(values = c("#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#3F4921", "#C0717C", "#CBD588", "#5F7FC7",
                                        "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD",
                                        "#D14285", "#6DDE88", "#652926", "#D1A33D", "#C84248",
                                        "azure2","chartreuse","cyan","blue4"),
                             drop=T)+
          scale_y_discrete(name=NULL, drop =F) +
          theme(axis.text.y  = element_text(angle=0, vjust = 0, hjust = 1,
                                            colour = c(rep("red",7),rep("steelblue",11),rep("#89C5DA",4),rep("#673770",3))))+
          scale_x_discrete(name=element_blank(), drop =FALSE) + 
          theme(axis.text.x  = element_text(angle=45, 
                                            vjust = 1, hjust = 1))+
          theme(legend.key=element_blank(), legend.key.size=unit(1,"point"), legend.title = element_blank(),
                legend.position="bottom")+
          geom_hline(yintercept = c(3,6,7,8,9,11,16), linetype="dotted", color="black")+
          geom_vline(xintercept = c(5,10,15,20,25,30,40,50,60,70,80,90,100,110,120), linetype="dotted", color="black")+
          ggtitle("Anchored GWAS disease variants miscellaneous")
        
        g <-g + aes(reorder(VAR,pos37),stringr::str_wrap(Disease_or_phenotype, 50)) + ylab(NULL)
        
        print(g)
        
        
        list_miscellaneous[[i]]<-miscellaneous.df
        
      }
    }
    
    
    
    
    #quit(status=1)
    
    #### KEYS ----
    
    for(k in 1:length(KEYS))
    {
      KEYS_sel<-KEYS[k]
      
      cat("--------------------------------------------------------------------->KEYS_sel:\t")
      cat(sprintf(as.character(k)))
      cat("\t")
      cat(sprintf(as.character(KEYS_sel)))
      cat("\n")
      
      Rosetta_subset_sel_KEY<-Rosetta_subset_sel[which(Rosetta_subset_sel$KEY == KEYS_sel),]
      
      cat("Rosetta_subset_sel_KEY\n")
      cat(str(Rosetta_subset_sel_KEY))
      cat("\n")
      
      # setwd(out)
      # 
      # write.table(Rosetta_subset_sel_KEY,file="test.tsv", sep="\t", quote=F, row.names = F)
      carried_variants_sel<-unique(Rosetta_subset_sel_KEY$carried_variants)
      
      cat("carried_variants_sel\n")
      cat(str(carried_variants_sel))
      cat("\n")
      
      POS_array<-c(Rosetta_subset_sel_KEY$start,Rosetta_subset_sel_KEY$end)
      
      min_pos<-min(POS_array)-100
      max_pos<-max(POS_array)+100
      
      chrom2 = unique(Rosetta_subset_sel_KEY$chr)
      chromstart2 = min_pos
      chromend2 = max_pos
      
      cat("chrom2,chromstart2,chromend2\t")
      cat(str(chrom2))
      cat("\t")
      cat(str(chromstart2))
      cat("\t")
      cat(str(chromend2))
      cat("\n")
      
      for(h in 1:length(carried_variants_sel))
      {
        carried_variants_sel_h<-carried_variants_sel[h]
        
        cat("--------------------------------------------------------------------->carried_variants_sel_h:\t")
        cat(sprintf(as.character(h)))
        cat("\t")
        cat(sprintf(as.character(carried_variants_sel_h)))
        cat("\n")
        
        MPRA_sel_Carried_variants_sel<-Rosetta_subset_sel_KEY[which(Rosetta_subset_sel_KEY$carried_variants == carried_variants_sel_h),]
        
        MPRA_sel_Carried_variants_sel$Tile<-factor(MPRA_sel_Carried_variants_sel$Tile,
                                                   levels=c("TILE_1","TILE_2","TILE_3","TILE_4","TILE_5"),
                                                   ordered=T)
        MPRA_sel_Carried_variants_sel<-MPRA_sel_Carried_variants_sel[order(MPRA_sel_Carried_variants_sel$Tile),]
        
        cat("MPRA_sel_Carried_variants_sel\n")
        cat(str(MPRA_sel_Carried_variants_sel))
        cat("\n")
        
        
        #### BED GC content ----
        
        BED<- data.frame(matrix(vector(), dim(MPRA_sel_Carried_variants_sel)[1], 7,
                                dimnames=list(c(),
                                              c("chrom","start","end","name","score","strand","type"))),
                         stringsAsFactors=F)
        
        
        BED$chrom<-unique(MPRA_sel_Carried_variants_sel$chr)
        BED$start<-MPRA_sel_Carried_variants_sel$start
        BED$end<-MPRA_sel_Carried_variants_sel$end
        BED$name<-MPRA_sel_Carried_variants_sel$REAL_TILE
        BED$score<-MPRA_sel_Carried_variants_sel$GC_content
        BED$strand<-'1'
        BED$type<-'utr'
        
        
        cat("BED\n")
        cat(str(BED))
        cat("\n")
        
        pg = plotGenes(BED,chrom2,chromstart2,chromend2,
                       types = BED$type,
                       colorby=BED$score,
                       colorbycol= SushiColors(5),colorbyrange=c(GC_min,GC_max),
                       labeltext=TRUE,maxrows=50,height=0.4,plotgenetype="box")
        labelgenome(chrom2,chromstart2,chromend2,n=10,scale="bp")
        addlegend(pg[[1]],palette=pg[[2]],title=paste("GC_content"),side="right",
                  bottominset=0.4,topinset=0,xoffset=-.035,labelside="left",
                  width=0.025,title.offset=0.055)
        labelplot(paste("C) Tiles and GC content","in",carried_variants_sel_h, sep=" "))
        
        #### BED mA content ----
        
        # breaks.mA<-sort(unique(pointsmA))
        # log_breaks.mA<-log(breaks.mA+1)
        # labels.mA<-as.character(breaks.mA)
        
        BED<- data.frame(matrix(vector(), dim(MPRA_sel_Carried_variants_sel)[1], 7,
                                dimnames=list(c(),
                                              c("chrom","start","end","name","score","strand","type"))),
                         stringsAsFactors=F)


        BED$chrom<-unique(MPRA_sel_Carried_variants_sel$chr)
        BED$start<-MPRA_sel_Carried_variants_sel$start
        BED$end<-MPRA_sel_Carried_variants_sel$end
        BED$name<-MPRA_sel_Carried_variants_sel$REAL_TILE
        BED$score<-MPRA_sel_Carried_variants_sel$mA
        BED$strand<-'1'
        BED$type<-'utr'



        pg = plotGenes(BED,chrom2,chromstart2,chromend2,
                       types = BED$type,
                       colorby=log(BED$score +1),
                       colorbycol= SushiColors(5),colorbyrange=c(log_breaks.mA[1],log_breaks.mA[length(log_breaks.mA)]),
                       labeltext=TRUE,maxrows=50,height=0.4,plotgenetype="box")


        labelgenome( chrom2, chromstart2,chromend2,n=10,scale="bp")
        addlegend(pg[[1]],palette=pg[[2]],title=paste("multiple alignments"),side="right",
                  bottominset=0.4,topinset=0,xoffset=-.035,labelside="left",
                  width=0.025,title.offset=0.055)
        labelplot(paste("C) Tiles and multiple alignments","in",carried_variants_sel_h, sep=" "))
        
      }#h  carried_variants_sel
      
    }# k KEYS
    
    #### Finish printing ----
    
    if (makepdf == TRUE)
    {
      dev.off()
    }
   
    # if(dim(GWAS_anchored)[1] > 0)
    # {
    #  quit(status=1) 
    # }
  }#i

  #### SAVE RDS DEF ----
  
  setwd(out)
  
  filename18<-paste(type,"_","_charac_tiles",".rds",sep='')
  
  saveRDS(list_charac_tiles, file = filename18)
  
  filename18<-paste(type,"_","_hmap_Sel",".rds",sep='')
  
  saveRDS(list_hmap_Sel, file = filename18)
  
  filename18<-paste(type,"_","_hmap_Rosetta",".rds",sep='')
  
  saveRDS(list_hmap_Rosetta, file = filename18)
  
  filename18<-paste(type,"_","_blood_disease",".rds",sep='')
  
  saveRDS(list_blood_disease, file = filename18)
  
  filename18<-paste(type,"_","_inflammatory_diseases",".rds",sep='')
  
  saveRDS(list_inflammatory_diseases, file = filename18)
  
  filename18<-paste(type,"_","_immunodeficiencies",".rds",sep='')
  
  saveRDS(list_immunodeficiencies, file = filename18)
  
  filename18<-paste(type,"_","_leukemias",".rds",sep='')
  
  saveRDS(list_leukemias, file = filename18)
  
  filename18<-paste(type,"_","_miscellaneous",".rds",sep='')
  
  saveRDS(list_miscellaneous, file = filename18)
  
}

Global_report = function(option_list)
{
  detach("package:Sushi", unload=TRUE)
  library("plyr", lib.loc = "/nfs/users/nfs_m/mt19/sOFTWARE/R_libs/")
  library("dplyr", lib.loc = "/nfs/users/nfs_m/mt19/sOFTWARE/R_libs/")
  library("lemon")
  
  
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
  
  #### Selection_file ----
  
  Selection_file<-as.data.frame(fread(file = opt$Selection_file,
                                      sep="\t",
                                      header=T), stringsAsFactors=F)
  
  Selection_file$Allelic_Series_ID<-paste(Selection_file$phenotype,Selection_file$block_no, sep="__")
  
  
  cat("Selection_file\n")
  cat(str(Selection_file))
  cat("\n")
  
  #### read RDS ----
  
  setwd(out)
  
  filename18<-paste(type,"_","_charac_tiles",".rds",sep='')
  
  List_charac_tiles<-readRDS(file = filename18)
  
  # cat("List_charac_tiles\n")
  # cat(str(List_charac_tiles))
  # cat("\n")
  
  charac_tiles = as.data.frame(data.table::rbindlist(List_charac_tiles))
  
  cat("charac_tiles\n")
  cat(str(charac_tiles))
  cat("\n")
  
 # quit(status=1)
  
  
  filename18<-paste(type,"_","_hmap_Sel",".rds",sep='')
  
  List_hmap_Sel<-readRDS(file = filename18)
  
  # cat("List_hmap_Sel\n")
  # cat(str(List_hmap_Sel))
  # cat("\n")
  
  hmap_Sel = as.data.frame(data.table::rbindlist(List_hmap_Sel))
  
  cat("hmap_Sel\n")
  cat(str(hmap_Sel))
  cat("\n")
  
  
  
  filename18<-paste(type,"_","_hmap_Rosetta",".rds",sep='')
  
  List_hmap_Rosetta<-readRDS(file = filename18)

  hmap_Rosetta = as.data.frame(data.table::rbindlist(List_hmap_Rosetta))
  
  
  cat("hmap_Rosetta\n")
  cat(str(hmap_Rosetta))
  cat("\n")
  
  
  filename18<-paste(type,"_","_blood_disease",".rds",sep='')
  
  List_blood_disease<-readRDS(file = filename18)
  
  # cat("List_blood_disease\n")
  # cat(str(List_blood_disease))
  # cat("\n")
  
  if(length(List_blood_disease) >0 )
  {
    
    List_blood_disease<-List_blood_disease[-which(sapply(List_blood_disease, is.null))]
    
    # cat("List_blood_disease_not_null\n")
    # cat(str(List_blood_disease))
    # cat("\n")
    
    blood_disease = as.data.frame(data.table::rbindlist(List_blood_disease))
    
    cat("blood_disease\n")
    cat(str(blood_disease))
    cat("\n")
  }else{
    
    blood_disease <- data.frame(matrix(vector(), 0, 15,
                                               dimnames=list(c(),
                                                             c("ensembl_gene_id","hgnc","activity","Disease_or_phenotype","Classif",
                                                               "predicted_disease_association_id","rs_id","VAR_38","VAR","G2V_RS_value",
                                                               "G2V_consequence","pos37","color","VAR_IN_AS","Allelic_Series_ID"))),
                                        stringsAsFactors=F) 
  }
  
  
  
  filename18<-paste(type,"_","_inflammatory_diseases",".rds",sep='')
  
  List_inflammatory_diseases<-readRDS(file = filename18)
  
  # cat("List_inflammatory_diseases\n")
  # cat(str(List_inflammatory_diseases))
  # cat("\n")
  
      if(length(List_inflammatory_diseases) >0 )
      {
        
        List_inflammatory_diseases<-List_inflammatory_diseases[-which(sapply(List_inflammatory_diseases, is.null))]
        
        # cat("List_inflammatory_diseases_not_null\n")
        # cat(str(List_inflammatory_diseases))
        # cat("\n")
        
        inflammatory_diseases = as.data.frame(data.table::rbindlist(List_inflammatory_diseases))
        
        cat("inflammatory_diseases\n")
        cat(str(inflammatory_diseases))
        cat("\n")
      }else{
        
        inflammatory_diseases <- data.frame(matrix(vector(), 0, 15,
                                           dimnames=list(c(),
                                                         c("ensembl_gene_id","hgnc","activity","Disease_or_phenotype","Classif",
                                                           "predicted_disease_association_id","rs_id","VAR_38","VAR","G2V_RS_value",
                                                           "G2V_consequence","pos37","color","VAR_IN_AS","Allelic_Series_ID"))),
                                    stringsAsFactors=F) 
      }
  
  
  
  
  filename18<-paste(type,"_","_immunodeficiencies",".rds",sep='')
  
  List_immunodeficiencies<-readRDS(file = filename18)
  
  # cat("List_immunodeficiencies\n")
  # cat(str(List_immunodeficiencies))
  # cat("\n")
  
  if(length(List_immunodeficiencies) >0 )
  {
    
    List_immunodeficiencies<-List_immunodeficiencies[-which(sapply(List_immunodeficiencies, is.null))]
    
    # cat("List_immunodeficiencies_not_null\n")
    # cat(str(List_immunodeficiencies))
    # cat("\n")
    
    immunodeficiencies = as.data.frame(data.table::rbindlist(List_immunodeficiencies))
    
    cat("immunodeficiencies\n")
    cat(str(immunodeficiencies))
    cat("\n")
  }else{
    
    immunodeficiencies <- data.frame(matrix(vector(), 0, 15,
                                               dimnames=list(c(),
                                                             c("ensembl_gene_id","hgnc","activity","Disease_or_phenotype","Classif",
                                                               "predicted_disease_association_id","rs_id","VAR_38","VAR","G2V_RS_value",
                                                               "G2V_consequence","pos37","color","VAR_IN_AS","Allelic_Series_ID"))),
                                        stringsAsFactors=F) 
  }
  
  
  
  
  filename18<-paste(type,"_","_leukemias",".rds",sep='')
  
  List_leukemias<-readRDS(file = filename18)
  
  # cat("List_leukemias\n")
  # cat(str(List_leukemias))
  # cat("\n")
  
  if(length(List_leukemias) >0 )
  {
    
    List_leukemias<-List_leukemias[-which(sapply(List_leukemias, is.null))]
    
    # cat("List_leukemias_not_null\n")
    # cat(str(List_leukemias))
    # cat("\n")
    
    leukemias = as.data.frame(data.table::rbindlist(List_leukemias))
    
    cat("leukemias\n")
    cat(str(leukemias))
    cat("\n")
  }else{
    
    leukemias <- data.frame(matrix(vector(), 0, 15,
                                               dimnames=list(c(),
                                                             c("ensembl_gene_id","hgnc","activity","Disease_or_phenotype","Classif",
                                                               "predicted_disease_association_id","rs_id","VAR_38","VAR","G2V_RS_value",
                                                               "G2V_consequence","pos37","color","VAR_IN_AS","Allelic_Series_ID"))),
                                        stringsAsFactors=F) 
  }
  
  
  
  
  filename18<-paste(type,"_","_miscellaneous",".rds",sep='')
  
  List_miscellaneous<-readRDS(file = filename18)
  
  # cat("List_miscellaneous\n")
  # cat(str(List_miscellaneous))
  # cat("\n")
  
  if(length(List_miscellaneous) >0 )
  {
    
    List_miscellaneous<-List_miscellaneous[-which(sapply(List_miscellaneous, is.null))]
    
    # cat("List_miscellaneous_not_null\n")
    # cat(str(List_miscellaneous))
    # cat("\n")
    
    miscellaneous = as.data.frame(data.table::rbindlist(List_miscellaneous))
    
    cat("miscellaneous\n")
    cat(str(miscellaneous))
    cat("\n")
  }else{
    
    miscellaneous <- data.frame(matrix(vector(), 0, 15,
                                               dimnames=list(c(),
                                                             c("ensembl_gene_id","hgnc","activity","Disease_or_phenotype","Classif",
                                                               "predicted_disease_association_id","rs_id","VAR_38","VAR","G2V_RS_value",
                                                               "G2V_consequence","pos37","color","VAR_IN_AS","Allelic_Series_ID"))),
                                        stringsAsFactors=F) 
  }
  
  
  ### X and Y points ----
  
  pointsGC<-round(as.numeric(summary(hmap_Rosetta$GC_content)),1)
  GC_min<-0
  GC_max<-pointsGC[length(pointsGC)]+10
  
  pointsGC<-c(GC_min,round(as.numeric(summary(hmap_Rosetta$GC_content)[-4]),1),GC_max)
  
  breaks.GC<-sort(unique(c(pointsGC)))
  labels.GC<-as.character(breaks.GC)
  
  cat("labels.GC\n")
  cat(sprintf(as.character(labels.GC)))
  cat("\n")
  
  ###  mA
  
  
  pointsmA<-round(as.numeric(summary(hmap_Rosetta$mA)),0)
  mA_min<-0
  mA_max<-pointsmA[length(pointsmA)]+10000
  
  pointsmA<-c(mA_min,round(as.numeric(summary(hmap_Rosetta$mA)[-4]),0),mA_max)
  
  cat("pointsmA\n")
  cat(sprintf(as.character(pointsmA)))
  cat("\n")
  
  breaks.mA<-sort(unique(pointsmA,10,100,1000,10000,25000,50000))
  log_breaks.mA<-log(breaks.mA+1)
  labels.mA<-as.character(breaks.mA)
  
  cat("log_breaks.mA\n")
  cat(sprintf(as.character(log_breaks.mA)))
  cat("\n")
  
  cat("labels.mA\n")
  cat(sprintf(as.character(labels.mA)))
  cat("\n")
  
 
  ####  PLOT GC content Label and Factor 4 ----
  
  Total_Label<-as.data.frame(table(hmap_Rosetta$Label))
 
  colnames(Total_Label)<-c("Label","Total_Label")
  
  cat("Total_Label\n")
  cat(str(Total_Label))
  cat("\n")
  
  Total_factor4<-as.data.frame(table(hmap_Rosetta$factor4))
  
  colnames(Total_factor4)<-c("factor4","Total_factor4")
  
  cat("Total_factor4\n")
  cat(str(Total_factor4))
  cat("\n")
  
  hmap_Rosetta<-merge(hmap_Rosetta,
                      Total_Label,
                      by="Label",
                      all.x=T)
  
  cat("hmap_Rosetta_1\n")
  cat(str(hmap_Rosetta))
  cat("\n")
  
  hmap_Rosetta<-merge(hmap_Rosetta,
                      Total_factor4,
                      by="factor4",
                      all.x=T)
  
  cat("hmap_Rosetta_2\n")
  cat(str(hmap_Rosetta))
  cat("\n")
  
 
  
  g<-hmap_Rosetta  %>%
    mutate(myaxis = paste0(Label, "\n", "n=", Total_Label),drop=F) %>%
    mutate(myaxis=fct_reorder(myaxis,as.numeric(Label)),drop=F) %>%
    ggplot(aes(x=myaxis, y=GC_content, color=factor4)) +
    geom_sina()+
    scale_color_manual(values = c("#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#3F4921", "#C0717C", "#CBD588", "#5F7FC7",
                                  "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD",
                                  "#D14285", "#6DDE88", "#652926", "#D1A33D", "#C84248","#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#3F4921", "#C0717C", "#CBD588","#673770", "#D3D93E", "#38333E"),
                       drop=F)+
    theme(plot.title = element_text(size=11)) +
    ggtitle("GC content per label") +
    scale_x_discrete(name=NULL, drop=F)+
    theme(axis.text.x  = element_text(angle=45, vjust = 1, hjust = 1))+
    scale_y_continuous(name="% GC", breaks=breaks.GC, labels=labels.GC, limits=c(breaks.GC[1],breaks.GC[length(breaks.GC)]))+
    theme(legend.key=element_blank(), legend.position="bottom")+
    guides(fill=guide_legend(nrow=2,byrow=TRUE))
  
  g2<-hmap_Rosetta  %>%
    mutate(myaxis = paste0(factor4, "\n", "n=", Total_factor4),drop=F) %>%
    mutate(myaxis=fct_reorder(myaxis,as.numeric(factor4)),drop=F) %>%
    ggplot(aes(x=myaxis, y=GC_content, color=Label)) +
    geom_sina()+
    scale_color_manual(values = c("#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#3F4921", "#C0717C", "#CBD588", "#5F7FC7",
                                  "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD",
                                  "#D14285", "#6DDE88", "#652926", "#D1A33D", "#C84248","#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#3F4921", "#C0717C", "#CBD588","#673770", "#D3D93E", "#38333E"),
                       drop=F)+
    theme(plot.title = element_text(size=11)) +
    ggtitle("GC content per GC content bin") +
    scale_x_discrete(name=NULL, drop=F)+
    theme(axis.text.x  = element_text(angle=45, vjust = 1, hjust = 1))+
    scale_y_continuous(name="% GC", breaks=breaks.GC, labels=labels.GC, limits=c(breaks.GC[1],breaks.GC[length(breaks.GC)]))+
    theme(legend.key=element_blank(), legend.position="bottom")+
    guides(fill=guide_legend(nrow=2,byrow=TRUE))
  
  
  g.mA<-hmap_Rosetta  %>%
    mutate(myaxis = paste0(Label, "\n", "n=", Total_Label),drop=F) %>%
    mutate(myaxis=fct_reorder(myaxis,as.numeric(Label)),drop=F) %>%
    ggplot(aes(x=myaxis, y=log(mA+1), color=factor4)) +
    geom_boxplot()+
    scale_color_manual(values = c("#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#3F4921", "#C0717C", "#CBD588", "#5F7FC7",
                                  "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD",
                                  "#D14285", "#6DDE88", "#652926", "#D1A33D", "#C84248","#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#3F4921", "#C0717C", "#CBD588","#673770", "#D3D93E", "#38333E"),
                       drop=F)+
    theme(plot.title = element_text(size=11)) +
    ggtitle("Hits in GRCh37") +
    scale_x_discrete(name=NULL, drop=F)+
    theme(axis.text.x  = element_text(angle=45, vjust = 1, hjust = 1))+
    scale_y_continuous(name="# Hits", breaks=log_breaks.mA, labels=labels.mA, limits=c(log_breaks.mA[1],log_breaks.mA[length(log_breaks.mA)]))+
    theme(legend.key=element_blank(), legend.position="bottom")+
    guides(fill=guide_legend(nrow=2,byrow=TRUE))
  
  
  #### Frequency table Flag ----
    
  cat("levels(charac_tiles$variable\n")
  cat(sprintf(as.character(levels(charac_tiles$variable))))
  cat("\n") 
  
  Freq_table_Label<-unique(charac_tiles[which(charac_tiles$variable == "Label"),c(which(colnames(charac_tiles) == "value"),
                                                         which(colnames(charac_tiles) == "seq_name"))])
  
  colnames(Freq_table_Label)[which(colnames(Freq_table_Label) == "value")]<-"Label"
  
  cat("Freq_table_Label_\n")
  cat(str(Freq_table_Label))
  cat("\n") 
 
  Freq_table_Flag_factor<-unique(charac_tiles[which(charac_tiles$variable == "Flag_factor"),c(which(colnames(charac_tiles) == "value"),
                                                                                  which(colnames(charac_tiles) == "seq_name"))])
  
  cat("Freq_table_Flag_factor_\n")
  cat(str(Freq_table_Flag_factor))
  cat("\n") 
  
  colnames(Freq_table_Flag_factor)[which(colnames(Freq_table_Flag_factor) == "value")]<-"Flag_factor"
  
  
  Freq_Flag<-merge(Freq_table_Label,
                   Freq_table_Flag_factor,
                   by="seq_name",
                   all=T)
  
  
  cat("Freq_Flag_\n")
  cat(str(Freq_Flag))
  cat("\n") 
  
  DEF.Freq_Flag<-as.data.frame(t(table(Freq_Flag$Label,Freq_Flag$Flag_factor)))
  
  colnames(DEF.Freq_Flag)<-c("Flag_factor","Label","Freq")
  
  cat("DEF.Freq_Flag_\n")
  cat(str(DEF.Freq_Flag))
  cat("\n") 
  
  DEF.Freq_Flag.dt<-data.table(DEF.Freq_Flag, key="Label")

  cat("DEF.Freq_Flag.dt_\n")
  cat(str(DEF.Freq_Flag.dt))
  cat("\n")
  
  DEF.Freq_Flag.total<-as.data.frame(DEF.Freq_Flag.dt[, Total<-sum(Freq), by=Label])
  
  cat("DEF.Freq_Flag.total_\n")
  cat(str(DEF.Freq_Flag.total))
  cat("\n")
  
  colnames(DEF.Freq_Flag.total)[2]<-"Total"
  
  DEF.Freq_Flag<-merge(DEF.Freq_Flag,
                   DEF.Freq_Flag.total,
                   by="Label")
  
  #DEF.Freq_Flag<-DEF.Freq_Flag[-which(DEF.Freq_Flag$Total == 0),]
  
  DEF.Freq_Flag$Perc<-round(100*(DEF.Freq_Flag$Freq/DEF.Freq_Flag$Total),2)
  
  cat("DEF.Freq_Flag_2\n")
  cat(str(DEF.Freq_Flag))
  cat("\n") 
  
  DEF.Freq_Flag$Label<-factor(DEF.Freq_Flag$Label,
                              levels=c("Negative_Control_Genomic_Regions","UNDETERMINED_CTRL","ASE_CTRL","ASSAYED_VARIANT"),
                          ordered=T)
  
  DEF.Freq_Flag$Flag_factor<-factor(DEF.Freq_Flag$Flag_factor,
                              levels=levels(as.factor(DEF.Freq_Flag$Flag_factor)),
                              ordered=T)
  
  cat("DEF.Freq_Flag_LAST\n")
  cat(str(DEF.Freq_Flag))
  cat("\n")
  
  break.y=seq(0,110,by=10)
  labels.y=as.character(break.y)
  
  g.Flag_factor<-DEF.Freq_Flag %>%
    mutate(myaxis = paste0(Label, "\n", "n=", Total)) %>%
    mutate(myaxis=fct_reorder(myaxis,as.numeric(Label))) %>%
    ggplot(aes(x=myaxis, y=Perc, fill=Flag_factor)) +
    geom_bar(stat="identity")+
    scale_fill_viridis(discrete = TRUE, drop=F) +
    theme(plot.title = element_text(size=11)) +
    ggtitle("Flags in seq-name tiles") +
    scale_x_discrete(name=NULL, drop=F)+
    theme(axis.text.x  = element_text(angle=45, vjust = 1, hjust = 1))+
    scale_y_continuous(name="Percentage of the total seq-names in each class", breaks=break.y, labels=labels.y, limits=c(break.y[1],
                                                                                                                        break.y[length(break.y)]))+
    theme(legend.key=element_blank(), legend.position="bottom")+
    scale_fill_discrete(name="\n", drop=T)+
    guides(fill=guide_legend(nrow=2,byrow=TRUE))
  
  #### Frequency table fmt3 ----
  
  
  
  Freq_table_Label<-unique(charac_tiles[which(charac_tiles$variable == "Label"),c(which(colnames(charac_tiles) == "value"),
                                                                                  which(colnames(charac_tiles) == "seq_name"))])
  
  colnames(Freq_table_Label)[which(colnames(Freq_table_Label) == "value")]<-"Label"
  
  cat("Freq_table_Label_\n")
  cat(str(Freq_table_Label))
  cat("\n") 
  
  Freq_table_fmt3_factor<-unique(charac_tiles[which(charac_tiles$variable == "fmt3_factor"),c(which(colnames(charac_tiles) == "value"),
                                                                                              which(colnames(charac_tiles) == "seq_name"))])
  
  cat("Freq_table_fmt3_factor_\n")
  cat(str(Freq_table_fmt3_factor))
  cat("\n") 
  
  colnames(Freq_table_fmt3_factor)[which(colnames(Freq_table_fmt3_factor) == "value")]<-"fmt3_factor"
  
  
  Freq_fmt3<-merge(Freq_table_Label,
                   Freq_table_fmt3_factor,
                   by="seq_name",
                   all=T)
  
  
  cat("Freq_fmt3_\n")
  cat(str(Freq_fmt3))
  cat("\n") 
  
  DEF.Freq_fmt3<-as.data.frame(t(table(Freq_fmt3$Label,Freq_fmt3$fmt3_factor)))
  
  colnames(DEF.Freq_fmt3)<-c("fmt3_factor","Label","Freq")
  
  cat("DEF.Freq_fmt3_\n")
  cat(str(DEF.Freq_fmt3))
  cat("\n") 
  
  DEF.Freq_fmt3.dt<-data.table(DEF.Freq_fmt3, key="Label")
  
  cat("DEF.Freq_fmt3.dt_\n")
  cat(str(DEF.Freq_fmt3.dt))
  cat("\n")
  
  DEF.Freq_fmt3.total<-as.data.frame(DEF.Freq_fmt3.dt[, Total<-sum(Freq), by=Label])
  
  cat("DEF.Freq_fmt3.total_\n")
  cat(str(DEF.Freq_fmt3.total))
  cat("\n")
  
  colnames(DEF.Freq_fmt3.total)[2]<-"Total"
  
  DEF.Freq_fmt3<-merge(DEF.Freq_fmt3,
                       DEF.Freq_fmt3.total,
                       by="Label")
  
  #DEF.Freq_fmt3<-DEF.Freq_fmt3[-which(DEF.Freq_fmt3$Total == 0),]
  
  DEF.Freq_fmt3$Perc<-round(100*(DEF.Freq_fmt3$Freq/DEF.Freq_fmt3$Total),2)
  
  cat("DEF.Freq_fmt3_2\n")
  cat(str(DEF.Freq_fmt3))
  cat("\n") 
  
  DEF.Freq_fmt3$Label<-factor(DEF.Freq_fmt3$Label,
                              levels=c("Negative_Control_Genomic_Regions","UNDETERMINED_CTRL","ASE_CTRL","ASSAYED_VARIANT"),
                              ordered=T)
  
  DEF.Freq_fmt3$fmt3_factor<-factor(DEF.Freq_fmt3$fmt3_factor,
                                    levels=levels(as.factor(DEF.Freq_fmt3$fmt3_factor)),
                                    ordered=T)
  
  cat("DEF.Freq_fmt3_LAST\n")
  cat(str(DEF.Freq_fmt3))
  cat("\n")
  
  break.y=seq(0,110,by=10)
  labels.y=as.character(break.y)
  
  g.fmt3_factor<-DEF.Freq_fmt3 %>%
    mutate(myaxis = paste0(Label, "\n", "n=", Total)) %>%
    mutate(myaxis=fct_reorder(myaxis,as.numeric(Label))) %>%
    ggplot(aes(x=myaxis, y=Perc, fill=fmt3_factor)) +
    geom_bar(stat="identity")+
    scale_fill_viridis(discrete = TRUE, drop=F) +
    theme(plot.title = element_text(size=11)) +
    ggtitle("fmt3 reconstruction from BLAST in seq-name tiles") +
    scale_x_discrete(name=NULL, drop=F)+
    theme(axis.text.x  = element_text(angle=45, vjust = 1, hjust = 1))+
    scale_y_continuous(name="Percentage of the total seq-names in each class", breaks=break.y, labels=labels.y, limits=c(break.y[1],
                                                                                                                         break.y[length(break.y)]))+
    theme(legend.key=element_blank(), legend.position="bottom")+
    guides(fill=guide_legend(nrow=2,byrow=TRUE))
  
  
  
  #### DISEASES ----
  

  Gather<- data.frame(matrix(vector(), 0, 3,
                             dimnames=list(c(),
                                           c("Allelic_Series_ID","Disease_or_phenotype","group"))),
                      stringsAsFactors=F)
  
  cat("Gather_0\n")
  cat(str(Gather))
  cat("\n")
  
  if(dim(blood_disease)[1] >0)
  {
    
    blood_disease_subset<-unique(blood_disease[,which(colnames(blood_disease)%in%colnames(Gather))])
    blood_disease_subset$group<-"blood_disease"
    Gather<-rbind(Gather,blood_disease_subset)
    
  }# dim(blood_disease)[1]
  if(dim(inflammatory_diseases)[1] >0)
  {
    
    inflammatory_diseases_subset<-unique(inflammatory_diseases[,which(colnames(inflammatory_diseases)%in%colnames(Gather))])
    inflammatory_diseases_subset$group<-"inflammatory_diseases"
    Gather<-rbind(Gather,inflammatory_diseases_subset)
    
  }# dim(inflammatory_diseases)[1]
  if(dim(immunodeficiencies)[1] >0)
  {
    
    immunodeficiencies_subset<-unique(immunodeficiencies[,which(colnames(immunodeficiencies)%in%colnames(Gather))])
    immunodeficiencies_subset$group<-"immunodeficiencies"
    Gather<-rbind(Gather,immunodeficiencies_subset)
    
  }# dim(immunodeficiencies)[1]
  if(dim(leukemias)[1] >0)
  {
    
    leukemias_subset<-unique(leukemias[,which(colnames(leukemias)%in%colnames(Gather))])
    leukemias_subset$group<-"leukemias"
    Gather<-rbind(Gather,leukemias_subset)
    
  }# dim(leukemias)[1]
  if(dim(miscellaneous)[1] >0)
  {
    
    miscellaneous_subset<-unique(miscellaneous[,which(colnames(miscellaneous)%in%colnames(Gather))])
    miscellaneous_subset$group<-"miscellaneous"
    Gather<-rbind(Gather,miscellaneous_subset)
    
  }# dim(miscellaneous)[1]
  
  cat("Gather\n")
  cat(str(Gather))
  cat("\n")
  
  if(dim(Gather)[1]>0)
  {
    DEF.Gather<-as.data.frame(t(table(Gather$group,Gather$Disease_or_phenotype)))
    
    colnames(DEF.Gather)<-c("Disease_or_phenotype","group","Freq")
    
    cat("DEF.Gather_\n")
    cat(str(DEF.Gather))
    cat("\n") 
    
    DEF.Gather.dt<-data.table(DEF.Gather, key="group")
    
    cat("DEF.Gather.dt_\n")
    cat(str(DEF.Gather.dt))
    cat("\n")
    
    DEF.Gather.total<-as.data.frame(DEF.Gather.dt[, Total<-sum(Freq), by=group])
    
    cat("DEF.Gather.total_\n")
    cat(str(DEF.Gather.total))
    cat("\n")
    
    colnames(DEF.Gather.total)[2]<-"Total"
    
    DEF.Gather<-merge(DEF.Gather,
                      DEF.Gather.total,
                      by="group")
    
    #DEF.Gather<-DEF.Gather[-which(DEF.Gather$Total == 0),]
    
    DEF.Gather$Perc<-round(100*(DEF.Gather$Freq/DEF.Gather$Total),2)
    
    cat("DEF.Gather_2\n")
    cat(str(DEF.Gather))
    cat("\n") 
    
    DEF.Gather$group<-factor(DEF.Gather$group,
                             levels=c("blood_disease","inflammatory_diseases","immunodeficiencies","leukemias","miscellaneous"),
                             ordered=T)
    
    DEF.Gather$Disease_or_phenotype<-factor(DEF.Gather$Disease_or_phenotype,
                                            levels=levels(as.factor(DEF.Gather$Disease_or_phenotype)),
                                            ordered=T)
    
    cat("DEF.Gather_LAST\n")
    cat(str(DEF.Gather))
    cat("\n")
    
    break.y=seq(0,110,by=10)
    groups.y=as.character(break.y)
    
    # setwd(out)
    # 
    # write.table(DEF.Gather,file="test.tsv", sep="\t", quote=F, row.names = F)
    
    g.Disease_or_phenotype<-DEF.Gather %>%
      mutate(myaxis = paste0(group, "\n", "n=", Total)) %>%
      mutate(myaxis=fct_reorder(myaxis,as.numeric(group))) %>%
      ggplot(aes(x=myaxis, y=Perc, fill=Disease_or_phenotype)) +
      geom_bar(stat="identity")+
      theme(plot.title = element_text(size=11)) +
      ggtitle("Allelic Series anchored by at least 1 GWAS disease variant from OT") +
      scale_x_discrete(name=NULL, drop=F)+
      theme(axis.text.x  = element_text(angle=45, vjust = 1, hjust = 1))+
      scale_y_continuous(name="Percentage of the total Allelic Series in each class", breaks=break.y, labels=labels.y, limits=c(break.y[1],
                                                                                                                                break.y[length(break.y)]))+
      scale_fill_manual(values = c("#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#3F4921", "#C0717C", "#CBD588", "#5F7FC7",
                                   "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD",
                                   "#D14285", "#6DDE88", "#652926", "#D1A33D", "#C84248",
                                   "azure2","chartreuse","cyan","blue4",
                                   "#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#3F4921", "#C0717C", "#CBD588", "#5F7FC7",
                                   "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD",
                                   "#D14285", "#6DDE88", "#652926", "#D1A33D", "#C84248",
                                   "#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#3F4921", "#C0717C", "#CBD588", "#5F7FC7",
                                   "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD",
                                   "#D14285", "#6DDE88", "#652926", "#D1A33D", "#C84248",
                                   "azure2","chartreuse","cyan","blue4",
                                   "#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#3F4921", "#C0717C", "#CBD588", "#5F7FC7",
                                   "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD",
                                   "#D14285", "#6DDE88", "#652926", "#D1A33D", "#C84248",
                                   "#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#3F4921", "#C0717C", "#CBD588", "#5F7FC7",
                                   "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD",
                                   "#D14285", "#6DDE88", "#652926", "#D1A33D", "#C84248",
                                   "azure2","chartreuse","cyan","blue4",
                                   "#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#3F4921", "#C0717C", "#CBD588", "#5F7FC7",
                                   "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD",
                                   "#D14285", "#6DDE88", "#652926", "#D1A33D", "#C84248",
                                   "#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#3F4921", "#C0717C", "#CBD588", "#5F7FC7",
                                   "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD",
                                   "#D14285", "#6DDE88", "#652926", "#D1A33D", "#C84248",
                                   "azure2","chartreuse","cyan","blue4",
                                   "#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#3F4921", "#C0717C", "#CBD588", "#5F7FC7",
                                   "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD",
                                   "#D14285", "#6DDE88", "#652926", "#D1A33D", "#C84248"),
                        drop=T)+
      theme(legend.key=element_blank(), legend.position="hidden")+
      guides(fill=guide_legend(nrow=20,byrow=TRUE))
    
    #g <-g + aes(myaxis,Perc,stringr::str_wrap(Disease_or_phenotype, 50))
    
    legend <- g_legend(g.Disease_or_phenotype + theme(legend.position='bottom', legend.title=element_blank(), legend.text = element_text(size=5)))
    legend2 <- g_legend(g.Disease_or_phenotype+guides(fill=guide_legend(reverse = TRUE))  + theme(legend.position='bottom', legend.title=element_blank(), legend.text = element_text(size=5)))
    
    
  }
  
   
  setwd(out)
  
  pdf(file=paste(type,'_GLOBAL_REPORT','.pdf', sep=''))
  print(g)
  print(g2)
  print(g.mA)
  print(g.Flag_factor)
  print(g.fmt3_factor)
  if(dim(Gather)[1]>0)
  {
    print(g.Disease_or_phenotype)
    grid.newpage()
    grid.draw(legend2)
    grid.newpage()
    grid.draw(legend)
  }
  dev.off
  
  
   #### SAVE ----
  
  filename<-paste(type,"_","global_GC_and_multi_align",".tsv",sep='')
  
  write.table(hmap_Rosetta,file=filename, sep="\t",quote=F,row.names = F)
  
  filename<-paste(type,"_","global_charac_tile",".tsv",sep='')
  
  write.table(charac_tiles,file=filename, sep="\t",quote=F,row.names = F)
  
  filename<-paste(type,"_","global_diseases",".tsv",sep='')
  
  write.table(Gather,file=filename, sep="\t",quote=F,row.names = F)
  
  
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
    make_option(c("--Open_targets_GLOBAL_file"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Selection_file"), type="character", default=NULL, 
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
  
  
 Analysis_fmt3(opt)
 Analysis_fmt7(opt)
 Per_AS_report(opt)
 Global_report(opt)
  
}


###########################################################################

system.time( main() )
