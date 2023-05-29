
suppressMessages(library("dplyr", lib.loc="/nfs/users/nfs_m/mt19/sOFTWARE/R_libs/"))
suppressMessages(library("data.table", lib.loc="/nfs/users/nfs_m/mt19/sOFTWARE/R_libs/"))
suppressMessages(library("ggplot2", lib.loc = "/nfs/users/nfs_m/mt19/sOFTWARE/R_libs"))
suppressMessages(library("optparse", lib.loc = "/nfs/users/nfs_m/mt19/sOFTWARE/R_libs"))
suppressMessages(library("cowplot", lib.loc = "/nfs/users/nfs_m/mt19/sOFTWARE/R_libs/"))
suppressMessages(library("plyr", lib.loc="/nfs/users/nfs_m/mt19/sOFTWARE/R_libs/"))
suppressMessages(library("GenomicRanges"))
suppressMessages(library("rtracklayer"))

opt = NULL


Function_1_Assemble_elements = function(option_list)
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
  
  #### READ and transform subset ----
  
  subset = read.table(opt$subset, sep="\t", stringsAsFactors = F, header = T)  
  
  cat("subset\n")
  cat(str(subset))
  cat("\n")
  
 
  
  
  ### Read the transformed dB_ALL file ----
  
  dB_ALL = as.data.frame(fread(opt$dB_ALL, sep="\t", header = T), stringsAsFactors = F)
  
  cat("dB_ALL\n")
  str(dB_ALL)
  cat("\n")
  
  ### Read the transformed VJ_MPRA_1 file ----
  
  VJ_MPRA_1 = read.table(opt$VJ_MPRA_1, sep=",", stringsAsFactors = F, header = T)
  
  VJ_MPRA_1$variant<-paste('chr',VJ_MPRA_1$variant,sep='')
  
  VJ_MPRA_1$chr<-gsub(":.+$","",VJ_MPRA_1$variant)
  VJ_MPRA_1$pos<-gsub("^[^:]+:","",VJ_MPRA_1$variant)
  VJ_MPRA_1$pos<-gsub(":.+$","",VJ_MPRA_1$pos)
  VJ_MPRA_1$ref<-gsub("^[^:]+:[^:]+:","",VJ_MPRA_1$variant)
  VJ_MPRA_1$ref<-gsub(":.+$","",VJ_MPRA_1$ref)
  VJ_MPRA_1$alt<-gsub("^[^:]+:[^:]+:[^:]+:","",VJ_MPRA_1$variant)
  
  cat("VJ_MPRA_1\n")
  str(VJ_MPRA_1)
  cat("\n")
  
  ### Read the transformed VJ_MPRA_2 file ----
  
  VJ_MPRA_2 = read.table(opt$VJ_MPRA_2, sep=",", stringsAsFactors = F, header = T)
  
  VJ_MPRA_2<-VJ_MPRA_2[,-grep("GATA1",colnames(VJ_MPRA_2))]
  
  Q1<-summary(VJ_MPRA_2$CTRL.fc)[2]
  Q3<-summary(VJ_MPRA_2$CTRL.fc)[5]
  
  VJ_MPRA_2$chr<-paste('chr',VJ_MPRA_2$chr,sep='')
  
  VJ_MPRA_2$variant<-paste(as.character(VJ_MPRA_2$chr),
                           as.character(VJ_MPRA_2$pos),as.character(VJ_MPRA_2$ref),
                           as.character(VJ_MPRA_2$alt), sep=':')
  
  VJ_MPRA_2$name.and.window<-paste(VJ_MPRA_2$type,"Bot",VJ_MPRA_2$bot,"Top",
                                     VJ_MPRA_2$top,sep="_")
  
 
  cat("VJ_MPRA_2\n")
  str(VJ_MPRA_2)
  cat("\n")
  
  #### READ and transform Sankaran_significance_threshold ----
  
  Sankaran_significance_threshold = opt$Sankaran_significance_threshold
  
  cat("Sankaran_significance_threshold\n")
  cat(sprintf(as.character(Sankaran_significance_threshold)))
  cat("\n")
  
  
  ### Read the transformed Fulco_1 file ----
  
  Fulco_1 = read.table(opt$Fulco_1, sep=",", stringsAsFactors = F, header = T)
  
  cat("Fulco_1\n")
  str(Fulco_1)
  cat("\n")
  
  Fulco_1.No.NA<-Fulco_1[!is.na(Fulco_1$start),]
  
  cat("Fulco_1.No.NA\n")
  str(Fulco_1.No.NA)
  cat("\n")
  
  #### READ and transform Fulco_1_CRISPRi_threshold_UP ----
  
  Fulco_1_CRISPRi_threshold_UP = opt$Fulco_1_CRISPRi_threshold_UP
  
  cat("Fulco_1_CRISPRi_threshold_UP\n")
  cat(sprintf(as.character(Fulco_1_CRISPRi_threshold_UP)))
  cat("\n")
  
  #### READ and transform Fulco_1_CRISPRi_threshold_DOWN ----
  
  Fulco_1_CRISPRi_threshold_DOWN = opt$Fulco_1_CRISPRi_threshold_DOWN
  
  cat("Fulco_1_CRISPRi_threshold_DOWN\n")
  cat(sprintf(as.character(Fulco_1_CRISPRi_threshold_DOWN)))
  cat("\n")
  
  
  ### VJ Select by ASE effect sizes ----
  
  VJ_MPRA_2.Q1<-VJ_MPRA_2[(VJ_MPRA_2$CTRL.fc <= Q1),]
  VJ_MPRA_2.Q3<-VJ_MPRA_2[(VJ_MPRA_2$CTRL.fc >= Q3),]
  VJ_MPRA_2.BALANCED<-VJ_MPRA_2[(VJ_MPRA_2$CTRL.fc >= 0.95 &
                                   VJ_MPRA_2$CTRL.fc <= 1.05),]
  
  VJ_MPRA_2.DEF<-rbind(VJ_MPRA_2.Q1,
                       VJ_MPRA_2.Q3,
                       VJ_MPRA_2.BALANCED)
  
  VJ_MPRA_2_subset<-VJ_MPRA_2.DEF[(VJ_MPRA_2.DEF$CTRL.padj < Sankaran_significance_threshold),]
  
  
  cat("VJ_MPRA_2_subset\n")
  str(VJ_MPRA_2_subset)
  cat("\n")
  
  Mut_and_Ref.Rsid<-as.character(VJ_MPRA_2_subset$dbSNP[which(duplicated(VJ_MPRA_2_subset$dbSNP))])
  
  cat("Mut_and_Ref.Rsid\n")
  str(Mut_and_Ref.Rsid)
  cat("\n")
  
  selected_constructs_1<-unique(as.character(VJ_MPRA_2_subset$construct[VJ_MPRA_2_subset$dbSNP%in%Mut_and_Ref.Rsid]))
  
  cat("selected_constructs_1\n")
  str(selected_constructs_1)
  cat("\n")
  
  #### CRISPR VALIDATED VJ ---
  
  SK.CRISPR<-c("rs737092","rs1175550","rs1546723")
  
  VJ_MPRA_2_subset_CRISPR<-VJ_MPRA_2[which(VJ_MPRA_2$dbSNP%in%SK.CRISPR),]
  
  cat("VJ_MPRA_2_subset_CRISPR\n")
  str(VJ_MPRA_2_subset_CRISPR)
  cat("\n")
  
  VJ_MPRA_2_subset_CRISPR_MAX<-as.data.frame(setDT(VJ_MPRA_2_subset_CRISPR, key="dbSNP")[, .SD[CTRL.fc %in% max(CTRL.fc)], by=.(dbSNP)])
  
  cat("VJ_MPRA_2_subset_CRISPR_MAX\n")
  str(VJ_MPRA_2_subset_CRISPR_MAX)
  cat("\n")
  
  Mut_and_Ref.Rsid<-as.character(VJ_MPRA_2_subset_CRISPR_MAX$dbSNP[which(duplicated(VJ_MPRA_2_subset_CRISPR_MAX$dbSNP))])
  
  cat("Mut_and_Ref.Rsid\n")
  str(Mut_and_Ref.Rsid)
  cat("\n")
  
  selected_constructs_2<-unique(as.character(VJ_MPRA_2_subset$construct[VJ_MPRA_2_subset$dbSNP%in%Mut_and_Ref.Rsid]))
  
  cat("selected_constructs_2\n")
  str(selected_constructs_2)
  cat("\n")
  
  #### VJ positive final pick ----
  
  VJ_ASE_POSITIVE<-VJ_MPRA_2[which(VJ_MPRA_2$construct%in%selected_constructs_1 |
                                 VJ_MPRA_2$construct%in%selected_constructs_2),]
  
  VJ_ASE_POSITIVE$variant<-paste(VJ_ASE_POSITIVE$chr,VJ_ASE_POSITIVE$pos,
                              VJ_ASE_POSITIVE$ref,VJ_ASE_POSITIVE$alt, sep=':')
  
 
  
  cat("VJ_ASE_POSITIVE_0\n")
  str(VJ_ASE_POSITIVE)
  cat("\n")
  
  VJ_MPRA_1_subset<-VJ_MPRA_1[(VJ_MPRA_1$variant%in%VJ_ASE_POSITIVE$variant),]
  
  VJ_ASE_POSITIVE<-merge(VJ_MPRA_1_subset,VJ_ASE_POSITIVE,
                     by=c("chr","pos","ref","alt","variant","name.and.window"),all.y=T)
  
  cat("VJ_ASE_POSITIVE_1\n")
  str(VJ_ASE_POSITIVE)
  cat("\n")
  
  #### Random Pick of Sankaran undetermined ctrls ----
  
  VJ_MPRA_1_negative_ASE_subset<-VJ_MPRA_1[-which(VJ_MPRA_2$variant%in%VJ_MPRA_1$variant),]
  
  cat("VJ_MPRA_1_negative_ASE_subset_1\n")
  str(VJ_MPRA_1_negative_ASE_subset)
  cat("\n")
  
  x2 <- runif(10, 1, dim(VJ_MPRA_1_negative_ASE_subset)[1])
  
  x3<-round(x2,0)
  
  Random_pick<-VJ_MPRA_1_negative_ASE_subset[(x3),]
  #Random_pick$fake_construct<-gsub("^Ref_|^Mut_","",Random_pick$name.and.window)
  #Random_pick$fake_construct<-gsub("","",Random_pick$fake_construct)
  
  cat("Random_pick_1\n")
  str(Random_pick)
  cat("\n")
  
  #### VJ_FINAL ----
  
  VJ_FINAL<-merge(unique(VJ_ASE_POSITIVE[,-c(which(colnames(VJ_ASE_POSITIVE) == "X"),
                                      which(colnames(VJ_ASE_POSITIVE) == "X.1"))],),
                  unique(Random_pick[,-c(which(colnames(Random_pick) == "X"),
                                             which(colnames(Random_pick) == "X.1"))],),
                  by=c("variant","name.and.window","oligonucleotide",
                       "chr","pos","ref","alt"),
                  all=T)
  VJ_FINAL$Label<-"NA"
  
  VJ_FINAL$Label[is.na(VJ_FINAL$CTRL.fc)]<-"UNDETERMINED_CTRL"
  VJ_FINAL$Label[!is.na(VJ_FINAL$CTRL.fc)]<-"ASE_CTRL"
  
  VJ_FINAL$VAR<-paste(VJ_FINAL$chr,
                      VJ_FINAL$pos,
                      VJ_FINAL$ref,
                      VJ_FINAL$alt,
                      sep='_')
  
  cat("VJ_FINAL_1\n")
  str(VJ_FINAL)
  cat("\n")
  
  indx.int<-c(which(colnames(VJ_FINAL) =="chr"),which(colnames(VJ_FINAL) =="pos"),
              which(colnames(VJ_FINAL) =="ref"),which(colnames(VJ_FINAL) =="alt"),
              which(colnames(VJ_FINAL) =="VAR"),which(colnames(VJ_FINAL) =="Label"))
  
  VJ_FINAL_subset<-unique(VJ_FINAL[,indx.int])
  
  #### Fulco select ----
  
  Fulco_1.No.NA_subset<-Fulco_1.No.NA[(Fulco_1.No.NA$CRISPRi.Score > Fulco_1_CRISPRi_threshold_DOWN &
                                  Fulco_1.No.NA$CRISPRi.Score < Fulco_1_CRISPRi_threshold_UP),]
  
  cat("Fulco_1.No.NA_subset_0\n")
  str(Fulco_1.No.NA_subset)
  cat("\n")
  
  Fulco_1.No.NA_subset<-Fulco_1.No.NA_subset[which(Fulco_1.No.NA_subset$Set == "Negative Control Genomic Regions"),]
  
  cat("Fulco_1.No.NA_subset_1\n")
  str(Fulco_1.No.NA_subset)
  cat("\n")
  
  Fulco_1.No.NA_subset$pos<-Fulco_1.No.NA_subset$start+((Fulco_1.No.NA_subset$end - Fulco_1.No.NA_subset$start+1)/2)
  
  Fulco_1.No.NA_subset$Label<-"Negative_Control_Genomic_Regions"
  
  cat("Fulco_1.No.NA_subset_2\n")
  str(Fulco_1.No.NA_subset)
  cat("\n")
  
  ind.int<-c(which(colnames(Fulco_1.No.NA_subset) =="chr"),which(colnames(Fulco_1.No.NA_subset) =="pos"),which(colnames(Fulco_1.No.NA_subset) =="Label"))
  
  Fulco_1.No.NA_subset_2<-unique(Fulco_1.No.NA_subset[,ind.int])
  Fulco_1.No.NA_subset_2$pos<-as.integer(Fulco_1.No.NA_subset_2$pos)
  
  Fulco_1.No.NA_subset_2$ref<-"NA"
  Fulco_1.No.NA_subset_2$alt<-"NA"
  Fulco_1.No.NA_subset_2$VAR<-paste(Fulco_1.No.NA_subset_2$chr,Fulco_1.No.NA_subset_2$pos,Fulco_1.No.NA_subset_2$ref,Fulco_1.No.NA_subset_2$alt,sep="_")
  
  cat("Fulco_1.No.NA_subset_2\n")
  str(Fulco_1.No.NA_subset_2)
  cat("\n")
  
  #### candidates from subset ----
  
  
  indx.int<-c(which(colnames(dB_ALL) =="chr"),which(colnames(dB_ALL) =="pos37"),
              which(colnames(dB_ALL) =="ref"),which(colnames(dB_ALL) =="alt"),
              which(colnames(dB_ALL) =="VAR"))
  
  
  dB_ALL_subset<-unique(dB_ALL[which(dB_ALL$VAR%in%subset$VAR),indx.int])
  
  dB_ALL_subset$Label<-"ASSAYED_VARIANT"
  
  colnames(dB_ALL_subset)[which(colnames(dB_ALL_subset) =="pos37")]<-"pos"
  
  
  cat("dB_ALL_subset_1\n")
  str(dB_ALL_subset)
  cat("\n")
  
  #### LIBRARY_CANDIDATES ----
  
  indx.dep1<-which(dB_ALL_subset$VAR%in%VJ_FINAL_subset$VAR)
  
  cat("indx.dep1\n")
  cat(sprintf(as.character(indx.dep1)))
  cat("\n")
  
  indx.dep2<-which(dB_ALL_subset$VAR%in%VJ_FINAL_subset$VAR)
  
  cat("indx.dep2\n")
  cat(sprintf(as.character(indx.dep2)))
  cat("\n")
  
  indx.total<-unique(c(indx.dep1,indx.dep2))
  
  cat("indx.total\n")
  cat(sprintf(as.character(indx.total)))
  cat("\n")
  
  if(length(indx.total)>0)
  {
    dB_ALL_subset<-dB_ALL_subset[-indx.total,]
    
  }else{
    
    dB_ALL_subset<-dB_ALL_subset
  }
  
  
  
  
  cat("dB_ALL_subset_2\n")
  str(dB_ALL_subset)
  cat("\n")
  
  LIBRARY_CANDIDATES<-rbind(dB_ALL_subset,
                            VJ_FINAL_subset,
                            Fulco_1.No.NA_subset_2)
  
  cat("LIBRARY_CANDIDATES_1\n")
  str(LIBRARY_CANDIDATES)
  cat("\n")
  
  
  
  
  LIBRARY_CANDIDATES$KEY<-paste("Element_",seq(1,dim(LIBRARY_CANDIDATES)[1],by=1), sep='')
  
  cat("LIBRARY_CANDIDATES_2\n")
  str(LIBRARY_CANDIDATES)
  cat("\n")
    
  #### SAVE ----
  
  setwd(out)
  
  filename1=paste("LIBRARY_CANDIDATES_",type,".tsv", sep='')
  
  write.table(LIBRARY_CANDIDATES, file = filename1, row.names = F, quote = F, sep="\t")

}

Function_2_generate_VAR_carried_Var_bed_file = function(option_list)
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
  
  #### Read table ----
  
  setwd(out)
  
  filename1=paste("LIBRARY_CANDIDATES_",type,".tsv", sep='')
  
  LIBRARY_CANDIDATES<-read.table(file = filename1, sep="\t", header=T, stringsAsFactors = F)
  
  cat("LIBRARY_CANDIDATES\n")
  str(LIBRARY_CANDIDATES)
  cat("\n")
  
  #### Fulco_1_CRISPRi_threshold_UP ----
  
  slidding_windows = opt$slidding_windows
  
  cat("slidding_windows\n")
  cat(sprintf(as.character(slidding_windows)))
  cat("\n")
  
  
  #### span ----
  
  span = opt$span
  
  cat("span\n")
  cat(sprintf(as.character(span)))
  cat("\n")
  
  
  #### sunk_cost ----
  
  sunk_cost = opt$sunk_cost
  
  cat("sunk_cost\n")
  cat(sprintf(as.character(sunk_cost)))
  cat("\n")
  
  sunk_cost_unstring<-unlist(strsplit(sunk_cost, split=","))
  
  cat("sunk_cost_unstring\n")
  cat(sprintf(as.character(sunk_cost_unstring)))
  cat("\n")
  
  TOTAL_sunk_cost<-sum(as.numeric(sunk_cost_unstring))
  
  cat("TOTAL_sunk_cost\n")
  cat(sprintf(as.character(TOTAL_sunk_cost)))
  cat("\n")
  
  #### net_span ---- 
  
  net_span<-span-TOTAL_sunk_cost
  
  cat("net_span\n")
  cat(sprintf(as.character(net_span)))
  cat("\n")
  
  #### sunk_cost ----
  
  fractions = opt$fractions
  
  cat("fractions\n")
  cat(sprintf(as.character(fractions)))
  cat("\n")
  
  fractions_unstring<-as.numeric(unlist(strsplit(fractions, split=",")))
  
  cat("fractions_unstring\n")
  cat(sprintf(as.character(fractions_unstring)))
  cat("\n")
  
  
  ##### Loop to generate the Tiles ----
  
  gr_windows <- GRanges(seqnames=NULL,ranges=NULL,strand=NULL, names=NULL, names2=NULL)
  
  list_df<-list()
  
  for(i in 1:slidding_windows)
  {
    cat("window\t")
    cat(sprintf(as.character(i)))
    cat("\t")
    
    
    if(i == 1 | i == slidding_windows)
    {
      Fraction_sel<-fractions_unstring[1]
      
      cat("Fraction_sel\t")
      cat(sprintf(as.character(Fraction_sel)))
      cat("\t")
      
      step<-round(net_span/Fraction_sel,0)
      
      cat("step\t")
      cat(sprintf(as.character(step)))
      cat("\n")
      
      if(i == 1)
      {
        start_pos<-LIBRARY_CANDIDATES$pos-9*step
        start_pos[start_pos <= 0]<-1
        
        end_pos<-LIBRARY_CANDIDATES$pos+step
        end_pos[end_pos <= 0]<-1
        
        distances<-end_pos-start_pos+1
        distance_constraint<-rep(net_span,length(distances))
        distances_round_up_error<-distances-distance_constraint
        
        new_end_pos<-end_pos-distances_round_up_error
        
        def_start_pos<-start_pos
        def_end_pos<-new_end_pos
        
        
      }else{
        
        start_pos<-LIBRARY_CANDIDATES$pos-step
        start_pos[start_pos <= 0]<-1
        
        end_pos<-LIBRARY_CANDIDATES$pos+9*step
        end_pos[end_pos <= 0]<-1
        
        distances<-end_pos-start_pos+1
        distance_constraint<-rep(net_span,length(distances))
        distances_round_up_error<-distances-distance_constraint
        
        new_start_pos<-start_pos+distances_round_up_error
        
        def_start_pos<-new_start_pos
        def_end_pos<-end_pos
        
      }
      
      def_start_pos[def_start_pos <= 0]<-1
      def_end_pos[def_end_pos <= 0]<-1
      
      # cat("def_start_pos\n")
      # cat(sprintf(as.character(length(def_start_pos))))
      # cat("\n")
      # 
      # cat("def_end_pos\n")
      # cat(sprintf(as.character(length(def_end_pos))))
      # cat("\n")
      
      
      df<-data.frame(chr=LIBRARY_CANDIDATES$chr,
                     start=def_start_pos,
                     end=def_end_pos,
                     KEY=LIBRARY_CANDIDATES$KEY,
                     Tile=paste("TILE_",i,sep=''),
                     VAR=LIBRARY_CANDIDATES$VAR, stringsAsFactors = F)
      
      # cat("df\n")
      # str(df)
      # cat("\n")
      
      
      list_df[[i]]<-df
    }# i ==1 or
    
    if(i == 2 | i == (slidding_windows-1))
    {
      Fraction_sel<-fractions_unstring[2]
      
      cat("Fraction_sel\t")
      cat(sprintf(as.character(Fraction_sel)))
      cat("\t")
      
      step<-round(net_span/Fraction_sel,0)
      
      cat("step\t")
      cat(sprintf(as.character(step)))
      cat("\n")
      
      if(i == 2)
      {
        start_pos<-LIBRARY_CANDIDATES$pos-2*step
        start_pos[start_pos <= 0]<-1
        
        end_pos<-LIBRARY_CANDIDATES$pos+step
        end_pos[end_pos <= 0]<-1
        
        distances<-end_pos-start_pos+1
        distance_constraint<-rep(net_span,length(distances))
        distances_round_up_error<-distances-distance_constraint
        
        new_end_pos<-end_pos-distances_round_up_error
        
        def_start_pos<-start_pos
        def_end_pos<-new_end_pos
        
        
      }else{
        
        start_pos<-LIBRARY_CANDIDATES$pos-step
        start_pos[start_pos <= 0]<-1
        
        end_pos<-LIBRARY_CANDIDATES$pos+2*step
        end_pos[end_pos <= 0]<-1
        
        distances<-end_pos-start_pos+1
        distance_constraint<-rep(net_span,length(distances))
        distances_round_up_error<-distances-distance_constraint
        
        new_start_pos<-start_pos+distances_round_up_error
        
        def_start_pos<-new_start_pos
        def_end_pos<-end_pos
        
      }
      
      def_start_pos[def_start_pos <= 0]<-1
      def_end_pos[def_end_pos <= 0]<-1
      
      # cat("def_start_pos\n")
      # cat(sprintf(as.character(length(def_start_pos))))
      # cat("\n")
      # 
      # cat("def_end_pos\n")
      # cat(sprintf(as.character(length(def_end_pos))))
      # cat("\n")
      
      
      df<-data.frame(chr=LIBRARY_CANDIDATES$chr,
                     start=def_start_pos,
                     end=def_end_pos,
                     KEY=LIBRARY_CANDIDATES$KEY,
                     Tile=paste("TILE_",i,sep=''),
                     VAR=LIBRARY_CANDIDATES$VAR, stringsAsFactors = F)
      
      # cat("df\n")
      # str(df)
      # cat("\n")
      
      
      list_df[[i]]<-df
    }# i ==1 or 
    
    if(i == 3)
    {
      Fraction_sel<-fractions_unstring[3]
      
      cat("Fraction_sel\t")
      cat(sprintf(as.character(Fraction_sel)))
      cat("\t")
      
      step<-round(net_span/Fraction_sel,0)
      
      cat("step\t")
      cat(sprintf(as.character(step)))
      cat("\n")
      
      if(i == 3)
      {
        start_pos<-LIBRARY_CANDIDATES$pos-step
        start_pos[start_pos <= 0]<-1
        
        end_pos<-LIBRARY_CANDIDATES$pos+step
        end_pos[end_pos <= 0]<-1
        
        distances<-end_pos-start_pos+1
        distance_constraint<-rep(net_span,length(distances))
        distances_round_up_error<-distances-distance_constraint
        
        new_end_pos<-end_pos-distances_round_up_error
        
        def_start_pos<-start_pos
        def_end_pos<-new_end_pos
        
        
      }
      
      def_start_pos[def_start_pos <= 0]<-1
      def_end_pos[def_end_pos <= 0]<-1
      
      # cat("def_start_pos\n")
      # cat(sprintf(as.character(length(def_start_pos))))
      # cat("\n")
      # 
      # cat("def_end_pos\n")
      # cat(sprintf(as.character(length(def_end_pos))))
      # cat("\n")
      
      
      df<-data.frame(chr=LIBRARY_CANDIDATES$chr,
                     start=def_start_pos,
                     end=def_end_pos,
                     KEY=LIBRARY_CANDIDATES$KEY,
                     Tile=paste("TILE_",i,sep=''),
                     VAR=LIBRARY_CANDIDATES$VAR, stringsAsFactors = F)
      
      # cat("df\n")
      # str(df)
      # cat("\n")
      
      
      list_df[[i]]<-df
    }# i ==1 or 
    
  }#i slidding_windows i == slidding_windows
  
  TILES_PRE_OV = unique(as.data.frame(data.table::rbindlist(list_df)))
  
  TILES_PRE_OV$carried_variants<-TILES_PRE_OV$VAR
    
  cat("TILES_PRE_OV\n")
  str(TILES_PRE_OV)
  cat("\n")
  
  #### Partial overlaps in other tiles ----
  
  VARS<-unique(TILES_PRE_OV$VAR)
  
  # cat("VARS\n")
  # str(VARS)
  # cat("\n")
  
  Gather<- data.frame(matrix(vector(), 0, 7,
                             dimnames=list(c(),
                                           c("chr","start","end","KEY","Tile","VAR","carried_variants"))),
                      stringsAsFactors=F)

  for(i in 1:length(VARS))
  {
    # cat("VAR_sel: ")
    # cat(sprintf(as.character(i)))
    # cat("\t")
    
    VAR_sel<-VARS[i]
    
    # cat(sprintf(as.character(VAR_sel)))
    # cat("\n")
    
    chr_VAR<-gsub("_.+$","",VAR_sel)
    pos_VAR<-gsub("^[^_]+_","",VAR_sel)
    pos_VAR<-gsub("_.+$","",pos_VAR)
    
    
    gr_VAR <-GRanges(
      seqnames = as.character(chr_VAR),
      ranges=IRanges(
        start=as.numeric(pos_VAR),
        end=as.numeric(pos_VAR),
        names = VAR_sel))
    
    
    # cat("gr_VAR_PRE\n")
    # str(gr_VAR)
    # cat("\n")
    
  
    gr_SUBJECT <-GRanges(
      seqnames = as.character(TILES_PRE_OV$chr),
      names2=paste(TILES_PRE_OV$KEY, TILES_PRE_OV$Tile,sep='__'),
      ranges=IRanges(
        start=as.numeric(TILES_PRE_OV$start),
        end=as.numeric(TILES_PRE_OV$end),
        names = TILES_PRE_OV$VAR))
    
    # cat("gr_SUBJECT\n")
    # str(gr_SUBJECT)
    # cat("\n")
    
    hits <- findOverlaps(gr_VAR,gr_SUBJECT, type="any", minoverlap=1,)
    
    p <- Pairs(gr_VAR, gr_SUBJECT, hits=hits)
    
    # cat("p\n")
    # str(p)
    # cat("\n")
    
    TILES_POST_OV_sel <- data.frame(chr=seqnames(p@second),
                                    start=start(p@second),
                                    end=end(p@second),
                                    VAR=names(p@second),
                                    KEY_TILE=p@second@elementMetadata@listData$names2,
                                    stringsAsFactors = F)


    # cat("TILES_POST_OV_sel\n")
    # cat(str(TILES_POST_OV_sel))
    # cat("\n")
    
    TILES_POST_OV_sel_NO_SELF_OV<-TILES_POST_OV_sel[-which(TILES_POST_OV_sel$VAR==VAR_sel),]
    
    if(dim(TILES_POST_OV_sel_NO_SELF_OV)[1] >0)
    {
      
      TILES_POST_OV_sel_NO_SELF_OV$query_VAR<-VAR_sel
      
      
      # cat("TILES_POST_OV_sel_NO_SELF_OV\n")
      # cat(str(TILES_POST_OV_sel_NO_SELF_OV))
      # cat("\n")
      
      
      KEY_TILE_space<-unique(TILES_POST_OV_sel_NO_SELF_OV$KEY_TILE)
      
      for(k in 1:length(KEY_TILE_space))
      {
        KEY_TILE_sel<-KEY_TILE_space[k]
        
        # cat("KEY_TILE_sel: ")
        # cat(sprintf(as.character(KEY_TILE_sel)))
        # cat("\n")
        
        NO_SELF_OV_sel<-TILES_POST_OV_sel_NO_SELF_OV[which(TILES_POST_OV_sel_NO_SELF_OV$KEY_TILE == KEY_TILE_sel),]
        
        # cat("NO_SELF_OV_sel\n")
        # cat(str(NO_SELF_OV_sel))
        # cat("\n")
        
        NEW_ALLELES<-data.frame(chr=rep(NO_SELF_OV_sel$chr,2),
                                start=rep(NO_SELF_OV_sel$start,2),
                                end=rep(NO_SELF_OV_sel$end,2),
                                VAR=rep(NO_SELF_OV_sel$VAR,2),
                                KEY=rep(gsub("__.+$","",NO_SELF_OV_sel$KEY_TILE),2),
                                Tile=rep(gsub("^.+__","",NO_SELF_OV_sel$KEY_TILE),2),
                                carried_variants=c(NO_SELF_OV_sel$query_VAR, paste(NO_SELF_OV_sel$VAR,NO_SELF_OV_sel$query_VAR, sep="|")),
                                stringsAsFactors = F)
        
        # cat("NEW_ALLELES\n")
        # cat(str(NEW_ALLELES))
        # cat("\n")
        
        Gather<-rbind(Gather,NEW_ALLELES)
        
        ## Create a Gather new alleles and rbind all & at the end add all of them together
        
      }# k KEY_TILE_space
      
      #quit(status=1)
    }
    
  }#i VARS overlap other tiles
  
  Gather$chr<-as.character(Gather$chr)
  cat("Gather\n")
  str(Gather)
  cat("\n")
  
  DEF_TILES<-unique(rbind(TILES_PRE_OV,Gather))
  
  cat("DEF_TILES\n")
  str(DEF_TILES)
  cat("\n")
  
 
  
  ### export bed of tiles ----
  
  gr_DEF_TILES <- GRanges(
    seqnames = gsub("^chr","Chr",DEF_TILES$chr),
    ranges=IRanges(
      start=DEF_TILES$start, # note the difference
      end=DEF_TILES$end,
      name=paste(DEF_TILES$KEY,DEF_TILES$Tile, sep="__")))
  
  gr_DEF_TILES_unique<-unique(gr_DEF_TILES)
  
  setwd(out)
  
  filename_20<-paste(type,"_TILES",".bed", sep='')
  
  
  export.bed(gr_DEF_TILES_unique,con=filename_20)
  
  
  #### SAVE MASTER TABLE ----
  
  setwd(out)
  
  write.table(DEF_TILES, file=paste(type,"_MASTER_TILES",".tsv", sep=''), sep="\t", row.names = F, quote=F)
 
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
    make_option(c("--subset"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--dB_ALL"), type="character", default=NULL,
                metavar="FILE.txt",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--VJ_MPRA_1"), type="character", default=NULL,
                metavar="FILE.txt",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--VJ_MPRA_2"), type="character", default=NULL,
                metavar="FILE.txt",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Fulco_1"), type="character", default=NULL,
                metavar="FILE.txt",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Sankaran_significance_threshold"), type="numeric", default=NULL,
                metavar="FILE.txt",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Fulco_1_CRISPRi_threshold_UP"), type="numeric", default=NULL,
                metavar="FILE.txt",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Fulco_1_CRISPRi_threshold_DOWN"), type="numeric", default=NULL,
                metavar="FILE.txt",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--slidding_windows"), type="numeric", default=NULL,
                metavar="FILE.txt",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--span"), type="numeric", default=NULL,
                metavar="FILE.txt",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--sunk_cost"), type="character", default=NULL,
                metavar="FILE.txt",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--fractions"), type="character", default=NULL,
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
  
  
  Function_1_Assemble_elements(opt)
  Function_2_generate_VAR_carried_Var_bed_file(opt)
  
}


###########################################################################

system.time( main() )
