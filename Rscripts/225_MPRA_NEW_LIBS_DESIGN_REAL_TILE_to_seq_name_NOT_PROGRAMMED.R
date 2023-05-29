
suppressMessages(library("dplyr", lib.loc="/nfs/users/nfs_m/mt19/sOFTWARE/R_libs/"))
suppressMessages(library("data.table", lib.loc="/nfs/users/nfs_m/mt19/sOFTWARE/R_libs/"))
suppressMessages(library("ggplot2", lib.loc = "/nfs/users/nfs_m/mt19/sOFTWARE/R_libs"))
suppressMessages(library("optparse", lib.loc = "/nfs/users/nfs_m/mt19/sOFTWARE/R_libs"))
suppressMessages(library("plyr", lib.loc="/nfs/users/nfs_m/mt19/sOFTWARE/R_libs/"))
suppressMessages(library("BiocGenerics", lib.loc="/nfs/users/nfs_m/mt19/sOFTWARE/R_libs/"))
suppressMessages(library("S4Vectors", lib.loc="/nfs/users/nfs_m/mt19/sOFTWARE/R_libs/"))
suppressMessages(library("IRanges", lib.loc="/nfs/users/nfs_m/mt19/sOFTWARE/R_libs/"))
suppressMessages(library("XVector", lib.loc="/nfs/users/nfs_m/mt19/sOFTWARE/R_libs/"))
suppressMessages(library("Biostrings", lib.loc="/nfs/users/nfs_m/mt19/sOFTWARE/R_libs/"))
suppressMessages(library("seqinr", lib.loc="/nfs/users/nfs_m/mt19/sOFTWARE/R_libs/"))
suppressMessages(library("stringr", lib.loc="/nfs/users/nfs_m/mt19/sOFTWARE/R_libs/"))
suppressMessages(library("DNABarcodes", lib.loc="/nfs/users/nfs_m/mt19/sOFTWARE/R_libs/"))
suppressMessages(library("seqRFLP", lib.loc="/nfs/users/nfs_m/mt19/sOFTWARE/R_libs/"))



opt = NULL

Function_3_generate_long_matrix = function(option_list)
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
  cat(str(net_span))
  cat("\n")
  
  #### Read master tiles ----
  
  setwd(out)
  
  filename<-paste(type,"_MASTER_TILES_PLUS_REF_AND_NCGR_SUBS_PLUS_ALT_fa",".tsv", sep='')
  
  DEF_TILES<-read.table(file=filename, sep=" ", header=T, stringsAsFactors = F)
  
  cat("DEF_TILES\n")
  cat(str(DEF_TILES))
  cat("\n")
  
  # index.debug<-grep("Element_1__TILE_1",DEF_TILES$REAL_TILE)
  # 
  # cat("index.debug\n")
  # cat(str(index.debug))
  # cat("\n")
  
  #quit(status=1)
  #### LOOP to read the ALT fasta files ----
  
  Gather<- data.frame(matrix(vector(), 0, dim(DEF_TILES)[2]+1,
                             dimnames=list(c(),
                                           c(colnames(DEF_TILES),"ALT_sequence"))),
                      stringsAsFactors=F)
  
 for(i in 1:dim(DEF_TILES)[1])
 {
   DEF_TILES_sel<-DEF_TILES[i,]
   
   # cat("DEF_TILES_sel\n")
   # cat(str(DEF_TILES_sel))
   # cat("\n")
   
   ALT_fa_file<-DEF_TILES_sel$ALT_fa_file
   
   fastaFile<-readDNAStringSet(file=ALT_fa_file)
   
   # cat("fastaFile\n")
   # cat(str(fastaFile))
   # cat("\n")
   
   seq_name = names(fastaFile)
   sequence = paste(fastaFile)
   ALT_fasta <- data.frame(seq_name, sequence, stringsAsFactors = F)
   
   ALT_fasta$chr<-gsub(":.+$","",ALT_fasta$seq_name)
   ALT_fasta$chr<-gsub("^[^C]+Chr","chr",ALT_fasta$chr)
   ALT_fasta$start<-as.integer(gsub("^[^:]+:","",ALT_fasta$seq_name))
   ALT_fasta$ALT_nchar<-nchar(ALT_fasta$sequence)
   
   ### Correct big insertions or deletions that alter the net span
   # cat("ALT_fasta\n")
   # cat(str(ALT_fasta))
   # cat("\n")
   
   if(ALT_fasta$ALT_nchar > net_span)
   {
     # cat("ALT_fasta$ALT_nchar\n")
     # cat(sprintf(as.character(ALT_fasta$ALT_nchar)))
     # cat("\n")
     # cat(sprintf(as.character(ALT_fasta$sequence)))
     # cat("\n")
     # cat("net_span\n")
     # cat(sprintf(as.character(net_span)))
     # cat("\n")
     
     difference<-ALT_fasta$ALT_nchar-net_span
     
     # cat("difference\n")
     # cat(sprintf(as.character(difference)))
     # cat("\n")
     
     #truncated_ALT_fasta<-str_trunc(ALT_fasta$sequence, net_span, side = c("right"))
     truncated_ALT_fasta<-strtrim(ALT_fasta$sequence, net_span)
     
     # cat("truncated_ALT_fasta\n")
     # cat(sprintf(as.character(nchar(truncated_ALT_fasta))))
     # cat("\n")
     # cat(sprintf(as.character(truncated_ALT_fasta)))
     # cat("\n")
     
     # truncated_ALT_fasta<-gsub("\\.","",truncated_ALT_fasta)
     # cat(sprintf(as.character(nchar(truncated_ALT_fasta))))
     # cat("\n")
     # cat(sprintf(as.character(truncated_ALT_fasta)))
     # cat("\n")
     
     df<-cbind(DEF_TILES_sel,truncated_ALT_fasta)
     
     # cat("df\n")
     # cat(str(df))
     # cat("\n")
     
     colnames(df)<-colnames(Gather)
     
     Gather<-rbind(Gather,df)
     
    # quit(status = 1)
   }else{
    
    # cat("Hello_world\n")
     
     df<-cbind(DEF_TILES_sel,ALT_fasta$sequence)
     
     colnames(df)<-colnames(Gather)
     
     Gather<-rbind(Gather,df)
    
     # cat("Gather\n")
     # cat(str(Gather))
     # cat("\n")
    }
   
   # if(i == index.debug){quit(status = 1)}
 }#i
  
  
  ##### From REAL_TILE to seq-name ----
  
  cat("Gather\n")
  cat(str(Gather))
  cat("\n")
  
  
  REF_ind<-c(which(colnames(Gather) == "VAR"),which(colnames(Gather) == "chr"),which(colnames(Gather) == "start"),which(colnames(Gather) == "end"),
             which(colnames(Gather) == "KEY"),which(colnames(Gather) == "Tile"),which(colnames(Gather) == "REAL_TILE"),which(colnames(Gather) == "Label"),
             which(colnames(Gather) == "REF_sequence"))
  
  REF_df<-Gather[,REF_ind]
  
  cat("REF_df_0\n")
  cat(str(REF_df))
  cat("\n")
  
  REF_df<-unique(REF_df)
  REF_df$carried_variants<-"REF"
  colnames(REF_df)[which(colnames(REF_df) == "REF_sequence")]<-"sequence"
  REF_df$sequence<-as.character(REF_df$sequence)
  
  cat("REF_df_1\n")
  cat(str(REF_df))
  cat("\n")
  
 
  
  ALT_ind<-c(which(colnames(Gather) == "VAR"),which(colnames(Gather) == "chr"),which(colnames(Gather) == "start"),which(colnames(Gather) == "end"),
             which(colnames(Gather) == "KEY"),which(colnames(Gather) == "Tile"),which(colnames(Gather) == "REAL_TILE"),which(colnames(Gather) == "Label"),
             which(colnames(Gather) == "ALT_sequence"),which(colnames(Gather) == "carried_variants"))
  
  ALT_df<-Gather[,ALT_ind]
  
  cat("ALT_df_0\n")
  cat(str(ALT_df))
  cat("\n")
  
  ALT_df<-unique(ALT_df)
  
  colnames(ALT_df)[which(colnames(ALT_df) == "ALT_sequence")]<-"sequence"
  ALT_df$sequence<-as.character(ALT_df$sequence)
  
  cat("ALT_df_1\n")
  cat(str(ALT_df))
  cat("\n")
  
  
  Long_df<-rbind(REF_df,ALT_df)
  
  cat("Long_df_1\n")
  cat(str(Long_df))
  cat("\n")
  
  Long_df_ordered<-Long_df[order(Long_df$KEY,Long_df$Tile,Long_df$carried_variants),]
  
  cat("Long_df_ordered_1\n")
  cat(str(Long_df_ordered))
  cat("\n")
  
  #### GC content CLASSIFICATION ----
  
  
  GC<-NULL

  for(i in 1: dim(Long_df_ordered)[1])
  {
    dnastring = DNAString(Long_df_ordered$sequence[i])
    GC_content<-100*(letterFrequency(dnastring, "G", as.prob=TRUE) + letterFrequency(dnastring, "C", as.prob=TRUE))
    GC[i]<-GC_content

  }
  
  Long_df_ordered$GC_content<-GC
  
  cat("Long_df_ordered_2\n")
  cat(str(Long_df_ordered))
  cat("\n")
  
  A<-summary(Long_df_ordered$GC_content)[c(1,2,5,6)]
  A[1]<-A[1]-1
  A[4]<-A[4]+1
  A<-round(A,2)
  
  cat("A\n")
  cat(sprintf(as.character(A)))
  cat("\n")
  
  
  Long_df_ordered$GC_content_intervals<-cut(Long_df_ordered$GC_content,breaks = A,right = FALSE)
  
  B<-summary(as.factor(Long_df_ordered$GC_content_intervals))
  B_levels<-levels(as.factor(Long_df_ordered$GC_content_intervals))
  
  cat("B\n")
  cat(sprintf(as.character(B_levels)))
  cat("\n")
  cat(sprintf(as.character(B)))
  cat("\n")
  
  Long_df_ordered$factor4<-"NA"
  
  Long_df_ordered$factor4[Long_df_ordered$GC_content_intervals == B_levels[1]]<-"High_AT"
  Long_df_ordered$factor4[Long_df_ordered$GC_content_intervals == B_levels[2]]<-"Medium"
  Long_df_ordered$factor4[Long_df_ordered$GC_content_intervals == B_levels[3]]<-"High_GC"
  
  Long_df_ordered$factor4<-factor(Long_df_ordered$factor4,
                                     levels=c("High_AT","Medium","High_GC"),
                                     ordered=T)
  
  B<-summary(Long_df_ordered$factor4)
  B_levels<-levels(Long_df_ordered$factor4)
  
  cat("B\n")
  cat(sprintf(as.character(B_levels)))
  cat("\n")
  cat(sprintf(as.character(B)))
  cat("\n")
  
  #### seq_name ----
  
  Long_df_ordered$seq_name<-paste(Long_df_ordered$REAL_TILE,Long_df_ordered$carried_variants, sep="__")
  
  # indx.deplete<-c(which(colnames(Long_df_ordered) == "VAR"),
  #                 which(colnames(Long_df_ordered) == "chr"),
  #                 which(colnames(Long_df_ordered) == "start"))
  
  cat("Long_df_ordered_3\n")
  cat(str(Long_df_ordered))
  cat("\n")
  
  #### SAVE MASTER TABLE ----
  
  # quit(status = 1)
  
  setwd(out)
  
  filename<-paste(type,"_MASTER_TILES_PLUS_REF_AND_NCGR_SUBS_PLUS_ALT_fa_PLUS_ALL_FILES",".tsv", sep='')
  write.table(Gather, file=filename, sep="\t", row.names = F, quote=F)

  filename<-paste(type,"_LONG_matrix_seq_names",".tsv", sep='')
  write.table(Long_df_ordered, file=filename, sep="\t", row.names = F, quote=F)
  
  
  
}

Function_4_add_constant_sequences_and_index = function(option_list)
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
  
  #### RE_sites ----
  
  RE_sites = opt$RE_sites
  
  # cat("RE_sites\n")
  # cat(sprintf(as.character(RE_sites)))
  # cat("\n")
  
  RE_sites_unstring<-toupper(unlist(strsplit(RE_sites, split=",")))
  
  cat("RE_sites_unstring\n")
  cat(sprintf(as.character(RE_sites_unstring)))
  cat("\n")
  
  #### PCR_arms ----
  
  PCR_arms = opt$PCR_arms
  
  # cat("PCR_arms\n")
  # cat(sprintf(as.character(PCR_arms)))
  # cat("\n")
  
  PCR_arms_unstring<-toupper(unlist(strsplit(PCR_arms, split=",")))
  
  cat("PCR_arms_unstring\n")
  cat(sprintf(as.character(PCR_arms_unstring)))
  cat("\n")
  
  #### Blacklisted_sequences ----
  
  Blacklisted_sequences = opt$Blacklisted_sequences
  
  # cat("Blacklisted_sequences\n")
  # cat(sprintf(as.character(Blacklisted_sequences)))
  # cat("\n")
  
  Blacklisted_sequences_unstring<-toupper(unlist(strsplit(Blacklisted_sequences, split=",")))
  
  cat("Blacklisted_sequences_unstring\n")
  cat(sprintf(as.character(Blacklisted_sequences_unstring)))
  cat("\n")
  
  #### Read LONG matrix ----
  
  setwd(out)
  
  filename<-paste(type,"_LONG_matrix_seq_names",".tsv", sep='')
  
  LONG_matrix<-read.table(file=filename, sep="\t", header=T, stringsAsFactors = F)
  
  cat("LONG_matrix\n")
  cat(str(LONG_matrix))
  cat("\n")
  
  #### Add constant sequences to table -----
  
  LONG_matrix$UP_FWD_Primer<-"NA"
  LONG_matrix$UP_REV_Primer<-"NA"
  
  cat("LONG_matrix\n")
  cat(str(LONG_matrix))
  cat("\n")
  
  
  indx.GC_content_Medium<-which(LONG_matrix$factor4 == "Medium")
  
  cat("indx.GC_content_Medium\n")
  cat(sprintf(as.character(indx.GC_content_Medium)))
  cat("\n")
  cat(sprintf(as.character(length(indx.GC_content_Medium))))
  cat("\n")
  
  LONG_matrix$UP_FWD_Primer[indx.GC_content_Medium]<-rep(PCR_arms_unstring[3],length(indx.GC_content_Medium))
  LONG_matrix$UP_REV_Primer[indx.GC_content_Medium]<-rep(PCR_arms_unstring[4],length(indx.GC_content_Medium))
  
  
  indx.GC_content_High_GC<-which(LONG_matrix$factor4 == "High_GC")
  
  LONG_matrix$UP_FWD_Primer[indx.GC_content_High_GC]<-rep(PCR_arms_unstring[5],length(indx.GC_content_High_GC))
  LONG_matrix$UP_REV_Primer[indx.GC_content_High_GC]<-rep(PCR_arms_unstring[6],length(indx.GC_content_High_GC))
  
  
  indx.GC_content_High_AT<-which(LONG_matrix$factor4 == "High_AT")
  
  LONG_matrix$UP_FWD_Primer[indx.GC_content_High_AT]<-rep(PCR_arms_unstring[1],length(indx.GC_content_High_AT))
  LONG_matrix$UP_REV_Primer[indx.GC_content_High_AT]<-rep(PCR_arms_unstring[2],length(indx.GC_content_High_AT))
  
  
  LONG_matrix$BamHI<-RE_sites_unstring[1]
  LONG_matrix$KpnI<-RE_sites_unstring[2]
  
  cat("LONG_matrix\n")
  cat(str(LONG_matrix))
  cat("\n")
  
  ##### CREATE THE BLACKLIST ----
  
  Black.list<-unique(c(Blacklisted_sequences_unstring))


  
  cat("Black.list\n")
  cat(sprintf(as.character(Black.list)))
  cat("\n")
  cat(str(Black.list))
  cat("\n")

  comp_v<-NULL
  rev_comp_v<-NULL
  comp_rev_comp_v<-NULL
  
  for(i in 1:length(Black.list))
  {
      Black.list_sel<-Black.list[i]

      ## cat("Black.list_sel\n")
      ## cat(sprintf(as.character(Black.list_sel)))
      ## cat("\n")

      Black.list_sel_split<-unlist(strsplit(Black.list_sel, split=''))

      ## cat("Black.list_sel_split\n")
      ## cat(sprintf(as.character(Black.list_sel_split)))
      ## cat("\n")

      comp_Black.list_sel_split <- toupper(comp(Black.list_sel_split))



      rev_comp_Black.list_sel_split <- toupper(rev(comp_Black.list_sel_split))



      comp_rev_comp_Black.list_sel_split <- toupper(comp(rev_comp_Black.list_sel_split))


      comp_Black.list_sel_split<-paste(comp_Black.list_sel_split, collapse='')
      rev_comp_Black.list_sel_split<-paste(rev_comp_Black.list_sel_split, collapse='')
      comp_rev_comp_Black.list_sel_split<-paste(comp_rev_comp_Black.list_sel_split, collapse='')


      comp_v[i]<-comp_Black.list_sel_split
      rev_comp_v[i]<-rev_comp_Black.list_sel_split
      comp_rev_comp_v[i]<-comp_rev_comp_Black.list_sel_split

      ## cat("comp_Black.list_sel_split\n")
      ## cat(sprintf(as.character(comp_Black.list_sel_split)))
      ## cat("\n")

      ## cat("rev_comp_Black.list_sel_split\n")
      ## cat(sprintf(as.character(rev_comp_Black.list_sel_split)))
      ## cat("\n")

      ## cat("comp_rev_comp_Black.list_sel_split\n")
      ## cat(sprintf(as.character(comp_rev_comp_Black.list_sel_split)))
      ## cat("\n")
      

      ## quit(status=1)

  }#i Black.list

  
  Black.list <- unique(c(Black.list,comp_v,rev_comp_v,comp_rev_comp_v,PCR_arms_unstring,RE_sites_unstring))

  cat("Black.list\n")
  cat(sprintf(as.character(Black.list)))
  cat("\n")
  cat(str(Black.list))
  cat("\n")
  
#### FLAG patterns in the candidate sequences -----
  FLAG_v <- NULL
  
  for(i in 1:dim(LONG_matrix)[1])
  {
      LONG_matrix_sel<-LONG_matrix[i,]

      ## cat("LONG_matrix_sel\n")
      ## cat(str(LONG_matrix_sel))
      ## cat("\n")

      list_Flags<-list()

      for(k in 1:length(Black.list))
      {
           Black.list_sel<-Black.list[k]

 
           ME<-matchPattern(Black.list_sel, as.character(LONG_matrix_sel$sequence),
                   max.mismatch=0, min.mismatch=0,
                   with.indels=FALSE, fixed=TRUE,
                   algorithm="auto")


           Flag<-ME@ranges@start
           
           ## cat("ME\n")
           ## cat(str(ME))
           ## cat("\n")

           if(length(Flag)>0)
           {
               cat("Flag\n")
               cat(sprintf(as.character(Flag)))
               cat("\n")
               
               cat("Black.list_sel\n")
               cat(sprintf(as.character(Black.list_sel)))
               cat("\n")


               list_Flags[[k]]<-Black.list_sel

#               quit(status=1)

           }
           else
           {
               list_Flags[[k]]<-"NO_FLAG"
           } 
        
      }#k

     FLAG_content<-paste(unique(as.character(unlist(list_Flags))),collapse="\\|")

     FLAG_v[i]<-FLAG_content


#      quit(status=1)

  }#i pattern

  LONG_matrix$Flag_intra_candidate_pattern <- FLAG_v


  cat("LONG_matrix\n")
  cat(str(LONG_matrix))
  cat("\n")

  
#### SAVE MASTER TABLE ----
  
#  quit(status = 1)
  
  setwd(out)
  
  filename<-paste(type,"_LONG_matrix_seq_names_Plus_Flags",".tsv", sep='')
  write.table(LONG_matrix, file=filename, sep="\t", row.names = F, quote=F)

  filename<-paste(type,"_Black_list",".tsv", sep='')
  write.table(Black.list, file=filename, sep="\t", row.names = F, quote=F)
  
  
  
}

Function_6_add_index = function(option_list)
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



  #### Read LONG matrix ----

  setwd(out)


  filename<-paste(type,"_LONG_matrix_seq_names_Plus_Flags",".tsv", sep='')
  LONG_matrix<-read.table(file=filename, sep="\t", header=T, stringsAsFactors = F)

  cat("LONG_matrix\n")
  cat(str(LONG_matrix))
  cat("\n")



  LONG_matrix$proposed<-paste(LONG_matrix$UP_FWD_Primer,LONG_matrix$sequence,LONG_matrix$UP_REV_Primer,sep='')

#  LONG_matrix$proposed2<-paste(LONG_matrix$UP_FWD_Primer,LONG_matrix$barcode,LONG_matrix$BamHI,LONG_matrix$KpnI,LONG_matrix$sequence,LONG_matrix$UP_REV_Primer,sep='')
  
  LONG_matrix$seq_name<-paste(LONG_matrix$REAL_TILE,LONG_matrix$carried_variants,LONG_matrix$Flag_intra_candidate_pattern,sep='|')
 # LONG_matrix$seq_name_bc<-paste(LONG_matrix$seq_name,LONG_matrix$barcode,sep=':')

  cat("LONG_matrix2\n")
  cat(str(LONG_matrix))
  cat("\n")

                                        # Here we convert the data frame to a fasta object

  indx.fasta<-c(which(colnames(LONG_matrix) == "seq_name"),which(colnames(LONG_matrix) == "proposed"))
      
  Merge2_TO_FASTA<-unique(LONG_matrix[,indx.fasta])

  df.fasta = dataframe2fas(Merge2_TO_FASTA, file="df.fasta")

  # Here we print the fasta file

  writeLines(df.fasta, sep ="\n", paste(type,'_TWIST_LIBRARY','.fasta',sep=''))

  filename<-paste(type,"_MPRA_Rosetta",".tsv", sep='')
  write.table(LONG_matrix, file=filename, sep="\t", row.names = F, quote=F)
  
  #### Print copy for alignment
  
  indx.fasta<-c(which(colnames(LONG_matrix) == "seq_name"),which(colnames(LONG_matrix) == "sequence"))
  
  
  Merge2_TO_FASTA<-unique(LONG_matrix[which(LONG_matrix$carried_variants != "REF"),indx.fasta])
  
  df.fasta = dataframe2fas(Merge2_TO_FASTA, file="df.fasta")
  
  # Here we print the fasta file
  
  writeLines(df.fasta, sep ="\n", paste(type,'_alignment_check','.fasta',sep=''))
  
  Merge2_TO_FASTA<-unique(LONG_matrix[,indx.fasta])
  
  df.fasta = dataframe2fas(Merge2_TO_FASTA, file="df.fasta")
  
  # Here we print the fasta file
  
  writeLines(df.fasta, sep ="\n", paste(type,'_reference','.fasta',sep=''))

#### Parameters

  param_elements<-unique(LONG_matrix$KEY)

  cat("param_elements:\t")
  cat(sprintf(as.character(length(param_elements))))
cat("\n")

  param_REAL_TILES<-unique(LONG_matrix$REAL_TILE)

  cat("param_REAL_TILES:\t")
  cat(sprintf(as.character(length(param_REAL_TILES))))
cat("\n")

  param_seq_name<-unique(LONG_matrix$seq_name)

  cat("param_seq_name:\t")
  cat(sprintf(as.character(length(param_seq_name))))
cat("\n")

  param_seq_name_bc<-unique(LONG_matrix$seq_name_bc)

  cat("param_seq_name_bc:\t")
  cat(sprintf(as.character(length(param_seq_name_bc))))
  cat("\n")
  

#  quit(status=1)

  
}

Function_7_printing_primers = function(option_list)
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
  
  
  
  #### Read LONG matrix ----
  
  setwd(out)
  
  
  filename<-paste(type,"_MPRA_Rosetta",".tsv", sep='')
  LONG_matrix<-read.table(file=filename, sep="\t", header=T, stringsAsFactors = F)
  
  cat("LONG_matrix\n")
  cat(str(LONG_matrix))
  cat("\n")
  
  
  LA_df<-data.frame(matrix(vector(), 0, 2,
                                      dimnames=list(c(),c("seq_name","sequence"))),
                               stringsAsFactors=F)
  
  #### LIBRARY AMPLIFICATION PCR #1 ----
  
  indx.int<-c(which(colnames(LONG_matrix) == "factor4"),
              which(colnames(LONG_matrix) == "UP_FWD_Primer"),
              which(colnames(LONG_matrix) == "UP_REV_Primer"),
              which(colnames(LONG_matrix) == "BamHI"),
              which(colnames(LONG_matrix) == "KpnI"))
  
  PCR_arms_df<-unique(LONG_matrix[,indx.int])
  
  REV_Primer_revComp<-NULL
  
  for(i in 1:dim(PCR_arms_df)[1])
  {
    REV_sel<-PCR_arms_df$UP_REV_Primer[i]
    
    cat("--->\t")
    cat(sprintf(as.character(i)))
    cat("\t")
    cat("REV_sel_\n")
    cat(sprintf(as.character(REV_sel)))
    cat("\n")
    
    REV_Primer_revComp[i]<-revComp(REV_sel)
  }# i
  
  PCR_arms_df$UP_REV_Primer_revComp<-REV_Primer_revComp
  
  cat("PCR_arms_df\n")
  cat(str(PCR_arms_df))
  cat("\n")
  
  FWD_PRIMERS_1_df<-data.frame(matrix(vector(), length(PCR_arms_df$factor4), 4,
                                      dimnames=list(c(),c("seq_name","BamHI","KpnI","UP_FWD_Primer"))),
                               stringsAsFactors=F)
  
  # cat("FWD_PRIMERS_0_df\n")
  # cat(str(FWD_PRIMERS_1_df))
  # cat("\n")
  
  FWD_PRIMERS_1_df$UP_FWD_Primer<-PCR_arms_df$UP_FWD_Primer
  FWD_PRIMERS_1_df$KpnI<-PCR_arms_df$KpnI
  FWD_PRIMERS_1_df$BamHI<-PCR_arms_df$BamHI
  FWD_PRIMERS_1_df$seq_name<-paste("LA_1_FWD",PCR_arms_df$factor4,sep='_')
  FWD_PRIMERS_1_df$sequence<-paste(FWD_PRIMERS_1_df$BamHI,
                                   FWD_PRIMERS_1_df$KpnI,
                                   FWD_PRIMERS_1_df$UP_FWD_Primer,
                                   sep='')
  
  cat("FWD_PRIMERS_1_df\n")
  cat(str(FWD_PRIMERS_1_df))
  cat("\n")
  
  indx.int<-c(which(colnames(FWD_PRIMERS_1_df) == "seq_name"),which(colnames(FWD_PRIMERS_1_df) == "sequence"))
  
  FWD_PRIMERS_1_df_subset<-unique(FWD_PRIMERS_1_df[,indx.int])
  
  cat("FWD_PRIMERS_1_df_subset\n")
  cat(str(FWD_PRIMERS_1_df_subset))
  cat("\n")
  
  LA_df<-rbind(LA_df,FWD_PRIMERS_1_df_subset)
  
  cat("LA_df\n")
  cat(str(LA_df))
  cat("\n")
  
  REV_PRIMERS_1_df<-data.frame(matrix(vector(), length(PCR_arms_df$factor4), 3,
                                      dimnames=list(c(),c("seq_name","UP_REV_Primer_revComp","Priming_REV_1"))),
                               stringsAsFactors=F)
  
  REV_PRIMERS_1_df$UP_REV_Primer_revComp<-PCR_arms_df$UP_REV_Primer_revComp
  REV_PRIMERS_1_df$Priming_REV_1<-toupper("gtcgaTCCTGGCCTAGTTG")
  REV_PRIMERS_1_df$seq_name<-paste("LA_1_REV",PCR_arms_df$factor4,sep='_')
  REV_PRIMERS_1_df$sequence<-paste(REV_PRIMERS_1_df$Priming_REV_1,
                                   REV_PRIMERS_1_df$UP_REV_Primer_revComp,
                                   sep='')
  
  
  
  
  cat("REV_PRIMERS_1_df\n")
  cat(str(REV_PRIMERS_1_df))
  cat("\n")
  
  indx.int<-c(which(colnames(REV_PRIMERS_1_df) == "seq_name"),which(colnames(REV_PRIMERS_1_df) == "sequence"))
  
  REV_PRIMERS_1_df_subset<-unique(REV_PRIMERS_1_df[,indx.int])
  
  cat("REV_PRIMERS_1_df_subset\n")
  cat(str(REV_PRIMERS_1_df_subset))
  cat("\n")
  
  LA_df<-rbind(LA_df,REV_PRIMERS_1_df_subset)
  
  cat("LA_df\n")
  cat(str(LA_df))
  cat("\n")
  
  #### LIBRARY AMPLIFICATION PCR #2 ----
  
  
  FWD_PRIMERS_2_df<-data.frame(matrix(vector(), length(PCR_arms_df$factor4), 5,
                                      dimnames=list(c(),c("seq_name","Gibson_5_prime","Spacer","BC","Priming"))),
                               stringsAsFactors=F)
  
  # cat("FWD_PRIMERS_0_df\n")
  # cat(str(FWD_PRIMERS_2_df))
  # cat("\n")
  
  FWD_PRIMERS_2_df$Priming<-FWD_PRIMERS_1_df_subset$sequence
  FWD_PRIMERS_2_df$seq_name<-FWD_PRIMERS_1_df_subset$seq_name
  FWD_PRIMERS_2_df$seq_name<-gsub("^LA_1_FWD","LA_2_FWD",FWD_PRIMERS_2_df$seq_name)
  FWD_PRIMERS_2_df$BC<-paste(rep("N",15),collapse='')
  FWD_PRIMERS_2_df$Spacer<-toupper("GCAAAGTGAACACATCGCTAAGCGAAAGCTAAG")
  FWD_PRIMERS_2_df$Gibson_5_prime<-toupper("tagagcatgcaccggtgata")
  
  FWD_PRIMERS_2_df$sequence<-paste(FWD_PRIMERS_2_df$Gibson_5_prime,
                                   FWD_PRIMERS_2_df$Spacer,
                                   FWD_PRIMERS_2_df$BC,
                                   FWD_PRIMERS_2_df$Priming,
                                   sep='')
  
  cat("FWD_PRIMERS_2_df\n")
  cat(str(FWD_PRIMERS_2_df))
  cat("\n")
  
  indx.int<-c(which(colnames(FWD_PRIMERS_2_df) == "seq_name"),which(colnames(FWD_PRIMERS_2_df) == "sequence"))
  
  FWD_PRIMERS_2_df_subset<-unique(FWD_PRIMERS_2_df[,indx.int])
  
  cat("FWD_PRIMERS_2_df_subset\n")
  cat(str(FWD_PRIMERS_2_df_subset))
  cat("\n")
  
  LA_df<-rbind(LA_df,FWD_PRIMERS_2_df_subset)
  
  cat("LA_df\n")
  cat(str(LA_df))
  cat("\n")
  
  REV_PRIMERS_2_df<-data.frame(matrix(vector(), length(PCR_arms_df$factor4), 3,
                                      dimnames=list(c(),c("seq_name","Gibson_3_prime","Priming"))),
                               stringsAsFactors=F)
  
  # cat("FWD_PRIMERS_0_df\n")
  # cat(str(FWD_PRIMERS_2_df))
  # cat("\n")
  
  REV_PRIMERS_2_df$Priming<-REV_PRIMERS_1_df_subset$sequence
  REV_PRIMERS_2_df$seq_name<-REV_PRIMERS_1_df_subset$seq_name
  REV_PRIMERS_2_df$seq_name<-gsub("^LA_1_REV","LA_2_REV",REV_PRIMERS_2_df$seq_name)
  REV_PRIMERS_2_df$Gibson_3_prime<-revComp(toupper("gaattcggccggcag"))
  
  REV_PRIMERS_2_df$sequence<-paste(REV_PRIMERS_2_df$Gibson_3_prime,
                                   REV_PRIMERS_2_df$Priming,
                                   sep='')
  
  
  
  
  cat("REV_PRIMERS_2_df\n")
  cat(str(REV_PRIMERS_2_df))
  cat("\n")
  
  indx.int<-c(which(colnames(REV_PRIMERS_2_df) == "seq_name"),which(colnames(REV_PRIMERS_2_df) == "sequence"))
  
  REV_PRIMERS_2_df_subset<-unique(REV_PRIMERS_2_df[,indx.int])
  
  cat("REV_PRIMERS_2_df_subset\n")
  cat(str(REV_PRIMERS_2_df_subset))
  cat("\n")
  
  LA_df<-rbind(LA_df,REV_PRIMERS_2_df_subset)
  
  cat("LA_df\n")
  cat(str(LA_df))
  cat("\n")
  
  #### Print copy for alignment
  
  # indx.fasta<-c(which(colnames(LONG_matrix) == "seq_name"),which(colnames(LONG_matrix) == "sequence"))
  # 
  # 
  # Merge2_TO_FASTA<-unique(LONG_matrix[which(LONG_matrix$carried_variants != "REF"),indx.fasta])
  # 
  # df.fasta = dataframe2fas(Merge2_TO_FASTA, file="df.fasta")
  # 
  # # Here we print the fasta file
  # 
  # writeLines(df.fasta, sep ="\n", paste(type,'_alignment_check','.fasta',sep=''))
  # 
  
  
  
  #  quit(status=1)
  
  
  
  #### LIBRARY PREPARATION 150 ----
  
  LP_150_df<-data.frame(matrix(vector(), 0, 2,
                           dimnames=list(c(),c("seq_name","sequence"))),
                    stringsAsFactors=F)
  
  FWD_PRIMERS_150_df<-data.frame(matrix(vector(), 2, 3,
                                      dimnames=list(c(),c("seq_name","P5","Priming"))),
                               stringsAsFactors=F)
  
  FWD_PRIMERS_150_df$Priming<-rep(toupper("GCAAAGTGAACACATCGCTAAGCGAAAGCTAAG"),2)
  FWD_PRIMERS_150_df$P5<-rep(toupper("AATGATACGGCGACCACCGAGATCTACAC"),2)
  
  FWD_PRIMERS_150_df$seq_name<-c("LP_PE150_FWD","P5")
  
  FWD_PRIMERS_150_df$sequence<-c(paste(FWD_PRIMERS_150_df$P5[1],FWD_PRIMERS_150_df$Priming[1], sep=''),
                                 FWD_PRIMERS_150_df$P5[1])
  
  cat("FWD_PRIMERS_150_df\n")
  cat(str(FWD_PRIMERS_150_df))
  cat("\n")
  
  indx.int<-c(which(colnames(FWD_PRIMERS_150_df) == "seq_name"),which(colnames(FWD_PRIMERS_150_df) == "sequence"))
  
  FWD_PRIMERS_150_df_subset<-unique(FWD_PRIMERS_150_df[,indx.int])
  
  cat("FWD_PRIMERS_150_df_subset\n")
  cat(str(FWD_PRIMERS_150_df_subset))
  cat("\n")
  
  LP_150_df<-rbind(LP_150_df,FWD_PRIMERS_150_df_subset)
  
  cat("LP_150_df\n")
  cat(str(LP_150_df))
  cat("\n")
  
  REV_PRIMERS_150_df<-data.frame(matrix(vector(), 2, 3,
                                        dimnames=list(c(),c("seq_name","P7","Priming"))),
                                 stringsAsFactors=F)
  
  REV_PRIMERS_150_df$Priming<-rep(toupper("CAACTAGGCCAGGAtcgac"),2)
  REV_PRIMERS_150_df$P7<-rep(toupper("ATCTCGTATGCCGTCTTCTGCTTG"),2)
  
  REV_PRIMERS_150_df$seq_name<-c("LP_PE150_REV","P7")
  
  REV_PRIMERS_150_df$revcomp_Priming<-sapply(REV_PRIMERS_150_df$Priming, revComp)
  REV_PRIMERS_150_df$revcomp_P7<-sapply(REV_PRIMERS_150_df$P7, revComp)
  
  
  
  REV_PRIMERS_150_df$sequence<-c(paste(REV_PRIMERS_150_df$revcomp_P7[1],REV_PRIMERS_150_df$revcomp_Priming[1], sep=''),
                                 REV_PRIMERS_150_df$revcomp_P7[1])
  
  cat("REV_PRIMERS_150_df\n")
  cat(str(REV_PRIMERS_150_df))
  cat("\n")
  
  indx.int<-c(which(colnames(REV_PRIMERS_150_df) == "seq_name"),which(colnames(REV_PRIMERS_150_df) == "sequence"))
  
  REV_PRIMERS_150_df_subset<-unique(REV_PRIMERS_150_df[,indx.int])
  
  cat("REV_PRIMERS_150_df_subset\n")
  cat(str(REV_PRIMERS_150_df_subset))
  cat("\n")
  
  LP_150_df<-rbind(LP_150_df,REV_PRIMERS_150_df_subset)
  
  cat("LP_150_df\n")
  cat(str(LP_150_df))
  cat("\n")
  
  
  #### LIBRARY PREPARATION 25 PCR #1 including sample indexes ----
  
  fastaFile <- read.fasta(file=opt$indexes_fasta)
  
  seq_name = fastaFile[seq(1, length(fastaFile), by=2)]
  sequence = fastaFile[seq(2, length(fastaFile), by=2)]
  
  FASTA.indexes <- data.frame(seq_name, sequence, stringsAsFactors = F)
  
  FASTA.indexes$seq_name<-gsub("^>","",FASTA.indexes$seq_name)
  
  cat("FASTA.indexes\n")
  cat(str(FASTA.indexes))
  cat("\n")
  
  
  
  LP_25_df<-data.frame(matrix(vector(), 0, 2,
                               dimnames=list(c(),c("seq_name","sequence"))),
                        stringsAsFactors=F)
  
  FWD_PRIMERS_25_df<-data.frame(matrix(vector(), dim(FASTA.indexes)[1], 4,
                                        dimnames=list(c(),c("seq_name","P5","index","Priming"))),
                                 stringsAsFactors=F)
  
  cat("FWD_PRIMERS_25_0_df\n")
  cat(str(FWD_PRIMERS_25_df))
  cat("\n")
  
  
  #sapply(FASTA.indexes$seq_name, paste())
  
  FWD_PRIMERS_25_df$seq_name<-paste("PE25_FWD",FASTA.indexes$seq_name, sep='_')
  cat("FWD_PRIMERS_25_0_df\n")
  cat(str(FWD_PRIMERS_25_df))
  cat("\n")
  
  FWD_PRIMERS_25_df$index<-FASTA.indexes$sequence
  cat("FWD_PRIMERS_25_1_df\n")
  cat(str(FWD_PRIMERS_25_df))
  cat("\n")
  FWD_PRIMERS_25_df$P5<-toupper("AATGATACGGCGACCACCGAGATCTACAC")
  cat("FWD_PRIMERS_25_2_df\n")
  cat(str(FWD_PRIMERS_25_df))
  cat("\n")
  FWD_PRIMERS_25_df$Priming<-toupper("GCAAAGTGAACACATCGCTAAGCGAAAGCTAAG")
  cat("FWD_PRIMERS_25_3_df\n")
  cat(str(FWD_PRIMERS_25_df))
  cat("\n")
  
  
  FWD_PRIMERS_25_df$sequence<-paste(FWD_PRIMERS_25_df$P5,
                                   FWD_PRIMERS_25_df$index,
                                   FWD_PRIMERS_25_df$Priming,
                                   sep='')
  
  cat("FWD_PRIMERS_25_df\n")
  cat(str(FWD_PRIMERS_25_df))
  cat("\n")
  
  indx.int<-c(which(colnames(FWD_PRIMERS_25_df) == "seq_name"),which(colnames(FWD_PRIMERS_25_df) == "sequence"))
  
  FWD_PRIMERS_25_df_subset<-unique(FWD_PRIMERS_25_df[,indx.int])
  
  cat("FWD_PRIMERS_25_df_subset\n")
  cat(str(FWD_PRIMERS_25_df_subset))
  cat("\n")
  
  LP_25_df<-rbind(LP_25_df,FWD_PRIMERS_25_df_subset)
  
  cat("LP_25_df\n")
  cat(str(LP_25_df))
  cat("\n")
  
  REV_PRIMERS_25_df<-data.frame(matrix(vector(), 1, 4,
                                        dimnames=list(c(),c("seq_name","P7","UMI","Priming"))),
                                 stringsAsFactors=F)
  
  REV_PRIMERS_25_df$Priming<-toupper("TCGACGAATTCGGCCGGCCGCTTCG")
  REV_PRIMERS_25_df$P7<-toupper("ATCTCGTATGCCGTCTTCTGCTTG")
  REV_PRIMERS_25_df$UMI<-toupper(paste(rep("N",10), collapse=''))
  
  REV_PRIMERS_25_df$seq_name<-c("LP_PE25_REV")
  
  REV_PRIMERS_25_df$revcomp_Priming<-sapply(REV_PRIMERS_25_df$Priming, revComp)
  REV_PRIMERS_25_df$revcomp_P7<-sapply(REV_PRIMERS_25_df$P7, revComp)
  
  
  
  REV_PRIMERS_25_df$sequence<-paste(REV_PRIMERS_25_df$revcomp_P7,REV_PRIMERS_25_df$UMI,REV_PRIMERS_25_df$revcomp_Priming, sep='')
  
  cat("REV_PRIMERS_25_df\n")
  cat(str(REV_PRIMERS_25_df))
  cat("\n")
  
  indx.int<-c(which(colnames(REV_PRIMERS_25_df) == "seq_name"),which(colnames(REV_PRIMERS_25_df) == "sequence"))
  
  REV_PRIMERS_25_df_subset<-unique(REV_PRIMERS_25_df[,indx.int])
  
  cat("REV_PRIMERS_25_df_subset\n")
  cat(str(REV_PRIMERS_25_df_subset))
  cat("\n")
  
  LP_25_df<-rbind(LP_25_df,REV_PRIMERS_25_df_subset)
  
  cat("LP_25_df\n")
  cat(str(LP_25_df))
  cat("\n")
  
  #### SEQUENCING ----
  
  
  SEQ_150_df<-data.frame(matrix(vector(), 0, 2,
                              dimnames=list(c(),c("seq_name","sequence"))),
                       stringsAsFactors=F)
  
  r1_150_df<-data.frame(matrix(vector(), dim(PCR_arms_df)[1], 4,
                                       dimnames=list(c(),c("seq_name","KpnI","BamHI","UP_FWD_Primer"))),
                                stringsAsFactors=F)
  
  r1_150_df$seq_name<-paste("PE_150_r1",PCR_arms_df$factor4,sep='_')
    r1_150_df$UP_FWD_Primer<-PCR_arms_df$UP_FWD_Primer
  r1_150_df$KpnI<-PCR_arms_df$KpnI
  r1_150_df$BamHI<-PCR_arms_df$BamHI
  r1_150_df$sequence<-paste(r1_150_df$BamHI,
                                   r1_150_df$KpnI,
                                   r1_150_df$UP_FWD_Primer,
                                   sep='')
  
  cat("r1_150_0_df\n")
  cat(str(r1_150_df))
  cat("\n")
  
  indx.int<-c(which(colnames(r1_150_df) == "seq_name"),which(colnames(r1_150_df) == "sequence"))
  
  r1_150_df_subset<-unique(r1_150_df[,indx.int])
  
  cat("r1_150_df_subset\n")
  cat(str(r1_150_df_subset))
  cat("\n")
  
  
  SEQ_150_df<-rbind(SEQ_150_df,r1_150_df_subset)
  
  cat("SEQ_150_df\n")
  cat(str(SEQ_150_df))
  cat("\n")
  
  index1_150_df<-r1_150_df
  
  index1_150_df$seq_name<-gsub("PE_150_r1","PE_150_index1",index1_150_df$seq_name)
  
  colnames(index1_150_df)[which(colnames(index1_150_df) == "sequence")]<-"FWD_sequence"
  
  
  index1_150_df$sequence<-sapply(index1_150_df$FWD_sequence, revComp)
  
  
  cat("index1_150_0_df\n")
  cat(str(index1_150_df))
  cat("\n")
  
  indx.int<-c(which(colnames(index1_150_df) == "seq_name"),which(colnames(index1_150_df) == "sequence"))
  
  index1_150_df_subset<-unique(index1_150_df[,indx.int])
  
  cat("index1_150_df_subset\n")
  cat(str(index1_150_df_subset))
  cat("\n")
  
  
  SEQ_150_df<-rbind(SEQ_150_df,index1_150_df_subset)
  
  cat("SEQ_150_df\n")
  cat(str(SEQ_150_df))
  cat("\n")
  
  r2_150_df<-data.frame(matrix(vector(), length(PCR_arms_df$factor4), 3,
                                      dimnames=list(c(),c("seq_name","UP_REV_Primer_revComp","Priming_REV_1"))),
                               stringsAsFactors=F)
  
  r2_150_df$UP_REV_Primer_revComp<-PCR_arms_df$UP_REV_Primer_revComp
  r2_150_df$Priming_REV_1<-toupper("CTAGTTG")
  r2_150_df$seq_name<-paste("PE_150_r2",PCR_arms_df$factor4,sep='_')
  r2_150_df$sequence<-paste(r2_150_df$Priming_REV_1,
                                   r2_150_df$UP_REV_Primer_revComp,
                                   sep='')
  
  
  cat("r2_150_df\n")
  cat(str(r2_150_df))
  cat("\n")
  
  
  indx.int<-c(which(colnames(r2_150_df) == "seq_name"),which(colnames(r2_150_df) == "sequence"))
  
  r2_150_df_subset<-unique(r2_150_df[,indx.int])
  
  cat("r2_150_df_subset\n")
  cat(str(r2_150_df_subset))
  cat("\n")
  
  
  SEQ_150_df<-rbind(SEQ_150_df,r2_150_df_subset)
  
  cat("SEQ_150_df\n")
  cat(str(SEQ_150_df))
  cat("\n")
  
  SEQ_25_df<-data.frame(matrix(vector(), 4, 2,
                                dimnames=list(c(),c("seq_name","sequence"))),
                         stringsAsFactors=F)
  
  # "GCAAAGTGAACACATCGCTAAGCGAAAGCTAAG" #r1
  # 
  # "TGAACCGATCGACGAATTCGGCCGGCCGCTTCG" #UMI
  # 
  # "GCCGGCCGAATTCGTCGATCGGTTCACGCAATG" #r2
  
  SEQ_25_df$seq_name<-c("PE25_r1","PE_25_index","PE_25_r2","PE_25_UMI")
  SEQ_25_df$sequence<-c(toupper("GCAAAGTGAACACATCGCTAAGCGAAAGCTAAG"),
                        revComp(toupper("GCAAAGTGAACACATCGCTAAGCGAAAGCTAAG")),
                        toupper("AATTCGTCGATCGGTTCACGCAATGGGATCC"),
                        toupper("TGAACCGATCGACGAATTCGGCCGGCCGCTTCG"))
  
  
  cat("SEQ_25_df\n")
  cat(str(SEQ_25_df))
  cat("\n")
  
  
  #### OTHER_primers ----
  
  OTHER_primers_df<-data.frame(matrix(vector(), 4, 2,
                               dimnames=list(c(),c("seq_name","sequence"))),
                        stringsAsFactors=F)
  
  # "GCAAAGTGAACACATCGCTAAGCGAAAGCTAAG" #r1
  # 
  # "TGAACCGATCGACGAATTCGGCCGGCCGCTTCG" #UMI
  # 
  # "GCCGGCCGAATTCGTCGATCGGTTCACGCAATG" #r2
  
  OTHER_primers_df$seq_name<-c("Gibson FWD","Gibson REV")
  OTHER_primers_df$sequence<-c(toupper("TCGACGAATTCGGCCGGCAG"),
                        toupper("TATCACCGGTGCATGCTCTA"))
  
  
  cat("OTHER_primers_df\n")
  cat(str(OTHER_primers_df))
  cat("\n")
  
  # >Gibson FWD
  # TCGACGAATTCGGCCGGCAG
  # >Gibson REV
  # TATCACCGGTGCATGCTCTA
  
  
  ##### Global_primers_fasta -----
  
  GLOBAL<-rbind(LA_df,LP_150_df,LP_25_df,SEQ_150_df,SEQ_25_df,OTHER_primers_df)
  
  cat("GLOBAL\n")
  cat(str(GLOBAL))
  cat("\n")
  
  GLOBAL$nchar<-nchar(GLOBAL$sequence)
  
  GLOBAL$seq_name<-paste(GLOBAL$seq_name,GLOBAL$nchar,sep='_')
  
  cat("GLOBAL\n")
  cat(str(GLOBAL))
  cat("\n")
  
  #### SAVE FASTA ----
  
  setwd(out)
  
  indx.fasta<-c(which(colnames(GLOBAL) == "seq_name"),which(colnames(GLOBAL) == "sequence"))
  
  Merge2_TO_FASTA<-unique(GLOBAL[,indx.fasta])
  
  df.fasta = dataframe2fas(Merge2_TO_FASTA, file="df.fasta")
  
  # Here we print the fasta file
  
  writeLines(df.fasta, sep ="\n", paste(type,'_Primers','.fasta',sep=''))
}


Function_7_printing_primers_second_part = function(option_list)
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
  
  #### indexes 96 ----
  
  indexes_96<-as.data.frame(fread(file=opt$indexes_96, sep="\t", header=T), stringsAsFactors = F)
  
  # cat("indexes_96\n")
  # cat(str(indexes_96))
  # cat("\n")
  # 
  colnames(indexes_96)<-c("ID","sequence")
  
  cat("indexes_96\n")
  cat(str(indexes_96))
  cat("\n")
  
  # quit(status = 1)
  
  
  #### LIBRARY PREPARATION 25 PCR #1 including sample indexes ----
  
  fastaFile <- read.fasta(file=opt$indexes_fasta)
  
  seq_name = fastaFile[seq(1, length(fastaFile), by=2)]
  sequence = fastaFile[seq(2, length(fastaFile), by=2)]
  
  FASTA.indexes <- data.frame(seq_name, sequence, stringsAsFactors = F)
  
  FASTA.indexes$seq_name<-gsub("^>","",FASTA.indexes$seq_name)
  
  cat("FASTA.indexes\n")
  cat(str(FASTA.indexes))
  cat("\n")
  
  
  #### indexes_96 minus indexes already used ----
  
  indx.dep<-which(indexes_96$i5_sequence%in%FASTA.indexes$sequence)
  
  if(length(indx.dep) >0)
  {
    indexes_96_RMV<-indexes_96[-indx.dep,]
    
  }else{
    
    indexes_96_RMV<-indexes_96
    
  }
  
  cat("indexes_96_RMV\n")
  cat(str(indexes_96_RMV))
  cat("\n")
  
  vector_id<-seq(11,dim(indexes_96_RMV)[1]+10, by=1)
  
  cat("vector_id\n")
  cat(str(vector_id))
  cat("\n")
  
  indexes_96_RMV$seq_name<-paste(rep("i",length(vector_id)),vector_id , sep='')
  
  cat("indexes_96_RMV\n")
  cat(str(indexes_96_RMV))
  cat("\n")
  
  LP_25_df<-data.frame(matrix(vector(), 0, 2,
                              dimnames=list(c(),c("seq_name","sequence"))),
                       stringsAsFactors=F)
  
  FWD_PRIMERS_25_df<-data.frame(matrix(vector(), dim(indexes_96_RMV)[1], 4,
                                       dimnames=list(c(),c("seq_name","P5","index","Priming"))),
                                stringsAsFactors=F)
  
  cat("FWD_PRIMERS_25_0_df\n")
  cat(str(FWD_PRIMERS_25_df))
  cat("\n")
  
  #quit(status=1)
  
  
  #sapply(indexes_96_RMV$seq_name, paste())
  
  FWD_PRIMERS_25_df$seq_name<-paste("PE25_FWD",indexes_96_RMV$seq_name, sep='_')
  cat("FWD_PRIMERS_25_0_df\n")
  cat(str(FWD_PRIMERS_25_df))
  cat("\n")
  
  FWD_PRIMERS_25_df$index<-indexes_96_RMV$sequence
  cat("FWD_PRIMERS_25_1_df\n")
  cat(str(FWD_PRIMERS_25_df))
  cat("\n")
  FWD_PRIMERS_25_df$P5<-toupper("AATGATACGGCGACCACCGAGATCTACAC")
  cat("FWD_PRIMERS_25_2_df\n")
  cat(str(FWD_PRIMERS_25_df))
  cat("\n")
  FWD_PRIMERS_25_df$Priming<-toupper("GCAAAGTGAACACATCGCTAAGCGAAAGCTAAG")
  cat("FWD_PRIMERS_25_3_df\n")
  cat(str(FWD_PRIMERS_25_df))
  cat("\n")
  
  
  FWD_PRIMERS_25_df$sequence<-paste(FWD_PRIMERS_25_df$P5,
                                    FWD_PRIMERS_25_df$index,
                                    FWD_PRIMERS_25_df$Priming,
                                    sep='')
  
  cat("FWD_PRIMERS_25_df\n")
  cat(str(FWD_PRIMERS_25_df))
  cat("\n")
  
  indx.int<-c(which(colnames(FWD_PRIMERS_25_df) == "seq_name"),which(colnames(FWD_PRIMERS_25_df) == "sequence"))
  
  FWD_PRIMERS_25_df_subset<-unique(FWD_PRIMERS_25_df[,indx.int])
  
  cat("FWD_PRIMERS_25_df_subset\n")
  cat(str(FWD_PRIMERS_25_df_subset))
  cat("\n")
  
  LP_25_df<-rbind(LP_25_df,FWD_PRIMERS_25_df_subset)
  
  cat("LP_25_df\n")
  cat(str(LP_25_df))
  
  indx.int2<-c(which(colnames(FWD_PRIMERS_25_df) == "seq_name"),which(colnames(FWD_PRIMERS_25_df) == "index"))
  
  FWD_PRIMERS_25_df_subset2<-unique(FWD_PRIMERS_25_df[,indx.int2])
  
  cat("FWD_PRIMERS_25_df_subset2\n")
  cat(str(FWD_PRIMERS_25_df_subset2))
  cat("\n")
  
  
  #### SAVE FASTA ----
  
  setwd(out)
  
  #indx.fasta<-c(which(colnames(GLOBAL) == "seq_name"),which(colnames(GLOBAL) == "sequence"))
  
  Merge2_TO_FASTA<-LP_25_df
  
  df.fasta = dataframe2fas(Merge2_TO_FASTA, file="df.fasta")
  
  # Here we print the fasta file
  
  writeLines(df.fasta, sep ="\n", paste(type,'_New_i5_Primers','.fasta',sep=''))
  
  df.fasta = dataframe2fas(FWD_PRIMERS_25_df_subset2, file="df.fasta")
  
  # Here we print the fasta file
  
  writeLines(df.fasta, sep ="\n", paste(type,'_New_i5_indexes','.fasta',sep=''))
  
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
    make_option(c("--span"), type="numeric", default=NULL,
                metavar="FILE.txt",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--sunk_cost"), type="character", default=NULL,
                metavar="FILE.txt",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--RE_sites_unstring"), type="character", default=NULL,
                metavar="FILE.txt",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--PCR_arms_unstring"), type="character", default=NULL,
                metavar="FILE.txt",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Blacklisted_sequences"), type="character", default=NULL,
                metavar="FILE.txt",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Barcode_source"), type="character", default=NULL,
                metavar="FILE.txt",
                help="Path to tab-separated input file listing regions to analyze. Required."),
     make_option(c("--Vetted_barcodes"), type="character", default=NULL,
                metavar="FILE.txt",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--indexes_fasta"), type="character", default=NULL,
                metavar="FILE.txt",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--indexes_96"), type="character", default=NULL,
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
  
  
 # Function_3_generate_long_matrix(opt)
 # Function_4_add_constant_sequences_and_index(opt)
 # Function_6_add_index(opt)
 # Function_7_printing_primers(opt)
 
 Function_7_printing_primers_second_part(opt)
 
  
}


###########################################################################

system.time( main() )
