
suppressMessages(library("dplyr", lib.loc="/nfs/users/nfs_m/mt19/sOFTWARE/R_libs/"))
suppressMessages(library("data.table", lib.loc="/nfs/users/nfs_m/mt19/sOFTWARE/R_libs/"))
suppressMessages(library("ggplot2", lib.loc = "/nfs/users/nfs_m/mt19/sOFTWARE/R_libs"))
suppressMessages(library("optparse", lib.loc = "/nfs/users/nfs_m/mt19/sOFTWARE/R_libs"))
suppressMessages(library("plyr", lib.loc="/nfs/users/nfs_m/mt19/sOFTWARE/R_libs/"))
suppressMessages(library("Biostrings"))
suppressMessages(library("vcfR"))


opt = NULL

Function_3_generate_fasta_auxiliary_files = function(option_list)
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
  
  filename<-paste(type,"_MASTER_TILES",".tsv", sep='')
  
  DEF_TILES<-read.table(file=filename, sep="\t", header=T, stringsAsFactors = F)
  
  DEF_TILES$REAL_TILE<-paste(DEF_TILES$KEY,DEF_TILES$Tile, sep="__")
  
  cat("DEF_TILES_0\n")
  cat(str(DEF_TILES))
  cat("\n")
  
  indx.NCGR<-grep("NA",DEF_TILES$VAR)
  
  DEF_TILES_NCGR<-DEF_TILES[indx.NCGR,]
  
  
  cat("DEF_TILES_NCGR_0\n")
  cat(str(DEF_TILES_NCGR))
  cat("\n")
  
  #### Read master tiles ----
  
  
  filename1=paste("LIBRARY_CANDIDATES_",type,".tsv", sep='')

  LIBRARY_CANDIDATES<-read.table(file=filename1, sep="\t", header=T, stringsAsFactors = F)
  
  cat("LIBRARY_CANDIDATES_0\n")
  cat(str(LIBRARY_CANDIDATES))
  cat("\n")
  
  indx.int<-c(which(colnames(LIBRARY_CANDIDATES) == "Label"),which(colnames(LIBRARY_CANDIDATES) == "KEY"))
  
  LIBRARY_CANDIDATES_subset<-unique(LIBRARY_CANDIDATES[,indx.int])
  
  cat("LIBRARY_CANDIDATES_subset_0\n")
  cat(str(LIBRARY_CANDIDATES_subset))
  cat("\n")
  
  #### Reference_fasta ----
  
  filename<-paste(type,"_REF",".fasta", sep='')
  
  fastaFile<-readDNAStringSet(file=filename)
  
  seq_name = names(fastaFile)
  sequence = paste(fastaFile)
  REF_fasta <- data.frame(seq_name, sequence, stringsAsFactors = F)
  
  REF_fasta$chr<-gsub(":.+$","",REF_fasta$seq_name)
  REF_fasta$chr<-gsub("^Chr","chr",REF_fasta$chr)
  
  REF_fasta$start<-gsub("^[^:]+:","",REF_fasta$seq_name)
  REF_fasta$start<-gsub("-.+$","",REF_fasta$start)
  REF_fasta$start<-as.integer(REF_fasta$start)
  REF_fasta$start<-as.integer(REF_fasta$start+1)
  
  REF_fasta$end<-gsub("^[^:]+:[0-9]+-","",REF_fasta$seq_name)
  REF_fasta$end<-gsub("-.+$","",REF_fasta$end)
  REF_fasta$end<-as.integer(REF_fasta$end)
  
  REF_fasta$distance<-REF_fasta$end-REF_fasta$start+1
  REF_fasta$nchar_seq<-nchar(REF_fasta$sequence)
  
  
  cat("REF_fasta\n")
  cat(str(REF_fasta))
  cat("\n")
  
  #quit(status = 1)
  
  
  indx.REF<-c(which(colnames(DEF_TILES) == "REAL_TILE"),which(colnames(DEF_TILES) == "chr"),which(colnames(DEF_TILES) == "start"),which(colnames(DEF_TILES) == "end"))
  
  REF_TILES<-unique(DEF_TILES[,indx.REF])
  
  cat("REF_TILES_0\n")
  cat(str(REF_TILES))
  cat("\n")
  
  #### merge fasta ref and REF_TILES ----
  
  REF_TILES<-merge(REF_TILES,
                   REF_fasta,
                   by=c("chr","start","end"),
                   all=T)
  
  cat("REF_TILES_1\n")
  cat(str(REF_TILES))
  cat("\n")
  
  
  
  
  #### vcf and intervals file ----
  
  path_vcf<-paste(out,'vcf_file', sep='')
  
  if (file.exists(path_vcf)){
    
  } else {
    dir.create(file.path(path_vcf))
    
    
  }
  
  path_intervals<-paste(out,'intervals_file', sep='')
  
  if (file.exists(path_intervals)){
    
  } else {
    dir.create(file.path(path_intervals))
    
    
  }
  
  vector_files_vcf<-NULL
  vector_files_interval<-NULL
  
  Gather<- data.frame(matrix(vector(), 0, dim(DEF_TILES)[2]+1,
                             dimnames=list(c(),
                                           c(colnames(DEF_TILES),"REF_sequence"))),
                      stringsAsFactors=F)
  
  for(i in 1:dim(DEF_TILES)[1])
  {
    
    DEF_TILES_sel<-DEF_TILES[i,]
    
    # cat("DEF_TILES_sel\n")
    # str(DEF_TILES_sel)
    # cat("\n")
    
    # if(i%in%indx.NCGR)
    # {
    #   cat("DEF_TILES_sel\n")
    #   str(DEF_TILES_sel)
    #   cat("\n")
    #   
    # }
    
    
    REF_TILES_sel<-REF_TILES[which(REF_TILES$REAL_TILE%in%DEF_TILES_sel$REAL_TILE),]
    
    # cat("REF_TILES_sel\n")
    # str(REF_TILES_sel)
    # cat("\n")
    
   
    carried_variants_sel<-unlist(strsplit(DEF_TILES_sel$carried_variants, split="\\|"))
    
    # cat("carried_variants_sel:")
    # cat(sprintf(as.character(carried_variants_sel)))
    # cat("\n")
    
    
    list_vcf<-list()
   # list_NCGR<-list()
    
    for(k in 1:length(carried_variants_sel))
    {
      carried_variants_sel_k<-carried_variants_sel[k]
      
      # cat("carried_variants_sel_k:")
      # cat(sprintf(as.character(carried_variants_sel_k)))
      # cat("\n")
      
      chr<-gsub("_.+$","",carried_variants_sel_k)
      pos<-gsub("^[^_]+_","",carried_variants_sel_k)
      pos<-as.integer(gsub("_.+$","",pos))
      ref<-gsub("^[^_]+_[^_]+_","",carried_variants_sel_k)
      ref<-gsub("_.+$","",ref)
      alt<-gsub("^[^_]+_[^_]+_[^_]+_","",carried_variants_sel_k)
      
      if(ref == "NA")
      {
        # cat("A_0\n")
        # str(A)
        # cat("\n")
        
        # cat("REF_TILES_sel\n")
        # str(REF_TILES_sel)
        # cat("\n")
        
        REL_POS_VAR<-pos-REF_TILES_sel$start
        
        # cat("REL_POS_VAR:")
        # cat(sprintf(as.character(REL_POS_VAR)))
        # cat("\n")
        
        ref_NCGR<-substr(REF_TILES_sel$sequence, REL_POS_VAR,REL_POS_VAR)
        alt_NCGR<-"NA"
        
        
        # cat("ref_NCGR:")
        # cat(sprintf(as.character(ref_NCGR)))
        # cat("\n")
        
        if(ref_NCGR == "C")
        {
          alt_NCGR<-"T"
        }
        if(ref_NCGR == "T")
        {
          alt_NCGR<-"C"
        }
        if(ref_NCGR == "G")
        {
          alt_NCGR<-"A"
        }
        if(ref_NCGR == "A")
        {
          alt_NCGR<-"G"
        }
        
        ref<-ref_NCGR
        alt<-alt_NCGR
        carried_variants_sel_k<-paste(chr,pos,ref,alt, sep="_")
        
        # 
        # quit(status = 1)
        
      }
      
      A<-as.data.frame(cbind(chr,pos,carried_variants_sel_k,ref,alt), stringsAsFactors=F)
      A$QUAL<-'100'
      A$FILTER<-'PASS'
      A$INFO<-'.'
      
      # cat("A_FIN\n")
      # str(A)
      # cat("\n")
      
      list_vcf[[k]]<-A
      
      
    }# k carried_variants_sel
    
    vcf_df = unique(as.data.frame(data.table::rbindlist(list_vcf)))
    
    if(i%in%indx.NCGR)
    {
      substituted_NCGR<-DEF_TILES_sel
      VAR_sel<-unique(DEF_TILES_sel$VAR)
      
      chr<-gsub("_.+$","",VAR_sel)
      pos<-gsub("^[^_]+_","",VAR_sel)
      pos<-as.integer(gsub("_.+$","",pos))
      ref<-gsub("^[^_]+_[^_]+_","",VAR_sel)
      ref<-gsub("_.+$","",ref)
      alt<-gsub("^[^_]+_[^_]+_[^_]+_","",VAR_sel)
      
        
        REL_POS_VAR<-pos-REF_TILES_sel$start
        
        # cat("REL_POS_VAR:")
        # cat(sprintf(as.character(REL_POS_VAR)))
        # cat("\n")
        
        ref_NCGR<-substr(REF_TILES_sel$sequence, REL_POS_VAR,REL_POS_VAR)
        alt_NCGR<-"NA"
        
        
        # cat("ref_NCGR:")
        # cat(sprintf(as.character(ref_NCGR)))
        # cat("\n")
        
        if(ref_NCGR == "C")
        {
          alt_NCGR<-"T"
        }
        if(ref_NCGR == "T")
        {
          alt_NCGR<-"C"
        }
        if(ref_NCGR == "G")
        {
          alt_NCGR<-"A"
        }
        if(ref_NCGR == "A")
        {
          alt_NCGR<-"G"
        }
        
        ref<-ref_NCGR
        alt<-alt_NCGR
        VAR_sel_NEW<-paste(chr,pos,ref,alt, sep="_")
     
      
      substituted_NCGR$VAR<-VAR_sel_NEW
      substituted_NCGR$carried_variants<-paste(vcf_df$carried_variants, collapse="\\|")

      substituted_NCGR<-cbind(substituted_NCGR,REF_TILES_sel$sequence)
      colnames(substituted_NCGR)<-colnames(Gather)
      
      Gather<-rbind(Gather,substituted_NCGR)
      
      # cat("substituted_NCGR\n")
      # str(substituted_NCGR)
      # cat("\n")
      # 
      # cat("Gather\n")
      # str(Gather)
      # cat("\n")
      # # 
      # # cat("vcf_df\n")
      # # str(vcf_df)
      # # cat("\n")
      # # 
      # quit(status = 1)
      
    }else{
      
      DEF_TILES_sel<-cbind(DEF_TILES_sel,REF_TILES_sel$sequence)
      colnames(DEF_TILES_sel)<-colnames(Gather)
      
      Gather<-rbind(Gather,DEF_TILES_sel)
      
    }
    
    # cat("Gather\n")
    # str(Gather)
    # cat("\n")
    # 
    # quit(status = 1)
    
    vcf_df$chr<-gsub("^chr","Chr",vcf_df$chr)
    
    vcf_df_ordered<-vcf_df[order(vcf_df$pos, decreasing = F),]
    
    filename_vcf<-paste(paste(DEF_TILES_sel$KEY,DEF_TILES_sel$Tile,i,sep="__"),".vcf",sep='')
    
    vector_files_vcf[i]<-filename_vcf
    
    # cat("filename_vcf:")
    # cat(sprintf(as.character(filename_vcf)))
    # cat("\n")
    
    setwd(path_vcf)
    
    cat(c('##fileformat=VCFv4.2',
          '##FILTER=<ID=PASS,Description="All filters passed">',
          '##ALT=<ID=NON_REF,Description="Represents any possible alternative allele at this location">',
          '##FILTER=<ID=LowQual,Description="Low quality">',
          '##GVCFBlock0-1=minGQ=0(inclusive),maxGQ=1(exclusive)',
          '##GVCFBlock1-2=minGQ=1(inclusive),maxGQ=2(exclusive)',
          '##GVCFBlock10-11=minGQ=10(inclusive),maxGQ=11(exclusive)',
          '##GVCFBlock11-12=minGQ=11(inclusive),maxGQ=12(exclusive)',
          '##GVCFBlock12-13=minGQ=12(inclusive),maxGQ=13(exclusive)',
          '##GVCFBlock13-14=minGQ=13(inclusive),maxGQ=14(exclusive)',
          '##GVCFBlock14-15=minGQ=14(inclusive),maxGQ=15(exclusive)',
          '##GVCFBlock15-16=minGQ=15(inclusive),maxGQ=16(exclusive)',
          '##GVCFBlock16-17=minGQ=16(inclusive),maxGQ=17(exclusive)',
          '##GVCFBlock17-18=minGQ=17(inclusive),maxGQ=18(exclusive)',
          '##GVCFBlock18-19=minGQ=18(inclusive),maxGQ=19(exclusive)',
          '##GVCFBlock19-20=minGQ=19(inclusive),maxGQ=20(exclusive)',
          '##GVCFBlock2-3=minGQ=2(inclusive),maxGQ=3(exclusive)',
          '##GVCFBlock20-21=minGQ=20(inclusive),maxGQ=21(exclusive)',
          '##GVCFBlock21-22=minGQ=21(inclusive),maxGQ=22(exclusive)',
          '##GVCFBlock22-23=minGQ=22(inclusive),maxGQ=23(exclusive)',
          '##GVCFBlock23-24=minGQ=23(inclusive),maxGQ=24(exclusive)',
          '##GVCFBlock24-25=minGQ=24(inclusive),maxGQ=25(exclusive)',
          '##GVCFBlock25-26=minGQ=25(inclusive),maxGQ=26(exclusive)',
          '##GVCFBlock26-27=minGQ=26(inclusive),maxGQ=27(exclusive)',
          '##GVCFBlock27-28=minGQ=27(inclusive),maxGQ=28(exclusive)',
          '##GVCFBlock28-29=minGQ=28(inclusive),maxGQ=29(exclusive)',
          '##GVCFBlock29-30=minGQ=29(inclusive),maxGQ=30(exclusive)',
          '##GVCFBlock3-4=minGQ=3(inclusive),maxGQ=4(exclusive)',
          '##GVCFBlock30-31=minGQ=30(inclusive),maxGQ=31(exclusive)',
          '##GVCFBlock31-32=minGQ=31(inclusive),maxGQ=32(exclusive)',
          '##GVCFBlock32-33=minGQ=32(inclusive),maxGQ=33(exclusive)',
          '##GVCFBlock33-34=minGQ=33(inclusive),maxGQ=34(exclusive)',
          '##GVCFBlock34-35=minGQ=34(inclusive),maxGQ=35(exclusive)',
          '##GVCFBlock35-36=minGQ=35(inclusive),maxGQ=36(exclusive)',
          '##GVCFBlock36-37=minGQ=36(inclusive),maxGQ=37(exclusive)',
          '##GVCFBlock37-38=minGQ=37(inclusive),maxGQ=38(exclusive)',
          '##GVCFBlock38-39=minGQ=38(inclusive),maxGQ=39(exclusive)',
          '##GVCFBlock39-40=minGQ=39(inclusive),maxGQ=40(exclusive)',
          '##GVCFBlock4-5=minGQ=4(inclusive),maxGQ=5(exclusive)',
          '##GVCFBlock40-41=minGQ=40(inclusive),maxGQ=41(exclusive)',
          '##GVCFBlock41-42=minGQ=41(inclusive),maxGQ=42(exclusive)',
          '##reference=file:///lustre/scratch116/vr/projects/hipsci/vrpipe/refs/human/ncbi37/hs37d5.fa',
          '##bcftools_concatVersion=1.3+htslib-1.3',
          '##bcftools_concatCommand=concat -a -D -f /lustre/scratch116/vr/projects/hipsci/vrpipe/a/8/4/3e0e3e374e1265c870395af411fa3/3053477/1_bcftools_concat/merge_list.txt',
          '##bcftools_annotateVersion=1.3+htslib-1.3',
          '##bcftools_annotateCommand=annotate -c CHROM,POS,ID,REF,ALT -a /lustre/scratch116/vr/projects/hipsci/vrpipe/refs/human/ncbi37/resources/annots-rsIDs-dbSNPv138.2014-04-02.tab.gz',
          '##bcftools_viewVersion=1.3.1+htslib-1.3.2',
          '##bcftools_viewCommand=view -Oz -o /lustre/scratch116/vr/projects/hipsci/releases/data/wgs/gatk_calls/HPSI0114i-kolf_2/HPSI0114i-kolf_2.wgs.gatk.haplotype_caller.20170425.genotypes.vcf.gz',
          '##bcftools_viewVersion=1.6+htslib-1.6',
          '##bcftools_viewCommand=view -R genERA_VAR_POST_HF1_2_3_region_file_FROM_TO.txt /lustre/scratch115/teams/soranzo/projects/genIE_analysis/reference_files/HPSI0114i-kolf_2.wgs.gatk.haplotype_caller.20170425.genotypes.vcf.gz; Date=Tue Oct 1 17:20:56 2019',
          paste('#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO', sep="\t")),
        file=filename_vcf,sep="\n")
    write.table(vcf_df_ordered,
                file=filename_vcf, sep="\t",append = T,
                quote=F,col.names = F, row.names = F, eol="\n")
    
    
    # TEST<-read.vcfR(
    #   filename_vcf)
    # 
    # cat("TEST\n")
    # str(TEST)
    # cat("\n")
    
    
    
   # write.vcf(vcf_df, file = filename_vcf, mask = FALSE, APPEND = FALSE)
    
    
    #quit(status = 1)
    
    
    ### interval
    
    filename_interval<-paste(paste(DEF_TILES_sel$KEY,DEF_TILES_sel$Tile,i,sep="__"),".intervals",sep='')
    
    vector_files_interval[i]<-filename_interval
    
    # cat("filename_interval:")
    # cat(sprintf(as.character(filename_interval)))
    # cat("\n")
    
    setwd(path_intervals)
    
    
    cat(paste(gsub("^chr","Chr",DEF_TILES_sel$chr),
              paste(DEF_TILES_sel$start,DEF_TILES_sel$end,sep="-"),
              sep=":"),file=filename_interval,sep="\n")
    
    # cat(paste("CHROM","POS","POS_TO", collapse="\t"),
    #     file=filename_interval,sep="\n")
    # 
    # cat(paste(gsub("^chr","Chr",DEF_TILES_sel$chr),
    #           DEF_TILES_sel$start,DEF_TILES_sel$end,collapse="\t"),
    #     file=filename_interval,sep="\n", append = T)
    # 
    
    
    
    # if(length(carried_variants_sel) >1)
    # {
    #   
    #   quit(status = 1)
    # }
    
    #Gather$REF_sequence<-as.character(Gather$REF_sequence)
    
    
  }#i DEF_TILES
  
  
  
  #### Add the two file columns to the master table
  
  Gather$vcf_file<-vector_files_vcf
  Gather$intervals_file<-vector_files_interval
  
  Gather<-merge(Gather,
                LIBRARY_CANDIDATES_subset,
                by="KEY",
                all=T)
  
  
  #### SAVE MASTER TABLE ----
  
  setwd(out)
  
  write.table(Gather, file=paste(type,"_MASTER_TILES_PLUS_REF_AND_NCGR_SUBS",".tsv", sep=''), sep="\t", row.names = F, quote=F)
  
  
  
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
  
  
 Function_3_generate_fasta_auxiliary_files(opt)
  
}


###########################################################################

system.time( main() )
