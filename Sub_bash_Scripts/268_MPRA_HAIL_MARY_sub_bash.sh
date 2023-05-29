#!/usr/bin/env bash
   
 
MASTER_ROUTE=$1
mem=$2
pc=$3
queue=$4
lane=$5
Bc_per_tile=$6
Threshold_reads_per_bc=$7
DEF_type=$8
base_quality_threshold=$(echo "6")
base_quality_threshold=$(echo "5")

umi_tools="/nfs/team151/software/umi-tools/bin/umi_tools"


source /software/hgi/installs/anaconda3/etc/profile.d/conda.sh


conda_bamtools=$(echo "/nfs/team151/software/bamtools/")

conda deactivate

conda activate $conda_bamtools


output_dir=$(echo "$MASTER_ROUTE")



 
type=$(echo "converting_and_filtering")
outfile_converting=$(echo "$output_dir""outfile""_""$type""_""$lane"".out")


touch $outfile_converting
echo -n "" > $outfile_converting



	
#### FILTER 1 base quality (bq) 

type=$(echo "filtering_1")
name_filtering_1=$(echo "$type""_""$lane""_job")


output_bwa_bam_sorted=$(echo "$output_dir""$lane""_sorted.bam")

echo "$output_bwa_bam_sorted"

output_bwa_bam_sorted_filtered=$(echo "$output_dir""$lane""_filtered.bam")

bsub -G team151 -o $outfile_converting -q normal -n$pc -J $name_filtering_1 -M $mem -R"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]" -- \
"samtools view -bq $base_quality_threshold $output_bwa_bam_sorted > $output_bwa_bam_sorted_filtered"

#### BLACKLIST EXTRACTION READ IDs unfiltered file
 
 
type=$(echo "BLACKLIST_EXTRACTION_pre_filter")
BLACKLIST_EXTRACTION_pre_filter=$(echo "$output_dir""outfile""_""$type""_""$lane"".out")
touch $BLACKLIST_EXTRACTION_pre_filter
echo -n "" > $BLACKLIST_EXTRACTION_pre_filter
name_BLACKLIST_EXTRACTION_pre_filter=$(echo "$type""_""$lane""_job")

output_readID_pre_filter=$(echo "$output_dir""$lane""_sorted_readID.txt")

bsub -G team151 -o $BLACKLIST_EXTRACTION_pre_filter -q normal -n$pc -J $name_BLACKLIST_EXTRACTION_pre_filter -M $mem -R"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]" -- \
"samtools view $output_bwa_bam_sorted|cut -f1 > $output_readID_pre_filter"

#### FILTER 2 NM

type=$(echo "filtering_2")
name_filtering_2=$(echo "$type""_""$lane""_job")

output_bwa_bam_filtered_NM_filtered=$(echo "$output_dir""$lane""_filtered_NM_filtered.bam")

conda activate $conda_bamtools

bsub -G team151 -o $outfile_converting -q normal -n$pc -w"done($name_filtering_1)"  -J $name_filtering_2 -M $mem -R"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]" -- \
"bamtools filter -tag NM:0 -in $output_bwa_bam_sorted_filtered -out $output_bwa_bam_filtered_NM_filtered"


#### LOOP To perform both things in filtering

Filter_files=$(echo "$output_bwa_bam_sorted_filtered"';'"$output_bwa_bam_filtered_NM_filtered"';')
echo $Filter_files
declare -a b=($(echo "$Filter_files" | tr ";" '\n'))

arraylength=${#b[@]}

Filter_options=$(echo "NoNM"';'"YesNM"';')
echo $Filter_options
declare -a c=($(echo "$Filter_options" | tr ";" '\n'))


echo "array_length: $arraylength"


for (( i=0; i<${arraylength}; i++ ));
do



      echo "index: $i, value: ${c[$i]}"

      
      Filter_file=${b[$i]}
      echo "Filter option: $Filter_file"

       Filter_option=${c[$i]}
       echo "Filter option: $Filter_option"


       DEF_type_2=$(echo "$Filter_option""_MAPQ_""$base_quality_threshold""_""$DEF_type")
       echo "------------------------------------------------->""$DEF_type_2"

       #### 1 sorting filtered results

       type=$(echo "sorting")
       name_sorting=$(echo "$type""_""$DEF_type_2""_job")

       output_bwa_bam_filtered_sorted=$(echo "$output_dir""$lane""_filtered_sorted_""$DEF_type_2"".bam")


       bsub -G team151 -o $outfile_converting -q normal -n$pc -w"done($name_filtering_2)" -J $name_sorting -M $mem -R"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]" -- \
       "samtools sort -@$pc $Filter_file -o $output_bwa_bam_filtered_sorted"

       #### 1 indexing filtered results
       
       type=$(echo "indexing")
       name_indexing=$(echo "$type""_""$DEF_type_2""_job")


       bsub -G team151 -o $outfile_converting -q normal -n$pc -w"done($name_sorting)" -J $name_indexing -M $mem -R"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]" -- \
       "samtools index $output_bwa_bam_filtered_sorted"
 

       #### UMI TOOLS GROUPING

       
       type=$(echo "umi_tools_grouping")
       outfile_umi_tools_grouping=$(echo "$output_dir""outfile""_""$type""_""$DEF_type_2"".out")
       touch $outfile_umi_tools_grouping
       echo -n "" > $outfile_umi_tools_grouping
       name_umi_tools_grouping=$(echo "$type""_""$DEF_type_2""_job")



       output_umi_tools_group=$(echo "$output_dir""$DEF_type_2""_group.tsv")
       output_umi_tools_log=$(echo "$output_dir""$DEF_type_2""_group.log")

       type=$(echo "BASHING")
       outfile_BASHING=$(echo "$output_dir""outfile""_""$type""_""$DEF_type_2"".out")
       touch $outfile_BASHING
       echo -n "" > $outfile_BASHING


       bsub -G team151 -o $outfile_umi_tools_grouping -M $mem -J $name_umi_tools_grouping -w"done($name_indexing)" -R"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]" -n$pc -q $queue -- \
       "$umi_tools group -I $output_bwa_bam_filtered_sorted --group-out=$output_umi_tools_group --log=$output_umi_tools_log"

       type=$(echo "UMI_sorting")
       name_UMI_sorting=$(echo  "$type""_""$DEF_type_2""_job")
       output_umi_summary=$(echo "$output_dir""$DEF_type_2""_UMI.tsv")


       bsub -G team151 -o $outfile_BASHING -M $mem  -w"done($name_umi_tools_grouping)" -J $name_UMI_sorting  -R"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]" -n$pc -q $queue -- \
     	    "cut -f2,7 $output_umi_tools_group|sort --parallel=$pc -S 20G > $output_umi_summary"


       type=$(echo "UMI_unique")
       name_UMI_unique=$(echo "$type""_""$DEF_type_2""_job")
       output_umi_summary_unique=$(echo "$output_dir""$DEF_type_2""_UMI_unique.tsv")

       bsub -G team151 -o $outfile_BASHING -M $mem  -w"done($name_UMI_sorting)" -J $name_UMI_unique  -R"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]" -n$pc -q $queue -- \
      "uniq -c $output_umi_summary > $output_umi_summary_unique"

          
       type=$(echo "Deduplicated_sorting")
       name_Deduplicated_sorting=$(echo  "$type""_""$DEF_type_2""_job")
       output_umi_tools_deduplicated_sorted=$(echo "$output_dir""$DEF_type_2""_deduplicated_sorted.tsv")


       bsub -G team151 -o $outfile_BASHING -M $mem  -w"done($name_UMI_unique)" -J $name_Deduplicated_sorting  -R"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]" -n$pc -q $queue -- \
       "cat $output_umi_summary_unique|tr -s ' '|cut -d ' ' -f 3|cut -f1|sort --parallel=$pc -S 20G > $output_umi_tools_deduplicated_sorted"


       type=$(echo "Deduplicated_unique")
       name_Deduplicated_unique=$(echo "$type""_""$DEF_type_2""_job")
       output_umi_tools_deduplicated_sorted_unique=$(echo "$output_dir""$DEF_type_2""_deduplicated_sorted_unique.tsv")

       bsub -G team151 -o $outfile_BASHING -M $mem  -w"done($name_Deduplicated_sorting)" -J $name_Deduplicated_unique  -R"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]" -n$pc -q $queue -- \
      "uniq -c $output_umi_tools_deduplicated_sorted > $output_umi_tools_deduplicated_sorted_unique"


       type=$(echo "Deduplicated_compress")
       name_Deduplicated_compress=$(echo "$type""_""$DEF_type_2""_job")
       output_umi_tools_group_collapsed=$(echo "$output_umi_tools_group_unique""_""$DEF_type_2"".gz")

       bsub -G team151 -o $outfile_BASHING -M $mem  -w"done($name_Deduplicated_unique)" -J $name_Deduplicated_compress  -R"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]" -n$pc -q $queue -- \
      "gzip -f $output_umi_tools_deduplicated_sorted_unique"



       type=$(echo "UMI_compress")
       name_UMI_compress=$(echo "$type""_""$DEF_type_2""_job")
       output_umi_summary_collapsed=$(echo "$output_umi_summary_unique""_""$DEF_type_2"".gz")

       bsub -G team151 -o $outfile_BASHING -M $mem  -w"done($name_Deduplicated_unique)" -J $name_UMI_compress  -R"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]" -n$pc -q $queue -- \
     	    "gzip -f $output_umi_summary_unique"

       #### BLACKLIST EXTRACTION READ IDs HERE HERE HERE HERE HERE

       type=$(echo "BLACKLIST_EXTRACTION")
       outfile_BLACKLIST_EXTRACTION=$(echo "$output_dir""outfile""_""$type""_""$DEF_type_2"".out")
       touch $outfile_BLACKLIST_EXTRACTION
       echo -n "" > $outfile_BLACKLIST_EXTRACTION
       name_BLACKLIST_EXTRACTION=$(echo "$type""_""$DEF_type_2""_job")
       output_readID_post_filter=$(echo "$output_dir""$DEF_type_2""_filtered_sorted_readID.txt")


       bsub -G team151 -o $outfile_BLACKLIST_EXTRACTION -q normal -n$pc -w"done($name_indexing)" -J $name_BLACKLIST_EXTRACTION -M $mem -R"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]" -- \
       "samtools view $output_bwa_bam_filtered_sorted|cut -f1 > $output_readID_post_filter"

       #### BLACKLIST R COMPARISON        ##### MEMORY ATTRITION STEP ----
       
       output_readID_pre_filter=$(echo "$output_dir""$lane""_sorted_readID.txt")
       output_readID_post_filter=$(echo "$output_dir""$DEF_type_2""_filtered_sorted_readID.txt")

       Rscript=/software/R-4.1.0/bin/Rscript
       Rscript_Blacklist=/nfs/users/nfs_m/mt19/Scripts/R/236_Read_PE150_Non_Programmed_BLACKLIST.R

       type=$(echo "Rscript_Blacklist")
       name_Rscript_Blacklist=$(echo "$type""_""$DEF_type_2""_job")
       Blacklist_type=$(echo "Blacklist""_""$DEF_type_2")

       Blacklist_file=$(echo "$output_dir""$Blacklist_type"".txt")

     #  echo "MEM:$mem:MEM"
       # echo "PC:$pc:PC"
       # echo "QUEUE:$queue:QUEUE"
       # echo "bsub -G team151 -o $outfile_BLACKLIST_EXTRACTION -M $mem -J $name_Rscript_Blacklist -R\"select[model==Intel_Platinum]\"  -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -n$pc -q $queue --"
       bsub -G team151 -o $outfile_BLACKLIST_EXTRACTION -M $mem -J $name_Rscript_Blacklist -w"done($name_BLACKLIST_EXTRACTION) && done($name_BLACKLIST_EXTRACTION_pre_filter)" -R"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]" -n$pc -q $queue -- \
       "$Rscript $Rscript_Blacklist \
        --TOTAL_ALIGNED_readIDs $output_readID_pre_filter \
        --FILTERED_ALIGNED_readIDs $output_readID_post_filter \
        --type $Blacklist_type --out $output_dir"

       #### REMOVE FILES from BLACKLIST EXTRACTION
       
       type=$(echo "Clear_ID_files")
       name_Clear_ID_files=$(echo "$type""_""$DEF_type_2""_job")

       bsub -G team151 -o $outfile_BLACKLIST_EXTRACTION -M $mem -J $name_Clear_ID_files -w"done($name_Rscript_Blacklist)" -R"select[model==Intel_Platinum]" -R"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]" -n$pc -q $queue -- \
     #   "rm -rf $output_readID_pre_filter $output_readID_post_filter"


       #### FINAL RDS OBJECT


       Rscript=/software/R-4.1.0/bin/Rscript
       Rscript_filtering_per_fraction=/nfs/users/nfs_m/mt19/Scripts/R/236_Read_PE150_Non_Programmed_MySeq.R

       type=$(echo "Rscript_part")
       name_Rscript_part=$(echo "$type""_""$DEF_type_2""_job")
       outfile_Rscript_part=$(echo "$output_dir""outfile""_""$type""_""$DEF_type_2"".out")
       touch $outfile_Rscript_part
       echo -n "" > $outfile_Rscript_part

       output_umi_tools_deduplicated_sorted_unique=$(echo "$DEF_type_2""_UMI_unique.tsv"".gz")
       output_umi_summary_collapsed=$(echo "$output_dir""$DEF_type_2""_UMI_unique.tsv"".gz")
       Blacklist_type=$(echo "Blacklist""_""$DEF_type_2")
       Blacklist_file=$(echo "$output_dir""$Blacklist_type"".txt")

       Bc_per_tile=$6
       Threshold_reads_per_bc=$7
       MPRA_Rosetta=$(echo "/lustre/scratch119/realdata/mdt2/teams/soranzo/users/mt19/NEW_TWIST_LIBS/NON_PROGRAMMED/NON_PROGRAMMED_MPRA_Rosetta.tsv")
       Threshold_min_length_barcode=$(echo "15")
       Threshold_min_bc_per_tile=$(echo "$Bc_per_tile")
       FINAL_type=$(echo "$DEF_type_2")
    


#       bsub -G team151 -o $outfile_Rscript_part -M $mem -J $name_Rscript_part -R"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]" -n$pc -q $queue -- \
       bsub -G team151 -o $outfile_Rscript_part -M $mem -J $name_Rscript_part -w"done($name_Rscript_Blacklist) && done($name_UMI_compress)" -R"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]" -n$pc -q $queue -- \
       "$Rscript $Rscript_filtering_per_fraction \
	--input_sequencing_seqname_collapsed $output_umi_summary_collapsed \
	--input_sequencing_Barcode_collapsed $output_umi_tools_deduplicated_sorted_unique \
	--MPRA_Rosetta $MPRA_Rosetta \
	--Blacklist $Blacklist_file \
	--Threshold_min_length_barcode $Threshold_min_length_barcode \
	--Threshold_reads_per_bc $Threshold_reads_per_bc \
	--Threshold_min_bc_per_tile $Threshold_min_bc_per_tile \
	--type $FINAL_type --out $output_dir"

       
done
