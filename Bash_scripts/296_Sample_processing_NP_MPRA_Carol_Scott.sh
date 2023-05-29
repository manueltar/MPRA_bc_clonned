#!/bin/bash>

 
#### Rscript

input_file=$1
MASTER_ROUTE=$2
batch=$3
mem=$4
pc=$5
queue=$6
input_file_2=$7
scratch_route=$8

REFERENCE_1=$(echo "/lustre/scratch119/realdata/mdt2/teams/soranzo/users/mt19/NEW_TWIST_LIBS/NON_PROGRAMMED/reference.fasta")
REFERENCE_2=$(echo "/lustre/scratch119/realdata/mdt2/teams/soranzo/users/mt19/NEW_TWIST_LIBS/NON_PROGRAMMED/SE_330_NovaSeq/reference.fasta")
LONG_MATRIX=$(echo "/lustre/scratch119/realdata/mdt2/teams/soranzo/users/mt19/NEW_TWIST_LIBS/NON_PROGRAMMED/NON_PROGRAMMED_LONG_matrix_seq_names_Plus_Flags.tsv")

GLOBAL_REFERENCE=$(echo "/lustre/scratch119/realdata/mdt2/teams/soranzo/users/mt19/NEW_TWIST_LIBS/NON_PROGRAMMED/combined_reference.fasta")




#### MASTER reference


# #### Rscript

# Rscript=/software/R-4.1.0/bin/Rscript


# output_dir=$(echo "/lustre/scratch119/realdata/mdt2/teams/soranzo/users/mt19/NEW_TWIST_LIBS/NON_PROGRAMMED/")

# Rscript_merge_references=/nfs/users/nfs_m/mt19/Scripts/R/249_0_NEW_REFERENCE_MERGER_v2.R

# type=$(echo "merge_references")
# outfile_merge_references=$(echo "$output_dir""outfile""_""$type"".out")
# touch $outfile_merge_references
# echo -n "" > $outfile_merge_references
# name_merge_references=$(echo "$type""_job")


# output=$(echo "/nfs/users/nfs_m/mt19/Scripts/Wraper_scripts/296_CREATION.sh")


# echo "bsub -G team151 -o $outfile_merge_references -M $mem  -J $name_merge_references -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -n$pc -q $queue -- \\" >> $output
# echo "\"$Rscript $Rscript_merge_references \\" >> $output
# echo "--REFERENCE_1 $REFERENCE_1 \\" >> $output
# echo "--REFERENCE_2 $REFERENCE_2 \\" >> $output
# echo "--LONG_MATRIX $LONG_MATRIX \\" >> $output
# echo "--type $type --out $output_dir\"" >> $output

# bash $output

# exit



# bwa align the reference




type_for_R=$(echo "100_percent")
declare -a arr


output_dir2=$scratch_route
echo "------------------------->$output_dir2"
output_dir=$MASTER_ROUTE

 
end_of_file=0
counter=0
array=()
array_Merge=()                  #



# while [[ $end_of_file == 0 ]]
# do
#   read -r line
#   end_of_file=$?

# #echo "LINE:$line"

#     if [[ "$line" =~ [^[:space:]] ]]; then

#     ((counter++))
#     a=($(echo "$line" | tr "\t" '\n'))
#     Lane_1=${a[0]}
#     Lane_2=${a[1]}
#     Merge=${a[2]}
  
#     echo "PROCESSED: $Lane_1 $Lane_2 > $Merge COUNTER $counter"

#     name_sample=$(echo "$Merge"|sed -r 's/_[RI][0-9]+\.fastq\.gz//g')
#     # name_job=$(echo "$batch""_""$name_sample""_""$counter")
#     # Type=$(echo "$name_sample"|sed -r 's/R[0-9]+_//g')
#     # Replicate=$(echo "$name_sample"|sed -r 's/_.+//g')

#     # echo "------------------------>name_job: $name_job"
#     echo "------------------------>$name_sample"
#     # echo "------------------------>$Lane_1 $Lane_2 $Merge"


#      type=$(echo "merging_Lanes""_""$Merge")
#      outfile_merging_Lanes=$(echo "$output_dir""outfile_""$type"".out")
#      touch $outfile_merging_Lanes
#      echo -n "" > $outfile_merging_Lanes
#      name_merging_Lanes=$(echo "$type""_job")
#      Sample_index_file=$(echo "/lustre/scratch119/realdata/mdt2/teams/soranzo/users/mt19/NEW_TWIST_LIBS/NON_PROGRAMMED/NovaSeq_Nov_2021_Demultiplex_Nov2021.txt")


    
#      bsub -G team151 -o $outfile_merging_Lanes -q $queue -n$pc -J $name_merging_Lanes -M $mem -R"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]" -- \
#      "cat $Lane_1 $Lane_2  > $Merge"



#     wait_done=$(echo "done(""$name_merging_Lanes"") &&")
#     array+=($wait_done)
#     array_Merge+=($Merge)

#     echo "$wait_done"


    
#     echo "---------------------------------------------------------------------------$Merge------------------------------------------>$counter\n"


#     if [ $counter == "4" ]; then

# 	#### CONDA

# 	source /software/hgi/installs/anaconda3/etc/profile.d/conda.sh


# 	conda_fastq_multx=$(echo "/nfs/users/nfs_m/mt19/sOFTWARE/fastq_multx")

# 	conda deactivate

# 	conda activate $conda_fastq_multx

# 	done_string=$(echo ""${array[@]}"""\"")

# 	echo "$done_string"

# 	trimmed_done_string=$(echo $done_string|sed -e 's/ \&\&\"$//g')
# #	trimmed_done_string=$(echo "\"""$trimmed_done_string""\"")
# 	echo "$trimmed_done_string"
# #	exit


# 	type=$(echo "demultiplexing_samples_R1_R2")
# 	outfile_demultiplexing_samples_R1_R2=$(echo "$output_dir""outfile_""$type"".out")
# 	touch $outfile_demultiplexing_samples_R1_R2
# 	echo -n "" > $outfile_demultiplexing_samples_R1_R2
# 	name_demultiplexing_samples_R1_R2=$(echo "$type""_job")
# 	Deconvolving_table=$(echo "/lustre/scratch119/realdata/mdt2/teams/soranzo/users/mt19/NEW_TWIST_LIBS/NON_PROGRAMMED/NovaSeq_Nov_2021_Demultiplex_Nov2021.txt")

# 	R1_file=$(echo "$output_dir""$name_sample""_R1.fastq.gz")
# 	R2_file=$(echo "$output_dir""$name_sample""_R2.fastq.gz")
# 	sample_index_file=$(echo "$output_dir""$name_sample""_I2.fastq.gz")
# 	UMI_file=$(echo "$output_dir""$name_sample""_I1.fastq.gz")

# 	output_R1_file=$(echo "$output_dir""R1"".%.fq.gz")
# 	output_R2_file=$(echo "$output_dir""R2"".%.fq.gz")
# 	output_sample_index_file=$(echo "$output_dir""I2"".%.fq.gz")

# 	echo "fastq-multx -m1 -d1 -x -B $Deconvolving_table $sample_index_file $R1_file $R2_file -o $output_sample_index_file $output_R1_file $output_R2_file"


#      	bsub -G team151 -o $outfile_demultiplexing_samples_R1_R2 -w"$trimmed_done_string" -q $queue -n$pc -J $name_demultiplexing_samples_R1_R2 -M $mem -R"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]" -- \
# 	"fastq-multx -m1 -d1 -x -B $Deconvolving_table $sample_index_file $R1_file $R2_file -o $output_sample_index_file $output_R1_file $output_R2_file"




# 	type=$(echo "demultiplexing_samples_UMI")
# 	outfile_demultiplexing_samples_UMI=$(echo "$output_dir""outfile_""$type"".out")
# 	touch $outfile_demultiplexing_samples_UMI
# 	echo -n "" > $outfile_demultiplexing_samples_UMI
# 	name_demultiplexing_samples_UMI=$(echo "$type""_job")
# 	Deconvolving_table=$(echo "/lustre/scratch119/realdata/mdt2/teams/soranzo/users/mt19/NEW_TWIST_LIBS/NON_PROGRAMMED/NovaSeq_Nov_2021_Demultiplex_Nov2021.txt")

# 	R1_file=$(echo "$output_dir""$name_sample""_R1.fastq.gz")
# 	R2_file=$(echo "$output_dir""$name_sample""_R2.fastq.gz")
# 	sample_index_file=$(echo "$output_dir""$name_sample""_I2.fastq.gz")
# 	UMI_file=$(echo "$output_dir""$name_sample""_I1.fastq.gz")

# 	output_UMI_file=$(echo "$output_dir""I1"".%.fq.gz")
# 	output_sample_index_file=$(echo "$output_dir""I2_NEW"".%.fq.gz")



# 	bsub -G team151 -o $outfile_demultiplexing_samples_UMI -w"$trimmed_done_string" -q $queue -n$pc -J $name_demultiplexing_samples_UMI -M $mem -R"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]" -- \
#         "fastq-multx -m1 -d1 -x -B $Deconvolving_table $sample_index_file $UMI_file -o $output_sample_index_file $output_UMI_file"

#         bsub -G team151 -o $outfile_demultiplexing_samples_UMI -w"$name_demultiplexing_samples_UMI" -q $queue -n$pc -J $name_demultiplexing_samples_UMI -M $mem -R"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]" -- \
# 	"rm $sample_index_file"

# #	exit
	
#        fi



#   fi #no spaces
# done < "$input_file"




######


	##### After this is finished:
	# extract UMI R1
	# extract UMI R2
	# merge
	# align against new reference



	
	


end_of_file=0

while [[ $end_of_file == 0 ]]
do
  read -r line
  end_of_file=$?

#echo "LINE:$line"

    if [[ "$line" =~ [^[:space:]] ]]; then

    ((counter++))
    a=($(echo "$line" | tr "\t" '\n'))
    R1_file=${a[0]}
    R2_file=${a[2]}
    UMI_file=${a[1]}
    Cell_Type=$(echo $R1_file|sed -e 's/^R1\.//g')
    Cell_Type=$(echo $Cell_Type|sed -e 's/_.*$//g')
    Rep=$(echo $R1_file|sed -e 's/^R1\.//g')
    Rep=$(echo $Rep|sed -e 's/^[^_]*_//g')
    Rep=$(echo $Rep|sed -e 's/_.*$//g')
    Type=$(echo $R1_file|sed -e 's/^R1\.//g')
    Type=$(echo $Type|sed -e 's/^[^_]*_[^_]*_//g')
    Type=$(echo $Type|sed -e 's/\..*$//g')
  

    echo "---------->R1:$R1_file R2:$R2_file UMI:$UMI_file<-"
    echo "---------->Cell_Type:$Cell_Type Rep:$Rep Type:$Type<-"
    master_prefix=$(echo "$Cell_Type""_""$Rep""_""$Type")
    echo "master_prefix: $master_prefix"

    ################################################################## UMI_LABELLING_R1 ####################################################################################################


    type=$(echo "UMI_LABELLING_R1""_""$master_prefix")
    outfile_UMI_LABELLING_R1=$(echo "$output_dir""outfile_""$type"".out")
    touch $outfile_UMI_LABELLING_R1
    echo -n "" > $outfile_UMI_LABELLING_R1
    name_UMI_LABELLING_R1=$(echo "$type""_job")


    output_processed_R1=$(echo "$output_dir""R1_UMI_added_""$Cell_Type""_""$Rep""_""$Type"".fq.gz")




    bsub -G team151 -o $outfile_UMI_LABELLING_R1 -q $queue -n$pc -J $name_UMI_LABELLING_R1 -M $mem -R"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]" -- \
    "/nfs/team151/software/umi-tools/bin/umi_tools extract --bc-pattern=NNNNNNNNNN --stdin=$UMI_file --read2-in=$R1_file --stdout=$output_processed_R1 --read2-stdout"

    ################################################################## UMI_LABELLING_R2 ####################################################################################################


    type=$(echo "UMI_LABELLING_R2""_""$master_prefix")
    outfile_UMI_LABELLING_R2=$(echo "$output_dir""outfile_""$type"".out")
    touch $outfile_UMI_LABELLING_R2
    echo -n "" > $outfile_UMI_LABELLING_R2
    name_UMI_LABELLING_R2=$(echo "$type""_job")


    output_processed_R2=$(echo "$output_dir""R2_UMI_added_""$Cell_Type""_""$Rep""_""$Type"".fq.gz")



    bsub -G team151 -o $outfile_UMI_LABELLING_R2 -q $queue -n$pc -J $name_UMI_LABELLING_R2 -M $mem -R"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]" -- \
    "/nfs/team151/software/umi-tools/bin/umi_tools extract --bc-pattern=NNNNNNNNNN --stdin=$UMI_file --read2-in=$R2_file --stdout=$output_processed_R2 --read2-stdout"

    ####################   flash_merge ####################################################################################################

    type=$(echo "flash_merge""_""$master_prefix")
    outfile_flash_merge=$(echo "$output_dir2""outfile_""$type"".out")
    touch $outfile_flash_merge
    echo -n "" > $outfile_flash_merge
    name_flash_merge=$(echo "$type""_job")

    master_prefix=$(echo "$Cell_Type""_""$Rep""_""$Type")
    echo "master_prefix: $master_prefix"
    output_processed_R1=$(echo "$output_dir""R1_UMI_added_""$Cell_Type""_""$Rep""_""$Type"".fq.gz")
    output_processed_R2=$(echo "$output_dir""R2_UMI_added_""$Cell_Type""_""$Rep""_""$Type"".fq.gz")

    

    
    flash2_output=$(echo "$output_dir2""$master_prefix"".extendedFrags.fastq.gz")

    echo "$flash2_output"
    flash2_not_combined_1=$(echo "$output_dir2""$master_prefix"".notCombined_1.fastq.gz")
    flash2_not_combined_2=$(echo "$output_dir2""$master_prefix"".notCombined_2.fastq.gz")


    #    bsub -G team151 -o $outfile_flash_merge -q $queue -n$pc -J $name_flash_merge -M $mem -R"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]" -- \
    bsub -G team151 -o $outfile_flash_merge -q $queue -n$pc -w"done($name_UMI_LABELLING_R1) && done($name_UMI_LABELLING_R2)" -J $name_flash_merge -M $mem -R"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]" -- \
    "flash2 -z --output-prefix=$master_prefix -t $pc -r 15 -f 15 -s 2 -m 12 -x 0.122 --output-directory $output_dir2 $output_processed_R1 $output_processed_R2"

    ################################################################## clean_1 ####################################################################################################


    type=$(echo "clean_1""_""$master_prefix")
    outfile_clean_1=$(echo "$output_dir2""outfile_""$type"".out")
    touch $outfile_clean_1
    echo -n "" > $outfile_clean_1
    name_clean_1=$(echo "$type""_job")

    master_prefix=$(echo "$Cell_Type""_""$Rep""_""$Type")
    echo "master_prefix: $master_prefix"
    output_processed_R1=$(echo "$output_dir""R1_UMI_added_""$Cell_Type""_""$Rep""_""$Type"".fq.gz")
    output_processed_R2=$(echo "$output_dir""R2_UMI_added_""$Cell_Type""_""$Rep""_""$Type"".fq.gz")

     
    
#    bsub -G team151 -o $outfile_clean_1 -q $queue -n$pc -J $name_clean_1 -M $mem -R"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]" -- \
    bsub -G team151 -o $outfile_clean_1 -q $queue -n$pc -w"done($name_flash_merge)" -J $name_clean_1 -M $mem -R"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]" -- \
    "rm $output_processed_R1 $output_processed_R2 $flash2_not_combined_1 $flash2_not_combined_2"

# #    bwa_mem ####################################################################################################


    type=$(echo "bwa_mem""_""$master_prefix")
    outfile_bwa_mem=$(echo "$output_dir2""outfile_""$type"".out")
    touch $outfile_bwa_mem
    echo -n "" > $outfile_bwa_mem
    name_bwa_mem=$(echo "$type""_job")

    master_prefix=$(echo "$Cell_Type""_""$Rep""_""$Type")
    echo "master_prefix: $master_prefix"

     flash2_output=$(echo "$output_dir2""$master_prefix"".extendedFrags.fastq.gz")

    sai_output=$(echo "$output_dir2""$master_prefix"".sai")
    output_processed_R1=$(echo "$output_dir2""R1_UMI_added_""$Cell_Type""_""$Rep""_""$Type"".fq.gz")
    output_processed_R2=$(echo "$output_dir2""R2_UMI_added_""$Cell_Type""_""$Rep""_""$Type"".fq.gz")

     



    bsub -G team151 -o $outfile_bwa_mem -q $queue -n$pc -w"done($name_flash_merge)" -J $name_bwa_mem -M $mem -R"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]" -- \
    "bwa aln -l 15 -O 100 -E 100 $GLOBAL_REFERENCE $flash2_output > $sai_output"
#    bsub -G team151 -o $outfile_bwa_mem -q $queue -n$pc -J $name_bwa_mem -M $mem -R"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]" -- \

#    bsub -G team151 -o $outfile_bwa_mem -q $queue -n$pc -J $name_bwa_mem -M $mem -R"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]" -- \
#    "bwa aln -l 15 -O 100 -E 100 $GLOBAL_REFERENCE $output_processed_R1 > $sai_output"
#    bsub -G team151 -o $outfile_bwa_mem -q $queue -n$pc -w"done($name_UMI_LABELLING_R2) && done($name_UMI_LABELLING_R1)" -J $name_bwa_mem -M $mem -R"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]" -- \
#    "bwa aln -l 15 -O 100 -E 100 $GLOBAL_REFERENCE $flash2_output > $sai_output"

    ################################################################## sai_to_sam ####################################################################################################


    type=$(echo "sai_to_sam""_""$master_prefix")
    outfile_sai_to_sam=$(echo "$output_dir2""outfile_""$type"".out")
    touch $outfile_sai_to_sam
    echo -n "" > $outfile_sai_to_sam
    name_sai_to_sam=$(echo "$type""_job")



    sai_output=$(echo "$output_dir2""$master_prefix"".sai")
    sam_output=$(echo "$output_dir2""$master_prefix"".sam")
    
#    bsub -G team151 -o $outfile_sai_to_sam -q $queue -n$pc -J $name_sai_to_sam -M $mem -R"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]" -- \
#    bsub -G team151 -o $outfile_sai_to_sam -q $queue -n$pc -w"done($name_bwa_mem)" -J $name_sai_to_sam -M $mem -R"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]" -- \
#     "bwa samse $GLOBAL_REFERENCE $sai_output $output_processed_R1 > $sam_output"

     bsub -G team151 -o $outfile_sai_to_sam -q $queue -n$pc -w"done($name_bwa_mem)" -J $name_sai_to_sam -M $mem -R"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]" -- \
    "bwa samse $GLOBAL_REFERENCE $sai_output $flash2_output > $sam_output"

     ################################################################## clean_2 ####################################################################################################


    type=$(echo "clean_2""_""$master_prefix")
    outfile_clean_2=$(echo "$output_dir2""outfile_""$type"".out")
    touch $outfile_clean_2
    echo -n "" > $outfile_clean_2
    name_clean_2=$(echo "$type""_job")


    bsub -G team151 -o $outfile_clean_2 -q $queue -n$pc -w"done($name_sai_to_sam)" -J $name_clean_2 -M $mem -R"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]" -- \
    "rm $sai_output"
    
    ################################################################## sam_to_bam_and_filter_first_alignment ####################################################################################################


    type=$(echo "sam_to_bam_and_filter_first_alignment""_""$master_prefix")
    outfile_sam_to_bam_and_filter_first_alignment=$(echo "$output_dir2""outfile_""$type"".out")
    touch $outfile_sam_to_bam_and_filter_first_alignment
    echo -n "" > $outfile_sam_to_bam_and_filter_first_alignment
    name_sam_to_bam_and_filter_first_alignment=$(echo "$type""_job")


    sam_output=$(echo "$output_dir2""$master_prefix"".sam")
    output_bwa_bam=$(echo "$output_dir2""$master_prefix"".bam")
    

    bsub -G team151 -o $outfile_sam_to_bam_and_filter_first_alignment -q $queue -n$pc -w"done($name_sai_to_sam)" -J $name_sam_to_bam_and_filter_first_alignment -M $mem -R"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]" -- \
    "samtools view -@$pc -h -F 20,4,256 $sam_output|samtools sort -@$pc -o $output_bwa_bam"


    ####### Get unalagined reads
    
    type=$(echo "unaligned_reads""_""$master_prefix")
    outfile_unaligned_reads=$(echo "$output_dir2""outfile_""$type"".out")
    touch $outfile_unaligned_reads
    echo -n "" > $outfile_unaligned_reads
    name_unaligned_reads=$(echo "$type""_job")


    sam_output=$(echo "$output_dir2""$master_prefix"".sam")
    unaligned_bwa_bam=$(echo "$output_dir2""Unaligned_reads_""$master_prefix"".bam")
    

    bsub -G team151 -o $outfile_unaligned_reads -q $queue -n$pc -w"done($name_sai_to_sam)" -J $name_unaligned_reads -M $mem -R"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]" -- \
    "samtools view -@$pc -h -f 20,4 -F 256 $sam_output|samtools sort -@$pc -o $unaligned_bwa_bam"

    ################################################################## clean_2 ####################################################################################################


    type=$(echo "clean_2_5""_""$master_prefix")
    name_clean_2=$(echo "$type""_job")


    bsub -G team151 -o $outfile_clean_2 -q $queue -n$pc -w"done($name_unaligned_reads)" -J $name_clean_2 -M $mem -R"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]" -- \
    "rm $sam_output"

    ################################################################## indexing_1 ####################################################################################################


    type=$(echo "indexing_1""_""$master_prefix")
    outfile_indexing_1=$(echo "$output_dir2""outfile_""$type"".out")
    touch $outfile_indexing_1
    echo -n "" > $outfile_indexing_1
    name_indexing_1=$(echo "$type""_job")


    sam_output=$(echo "$output_dir2""$master_prefix"".sam")
    output_bwa_bam=$(echo "$output_dir2""$master_prefix"".bam")
    

    bsub -G team151 -o $outfile_indexing_1 -q $queue -n$pc -w"done($name_sam_to_bam_and_filter_first_alignment)" -J $name_indexing_1 -M $mem -R"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]" -- \
    "samtools index $output_bwa_bam"


    ################################################################## bamtools_filter ####################################################################################################
    
    source /software/hgi/installs/anaconda3/etc/profile.d/conda.sh


    conda_bamtools=$(echo "/nfs/team151/software/bamtools/")

    conda deactivate

    conda activate $conda_bamtools


    type=$(echo "bamtools_filter""_""$master_prefix")
    outfile_bamtools_filter=$(echo "$output_dir2""outfile_""$type"".out")
    touch $outfile_bamtools_filter
    echo -n "" > $outfile_bamtools_filter
    name_bamtools_filter=$(echo "$type""_job")


    output_bwa_bam=$(echo "$output_dir2""$master_prefix"".bam")
    output_bwa_bam_filtered=$(echo "$output_dir2""$master_prefix""_filtered"".bam")

    


    #     bsub -G team151 -o $outfile_bamtools_filter -q $queue -n$pc -J $name_bamtools_filter -M $mem -R"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]" -- \
    bsub -G team151 -o $outfile_bamtools_filter -q $queue -n$pc -w"done($name_indexing_1)" -J $name_bamtools_filter -M $mem -R"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]" -- \
    "bamtools filter -tag NM:0 -in $output_bwa_bam -out $output_bwa_bam_filtered"

    ################################################################## indexing_2 ####################################################################################################


    type=$(echo "indexing_2""_""$master_prefix")
    outfile_indexing_2=$(echo "$output_dir2""outfile_""$type"".out")
    touch $outfile_indexing_2
    echo -n "" > $outfile_indexing_2
    name_indexing_2=$(echo "$type""_job")

    output_bwa_bam_filtered=$(echo "$output_dir2""$master_prefix""_filtered"".bam")
    

    bsub -G team151 -o $outfile_indexing_2 -q $queue -n$pc -w"done($name_bamtools_filter)" -J $name_indexing_2 -M $mem -R"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]" -- \
    "samtools index $output_bwa_bam_filtered"

    ###############################    /nfs/team151/software/umi-tools/bin/umi_tools_grouping ####################################################################################################

    type=$(echo "umi_tools_grouping""_""$master_prefix")
    outfile_umi_tools_grouping=$(echo "$output_dir2""outfile""_""$type"".out")
    touch $outfile_umi_tools_grouping
    echo -n "" > $outfile_umi_tools_grouping
    name_umi_tools_grouping=$(echo "$type""_job")


#    output_bwa_bam=$(echo "$output_dir2""$master_prefix"".bam")
    output_bwa_bam_filtered=$(echo "$output_dir2""$master_prefix""_filtered"".bam")
    output_umi_tools_group=$(echo "$output_dir2""$master_prefix""_group.tsv")
    output_umi_tools_log=$(echo "$output_dir2""$master_prefix""_group.log")
    output_bwa_bam=$(echo "$output_dir2""$master_prefix"".bam")


#    bsub -G team151 -o $outfile_umi_tools_grouping -M $mem -J $name_umi_tools_grouping -R"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]" -n$pc -q $queue -- \
    #    bsub -G team151 -o $outfile_umi_tools_grouping -M 32000 -J $name_umi_tools_grouping -w"done($name_indexing_1)" -R"select[mem>=32000] rusage[mem=32000] span[hosts=1]" -n8 -q $queue -- \

    bsub -G team151 -o $outfile_umi_tools_grouping -M 32000 -J $name_umi_tools_grouping -w"done($name_indexing_2)" -R"select[mem>=32000] rusage[mem=32000] span[hosts=1]" -n8 -q $queue -- \
    "/nfs/team151/software/umi-tools/bin/umi_tools group -I $output_bwa_bam_filtered --group-out=$output_umi_tools_group --log=$output_umi_tools_log"

#    bsub -G team151 -o $outfile_umi_tools_grouping -M 64000 -J $name_umi_tools_grouping -w"done($name_indexing_2)" -R"select[mem>=64000] rusage[mem=64000] span[hosts=1]" -n8 -q $queue -- \
#    bsub -G team151 -o $outfile_umi_tools_grouping -M 64000 -J $name_umi_tools_grouping -R"select[mem>=64000] rusage[mem=64000] span[hosts=1]" -n8 -q $queue -- \
#    "/nfs/team151/software/umi-tools/bin/umi_tools group -I $output_bwa_bam --group-out=$output_umi_tools_group --log=$output_umi_tools_log"

    ################################################################## bash_cut_and_sort ####################################################################################################

    type=$(echo "bash_cut_and_sort""_""$master_prefix")
    outfile_bash_cut_and_sort=$(echo "$output_dir2""outfile""_""$type"".out")
    touch $outfile_bash_cut_and_sort
    echo -n "" > $outfile_bash_cut_and_sort
    name_bash_cut_and_sort=$(echo "$type""_job")


    output_umi_tools_group=$(echo "$output_dir2""$master_prefix""_group.tsv")
    output_umi_summary=$(echo "$output_dir2""$master_prefix""_UMI.tsv")
    max_memory_G=$(expr $mem / 1000)
    max_memory_G=$(echo "$max_memory_G""G")

    echo "max_memory_G_1: $max_memory_G"


    bsub -G team151 -o $outfile_bash_cut_and_sort -M $mem -J $name_bash_cut_and_sort -w"done($name_umi_tools_grouping)" -R"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]" -n$pc -q $queue -- \
    "cut -f2,7 $output_umi_tools_group|sort --parallel=$pc -S $max_memory_G > $output_umi_summary"


    
    ################################################################## bash_unique ####################################################################################################

    type=$(echo "bash_unique""_""$master_prefix")
    outfile_bash_unique=$(echo "$output_dir2""outfile""_""$type"".out")
    touch $outfile_bash_unique
    echo -n "" > $outfile_bash_unique
    name_bash_unique=$(echo "$type""_job")


    output_umi_summary=$(echo "$output_dir2""$master_prefix""_UMI.tsv")
    output_umi_summary_unique=$(echo "$output_dir2""$master_prefix""_UMI_unique.tsv")
    

#    bsub -G team151 -o $outfile_bash_unique -M $mem -J $name_bash_unique -R"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]" -n$pc -q $queue -- \
    bsub -G team151 -o $outfile_bash_unique -M $mem -J $name_bash_unique -w"done($name_bash_cut_and_sort)" -R"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]" -n$pc -q $queue -- \
    "uniq -c $output_umi_summary > $output_umi_summary_unique"

    ################################################################## deduplicate ####################################################################################################

    type=$(echo "deduplicate""_""$master_prefix")
    outfile_deduplicate=$(echo "$output_dir2""outfile""_""$type"".out")
    touch $outfile_deduplicate
    echo -n "" > $outfile_deduplicate
    name_deduplicate=$(echo "$type""_job")

    output_umi_summary_unique=$(echo "$output_dir2""$master_prefix""_UMI_unique.tsv")
    
    output_umi_tools_deduplicated_sorted=$(echo "$output_dir2""$master_prefix""_deduplicated_sorted.tsv")

    max_memory_G=$(expr $mem / 1000)
    max_memory_G=$(echo "$max_memory_G""G")

    echo "max_memory_G_2: $max_memory_G"


    bsub -G team151 -o $outfile_deduplicate -M $mem -J $name_deduplicate -w"done($name_bash_unique)" -R"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]" -n$pc -q $queue -- \
    "cat $output_umi_summary_unique|tr -s ' '|cut -d ' ' -f 3|cut -f1|sort --parallel=$pc -S $max_memory_G > $output_umi_tools_deduplicated_sorted"

    ################################################################## COMPRESS ####################################################################################################

    type=$(echo "bash_unique_COMPRESS""_""$master_prefix")
    name_bash_unique_COMPRESS=$(echo "$type""_job")

    bsub -G team151 -o $outfile_bash_unique -q $queue -M $mem -J $name_bash_unique_COMPRESS -w"done($name_deduplicate)" -R"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]" -n$pc -- \
    "gzip -f $output_umi_summary_unique"



    ################################################################## bash_unique_DEDUPLICATED ####################################################################################################

    type=$(echo "bash_unique_DEDUPLICATED""_""$master_prefix")
    outfile_bash_unique_DEDUPLICATED=$(echo "$output_dir2""outfile""_""$type"".out")
    touch $outfile_bash_unique_DEDUPLICATED
    echo -n "" > $outfile_bash_unique_DEDUPLICATED
    name_bash_unique_DEDUPLICATED=$(echo "$type""_job")

    output_umi_tools_deduplicated_sorted=$(echo "$output_dir2""$master_prefix""_deduplicated_sorted.tsv")
    output_umi_tools_deduplicated_sorted_unique=$(echo "$output_dir2""$master_prefix""_deduplicated_sorted_unique.tsv")
    
    
    bsub -G team151 -o $outfile_bash_unique_DEDUPLICATED -M $mem -J $name_bash_unique_DEDUPLICATED -w"done($name_deduplicate)" -R"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]" -n$pc -q $queue -- \
    "uniq -c $output_umi_tools_deduplicated_sorted > $output_umi_tools_deduplicated_sorted_unique"

    type=$(echo "bash_unique_DEDUPLICATED_COMPRESS""_""$master_prefix")
    name_bash_unique_DEDUPLICATED_COMPRESS=$(echo "$type""_job")

    bsub -G team151 -o $outfile_bash_unique_DEDUPLICATED -M $mem -J $name_bash_unique_DEDUPLICATED_COMPRESS -w"done($name_bash_unique_DEDUPLICATED)" -R"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]" -n$pc -q $queue -- \
    "gzip -f $output_umi_tools_deduplicated_sorted_unique"

	 

    #############################     clean_3 ####################################################################################################

    type=$(echo "clean_3""_""$master_prefix")
    outfile_clean_3=$(echo "$output_dir2""outfile_""$type"".out")
    touch $outfile_clean_3
    echo -n "" > $outfile_clean_3
    name_clean_3=$(echo "$type""_job")


    bsub -G team151 -o $outfile_clean_3 -q $queue -n$pc -w"done($name_bash_unique_DEDUPLICATED_COMPRESS) && done($name_bash_unique_COMPRESS)" -J $name_clean_3 -M $mem -R"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]" -- \
    "rm $output_umi_tools_deduplicated_sorted $output_umi_summary $output_umi_tools_log $output_umi_tools_group"
    

    
#    exit

  fi #no spaces
done < "$input_file_2"


