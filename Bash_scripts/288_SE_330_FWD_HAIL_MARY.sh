#!/bin/bash>

   
#### Rscript

MASTER_ROUTE=$1
mem=$2
pc=$3
queue=$4
output=$5
#R1=$6
R1=$6
lane=$7
 
#MAPQ=$(echo "6")

# DEF_type=$(echo "MAPQ_6_15_15")
# Bc_per_tile=$(echo "10")
# Threshold_reads_per_bc=$(echo "5")

# Bc_per_tile=$(echo "1")
# Threshold_reads_per_bc=$(echo "1")

Bc_per_tile=$(echo "1")
Threshold_reads_per_bc=$(echo "1")

Bc_per_tile=$(echo "5")
Threshold_reads_per_bc=$(echo "2")


# Bc_per_tile=$(echo "1")
# Threshold_reads_per_bc=$(echo "1")
#DEF_type=$(echo "MAPQ_6_""$Bc_per_tile""_""$Threshold_reads_per_bc")

bash_script=$(echo "/nfs/users/nfs_m/mt19/Scripts/Wraper_scripts/268_MPRA_HAIL_MARY_sub_bash.sh")
DEF_type=$(echo "$Bc_per_tile""_""$Threshold_reads_per_bc")

#bash_script=$(echo "/nfs/users/nfs_m/mt19/Scripts/Wraper_scripts/268_MPRA_Sample_procesing_UMI_TOOLS_v6_NM_mapq6_NEW_PRIMERS_sub_bash.sh")
#DEF_type=$(echo "YesNM_MAPQ_6_""$Bc_per_tile""_""$Threshold_reads_per_bc")



touch $output
echo -n "" > $output

echo "#!/bin/bash"  >> $output



output_dir=$(echo "$MASTER_ROUTE")


# echo "################################################# UMI-TOOLS EXTRACT ###################################################################################################################"  >> $output


type=$(echo "R1_barcode_extraction")
# outfile_R1_barcode_extraction=$(echo "$output_dir""outfile""_""$type""_""$lane"".out")
# touch $outfile_R1_barcode_extraction
# echo -n "" > $outfile_R1_barcode_extraction
# name_R1_barcode_extraction=$(echo "$type""_""$lane""_job")

# umi_tools="/nfs/team151/software/umi-tools/bin/umi_tools"

output_R1_processed=$(echo "$output_dir""$type""_""$lane"".fastq.gz")
# #output_trimmomatic_r2_partII=$(echo "$output_dir""$type""_""$lane"".fastq.gz")


# echo "bsub -G team151 -o $outfile_R1_barcode_extraction  -M $mem -J $name_R1_barcode_extraction -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -n$pc -q $queue -- \\" >> $output
# echo "\"$umi_tools extract --bc-pattern=NNNNNNNNNNNNNNN --stdin=$R1 --read2-in=$R1 --stdout=$output_R1_processed\"" >> $output


# echo  -e "\n" >> $output
# echo  -e "\n" >> $output



echo "#########################################################################################################################################################################"  >> $output
echo "################################################# TRIMMING CONSTANT SEQUENCE ###################################################################################################################"  >> $output

# type=$(echo "R1_trimming_1")
# outfile_trimmomatic_r1=$(echo "$output_dir""outfile""_""$type""_""$lane"".out")
# touch $outfile_trimmomatic_r1
# echo -n "" > $outfile_trimmomatic_r1

# name_trimmomatic_r1=$(echo "$type""_""$lane""_job")
# output_trimmomatic_r1=$(echo "$output_dir""$type""_""$lane"".fastq.gz")
# output_r2_unzip=$(echo "$output_dir""$type""_""$lane"".fastq")
# echo "->$name_trimmomatic $output_trimmomatic_r1<-"

# #echo "bsub -G team151 -o $outfile_trimmomatic_r1 -q $queue -n$pc -w\"done($name_R1_barcode_extraction)\"  -J $name_trimmomatic_r1 -M $mem -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -- \\" >> $output
# echo "bsub -G team151 -o $outfile_trimmomatic_r1 -q $queue -n$pc -J $name_trimmomatic_r1 -M $mem -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -- \\" >> $output
# echo "\"java -jar /nfs/users/nfs_m/mt19/sOFTWARE/Trimmomatic-0.39/trimmomatic-0.39.jar SE -threads $pc $output_R1_processed $output_trimmomatic_r1 HEADCROP:27\"" >> $output

# echo  -e "\n" >> $output
# echo  -e "\n" >> $output


type=$(echo "R1_trimming_2")
# outfile_trimmomatic_r1_PARTII=$(echo "$output_dir""outfile""_""$type""_""$lane"".out")
# touch $outfile_trimmomatic_r1_PARTII
# echo -n "" > $outfile_trimmomatic_r1_PARTII

# name_trimmomatic_r1_PARTII=$(echo "$type""_""$lane""_job")
 output_trimmomatic_r1_PARTII=$(echo "$output_dir""$type""_""$lane"".fastq.gz")
# output_r2_unzip=$(echo "$output_dir""$type""_""$lane"".fastq")
# echo "->$name_trimmomatic $output_trimmomatic_r1_PARTII<-"

# echo "bsub -G team151 -o $outfile_trimmomatic_r1_PARTII -q $queue -n$pc -w\"done($name_trimmomatic_r1)\"  -J $name_trimmomatic_r1_PARTII -M $mem -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -- \\" >> $output
# #echo "bsub -G team151 -o $outfile_trimmomatic_r1_PARTII -q $queue -n$pc -J $name_trimmomatic_r1_PARTII -M $mem -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -- \\" >> $output
# echo "\"java -jar /nfs/users/nfs_m/mt19/sOFTWARE/Trimmomatic-0.39/trimmomatic-0.39.jar SE -threads $pc $output_trimmomatic_r1 $output_trimmomatic_r1_PARTII CROP:271\"" >> $output

# echo  -e "\n" >> $output
# echo  -e "\n" >> $output



echo "################################################### ALIGNMENT ##################################################################################################################"  >> $output

#exit
# type=$(echo "bwa_aligning")
# outfile_bwa_aligning=$(echo "$output_dir""outfile""_""$type""_""$lane"".out")
# #touch $outfile_bwa_aligning
# #echo -n "" > $outfile_bwa_aligning
# name_bwa_aligning=$(echo "$type""_""$lane""_job")


# REFERENCE=/lustre/scratch119/realdata/mdt2/teams/soranzo/users/mt19/NEW_TWIST_LIBS/NON_PROGRAMMED/reference_files/NON_PROGRAMMED_reference.fasta
# output_bwa_bam=$(echo "$output_dir""$lane"".bam")


# #echo "bsub -G team151 -o $outfile_bwa_aligning -M $mem -w\"done($name_trimmomatic_r1_PARTII)\" -J $name_bwa_aligning -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -n$pc -q $queue -- \\" >> $output
#  echo "bsub -G team151 -o $outfile_bwa_aligning -M $mem -J $name_bwa_aligning -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -n$pc -q $queue -- \\" >> $output
#  echo "\"bwa mem -M -t $pc -A 4 -B 16 $REFERENCE $output_trimmomatic_r1_PARTII > $output_bwa_bam\"" >> $output



#  type=$(echo "exclude_secondary_alignments_and_sorting")
#  name_exclude_secondary_alignments_and_sorting=$(echo "$type""_""$lane""_job")
#  output_bwa_bam_sorted=$(echo "$output_dir""$lane""_sorted.bam")

#  echo "bsub -G team151 -o $outfile_bwa_aligning -M $mem -w\"done($name_bwa_aligning)\" -J $name_exclude_secondary_alignments_and_sorting -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -n$pc -q $queue -- \\" >> $output
#  echo "\"samtools view -@$pc -h -F 256 $output_bwa_bam|samtools sort -@$pc -o $output_bwa_bam_sorted\"" >> $output

# # echo "bash $output"
# # exit

# type=$(echo "clean_1")
# outfile_clean_1=$(echo "$output_dir""outfile""_""$type""_""$lane"".out")
# touch $outfile_clean_1
# echo -n "" > $outfile_clean_1

# name_clean_1=$(echo "$type""_""$lane""_job")

# echo "bsub -G team151 -o $outfile_clean_1 -q $queue -n$pc -w\"done($name_exclude_secondary_alignments_and_sorting)\" -J $name_clean_1 -M $mem -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -- \\" >> $output
# #echo "bsub -G team151 -o $outfile_clean_1 -q $queue -n$pc -J $name_clean_1 -M $mem -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -- \\" >> $output
# echo "\"rm -rf $output_bwa_bam\"" >> $output

# echo  -e "\n" >> $output
# echo  -e "\n" >> $output

# type=$(echo "indexing")
# name_indexing=$(echo "$type""_""$lane""_job")


# echo "bsub -G team151 -o $outfile_bwa_aligning -M $mem -J $name_indexing -w\"done($name_exclude_secondary_alignments_and_sorting)\" -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -n$pc -q $queue -- \\" >> $output
# #echo "bsub -G team151 -o $outfile_bwa_aligning -M $mem -J $name_indexing  -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -n$pc -q $queue -- \\" >> $output
# echo "\"samtools index $output_bwa_bam_sorted\"" >> $output





# echo "############################################# FILTER FINAL ############################################################################################################################"  >> $output



type=$(echo "filter_and_rest")
name_filter_and_rest=$(echo "$type""_""$lane""_job")


outfile_filter_and_rest=$(echo "$output_dir""outfile""_""$type""_""$lane"".out")
touch $outfile_filter_and_rest
echo -n "" > $outfile_filter_and_rest
name_filter_and_rest=$(echo "$type""_""$lane""_job")
DEF_type=$(echo "$DEF_type""_""$lane")


#echo "bsub -G team151 -o $outfile_filter_and_rest -q $queue -n$pc -w\"done($name_indexing)\" -J $name_filter_and_rest -M $mem -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -- \\"  >> $output
echo "bsub -G team151 -o $outfile_filter_and_rest -M $mem -J $name_filter_and_rest -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -n$pc -q $queue -- \\"  >> $output
echo "\"bash $bash_script $MASTER_ROUTE $mem $pc $queue $lane $Bc_per_tile $Threshold_reads_per_bc $DEF_type\"" >> $output

echo -e "\n"  >> $output

bash $output
exit

#subsampling_string=$(echo '0.01_''0.1_''0.25_''0.5_''0.75_''0.9_')
subsampling_string=$(echo '0.9_')

echo $subsampling_string

a=($(echo "$subsampling_string" | tr "_" '\n'))



for i  in "${a[@]}"
    do
        subsampling_fraction=${i}
        echo "$subsampling_fraction"
	if [ $subsampling_fraction == "0.01" ] || [ $subsampling_fraction == "0.1" ] || [ $subsampling_fraction == "0.25" ];

        then

            queue=$(echo "normal")

         else
            queue=$queue
        fi

# 	type=$(echo "subsampling_bam""_""$subsampling_fraction")
# 	name_subsampling_bam=$(echo "$type""_""$lane""_""$subsampling_fraction""_job")


# 	outfile_subsampling_bam=$(echo "$output_dir""outfile""_""$type""_""$lane""_""$subsampling_fraction"".out")
# 	touch $outfile_subsampling_bam
# 	echo -n "" > $outfile_subsampling_bam
# 	name_subsampling_bam=$(echo "$type""_""$lane""_""$subsampling_fraction""_job")
# 	output_bwa_bam_sorted_subsampled=$(echo "$output_dir""$lane""_""$subsampling_fraction""_sorted.bam")

	
# 	echo "bsub -G team151 -o $outfile_subsampling_bam -M $mem -J $name_subsampling_bam -w\"done($name_exclude_secondary_alignments_and_sorting)\" -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -n$pc -q $queue -- \\" >> $output
# #	echo "bsub -G team151 -o $outfile_subsampling_bam -M $mem -J $name_subsampling_bam  -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -n$pc -q $queue -- \\" >> $output

# 	echo "\"samtools view -s $subsampling_fraction -b $output_bwa_bam_sorted > $output_bwa_bam_sorted_subsampled\"" >> $output

	

# 	type=$(echo "filter_and_rest""_""$subsampling_fraction")
# 	name_filter_and_rest=$(echo "$type""_""$lane""_""$subsampling_fraction""_job")


	outfile_filter_and_rest=$(echo "$output_dir""outfile""_""$type""_""$lane""_""$subsampling_fraction"".out")
	touch $outfile_filter_and_rest
	echo -n "" > $outfile_filter_and_rest
	name_filter_and_rest=$(echo "$type""_""$lane""_""$subsampling_fraction""_job")

#	DEF_type_adapted=$(echo "MAPQ_6_""$Bc_per_tile""_""$Threshold_reads_per_bc""_""$subsampling_fraction")

	lane_adapted=$(echo "$lane""_""$subsampling_fraction")
	DEF_type_adapted=$(echo "$DEF_type""_""$lane_adapted")

#	echo "bsub -G team151 -o $outfile_filter_and_rest -q $queue -n$pc -w\"done($name_subsampling_bam)\" -J $name_filter_and_rest -M $mem -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -- \\"  >> $output
	echo "bsub -G team151 -o $outfile_filter_and_rest -q $queue -n$pc -J $name_filter_and_rest -M $mem -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -- \\"  >> $output
	echo "\"bash $bash_script $MASTER_ROUTE $mem $pc $queue $lane_adapted $Bc_per_tile $Threshold_reads_per_bc $DEF_type_adapted\"" >> $output

	echo -e "\n"  >> $output
 done



echo "bash $output"


