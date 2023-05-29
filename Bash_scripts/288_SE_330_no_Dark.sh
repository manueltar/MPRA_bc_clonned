#!/bin/bash>

   
#### Rscript

MASTER_ROUTE=$1
mem=$2
pc=$3
queue=$4
output=$5
#R1=$6
R2=$6
lane=$7
 
#MAPQ=$(echo "6")
#Bc_per_tile=$(echo "15")
# Threshold_reads_per_bc=$(echo "15")
# DEF_type=$(echo "MAPQ_6_15_15")
# Bc_per_tile=$(echo "10")
# Threshold_reads_per_bc=$(echo "5")

# Bc_per_tile=$(echo "1")
# Threshold_reads_per_bc=$(echo "1")

Bc_per_tile=$(echo "1")
Threshold_reads_per_bc=$(echo "1")

# Bc_per_tile=$(echo "1")
# Threshold_reads_per_bc=$(echo "1")
#DEF_type=$(echo "MAPQ_6_""$Bc_per_tile""_""$Threshold_reads_per_bc")

bash_script=$(echo "/nfs/users/nfs_m/mt19/Scripts/Wraper_scripts/268_MPRA_Sample_procesing_UMI_TOOLS_v5_mapq4_NEW_PRIMERS_sub_bash.sh")
DEF_type=$(echo "NoNM_MAPQ_6_""$Bc_per_tile""_""$Threshold_reads_per_bc")

#bash_script=$(echo "/nfs/users/nfs_m/mt19/Scripts/Wraper_scripts/268_MPRA_Sample_procesing_UMI_TOOLS_v6_NM_mapq6_NEW_PRIMERS_sub_bash.sh")
#DEF_type=$(echo "YesNM_MAPQ_6_""$Bc_per_tile""_""$Threshold_reads_per_bc")





touch $output
echo -n "" > $output

echo "#!/bin/bash"  >> $output



output_dir=$(echo "$MASTER_ROUTE")




echo "#########################################################################################################################################################################"  >> $output
echo "################################################# TRIMMING CONSTANT SEQUENCE ###################################################################################################################"  >> $output

type=$(echo "R2_trimming")
outfile_trimmomatic_r2=$(echo "$output_dir""outfile""_""$type""_""$lane"".out")
touch $outfile_trimmomatic_r2
echo -n "" > $outfile_trimmomatic_r2

name_trimmomatic_r2=$(echo "$type""_""$lane""_job")
output_trimmomatic_r2=$(echo "$output_dir""$type""_""$lane"".fastq.gz")
output_r2_unzip=$(echo "$output_dir""$type""_""$lane"".fastq")
echo "->$name_trimmomatic $output_trimmomatic_r2<-"

echo "bsub -G team151 -o $outfile_trimmomatic_r2 -q $queue -n$pc -J $name_trimmomatic_r2 -M $mem -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -- \\" >> $output
echo "\"java -jar /nfs/users/nfs_m/mt19/sOFTWARE/Trimmomatic-0.39/trimmomatic-0.39.jar SE -threads $pc $R2 $output_trimmomatic_r2 HEADCROP:14\"" >> $output

echo  -e "\n" >> $output
echo  -e "\n" >> $output



type=$(echo "r2_unzip")
outfile_r2_unzip=$(echo "$output_dir""outfile""_""$type""_""$lane"".out")
touch $outfile_r2_unzip
echo -n "" > $outfile_r2_unzip

name_r2_unzip=$(echo "$type""_""$lane""_job")


echo "bsub -G team151 -o $outfile_r2_unzip -q $queue -n$pc -w\"done($name_trimmomatic_r2)\" -J $name_r2_unzip -M $mem -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -- \\" >> $output
echo "\"gunzip $output_trimmomatic_r2\"" >> $output


#$ gunzip R2_trimming_lane_1.fastq.gz


echo "################################################# REVERSE COMPLEMENT ###################################################################################################################"  >> $output

hgi_base=$(echo "hgi_base")

echo "source /software/hgi/installs/anaconda3/etc/profile.d/conda.sh"  >> $output

echo "conda deactivate" >> $output

echo "conda activate $hgi_base" >> $output


type=$(echo "R2_REV_COMP")
outfile_rev_comp_r2=$(echo "$output_dir""outfile""_""$type""_""$lane"".out")
touch $outfile_rev_comp_r2
echo -n "" > $outfile_rev_comp_r2

name_rev_comp_r2=$(echo "$type""_""$lane""_job")
output_rev_comp_r2=$(echo "$output_dir""$type""_""$lane"".fastq.gz")
echo "->$name_trimmomatic $output_rev_comp_r2<-"

echo "bsub -G team151 -o $outfile_rev_comp_r2 -q $queue -n$pc -w\"done($name_r2_unzip)\" -J $name_rev_comp_r2 -M $mem -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -- \\" >> $output
echo "\"fastx_reverse_complement -z -i $output_r2_unzip -o $output_rev_comp_r2\"" >> $output

echo  -e "\n" >> $output
echo  -e "\n" >> $output

type=$(echo "clean_r2_1")
outfile_clean_r2_1=$(echo "$output_dir""outfile""_""$type""_""$lane"".out")
touch $outfile_clean_r2_1
echo -n "" > $outfile_clean_r2_1

name_clean_r2_1=$(echo "$type""_""$lane""_job")


echo "bsub -G team151 -o $outfile_clean_r2_1 -q $queue -n$pc -w\"done($name_rev_comp_r2)\" -J $name_clean_r2_1 -M $mem -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -- \\" >> $output
echo "\"rm -rf $output_r2_unzip\"" >> $output

echo  -e "\n" >> $output
echo  -e "\n" >> $output

#echo "bash $output"
#exit

echo "################################################# TRIMMING CONSTANT SEQUENCE PART II ###################################################################################################################"  >> $output

 type=$(echo "R2_trimming_partII")
outfile_trimmomatic_r2_partII=$(echo "$output_dir""outfile""_""$type""_""$lane"".out")
touch $outfile_trimmomatic_r2_partII
echo -n "" > $outfile_trimmomatic_r2_partII

name_trimmomatic_r2_partII=$(echo "$type""_""$lane""_job")
 output_trimmomatic_r2_partII=$(echo "$output_dir""$type""_""$lane"".fastq.gz")

echo "->$name_trimmomatic $output_trimmomatic_r2_partII<-"

echo "bsub -G team151 -o $outfile_trimmomatic_r2_partII -q $queue -n$pc -w\"done($name_rev_comp_r2)\" -J $name_trimmomatic_r2_partII -M $mem -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -- \\" >> $output
echo "\"java -jar /nfs/users/nfs_m/mt19/sOFTWARE/Trimmomatic-0.39/trimmomatic-0.39.jar SE -threads $pc $output_rev_comp_r2 $output_trimmomatic_r2_partII HEADCROP:4\"" >> $output


echo "################################################# UMI-TOOLS EXTRACT ###################################################################################################################"  >> $output
## problem HERE



type=$(echo "R2_barcode_extraction")
outfile_R2_barcode_extraction=$(echo "$output_dir""outfile""_""$type""_""$lane"".out")
touch $outfile_R2_barcode_extraction
echo -n "" > $outfile_R2_barcode_extraction
name_R2_barcode_extraction=$(echo "$type""_""$lane""_job")

umi_tools="/nfs/team151/software/umi-tools/bin/umi_tools"

output_R2_processed=$(echo "$output_dir""$type""_""$lane"".fastq.gz")
#output_trimmomatic_r2_partII=$(echo "$output_dir""$type""_""$lane"".fastq.gz")


echo "bsub -G team151 -o $outfile_R2_barcode_extraction  -M $mem -w\"done($name_trimmomatic_r2_partII)\" -J $name_R2_barcode_extraction -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -n$pc -q $queue -- \\" >> $output
#echo "bsub -G team151 -o $outfile_R2_barcode_extraction  -M $mem -J $name_R2_barcode_extraction -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -n$pc -q $queue -- \\" >> $output
echo "\"$umi_tools extract --bc-pattern=NNNNNNNNNNNNNNN --stdin=$output_trimmomatic_r2_partII --read2-in=$output_trimmomatic_r2_partII --stdout=$output_R2_processed --read2-stdout\"" >> $output


echo  -e "\n" >> $output
echo  -e "\n" >> $output

echo "################################################# TRIMMING CONSTANT SEQUENCE PART III ###################################################################################################################"  >> $output

 type=$(echo "R2_trimming_partIII")
outfile_trimmomatic_r2_partIII=$(echo "$output_dir""outfile""_""$type""_""$lane"".out")
touch $outfile_trimmomatic_r2_partIII
echo -n "" > $outfile_trimmomatic_r2_partIII

name_trimmomatic_r2_partIII=$(echo "$type""_""$lane""_job")
 output_trimmomatic_r2_partIII=$(echo "$output_dir""$type""_""$lane"".fastq.gz")
output_r2_unzip=$(echo "$output_dir""$type""_""$lane"".fastq")
echo "->$name_trimmomatic $output_trimmomatic_r2_partIII<-"

echo "bsub -G team151 -o $outfile_trimmomatic_r2_partIII -q $queue -n$pc -w\"done($name_R2_barcode_extraction)\" -J $name_trimmomatic_r2_partIII -M $mem -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -- \\" >> $output
echo "\"java -jar /nfs/users/nfs_m/mt19/sOFTWARE/Trimmomatic-0.39/trimmomatic-0.39.jar SE -threads $pc $output_R2_processed $output_trimmomatic_r2_partIII HEADCROP:42\"" >> $output

echo  -e "\n" >> $output
echo  -e "\n" >> $output

echo "################################################### ALIGNMENT ##################################################################################################################"  >> $output

#exit
type=$(echo "bwa_aligning")
outfile_bwa_aligning=$(echo "$output_dir""outfile""_""$type""_""$lane"".out")
#touch $outfile_bwa_aligning
#echo -n "" > $outfile_bwa_aligning
name_bwa_aligning=$(echo "$type""_""$lane""_job")



REFERENCE=/lustre/scratch119/realdata/mdt2/teams/soranzo/users/mt19/NEW_TWIST_LIBS/NON_PROGRAMMED/reference_files/NON_PROGRAMMED_reference.fasta
output_bwa_bam=$(echo "$output_dir""$lane"".bam")




echo "bsub -G team151 -o $outfile_bwa_aligning -M $mem -w\"done($name_trimmomatic_r2_partIII)\" -J $name_bwa_aligning -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -n$pc -q $queue -- \\" >> $output
# echo "bsub -G team151 -o $outfile_bwa_aligning -M $mem -J $name_bwa_aligning -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -n$pc -q $queue -- \\" >> $output
 echo "\"bwa mem -M -t $pc -A 4 -B 16 $REFERENCE $output_trimmomatic_r2_partIII > $output_bwa_bam\"" >> $output



 type=$(echo "exclude_secondary_alignments_and_sorting")
 name_exclude_secondary_alignments_and_sorting=$(echo "$type""_""$lane""_job")
 output_bwa_bam_sorted=$(echo "$output_dir""$lane""_sorted.bam")

 echo "bsub -G team151 -o $outfile_bwa_aligning -M $mem -w\"done($name_bwa_aligning)\" -J $name_exclude_secondary_alignments_and_sorting -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -n$pc -q $queue -- \\" >> $output
 echo "\"samtools view -@$pc -h -F 256 $output_bwa_bam|samtools sort -@$pc -o $output_bwa_bam_sorted\"" >> $output

# echo "bash $output"
# exit

type=$(echo "clean_r2_2")
outfile_clean_r2_2=$(echo "$output_dir""outfile""_""$type""_""$lane"".out")
touch $outfile_clean_r2_2
echo -n "" > $outfile_clean_r2_2

name_clean_r2_2=$(echo "$type""_""$lane""_job")

echo "bsub -G team151 -o $outfile_clean_r2_2 -q $queue -n$pc -w\"done($name_exclude_secondary_alignments_and_sorting)\" -J $name_clean_r2_2 -M $mem -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -- \\" >> $output
#echo "bsub -G team151 -o $outfile_clean_r2_2 -q $queue -n$pc -J $name_clean_r2_2 -M $mem -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -- \\" >> $output
echo "\"rm -rf $output_R2_processed $output_bwa_bam\"" >> $output

echo  -e "\n" >> $output
echo  -e "\n" >> $output

type=$(echo "indexing")
name_indexing=$(echo "$type""_""$lane""_job")


echo "bsub -G team151 -o $outfile_bwa_aligning -M $mem -J $name_indexing -w\"done($name_exclude_secondary_alignments_and_sorting)\" -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -n$pc -q $queue -- \\" >> $output
#echo "bsub -G team151 -o $outfile_bwa_aligning -M $mem -J $name_indexing  -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -n$pc -q $queue -- \\" >> $output
echo "\"samtools index $output_bwa_bam_sorted\"" >> $output


#exit


# echo "############################################# FILTER FINAL ############################################################################################################################"  >> $output



type=$(echo "filter_and_rest")
name_filter_and_rest=$(echo "$type""_""$lane""_job")
#bash_script=$(echo "/nfs/users/nfs_m/mt19/Scripts/Wraper_scripts/268_MPRA_Sample_procesing_UMI_TOOLS_v4_NEW_PRIMERS_sub_bash.sh")


outfile_filter_and_rest=$(echo "$output_dir""outfile""_""$type""_""$lane"".out")
touch $outfile_filter_and_rest
echo -n "" > $outfile_filter_and_rest
name_filter_and_rest=$(echo "$type""_""$lane""_job")
DEF_type=$(echo "$DEF_type""_""$lane")


#echo "bsub -G team151 -o $outfile_filter_and_rest -q $queue -n$pc -w\"done($name_indexing)\" -J $name_filter_and_rest -M $mem -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -- \\"  >> $output
echo "bsub -G team151 -o $outfile_filter_and_rest -M $mem -J $name_filter_and_rest -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -n$pc -q $queue -- \\"  >> $output
echo "\"bash $bash_script $MASTER_ROUTE $mem $pc $queue $lane $Bc_per_tile $Threshold_reads_per_bc $DEF_type\"" >> $output

echo -e "\n"  >> $output



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


