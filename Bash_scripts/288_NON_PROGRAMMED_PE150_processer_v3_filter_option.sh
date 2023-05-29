#!/bin/bash>

   
#### Rscript

MASTER_ROUTE=$1
mem=$2
pc=$3
queue=$4
output=$5
R1=$6
R2=$7
INDEX=$8
lane=$9
 
MAPQ=$(echo "6")
#Bc_per_tile=$(echo "15")
# Threshold_reads_per_bc=$(echo "15")
# DEF_type=$(echo "MAPQ_6_15_15")
 Bc_per_tile=$(echo "10")
 Threshold_reads_per_bc=$(echo "5")
# Bc_per_tile=$(echo "5")
# Threshold_reads_per_bc=$(echo "2")
# Bc_per_tile=$(echo "1")
# Threshold_reads_per_bc=$(echo "1")
DEF_type=$(echo "MAPQ_6_""$Bc_per_tile""_""$Threshold_reads_per_bc")


touch $output
echo -n "" > $output

echo "#!/bin/bash"  >> $output



output_dir=$(echo "$MASTER_ROUTE")
type_for_R=$(echo "100_percent")
declare -a arr

# echo "#########################################################################################################################################################################"  >> $output
# echo "#########################################################################################################################################################################"  >> $output


### With the new sequencing primers this is going to change. The source of the UMI is the R1. If we do dark cycling the pattern will be Nx15. i fwe don't then the pattern will have before ten bases of constant sequence
### Then it is UMI tools extract from R1 to R1 and UMI tools extract from R1 to R2

type=$(echo "R1_processing")
outfile_R1_processing=$(echo "$output_dir""outfile""_""$type""_""$lane"".out")
touch $outfile_R1_processing
echo -n "" > $outfile_R1_processing
name_R1_processing=$(echo "$type""_""$lane""_job")

output_R1_processed=$(echo "$output_dir""$type""_""$lane"".fastq.gz")


echo "bsub -G team151 -o $outfile_R1_processing -M $mem  -J $name_R1_processing -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -n$pc -q $queue -- \\" >> $output
echo "\"umi_tools extract --bc-pattern=NNNNNNNNNNNNNNN --stdin=$INDEX --read2-in=$R1 --stdout=$output_R1_processed --read2-stdout\"" >> $output

type=$(echo "R2_processing")
outfile_R2_processing=$(echo "$output_dir""outfile""_""$type""_""$lane"".out")
touch $outfile_R2_processing
echo -n "" > $outfile_R2_processing
name_R2_processing=$(echo "$type""_""$lane""_job")

output_R2_processed=$(echo "$output_dir""$type""_""$lane"".fastq.gz")


echo "bsub -G team151 -o $outfile_R2_processing -M $mem  -J $name_R2_processing -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -n$pc -q $queue -- \\" >> $output
echo "\"umi_tools extract --bc-pattern=NNNNNNNNNNNNNNN --stdin=$INDEX --read2-in=$R2 --stdout=$output_R2_processed --read2-stdout\"" >> $output

type=$(echo "flash_merging")
outfile_flash_merging=$(echo "$output_dir""outfile""_""$type""_""$lane"".out")
touch $outfile_flash_merging
echo -n "" > $outfile_flash_merging
name_flash_merging=$(echo "$type""_""$lane""_job")


echo "bsub -G team151 -o $outfile_flash_merging -M $mem -w\"done($name_R1_processing) && done($name_R2_processing)\" -J $name_flash_merging -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -n$pc -q $queue -- \\" >> $output
#echo "bsub -G team151 -o $outfile_flash_merging -M $mem -J $name_flash_merging -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -n$pc -q $queue -- \\" >> $output
echo "\"flash2 -z -t 1 -r 151 -f 270 -s 27 $output_R1_processed $output_R2_processed --output-prefix=$lane\"" >> $output

echo "#########################################################################################################################################################################"  >> $output
echo "#########################################################################################################################################################################"  >> $output

#ALIGNMENT #########

type=$(echo "bwa_aligning")
outfile_bwa_aligning=$(echo "$output_dir""outfile""_""$type""_""$lane"".out")
touch $outfile_bwa_aligning
echo -n "" > $outfile_bwa_aligning
name_bwa_aligning=$(echo "$type""_""$lane""_job")


flash_output=$(echo "$output_dir""$lane"".extendedFrags.fastq.gz")
REFERENCE="/lustre/scratch119/realdata/mdt2/teams/soranzo/users/mt19/NEW_TWIST_LIBS/NON_PROGRAMMED/reference_files/NON_PROGRAMMED_reference.fasta"
output_bwa_bam=$(echo "$output_dir""$lane""_sorted.bam")
output_bwa_bam_filtered_sorted=$(echo "$output_dir""$lane""_filtered_sorted.bam")



echo "bsub -G team151 -o $outfile_bwa_aligning -M $mem -w\"done($name_flash_merging)\" -J $name_bwa_aligning -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -n$pc -q $queue -- \\" >> $output
 #echo "bsub -G team151 -o $outfile_bwa_aligning -M $mem -J $name_bwa_aligning -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -n$pc -q $queue -- \\" >> $output
 echo "\"bwa mem -M -t $pc $REFERENCE $flash_output |samtools view -@$pc -q $MAPQ -h -F 256|samtools sort -@$pc -o $output_bwa_bam_filtered_sorted\"" >> $output



type=$(echo "indexing_filtered")
name_indexing_filtered=$(echo "$type""_""$lane""_job")


echo "bsub -G team151 -o $outfile_bwa_aligning -M $mem -J $name_indexing_filtered -w\"done($name_bwa_aligning)\" -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -n$pc -q $queue -- \\" >> $output
echo "\"samtools index $output_bwa_bam_filtered_sorted\"" >> $output


echo "#########################################################################################################################################################################"  >> $output
echo "#########################################################################################################################################################################"  >> $output


type=$(echo "umi_tools_grouping")
outfile_umi_tools_grouping=$(echo "$output_dir""outfile""_""$type""_""$lane"".out")
touch $outfile_umi_tools_grouping
echo -n "" > $outfile_umi_tools_grouping
name_umi_tools_grouping=$(echo "$type""_""$lane""_job")


output_umi_tools_group=$(echo "$output_dir""$lane""_group.tsv")
output_umi_tools_log=$(echo "$output_dir""$lane""_group.log")



echo "bsub -G team151 -o $outfile_umi_tools_grouping -M $mem -J $name_umi_tools_grouping -w\"done($name_indexing_filtered)\" -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -n$pc -q $queue -- \\" >> $output
#echo "bsub -G team151 -o $outfile_umi_tools_grouping -M $mem -J $name_umi_tools_grouping -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -n$pc -q $queue -- \\" >> $output
echo "\"umi_tools group -I $output_bwa_bam_filtered_sorted --group-out=$output_umi_tools_group --log=$output_umi_tools_log\"" >> $output
 


	type=$(echo "Read_sorting")
        name_Read_collapsing_1=$(echo "$type""_""$lane""_""_job")
	output_umi_tools_group_sorted=$(echo "$output_dir""$type""_""$lane""_""sorted"".tsv")

	echo "bsub -G team151 -o $outfile_umi_tools_grouping -M $mem  -w\"done($name_umi_tools_grouping)\" -J $name_Read_collapsing_1  -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -n$pc -q $queue -- \\" >> $output
#	echo "bsub -G team151 -o $outfile_umi_tools_grouping -M $mem  -J $name_Read_collapsing_1  -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -n$pc -q $queue -- \\" >> $output
	echo "\"cut -f2,7 $output_umi_tools_group|sort --parallel=$pc -S 20G > $output_umi_tools_group_sorted\"" >> $output

	type=$(echo "Read_unique")
        name_Read_collapsing_2=$(echo "$type""_""$lane""_""_job")
	output_umi_tools_group_unique=$(echo "$output_dir""$type""_""$lane""_""sorted_unique"".tsv")

	echo "bsub -G team151 -o $outfile_umi_tools_grouping -M $mem  -w\"done($name_Read_collapsing_1)\" -J $name_Read_collapsing_2  -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -n$pc -q $queue -- \\" >> $output
#	echo "bsub -G team151 -o $outfile_umi_tools_grouping -M $mem  -J $name_Read_collapsing_2  -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -n$pc -q $queue -- \\" >> $output
	echo "\"uniq -c $output_umi_tools_group_sorted > $output_umi_tools_group_unique\"" >> $output

	type=$(echo "Read_compress")
        name_Read_collapsing_3=$(echo "$type""_""$lane""_""_job")
	output_umi_tools_group_collapsed=$(echo "$output_umi_tools_group_unique"".gz")
#	output_umi_tools_group_collapsed=$(echo "$output_umi_tools_group"".gz")

	echo "bsub -G team151 -o $outfile_umi_tools_grouping -M $mem  -w\"done($name_Read_collapsing_2)\" -J $name_Read_collapsing_3  -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -n$pc -q $queue -- \\" >> $output
#	echo "bsub -G team151 -o $outfile_umi_tools_grouping -M $mem  -J $name_Read_collapsing_3  -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -n$pc -q $queue -- \\" >> $output
	echo "\"gzip -f $output_umi_tools_group_unique\"" >> $output
	

	type=$(echo "CLEANING_2")
        name_CLEANING=$(echo "$type""_""$lane""_""_job")

	echo "bsub -G team151 -o $outfile_umi_tools_grouping -M 4000  -w\"done($name_Read_collapsing_3)\" -J $name_CLEANING  -R\"select[mem>=4000] rusage[mem=4000] span[hosts=1]\" -n1 -q $queue -- \\" >> $output
	echo "\"rm $output_umi_tools_group $output_umi_tools_group_sorted $output_umi_tools_group_unique\"" >> $output



	
echo "#########################################################################################################################################################################"  >> $output
echo "#########################################################################################################################################################################"  >> $output

############### Rscript

type=$(echo "Read_unique")
output_umi_tools_group_unique=$(echo "$output_dir""$type""_""$lane""_""sorted_unique"".tsv")
output_umi_tools_group_collapsed=$(echo "$output_umi_tools_group_unique"".gz")

type=$(echo "R_processing")
outfile_R_processing=$(echo "$output_dir""outfile""_""$type""_""$lane"".out")
touch $outfile_R_processing
echo -n "" > $outfile_R_processing
name_R_processing=$(echo "$type""_""$lane""_job")



Rscript="/software/R-3.6.1/bin/Rscript"
Rscript_filtering_per_fraction="/nfs/users/nfs_m/mt19/Scripts/R/236_Read_PE150_Non_Programmed.R"


MPRA_Rosetta=$(echo "/lustre/scratch119/realdata/mdt2/teams/soranzo/users/mt19/NEW_TWIST_LIBS/NON_PROGRAMMED/NON_PROGRAMMED_MPRA_Rosetta.tsv")
Threshold_min_length_barcode=$(echo "15")
#Threshold_reads_per_bc=$(echo "2")
Threshold_min_bc_per_tile=$(echo $Bc_per_tile)
FINAL_type=$(echo "$lane""_""$type_for_R")



echo "bsub -G team151 -o $outfile_R_processing -M $mem -J $name_R_processing -w\"done($name_Read_collapsing_3)\" -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -n$pc -q $queue -- \\" >> $output
#echo "bsub -G team151 -o $outfile_R_processing -M $mem -J $name_R_processing -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -n$pc -q $queue -- \\" >> $output
echo "\"$Rscript $Rscript_filtering_per_fraction \\" >> $output
echo "--input $output_umi_tools_group_collapsed \\" >> $output
echo "--MPRA_Rosetta $MPRA_Rosetta \\" >> $output
echo "--Threshold_min_length_barcode $Threshold_min_length_barcode \\" >> $output
echo "--Threshold_reads_per_bc $Threshold_reads_per_bc \\" >> $output
echo "--Threshold_min_bc_per_tile $Threshold_min_bc_per_tile \\" >> $output
echo "--type $FINAL_type --out $output_dir\"" >> $output


name_R_processing_100_percent=$(echo "$name_R_processing")

# R_processing_string=$(echo "&& done($name_R_processing)")
# echo "->>>$R_processing_string"
# arr[${#arr[@]}]="$R_processing_string"

#echo "bash $output"
#exit


echo "#########################################################################################################################################################################"  >> $output
echo "#########################################################################################################################################################################"  >> $output


subsampling_string=$(echo '0.01_''0.1_''0.25_''0.5_''0.75_''0.9_')

echo $subsampling_string

a=($(echo "$subsampling_string" | tr "_" '\n'))



for i  in "${a[@]}"
    do
	subsampling_fraction=${i}
	echo "$subsampling_fraction"
	
	type=$(echo "seqtk_fractioning")

	if [ $subsampling_fraction == "0.01" ] || [ $subsampling_fraction == "0.1" ] || [ $subsampling_fraction == "0.25" ];

	then

	    queue=$(echo "normal")

 	 else
	    queue=$4
	fi

	outfile_seqtk_fractioning=$(echo "$output_dir""outfile_seqtk_fractioning""_""$type""_""$lane""_""$subsampling_fraction"".out")
	touch $outfile_seqtk_fractioning
	echo -n "" > $outfile_seqtk_fractioning
	name_seqtk_fractioning=$(echo "$type""_""$lane""_""$subsampling_fraction""_job")

	output_fraction=$(echo "$output_dir""$type""_""$lane""_""$subsampling_fraction""FLASH_output.fastq.gz")

	echo "bsub -G team151 -o $outfile_seqtk_fractioning -M $mem -w\"done($name_flash_merging)\"  -J $name_seqtk_fractioning  -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -n$pc -q $queue -- \\" >> $output
#	echo "bsub -G team151 -o $outfile_seqtk_fractioning -M $mem  -J $name_seqtk_fractioning  -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -n$pc -q $queue -- \\" >> $output
	echo "\"seqtk sample -s 100 \\" >> $output
	echo "$flash_output \\" >> $output
	echo "$subsampling_fraction|gzip -f > \\" >> $output
	echo "$output_fraction\"" >> $output

echo "#########################################################################################################################################################################"  >> $output
echo "#########################################################################################################################################################################"  >> $output


#### ALIGNMENT #########

type=$(echo "bwa_aligning")
outfile_bwa_aligning=$(echo "$output_dir""outfile""_""$type""_""$lane""_""$subsampling_fraction"".out")
touch $outfile_bwa_aligning
echo -n "" > $outfile_bwa_aligning
name_bwa_aligning=$(echo "$type""_""$lane""_""$subsampling_fraction""_job")



REFERENCE=/lustre/scratch119/realdata/mdt2/teams/soranzo/users/mt19/NEW_TWIST_LIBS/NON_PROGRAMMED/reference_files/NON_PROGRAMMED_reference.fasta
output_bwa_bam=$(echo "$output_dir""$lane""_""$subsampling_fraction""_sorted.bam")
output_bwa_bam_filtered_sorted=$(echo "$output_dir""$lane""_""$subsampling_fraction""_filtered_sorted.bam")



echo "bsub -G team151 -o $outfile_bwa_aligning -M $mem -w\"done($name_seqtk_fractioning)\" -J $name_bwa_aligning -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -n$pc -q $queue -- \\" >> $output
#echo "bsub -G team151 -o $outfile_bwa_aligning -M $mem -J $name_bwa_aligning -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -n$pc -q $queue -- \\" >> $output
echo "\"bwa mem -M -t $pc $REFERENCE $output_fraction |samtools view -@$pc -q $MAPQ -h -F 256|samtools sort -@$pc -o $output_bwa_bam_filtered_sorted\"" >> $output



type=$(echo "indexing_filtered")
name_indexing_filtered=$(echo "$type""_""$lane""_""$subsampling_fraction""_job")


echo "bsub -G team151 -o $outfile_bwa_aligning -M $mem -J $name_indexing_filtered -w\"done($name_bwa_aligning)\" -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -n$pc -q $queue -- \\" >> $output
echo "\"samtools index $output_bwa_bam_filtered_sorted\"" >> $output


echo "#########################################################################################################################################################################"  >> $output
echo "#########################################################################################################################################################################"  >> $output

################## umi tools group ################

type=$(echo "umi_tools_grouping")
outfile_umi_tools_grouping=$(echo "$output_dir""outfile""_""$type""_""$lane""_""$subsampling_fraction"".out")
touch $outfile_umi_tools_grouping
echo -n "" > $outfile_umi_tools_grouping
name_umi_tools_grouping=$(echo "$type""_""$lane""_""$subsampling_fraction""_job")


output_umi_tools_group=$(echo "$output_dir""$lane""_""$subsampling_fraction""_group.tsv")
output_umi_tools_log=$(echo "$output_dir""$lane""_""$subsampling_fraction""_group.log")



echo "bsub -G team151 -o $outfile_umi_tools_grouping -M $mem -J $name_umi_tools_grouping -w\"done($name_indexing_filtered)\" -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -n$pc -q $queue -- \\" >> $output
echo "\"umi_tools group -I $output_bwa_bam_filtered_sorted --group-out=$output_umi_tools_group --log=$output_umi_tools_log\"" >> $output



	type=$(echo "Read_sorting")
        name_Read_collapsing_1=$(echo "$type""_""$lane""_""$subsampling_fraction""_""_job")
	output_umi_tools_group_sorted=$(echo "$output_dir""$type""_""$lane""_""$subsampling_fraction""_""sorted"".tsv")

	echo "bsub -G team151 -o $outfile_umi_tools_grouping -M $mem  -w\"done($name_umi_tools_grouping)\" -J $name_Read_collapsing_1  -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -n$pc -q $queue -- \\" >> $output
#	echo "bsub -G team151 -o $outfile_umi_tools_grouping -M $mem  -J $name_Read_collapsing_1  -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -n$pc -q $queue -- \\" >> $output
	echo "\"cut -f2,7 $output_umi_tools_group|sort --parallel=$pc -S 20G > $output_umi_tools_group_sorted\"" >> $output

	type=$(echo "Read_unique")
        name_Read_collapsing_2=$(echo "$type""_""$lane""_""$subsampling_fraction""_""_job")
	output_umi_tools_group_unique=$(echo "$output_dir""$type""_""$lane""_""$subsampling_fraction""_""sorted_unique"".tsv")

	echo "bsub -G team151 -o $outfile_umi_tools_grouping -M $mem  -w\"done($name_Read_collapsing_1)\" -J $name_Read_collapsing_2  -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -n$pc -q $queue -- \\" >> $output
#	echo "bsub -G team151 -o $outfile_umi_tools_grouping -M $mem  -J $name_Read_collapsing_2  -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -n$pc -q $queue -- \\" >> $output
	echo "\"uniq -c $output_umi_tools_group_sorted > $output_umi_tools_group_unique\"" >> $output

	type=$(echo "Read_compress")
        name_Read_collapsing_3=$(echo "$type""_""$lane""_""$subsampling_fraction""_""_job")
	output_umi_tools_group_collapsed=$(echo "$output_umi_tools_group_unique"".gz")
#	output_umi_tools_group_collapsed=$(echo "$output_umi_tools_group"".gz")

	echo "bsub -G team151 -o $outfile_umi_tools_grouping -M $mem  -w\"done($name_Read_collapsing_2)\" -J $name_Read_collapsing_3  -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -n$pc -q $queue -- \\" >> $output
#	echo "bsub -G team151 -o $outfile_umi_tools_grouping -M $mem  -J $name_Read_collapsing_3  -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -n$pc -q $queue -- \\" >> $output
	echo "\"gzip -f $output_umi_tools_group_unique\"" >> $output
	

	type=$(echo "CLEANING_2")
        name_CLEANING=$(echo "$type""_""$lane""_""$subsampling_fraction""_""_job")

	echo "bsub -G team151 -o $outfile_umi_tools_grouping -M 4000  -w\"done($name_Read_collapsing_3)\" -J $name_CLEANING  -R\"select[mem>=4000] rusage[mem=4000] span[hosts=1]\" -n1 -q $queue -- \\" >> $output
	echo "\"rm $output_umi_tools_group $output_umi_tools_group_sorted $output_umi_tools_group_unique\"" >> $output
	




echo "#########################################################################################################################################################################"  >> $output
echo "#########################################################################################################################################################################"  >> $output

############### Rscript
type=$(echo "Read_unique")
output_umi_tools_group_unique=$(echo "$output_dir""$type""_""$lane""_""$subsampling_fraction""_""sorted_unique"".tsv")
output_umi_tools_group_collapsed=$(echo "$output_umi_tools_group_unique"".gz")

type=$(echo "R_processing")
outfile_R_processing=$(echo "$output_dir""outfile""_""$type""_""$lane""_""$subsampling_fraction"".out")
touch $outfile_R_processing
echo -n "" > $outfile_R_processing
name_R_processing=$(echo "$type""_""$lane""_""$subsampling_fraction""_job")


Rscript=/software/R-3.6.1/bin/Rscript
Rscript_filtering_per_fraction=/nfs/users/nfs_m/mt19/Scripts/R/236_Read_PE150_Non_Programmed.R


MPRA_Rosetta=$(echo "/lustre/scratch119/realdata/mdt2/teams/soranzo/users/mt19/NEW_TWIST_LIBS/NON_PROGRAMMED/NON_PROGRAMMED_MPRA_Rosetta.tsv")
Threshold_min_length_barcode=$(echo "15")

FINAL_type=$(echo "$lane""_""$subsampling_fraction")




echo "bsub -G team151 -o $outfile_R_processing -M $mem -J $name_R_processing -w\"done($name_Read_collapsing_3)\" -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -n$pc -q $queue -- \\" >> $output
#echo "bsub -G team151 -o $outfile_R_processing -M $mem -J $name_R_processing -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -n$pc -q $queue -- \\" >> $output
echo "\"$Rscript $Rscript_filtering_per_fraction \\" >> $output
echo "--input $output_umi_tools_group_collapsed \\" >> $output
echo "--MPRA_Rosetta $MPRA_Rosetta \\" >> $output
echo "--Threshold_min_length_barcode $Threshold_min_length_barcode \\" >> $output
echo "--Threshold_reads_per_bc $Threshold_reads_per_bc \\" >> $output
echo "--Threshold_min_bc_per_tile $Threshold_min_bc_per_tile \\" >> $output
echo "--type $FINAL_type --out $output_dir\"" >> $output

R_processing_string=$(echo "&& done($name_R_processing)")
#echo "->>>$R_processing_string"
arr[${#arr[@]}]="$R_processing_string"

done


done_string=$(echo ""${arr[@]}"""\"")
first_half_done_string=$(echo "\"""done($name_R_processing_100_percent) ")

complete_done_string=$(echo "$first_half_done_string""$done_string")

# echo "$complete_done_string"


echo "#########################################################################################################################################################################"  >> $output
echo "#########################################################################################################################################################################"  >> $output
echo "#########################################################################################################################################################################"  >> $output
echo "#########################################################################################################################################################################"  >> $output
echo "#########################################################################################################################################################################"  >> $output



Rscript_ALL_TOGETHER=/nfs/users/nfs_m/mt19/Scripts/R/237_NON_PROGRAMMED_MPRA_ALL_TOGETHER.R

Initial_Selection=$(echo '/lustre/scratch119/realdata/mdt2/teams/soranzo/users/mt19/V2F_paper/ER_Labelling_Initial_Selection.txt')
dB_ALL=$(echo '/lustre/scratch119/realdata/mdt2/teams/soranzo/users/mt19/ALL_db.tsv')

MPRA_Rosetta=$(echo "/lustre/scratch119/realdata/mdt2/teams/soranzo/users/mt19/NEW_TWIST_LIBS/NON_PROGRAMMED/NON_PROGRAMMED_MPRA_Rosetta.tsv")
type=$(echo "FINAL_processing")

outfile_ALL_TOGETHER=$(echo "$output_dir""outfile_ALL_TOGETHER""_""$type""_""$lane"".out")
touch $outfile_ALL_TOGETHER
echo -n "" > $outfile_ALL_TOGETHER
name_ALL_TOGETHER=$(echo "$type""_""$lane""_job")


echo "bsub -G team151 -o $outfile_ALL_TOGETHER -M $mem  -w$complete_done_string -J $name_ALL_TOGETHER -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -n$pc -q $queue -- \\" >> $output
echo "\"$Rscript $Rscript_ALL_TOGETHER \\" >> $output
echo "--MPRA_Rosetta $MPRA_Rosetta \\" >> $output
echo "--type $DEF_type --out $output_dir\"" >> $output

echo "#########################################################################################################################################################################"  >> $output

echo "bash $output"

