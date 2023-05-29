#!/bin/bash
 
Rscript=/software/R-3.6.1/bin/Rscript

output=$8

touch $output
echo -n "" > $output

echo -e "#!/bin/bash"  >> $output
echo -e "\n\n"  >> $output



Rscript_module_step_I_Library_bed=$(echo "/nfs/users/nfs_m/mt19/Scripts/R/223_MPRA_NEW_LIBS_DESIGN_v2.R")
dB_ALL=$(echo "/lustre/scratch119/realdata/mdt2/teams/soranzo/users/mt19/ALL_db.tsv")
subset=$1
type=$2
MASTER_ROUTE=$3
mem=$4
pc=$5
queue=$6


output_dir=$(echo "$MASTER_ROUTE""/""$type""/")

rm -rf $output_dir
mkdir -p $output_dir

### step_I_Library_bed

outfile_step_I_Library_bed=$(echo "$output_dir""outfile_step_I_Library_bed""_""$type"".out")
touch $outfile_step_I_Library_bed
echo -n "" > $outfile_step_I_Library_bed
name_step_I_Library_bed=$(echo "$type""_name_step_I_Library_bed")



VJ_MPRA_1=$(echo "/lustre/scratch119/realdata/mdt2/teams/soranzo/users/mt19/NEW_TWIST_LIBS/Sankaran.MPRA.contructs.csv")
VJ_MPRA_2=$(echo "/lustre/scratch119/realdata/mdt2/teams/soranzo/users/mt19/NEW_TWIST_LIBS/MPRA_Sankaran.csv")
Fulco_1=$(echo "/lustre/scratch119/realdata/mdt2/teams/soranzo/users/mt19/NEW_TWIST_LIBS/aag2445_Table_S2_Engreitzz_CRISPR_i_qtls.csv")
Sankaran_significance_threshold=$(echo "0.05")
Fulco_1_CRISPRi_threshold_DOWN=$(echo "0")
Fulco_1_CRISPRi_threshold_UP=$(echo "0.001")
slidding_windows=$(echo "5")
span=$(echo "300")
#sunk_cost=$(echo "30,12,11")
sunk_cost=$(echo "30")
fractions=$(echo "10,3,2")

echo "bsub -G team151 -o $outfile_step_I_Library_bed -M $mem  -J $name_step_I_Library_bed -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -n$pc -q $queue -- \\" >> $output
echo "\"$Rscript $Rscript_module_step_I_Library_bed \\" >> $output
echo "--dB_ALL $dB_ALL \\" >> $output
echo "--subset $subset \\" >> $output
echo "--VJ_MPRA_1 $VJ_MPRA_1 \\" >> $output
echo "--VJ_MPRA_2 $VJ_MPRA_2 \\" >> $output
echo "--Fulco_1 $Fulco_1 \\" >> $output
echo "--Fulco_1_CRISPRi_threshold_DOWN $Fulco_1_CRISPRi_threshold_DOWN \\" >> $output
echo "--Fulco_1_CRISPRi_threshold_UP $Fulco_1_CRISPRi_threshold_UP \\" >> $output
echo "--Sankaran_significance_threshold $Sankaran_significance_threshold \\" >> $output
echo "--slidding_windows $slidding_windows \\" >> $output
echo "--span $span \\" >> $output
echo "--sunk_cost $sunk_cost \\" >> $output
echo "--fractions $fractions \\" >> $output
echo "--type $type --out $output_dir\"" >> $output



##### step IIa REF fasta

REFERENCE=$7

outfile_step_II_Library_bed_to_fasta=$(echo "$output_dir""outfile_step_II_Library_bed_to_fasta""_""$type"".out")
touch $outfile_step_II_Library_bed_to_fasta
echo -n "" > $outfile_step_II_Library_bed_to_fasta
name_step_II_Library_bed_to_fasta=$(echo "$type""_name_step_II_Library_bed_to_fasta")


bedfile=$(echo "$output_dir""$type""_TILES"".bed")

REF_fasta=$(echo "$output_dir""$type""_REF"".fasta")


echo "bsub -G team151 -o $outfile_step_II_Library_bed_to_fasta -M $mem  -w\"done($name_step_I_Library_bed)\" -J $name_step_II_Library_bed_to_fasta -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -n$pc -q $queue -- \\" >> $output
echo "\"bedtools getfasta -fi $REFERENCE -fo $REF_fasta -bed $bedfile\"" >> $output


##### step IIb ALT fasta


name_step_II_Library_bed_to_fasta_B_SIDE=$(echo "$type""_name_step_II_Library_bed_to_fasta_B_SIDE")
Rscript_module_step_II_Library_bed_to_fasta_B_SIDE=$(echo "/nfs/users/nfs_m/mt19/Scripts/R/224_MPRA_NEW_LIBS_DESIGN_step_bed_to_fasta.R")


echo "bsub -G team151 -o $outfile_step_II_Library_bed_to_fasta -M $mem -w\"done($name_step_II_Library_bed_to_fasta)\"  -J $name_step_II_Library_bed_to_fasta_B_SIDE -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -n$pc -q $queue -- \\" >> $output
echo "\"$Rscript $Rscript_module_step_II_Library_bed_to_fasta_B_SIDE \\" >> $output
echo "--type $type --out $output_dir\"" >> $output



##### step IIc ALT fasta

name_step_II_Library_bed_to_fasta_C_SIDE=$(echo "$type""_name_step_II_Library_bed_to_fasta_C_SIDE")
FASTA_ALTERNATE=$(echo "/nfs/users/nfs_m/mt19/Scripts/Wraper_scripts/91_GATK_ALTERNATE_genIE.sh")
master_file=$(echo "$output_dir""$type""_MASTER_TILES_PLUS_REF_AND_NCGR_SUBS.tsv")


echo "bsub -G team151 -o $outfile_step_II_Library_bed_to_fasta -w\"done($name_step_II_Library_bed_to_fasta_B_SIDE)\"  -J $name_step_II_Library_bed_to_fasta_C_SIDE -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -M $mem -n$pc -q $queue -- \\" >> $output
echo "\"bash $FASTA_ALTERNATE $master_file $MASTER_ROUTE $type $REFERENCE $mem $pc $queue\"" >> $output





