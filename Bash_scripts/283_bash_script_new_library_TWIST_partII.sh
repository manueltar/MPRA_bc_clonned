#!/bin/bash
  
Rscript=/software/R-3.6.1/bin/Rscript

output=$8

touch $output
echo -n "" > $output

echo -e "#!/bin/bash"  >> $output
echo -e "\n\n"  >> $output

### Library_CREATION



#Rscript_module_Library_CREATION=$(echo "/nfs/users/nfs_m/mt19/Scripts/R/225_MPRA_NEW_LIBS_DESIGN_REAL_TILE_to_seq_name.R")
Rscript_module_Library_CREATION=$(echo "/nfs/users/nfs_m/mt19/Scripts/R/225_MPRA_NEW_LIBS_DESIGN_REAL_TILE_to_seq_name_NOT_PROGRAMMED.R")
dB_ALL=$(echo "/lustre/scratch119/realdata/mdt2/teams/soranzo/users/mt19/ALL_db.tsv")
subset=$1
type=$2
MASTER_ROUTE=$3
mem=$4
pc=$5
queue=$6
REFERENCE=$7

output_dir=$(echo "$MASTER_ROUTE""/""$type""/")

##### step Library_CREATION


outfile_Library_CREATION=$(echo "$output_dir""outfile_Library_CREATION""_""$type"".out")
touch $outfile_Library_CREATION
echo -n "" > $outfile_Library_CREATION
name_Library_CREATION=$(echo "$type""_name_Library_CREATION")


span=$(echo "300")
#sunk_cost=$(echo "30,12,11")
sunk_cost=$(echo "30")
RE_sites_unstring=$(echo "GGATCC,GGTACC")
PCR_arms_unstring=$(echo "CTGCGCCTGATGCAG,GGTGCTCGCTATCG,AGGACCGGATCAACT,CCTGCAGGGAATTC,ACTGGCCGCTTGACG,CACTGCGGCTCCTG")
Blacklisted_sequences=$(echo "TAGAGCATGCACCGGtgata,ctgccggccgaattcgtcga")
Barcode_source=$(echo "/nfs/users/nfs_m/mt19/RareVar_2019/MIPRA_design/publications-data/BC_1.txt")
Vetted_barcodes=$(echo "/lustre/scratch119/realdata/mdt2/teams/soranzo/users/mt19/NEW_TWIST_LIBS/MPRA_E_A_Vetted_barcodes.tsv")
indexes_fasta=$(echo "/nfs/users/nfs_m/mt19/RareVar_2019/DEF_design/Resynthesis/indexes.fasta")
indexes_96=$(echo "/nfs/users/nfs_m/mt19/RareVar_2019/DEF_design/Resynthesis/New_indexes_Illumina.tsv")


echo "bsub -G team151 -o $outfile_Library_CREATION -M $mem -J $name_Library_CREATION -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -n$pc -q $queue -- \\" >> $output
echo "\"$Rscript $Rscript_module_Library_CREATION \\" >> $output
echo "--span $span \\" >> $output
echo "--sunk_cost $sunk_cost \\" >> $output
echo "--RE_sites_unstring $RE_sites_unstring \\" >> $output
echo "--PCR_arms_unstring $PCR_arms_unstring \\" >> $output
echo "--Blacklisted_sequences $Blacklisted_sequences \\" >> $output
echo "--Barcode_source $Barcode_source \\" >> $output
echo "--Vetted_barcodes $Vetted_barcodes \\" >> $output
echo "--indexes_fasta $indexes_fasta \\" >> $output
echo "--indexes_96 $indexes_96 \\" >> $output
echo "--type $type --out $output_dir\"" >> $output

echo "bash $output"
