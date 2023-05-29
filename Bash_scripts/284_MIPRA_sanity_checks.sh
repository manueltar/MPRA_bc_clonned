#!/bin/bash
 
Rscript=/software/R-3.6.1/bin/Rscript
perl=perl

output=$7

touch $output
echo -n "" > $output

echo -e "#!/bin/bash"  >> $output
echo -e "\n\n"  >> $output

### Library_CREATION


MASTER_ROUTE=$1
Aligner_fmt3=$(echo "/nfs/users/nfs_m/mt19/Scripts/Perl/16_Blaster_GRCh37_MIPRA_edition_restricted.pl")
input=$2
reference="/lustre/scratch119/realdata/mdt2/teams/soranzo/users/mt19/NEW_TWIST_LIBS/Homo_sapiens.GRCh37.dna.all.fa"
global_file_fmt3=$(echo "$MASTER_ROUTE""/""global_file_fmt3.txt")
type=$3
mem=$4
pc=$5
queue=$6


output_dir=$(echo "$MASTER_ROUTE")

##### step BLAST_fmt3


outfile_BLAST_fmt3=$(echo "$output_dir""outfile_BLAST_fmt3""_""$type"".out")
touch $outfile_BLAST_fmt3
echo -n "" > $outfile_BLAST_fmt3
name_BLAST_fmt3=$(echo "$type""_name_BLAST_fmt3")


echo "bsub -G team151 -o $outfile_BLAST_fmt3 -M $mem -J $name_BLAST_fmt3 -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -n$pc -q $queue -- \\" >> $output
echo "\"$perl $Aligner_fmt3 $input $reference $global_file_fmt3\"" >> $output

##### step BLAST_fmt7

Aligner_fmt7=$(echo "/nfs/users/nfs_m/mt19/Scripts/Perl/16_Blaster_GRCh37_fmt7.sh")
global_file_fmt7=$(echo "$MASTER_ROUTE""/""global_file_fmt7.txt")


outfile_BLAST_fmt7=$(echo "$output_dir""outfile_BLAST_fmt7""_""$type"".out")
touch $outfile_BLAST_fmt7
echo -n "" > $outfile_BLAST_fmt7
name_BLAST_fmt7=$(echo "$type""_name_BLAST_fmt7")


echo "bsub -G team151 -o $outfile_BLAST_fmt7 -M $mem -J $name_BLAST_fmt7 -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -n$pc -q $queue -- \\" >> $output
echo "\"$perl $Aligner_fmt7 $input $reference $global_file_fmt7\"" >> $output

##### step parse_fmt3

parser_fmt3=$(echo "/nfs/users/nfs_m/mt19/Scripts/Perl/17_Parsers_GRCh38_5.pl")

outfile_parse_fmt3=$(echo "$output_dir""outfile_parse_fmt3""_""$type"".out")
touch $outfile_parse_fmt3
echo -n "" > $outfile_parse_fmt3
name_parse_fmt3=$(echo "$type""_name_parse_fmt3")

parse_file_fmt3=$(echo "$MASTER_ROUTE""/""global_file_fmt3_parsed.tsv")


echo "bsub -G team151 -o $outfile_parse_fmt3 -M $mem -w\"done($name_BLAST_fmt3)\"  -J $name_parse_fmt3 -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -n$pc -q $queue -- \\" >> $output
echo "\"$perl $parser_fmt3 $global_file_fmt3 $parse_file_fmt3\"" >> $output

##### step R_check_and_report

R_check_and_report=$(echo "/nfs/users/nfs_m/mt19/Scripts/R/227_MIPRA_SANITY_CHECK.R")

outfile_R_check_and_report=$(echo "$output_dir""outfile_R_check_and_report""_""$type"".out")
touch $outfile_R_check_and_report
echo -n "" > $outfile_R_check_and_report
name_R_check_and_report=$(echo "$type""_name_R_check_and_report")

parse_file_fmt3=$(echo "$MASTER_ROUTE""/""global_file_fmt3_parsed.tsv")


Open_targets_GLOBAL_file=$(echo "/lustre/scratch119/realdata/mdt2/teams/soranzo/users/ALL_dB/Allelic_Series_csv_tables_AS_DEFINITIVE/AS_DEFINITIVE_Open_Targets_GWAS_DISEASES_variants.tsv")
Selection_file=$(echo "/lustre/scratch119/realdata/mdt2/teams/soranzo/users/mt19/NEW_TWIST_LIBS/Selection_""$type"".tsv")

echo "bsub -G team151 -o $outfile_R_check_and_report -M $mem -w\"done($name_parse_fmt3) && done($name_BLAST_fmt7)\"  -J $name_R_check_and_report -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -n$pc -q $queue -- \\" >> $output
echo "\"$Rscript $R_check_and_report --Selection_file $Selection_file --Open_targets_GLOBAL_file $Open_targets_GLOBAL_file --type $type --out $output_dir\"" >> $output

#echo "bsub -G team151 -o $outfile_R_check_and_report -M $mem  -J $name_R_check_and_report -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -n$pc -q $queue -- \\" >> $output
