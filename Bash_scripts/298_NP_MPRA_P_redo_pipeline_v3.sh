#!/bin/bash>
 
 
#### Rscript
 
Rscript=/software/R-4.1.0/bin/Rscript



MASTER_ROUTE=$1
enhancer_logval_Threshold=$2
mem=$3
pc=$4
queue=$5 





fdr_Threshold=$(echo '0.05')
enhancer_empirical_log_pval_Threshold=$(echo '1.3')
enhancer_pval_Threshold=$(echo '1.3')
ASE_log_pval_Threshold=$(echo '1.3')
output_dir=$MASTER_ROUTE

output="/nfs/users/nfs_m/mt19/Scripts/Wraper_scripts/298_CREATION_NP""_""$enhancer_logval_Threshold"".sh"
output_dir_2=$(echo "/lustre/scratch126/humgen/teams/soranzo/users/mt19/NEW_TWIST_LIBS/NON_PROGRAMMED/""NEW_LOOK""_""$enhancer_logval_Threshold""/")

#rm -rf $output_dir_2
#mkdir -p $output_dir_2

Log_files_path=$(echo "$output_dir_2""Log_files""/")

 
#rm -rf $Log_files_path
#mkdir -p $Log_files_path


Rscript_folder=$(echo "/nfs/users/nfs_m/mt19/Scripts/R""/")
dependencies_folder=$(echo "/lustre/scratch126/humgen/teams/soranzo/users/mt19/NEW_TWIST_LIBS/NON_PROGRAMMED/""Dependencies""/")
ALL_dB=$(echo '/lustre/scratch123/hgi/mdt1/teams/soranzo/projects/Explore_Patricks_tables/ALL_db.tsv')
dB=$(echo "$ALL_dB")
LONG_MATRIX=$(echo "/lustre/scratch126/humgen/teams/soranzo/users/mt19/NEW_TWIST_LIBS/NON_PROGRAMMED/NON_PROGRAMMED_LONG_matrix_seq_names_Plus_Flags.tsv")


FC_Threshold=$(echo "0.05")
ASE_Threshold=$(echo "0.95,1.05")
finemap_prob_Threshold=$(echo '0.1')


K562_replicates_QC_PASS=$(echo 'K562_Rep1,K562_Rep2,K562_Rep3,K562_Rep4,K562_Rep5,K562_Rep7,K562_Rep8')
CHRF_replicates_QC_PASS=$(echo 'CHRF_Rep2,CHRF_Rep3,CHRF_Rep4,CHRF_Rep5,CHRF_Rep6,CHRF_Rep7,CHRF_Rep8')
HL60_replicates_QC_PASS=$(echo 'HL60_Rep1,HL60_Rep2,HL60_Rep4,HL60_Rep5,HL60_Rep6,HL60_Rep7,HL60_Rep8')
THP1_replicates_QC_PASS=$(echo 'THP1_Rep2,THP1_Rep3,THP1_Rep4,THP1_Rep5,THP1_Rep6,THP1_Rep7,THP1_Rep8')

# CAREFUL!!!!



touch $output
echo -n "" > $output

echo "#!/bin/bash"  >> $output
 




echo "#########################################################################################################################################################################"  >> $output
echo "#########################################################################################################################################################################"  >> $output
echo "#####################################################################-----> NORMALISATION <-----###############################################################################"  >> $output
echo "#########################################################################################################################################################################"  >> $output


Rscript_data_wrangling_normalization_metrics=/nfs/users/nfs_m/mt19/Scripts/R/249_MPRA_NP_1_DW_PLUS_NOR_v3.R



type=$(echo "data_wrangling_normalization_metrics""_""$enhancer_logval_Threshold")
outfile_data_wrangling_normalization_metrics=$(echo "$Log_files_path""outfile""_""$type"".out")
touch $outfile_data_wrangling_normalization_metrics
echo -n "" > $outfile_data_wrangling_normalization_metrics
name_data_wrangling_normalization_metrics=$(echo "$type""_job")



#LONG_MATRIX=$(echo "/lustre/scratch126/humgen/teams/soranzo/users/mt19/NEW_TWIST_LIBS/NON_PROGRAMMED/NON_PROGRAMMED_LONG_matrix_seq_names_Plus_Flags_UPDATED_BY_REFERENCE.tsv")
#indir=$(echo "/lustre/scratch126/humgen/teams/soranzo/users/mt19/NEW_TWIST_LIBS/NON_PROGRAMMED/PE_15_10/")

indir=$(echo " /lustre/scratch126/humgen/teams/soranzo/users/mt19/NEW_TWIST_LIBS/NON_PROGRAMMED/NovaSeq_HAIL_MARY_results/")

step_mem=$(expr $mem \* 2)
step_pc=$(expr $pc \* 2)

echo "$mem""->""$step_mem"
echo "$pc""->""$step_pc"




echo "$step_mem"
echo "$step_pc"


echo "bsub -G hematopoiesis -o $outfile_data_wrangling_normalization_metrics -M $step_mem  -J $name_data_wrangling_normalization_metrics -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$step_mem] rusage[mem=$step_mem] span[hosts=1]\" -n$step_pc -q $queue -- \\" >> $output
echo "\"$Rscript $Rscript_data_wrangling_normalization_metrics \\" >> $output
echo "--LONG_MATRIX $LONG_MATRIX \\" >> $output
echo "--indir $indir \\" >> $output
echo "--type $type --out $output_dir\"" >> $output




Rscript_QC_graphs=/nfs/users/nfs_m/mt19/Scripts/R/251_MPRA_NP_2_QC_GRAPHS.R

type=$(echo "QC_graphs""_""$enhancer_logval_Threshold")
outfile_QC_graphs=$(echo "$Log_files_path""outfile""_""$type"".out")
touch $outfile_QC_graphs
echo -n "" > $outfile_QC_graphs
name_QC_graphs=$(echo "$type""_job")




#echo "bsub -G team151 -o $outfile_QC_graphs -M $mem  -J $name_QC_graphs -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -n$pc -q $queue -- \\" >> $output
echo "bsub -G team151 -o $outfile_QC_graphs -M $mem  -w\"done($name_data_wrangling_normalization_metrics)\" -J $name_QC_graphs -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -n$pc -q $queue -- \\" >> $output
echo "\"$Rscript $Rscript_QC_graphs \\" >> $output
echo "--type $type --out $output_dir\"" >> $output




echo "#########################################################################################################################################################################"  >> $output
echo "#########################################################################################################################################################################"  >> $output
echo "#####################################################################-----> POST QC <-----###############################################################################"  >> $output
echo "#########################################################################################################################################################################"  >> $output
echo "#########################################################################################################################################################################"  >> $output

echo "#####################################################################-----> MPRAnalyze <-----###############################################################################"  >> $output




output_dir=$MASTER_ROUTE


Rscript_QC_pass=/nfs/users/nfs_m/mt19/Scripts/R/251_MPRA_NP_3_QC_PASS.R

type=$(echo "QC_pass""_""$enhancer_logval_Threshold")
outfile_QC_pass=$(echo "$Log_files_path""outfile""_""$type"".out")
touch $outfile_QC_pass
echo -n "" > $outfile_QC_pass
name_QC_pass=$(echo "$type""_job")



echo "bsub -G team151 -o $outfile_QC_pass -M $mem -w\"done($name_data_wrangling_normalization_metrics)\"  -J $name_QC_pass -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -n$pc -q $queue -- \\" >> $output
#echo "bsub -G team151 -o $outfile_QC_pass -M $mem -J $name_QC_pass -R\"select[model==Intel_Platinum]\"  -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -n$pc -q $queue -- \\" >> $output
echo "\"$Rscript $Rscript_QC_pass \\" >> $output
echo "--K562_replicates_QC_PASS $K562_replicates_QC_PASS \\" >> $output
echo "--CHRF_replicates_QC_PASS $CHRF_replicates_QC_PASS \\" >> $output
echo "--HL60_replicates_QC_PASS $HL60_replicates_QC_PASS \\" >> $output
echo "--THP1_replicates_QC_PASS $THP1_replicates_QC_PASS \\" >> $output
echo "--type $type --out $output_dir\"" >> $output





Cell_Type_string=$(echo 'K562;''CHRF;''HL60;''THP1;')
#Cell_Type_string=$(echo 'HL60;''THP1;')


echo $Cell_Type_string

a=($(echo "$Cell_Type_string" | tr ";" '\n'))

declare -a arr

for i  in "${a[@]}"
    do
        Cell_Type=${i}
        echo "$Cell_Type"

	Rscript_QC_pass_MPRAnalyze=/nfs/users/nfs_m/mt19/Scripts/R/251_MPRA_NP_3_5_MPRAnalyze_partII.R

	type=$(echo "QC_pass_MPRAnalyze""_""$Cell_Type""_""$enhancer_logval_Threshold")
	outfile_QC_pass_MPRAnalyze=$(echo "$Log_files_path""outfile""_""$type"".out")
	touch $outfile_QC_pass_MPRAnalyze
	echo -n "" > $outfile_QC_pass_MPRAnalyze
	name_QC_pass_MPRAnalyze=$(echo "$type""_job")



	echo "bsub -G team151 -o $outfile_QC_pass_MPRAnalyze -M $mem -w\"done($name_QC_pass)\"  -J $name_QC_pass_MPRAnalyze -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -n$pc -q $queue -- \\" >> $output
#	echo "bsub -G team151 -o $outfile_QC_pass_MPRAnalyze -M $mem  -J $name_QC_pass_MPRAnalyze -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -n$pc -q $queue -- \\" >> $output
	echo "\"$Rscript $Rscript_QC_pass_MPRAnalyze \\" >> $output
	echo "--Cell_Type $Cell_Type \\" >> $output
	echo "--type $type --out $output_dir\"" >> $output

	QC_pass_MPRAnalyze_string=$(echo "&& done($name_QC_pass_MPRAnalyze)")
        echo "->>>$QC_pass_MPRAnalyze_string"
	arr[${#arr[@]}]="$QC_pass_MPRAnalyze_string"
	
done
 



 echo "#####################################################################-----> parameters_recalculation.R  <-----###############################################################################"  >> $output



Rscript_Parameter_recalculation=/nfs/users/nfs_m/mt19/Scripts/R/251_MPRA_4_QC_parameters_recalculation_v2_NP.R

type=$(echo "Parameter_recalculation""_""$enhancer_logval_Threshold")
outfile_Parameter_recalculation=$(echo "$Log_files_path""outfile""_""$type"".out")
touch $outfile_Parameter_recalculation
echo -n "" > $outfile_Parameter_recalculation
name_Parameter_recalculation=$(echo "$type""_job")

 


echo "bsub -G team151 -o $outfile_Parameter_recalculation -M $mem -w\"done($name_data_wrangling_normalization_metrics)\"  -J $name_Parameter_recalculation -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -n$pc -q $queue -- \\" >> $output
#echo "bsub -G team151 -o $outfile_Parameter_recalculation -M $mem  -J $name_Parameter_recalculation -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -n$pc -q $queue -- \\" >> $output
echo "\"$Rscript $Rscript_Parameter_recalculation \\" >> $output
echo "--K562_replicates_QC_PASS $K562_replicates_QC_PASS \\" >> $output
echo "--CHRF_replicates_QC_PASS $CHRF_replicates_QC_PASS \\" >> $output
echo "--HL60_replicates_QC_PASS $HL60_replicates_QC_PASS \\" >> $output
echo "--THP1_replicates_QC_PASS $THP1_replicates_QC_PASS \\" >> $output
echo "--type $type --out $output_dir --out2 $output_dir_2\"" >> $output

first_half_done_string=$(echo "\"""done($name_Parameter_recalculation) ")
#echo "$first_half_done_string"

done_string=$(echo ""${arr[@]}"""\"")
#echo "$done_string"

complete_done_string=$(echo "$first_half_done_string""$done_string")

echo "$complete_done_string"




#  # echo "#########################################################################################################################################################################"  >> $output
# # echo "#####################################################################-----> Volcanos overview <-----###############################################################################"  >> $output






Rscript_volcano_overview=$(echo "/nfs/users/nfs_m/mt19/Scripts/R/251_MPRA_5_Volcanos_and_ACTIVE_tile_definition_v3_NP.R")


type=$(echo "volcano_overview""_""$enhancer_logval_Threshold")
outfile_volcano_overview=$(echo "$Log_files_path""outfile""_""$type"".out")
touch $outfile_volcano_overview
echo -n "" > $outfile_volcano_overview
name_volcano_overview=$(echo "$type""_job")

enhancer_result_K562=$(echo "$output_dir""MPRAnalyze_enhancer_results_K562.txt")
ASE_result_K562=$(echo "$output_dir""MPRAnalyze_ASE_results_ASE_K562.txt")

enhancer_result_CHRF=$(echo "$output_dir""MPRAnalyze_enhancer_results_CHRF.txt")
ASE_result_CHRF=$(echo "$output_dir""MPRAnalyze_ASE_results_ASE_CHRF.txt")

enhancer_result_HL60=$(echo "$output_dir""MPRAnalyze_enhancer_results_HL60.txt")
ASE_result_HL60=$(echo "$output_dir""MPRAnalyze_ASE_results_ASE_HL60.txt")

enhancer_result_THP1=$(echo "$output_dir""MPRAnalyze_enhancer_results_THP1.txt")
ASE_result_THP1=$(echo "$output_dir""MPRAnalyze_ASE_results_ASE_THP1.txt")

VEP_CSQ=$(echo '/lustre/scratch126/humgen/teams/soranzo/users/ALL_dB/Variant_csv_tables/VEP_consequence_graphs.csv')


echo "bsub -G team151 -o $outfile_volcano_overview -M $mem -w$complete_done_string  -J $name_volcano_overview -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -n$pc -q $queue -- \\" >> $output
#echo "bsub -G team151 -o $outfile_volcano_overview -M $mem -w\"($name_Parameter_recalculation)\"  -J $name_volcano_overview -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -n$pc -q $queue -- \\" >> $output
#echo "bsub -G team151 -o $outfile_volcano_overview -M $mem -J $name_volcano_overview -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -n$pc -q $queue -- \\" >> $output
echo "\"$Rscript $Rscript_volcano_overview \\" >> $output
echo "--fdr_Threshold $fdr_Threshold \\" >> $output
echo "--FC_Threshold $FC_Threshold \\" >> $output
echo "--ASE_Threshold $ASE_Threshold \\" >> $output
echo "--enhancer_logval_Threshold $enhancer_logval_Threshold \\" >> $output
echo "--ASE_log_pval_Threshold $ASE_log_pval_Threshold \\" >> $output
echo "--enhancer_result_K562 $enhancer_result_K562 \\" >> $output
echo "--ASE_result_K562 $ASE_result_K562 \\" >> $output
echo "--enhancer_result_CHRF $enhancer_result_CHRF \\" >> $output
echo "--ASE_result_CHRF $ASE_result_CHRF \\" >> $output
echo "--enhancer_result_HL60 $enhancer_result_HL60 \\" >> $output
echo "--ASE_result_HL60 $ASE_result_HL60 \\" >> $output
echo "--enhancer_result_THP1 $enhancer_result_THP1 \\" >> $output
echo "--ASE_result_THP1 $ASE_result_THP1 \\" >> $output
echo "--VEP_CSQ $VEP_CSQ \\" >> $output
echo "--type $type --out $output_dir --out2 $output_dir_2\"" >> $output


echo "#####################################################################-----> VJ_directionality  <-----###############################################################################"  >> $output


Rscript_VJ_directionality=$(echo "$Rscript_folder""251_MPRA_12_VJ_check_directionality_v3.R")

type=$(echo "VJ_directionality""_""$enhancer_logval_Threshold")
outfile_VJ_directionality=$(echo "$Log_files_path""outfile""_""$type"".out")
touch $outfile_VJ_directionality
echo -n "" > $outfile_VJ_directionality
name_VJ_directionality=$(echo "$type""_job")

Sankaran_MPRA=$(echo $dependencies_folder"MPRA_Sankaran.csv")
MPRA_Real_tile_QC2_PASS=$(echo "$output_dir_2""MPRA_Real_Tile_QC2_PASS.rds")



step_mem=$(echo  "4000")
step_pc=$(echo "1")


echo "bsub -G team151 -o $outfile_VJ_directionality -M $step_mem -w\"done($name_volcano_overview)\" -J $name_VJ_directionality -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$step_mem] rusage[mem=$step_mem] span[hosts=1]\" -n$step_pc -q $queue -- \\" >> $output
#echo "bsub -G team151 -o $outfile_VJ_directionality -M $step_mem -J $name_VJ_directionality -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$step_mem] rusage[mem=$step_mem] span[hosts=1]\" -n$step_pc -q $queue -- \\" >> $output
echo "\"$Rscript $Rscript_VJ_directionality \\" >> $output
echo "--Sankaran_MPRA $Sankaran_MPRA \\" >> $output
echo "--MPRA_Real_tile_QC2_PASS $MPRA_Real_tile_QC2_PASS \\" >> $output
echo "--type $type --out $output_dir --out2 $output_dir_2\"" >> $output





echo "#####################################################################-----> POST_QC_PLOTS  <-----###############################################################################"  >> $output


Rscript_POST_QC_graphs=$(echo "$Rscript_folder""251_MPRA_13_POST_QC_PLOTS_v2_NP.R")

type=$(echo "POST_QC_graphs""_""$enhancer_logval_Threshold")
outfile_POST_QC_graphs=$(echo "$Log_files_path""outfile""_""$type"".out")
touch $outfile_POST_QC_graphs
echo -n "" > $outfile_POST_QC_graphs
name_POST_QC_graphs=$(echo "$type""_job")

MPRA_Real_tile_QC2_PASS=$(echo "$output_dir_2""MPRA_Real_Tile_QC2_PASS.rds")


step_mem=$(echo  "4000")
step_pc=$(echo "1")

echo "bsub -G team151 -o $outfile_POST_QC_graphs -M $step_mem -w\"done($name_volcano_overview)\" -J $name_POST_QC_graphs -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$step_mem] rusage[mem=$step_mem] span[hosts=1]\" -n$step_pc -q $queue -- \\" >> $output
#echo "bsub -G team151 -o $outfile_POST_QC_graphs -M $step_mem -J $name_POST_QC_graphs -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$step_mem] rusage[mem=$step_mem] span[hosts=1]\" -n$step_pc -q $queue -- \\" >> $output
echo "\"$Rscript $Rscript_POST_QC_graphs \\" >> $output
echo "--MPRA_Real_tile_QC2_PASS $MPRA_Real_tile_QC2_PASS \\" >> $output
echo "--type $type --out $output_dir --out2 $output_dir_2\"" >> $output



echo "#####################################################################-----> KTP_collapse <-----###############################################################################"  >> $output

Rscript_KTP_collapse=$(echo "$Rscript_folder""251_MPRA_6_KTP_collapse_v2.R")

type=$(echo "KTP_collapse""_""$enhancer_logval_Threshold")
outfile_KTP_collapse=$(echo "$Log_files_path""outfile""_""$type"".out")
touch $outfile_KTP_collapse
echo -n "" > $outfile_KTP_collapse
name_KTP_collapse=$(echo "$type""_""job")


ACTIVE_TILES=$(echo "$output_dir_2""ACTIVE_TILES.txt")
MPRA_Real_tile_QC2_PASS=$(echo "$output_dir_2""MPRA_Real_Tile_QC2_PASS.rds")
CSQ_colors=$(echo $dependencies_folder"df_CSQ_colors.rds")

step_mem=$(echo  "4000")
step_pc=$(echo "1")

echo "$step_mem"
echo "$step_pc"

echo "bsub -G team151 -o $outfile_KTP_collapse -M $step_mem -w\"done($name_volcano_overview)\" -J $name_KTP_collapse -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$step_mem] rusage[mem=$step_mem] span[hosts=1]\" -n$step_pc -q $queue -- \\" >> $output
#echo "bsub -G team151 -o $outfile_KTP_collapse -M $step_mem -J $name_KTP_collapse -R\"select[mem>=$step_mem] rusage[mem=$step_mem] span[hosts=1]\" -n$step_pc -q $queue -- \\" >> $output
echo "\"$Rscript $Rscript_KTP_collapse \\" >> $output
echo "--MPRA_Real_tile_QC2_PASS $MPRA_Real_tile_QC2_PASS \\" >> $output
echo "--ACTIVE_TILES $ACTIVE_TILES \\" >> $output
echo "--CSQ_colors $CSQ_colors \\" >> $output
echo "--type $type --out $output_dir --out2 $output_dir_2\"" >> $output

echo "#########################################################################################################################################################################"  >> $output
echo "#####################################################################-----> Per variant graphs <-----###############################################################################"  >> $output


Rscript_dot_plot_variant=$(echo "$Rscript_folder""251_MPRA_16_Per_variant_dot_plot_v2_NP.R")

type=$(echo "dot_plot_variant""_""$enhancer_logval_Threshold")
outfile_dot_plot_variant=$(echo "$Log_files_path""outfile""_""$type"".out")
touch $outfile_dot_plot_variant
echo -n "" > $outfile_dot_plot_variant
name_dot_plot_variant=$(echo "$type""_job")

MPRA_Real_tile_QC2_PASS=$(echo "$output_dir_2""MPRA_Real_Tile_QC2_PASS.rds")
echo "$LONG_MATRIX"

step_mem=$(echo  "4000")
step_pc=$(echo "1")

echo "$step_mem"
echo "$step_pc"

echo "bsub -G team151 -o $outfile_dot_plot_variant -M $step_mem -w\"done($name_KTP_collapse)\" -J $name_dot_plot_variant -R\"select[mem>=$step_mem] rusage[mem=$step_mem] span[hosts=1]\" -n$step_pc -q $queue -- \\" >> $output
#echo "bsub -G team151 -o $outfile_dot_plot_variant -M $step_mem -J $name_dot_plot_variant -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$step_mem] rusage[mem=$step_mem] span[hosts=1]\" -n$step_pc -q $queue -- \\" >> $output
echo "\"$Rscript $Rscript_dot_plot_variant \\" >> $output
echo "--LONG_MATRIX $LONG_MATRIX \\" >> $output
echo "--ALL_dB $ALL_dB \\" >> $output
echo "--MPRA_Real_tile_QC2_PASS $MPRA_Real_tile_QC2_PASS \\" >> $output
echo "--type $type --out $output_dir --out2 $output_dir_2\"" >> $output


echo "#####################################################################-----> UPSETR <-----###############################################################################"  >> $output
echo "#####################################################################-----> Tier Printer <-----###############################################################################"  >> $output



Rscript_UPSETR=$(echo "$Rscript_folder""251_MPRA_8_UPSETR_PLOTS_v3.R")

type=$(echo "UPSETR""_""$enhancer_logval_Threshold")
outfile_UPSETR=$(echo "$Log_files_path""outfile""_""$type"".out")
touch $outfile_UPSETR
echo -n "" > $outfile_UPSETR
name_UPSETR=$(echo "$type""_""job")

step_mem=$(echo  "4000")
step_pc=$(echo "1")

echo "$step_mem"
echo "$step_pc"

ACTIVE_TILES=$(echo "$output_dir_2""ACTIVE_TILES.txt")
MPRA_Real_tile_QC2_PASS=$(echo "$output_dir_2""MPRA_Real_Tile_QC2_PASS.rds")
CSQ_colors=$(echo $dependencies_folder"df_CSQ_colors.rds")


echo "bsub -G team151 -o $outfile_UPSETR -M $step_mem -w\"done($name_KTP_collapse)\" -J $name_UPSETR -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$step_mem] rusage[mem=$step_mem] span[hosts=1]\" -n$step_pc -q $queue -- \\" >> $output
#echo "bsub -G team151 -o $outfile_UPSETR -M $step_mem  -J $name_UPSETR -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$step_mem] rusage[mem=$step_mem] span[hosts=1]\" -n$step_pc -q $queue -- \\" >> $output
echo "\"$Rscript $Rscript_UPSETR \\" >> $output
echo "--CSQ_colors $CSQ_colors \\" >> $output
echo "--type $type --out $output_dir --out2 $output_dir_2\"" >> $output


echo "#####################################################################-----> CUMULATIVE_CURVES <-----###############################################################################"  >> $output

Rscript_CUMULATIVE_CURVES=$(echo "$Rscript_folder""251_MPRA_7_5_CUMMULATIVE_CURVES_v2_NP.R")

type=$(echo "CUMULATIVE_CURVES""_""$enhancer_logval_Threshold")
outfile_CUMULATIVE_CURVES=$(echo "$Log_files_path""outfile""_""$type"".out")
touch $outfile_CUMULATIVE_CURVES
echo -n "" > $outfile_CUMULATIVE_CURVES
name_CUMULATIVE_CURVES=$(echo "$type""_""job")

ACTIVE_TILES=$(echo "$output_dir_2""ACTIVE_TILES_for_cummulative.txt")
CSQ_colors=$(echo $dependencies_folder"df_CSQ_colors.rds")


echo "bsub -G team151 -o $outfile_CUMULATIVE_CURVES -M $mem -w\"done($name_KTP_collapse)\" -J $name_CUMULATIVE_CURVES -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -n$pc -q $queue -- \\" >> $output
#echo "bsub -G team151 -o $outfile_CUMULATIVE_CURVES -M $mem -J $name_CUMULATIVE_CURVES -R\"select[model==Intel_Platinum]\"  -R\"select[mem>=$mem] rusage[mem=$mem] span[hosts=1]\" -n$pc -q $queue -- \\" >> $output
echo "\"$Rscript $Rscript_CUMULATIVE_CURVES \\" >> $output
echo "--CSQ_colors $CSQ_colors \\" >> $output
echo "--ACTIVE_TILES $ACTIVE_TILES \\" >> $output
echo "--type $type --out $output_dir --out2 $output_dir_2\"" >> $output




echo "#####################################################################-----> CUMULATIVE_CURVES_PER_CT_Prioritisation_bins <-----###############################################################################"  >> $output

Rscript_CT_Prioritisation_bins_CUMULATIVE_CURVES=$(echo "$Rscript_folder""251_MPRA_23_CUMMULATIVE_CURVES_Prioritisation_bins_NP.R")

type=$(echo "CT_Prioritisation_bins_CUMULATIVE_CURVES""_""$enhancer_logval_Threshold")
outfile_CT_Prioritisation_bins_CUMULATIVE_CURVES=$(echo "$Log_files_path""outfile""_""$type"".out")
touch $outfile_CT_Prioritisation_bins_CUMULATIVE_CURVES
echo -n "" > $outfile_CT_Prioritisation_bins_CUMULATIVE_CURVES
name_CT_Prioritisation_bins_CUMULATIVE_CURVES=$(echo "$type""_""job")

CUMMULATIVE_CLASSES=$(echo "$output_dir_2""cummulative_plots/cummulative_classification_with_CTRLS.rds")
CSQ_colors=$(echo $dependencies_folder"df_CSQ_colors.rds")
VAR_Prioritization_dB=$(echo "/lustre/scratch126/humgen/teams/soranzo/users/mt19/V2F_paper/MPRA_E_Plus_ASE_ACTIVE_AT_LEAST_1_genIE_Tier_del_0.1/ER_Labelling_Initial_Selection.rds")

step_mem=$(echo  "4000")
step_pc=$(echo "1")

echo "$step_mem"
echo "$step_pc"
echo "bsub -G team151 -o $outfile_CT_Prioritisation_bins_CUMULATIVE_CURVES -M $step_mem -w\"done($name_CUMULATIVE_CURVES)\" -J $name_CT_Prioritisation_bins_CUMULATIVE_CURVES -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$step_mem] rusage[mem=$step_mem] span[hosts=1]\" -n$step_pc -q $queue -- \\" >> $output
#echo "bsub -G team151 -o $outfile_CT_Prioritisation_bins_CUMULATIVE_CURVES -M $step_mem -J $name_CT_Prioritisation_bins_CUMULATIVE_CURVES -R\"select[model==Intel_Platinum]\"  -R\"select[mem>=$step_mem] rusage[mem=$step_mem] span[hosts=1]\" -n$step_pc -q $queue -- \\" >> $output
echo "\"$Rscript $Rscript_CT_Prioritisation_bins_CUMULATIVE_CURVES \\" >> $output
echo "--CUMMULATIVE_CLASSES $CUMMULATIVE_CLASSES \\" >> $output
echo "--VAR_Prioritization_dB $VAR_Prioritization_dB \\" >> $output
echo "--type $type --out $output_dir --out2 $output_dir_2\"" >> $output

echo "#####################################################################-----> CUMULATIVE_CURVES_PER_CT_PP_bins <-----###############################################################################"  >> $output

Rscript_CT_PP_bins_CUMULATIVE_CURVES=$(echo "$Rscript_folder""251_MPRA_24_CUMMULATIVE_CURVES_PP_bins_NP.R")

type=$(echo "CT_PP_bins_CUMULATIVE_CURVES""_""$enhancer_logval_Threshold")
outfile_CT_PP_bins_CUMULATIVE_CURVES=$(echo "$Log_files_path""outfile""_""$type"".out")
touch $outfile_CT_PP_bins_CUMULATIVE_CURVES
echo -n "" > $outfile_CT_PP_bins_CUMULATIVE_CURVES
name_CT_PP_bins_CUMULATIVE_CURVES=$(echo "$type""_""job")

CUMMULATIVE_CLASSES=$(echo "$output_dir_2""cummulative_plots/cummulative_classification_with_CTRLS.rds")
CSQ_colors=$(echo $dependencies_folder"df_CSQ_colors.rds")
ALL_dB=$(echo '/lustre/scratch123/hgi/mdt1/teams/soranzo/projects/Explore_Patricks_tables/ALL_db.tsv')
VAR_Prioritization_dB=$(echo "/lustre/scratch126/humgen/teams/soranzo/users/mt19/V2F_paper/MPRA_E_Plus_ASE_ACTIVE_AT_LEAST_1_genIE_Tier_del_0.1/ER_Labelling_Initial_Selection_with_CSQ_labels.rds")

step_mem=$(echo  "4000")
step_pc=$(echo "1")

echo "$step_mem"
echo "$step_pc"
echo "bsub -G team151 -o $outfile_CT_PP_bins_CUMULATIVE_CURVES -M $step_mem -w\"done($name_CUMULATIVE_CURVES)\" -J $name_CT_PP_bins_CUMULATIVE_CURVES -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$step_mem] rusage[mem=$step_mem] span[hosts=1]\" -n$step_pc -q $queue -- \\" >> $output
#echo "bsub -G team151 -o $outfile_CT_PP_bins_CUMULATIVE_CURVES -M $step_mem -J $name_CT_PP_bins_CUMULATIVE_CURVES -R\"select[model==Intel_Platinum]\"  -R\"select[mem>=$step_mem] rusage[mem=$step_mem] span[hosts=1]\" -n$step_pc -q $queue -- \\" >> $output
echo "\"$Rscript $Rscript_CT_PP_bins_CUMULATIVE_CURVES \\" >> $output
echo "--CUMMULATIVE_CLASSES $CUMMULATIVE_CLASSES \\" >> $output
echo "--VAR_Prioritization_dB $VAR_Prioritization_dB \\" >> $output
echo "--ALL_dB $ALL_dB \\" >> $output
echo "--type $type --out $output_dir --out2 $output_dir_2\"" >> $output





echo "#####################################################################-----> CUMULATIVE_CURVES_PER_CT_LABEL_2 <-----###############################################################################"  >> $output

Rscript_CT_Label_2_CUMULATIVE_CURVES=$(echo "$Rscript_folder""251_MPRA_17_CUMMULATIVE_CURVES_Label2.R")

type=$(echo "CT_Label_2_CUMULATIVE_CURVES""_""$enhancer_logval_Threshold")
outfile_CT_Label_2_CUMULATIVE_CURVES=$(echo "$Log_files_path""outfile""_""$type"".out")
touch $outfile_CT_Label_2_CUMULATIVE_CURVES
echo -n "" > $outfile_CT_Label_2_CUMULATIVE_CURVES
name_CT_Label_2_CUMULATIVE_CURVES=$(echo "$type""_""job")

CUMMULATIVE_CLASSES=$(echo "$output_dir_2""cummulative_plots/cummulative_classification_with_CTRLS.rds")
CSQ_colors=$(echo $dependencies_folder"df_CSQ_colors.rds")

step_mem=$(echo  "4000")
step_pc=$(echo "1")

echo "$step_mem"
echo "$step_pc"
echo "bsub -G team151 -o $outfile_CT_Label_2_CUMULATIVE_CURVES -M $step_mem -w\"done($name_CUMULATIVE_CURVES)\" -J $name_CT_Label_2_CUMULATIVE_CURVES -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$step_mem] rusage[mem=$step_mem] span[hosts=1]\" -n$step_pc -q $queue -- \\" >> $output
#echo "bsub -G team151 -o $outfile_CT_Label_2_CUMULATIVE_CURVES -M $step_mem -J $name_CT_Label_2_CUMULATIVE_CURVES -R\"select[model==Intel_Platinum]\"  -R\"select[mem>=$step_mem] rusage[mem=$step_mem] span[hosts=1]\" -n$step_pc -q $queue -- \\" >> $output
echo "\"$Rscript $Rscript_CT_Label_2_CUMULATIVE_CURVES \\" >> $output
echo "--CUMMULATIVE_CLASSES $CUMMULATIVE_CLASSES \\" >> $output
echo "--type $type --out $output_dir --out2 $output_dir_2\"" >> $output
 
echo "#####################################################################-----> CUMMULATIVE_CURVES_VEP_DEF_LABELS <-----###############################################################################"  >> $output

Rscript_CT_VEP_DEF_LABELS_CUMULATIVE_CURVES=$(echo "$Rscript_folder""251_MPRA_18_CUMMULATIVE_CURVES_VEP_DEF_LABELS.R")

type=$(echo "CT_VEP_DEF_LABELS_CUMULATIVE_CURVES""_""$enhancer_logval_Threshold")
outfile_CT_VEP_DEF_LABELS_CUMULATIVE_CURVES=$(echo "$Log_files_path""outfile""_""$type"".out")
touch $outfile_CT_VEP_DEF_LABELS_CUMULATIVE_CURVES
echo -n "" > $outfile_CT_VEP_DEF_LABELS_CUMULATIVE_CURVES
name_CT_VEP_DEF_LABELS_CUMULATIVE_CURVES=$(echo "$type""_""job")

CUMMULATIVE_CLASSES=$(echo "$output_dir_2""cummulative_plots/cummulative_classification_with_CTRLS.rds")
CSQ_colors=$(echo $dependencies_folder"df_CSQ_colors.rds")


step_mem=$(echo  "4000")
step_pc=$(echo "1")

echo "$step_mem"
echo "$step_pc"


echo "bsub -G team151 -o $outfile_CT_VEP_DEF_LABELS_CUMULATIVE_CURVES -M $step_mem -w\"done($name_CUMULATIVE_CURVES)\" -J $name_CT_VEP_DEF_LABELS_CUMULATIVE_CURVES -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$step_mem] rusage[mem=$step_mem] span[hosts=1]\" -n$step_pc -q $queue -- \\" >> $output
#echo "bsub -G team151 -o $outfile_CT_VEP_DEF_LABELS_CUMULATIVE_CURVES -M $step_mem -J $name_CT_VEP_DEF_LABELS_CUMULATIVE_CURVES -R\"select[model==Intel_Platinum]\"  -R\"select[mem>=$step_mem] rusage[mem=$step_mem] span[hosts=1]\" -n$step_pc -q $queue -- \\" >> $output
echo "\"$Rscript $Rscript_CT_VEP_DEF_LABELS_CUMULATIVE_CURVES \\" >> $output
echo "--CUMMULATIVE_CLASSES $CUMMULATIVE_CLASSES \\" >> $output
echo "--CSQ_colors $CSQ_colors \\" >> $output
echo "--type $type --out $output_dir --out2 $output_dir_2\"" >> $output

echo "#####################################################################-----> CUMMULATIVE_CURVES_Factor4 <-----###############################################################################"  >> $output

Rscript_CT_Factor4_CUMULATIVE_CURVES=$(echo "$Rscript_folder""251_MPRA_20_CUMMULATIVE_CURVES_Factor4.R")

type=$(echo "CT_Factor4_CUMULATIVE_CURVES""_""$enhancer_logval_Threshold")
outfile_CT_Factor4_CUMULATIVE_CURVES=$(echo "$Log_files_path""outfile""_""$type"".out")
touch $outfile_CT_Factor4_CUMULATIVE_CURVES
echo -n "" > $outfile_CT_Factor4_CUMULATIVE_CURVES
name_CT_Factor4_CUMULATIVE_CURVES=$(echo "$type""_""job")

CUMMULATIVE_CLASSES=$(echo "$output_dir_2""cummulative_plots/cummulative_classification_with_CTRLS.rds")
CSQ_colors=$(echo $dependencies_folder"df_CSQ_colors.rds")
MPRA_Real_tile_QC2_PASS=$(echo "$output_dir_2""MPRA_Real_Tile_QC2_PASS.rds")

step_mem=$(echo  "4000")
step_pc=$(echo "1")
echo "$step_mem"
echo "$step_pc"


echo "bsub -G team151 -o $outfile_CT_Factor4_CUMULATIVE_CURVES -M $step_mem -w\"done($name_CUMULATIVE_CURVES)\" -J $name_CT_Factor4_CUMULATIVE_CURVES -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$step_mem] rusage[mem=$step_mem] span[hosts=1]\" -n$step_pc -q $queue -- \\" >> $output
#echo "bsub -G team151 -o $outfile_CT_Factor4_CUMULATIVE_CURVES -M $step_mem -J $name_CT_Factor4_CUMULATIVE_CURVES -R\"select[model==Intel_Platinum]\"  -R\"select[mem>=$step_mem] rusage[mem=$step_mem] span[hosts=1]\" -n$step_pc -q $queue -- \\" >> $output
echo "\"$Rscript $Rscript_CT_Factor4_CUMULATIVE_CURVES \\" >> $output
echo "--CUMMULATIVE_CLASSES $CUMMULATIVE_CLASSES \\" >> $output
echo "--CSQ_colors $CSQ_colors \\" >> $output
echo "--type $type --out $output_dir --out2 $output_dir_2\"" >> $output

echo "#####################################################################-----> CUMMULATIVE_CURVES_LINEAGES <-----###############################################################################"  >> $output

Rscript_CT_LINEAGES_CUMULATIVE_CURVES=$(echo "$Rscript_folder""251_MPRA_19_CUMMULATIVE_CURVES_Lineages.R")

type=$(echo "CT_LINEAGES_CUMULATIVE_CURVES""_""$enhancer_logval_Threshold")
outfile_CT_LINEAGES_CUMULATIVE_CURVES=$(echo "$Log_files_path""outfile""_""$type"".out")
touch $outfile_CT_LINEAGES_CUMULATIVE_CURVES
echo -n "" > $outfile_CT_LINEAGES_CUMULATIVE_CURVES
name_CT_LINEAGES_CUMULATIVE_CURVES=$(echo "$type""_""job")



echo "--1-->""$ALL_dB"

CUMMULATIVE_CLASSES=$(echo "$output_dir_2""cummulative_plots/cummulative_classification_with_CTRLS.rds")


echo "--2-->""$dB"

TOME_correspondence=$(echo $dependencies_folder"Correspondence_phenotype_TOME.txt")


step_mem=$(echo  "4000")
step_pc=$(echo "1")


echo "$step_mem"
echo "$step_pc"
echo "bsub -G team151 -o $outfile_CT_LINEAGES_CUMULATIVE_CURVES -M $step_mem -w\"done($name_CUMULATIVE_CURVES)\" -J $name_CT_LINEAGES_CUMULATIVE_CURVES -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$step_mem] rusage[mem=$step_mem] span[hosts=1]\" -n$step_pc -q $queue -- \\" >> $output
#echo "bsub -G team151 -o $outfile_CT_LINEAGES_CUMULATIVE_CURVES -M $step_mem -J $name_CT_LINEAGES_CUMULATIVE_CURVES -R\"select[model==Intel_Platinum]\"  -R\"select[mem>=$step_mem] rusage[mem=$step_mem] span[hosts=1]\" -n$step_pc -q $queue -- \\" >> $output
echo "\"$Rscript $Rscript_CT_LINEAGES_CUMULATIVE_CURVES \\" >> $output
echo "--dB $dB \\" >> $output
echo "--TOME_correspondence $TOME_correspondence \\" >> $output
echo "--finemap_prob_Threshold $finemap_prob_Threshold \\" >> $output
echo "--CUMMULATIVE_CLASSES $CUMMULATIVE_CLASSES \\" >> $output
echo "--type $type --out $output_dir --out2 $output_dir_2\"" >> $output




echo "#####################################################################-----> EFFECT_SIZE_STATS_and_LM <-----###############################################################################"  >> $output

Rscript_EFFECT_SIZE_STATS_and_LM=$(echo "$Rscript_folder""251_MPRA_10_Effect_sizes_v2_NP.R")

type=$(echo "EFFECT_SIZE_STATS_and_LM""_""$enhancer_logval_Threshold")
outfile_EFFECT_SIZE_STATS_and_LM=$(echo "$Log_files_path""outfile""_""$type"".out")
touch $outfile_EFFECT_SIZE_STATS_and_LM
echo -n "" > $outfile_EFFECT_SIZE_STATS_and_LM
name_EFFECT_SIZE_STATS_and_LM=$(echo "$type""_""job")



TOME_correspondence=$(echo $dependencies_folder"Correspondence_phenotype_TOME.txt")
MPRA_Real_tile_QC2_PASS=$(echo "$output_dir_2""MPRA_Real_Tile_QC2_PASS.rds")
CUMMULATIVE_CLASSES=$(echo "$output_dir_2""LINEAGE_CUMULATIVE_CURVES/cummulative_classification_ONLY_ASSAYED_with_Lineage_classification.rds")

step_mem=$(echo  "8000")
step_pc=$(echo "2")

echo "$step_mem"
echo "$step_pc"

echo "bsub -G team151 -o $outfile_EFFECT_SIZE_STATS_and_LM -M $step_mem -w\"done($name_CT_LINEAGES_CUMULATIVE_CURVES)\" -J $name_EFFECT_SIZE_STATS_and_LM -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$step_mem] rusage[mem=$step_mem] span[hosts=1]\" -n$step_pc -q $queue -- \\" >> $output
#echo "bsub -G team151 -o $outfile_EFFECT_SIZE_STATS_and_LM -M $step_mem -J $name_EFFECT_SIZE_STATS_and_LM -R\"select[model==Intel_Platinum]\"  -R\"select[mem>=$step_mem] rusage[mem=$step_mem] span[hosts=1]\" -n$step_pc -q $queue -- \\" >> $output
echo "\"$Rscript $Rscript_EFFECT_SIZE_STATS_and_LM \\" >> $output
echo "--CUMMULATIVE_CLASSES $CUMMULATIVE_CLASSES \\" >> $output
echo "--dB $dB \\" >> $output
echo "--TOME_correspondence $TOME_correspondence \\" >> $output
echo "--MPRA_Real_tile_QC2_PASS $MPRA_Real_tile_QC2_PASS \\" >> $output
echo "--finemap_prob_Threshold $finemap_prob_Threshold \\" >> $output
echo "--type $type --out $output_dir --out2 $output_dir_2\"" >> $output


#exit

echo "#####################################################################-----> EFFECT_SIZE_PLOTS <-----###############################################################################"  >> $output

Rscript_EFFECT_SIZE_PLOTS=$(echo "$Rscript_folder""251_MPRA_11_Effect_sizes_selected_plots_v3.R")

type=$(echo "EFFECT_SIZE_PLOTS""_""$enhancer_logval_Threshold")
outfile_EFFECT_SIZE_PLOTS=$(echo "$Log_files_path""outfile""_""$type"".out")
touch $outfile_EFFECT_SIZE_PLOTS
echo -n "" > $outfile_EFFECT_SIZE_PLOTS
name_EFFECT_SIZE_PLOTS=$(echo "$type""_""job")

step_mem=$(echo  "4000")
step_pc=$(echo "1")

echo "$step_mem"
echo "$step_pc"

echo "bsub -G team151 -o $outfile_EFFECT_SIZE_PLOTS -M $step_mem -w\"done($name_EFFECT_SIZE_STATS_and_LM)\" -J $name_EFFECT_SIZE_PLOTS -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$step_mem] rusage[mem=$step_mem] span[hosts=1]\" -n$step_pc -q $queue -- \\" >> $output
#echo "bsub -G team151 -o $outfile_EFFECT_SIZE_PLOTS -M $step_mem -J $name_EFFECT_SIZE_PLOTS -R\"select[model==Intel_Platinum]\"  -R\"select[mem>=$step_mem] rusage[mem=$step_mem] span[hosts=1]\" -n$step_pc -q $queue -- \\" >> $output
echo "\"$Rscript $Rscript_EFFECT_SIZE_PLOTS \\" >> $output
echo "--type $type --out $output_dir --out2 $output_dir_2\"" >> $output


echo "#####################################################################-----> AS_printing  <-----###############################################################################"  >> $output


Rscript_CT_AS_printing_CUMULATIVE_CURVES=$(echo "$Rscript_folder""251_MPRA_25_CUMMULATIVE_CURVES_AS_PRINTER_NP.R")

type=$(echo "CT_AS_printing_CUMULATIVE_CURVES""_""$enhancer_logval_Threshold")
outfile_CT_AS_printing_CUMULATIVE_CURVES=$(echo "$Log_files_path""outfile""_""$type"".out")
touch $outfile_CT_AS_printing_CUMULATIVE_CURVES
echo -n "" > $outfile_CT_AS_printing_CUMULATIVE_CURVES
name_CT_AS_printing_CUMULATIVE_CURVES=$(echo "$type""_""job")

CUMMULATIVE_CLASSES=$(echo "$output_dir_2""cummulative_plots/cummulative_classification_with_CTRLS.rds")
CSQ_colors=$(echo $dependencies_folder"df_CSQ_colors.rds")
ALL_dB=$(echo '/lustre/scratch123/hgi/mdt1/teams/soranzo/projects/Explore_Patricks_tables/ALL_db.tsv')
Supp4_Table_CURATED_PLUS_PHENOTYPES=$(echo "/lustre/scratch123/hgi/mdt1/teams/soranzo/projects/VARIANT_INTERPRETATION/OVERHAUL_INTERVAL/FINAL_RESULTS/Fig4_pannels/""Supp_Table_4_CURATED_Plus_phenotypes.rds")
VAR_Prioritization_dB=$(echo "/lustre/scratch126/humgen/teams/soranzo/users/mt19/V2F_paper/MPRA_E_Plus_ASE_ACTIVE_AT_LEAST_1_genIE_Tier_del_0.1/ER_Labelling_Initial_Selection_with_CSQ_labels.rds")

step_mem=$(echo  "4000")
step_pc=$(echo "1")

echo "$step_mem"
echo "$step_pc"

echo "bsub -G team151 -o $outfile_CT_AS_printing_CUMULATIVE_CURVES -M $step_mem -w\"done($name_CUMULATIVE_CURVES)\" -J $name_CT_AS_printing_CUMULATIVE_CURVES -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$step_mem] rusage[mem=$step_mem] span[hosts=1]\" -n$step_pc -q $queue -- \\" >> $output
#echo "bsub -G team151 -o $outfile_CT_AS_printing_CUMULATIVE_CURVES -M $step_mem -J $name_CT_AS_printing_CUMULATIVE_CURVES -R\"select[model==Intel_Platinum]\"  -R\"select[mem>=$step_mem] rusage[mem=$step_mem] span[hosts=1]\" -n$step_pc -q $queue -- \\" >> $output
echo "\"$Rscript $Rscript_CT_AS_printing_CUMULATIVE_CURVES \\" >> $output
echo "--CUMMULATIVE_CLASSES $CUMMULATIVE_CLASSES \\" >> $output
echo "--Supp4_Table_CURATED_PLUS_PHENOTYPES $Supp4_Table_CURATED_PLUS_PHENOTYPES \\" >> $output
echo "--VAR_Prioritization_dB $VAR_Prioritization_dB \\" >> $output
echo "--ALL_dB $ALL_dB \\" >> $output
echo "--type $type --out $output_dir --out2 $output_dir_2\"" >> $output

echo "#####################################################################-----> explore_other_variants_of_design  <-----###############################################################################"  >> $output


Rscript_CT_explore_other_variants_of_design_CUMULATIVE_CURVES=$(echo "$Rscript_folder""251_MPRA_26_CUMULATIVE_CURVES_Explore_non_screened_variants_NP.R")

type=$(echo "CT_explore_other_variants_of_design_CUMULATIVE_CURVES""_""$enhancer_logval_Threshold")
outfile_CT_explore_other_variants_of_design_CUMULATIVE_CURVES=$(echo "$Log_files_path""outfile""_""$type"".out")
touch $outfile_CT_explore_other_variants_of_design_CUMULATIVE_CURVES
echo -n "" > $outfile_CT_explore_other_variants_of_design_CUMULATIVE_CURVES
name_CT_explore_other_variants_of_design_CUMULATIVE_CURVES=$(echo "$type""_""job")

CUMMULATIVE_CLASSES=$(echo "$output_dir_2""cummulative_plots/cummulative_classification_with_CTRLS.rds")
CSQ_colors=$(echo $dependencies_folder"df_CSQ_colors.rds")
ALL_dB=$(echo '/lustre/scratch123/hgi/mdt1/teams/soranzo/projects/Explore_Patricks_tables/ALL_db.tsv')
Supp4_Table_CURATED_PLUS_PHENOTYPES=$(echo "/lustre/scratch123/hgi/mdt1/teams/soranzo/projects/VARIANT_INTERPRETATION/OVERHAUL_INTERVAL/FINAL_RESULTS/Fig4_pannels/""Supp_Table_4_CURATED_Plus_phenotypes.rds")
VAR_Prioritization_dB=$(echo "/lustre/scratch126/humgen/teams/soranzo/users/mt19/V2F_paper/MPRA_E_Plus_ASE_ACTIVE_AT_LEAST_1_genIE_Tier_del_0.1/ER_Labelling_Initial_Selection_with_CSQ_labels.rds")
Selection_top_500=$(echo "/lustre/scratch126/humgen/teams/soranzo/users/mt19/NEW_TWIST_LIBS/Selection_top_500.tsv")
Selection_top_restricted_500=$(echo "/lustre/scratch126/humgen/teams/soranzo/users/mt19/NEW_TWIST_LIBS/Selection_top_restricted_500.tsv")
Selection_AS_Druggable_2=$(echo "/lustre/scratch126/humgen/teams/soranzo/users/mt19/NEW_TWIST_LIBS/Selection_AS_Druggable_2.tsv")


step_mem=$(echo  "4000")
step_pc=$(echo "1")

echo "$step_mem"
echo "$step_pc"

echo "bsub -G team151 -o $outfile_CT_explore_other_variants_of_design_CUMULATIVE_CURVES -M $step_mem -w\"done($name_CUMULATIVE_CURVES)\" -J $name_CT_explore_other_variants_of_design_CUMULATIVE_CURVES -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$step_mem] rusage[mem=$step_mem] span[hosts=1]\" -n$step_pc -q $queue -- \\" >> $output
#echo "bsub -G team151 -o $outfile_CT_explore_other_variants_of_design_CUMULATIVE_CURVES -M $step_mem -J $name_CT_explore_other_variants_of_design_CUMULATIVE_CURVES -R\"select[model==Intel_Platinum]\"  -R\"select[mem>=$step_mem] rusage[mem=$step_mem] span[hosts=1]\" -n$step_pc -q $queue -- \\" >> $output
echo "\"$Rscript $Rscript_CT_explore_other_variants_of_design_CUMULATIVE_CURVES \\" >> $output
echo "--CUMMULATIVE_CLASSES $CUMMULATIVE_CLASSES \\" >> $output
echo "--Supp4_Table_CURATED_PLUS_PHENOTYPES $Supp4_Table_CURATED_PLUS_PHENOTYPES \\" >> $output
echo "--VAR_Prioritization_dB $VAR_Prioritization_dB \\" >> $output
echo "--Selection_top_500 $Selection_top_500 \\" >> $output
echo "--Selection_top_restricted_500 $Selection_top_restricted_500 \\" >> $output
echo "--Selection_AS_Druggable_2 $Selection_AS_Druggable_2 \\" >> $output
echo "--ALL_dB $ALL_dB \\" >> $output
echo "--type $type --out $output_dir --out2 $output_dir_2\"" >> $output

echo "#####################################################################-----> disease_variants  <-----###############################################################################"  >> $output


Rscript_CT_disease_variants_CUMULATIVE_CURVES=$(echo "$Rscript_folder""251_MPRA_27_CUMMULATIVE_CURVES_disease_vs_01PP_NP.R")

type=$(echo "CT_disease_variants_CUMULATIVE_CURVES""_""$enhancer_logval_Threshold")
outfile_CT_disease_variants_CUMULATIVE_CURVES=$(echo "$Log_files_path""outfile""_""$type"".out")
touch $outfile_CT_disease_variants_CUMULATIVE_CURVES
echo -n "" > $outfile_CT_disease_variants_CUMULATIVE_CURVES
name_CT_disease_variants_CUMULATIVE_CURVES=$(echo "$type""_""job")

CUMULATIVE_CLASSES_PREPARED_leukemias=$(echo "$output_dir_2""Cumulative_frequencies_for_leukemias.rds")
CUMULATIVE_CLASSES_PREPARED_inflammatory_disease=$(echo "$output_dir_2""Cumulative_frequencies_for_inflammatory_diseases.rds")
CUMULATIVE_CLASSES_PREPARED_immunodeficiencies=$(echo "$output_dir_2""Cumulative_frequencies_for_immunodeficiencies.rds")
CUMULATIVE_CLASSES_PREPARED_blood_disease=$(echo "$output_dir_2""Cumulative_frequencies_for_blood_disease.rds")
CUMULATIVE_CLASSES_PREPARED_miscellaneous=$(echo "$output_dir_2""Cumulative_frequencies_for_miscellaneous.rds")



CSQ_colors=$(echo $dependencies_folder"df_CSQ_colors.rds")
ALL_dB=$(echo '/lustre/scratch123/hgi/mdt1/teams/soranzo/projects/Explore_Patricks_tables/ALL_db.tsv')
Supp4_Table_CURATED_PLUS_PHENOTYPES=$(echo "/lustre/scratch123/hgi/mdt1/teams/soranzo/projects/VARIANT_INTERPRETATION/OVERHAUL_INTERVAL/FINAL_RESULTS/Fig4_pannels/""Supp_Table_4_CURATED_Plus_phenotypes.rds")
VAR_Prioritization_dB=$(echo "/lustre/scratch126/humgen/teams/soranzo/users/mt19/V2F_paper/MPRA_E_Plus_ASE_ACTIVE_AT_LEAST_1_genIE_Tier_del_0.1/ER_Labelling_Initial_Selection_with_CSQ_labels.rds")
Selection_top_500=$(echo "/lustre/scratch126/humgen/teams/soranzo/users/mt19/NEW_TWIST_LIBS/Selection_top_500.tsv")
Selection_top_restricted_500=$(echo "/lustre/scratch126/humgen/teams/soranzo/users/mt19/NEW_TWIST_LIBS/Selection_top_restricted_500.tsv")
Selection_AS_Druggable_2=$(echo "/lustre/scratch126/humgen/teams/soranzo/users/mt19/NEW_TWIST_LIBS/Selection_AS_Druggable_2.tsv")


step_mem=$(echo  "4000")
step_pc=$(echo "1")

echo "$step_mem"
echo "$step_pc"

echo "bsub -G team151 -o $outfile_CT_disease_variants_CUMULATIVE_CURVES -M $step_mem -w\"done($name_CT_explore_other_variants_of_design_CUMULATIVE_CURVES)\" -J $name_CT_disease_variants_CUMULATIVE_CURVES -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$step_mem] rusage[mem=$step_mem] span[hosts=1]\" -n$step_pc -q $queue -- \\" >> $output
#echo "bsub -G team151 -o $outfile_CT_disease_variants_CUMULATIVE_CURVES -M $step_mem -J $name_CT_disease_variants_CUMULATIVE_CURVES -R\"select[model==Intel_Platinum]\"  -R\"select[mem>=$step_mem] rusage[mem=$step_mem] span[hosts=1]\" -n$step_pc -q $queue -- \\" >> $output
echo "\"$Rscript $Rscript_CT_disease_variants_CUMULATIVE_CURVES \\" >> $output
echo "--CUMULATIVE_CLASSES_PREPARED_inflammatory_disease $CUMULATIVE_CLASSES_PREPARED_inflammatory_disease \\" >> $output
echo "--CUMULATIVE_CLASSES_PREPARED_leukemias $CUMULATIVE_CLASSES_PREPARED_leukemias \\" >> $output
echo "--CUMULATIVE_CLASSES_PREPARED_immunodeficiencies $CUMULATIVE_CLASSES_PREPARED_immunodeficiencies \\" >> $output
echo "--CUMULATIVE_CLASSES_PREPARED_blood_disease $CUMULATIVE_CLASSES_PREPARED_blood_disease \\" >> $output
echo "--CUMULATIVE_CLASSES_PREPARED_miscellaneous $CUMULATIVE_CLASSES_PREPARED_miscellaneous \\" >> $output
echo "--VAR_Prioritization_dB $VAR_Prioritization_dB \\" >> $output
echo "--ALL_dB $ALL_dB \\" >> $output
echo "--type $type --out $output_dir --out2 $output_dir_2\"" >> $output

echo "#####################################################################-----> PRIORITISED_variants_CUMULATIVE_CURVES  <-----###############################################################################"  >> $output


Rscript_CT_PRIORITISED_variants_CUMULATIVE_CURVES=$(echo "$Rscript_folder""251_MPRA_28_CUMMULATIVE_CURVES_selected_bins_vs_01PP_NP.R")

type=$(echo "CT_PRIORITISED_variants_CUMULATIVE_CURVES""_""$enhancer_logval_Threshold")
outfile_CT_PRIORITISED_variants_CUMULATIVE_CURVES=$(echo "$Log_files_path""outfile""_""$type"".out")
touch $outfile_CT_PRIORITISED_variants_CUMULATIVE_CURVES
echo -n "" > $outfile_CT_PRIORITISED_variants_CUMULATIVE_CURVES
name_CT_PRIORITISED_variants_CUMULATIVE_CURVES=$(echo "$type""_""job")

CUMULATIVE_CLASSES_PREPARED_Druggable_target=$(echo "$output_dir_2""Cumulative_frequencies_for_Druggable_target.rds")
CUMULATIVE_CLASSES_PREPARED_top_500=$(echo "$output_dir_2""Cumulative_frequencies_for_top_500.rds")
CUMULATIVE_CLASSES_PREPARED_top_restricted_500=$(echo "$output_dir_2""Cumulative_frequencies_for_top_restricted_500.rds")




CSQ_colors=$(echo $dependencies_folder"df_CSQ_colors.rds")
ALL_dB=$(echo '/lustre/scratch123/hgi/mdt1/teams/soranzo/projects/Explore_Patricks_tables/ALL_db.tsv')
Supp4_Table_CURATED_PLUS_PHENOTYPES=$(echo "/lustre/scratch123/hgi/mdt1/teams/soranzo/projects/VARIANT_INTERPRETATION/OVERHAUL_INTERVAL/FINAL_RESULTS/Fig4_pannels/""Supp_Table_4_CURATED_Plus_phenotypes.rds")
VAR_Prioritization_dB=$(echo "/lustre/scratch126/humgen/teams/soranzo/users/mt19/V2F_paper/MPRA_E_Plus_ASE_ACTIVE_AT_LEAST_1_genIE_Tier_del_0.1/ER_Labelling_Initial_Selection_with_CSQ_labels.rds")
Selection_top_500=$(echo "/lustre/scratch126/humgen/teams/soranzo/users/mt19/NEW_TWIST_LIBS/Selection_top_500.tsv")
Selection_top_restricted_500=$(echo "/lustre/scratch126/humgen/teams/soranzo/users/mt19/NEW_TWIST_LIBS/Selection_top_restricted_500.tsv")
Selection_AS_Druggable_2=$(echo "/lustre/scratch126/humgen/teams/soranzo/users/mt19/NEW_TWIST_LIBS/Selection_AS_Druggable_2.tsv")


step_mem=$(echo  "4000")
step_pc=$(echo "1")

echo "$step_mem"
echo "$step_pc"

echo "bsub -G team151 -o $outfile_CT_PRIORITISED_variants_CUMULATIVE_CURVES -M $step_mem -w\"done($name_CT_explore_other_variants_of_design_CUMULATIVE_CURVES)\" -J $name_CT_PRIORITISED_variants_CUMULATIVE_CURVES -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$step_mem] rusage[mem=$step_mem] span[hosts=1]\" -n$step_pc -q $queue -- \\" >> $output
#echo "bsub -G team151 -o $outfile_CT_PRIORITISED_variants_CUMULATIVE_CURVES -M $step_mem -J $name_CT_PRIORITISED_variants_CUMULATIVE_CURVES -R\"select[model==Intel_Platinum]\"  -R\"select[mem>=$step_mem] rusage[mem=$step_mem] span[hosts=1]\" -n$step_pc -q $queue -- \\" >> $output
echo "\"$Rscript $Rscript_CT_PRIORITISED_variants_CUMULATIVE_CURVES \\" >> $output
echo "--CUMULATIVE_CLASSES_PREPARED_top_500 $CUMULATIVE_CLASSES_PREPARED_top_500 \\" >> $output
echo "--CUMULATIVE_CLASSES_PREPARED_Druggable_target $CUMULATIVE_CLASSES_PREPARED_Druggable_target \\" >> $output
echo "--CUMULATIVE_CLASSES_PREPARED_top_restricted_500 $CUMULATIVE_CLASSES_PREPARED_top_restricted_500 \\" >> $output
echo "--VAR_Prioritization_dB $VAR_Prioritization_dB \\" >> $output
echo "--ALL_dB $ALL_dB \\" >> $output
echo "--type $type --out $output_dir --out2 $output_dir_2\"" >> $output

bash $output
