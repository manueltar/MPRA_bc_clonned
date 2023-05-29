#0;136;0c#!/bin/bash>

 
#### Rscript
 
Rscript=/software/R-4.1.0/bin/Rscript



mem=$1
pc=$2
queue=$3 

output="/nfs/users/nfs_m/mt19/Scripts/Wraper_scripts/315_CREATION.sh"



output_dir=$(echo "/lustre/scratch126/humgen/teams/soranzo/users/mt19/NEW_TWIST_LIBS/NON_PROGRAMMED/""P_vs_NP""/")

#rm -rf $output_dir
#mkdir -p $output_dir

Log_files_path=$(echo "$output_dir""Log_files""/")


#rm -rf $Log_files_path
#mkdir -p $Log_files_path


Rscript_folder=$(echo "/nfs/users/nfs_m/mt19/Scripts/R""/")
dependencies_folder=$(echo "/lustre/scratch126/humgen/teams/soranzo/users/mt19/NEW_TWIST_LIBS/NON_PROGRAMMED/""Dependencies""/")
ALL_dB=$(echo '/lustre/scratch123/hgi/mdt1/teams/soranzo/projects/Explore_Patricks_tables/ALL_db.tsv')
dB=$(echo "$ALL_dB")
LONG_MATRIX=$(echo "/lustre/scratch126/humgen/teams/soranzo/users/mt19/NEW_TWIST_LIBS/NON_PROGRAMMED/NON_PROGRAMMED_LONG_matrix_seq_names_Plus_Flags.tsv")



touch $output
echo -n "" > $output

echo "#!/bin/bash"  >> $output
 

echo "#####################################################################-----> UPSETR_P_vs_NP <-----###############################################################################"  >> $output
echo "#####################################################################-----> Tier Printer <-----###############################################################################"  >> $output



Rscript_UPSETR_P_vs_NP=$(echo "$Rscript_folder""251_MPRA_22_P_vs_NP.R")

type=$(echo "UPSETR_P_vs_NP""_""$enhancer_logval_Threshold")
outfile_UPSETR_P_vs_NP=$(echo "$Log_files_path""outfile""_""$type"".out")
touch $outfile_UPSETR_P_vs_NP
echo -n "" > $outfile_UPSETR_P_vs_NP
name_UPSETR_P_vs_NP=$(echo "$type""_""job")

step_mem=$(echo  "4000")
step_pc=$(echo "1")

echo "$step_mem"
echo "$step_pc"


CSQ_colors=$(echo $dependencies_folder"df_CSQ_colors.rds")
KEY_collapse_P=$(echo "/lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/ASE_CHANGE/Element_collapse.rds")
KEY_collapse_NP_1_3=$(echo "/lustre/scratch126/humgen/teams/soranzo/users/mt19/NEW_TWIST_LIBS/NON_PROGRAMMED/NEW_LOOK_1.3/Element_collapse.rds")
KEY_collapse_NP_3=$(echo "/lustre/scratch126/humgen/teams/soranzo/users/mt19/NEW_TWIST_LIBS/NON_PROGRAMMED/NEW_LOOK_3/Element_collapse.rds")
df_Cell_colors=$(echo $dependencies_folder"df_Cell_colors.rds")


echo "bsub -G team151 -o $outfile_UPSETR_P_vs_NP -M $step_mem  -J $name_UPSETR_P_vs_NP -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$step_mem] rusage[mem=$step_mem] span[hosts=1]\" -n$step_pc -q $queue -- \\" >> $output
echo "\"$Rscript $Rscript_UPSETR_P_vs_NP \\" >> $output
echo "--KEY_collapse_P $KEY_collapse_P \\" >> $output
echo "--KEY_collapse_NP_1_3 $KEY_collapse_NP_1_3 \\" >> $output
echo "--KEY_collapse_NP_3 $KEY_collapse_NP_3 \\" >> $output
echo "--df_Cell_colors $df_Cell_colors \\" >> $output
echo "--CSQ_colors $CSQ_colors \\" >> $output
echo "--type $type --out $output_dir\"" >> $output


bash $output

#exit
