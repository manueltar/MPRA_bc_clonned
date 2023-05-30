# All of the sequencing files are stored in /nfs/team151_data02/MPRA_Nova_S2_May_2022/ . Also: the design fasta NON_PROGRAMMED_reference.fasta, the barcode reference /lustre/scratch126/humgen/teams/soranzo/users/mt19/NEW_TWIST_LIBS/NON_PROGRAMMED/combined_reference.fasta and the *deduplicated_sorted_unique.tsv.gz counting estimates

# MASTER SELECTOR

$ /software/R-3.6.1/bin/Rscript /nfs/users/nfs_m/mt19/Scripts/R/226_Prioritizing_subsets_v2.R --Rank_file /lustre/scratch119/humgen/teams/soranzo/users/ALL_dB/Allelic_Series/AS_MPRA_simplified_E_A_Allelic_Series_Generation.tsv --Rank_file_3 /lustre/scratch119/humgen/teams/soranzo/users/ALL_dB/Allelic_Series/AS_MPRA_simplified_E_A_Allelic_Series_Generation.tsv --Rank_file_4 /lustre/scratch126/humgen/teams/soranzo/users/ALL_dB/Allelic_Series/AS_Druggable_2_Allelic_Series_Generation.tsv --Rank_file_5 /lustre/scratch126/humgen/teams/soranzo/users/ALL_dB/Allelic_Series/top_500.tsv --Rank_file_6/lustre/scratch126/humgen/teams/soranzo/users/ALL_dB/Allelic_Series/top_restricted_500.tsv --Rank_file_7 /lustre/scratch119/realdata/mdt2/teams/soranzo/users/ALL_dB/Allelic_Series/AS_DEFINITIVE_Allelic_Series_Generation.tsv --Open_targets_GLOBAL_file /lustre/scratch126/humgen/teams/soranzo/users/ALL_dB/Allelic_Series_csv_tables_AS_DEFINITIVE/AS_DEFINITIVE_Open_Targets_GWAS_DISEASES_variants.tsv --dB_ALL /lustre/scratch126/humgen/teams/soranzo/users/mt19/ALL_db.tsv --PP_threshold 0.1 --cut_3 50 --Original_selection_Grouped_variants /nfs/users/nfs_m/mt19/HumGEN_seminar/Grouped_Variants.txt --categories_of_RRVV G3_RRVV_09,G2_RRVV_InfoScore --type_1 MPRA_E_A_characterization --type_2 RareVar_bins --type_3 MPRA_E_A_PROGRAMMED --type_4 AS_Druggable_2 --type_5 top_500 --type_6 top_restricted_500 --type_7 AS_DEFINITIVE --out /lustre/scratch126/humgen/teams/soranzo/users/mt19/NEW_TWIST_LIBS/ > Report_libraries.out


# Design NP library

$ bash ~/Scripts/Wraper_scripts/282_bash_script_new_library_TWIST.sh /lustre/scratch126/humgen/teams/soranzo/users/mt19/NEW_TWIST_LIBS/ Selection_NON_PROGRAMMED.tsv NON_PROGRAMMED /lustre/scratch126/humgen/teams/soranzo/users/mt19/NEW_TWIST_LIBS 4000 1 normal /lustre/scratch115/resources/ref/Homo_\
sapiens/GRCh37_53/Homo_sapiens.GRCh37.dna.all.fa ~/Scripts/Wraper_scripts/282_CREATION_NON_PROGRAMMED.sh


$ bash ~/Scripts/Wraper_scripts/283_bash_script_new_library_TWIST_partII.sh Selection_NON_PROGRAMMED.tsv NON_PROGRAMMED /lustre/scratch126/humgen/teams/soranzo/users/mt19/NEW_TWIST_LIBS 4000 1 normal /lustre/scratch115/resources/ref/Homo_sapiens/GRCh37_53/Homo_sapiens.GRCh37.dna.all.fa ~/Scripts/\
Wraper_scripts/283_CREATION_NON_PROGRAMMED.sh

$ bash ~/Scripts/Wraper_scripts/284_MIPRA_sanity_checks.sh /lustre/scratch126/humgen/teams/soranzo/users/mt19/NEW_TWIST_LIBS/NON_PROGRAMMED/ NON_PROGRAMMED_alignment_check.fasta  NON_PROGRAMMED 8000 2 normal ~/Scripts/Wraper_scripts/284_CREATION_NON_PROGRAMMED.sh

# Size of NP library

param_elements: 1914
param_REAL_TILES:       9570
param_seq_name: 21304
param_seq_name_bc:      0

# PE150 

bash ~/Scripts/Wraper_scripts/288_NON_PROGRAMMED_PE150_processer_v3_filter_option.sh /lustre/scratch126/humgen/teams/soranzo/users/mt19/NEW_TWIST_LIBS/NON_PROGRAMMED/PE150/ 40000 10 long  ~/Scripts/Wraper_scripts/288_CREATION_lane_2.sh /lustre/scratch126/humgen/teams/soranzo/users/mt19/NEW_TWIST_LIBS/NON_PROGRAMMED/PE150/Undetermined_S0_L002_R1_001.fastq.gz /lustre/scratch126/humgen/teams/soranzo/users/mt19/NEW_TWIST_LIBS/NON_PROGRAMMED/PE150/Undetermined_S0_L002_R2_001.fastq.gz /lustre/scratch126/humgen/teams/soranzo/users/mt19/NEW_TWIST_LIBS/NON_PROGRAMMED/PE150/Undetermined_S0_L002_I1_001.fastq.gz lane_2

bash ~/Scripts/Wraper_scripts/288_NON_PROGRAMMED_PE150_processer_v3_filter_option.sh /lustre/scratch126/humgen/teams/soranzo/users/mt19/NEW_TWIST_LIBS/NON_PROGRAMMED/PE150/ 32000 8 normal  ~/Scripts/Wraper_scripts/288_CREATION_lane_1.sh /lustre/scratch126/humgen/teams/soranzo/users/mt19/NEW_TWIST_LIBS/NON_PROGRAMMED/PE150/Undetermined_S0_L001_R1_001.fastq.gz /lustre/scratch126/humgen/teams/soranzo/users/mt19/NEW_TWIST_LIBS/NON_PROGRAMMED/PE150/Undetermined_S0_L001_R2_001.fastq.gz /lustre/scratch126/humgen/teams/soranzo/users/mt19/NEW_TWIST_LIBS/NON_PROGRAMMED/PE150/Undetermined_S0_L001_I1_001.fastq.gz  lane_1


# PE150_Repeat

bash ~/Scripts/Wraper_scripts/288_NON_PROGRAMMED_PE150_processer_v3_filter_option.sh /lustre/scratch126/humgen/teams/soranzo/users/mt19/NEW_TWIST_LIBS/NON_PROGRAMMED/PE150_Repeat/ 40000 10 long  ~/Scripts/Wraper_scripts/288_CREATION_lane_2.sh /lustre/scratch126/humgen/teams/soranzo/users/mt19/NEW_TWIST_LIBS/NON\
_PROGRAMMED/PE150_Repeat/Undetermined_S0_L002_R1_001.fastq.gz /lustre/scratch126/humgen/teams/soranzo/users/mt19/NEW_TWIST_LIBS/NON_PROGRAMMED/PE150_Repeat/Undetermined_S0_L002_R2_001.fastq.gz /lustre/scratch126/humgen/teams/soranzo/users/mt19/NEW_TWIST_LIBS/NON_PROGRAMMED/PE150_Repeat/Undetermined_S0_L002_I1_001.fastq.gz lane_2

bash ~/Scripts/Wraper_scripts/288_NON_PROGRAMMED_PE150_processer_v3_filter_option.sh /lustre/scratch126/humgen/teams/soranzo/users/mt19/NEW_TWIST_LIBS/NON_PROGRAMMED/PE150_Repeat/ 32000 8 normal  ~/Scripts/Wraper_scripts/288_CREATION_lane_1.sh /lustre/scratch126/humgen/teams/soranzo/users/mt19/NEW_TWIST_LIBS/NO\
N_PROGRAMMED/PE150_Repeat/Undetermined_S0_L001_R1_001.fastq.gz /lustre/scratch126/humgen/teams/soranzo/users/mt19/NEW_TWIST_LIBS/NON_PROGRAMMED/PE150_Repeat/Undetermined_S0_L001_R2_001.fastq.gz /lustre/scratch126/humgen/teams/soranzo/users/mt19/NEW_TWIST_LIBS/NON_PROGRAMMED/PE150_Repeat/Undetermined_S0_L001_I1_001.fastq.gz  lane_1

# Put together PE150 and PE150_Repeat

/software/R-3.6.1/bin/Rscript /nfs/users/nfs_m/mt19/Scripts/R/238_NON_PROGRAMMED_Nova_Seq_overall_printer.R --MPRA_Rosetta /lustre/scratch126/humgen/teams/soranzo/users/mt19/NEW_TWIST_LIBS/NON_PROGRAMMED/NON_PROGRAMMED_MPRA_Rosetta.tsv --NovaSeq_1  /lustre/scratch126/humgen/teams/soranzo/users/mt19/NEW_TWIST_LIBS/NON_PROGRAMMED/PE150/SUMMARY_TABLE_MAPQ_6_10_5.tsv --NovaSeq_2  /lustre/scratch126/humgen/teams/soranzo/users/mt19/NEW_TWIST_LIBS/NON_PROGRAMMED/PE150_Repeat/SUMMARY_TABLE_MAPQ_6_10_5.tsv  --type GLOBAL --out /lustre/scratch126/humgen/teams/soranzo/users/mt19/NEW_TWIST_LIBS/NON_PROGRAMMED/


# SE_330_HAIL_MARY

bash ~/Scripts/Wraper_scripts/288_SE_330_FWD_HAIL_MARY.sh /lustre/scratch126/humgen/teams/soranzo/users/mt19/NEW_TWIST_LIBS/NON_PROGRAMMED/SE_330_HAIL_MARY/ 32000 8 normal  ~/Scripts/Wraper_scripts/288_CREATION_lane_1.sh Pool_M_Pool_A_Pool_G_S1_L001_R1_001.fastq.gz lane_1

bash ~/Scripts/Wraper_scripts/288_SE_330_FWD_HAIL_MARY.sh /lustre/scratch126/humgen/teams/soranzo/users/mt19/NEW_TWIST_LIBS/NON_PROGRAMMED/SE_330_HAIL_MARY/ 32000 8 normal  ~/Scripts/Wraper_scripts/288_CREATION_lane_2.sh Pool_A_Pool_G_S2_L002_R1_001.fastq.gz lane_2


echo -n "" > SE_330_HAIL_MARY.out

bsub -G team151 -o SE_330_HAIL_MARY.out -M 12000 -R"select[mem >12000] rusage[mem=12000] span[hosts=1]" -n3 -q normal "/software/R-4.1.0/bin/Rscript /nfs/users/nfs_m/mt19/Scripts/R/238_NON_PROGRAMMED_Nova_Seq_HAIL_MARY_printer.R --MPRA_Rosetta /lustre/scratch126/humgen/teams/soranzo/users/mt19/NEW_TWIST_LIBS/NO\
N_PROGRAMMED/NON_PROGRAMMED_MPRA_Rosetta.tsv --QC_dnaALT_REF_threshold 0.33,3.33 --NovaSeq_1 SUMMARY_TABLE_NovaSeq_SE_no_dark.rds --NovaSeq_1_TILES TAGGING_BARCODES_TABLE_NovaSeq_SE_no_dark.rds --selected_option NoNM_MAPQ_5_5_2  --type NovaSeq_SE_no_dark --out /lustre/scratch126/humgen/teams/soranzo/users/mt19/NEW_TWIST_LIBS/NON_PROGRAMMED/SE_330_HAIL_MARY/"

# SE_330_NovaSeq_Repeat

bash ~/Scripts/Wraper_scripts/288_SE_330_no_Dark.sh /lustre/scratch126/humgen/teams/soranzo/users/mt19/NEW_TWIST_LIBS/NON_PROGRAMMED/SE_330_NovaSeq_Repeat/ 32000 8 normal  ~/Scripts/Wraper_scripts/288_CREATION_lane_1.sh Undetermined_S0_L001_R1_001.fastq.gz lane_1

bash ~/Scripts/Wraper_scripts/288_SE_330_no_Dark.sh /lustre/scratch126/humgen/teams/soranzo/users/mt19/NEW_TWIST_LIBS/NON_PROGRAMMED/SE_330_NovaSeq_Repeat/ 48000 12 normal  ~/Scripts/Wraper_scripts/288_CREATION_lane_2.sh Undetermined_S0_L002_R1_001.fastq.gz lane_2


echo -n "" > SE_330_NovaSeq_Repeat.out

bsub -G team151 -o SE_330_NovaSeq_Repeat.out -M 12000 -R"select[mem >12000] rusage[mem=12000] span[hosts=1]" -n3 -q normal "/software/R-4.1.0/bin/Rscript /nfs/users/nfs_m/mt19/Scripts/R/238_NON_PROGRAMMED_Nova_Seq_overall_printer.R --MPRA_Rosetta /lustre/scratch126/humgen/teams/soranzo/users/mt19/NEW_TWIST_LIBS\
/NON_PROGRAMMED/NON_PROGRAMMED_MPRA_Rosetta.tsv --QC_dnaALT_REF_threshold 0.33,3.33 --NovaSeq_1 SUMMARY_TABLE_NovaSeq_SE_no_dark.rds --NovaSeq_1_TILES TAGGING_BARCODES_TABLE_NovaSeq_SE_no_dark.rds --type NovaSeq_SE_no_dark --out /lustre/scratch126/humgen/teams/soranzo/users/mt19/NEW_TWIST_LIBS/NON_PROGRAMMED/SE_330_NovaSeq_Repeat/"


# build combined reference

echo -n "" > combined_reference.out


bsub -G team151 -o combined_reference.out -M 32000 -R"select[mem >32000] rusage[mem=32000] span[hosts=1]" -n8 -q normal "/software/R-4.1.0/bin/Rscript /nfs/users/nfs_m/mt19/Scripts/R/257_NON_PROGRAMMED_combiner_sequencing_HAIL_MARY.R --Fasta_Reference_1 /lustre/scratch126/humgen/teams/soranzo/users/mt19/NEW_TWIST_LIBS/NON_PROGRAMMED/reference.fasta --Fasta_Reference_2 /lustre/scratch126/humgen/teams/soranzo/users/mt19/NEW_TWIST_LIBS/NON_PROGRAMMED/SE_330_HAIL_MARY/reference.fasta --Fasta_Reference_3 /lustre/scratch126/humgen/teams/soranzo/users/mt19/NEW_TWIST_LIBS/NON_PROGRAMMED/SE_330_NovaSeq_Repeat/reference.fasta \
 --MPRA_Rosetta /lustre/scratch126/humgen/teams/soranzo/users/mt19/NEW_TWIST_LIBS/NON_PROGRAMMED/NON_PROGRAMMED_MPRA_Rosetta.tsv --type test --out /lustre/scratch126/humgen/teams/soranzo/users/mt19/NEW_TWIST_LIBS/NON_PROGRAMMED/"


# Sample processing and aligning

nohup bash /nfs/users/nfs_m/mt19/Scripts/Wraper_scripts/296_Sample_processing_NP_MPRA_Carol_Scott.sh sample_file.txt /nfs/team151_data02/MPRA_Nova_S2_May_2022/ Nova_S2_May_2022 32000 8 normal sample_file_DEMULTIPLEXED.txt /lustre/scratch126/humgen/teams/soranzo/users/mt19/NEW_TWIST_LIBS/NON_PROGRAMMED/NovaSeq_HAIL_MARY_results/ & # bwa aln -l 15 -O 100 -E 100

# Analysis of the results

bash ~/Scripts/Wraper_scripts/298_NP_MPRA_P_redo_pipeline_v3.sh /lustre/scratch126/humgen/teams/soranzo/users/mt19/NEW_TWIST_LIBS/NON_PROGRAMMED/NovaSeq_HAIL_MARY_results/ 3 16000 4 normal # 3 logpvalenhancer threshold

bash ~/Scripts/Wraper_scripts/298_NP_MPRA_P_redo_pipeline_v3.sh /lustre/scratch126/humgen/teams/soranzo/users/mt19/NEW_TWIST_LIBS/NON_PROGRAMMED/NovaSeq_HAIL_MARY_results/ 1.3 16000 4 normal # 1.3 logpval enhancer threshold

bash ~/Scripts/Wraper_scripts/315_MPRA_P_vs_NP.sh 4000 1 normal # comparison P vs NP in overlapping variants

# Rclone copy

echo -n "" > menial_3.out

bsub -G team151 -o menial_3.out -M 4000 -R"select[model==Intel_Platinum]" -R"select[mem >4000] rusage[mem=4000] span[hosts=1]" -n1 -q normal "rclone copyto NEW_LOOK_3/ mt19_g_drive:/Project_WetLab_Projects/V2F_paper/RESULTS/NP_MPRA_3/ --drive-shared-with-me"


echo -n "" > menial_1_3.out

bsub -G team151 -o menial_1_3.out -M 4000 -R"select[model==Intel_Platinum]" -R"select[mem >4000] rusage[mem=4000] span[hosts=1]" -n1 -q normal "rclone copyto NEW_LOOK_1.3/ mt19_g_drive:/Project_WetLab_Projects/V2F_paper/RESULTS/NP_MPRA_1.3/ --drive-shared-with-me"

echo -n "" > menial_P_vs_NP.out

bsub -G team151 -o menial_P_vs_NP.out -M 4000 -R"select[model==Intel_Platinum]" -R"select[mem >4000] rusage[mem=4000] span[hosts=1]" -n1 -q normal "rclone copyto P_vs_NP/ mt19_g_drive:/Project_WetLab_Projects/V2F_paper/RESULTS/P_vs_NP --drive-shared-with-me"


# To find the results associated to Open Targets diseases look in this file. The release version of the OpenTargets database used was OpenTargets/20.11 (Nov 2020)

/lustre/scratch126/humgen/teams/soranzo/users/ALL_dB/Allelic_Series_csv_tables_AS_DEFINITIVE/Open_Target_DISEASES.tsv