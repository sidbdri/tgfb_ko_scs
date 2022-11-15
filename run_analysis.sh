#!/bin/bash

function cleanup {
   echo "Cleaning tmp..."
   rm -rf ${MAIN_DIR}/*.tmp
   echo "Killing all sub-processes..."
   kill -- -$$
}

trap exit INT
trap cleanup EXIT


export WORKON_HOME=${HOME}/.virtualenvs
export PROJECT_HOME=${HOME}/projects
source /usr/local/bin/virtualenvwrapper.sh

NUM_TOTAL_THREADS=2

MAIN_DIR=${HOME}/projects/tgfb_ko_scs
DATA_DIR=${MAIN_DIR}/data
RESULTS_DIR=${MAIN_DIR}/results
LOG_DIR=${RESULTS_DIR}/logs

mkdir -p ${RESULTS_DIR}

SPECIES=mouse
ENSEMBL_DIR=(${DATA_DIR}/mouse_ensembl_107)

# make sure we are in project ve https://github.com/sidbdri/cookiecutter-de_analysis_skeleton/issues/108
workon tgfb_ko_scs
ulimit -Sn 16000


cellranger-7.0.1 multi --id AAATHKMHV_WT --csv WT_config.csv
cellranger-7.0.1 multi --id AAATHKMHV_KO --csv KO_config.csv

cellranger-7.0.1 aggr --id AAATHKMHV_aggr --csv aggr.csv



# for sample in KO1 KO2 KO3 KO4 WT1 WT2 WT3 WT4;do
# # for sample in KO1;do
# 	mkdir -p /home/xinhe/projects/tgfb_ko_scs/igv/${sample}
# 	for cell in Microglia Endothelial Neurons pericytes Oligo Astrocytes OPC; do
# 		bash /home/xinhe/Projects/sidbdri-utils/scripts/bam4igv.sh -f \
# 		    --input-bam /home/xinhe/projects/tgfb_ko_scs/AAATHKMHV_${sample::2}/outs/per_sample_outs/${sample}/count/sample_alignments.bam \
# 		    --output-dir /home/xinhe/projects/tgfb_ko_scs/igv/${sample} \
# 		    --location-string "chr1:186622792-186705989" \
# 		    --cell-barcode-tsv /home/xinhe/projects/tgfb_ko_scs/cluster_files/${cell}.tsv \
# 		    --cores 5
# 	done
# done
#
# for type in KO WT;do
#      sambamba merge -t 10 ${type}.bam ./${type}[1234]/chr1:186622792-186705989_sample_alignments.bam
#     for cell in Microglia Endothelial Neurons pericytes Oligo Astrocytes OPC; do
#         sambamba merge -t 10 ${type}_${cell}.bam ./${type}[1234]/${cell}_chr1:186622792-186705989_sample_alignments.bam
#     done
# done
#
# for type in WT;do
#      sambamba merge -t 10 ${type}.bam ./${type}[134]/chr1:186622792-186705989_sample_alignments.bam
#     for cell in Microglia Endothelial Neurons pericytes Oligo Astrocytes OPC; do
#         sambamba merge -t 10 ${type}_${cell}.bam ./${type}[134]/${cell}_chr1:186622792-186705989_sample_alignments.bam
#     done
# done
#
# # for sample in KO1 KO2 KO3 KO4 WT1 WT2 WT3 WT4;do
# for sample in KO1 KO2 KO3 KO4 WT1 WT2 WT3 WT4;do
#     for file in ${sample}/*chr1:*;do
#         cp ${file} ./${sample}_`basename $file`
# 	done
# done
# for i in *.sql ; do
#     mv -v $i ${i%.sql}-AM.sql
# done
#
#
#
# ## appearantly this is how you pass parameter to snakemake, well...
# ## https://snakemake.readthedocs.io/en/stable/project_info/faq.html#is-it-possible-to-pass-variable-values-to-the-workflow-via-the-command-line
# export AGGR_CSV='aggr.csv'
# python3 -m snakemake -s Snakefile.cellranger aggr -j ${NUM_TOTAL_THREADS}
# ## for extra aggr runs, create the csv file and set it to AGGR_CSV and run snakemake again
# # export AGGR_CSV='aggr_run2.csv'
# # python3 -m snakemake -s Snakefile.cellranger aggr -j ${NUM_TOTAL_THREADS}
#
# # The cellranger moudle did not work for cellranger 7.0+,
# # disable for now till it updates
# # python3 -m snakemake -s Snakefile.common multiqc -j 8
#
# echo 'loupe files:'
# find . -type f -name "*cloupe*" | sort
#
# echo 'summary files:'
# find . -type f -name "web_summary.html" | sort
#
# exit;
#
# Rscript seurat_analysis_mouse.R



# to check read from cell
# samtools view sample_alignments.bam chr1:186,622,792-186,632,842 | grep