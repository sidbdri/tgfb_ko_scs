#!/bin/bash

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