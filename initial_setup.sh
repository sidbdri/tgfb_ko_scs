#!/bin/bash

export WORKON_HOME=${HOME}/.virtualenvs
export PROJECT_HOME=${HOME}/projects
source /usr/local/bin/virtualenvwrapper.sh

## we do some init setup such as recording versions
[[ -f "./project_setup.sh" ]] && bash ./project_setup.sh

mkproject -f -p python3.8 tgfb_ko_scs

source config.sh

pip install -I --requirement requirements.txt
pip install git+https://github.com/sidbdri/transcript-utils.git@${transcript_utils_hash}
pip install --upgrade --force-reinstall git+https://github.com/HTGenomeAnalysisUnit/MultiQC.git@cellranger-module

## Clone sidbdri-utils package
rm -rf sidbdri-utils && git clone https://github.com/sidbdri/sidbdri-utils.git
(cd sidbdri-utils && git checkout -q ${sidbdri_utils_hash})
source sidbdri-utils/includes.sh

# we check sample name for special character
# https://github.com/sidbdri/cookiecutter-de_analysis_skeleton/issues/98
if [[ `echo AAATHKMHV | grep -P '[\t\n.]'` != ''  ]];then
    echo "Error: Please make sure sample names don't contain specieal characters."
    exit 1
fi

DATA_DIR=data
RNASEQ_DIR=${DATA_DIR}/rnaseq
mkdir -p ${RNASEQ_DIR}

for sample in AAATHKMHV; do
    ln -s /srv/data/jqiu/tgfb_ko_scs/JQ_211022/outs/fastq_path//$sample ${RNASEQ_DIR}/$sample
done

# we create the csv file required for cellranger aggr
if [[ ! -s 'aggr.csv' ]];then
    echo "init agge.csv for cellranger aggr..."
    echo "sample_id,molecule_h5" > aggr.csv
    for sample in AAATHKMHV; do
        echo "${sample},results/${sample}/outs/molecule_info.h5" >> aggr.csv
    done
fi



ENSEMBL_DIR=${DATA_DIR}/mouse_ensembl_107
mkdir -p ${ENSEMBL_DIR}
ln -s /srv/data/genome/mouse/refdata-gex-mm10-2020-A  ${ENSEMBL_DIR}



HUMAN_ENSEMBL_DIR=${DATA_DIR}/human_ensembl_107
mkdir -p ${HUMAN_ENSEMBL_DIR}


echo "Init git"
git init
mv gitignore .gitignore
