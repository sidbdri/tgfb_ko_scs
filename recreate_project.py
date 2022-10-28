#!/usr/bin/env python

from cookiecutter.main import cookiecutter
import os

HOME = os.path.expanduser("~")
BRANCH='master'
OVERWRITE_IF_EXISTS=True
OUT=os.path.join(HOME, "projects")

PARAMETERS={
    "_template": "/home/xinhe/Projects/cookiecutter-sc_analysis_skeleton",
    "cellranger_references": {
        "human": "refdata-gex-GRCh38-2020-A",
        "mouse": "refdata-gex-mm10-2020-A",
        "rat": "10x_refseq",
        "zebrafish": "Danio.rerio_genome"
    },
    "cellranger_version": "7.0.1",
    "ensembl_version": "107",
    "number_total_threads": "2",
    "project_name": "tgfb_ko_scs",
    "projects_base": "projects",
    "samples": "AAATHKMHV",
    "samples_dir": "/srv/data/jqiu/tgfb_ko_scs/JQ_211022/outs/fastq_path/",
    "species": "mouse",
    "virtualenv_home": ".virtualenvs"
}

COOKIECUTTER_PATH='/home/xinhe/Projects/cookiecutter-sc_analysis_skeleton'

cookiecutter(COOKIECUTTER_PATH,
             no_input=True,
             overwrite_if_exists=OVERWRITE_IF_EXISTS,
             output_dir=OUT,
             extra_context=PARAMETERS
             )

# cookiecutter('git@github.com:sidbdri/cookiecutter-de_analysis_skeleton.git',
#              no_input=True,
#              checkout=BRANCH,
#              overwrite_if_exists=OVERWRITE_IF_EXISTS,
#              output_dir=OUT,
#              extra_context=PARAMETERS
#              )