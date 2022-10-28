import os

HOME_DIR = os.path.expanduser("~")
MAIN_DIR = os.path.join(HOME_DIR, 'projects/tgfb_ko_scs')
DATA_DIR = os.path.join(MAIN_DIR, 'data')
SAMPLES="AAATHKMHV"
SAMPLES=SAMPLES.split()
CELLRANGER_EXECUTABLE="/usr/local/bin/cellranger-7.0.1"
MAPPER_INDEX = ("%s/mouse_ensembl_107/refdata-gex-mm10-2020-A" % DATA_DIR)
SPECIES="mouse"
ENSEMBL_DIR=("%s/mouse_ensembl_107" % DATA_DIR)
AGGR_CSV = os.environ.get("AGGR_CSV", "aggr.csv")
AGGR_ID = os.path.splitext(os.path.basename(AGGR_CSV))[0]