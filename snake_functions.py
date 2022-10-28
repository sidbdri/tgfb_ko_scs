from snake_variables import *
from glob import glob
import os

def retrieve_fastqs(sample):
    fastqs = glob(os.path.join(DATA_DIR,'rnaseq',sample))
    return fastqs

