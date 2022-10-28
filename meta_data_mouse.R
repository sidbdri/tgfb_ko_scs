library(dplyr)
library(Seurat)
library(cowplot)
library(purrr)
library(tibble)
library(magrittr)
library(ggplot2)
library(readr)
library(SCINA)

library(aggregateBioVar)
library(SummarizedExperiment)
library(SingleCellExperiment)
library(DESeq2)

# For data transformation and visualization
library(magrittr)
library(cowplot)
library(ggtext)

qc.cutoff.mt.percent=5
qc.cutoff.nFeature.min = 10
qc.cutoff.nFeature.max = 99999

SPECIES='mouse'

SAMPLE_DATA <- tribble(
  ~sample_name, ~genotype, ~counts_dir,
  # "N1","KO","N1-SCI7T002-SCI5T002_HMTVVDSX3_count/outs/filtered_feature_bc_matrix",
)