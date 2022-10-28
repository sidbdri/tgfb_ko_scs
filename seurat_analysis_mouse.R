# https://satijalab.org/seurat/articles/pbmc3k_tutorial.html
# https://satijalab.org/seurat/articles/essential_commands.html


SPECIES <- "mouse"
source(paste0('meta_data_', SPECIES, '.R'))

load_10x_sample <- function(sample_id, counts_dir) {
  sample_data <- Read10X(data.dir = counts_dir)
  CreateSeuratObject(counts = sample_data,
                     project = sample_id,
                     min.cells = 0, min.features = 0)
}



###############################################
#  read in 10X data and create Seurat object
#  Merge
###############################################
seurat.objects <- map2(SAMPLE_DATA$sample_id,
                       SAMPLE_DATA$sample_id %>% file.path('results',.,'outs','filtered_feature_bc_matrix'),
                       load_10x_sample)
seurat.object <- merge(seurat.objects[[1]],
                       seurat.objects[-1],
                       add.cell.ids = (SAMPLE_DATA$sample_id))

###############################################
# we can also read the aggr result in some case
###############################################
# cell_aggr <- Read10X(data.dir = 'results/aggr_samples/outs/count/filtered_feature_bc_matrix')
# seurat.object <- CreateSeuratObject(counts = cell_aggr, min.cells = 5)
# aggr_csv<-read.csv("results/aggr.csv")
# SAMPLE_DATA %<>% left_join(aggr_csv,by=c('sample_id'='sample_id'))


# we add all meta data from sample_data table
addition_meta <- seurat.object@meta.data %>%
  rownames_to_column("barcodes") %>%
  dplyr::select(barcodes,orig.ident) %>%
  left_join(SAMPLE_DATA,by=c('orig.ident'='sample_id')) %>%
  column_to_rownames('barcodes')
seurat.object<-AddMetaData(seurat.object, metadata = addition_meta)


#########################################
#                  QC                   #
#########################################
## https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4758103/
## - The number of unique genes detected in each cell.
##   - Low-quality cells or empty droplets will often have very few genes
##   - Cell doublets or multiplets may exhibit an aberrantly high gene count
## - Similarly, the total number of molecules detected within a cell (correlates strongly with unique genes)
## - The percentage of reads that map to the mitochondrial genome
##   - Low-quality / dying cells often exhibit extensive mitochondrial contamination
##   - We calculate mitochondrial QC metrics with the PercentageFeatureSet() function, which calculates the percentage of counts originating from a set of features
##   - We use the set of all genes starting with MT- as a set of mitochondrial genes
## The [[ operator can add columns to object metadata. This is a great place to stash QC stats
mitochondrial_pattern = switch(SPECIES, "human" = "^MT-", "mouse" = "^mt-")
seurat.object[["percent.mt"]] <- PercentageFeatureSet(seurat.object, pattern = mitochondrial_pattern)
# VlnPlot(seurat.object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

## FeatureScatter is typically used to visualize feature-feature relationships,
# plot1 <- FeatureScatter(seurat.object, feature1 = "nCount_RNA", feature2 = "percent.mt")
# plot2 <- FeatureScatter(seurat.object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# plot1 + plot2

## plot number of cell per sample
# seurat.object@meta.data %>% group_by(orig.ident) %>% summarise(count=n()) %>%
#   ggplot(aes(x=orig.ident, y=count)) +
#   geom_bar(stat="identity", fill="steelblue")+
#   theme_minimal() + ggtitle('#cell per sample ')

## base on the plots above, make a cutoff
seurat.object <- subset(seurat.object,
                        subset = nFeature_RNA > qc.cutoff.nFeature.min &
                          nFeature_RNA < qc.cutoff.nFeature.max &
                          percent.mt < qc.cutoff.mt.percent)

#########################################
#          Normalizing the data         #
# @todo also check sctransform
#########################################
seurat.object <- NormalizeData(seurat.object, verbose = FALSE,
                               normalization.method = "LogNormalize",
                               scale.factor = 10000)

#########################################
#       highly variable features        #
#          (feature selection)          #
#########################################
seurat.object <- FindVariableFeatures(seurat.object,
                                      selection.method = "vst",
                                      nfeatures = 1000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(seurat.object), 10)

plot1 <- VariableFeaturePlot(seurat.object) %>%
  LabelPoints(plot = ., points = top10, repel = TRUE)
plot1


#########################################
#         Scaling the data              #
#########################################
## Next, we apply a linear transformation (‘scaling’) that is a standard pre-processing step prior to dimensional reduction techniques like PCA. The ScaleData() function:
##
## Shifts the expression of each gene, so that the mean expression across cells is 0
## Scales the expression of each gene, so that the variance across cells is 1
## This step gives equal weight in downstream analyses, so that highly-expressed genes do not dominate
## The results of this are stored in pbmc[["RNA"]]@scale.data
seurat.object <- ScaleData(seurat.object, features = rownames(seurat.object))


#########################################
#               PCA       ???          #
#########################################
seurat.object <- RunPCA(seurat.object, features = VariableFeatures(object = seurat.object))
# print(seurat.object[["pca"]], dims = 1:50, nfeatures = 50)
# DimPlot(seurat.object, reduction = "pca")

##### @todo not sure needed
# seurat.object <- JackStraw(seurat.object, num.replicate = 100)
# seurat.object <- ScoreJackStraw(seurat.object, dims = 1:30)
# JackStrawPlot(seurat.object, dims = 1:20)
# ElbowPlot(seurat.object)


#########################################
#               Clustering              #
#########################################
seurat.object <- FindNeighbors(seurat.object, dims = 1:10)
## The resolution parameter in FindClusters sets the cluster granularity
## variable name changed to seurat.object so that can change the resolution
## without having to rerun above code
seurat.object <- FindClusters(seurat.object, resolution = 0.7)
seurat.object <- RunTSNE(seurat.object, dims = 1:10, check_duplicates = FALSE)
seurat.object <- RunUMAP(seurat.object, dims = 1:10)


#########################################
#               Plotting                #
#########################################
DimPlot(seurat.object, reduction = "tsne", pt.size=0.7, order = TRUE, label = TRUE)
DimPlot(seurat.object, reduction = "tsne", pt.size=0.7, order = TRUE, label = TRUE, split.by="genotype")
DimPlot(seurat.object, reduction = "umap", group.by = 'orig.ident', shape.by='genotype')


#########################################
#               DE              #
#########################################
## find all markers of cluster 2
FindMarkers(seurat.object, ident.1 = 2, min.pct = 0.25) %>% head(5)
## find all markers distinguishing cluster 5 from clusters 0 and 3
FindMarkers(seurat.object, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25) %>% head(5)
## find markers for every cluster compared to all remaining cells, report only the positive ones
# seurat.markers <- FindAllMarkers(seurat.object, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# seurat.markers %>% group_by(cluster) %>% slice_max(n = 2, order_by = avg_log2FC)

VlnPlot(seurat.object, features = c("Aif1", "Trem2"))
FeaturePlot(seurat.object, features = c("Aif1", "Trem2"))

#### cell annotation
# signatures=preprocess.signatures('meta_data/marker.csv')
# results = SCINA(seurat.object[['RNA']]@data %>% as.matrix(), signatures, max_iter = 100, convergence_n = 10,
#                 convergence_rate = 0.999, sensitivity_cutoff = 0.9, rm_overlap=TRUE, allow_unknown=TRUE, log_file='SCINA.log')
# seurat.object$scina <- results$cell_labels
# p1 <- DimPlot(seurat.object, reduction = "umap",group.by = 'scina')
# p2 <- DimPlot(seurat.object, reduction = "umap",group.by = 'orig.ident')
# p3 <- DimPlot(seurat.object, reduction = "umap",group.by = 'genotype')
# p1  + p3
# DimPlot(seurat.object, reduction = "tsne",group.by = 'scina',split.by="genotype")

# we use online cell marker database
markerfile <- 'http://xteam.xbio.top/CellMarker/download/Single_cell_markers.txt'
markers_cellmarker<- url(markerfile) %>% read_delim(delim = '\t')
my_marker<-c('Oligodendrocyte','Microglial cell','Oligodendrocyte progenitor cell',
             'Endothelial cell','Neuron','Astrocyte','Pericy')
cmmarkers <- markers_cellmarker %>%
  filter(speciesType=='Mouse',tissueType=='Brain',cellName %in% my_marker) %>%
  dplyr::select(cellName,cellMarker) %>%
  group_by(cellName) %>%
  summarise(m=paste(cellMarker,collapse = ',')) %>%
  deframe() %>% as.list() %>%
  lapply(function(x){
    strsplit(x,',')%>%unlist() %>% sub(pattern = ' ',replacement = '') %>%unique()
  })
cmresults = SCINA(seurat.object[['RNA']]@data %>% as.matrix(), cmmarkers,
                  max_iter = 100, convergence_n = 10,convergence_rate = 0.999,
                  sensitivity_cutoff = 0.9, rm_overlap=TRUE, allow_unknown=TRUE,
                  log_file='SCINA.log')
seurat.object$cm <- cmresults$cell_labels
DimPlot(seurat.object, reduction = "umap", order = TRUE, pt.size=0.7,
        group.by = 'cm', split.by="genotype")

####### DE
seurat.object.sce <- as.SingleCellExperiment(seurat.object)
assays(seurat.object.sce)
# dim(counts(seurat.object.sce))
# counts(seurat.object.sce)[1:5, 1:30]
# table(seurat.object.sce$orig.ident)
# table(seurat.object.sce$genotype)
# colData(seurat.object.sce)

# countsBySubject(scExp = seurat.object.sce, subjectVar = "orig.ident")
# subjectMetaData(scExp = seurat.object.sce, subjectVar = "orig.ident")
# summarizedCounts(scExp = seurat.object.sce, subjectVar = "orig.ident")
aggregate_counts<-aggregateBioVar(scExp = seurat.object.sce,
                                  subjectVar = "orig.ident",
                                  cellVar = "seurat_clusters")
subj_dds_dataset <- DESeqDataSetFromMatrix(
  countData = assay(aggregate_counts$`1`, "counts"),
  colData = colData(aggregate_counts$`2`),
  design = ~ genotype
)
subj_dds <- DESeq(subj_dds_dataset)
subj_dds_results <- results(subj_dds, contrast = c("genotype", "KO", "WT"))
subj_dds_transf <- as.data.frame(subj_dds_results) %>%
  bind_cols(feature = rownames(subj_dds_results)) %>%
  mutate(log_padj = - log(.data$padj, base = 10))
fs <- subj_dds_results %>% as.data.frame() %>% tibble::rownames_to_column() %>%
  arrange(padj) %>% head(20) %>% pull(rowname)
FeaturePlot(seurat.object, features = fs[c(1:2)])


## Save the objects in the workspace for future analysis
rws <- "results/Rworkspace/"
if (!dir.exists(rws)) {
  dir.create(rws, recursive = TRUE)
}

save(list = ls(), file = paste0(rws, "sc_", SPECIES,".RData"))



# #########################################
# # integration   #
# #########################################
# # split the dataset into a list of seurat objects (stim and CTRL)
# seurat.object.list <- SplitObject(seurat.object, split.by = 'batch')
# # normalize and identify variable features for each dataset independently
# seurat.object.list <- lapply(X = seurat.object.list, FUN = function(x) {
#   x <- NormalizeData(x)
#   x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
# })
# features <- SelectIntegrationFeatures(object.list = seurat.object.list)
# batch.anchors <- FindIntegrationAnchors(object.list = seurat.object.list, anchor.features = features)
# batch.combined <- IntegrateData(anchorset = batch.anchors)
# DefaultAssay(batch.combined) <- "integrated"
# batch.combined <- ScaleData(batch.combined, verbose = FALSE)
# batch.combined <- RunPCA(batch.combined, npcs = 30, verbose = FALSE)
# batch.combined <- RunUMAP(batch.combined, reduction = "pca", dims = 1:30)
# batch.combined <- FindNeighbors(batch.combined, reduction = "pca", dims = 1:30)
# batch.combined <- FindClusters(batch.combined, resolution = 0.5)
# p1 <- DimPlot(batch.combined, reduction = "umap", group.by = "batch")
# p2 <- DimPlot(batch.combined, reduction = "umap", label = TRUE, repel = TRUE)
# p1 + p2
# DimPlot(batch.combined, reduction = "umap", split.by = "batch")