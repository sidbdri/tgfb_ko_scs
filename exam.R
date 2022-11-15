## we want to be able to find reads from a regions, and find out what cell expressing these reads.
library(GenomicAlignments)
library(dplyr)

bam_file <- 'igv/KO.bam'
roi <- 'chr1:186623186-186629036'
strand <- '-'
picard_executable <- '/opt/picard-tools-2.25.6/picard.jar'
cellranger_barcodes_tsv_gz <- 'AAATHKMHV_aggr/outs/count/filtered_feature_bc_matrix/barcodes.tsv.gz'
cell_barcode_postfix <- '-5'
loupe_import_file <- 'filtered_by_read_id_loupe_import.tsv'

gr = readGAlignments(bam_file,use.names=T)
# gr['VL00273:7:AAATHKMHV:1:2604:52486:1473']

roi_vector <- roi %>% strsplit('[:-]') %>% unlist()
q=GRanges(seqnames=roi_vector[1],ranges=IRanges(start = roi_vector[2], end = roi_vector[3]),strand=strand)
## this contains all reads overlap with roi
g <-subsetByOverlaps(gr, q)

## we sometimes want to check KO exon reads, thus requires futher filtering
# g <-g[start(g)>=186622787 & start(g)<=186706220 & end(g)>186628883 & end(g)<186629036  & strand(g)=='-' & njunc(g)>0]

# we produce a list of reads of interests
g %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column('read_name') %>% 
  select(read_name) %>% 
  write_delim('read_ids.txt',delim = '\n',col_names = F)

## we can create a bam file containing only these reads, for igv inspection
# cmd <- paste0('java -jar ',picard_executable,' FilterSamReads -I ',bam_file,' -O filtered_by_read_id.bam -READ_LIST_FILE read_ids.txt -FILTER includeReadList')
# system(cmd)
# system('samtools index filtered_by_read_id.bam')

## we find cells that express these reads, the CB field in the bam file contains the cell barcode.
## https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/bam
library(Rsamtools)
ScanBamParam(what=c("rname", "strand", "pos", "qwidth","CB"))
bam <- scanBam("filtered_by_read_id.bam",param = ScanBamParam(what=c("qname"),tag=c("CB")))
read2cell <- bam[[1]]$tag$CB 
names(read2cell)<-bam[[1]]$qname

## cellranger will filter out  some cell, we want to only include cells that exist in the cellranger output
existing_barcode<-readr::read_csv(cellranger_barcodes_tsv_gz,col_names = 'Barcode')%>%pull(Barcode)
## in a cellranger count output, the cell barcode ends with -1
## in a cellranger multi output, the cell barcode ends with the sample index, i.e -1, -2, -3...
## depend on which output using, this can changed. 
tb.out <- read2cell %>% unique() %>% 
  sub(pattern = '-1',replacement = cell_barcode_postfix) %>% 
  as.data.frame() %>% set_names('Barcode') %>%
  # we set a dummy so loupe browser can import it.
  dplyr::mutate(cell_of_interests=1) %>%
  mutate(valid=Barcode %in% existing_barcode)

## these are the cells expressing the reads, but filtered by cell ranger
tb.out %>% filter(!valid)

## output the cell ids
tb.out %>% filter(valid) %>% 
  write_csv(loupe_import_file)

message('Please see ',loupe_import_file, ' for a list of cell barcodes which contains the reads from the region of interests.')
message('This file can be directly imported to loupe browser')