library(infercnv)
library(Seurat)
library(class)
library(ggplot2)
library(matrixStats)
library(digest)

source('CNVcaller/CNVcaller_Algo.R')
source('Utils.R')

# Try to run CNVcaller on the scRNASeq dataset
dirs <- setDirectory()
outDir <- dirs[[1]]
dataDir <- dirs[[2]]
plottingDir <- dirs[[3]]
date <- Sys.Date()

# Load data, and the start stop locations of each transcript
data <- readRDS(paste0(outDir, 'dataFiltered.small.RDS'))
geneLocs.temp <- read.table(paste0(dataDir, 'GeneLocs.txt'),
                            sep = '\t',
                            row.names = 1)
# So InferCNV is using all the chromosomes (including mito etc, so trim them out, and leave chr 1-22
geneLocs.temp <- geneLocs.temp[grep('chr[0-9]+$', geneLocs.temp$V2),]

# The plotting function orders chromosomes by order that they appear in our gene Locs so reorder
chrOrder <- c('chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 
              'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22')

# Yes this method of sorting is hacky, but it only took me 2 min to implement and I'm a bit rushed
geneLocs <- geneLocs.temp[1,,drop = F]
for(j in 1:length(chrOrder)){
  geneLocs <- rbind(geneLocs, geneLocs.temp[geneLocs.temp$V2%in%chrOrder[j],])
}
geneLocs <- geneLocs[2:nrow(geneLocs),]
colnames(geneLocs) <- c('Chromosome', 'Start', 'End')

# Trimming dataset to genes available in GeneLocs (~500 genes different - all ribosomal variants)
data <- data[rownames(data)[(rownames(data) %in% rownames(geneLocs))]]
sampleAnnotation <- data@meta.data
sampleAnnotation <- data.frame(row.names = rownames(sampleAnnotation), 
                               cellType = sampleAnnotation$cell.type)

# Get the unique cell types and remove HCC so we have a list of reference cell types
ref_groups <- unique(sampleAnnotation$cellType)
ref_groups <- ref_groups[!grepl('HCC', ref_groups)]

# load HallmarkGenes and stick in a list
pathways <- GSA::GSA.read.gmt(paste0(dataDir, 'h.all.v7.2.symbols.gmt'))
names <- pathways[[2]]
pathways <- pathways[[1]]
names(pathways) <- names

# Extract the pca loadings from the seurat object
pca.mat <- data@reductions$pca@cell.embeddings


# Make an output directory for this run
pThresh <- 0.1
dynRange <- 1000

minRat <- 1/sqrt(dynRange)
maxRat <- sqrt(dynRange)

infercnv_obj <- infercnv::CreateInfercnvObject(raw_counts_matrix = as.matrix(data@assays$RNA@counts), 
                                               gene_order_file = geneLocs,
                                               annotations_file = sampleAnnotation,
                                               ref_group_names = ref_groups)

CNVcallerOut <- paste0(outDir, 'CNVcaller_pathwayNormOnly_pThresh', pThresh,'_dynRange',dynRange, '_', date)
dir.create(CNVcallerOut)

# Auto detect the number of cpus
# Otherwise it defaults to 4 and we only utilize 25% of our cloud server XD
nCores <- parallel::detectCores()

# Run infercnv
#infercnv_obj <- refactored_run(infercnv_obj,
infercnv_obj <- run_no_smoothing(infercnv_obj,
                                 cutoff = 0.1,
                                 out_dir = CNVcallerOut,
                                 cluster_by_groups = T,
                                 HMM = T,
                                 num_threads = nCores,
                                 denoise = T,
                                 pathways = pathways,
                                 pcaLoadings = pca.mat,
                                 minExprRatio = minRat,
                                 maxExprRatio = maxRat,
                                 pThresh = pThresh,
                                 diagnostics = T,
                                 resume_mode = F,
                                 skipChrSmooth = F)
print(paste("Finished P Thresh:", pThresh))
