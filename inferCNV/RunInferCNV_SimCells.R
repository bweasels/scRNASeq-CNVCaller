library(infercnv)

# Load function
source('Utils.R')

### This file needs the sim object generated in generateSplatterDatasets.R

# Try to run InferCNV on the scRNASeq dataset
dirs <- setDirectory()
outDir <- dirs[[1]]
dataDir <- dirs[[2]]
plottingDir <- dirs[[3]]
date <- Sys.Date()

# Load data, and the start stop locations of each transcript
#data_sim <- readRDS(paste0(outDir, 'data_sim.RDS'))
data_real <- readRDS(paste0(outDir, 'data_all.RDS'))
geneLocs <- read.table(paste0(outDir, 'GeneLocs.txt'),
                       sep = '\t',
                       row.names = 1)

# Trimming dataset to genes available in GeneLocs 
data_real <- data_real[rownames(data_real)[(rownames(data_real) %in% rownames(geneLocs))]]
# Dimensions: 15252 genes x 5278 cells

# Renaming gene names to match those of fake data
rownames(geneLocs) <- paste0("Gene", seq(1,nrow(geneLocs)))

sampleAnnotation <- data_real@meta.data
# Row names of sample Annotation should match cell names of fake data 
sampleAnnotation <- data.frame(row.names = colnames(sim), 
                               cellType = sampleAnnotation$cell.type)


# Get the unique cell types and remove HCC so we have a list of reference cell types
ref_groups <- unique(sampleAnnotation$cellType)
ref_groups <- ref_groups[!grepl('HCC', ref_groups)]

# Make the infercnv object
infercnv_obj <- CreateInfercnvObject(raw_counts_matrix = as.matrix(counts(sim)),
                                     annotations_file = sampleAnnotation,
                                     gene_order_file = geneLocs,
                                     ref_group_names = ref_groups)

# Make an output directory for this run
inferCNVOut <- paste0(outDir, 'InferCNV_allFake_', date)
dir.create(inferCNVOut)

# Auto detect the number of cpus
# Otherwise it defaults to 4 and we only utilize 25% of our cloud server XD
nCores <- parallel::detectCores()

# Run infercnv
infercnv_obj <- infercnv::run(infercnv_obj,
                              cutoff = 0.1,
                              out_dir = inferCNVOut,
                              cluster_by_groups = T,
                              HMM = T,
                              num_threads = nCores)
# Reference data size after inferCNV preprocessing: Cells= 5278 Genes= 4519
