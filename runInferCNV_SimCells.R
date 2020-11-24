library(infercnv)

# Load function
source('Utils.R')

# Try to run InferCNV on the scRNASeq dataset
dirs <- setDirectory()
outDir <- dirs[[1]]
dataDir <- dirs[[2]]
plottingDir <- dirs[[3]]
date <- Sys.Date()

# Load data, and the start stop locations of each transcript
data_sim <- readRDS(paste0(outDir, 'data_sim.RDS'))
data_real <- readRDS(paste0(outDir, 'data_all.RDS'))
geneLocs <- read.table(paste0(outDir, 'GeneLocs.txt'),
                       sep = '\t',
                       row.names = 1)

# Trim tx in data to match tx I was able to find (~500 genes different - all ribosomal variants)
data_sim <- data_sim[rownames(data_sim)%in%geneLocs[,1],]
data_real <- data_real[rownames(data_real)%in%geneLocs[,1],]
sampleAnnotation <- data_real@meta.data
sampleAnnotation <- data.frame(row.names = rownames(sampleAnnotation), 
                               cellType = sampleAnnotation$cell.type)

# Get the unique cell types and remove HCC so we have a list of reference cell types
ref_groups <- unique(sampleAnnotation$cellType)
ref_groups <- ref_groups[!grepl('HCC', ref_groups)]

# Make the infercnv object
#TODO: figure out expression in as.matrix() and adapt to sim data
infercnv_obj <- CreateInfercnvObject(raw_counts_matrix = as.matrix(data@assays$RNA@counts),
                                     annotations_file = sampleAnnotation,
                                     gene_order_file = geneLocs,
                                     ref_group_names = ref_groups)

# Make an output directory for this run
inferCNVOut <- paste0(outDir, 'InferCNV_allCells_', date)
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
