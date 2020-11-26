library(infercnv)

# Load function
source('Utils.R')

# Try to run InferCNV on the scRNASeq dataset
dirs <- setDirectory()
outDir <- dirs[[1]]
dataDir <- dirs[[2]]
plottingDir <- dirs[[3]]
date <- Sys.Date()

data.paths <- list.files(path = outDir,
                         pattern = '_UMM',
                         full.names = T)
sampleNames <- gsub('.*Melanoma_(.*).RDS$', '\\1', data.paths)

# To run stuff in parallel, us mclapply (multi-core list apply) - it operates like a python lambda, but applied on elements of a list
# There are equivalents for strings (sapply), and matrices (apply), and regular lists (lapply), but it needs the code packaged into a function
mel.samps <- list(list(data.paths[1], sampleNames[1], outDir),
                  list(data.paths[2], sampleNames[2], outDir))

RunInferCNV <- function(inputs){
  # Break out the input data
  input.dir <- inputs[[1]]
  input.name <- inputs[[2]]
  outDir <- inputs[[3]]
  
  # Load data, and the start stop locations of each transcript
  data <- readRDS(input.dir)
  geneLocs <- read.table(paste0(outDir, 'GeneLocs.txt'),
                         sep = '\t',
                         row.names = 1)
  
  # Trimming dataset to genes available in GeneLocs (~500 genes different - all ribosomal variants)
  data <- data[rownames(data)[(rownames(data) %in% rownames(geneLocs))]]
  sampleAnnotation <- data@meta.data
  sampleAnnotation <- data.frame(row.names = rownames(sampleAnnotation), 
                                 cellType = sampleAnnotation$cell.type)
  
  # Get the unique cell types and remove HCC so we have a list of reference cell types
  ref_groups <- unique(sampleAnnotation$cellType)
  ref_groups <- ref_groups[!grepl('Melanoma', ref_groups)]
  
  # Make the infercnv object
  infercnv_obj <- CreateInfercnvObject(raw_counts_matrix = as.matrix(data@assays$RNA@counts),
                                       annotations_file = sampleAnnotation,
                                       gene_order_file = geneLocs,
                                       ref_group_names = ref_groups)
  
  # Make an output directory for this run
  inferCNVOut <- paste0(outDir, 'InferCNV_', input.name, '_', date)
  dir.create(inferCNVOut)
  
  # Auto detect the number of cpus
  # Since we're running in parallel, only ask for half to cores to allow for some overhead
  nCores <- parallel::detectCores()/2
  
  # Run infercnv
  infercnv_obj <- infercnv::run(infercnv_obj,
                                cutoff = 0.1,
                                out_dir = inferCNVOut,
                                cluster_by_groups = T,
                                HMM = T,
                                num_threads = nCores)
}

mclapply(mel.samps, RunInferCNV)