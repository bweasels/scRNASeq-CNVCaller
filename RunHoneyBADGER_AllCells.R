# Load function
source('Utils.R')

# Prep for HoneyBADGER on the scRNASeq dataset
dirs <- setDirectory()
outDir <- dirs[[1]]
dataDir <- dirs[[2]]
plottingDir <- dirs[[3]]

date <- Sys.Date()

library(HoneyBADGER)
library(Seurat)
library(GenomicRanges)
library(VariantAnnotation)
library(biomaRt)

# Load data, and the start stop locations of each transcript
data <- readRDS(paste0(outDir, 'data_all.RDS'))
geneLocs <- read.table(paste0(outDir, 'GeneLocs.txt'),
                       sep = '\t',
                       row.names = 1)

# Trim tx in data to match tx I was able to find (~500 genes different - all ribosomal varants)
data <- data[rownames(data)%in%geneLocs[,1],]
sampleAnnotation <- data@meta.data
sampleAnnotation <- data.frame(row.names = rownames(sampleAnnotation), 
                               cellType = sampleAnnotation$cell.type)

# Get the unique cell types and remove HCC so we have a list of reference cell types
ref_groups <- unique(sampleAnnotation$cellType)
ref_groups <- ref_groups[!grepl('HCC', ref_groups)]

HCC <- data[,data$cell.type=='HCC']
HCC <- as.matrix(HCC@assays$RNA@counts)

ref <- data[,data$cell.type!='HCC']
ref <- as.matrix(ref@assays$RNA@counts)

# Get the gene positions
mart.obj <- useMart(biomart = 'ENSEMBL_MART_ENSEMBL', 
                    dataset = 'hsapiens_gene_ensembl')

# make honeybadger object
hb <- new('HoneyBADGER', name='Liver_Cancer')
hb$setGexpMats(gexp.sc.init = HCC, 
               gexp.ref.init = ref, 
               mart.obj = mart.obj, 
               filter = F, 
               scale = T, 
               verbose = T)

# Model variance and necessary exp deviance to ID CNV
hb$setMvFit(verbose = T)
hb$setGexpDev(verbose=T)
hb$calcGexpCnvBoundaries(init=T,
                         verbose = T)

# Output the identified CNVs
print(range(hb$genes[unlist(hb$bound.genes.final)]))

# Retest to ID posterior probability & get final results
hb$retestIdentifiedCnvs(retestBoundGenes = T,
                        retestBoundSnps = F,
                        verbose = T)

results <- hb$summarizeResults(geneBased = T,
                               alleleBased = F)

# Make result plot
pdf(paste0(plottingDir, 'HoneyBADGERRes_AllCells_', date, '.pdf'))
trees <- hb$visualizeResults(geneBased = T,
                             alleleBased = F,
                             details = T,
                             margins = c(25, 15))

hc <- trees$hc
order <- hc$labels[hc$order]

hb$plotGexpProfile()
hb$plotGexpProfile(region = hb$cnvs[['gene-based']][['amp']])
hb$plotGexpProfile(region = hb$cnvs[['gene-based']][['del']])
