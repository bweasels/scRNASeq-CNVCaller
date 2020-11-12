## Generate our test datasets with Splatter
#Ben's Local Directory
setwd('/OneDrive/PhD/Fall 2020/Computational Genomics/scRNASeq_CNV/')
outDir <- '/OneDrive/PhD/Fall 2020/Computational Genomics/data/'
dataDir <- '/OneDrive/PhD/Fall 2020/Computational Genomics/scRNASeq-CNVCaller/Liver Cancer/Pt13.a/'
plottingDir <- '/OneDrive/PhD/Fall 2020/Computational Genomics/scRNASeq-CNVCaller/plots/'

#Load libraries and set date
library(splatter)
library(Seurat)
library(ggplot2)
library(cowplot)
date <- Sys.Date()

#set.seed locks our random number indexer so we get repeatable results
set.seed(10000)

#Seurat's function to load 10X scRNA-Seq data
#Removes genes that are expressed in 3 or fewer cells, and cells that express fewer than 200 unique genes
#data from https://www.nature.com/articles/s41467-019-14050-z?proof=t
data <- Read10X(data.dir = '/OneDrive/PhD/Fall 2020/Computational Genomics/scRNASeq-CNVCaller/Liver Cancer/Pt13.a/')
data <- CreateSeuratObject(counts = data,
                           project = 'HCC',
                           min.cells = 3,
                           min.features = 200)

# Preprocessing follows this vignette from Seurat: https://satijalab.org/seurat/v3.2/pbmc3k_tutorial.html
# Counts stored as a sparse matrix object, so they need to be converted to a normal matrix before applying regular operations 
counts <- as.matrix(data@assays$RNA@counts)

# get percent mt and plot out statistics
data$percent.mt <- colSums(counts[grep('^MT\\-.*', rownames(counts)),])/colSums(counts)

#access metadata and plot out 
metadata <- data@meta.data

vlnPlts <- VlnPlot(data, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'))

p1 <- ggplot(data = metadata, aes(x = nCount_RNA, y = percent.mt)) + geom_point()
p1 <- p1 + labs(x = "Number of Reads", y = 'Percent MT RNA')

p2 <- ggplot(data = metadata, aes(x = nCount_RNA, y = nFeature_RNA)) + geom_point()
p2 <- p2 + labs(x = 'Number of Reads', y = 'Number of Expressed Genes')

#organize the plots into a grid
preFiltered <- plot_grid(plotlist=list(p1, p2), ncol = 2)
preFiltered <- plot_grid(plotlist = list(preFiltered, vlnPlts), ncol = 1)

#Data seems homogenous and of good quality so set limits at 
# nFeature_RNA < 3000 - seurat says 2500, but that's for immune cells which are under sequenced via 3' rnaseq
# percent.mt < 0.15 - clear ~poisson continuum of percent.mt, but .15 seems to cut off the upwards skew, leaving a normal curve behind
data <- data[,data$percent.mt<0.15]
data <- data[,data$nFeature_RNA<3000]

#Make new plots and write out for posterity
metadata <- data@meta.data

vlnPlts <- VlnPlot(data, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'))

p1 <- ggplot(data = metadata, aes(x = nCount_RNA, y = percent.mt)) + geom_point()
p1 <- p1 + labs(x = "Number of Reads", y = 'Percent MT RNA')

p2 <- ggplot(data = metadata, aes(x = nCount_RNA, y = nFeature_RNA)) + geom_point()
p2 <- p2 + labs(x = 'Number of Reads', y = 'Number of Expressed Genes')

#organize the plots into a grid
postFiltered <- plot_grid(plotlist=list(p1, p2), ncol = 2)
postFiltered <- plot_grid(plotlist = list(postFiltered, vlnPlts), ncol = 1)

#Default plot width and height is 7, so plotting 3 default plots wide, and 2 default plots tall
pdf(paste0(plottingDir, '/PrePostFilteringQCPlots_', date, '.pdf'), width = 21, height = 14)
plot(preFiltered)
plot(postFiltered)
dev.off()


# Normalize the data by Reads Per 10,000 and then take natural log of the data+1
data <- NormalizeData(data, 
                      normalization.method = 'LogNormalize',
                      scale.factor = 10000)

# Get essential stats for splatter simulation
# 11K Cells took ~1.5 hours on my laptop
if(!file.exists('InitialParams_2020-11-11.RDS')){
  params <- splatEstimate(as.matrix(data@assays$RNA@counts))
  saveRDS(params, paste0('InitialParams_', date, '.RDS'))
}else{
  params <- readRDS('InitialParams_2020-11-11.RDS')
}
sim <- splatSimulate(params)
