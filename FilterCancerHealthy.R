##To filter out all other weird cell types from the HCC dataset
# Use source for homebrew scripts

#Load libraries and set date
library(Seurat)
library(ggplot2)
library(cowplot)

# If working in RStudio we'll have to set our local working dir manually
# This allows our scripts to run headless
source('Utils.R')

dirs <- setDirectory()
outDir <- dirs[[1]]
dataDir <- dirs[[2]]
plottingDir <- dirs[[3]]

date <- Sys.Date()

#set.seed locks our random number indexer so we get repeatable results
set.seed(10000)

#Seurat's function to load 10X scRNA-Seq data
#Removes genes that are expressed in 3 or fewer cells, and cells that express fewer than 200 unique genes
#data from https://www.nature.com/articles/s41467-019-14050-z?proof=t
data <- Read10X(data.dir = paset0(dataDir, 'Pt13.a/')
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
#if(!file.exists('InitialParams_2020-11-15.RDS')){
#  params <- splatEstimate(as.matrix(data@assays$RNA@counts))
#  saveRDS(params, paste0('InitialParams_', date, '.RDS'))
#}else{
#  params <- readRDS('InitialParams_2020-11-15.RDS')
#}
#sim <- splatSimulate(params)

# plot out scrnaseq data and get clusters
data <- FindVariableFeatures(data,
                             selection.method = 'vst',
                             nfeatures = 2000)
genes <- rownames(data)
data <- ScaleData(data, features = genes)
data <- RunPCA(data,
               features = VariableFeatures(data))

# based on the pca, most variablility is within top 20 pcs
data <- FindNeighbors(data, dims = 1:20)

#Tried from 2 - 0.1 in intervals of 0.1, and found that 0.1 doesn't break the large clusters unhelpfully in half
data <- FindClusters(data, resolution = 0.1)

# run UMAP to cluster the data & find idents
data <- RunUMAP(data, dims = 1:20)

pdf(paste0(plottingDir, 'UMAPofClusters_', date, '.pdf'))
DimPlot(data, reduction = 'umap')
dev.off()

idents <- unique(Idents(data))
markers <- FindAllMarkers(data)
markers.list <- vector(mode = 'list', length = length(idents))
names(markers.list) <- idents
for(i in 1:length(markers.list)){
  markers.list[[i]] <- markers[markers$cluster==names(markers.list)[i],]
}

# print out top markers and compare to markers used in the paper
sapply(markers.list, function(x) print(head(x, 20)))

# Cluster 0 - HCC (Hepatocellular Carcinoma)
# Cluster 1 - myeloid derived  b cells (ignore)
# Cluster 2 - Endothelial (Regular Liver Cells)
# Cluster 3 - Cancer assoc Fibroblasts (maybe ignore - is prob v weird)
# Cluster 4 - Cancer assoc Fibroblasts (they're kinda heterogenous...)
# Cluster 5 - NK/T Cells
# Cluster 6 - T Helper cells
clusterNames <- c('HCC', 'Myeloid.Der', 'Healthy', 'CaF', 'CaF', 'NK-T.Cells', 'T-Helper')
clusterNo <- c(0, 1, 2, 3, 4, 5, 6)
data$cell.type <- clusterNames[match(Idents(data), clusterNo)]
#dir.create(paste0(outDir, 'TxToCloud/'))
#saveRDS(data, paste0(outDir, 'TxToCloud/data_all.RDS'))
saveRDS(data, paste0(outDir, 'data_all.RDS'))

# So we only want Clusters 0 & 2
data.filtered <- data[,Idents(data)%in%c(0, 2)]
data.filtered$cell.type <- ifelse(Idents(data.filtered)==0, 'HCC', 'Healthy')

# Save a full sized version to tranfer to cloud
# Has same name as other 
#saveRDS(data.filtered, paste0(outDir, 'TxToCloud/dataFiltered.RDS'))
saveRDS(data.filtered, paste0(outDir, 'dataFiltered.RDS'))

# Get HCC and Healthy Indices for proportional representation
indices.HCC <- grep('HCC', data.filtered$cell.type)
indices.Healthy <- grep('Healthy', data.filtered$cell.type)

# Sample 10% of the data to make our mini scRNASeq dataset
indices <- c(sample(indices.HCC, 
                   size = round(length(indices.HCC)*0.1), 
                   replace = F),
             sample(indices.Healthy,
                    size = round(length(indices.Healthy)*0.1),
                    replace = F))
data.small <- data.filtered[,indices]

saveRDS(data.small, paste0(outDir, 'dataFiltered.RDS'))
write.csv(data.small@meta.data, paste0(outDir, 'metadataFilt_', date, '.csv'))
