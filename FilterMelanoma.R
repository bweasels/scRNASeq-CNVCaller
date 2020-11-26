#Analyze the melanoma samples
source('Utils.R')
dirs <- setDirectory()
outDir <- dirs[[1]]
dataDir <- dirs[[2]]
plottingDir <- dirs[[3]]

library(Seurat)
library(ggplot2)
library(cowplot)

date <- Sys.Date()

# Get the directories with the melanoma samples
mel.dirs <- list.dirs(path = dataDir,
                      full.names = T)
mel.dirs <- grep('UMM', mel.dirs, value = T)
mel.samps <- gsub('.*(UMM.*)', '\\1', mel.dirs)

# These labels are based on the plots generated between lines 118 & 136
clusterLabels <- list(data.frame(Cluster = c(0, 1, 2, 3, 4, 5, 6), 
                                 Label = c('Immune', 'Immune', 'Melanoma', 'Immune', 'Healthy', 'Melanoma', 'Immune')),
                      data.frame(Cluster = c(0, 1, 2, 3, 4),
                                 Label = c('Melanoma', 'Melanoma', 'Immune', 'Immune', 'Healthy')))
names(clusterLabels) <- mel.samps

# Load the melanoma dataset and follow the paper's analysis methods
for(i in 1:2){
  data <- Read10X(data.dir = mel.dirs[i])
  data <- CreateSeuratObject(counts = data,
                             project = mel.samps[i],
                             min.cells = 3,
                             min.features = 400)
  
  # get percent mt and plot out statistics
  counts <- as.matrix(data@assays$RNA@counts)
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
  data <- data[,data$percent.mt<0.1]
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
  pdf(paste0(plottingDir, '/InitQCPlots_', mel.samps[i], '_', date, '.pdf'), width = 21, height = 14)
  plot(preFiltered)
  plot(postFiltered)
  dev.off()
  
  # Normalize data and scale for PCA
  # Normalize the data by Reads Per 10,000 and then take natural log of the data+1
  data <- NormalizeData(data, 
                        normalization.method = 'LogNormalize',
                        scale.factor = 10000)
  
  # plot out scrnaseq data and get clusters
  data <- FindVariableFeatures(data,
                               selection.method = 'vst',
                               nfeatures = 2000, )
  genes <- rownames(data)
  data <- ScaleData(data, features = genes)
  data <- RunPCA(data,
                 features = VariableFeatures(data))
  
  # Cluster using original Louvian algo (algorithm = 1)
  data <- FindNeighbors(data, dims = 1:20)
  
  # I'm sorry, but their resolution is kinda bupkus
  data <- FindClusters(data, resolution = 0.1, algorithm = 1)
  
  # Plot the tsne out - the paper uses a t-SNE
  data <- RunTSNE(data, dims = 1:20)
  data <- RunUMAP(data, dims = 1:20)
  
  # Get the markers and cluster names
  idents <- unique(Idents(data))
  markers <- FindAllMarkers(data)
  markers.list <- vector(mode = 'list', length = length(idents))
  names(markers.list) <- idents
  for(j in 1:length(markers.list)){
    markers.list[[j]] <- markers[markers$cluster==names(markers.list)[j],]
  }
  
  # print out top markers and compare to markers used in the paper
  sapply(markers.list, function(x) print(head(x, 20)))
  
  ## Marker Genes Lists
  Mel <- c('MLANA', 'MITF', 'DCT', 'PRAME', 'GEP')
  Immune <- c('CD3D', 'CD3E', 'CD8A', 'CD19', 'CD79A', 'MS4A1', 
              'CD20', 'IGHG1', 'MZB1', 'SDC1', 'CD79A', 'CD68', 
              'CD163', 'CD14', 'FGFBP2', 'FCG3RA', 'CX3CR1')
  Healthy <- c('RPE65', 'RCVRN', 'FGF7', 'PECAM1', 'VWF')
  
  # Make gene list scores
  reads <- as.matrix(data@assays$RNA@data)
  data$mel.score <- colSums(reads[rownames(reads)%in%Mel,])
  data$Imm.score <- colSums(reads[rownames(reads)%in%Immune,])
  data$Healthy.score <- colSums(reads[rownames(reads)%in%Healthy,])
  
  # So the paper does t-sne, but UMAP seperates the clusters singificantly better
  P.mel <- FeaturePlot(data, features = 'mel.score', reduction = 'umap')
  P.imm <- FeaturePlot(data, features = 'Imm.score', reduction = 'umap')
  P.healthy <- FeaturePlot(data, features = 'Healthy.score', reduction = 'umap')
  P.clusters <- DimPlot(data, reduction = 'umap')
  
  pdf(paste0(plottingDir, '/UMAPs_', mel.samps[i], '_', date, '.pdf'), width = 16, height = 14)
  plot_grid(plotlist = list(P.clusters, P.mel, P.imm, P.healthy), ncol = 2)
  dev.off()
  
  data$cell.type <- clusterLabels[[i]]$Label[match(Idents(data), clusterLabels[[i]]$Cluster)]
  saveRDS(data, paste0(outDir, 'Melanoma_', mel.samps[i], '.RDS'))
}