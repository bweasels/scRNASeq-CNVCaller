#First draft of the pathway caller

source('Utils.R')
dirs <- setDirectory()
outDir <- dirs[[1]]
dataDir <- dirs[[2]]
plottingDir <- dirs[[3]]

library(Seurat)
library(class)
library(pheatmap)
library(ggplot2)
library(matrixStats)
source('PathwayFunctions.R')

###Key Variables
binSize <- 5 #Number of cells to bin by
dist.method <- 'euclidean' # Options: "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski"
nIter <- 100000 #Number of iterations for the null pvalue distributions
pThresh <- 0.05 #Pvalue threshold for enriched pathways
numCores <- detectCores() # Num cores for my parallel hackery to solve my slow code
plotChrOverlap <- F # If trying new pathways, validate chromosomal coverage with a plot

# Load data and simulate InferCNV's preprocessing steps
data <- readRDS('data/TxToCloud/data_all.RDS')
data <- NormalizeData(data, scale.factor=10000, verbose = T, normalization.method = 'RC')
geneAvg <- Matrix::rowMeans(data@assays$RNA@counts)
print(paste('Number of Genes filtered out:', sum(geneAvg<0.1)))
data <- data[geneAvg>0.1,]

# load HallmarkGenes and stick in a list
pathways <- GSA::GSA.read.gmt('h.all.v7.2.symbols.gmt')
names <- pathways[[2]]
pathways <- pathways[[1]]
names(pathways) <- names

# Load gene positions and trim to the genes that are in both position and data
GenePosition <- read.table('data/GeneLocs.txt', sep = '\t')
colnames(GenePosition) <- c('Gene', 'Chromosome', 'Start', 'End')
GenePosition <- GenePosition[grep('chr[0-9]+', GenePosition$Chromosome),]
GenePosition <- GenePosition[GenePosition$Gene%in%rownames(data), ]
data <- data[GenePosition$Gene,]

# Calculate PCA on the trimmed data and extract embeddings
data <- RunPCA(data)
data.pca <- Embeddings(data, reduction='pca')

# Calculate the probability of each point being the neighbor of another
p.mat <- .CalculateDistance(PCAEmbeddings = data.pca,
                           nCores = numCores, 
                           binSize = binSize)

# Get top bin probability and then bin cells via the most probable neighbors
nearestCells <- apply(p.mat, 2, function(x) rownames(p.mat)[order(x, decreasing = T)[1:binSize]])
binnedCells <- as.matrix(data@assays$RNA@data)
for(i in 1:ncol(nearestCells)){
  binnedCells[,i] <- rowMeans(binnedCells[,nearestCells[,i]])
}

# Make chromosomal overlap plots to validate chr overlap
if(plotChrOverlap){
  .ChromosomeCoveragePlot(pathways)
}

# Calculate the ratio of pathway expression/chromosome to overall chromosome expression
perChrExp <- .AvgExpPerChromosome(pathways, as.matrix(data@assays$RNA@data))

# Get ranking of each gene in each cell
binnedCells.order <- apply(binnedCells, 2, function(x) rownames(binnedCells)[order(x, decreasing = T)])

# Calculate the sizes of pathways & Make a control p val population in increments of 50
path.sizes <- sapply(pathways, function(x) sum(x%in%rownames(data)))
nGenes <- (1:round(max(path.sizes/50)))*50
nGenes <- as.list(nGenes)

# Generate null enrichment distributions to calculate p value
randomControls <- .GenNullEnrich(pathways = pathways,
                                 data = binnedCells,
                                 nIter = nIter,
                                 nGenes = nGenes)

# Calculate Mean and SD for each of the nGenes null Distribution
rand.stats <- data.frame(apply(randomControls, 2, function(x) c(mean(x), sd(x))))

# Calculate the median binned rank for each pathway
medPathwayRank <- matrix(ncol = ncol(binnedCells),
                         nrow = length(pathways))
rownames(medPathwayRank) <- names(pathways)

for(i in 1:nrow(medPathwayRank)){
  genes <- pathways[[i]]
  for(j in 1:ncol(medPathwayRank)){
    medPathwayRank[i,j] <- median(match(genes, binnedCells.order[,j]), na.rm = T)
  }
}

# Get the index of the best random control for the pathway (is the pathway closest to the 50, 100, 150 ... gene random control?)
RandIndex <- sapply(pathways, function(x) which.min(abs(sum(x%in%rownames(data))-unlist(nGenes))))

# Calculate the pValue
for(i in 1:nrow(medPathwayRank)){
  ExprRanks <- medPathwayRank[i,]
  mean <- rand.stats[1, RandIndex[i]]
  sd <- rand.stats[2, RandIndex[i]]
  medPathwayRank[i,] <- pnorm((ExprRanks-mean)/sd)
}

# Correct for multiple hypothesis testing for each pathway (multiply by the number of cells)
medPathwayRank <- medPathwayRank*ncol(medPathwayRank)

# get the expression matri for normalizing
reads <- as.matrix(data@assays$RNA@data)

# normalize each chr expression for significant pathways
for(i in 1:ncol(medPathwayRank)){
  sigPathways <- rownames(medPathwayRank)[medPathwayRank[,i]<0.05]
  for(path in sigPathways){
    path.genes <- pathways[[path]]
    path.genes <- GenePosition[GenePosition$Gene%in%path.genes,]
    for(chr in unique(path.genes$Chromosome)){
      chr.path.genes <- path.genes[path.genes$Chromosome==chr,]
      reads.chr.path <- reads[rownames(reads)%in%chr.path.genes$Gene,i,drop = F]
      
      #Get the pathway expression ratio from the array and divide the binned counts by that expression
      reads.chr.path <- reads.chr.path/perChrExp[path,i,chr]
      
      #Save the altered reads
      reads[rownames(reads)%in%chr.path.genes$Gene, i] <- reads.chr.path
    }
  }
}
