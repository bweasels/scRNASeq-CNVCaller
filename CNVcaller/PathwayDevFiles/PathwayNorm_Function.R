#First draft of the pathway caller

normalize_by_pathway <- function(infercnv_obj, # Object passed via inferCNV
                                 numNeighbors, # Number of nearest cells to bin by
                                 dist.method, # How to calc dist - options: "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski"
                                 numIter, # Number of iterations for null pathway enrichment distributions
                                 pThresh, # Padjested threshold for enriched pathways
                                 numCores, # number of cores to analyze with
                                 plottingFlag, # If we should make plots of not
                                 validateClustering){ # Limit to only two PCs for interpretability of clustering 

  source('PathwayFunctions.R')
  if(Sys.info()[[1]]=='Windows'){
    numCores <- 1 # Windows handles parallelization differently than linux
  }

  if(validateClustering){
    warning('Will only Cluster using the first 2 PCs for interpretability w/ reduced performance')
    plottingFlag <- T
  }
  
  # Extract the necessary information from the inferCNV object
  data <- infercnv_obj@expr.data[,unlist(infercnv_obj@observation_grouped_cell_indices)]
  pathways <- infercnv_obj@norm_pathways
  
  # Trim the pca data to the cells in the actual data
  data.pca <- infercnv_obj@pca_loadings[colnames(data),]
  
  # Gene_order has rownames for genes, but I want them to be in their own column so quickly do that
  GenePosition <- infercnv_obj@gene_order
  GenePosition <- cbind(Gene = rownames(GenePosition), GenePosition)
  
  # Calculate the probability of each point being the neighbor of another
  p.mat <- .CalculateDistance(PCAEmbeddings = data.pca,
                              nCores = numCores, 
                              binSize = numNeighbors,
                              plot = plottingFlag,
                              valClust = validateClustering)
  
  # Get top bin probability and then bin cells via the most probable neighbors
  nearestCells <- apply(p.mat, 2, function(x) rownames(p.mat)[order(x, decreasing = F)[1:numNeighbors]])
  binnedCells <- data
  for(i in 1:ncol(nearestCells)){
    binnedCells[,i] <- rowMeans(binnedCells[,nearestCells[,i]])
  }
  
  # Make chromosomal overlap plots to validate chr overlap
  if(plottingFlag){
    .ChromosomeCoveragePlot(pathways)
  }
  
  # Calculate the ratio of pathway expression/chromosome to overall chromosome expression
  perChrExp <- .AvgExpPerChromosome(pathways, data)
  
  # Get ranking of each gene in each cell
  binnedCells.order <- apply(binnedCells, 2, function(x) rownames(binnedCells)[order(x, decreasing = T)])
  
  # Calculate the sizes of pathways & Make a control p val population in increments of 50
  path.sizes <- sapply(pathways, function(x) sum(x%in%rownames(data)))
  nGenes <- (1:round(max(path.sizes/50)))*50
  nGenes <- as.list(nGenes)
  
  # Generate null enrichment distributions to calculate p value
  randomControls <- .GenNullEnrich(pathways = pathways,
                                   data = binnedCells,
                                   nIter = numIter,
                                   nGenes = nGenes,
                                   numCores = numCores)
  
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
  
  # get the expression matrix for normalizing
  reads <- data
  
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
  infercnv_obj@expr.data[,unlist(infercnv_obj@observation_grouped_cell_indices)] <- reads
  return(infercnv_obj)
}