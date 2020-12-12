#' @include inferCNV.R
NULL

#' @title .calculateDistance()
#'
#' @description Calculate the distance between points in the pca reduction 
#' Also calculates Joccard's index to calculate probability of neighbors being neighbors
#' Joccards index helps reduce large distance neigbors
#' 
#' @param PCAEmbeddings PCA representation of the data
#' 
#' @param nCores Number of cores to use
#' 
#' @param binSize number of neighbors to cluster by
#' 
#' @param plot boolean to plot or not
#' 
#' @param valClust boolean to limit to first 2 pcs to validate clustering algos
#' 
#' @return p.mat
#' 
#' @keywords internal
#' @noRd
#' 

.CalculateDistance <- function(PCAEmbeddings, nCores, binSize, dist.method, plot=F, valClust = F){
  require(parallel)
  
  # Identify the distance from one point to another
  if(valClust){
    dists <- as.matrix(dist(PCAEmbeddings[,1:2], 
                            method = dist.method))
  } else{
    dists <- as.matrix(dist(PCAEmbeddings, 
                            method = dist.method))
  }
  # Convert unique names to numbers to keep memory low
  lookup <- data.frame(CellBarcode = rownames(dists),
                       Code = 1:nrow(dists))
  rownames(dists) <- lookup$Code
  
  dists <- as.list(as.data.frame(dists))
  
  # Sort each column & replace with the code for the closest to farthest cells
  dists <- mclapply(dists, function(x) as.numeric(lookup$Code[order(x)]), 
                       mc.cores = nCores)
  
  # Take 2N first frows of the dist mat and calculate prob belonging
  dists <- mclapply(dists, function(x) x[1:(1+binSize*2)], mc.cores = nCores)
  
  if(plot){
    plot.df <- as.data.frame(dists)
    reps <- paste0('rep', 1:10)
    # sample 10 sets of random cells and color the adjacent cells the same color
    plot.df <- matrix(ncol = 10, nrow = length(dists), dimnames = list(names(dists), paste0('rep', 1:10)))
    for(i in 1:ncol(plot.df)){
      cells <- sample(names(dists), length(dists)*0.1)
      for(j in 1:length(cells)){
        pairedCells <- dists[[grep(cells[j], names(dists))]][1:5]
        plot.df[pairedCells,i] <- lookup$Code[grep(cells[j], lookup$CellBarcode)]
      }
    }
    plot.df <- data.frame(plot.df, PCA.x = PCAEmbeddings[,1], PCA.y = PCAEmbeddings[,2])
    pdf(paste0(plottingDir, 'ClusteredCellsDist.pdf'))
    for(i in 1:length(reps)){
      p <- ggplot(plot.df, aes_string(x='PCA.x', y='PCA.y', color = reps[i])) + geom_point(size = 0.5)
      p <- p + labs(x = 'PC.1', y='PC.2', title = paste('Random 10% of Clusters |', reps[i]))
      plot(p)
    }
    dev.off()
  }
  
  # Calculate probability matrix (AUB)/(AnB) (Joccard Index)
  # This is an extremely slow step, so validate if it is necessary
  p.mat <- mclapply(dists, function(x) lapply(dists, function(y) sum(x%in%y)/length(unique(c(x,y)))),
                    mc.cores = nCores)
  p.mat <- mclapply(p.mat, unlist, mc.cores = nCores)
  p.mat <- do.call(rbind, p.mat)
  
  # Get top bin probability and then bin cells via the most probable neighbors
  nearestCells <- apply(p.mat, 2, function(x) rownames(p.mat)[order(x, decreasing = T)[1:binSize]])
  if(plot){
    plot.df <- as.data.frame(nearestCells)
    reps <- paste0('rep', 1:10)
    # sample 10 sets of random cells and color the adjacent cells the same color
    plot.df <- matrix(ncol = 10, nrow = ncol(nearestCells), dimnames = list(colnames(nearestCells), paste0('rep', 1:10)))
    for(i in 1:ncol(plot.df)){
      cells <- sample(colnames(nearestCells), ncol(nearestCells)*0.1)
      for(j in 1:length(cells)){
        pairedCells <- nearestCells[,grep(cells[j], colnames(nearestCells))]
        plot.df[pairedCells,i] <- lookup$Code[grep(cells[j], lookup$CellBarcode)]
      }
    }
    plot.df <- data.frame(plot.df, PCA.x = PCAEmbeddings[,1], PCA.y = PCAEmbeddings[,2])
    pdf(paste0(plottingDir, 'ClusteredCellsJoccardDist.pdf'))
    for(i in 1:length(reps)){
      p <- ggplot(plot.df, aes_string(x='PCA.x', y='PCA.y', color = reps[i])) + geom_point(size = 0.5)
      p <- p + labs(x = 'PC.1', y='PC.2', title = paste('Random 10% of Clusters |', reps[i]))
      plot(p)
    }
    dev.off()
  }
  # This is a for loop implementation of the above - mcapply should be faster than the nested for loops
  #for(i in 1:length(dists)){
  #  col <- dists[[i]]
  #  for(j in 1:length(dists)){
  #    intersection <- unique(c(col, dists[[j]]))
  #    p.mat[i,j] <- sum(col%in%dists[[j]])/length(intersection)
  #  }
  #}

  return(p.mat)
}


#' @title .ChromosomalCoveragePlot()
#'
#' @description Function to plot the distribution of gene sites across the chromosomes
#' 
#' @param pathways list of genes included in pathways
#' 
#' @param GenePosition list of gene start and stop sites and chr coverage
#' 
#' @return None
#' 
#' @keywords internal
#' @noRd
#' 
.ChromosomeCoveragePlot <- function(pathways,
                                    GenePosition){
  
  pathwayLocs <- c()
  for(i in 1:length(pathways)){
    # For each pathway, get the start/stop sites
    pathway <- pathways[[i]]
    name <- names(pathways)[i]
    for(j in 1:length(pathway)){
      
      # Get the genes in this pathway
      info <- unlist(GenePosition[GenePosition$Gene==pathway[j],])
      if(length(info)!=0){
        start <- info[3]
        end <- info[4]
        
        # If an exon is more than 1kbp, repeat an entry every 1kbp to create steaks of points
        pos <- seq(from = start, to = end, by = 1000)
        chunk <- data.frame(pathway = rep(name, length(pos)),
                            loci = pos,
                            chromosome = rep(info[2], length(pos)))
        pathwayLocs <- rbind(pathwayLocs, chunk)
      }
    }
    print(paste0('Finished Pathway: ', names(pathways)[i]))
  }
  
  # Make an overlap plot for each chromosome
  chromosomes <- paste0('chr', 1:22)
  
  options(scipen = 8)
  pdf(paste0(plottingDir, 'ChromosomalCoveragePerPathway.pdf'), width = 21, height = 10)
  for(chr in chromosomes){
    p <- ggplot(pathwayLocs[pathwayLocs$chromosome==chr,], aes(x = loci, y = pathway)) + geom_point(size = 0.8)
    p <- p + labs(x = "Position", title = chr)
    plot(p)
  }
  dev.off()
  
}

#' @title .AvgExpPerChromosome()
#'
#' @description calculate the average expression per chromsome for the 
#' cells in question. Calculates the ratio of the pathway expression to 
#' overall average expression per chromosome
#' 
#' @param pathways List of genes in each pathway
#' 
#' @param data expression matrix to use
#' 
#' @param GenePositions matrix of the genes and the start stop locations
#' 
#' @return perChrExp
#' 
#' @keywords internal
#' @noRd
#' 
.AvgExpPerChromosome <- function(pathways, data, GenePosition){
  
  # Create list of genes in each chromosome
  perChrExp <- array(0, dim = c(length(pathways),
                                 ncol(data),
                                 length(unique(GenePosition$Chromosome))),
                     dimnames = list(names(pathways),
                                     colnames(data),
                                     unique(GenePosition$Chromosome)))
  
  # populate each element in the list with the ratio of the pathway expression within chromosome to chromosomal average
  for(i in 1:dim(perChrExp)[3]){
    
    # Get the genes in each chromosome
    chrName <- dimnames(perChrExp)[[3]][i]
    genes <- GenePosition$Gene[GenePosition$Chromosome==chrName]
    data.subset <- data[genes,]
    
    # Get the average expression per cell per chromosome
    avgExpCell <- colMeans(data.subset)
    
    # For each pathway, make ratio | mean(pathwayExp/chr):mean(avgExp/chr)
    for(j in 1:length(pathways)){
      pway.genes <- pathways[[j]]
      
      # If there are any genes in that chromosome calc ratio and store in 3d array
      if(any(pway.genes%in%rownames(data.subset))){
        pway.expression <- colMeans(data.subset[rownames(data.subset)%in%pway.genes,,drop = F])
        ratio <- pway.expression/avgExpCell
        perChrExp[j,,i] <- ratio
      }
    }
    print(paste('Finished chromosome: ', dimnames(perChrExp)[[3]][i]))
  }
  return(perChrExp)
}

#' @title .GenNullEnrich
#'
#' @description Generates null pathway enrichment distributions for p value
#' 
#' @param pathways list of pathways to generate null dists for
#' 
#' @param data Expression Matrix
#' 
#' @param nIter number of iterations to perform this over
#' 
#' @param nGenes list of number of genes needed to generate null distributions for
#' 
#' @param numCores number of cpus to run on
#' 
#' @return randomControls
#' 
#' @keywords internal
#' @noRd
#' 

.GenNullEnrich <- function(pathways, data, nIter, nGenes, numCores){
  
  # Calculate non-zero median of each gene
  medianExpr <- apply(data, 1, median)
  names(medianExpr) <- rownames(data)
  medianExpr <- sort(medianExpr, decreasing = T)
  
  # Define a function to do the simulation so we can use mclapply
  .ctrlSim <- function(nGenes, nIter, medianExpr){
    out <- rep(0, times = nIter)
    for(i in 1:nIter){
      rand.genes <- sample(x = names(medianExpr), size = nGenes)
      out[i] <- median(match(rand.genes, names(medianExpr)))
    }
    return(out)
  }
  
  # run the simulations in parallel bcos 100,000K sims are a bit slow
  randomControls <- mclapply(nGenes, function(x) .ctrlSim(x, nIter, medianExpr), mc.cores = numCores)
  randomControls <- do.call(cbind, randomControls)
  colnames(randomControls) <- paste0('NumGenes', nGenes)
  return(randomControls)
}