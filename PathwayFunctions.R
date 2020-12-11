# Pathway Normalization fuctions

####################################################################
# Function to calculate the probability that a cell is its neighbor
# Output: nCells x nCells matrix of pairwise probabilities of overlap
.CalculateDistance <- function(PCAEmbeddings, nCores, binSize){
  require(parallel)
  
  # Identify the distance from one point to another
  dists <- as.matrix(dist(PCAEmbeddings, 
                          method = dist.method))
  
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
  
  # Calculate probability matrix (AUB)/(AnB) (Joccard Index)
  # This is an extremely slow step, so validate if it is necessary
  p.mat <- matrix(nrow = length(dists),
                  ncol = length(dists),
                  dimnames = list(names(dists), names(dists)))
  
  # For each point, calculate the Joccard Index for all other points
  p.mat <- mclapply(dists, function(x) lapply(dists, function(y) sum(x%in%y)/length(unique(c(x,)))))
  p.mat <- mclapply(p.mat, unlist)
  p.mat <- do.call(rbind, p.mat)

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

####################################################################
# Function to plot the chromosomal coverage of the genes in the pathway genes
# Output: pdf plot of chromosomal coverage
.ChromosomeCoveragePlot <- function(pathways){
  
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

####################################################################
# Function to calculate ratio of pathway expression to average expression per chromosome
# Output: nCells x nCells matrix of pairwise probabilities of overlap
.AvgExpPerChromosome <- function(pathways, data){
  
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

####################################################################
# Function to generate null distributions to calculate pathway enrichment pvalues against
# Output: list of pathway enrichment null distributions

.GenNullEnrich <- function(pathways, data, nIter, nGenes){
  
  # Calculate non-zero median of each gene
  medianExpr <- apply(binnedCells, 1, median)
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
  randomControls <- mclapply(nGenes, function(x) .ctrlSim(x, nIter, medianExpr))
  randomControls <- do.call(cbind, randomControls)
  colnames(randomControls) <- paste0('NumGenes', nGenes)
  return(randomControls)
}