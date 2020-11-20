# Put all utility functions here

setDirectory <- function(){
  # Returns info about the machine and User login
  info <- Sys.info()
  
  # Local directory on my PC
  if(info['user']=='benku'){
    # I keep my data on an external hdd for backup
    setwd('E:')
    
    setwd('/OneDrive/PhD/Fall 2020/Computational Genomics/scRNASeq-CNVCaller/')
    outDir <- '/OneDrive/PhD/Fall 2020/Computational Genomics/scRNASeq-CNVCaller/data/'
    dataDir <- '/OneDrive/PhD/Fall 2020/Computational Genomics/scRNASeq-CNVCaller/Liver Cancer/Pt13.a/'
    plottingDir <- '/OneDrive/PhD/Fall 2020/Computational Genomics/scRNASeq-CNVCaller/plots/'
    
    # This is my local directory on our google cloud server - you might have to set up your own if the user permissions don't play well
  }else if(info['user']=='bkw2118'){
    setwd('/home/bkw2118/scRNASeq-CNVCaller')
    outDir <- 'outDir/'
    dataDir <- 'dataDir/'
    plottingDir <- 'plots/'
  }
  
  return(list(outDir, dataDir, plottingDir))
}