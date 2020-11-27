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

    }else if(info['user']=='charlotterochereau'){
    #setwd('~/Code/scRNASeq-CNVCaller')
    outDir <- 'output/'
    dataDir <- 'raw_data/'
    plottingDir <- 'plots/'
    
    # Turns out adding storage adds an external disk, so I had to mount it and move everything there
    # I think we both should have permissions to mess with it so I'll set it as the default
    }else if(info['release']=='4.19.0-12-cloud-amd64'){
    setwd('/mnt/disks/ext/scRNASeq-CNVCaller')
    outDir <- '/mnt/disks/ext/scRNASeq-CNVCaller/output/'
    dataDir <- '/mnt/disks/ext/scRNASeq-CNVCaller/data/'
    plottingDir <- '/mnt/disks/ext/scRNASeq-CNVCaller/plots/'
  }
  
  return(list(outDir, dataDir, plottingDir))
}