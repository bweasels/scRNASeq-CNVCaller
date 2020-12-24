setDirectory <- function(){
  # Returns info about the machine and User login
  info <- Sys.info()
  
  # Local directory on my PC
  if(info['user']=='benku'){
    # Data kept on an external hdd for backup
    setwd('E:')
    
    setwd('/OneDrive/PhD/Fall 2020/Computational Genomics/scRNASeq-CNVCaller/')
    outDir <- '/OneDrive/PhD/Fall 2020/Computational Genomics/scRNASeq-CNVCaller/data/'
    dataDir <- '/OneDrive/PhD/Fall 2020/Computational Genomics/scRNASeq-CNVCaller/Liver Cancer/'
    plottingDir <- '/OneDrive/PhD/Fall 2020/Computational Genomics/scRNASeq-CNVCaller/plots/'

    }else if(info['user']=='charlotterochereau'){
    #setwd('~/Code/scRNASeq-CNVCaller')
    outDir <- 'output/'
    dataDir <- 'raw_data/'
    plottingDir <- 'plots/'
    
    # Add storage adds an external disk. Mount it and move everything there
    }else if(info['release']=='4.19.0-13-cloud-amd64'){
    setwd('/mnt/disks/ext/scRNASeq-CNVCaller')
    outDir <- '/mnt/disks/ext/scRNASeq-CNVCaller/output/'
    dataDir <- '/mnt/disks/ext/scRNASeq-CNVCaller/data/'
    plottingDir <- '/mnt/disks/ext/scRNASeq-CNVCaller/plots/'
  }
  
  return(list(outDir, dataDir, plottingDir))
}
