# Try to run CaSPeR
library(CaSpER)
library(Seurat)

source('Utils.R')
dirs <- setDirectory()
outputDir <- dirs[[1]]
dataDir <- dirs[[2]]
plottingDir <- dirs[[3]]

# Load the HCC dataset
data <- readRDS(paste0(outputDir, 'data_all.RDS'))

# Just for my janky file directory
data <- readRDS('data/TxToCloud/data_all.RDS')

# Code from the CaSpER Tutorial https://rpubs.com/akdes/673120
cytoband <- read.delim(paste0(dataDir, 'cytoBand.txt'), header = F)
cytoband <- data.frame(V1=gsub("chr", "", cytoband[,1]), V2=cytoband[,2], V3=cytoband[,3], V4=substring(cytoband$V4, 1, 1), stringsAsFactors=F)
start <- do.call(rbind, lapply(split(cytoband$V2, paste0(cytoband$V1, cytoband$V4)), min))
end <- do.call(rbind, lapply(split(cytoband$V3, paste0(cytoband$V1, cytoband$V4)), max))
cytoband <- data.frame(V1=gsub("p", "", gsub("q", "", rownames(start))), V2=start, V3=end, V4=rownames(start), stringsAsFactors=F)
cytoband <- cytoband [as.vector(unlist(sapply(c(1:22, "X"), function(x) which(cytoband$V1 %in% x)))), ]
cytoband$V4[grep("q", cytoband$V4)] <- "q"
cytoband$V4[grep("p", cytoband$V4)] <- "p"
rownames(cytoband) <- NULL

annotation <- generateAnnotation(id_type="hgnc_symbol", 
                                  genes=rownames(data), 
                                  ishg19=T, 
                                  centromere,
                                  host="useast.ensembl.org")

# Trim data to fit annotation
data <- data[match(annotation$Gene, rownames(data)),]

loh <- readBAFExtractOutput(dataDir,
                            sequencing.type = 'single-cell')

# Tells which loh file each cell is attributable to - aka only 1 file
names(loh) <- ('HCC')
loh.mapping <- data.frame(loh.name = 'HCC', sample.name = colnames(data))

# Make the CaSpER object
norm.reads <- as.matrix(data@assays$RNA@data)
ctrl.samps <- colnames(data)[data$cell.type!='HCC']

casper.obj <- CreateCasperObject(raw.data = norm.reads,
                                 annotation = annotation, # Cell type annotation
                                 matrix.type = 'normalized',
                                 expr.cutoff = 0.1, # Average expression below which you throw out a gene
                                 loh = loh,
                                 loh.name.mapping = loh.mapping,
                                 loh.scale = 3,
                                 cnv.scale = 3,
                                 filter = 'median',
                                 control.sample.ids = ctrl.samps,
                                 method = 'iterative',
                                 cytoband = cytoband,
                                 sequencing.type = 'single-cell')

gc()
final.objects <- runCaSpER(casper.obj, 
                           removeCentromere=F, 
                           method="iterative")
