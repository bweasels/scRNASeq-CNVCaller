# Load function
source('Utils.R')

# Prep for HoneyBADGER on the scRNASeq dataset
dirs <- setDirectory()
outDir <- dirs[[1]]
dataDir <- dirs[[2]]
plottingDir <- dirs[[3]]

date <- Sys.Date()

library(HoneyBADGER)
library(Seurat)
library(GenomicRanges)
library(VariantAnnotation)
library(biomaRt)

# Load data, and the start stop locations of each transcript
data <- readRDS(paste0(outDir, 'dataFiltered.RDS'))

# Trim tx in data to match tx I was able to find (~500 genes different - all ribosomal varants)
sampleAnnotation <- data@meta.data
sampleAnnotation <- data.frame(row.names = rownames(sampleAnnotation), 
                               cellType = sampleAnnotation$cell.type)

# Get the unique cell types and remove HCC so we have a list of reference cell types
ref_groups <- unique(sampleAnnotation$cellType)
ref_groups <- ref_groups[!grepl('HCC', ref_groups)]

HCC <- data[,data$cell.type=='HCC']
HCC <- as.matrix(HCC@assays$RNA@counts)

ref <- data[,data$cell.type!='HCC']
ref <- as.matrix(ref@assays$RNA@counts)

# Get the gene positions
mart.obj <- useMart(biomart = 'ENSEMBL_MART_ENSEMBL', 
                    dataset = 'hsapiens_gene_ensembl')

# make honeybadger object
hb <- new('HoneyBADGER', name='Liver_Cancer')
hb$setGexpMats(gexp.sc.init = HCC, 
               gexp.ref.init = ref, 
               mart.obj = mart.obj, 
               filter = F, 
               scale = T, 
               verbose = T)

# Model variance and necessary exp deviance to ID CNV
hb$setMvFit(verbose = T)
hb$setGexpDev(verbose=T)
hb$calcGexpCnvBoundaries(init=T,
                         verbose = T)

# Output the identified CNVs
print(range(hb$genes[unlist(hb$bound.genes.final)]))

# Retest to ID posterior probability & get final results
hb$retestIdentifiedCnvs(retestBoundGenes = T,
                        retestBoundSnps = F,
                        verbose = T)

results <- hb$summarizeResults(geneBased = T,
                               alleleBased = F)

# Make result plot
pdf(paste0(plottingDir, 'HoneyBADGERRes_AllCells_', date, '.pdf'))
trees <- hb$visualizeResults(geneBased = T,
                             alleleBased = F,
                             details = T,
                             margins = c(25, 15))

hc <- trees$hc
order <- hc$labels[hc$order]

hb$plotGexpProfile()
hb$plotGexpProfile(region = hb$cnvs[['gene-based']][['amp']])
hb$plotGexpProfile(region = hb$cnvs[['gene-based']][['del']])
dev.off()

# Add allele level information
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

##PREP SNP FILES
# Load vcfFile
vcfFile <- paste0(dataDir, 'pt14.d.vcf.filtAF.gz')
vcf <- readVcf(vcfFile, 'hg38')

# Get snps
snps <- rowRanges(vcf)

#Get allele frequency of each alternate allele
info <- info(vcf)
maf <- info[,'AF']

# Limit to common snps by MAF (>10% in pop)
vi <- sapply(maf, function(x) any(x>0.1))
snps <- snps[vi,]

# Get rid of non-snps
vi <- width(snps@elementMetadata$REF) == 1
snps <- snps[vi,]


## Make SNP Mats
bamFile <- paste0(dataDir, 'cellRangerOutput/outs/possorted_genome_bam.bam')
baiFile <- paste0(dataDir, 'cellRangerOutput/outs/possorted_genome_bam.bam.bai')
barcodes <- paste0(dataDir, 'cellRangerOutput/outs/filtered_feature_bc_matrix/barcodes.tsv.gz')

results <- getSnpMats10X(snps = snps,
                         bamFile = bamFile,
                         indexFile = baiFile,
                         cellBarcodes = barcodes)

# Get the values required for that
r <- results$refCount
cov <- results$refCount + results$altCount

# load transcripts to map SNPs to genes
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# Use the allele metrics
hb$setAlleleMats(r.init=r,
                 n.sc.ini=cov,
                 het.deviance.threshold=0.1,
                 n.cores = detectCores())
hb$setGeneFactors(txdb)

# Plot allele level snps
pdf(paste0(plottingDir, 'HoneyBADGERRes_AllCellsAllele_', date, '.pdf'))
hb$plotAlleleProfile()
hb$calcAlleleCnvBoundaries(init=T, 
                           verbose=F)

# double check what CNVs were identified
bsf <- get('bound.snps.final', slot(hb, '.xData'))
snps <- get('snps', slot(hb, '.xData'))
regions.snp <- range(snps[unlist(bsf)])
print(regions.snp)

# Retest for the posteriori probability
hb$retestIdentifiedCnvs(retestBoundGenes=FALSE, 
                        retestBoundSnps=TRUE, 
                        verbose=FALSE)
results <- hb$summarizeResults(geneBased=FALSE, 
                               alleleBased=TRUE)
print(head(results[,1:5]))

# Plot allele level heatmap
trees2 <- hb$visualizeResults(geneBased=FALSE, 
                              alleleBased=TRUE, 
                              details=TRUE)
hb$plotAlleleProfile()
hb$plotAlleleProfile(region=hb$cnvs[['allele-based']][['del.loh']])
res.allele <- res.gene <- hb$summarizeResults(geneBased=T,
                                              alleleBased=T)
saveRDS(paste0(outDir, 'HoneyBADGERAlleleres_', date, '.RDS'))
