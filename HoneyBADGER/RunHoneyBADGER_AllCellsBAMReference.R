# Load function
source('Utils.R')

# Prep for HoneyBADGER on the scRNASeq dataset
dirs <- setDirectory()
outDir <- dirs[[1]]
dataDir <- dirs[[2]]
plottingDir <- dirs[[3]]
date <- Sys.Date()

library(HoneyBADGER)
library(GenomicRanges)
library(VariantAnnotation)
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
cellBarcodes <- paste0(dataDir, 'cellRangerOutput/outs/filtered_feature_bc_matrix/barcodes.tsv.gz')

results <- getSnpMats10X(snps = snps,
                         bamFile = bamFile,
                         indexFile = baiFile,
                         barcodes = cellBarcodes)

# Get the values required for that
r <- results$refCount
cov <- results$refCount + results$altCount

# load transcripts to map SNPs to genes
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
