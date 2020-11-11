## Generate our test datasets with Splatter
#Ben's Local Directory
setwd('/OneDrive/PhD/Fall 2020/Computational Genomics/scRNASeq_CNV/')
outDir <- '/OneDrive/PhD/Fall 2020/Computational Genomics/data/'
dataDir <- '/OneDrive/PhD/Fall 2020/Computational Genomics/scRNASeq_CNV/Liver Cancer/Pt13.a/'

###For Charlotte###
#There are two package mangers in R world
#1) install.packages('your package') is in-built and generally my first go to
#2) Bioconductor is a bioinformatics focused repository of packages
#Most packages are listed on both, but there are many bioinformatic tools that can only
#be found on Bioconductor, so it is essential to have too

##To install bioconductor##
if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}

##To install packages with bioconductor
BiocManager::install('splatter')
#The :: operator asks R to temporarily load that package and use the function afterwards. 

##Other tips
#If you have imported two packages each with a function with identical names, 
#R will use your provided arguments to automatically select the correct funcion
##For example given functions A::multiply(number1, number2, decimalPlace), and B::multiply(number1, number2, returnMessage)
#If you do multiply(1, 2, 'done!'), R will automatically detect that you are using package B's multiply

#In addition, R will ~generally~ figure out your arguments without specifying the exact argument that your variable is passed to
#So, if you do multiply('done!', 1, 2), R will figure you want package B's multiply and will ~generally~ map the string 'done!' 
#to variable returnMessage


#But it is good practice to be explicit in your variable mapping
#For example
#multiply(number1 = 1,
#         number2 = 2,
#         returnMessage = 'Done!')

#R's version of import
library(splatter)
library(Seurat)

#set.seed locks our random number indexer so we get repeatable results
set.seed(10000)

#Seurat's function to load 10X scRNA-Seq data
#Removes genes that are expressed in 3 or fewer cells, and cells that express fewer than 200 unique genes
data <- Read10X(data.dir = dataDir)
data <- CreateSeuratObject(counts = data,
                           project = 'HCC',
                           min.cells = 3,
                           min.features = 200)
