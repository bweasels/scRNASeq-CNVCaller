---
title: "README"
output: html_document
---

This repository contains code for scRNASeq-CNVCaller and for our implentations of HoneyBADGER and inferCNV. scRNASeq-CNVCaller is a supplement to the existing inferCNV package which refines the program’s ability to identify copy number variations in scRNA-Seq datasets. This feature utilizes a simplified gene set enrichment algorithm to identify and normalize overexpressed gene expression pathways.

The scRNA-Seq datasets we have tested are available at:
- link to liver dataset
- link to melanoma dataset

**Usage**

As a prerequisite, you must have R installed to use this repository.

Input files for the scRNA-Seq dataset consist of: 

-
- 
- 
-

All input files should be added to the inputs folder. 
During inference, output files will be generated in the outputs folder.

Inference takes 15 min to run on a cloud server and ~ 90 min on a local computer.

You can run the different algorithms using the following commands.

First, preprocess the dataset:

Run scRNASeq-CNVCaller:
```
Rscript CNVCaller/RunCNVcaller_AllCells.R
```

Run inferCNV:
```
Rscript inferCNV/RunInferCNV_AllCells.R
```

Run HoneyBADGER:
```
Rscript HoneyBADGER/RunHoneyBADGER_AllCells.R
```

#TODO: add what each file is doing
#TODO: add inputs & outputs folders to git with a filtered example
#TODO: add steps to perform to obtain all input files
#TODO: point to dataDir, outDir, plottingDir in Utils   