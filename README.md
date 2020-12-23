# Pathway Normalization Based CNV Caller

This repository contains code for our pathway based CNVCaller and implentations of HoneyBADGER and inferCNV. The Pathway Based CNVCaller is a supplement to the existing inferCNV package which refines the program’s ability to identify copy number variations in scRNA-Seq datasets. This feature utilizes a simplified gene set enrichment algorithm to identify and normalize overexpressed gene expression pathways.


##Data and Package Availability##

Data is available via a google bucket with the following command: `gsutil cp gs://compgenomics_sccnvcaller/data.zip` (requires [gsutil](https://cloud.google.com/storage/docs/gsutil_install))

All of these files are run on R 4.0.3. Required packages and their versions can be found in the `sessionInfo.txt` file in the attached data

In order to make the program work on multiple systems, all paths are relative to:

- dataDir: directory with raw, unprocessed files
- outDir: directory with processed files, and outputs from scripts
- plottingDir: directory to place stand alone plots

We have provided a convience function, `setDirectory()` within `Utils.R` to automatically set those directories to suit your local filesystem. 
Please edit it such that it addresses the datDir and outDir in the zipped data file, and a convenient plotting directory of your choice.  

##Usage##


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
