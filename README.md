# Pathway Normalization Based CNV Caller

This repository contains code for our pathway based CNVCaller and implentations of HoneyBADGER and inferCNV. The Pathway Based CNVCaller is a supplement to the existing inferCNV package which refines the program’s ability to identify copy number variations in scRNA-Seq datasets. This feature utilizes a simplified gene set enrichment algorithm to identify and normalize overexpressed gene expression pathways.


## Data and Package Availability

Data is available via a google bucket: `gs://compgenomics_sccnvcaller/data.zip` (access requires [gsutil](https://cloud.google.com/storage/docs/gsutil_install))

All of these files are run on R 4.0.3. Required packages and their versions can be found in the `sessionInfo.txt` file in top level.

In order to make the program work on multiple systems, all paths are relative to:

- dataDir: directory with raw, unprocessed files (see example on small dataset)
- outDir: directory with processed files, and outputs from scripts (see example on small dataset)
- plottingDir: directory to place stand alone plots

We have provided a convience function, `setDirectory()`, within `Utils.R` to automatically set those directories to suit your local filesystem. 
There is a miniaturized version of the test data in the correct out and dataDirs in the github, and full sized ones in the zipped data file from the google bucket. 
Please edit `setDirectory()` to address the dataDirs and outDirs you want to use, as well as a convenient plotting directory of your choice.  

Like inferCNV, scRNASeq-CNVCaller require three inputs:
- a read count matrix 
- a gmt file of gene sets to analyze
- a gene ordering file 


## Usage


### Dataset preprocessing
In the preprocessing directory, we provide scripts to filter out bad quality data and to annotate the remaining data.
We use FilterMelanoma.R for the melanoma dataset (*quote*) and FilterLiverCancer.R for the liver dataset (*quote*).
The preprocessed data is saved as an RDS object (data_all.RDS) in outDir.
```
Rscript preprocessing/FilterLiverCancer.R
```

We also simulate a null dataset to provide a negative control. To generate it, run generateSplatterDataset.R.
The null dataset is also saved as a RDS object (data_sim.RDS) in outDir.
```
Rscript preprocessing/generateSplatterDataset.R
```

### Inference
Inference takes 15 min to run on a cloud server and ~ 90 min on a local computer.
You can run the different algorithms using the following commands.

#### scRNASeq-CNVCaller
We describe scRNASeq-CNVCaller algorithm in CNVcaller_Algo.R. Use RunCNVcaller_AllCells.R to run the algorithm on your preprocessed dataset.

```
Rscript CNVCaller/RunCNVcaller_AllCells.R
```
CNVcaller subdirectories include:

- inferCNVscripts - files from the original inferCNV package that are required to run scRNASeq-CNVCaller
- PathwayDevFiles - scripts used for developing the new pathway normalization feature
- CNVcaller_Refactoring - scripts for streamlining the code of the current scRNASeq-CNVCaller implementation (WIP)


#### Run inferCNV
```
Rscript inferCNV/RunInferCNV_AllCells.R
```
We use separate scripts for each of the datasets we tested (melanoma dataset, null dataset, smaller liver cancer dataset, etc). All are modified versions of RunInferCNV_AllCells.R.


#### Run HoneyBADGER
```
Rscript HoneyBADGER/RunHoneyBADGER_AllCells.R
```
We use separate scripts for each of the datasets we tested. All are modified versions of RunHoneyBADGER_AllCells.R.
RunHoneyBADGER_AllCellsBAMReference.R was implemented using BAM files generated from Fastq by CellRanger and VCF files generated using bcftools.


#### Run CaSpER
```
Rscript CaSpER/RunCaSpER_AllCells.R
```