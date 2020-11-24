#!/bin/bash
# This will dl the raw 10x data and break into fastqs

OUT_DIR=/home/bkw2118/extDrive/scRNASeq-CNVCaller/data/raw10XFiles/
GEO_URL=https://sra-pub-src-2.s3.amazonaws.com/SRR6885508/Pt14.d_possorted_genome_bam.bam.1
BAM_FILE=/home/bkw2118/extDrive/scRNASeq-CNVCaller/data/raw10XFiles/Pt14.d.raw.bam
FASTQ_DIR=/home/bkw2118/extDrive/scRNASeq-CNVCaller/data/raw10XFiles/fastqs/

# Download the bam file
curl $GEO_URL --output $BAM_FILE

# Run bcltofastq on it all
bamtofastq-1.3.2 --nthreads=16 $BAM_FILE $FASTQ_DIR

REFERENCE=/mnt/disks/ext/tools/cellranger-5.0.0/reference/refdata-gex-GRCh38-2020-A
ID=CNVCALLER
SAMPLE=PT14D

cellranger count --id=$ID --transcriptome=$REFERENCE --fastqs=$FASTQ_DIR
exit
