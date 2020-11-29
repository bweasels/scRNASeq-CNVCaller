#!/bin/bash
#Script to generate the vcf file for HoneyBADGER

# Variables for the scripts
FASTA=/home/bkw2118/extDrive/tools/cellranger-5.0.0/reference/refdata-gex-GRCh38-2020-A/fasta/genome.fa
BAM=/mnt/disks/ext/scRNASeq-CNVCaller/data/cellRangerOutput/outs/possorted_genome_bam.bam
OUT=/mnt/disks/ext/scRNASeq-CNVCaller/data/pt14.d.vcf.gz
OUT_AF=/mnt/disks/ext/scRNASeq-CNVCaller/data/pt14.d.vcf.filtAF.gz

# Create the raw vcf file
bcftools mpileup -Ou -f $FASTA $BAM | bcftools call -vmO z -o $OUT

# Filter the variant calls
bcftools filter -O z -o $OUT_FILTERED -s LOWQUAL -i'%QUAL>10' $OUT

# Add all the tags for HoneyBADGER
bcftools +fill-tags $OUT_FILTERED -o $OUT_AF -- -t AF
