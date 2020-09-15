#!/bin/bash

# This is a pipeline for SNVs analysis using GATK
# It is the first step: Data pre-processing for variant discovery
# It contains:
    # 0. Prepare required files
    # 1. Alignment
    # 2. Process aligned reads to analysis-ready reads

# Author: PJ Wu
# Last update: 2020-09-15


# Reference: GATK Resource Bundle
# Species: Human
GENOME_PATH="PATH-to-GENOME-FASTA"
KNOWN_VARIANTS_PATH="PATH-to-KNOWN-VARIANTS"

# Dataset
DATA_PATH="PATH-to-FASTQ"

mkdir mapped_reads


# 0. Prepare required files: ref.fai, ref.dict, vcf.tbi
# Indexing reference genome
samtools faidx $GENOME_PATH 

# Create reference dictionary
gatk CreateSequenceDictionary \
    -R $GENOME_FA/Homo_sapiens_assembly38.fasta \
    -O $GENOME_FA/Homo_sapiens_assembly38.dict 

## Indexing known-sites VCF file
gatk IndexFeatureFile \
    -I $KNOWN_VARIANTS_PATH/1000G_phase1.snps.high_confidence.hg38.vcf.gz 
gatk IndexFeatureFile \
    -I $KNOWN_VARIANTS_PATH/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz 
gatk IndexFeatureFile \
    -I $KNOWN_VARIANTS_PATH/Homo_sapiens_assembly38.dbsnp138.vcf 


# 1. Align reads to reference genome
    # Build index
    # Align reads

## Build index
bwa index $GENOME_PATH/Homo_sapiens_assembly38.fasta 

## Align reads
bwa mem -t 20 \
    -M \
    $GENOME_PATH/Homo_sapiens_assembly38.fasta \
    $DATA_PATH/tumor_R1.clean.fastq.gz \
    $DATA_PATH/tumor_R2.clean.fastq.gz \
    -R '@RG\tID:ID\tLB:LABEL\tSM:SAMPLE\tPL:ILLUMINA' \
    2>./mapped_reads/tumor_bwa.log | samtools view -S -b -@ 20 -o ./mapped_reads/tumor.bam


#2. Process aligned to analysis-ready reads
    # Sort by coordinate and mark duplicates
    # Recalibrate base
    # Apply base quality score recalibration (BQSR)

## Sort and mark dup
picard SortSam \
    I=./mapped_reads/tumor.bam \
    O=./mapped_reads/tumor.sorted.bam \
    SORT_ORDER=coordinate
picard MarkDuplicates \
    I=./mapped_reads/tumor.sorted.bam \
    O=./mapped_reads/tumor.sorted.markDup.bam \
    M=./mapped_reads/tumor.sorted.markDup.txt 

## Recalibrate base
gatk BaseRecalibrator \
    -I ./mapped_reads/tumor.sorted.markDup.bam \
    -R $GENOME_PATH/Homo_sapiens_assembly38.fasta \
    --known-sites $KNOWN_VARIANTS_PATH/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
    --known-sites $KNOWN_VARIANTS_PATH/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
    --known-sites $KNOWN_VARIANTS_PATH/Homo_sapiens_assembly38.dbsnp138.vcf \
    -O ./mapped_reads/tumor.recal.table

gatk ApplyBQSR \
    -R $GENOME_PATH/Homo_sapiens_assembly38.fasta \
    -I ./mapped_reads/tumor.sorted.markDup.bam \
    --bqsr-recal-file ./mapped_reads/tumor.recal.table \
    -O ./mapped_reads/tumor.sorted.markDup.recal.bam \

gatk AnalyzeCovariates \
    -bqsr ./mapped_reads/tumor.recal.table \
    -plots ./oecm1.AnalyzeCovariates.pdf
