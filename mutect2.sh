#!/bin/bash


# Author: Pin-Jou Wu
# Last update: 2021-01-26

# This script is for somatic variant calling (GATK-Mutect2)
# It takes the processed BAMs from dataPreprocess.sh

# The pipeline contains:
    # 1. Call candidate variants
    # 2. Filter variants 

# Reference (Human)
GENOME_PATH="PATH-to-GENOME-FASTA"
GERMLINE_RESOURCE="PATH-to-GERMLINE-RESOURCE"
PON="PATH-to-PON"

# Dataset (e.g. ./mapped_reads)
DATA_PATH="PATH-to-ANALYSIS-READY-BAM"
OUTPUT="./mutect2"
mkdir mutect2


#1. Somatic variants calling (Mutect2)
gatk Mutect2 \
    -R $GENOME_PATH/Homo_sapiens_assembly38.fasta \
    -I $DATA_PATH/tumor.sorted.markDup.recal.bam \
    --germline-resource $GERMLINE_RESOURCE/af-only-gnomad.hg38.vcf.gz \
    --panel-of-normals $PON/1000g_pon.hg38.vcf.gz \
    --f1r2-tar-gz $OUTPUT/f1r2.tar.gz \
    -O $OUTPUT/tumor.raw.vcf.gz

gatk LearnReadOrientationModel \
    -I $OUTPUT/f1r2.tar.gz \
    -O $OUTPUT/f1r2.priors.tar.gz

#2. Filter variants
gatk FilterMutectCalls \
    -R $GENOME_PATH/Homo_sapiens_assembly38.fasta \
    -V $OUTPUT/tumor.raw.vcf.gz \
    --ob-priors $OUTPUT/f1r2.priors.tar.gz \
    -O $OUTPUT/tumor.filtered.vcf.gz 

# Select variants which pass the filter (FILTER==PASS)
bcftools view -f PASS -Oz $OUTPUT/tumor.filtered.vcf.gz -o $OUTPUT/tumor.final.vcf.gz
