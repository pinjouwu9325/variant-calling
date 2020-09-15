#!/bin/bash

# This script is for somatic variant calling (GATK-Mutect2)
# It takes the processed BAMs from dataPreprocess.sh

# Author: PJ Wu
# Last update: 2020-09-15

# The pipeline contains:
    # 1. Call candidate variants
    # 2. Filter Variants 

# Reference (Human)
GENOME_PATH="PATH-to-GENOME-FASTA"
GERMLINE_RESOURCE="PATH-to-GERMLINE-RESOURCE"
PON="PATH-to-PON"

# Dataset
DATA_PATH="PATH-to-ANALYSIS-READY-BAM"
    

#1. Somatic variants calling (Mutect2)
gatk Mutect2 \
    -R $GENOME_PATH/Homo_sapiens_assembly38.fasta \
    -I ./mapped_reads/tumor.sorted.markDup.recal.bam \
    --germline-resource $GERMLINE_RESOURCE/af-only-gnomad.hg38.vcf.gz \
    --panel-of-normals $PON/1000g_pon.hg38.vcf.gz \
    --f1r2-tar-gz f1r2.tar.gz \
    -O tumor.vcf.gz

gatk LearnReadOrientationModel \
    -I f1r2.tar.gz \
    -O f1r2.priors.tar.gz

#2. Filter variants
gatk FilterMutectCalls \
    -R $GENOME_PATH/Homo_sapiens_assembly38.fasta \
    -V ./tumor.vcf.gz \
    --ob-priors f1r2.priors.tar.gz \
    -O tumor.filtered.vcf.gz 

# Filter vcf file with PASS and chr regions
bcftools view -f PASS tumor.filtered.vcf.gz -O z > tumor.PASS.vcf.gz
