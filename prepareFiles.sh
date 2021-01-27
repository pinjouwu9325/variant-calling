#!/bin/bash


# Author: Pin-Jou Wu
# Last update: 2021-01-26

# Before running GATK variant calling pipelines, there are several files need to be prepared. They are: 
# 1. Reference genome (.fasta)
# 2. Indexed reference genome (.fasta.fai)
# 3. Reference dictionary (.dict)
# 4. Indexed feature files (.vcf.gz.tbi)

# They can be downloaded from GATK Resource bundle or created using the following commands.
# GATK Resouce bundle:
# https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0/
# gs://gcp-public-data--broad-references (using Google Cloud SDK to download)

GENOME_PATH="PATH-to-GENOME-FASTA"
KNOWN_VARIANTS_PATH="PATH-to-KNOWN-VARIANTS"

# Index reference genome
samtools faidx $GENOME_PATH/Homo_sapiens_assembly38.fasta 

# Create reference dictionary
gatk CreateSequenceDictionary \
    -R $GENOME_PATH/Homo_sapiens_assembly38.fasta \
    -O $GENOME_PATH/Homo_sapiens_assembly38.dict 

# Index feature files
gatk IndexFeatureFile \
    -I $KNOWN_VARIANTS_PATH/1000G_phase1.snps.high_confidence.hg38.vcf.gz
gatk IndexFeatureFile \
    -I $KNOWN_VARIANTS_PATH/Homo_sapiens_assembly38.dbsnp138.vcf 
gatk IndexFeatureFile \
    -I $KNOWN_VARIANTS_PATH/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz 
gatk IndexFeatureFile \
    -I $KNOWN_VARIANTS_PATH/1000g_pon.hg38.vcf.gz 
gatk IndexFeatureFile \
    -I $KNOWN_VARIANTS_PATH/af-only-gnomad.hg38.vcf.gz 

# Build reference genome index for BWA
bwa index $GENOME_PATH/Homo_sapiens_assembly38.fasta 
