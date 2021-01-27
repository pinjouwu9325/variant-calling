#!/bin/basha


# Author: PJ Wu
# Last update: 2021-01-26

# This script is for single-sample data variant calling using GATK (HaplotypeCaller)

# Reference
GENOME_PATH="PATH-to-GENOME-FASTA"
KNOWN_VARIANTS_PATH="PATH-to-KNOWN-VARIANTS"

# Dataset
DATA_PATH="PATH-to-ANAYSIS-READY-BAM"
INTERVALS="PATH-to-INTERVALS"
OUTPUT="./hapCaller"
mkdir hapCaller

# Variant calling with snp reference, intervals and bamout to show realigned reads
# --dbsnp 1303/mgp.v3.snps.rsIDdbSNPv137.vcf.gz
gatk --java-options "-xmx4g" haplotypecaller \
    --dbsnp $KNOWN_VARIANTS_PATH/Homo_sapiens_assembly38.dbsnp138.vcf \
    --intervals $INTERVALS/intervals.list \
    -R $GENOME_PATH/Homo_sapiens_assembly38.fasta \
    -I $DATA_PATH/tumor_sorted_markdup_recal.bam \
    -O $OUTPUT/tumor.raw.vcf.gz \
    -bamout $OUTPUT/tumor_Hapbamout.bam
