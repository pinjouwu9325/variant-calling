#!/bin/bash


# Author: Pin-Jou Wu
# Last update: 2021-01-26

# This script is for Variant Quality Score Recalibration (VQSR) in GATK germline variant calling pipeline
# It contains:
#   1. Build the models: VariantRecalibrator
#   2. Apply a filtering threshold: ApplyRecalibration
#   3. Select variants which pass the filter: bcftools 
# Recalivrate variant types (SNP, INDEL) separately

# Reference
GENOME_PATH="PATH-to-GENOME-FASTA"
KNOWN_VARIANTS_PATH="PATH-to-KNOWN-VARIANTS"

# Dataset and output
DATA_PATH="./hapCaller"
OUTPUT_VQSR="./vqsr"
mkdir vqsr

# SNP
# 1. Build the models
gatk VariantRecalibrator \
    -R $GENOME_PATH/Homo_sapiens_assembly38.fasta \
    -V $DATA_PATH/tumor.raw.vcf.gz \
    -resource:hapmap,known=false,training=true,truth=true,prior=15.0 $KNOWN_VARIANTS_PATH/hapmap_3.3.hg38.vcf.gz \
    -resource:omni,known=false,training=true,truth=true,prior=12.0 $KNOWN_VARIANTS_PATH/1000G_omni2.5.hg38.vcf.gz \
    -resource:1000G,known=false,training=true,truth=false,prior=10.0 $KNOWN_VARIANTS_PATH/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
    -resource:dbsnp,known=true,training=false,truth=false,prior=7.0 $KNOWN_VARIANTS_PATH/Homo_sapiens_assembly38.dbsnp138.vcf \
    -an QD -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an SOR -an DP \
    -mode SNP \
    --max-gaussians 6 \
    -O $OUTPUT_VQSR/tumor.snp.recal \
    --tranches-file $OUTPUT_VQSR/tumor.snp.tranches \
    --rscript-file $OUTPUT_VQSR/tumor.snp.plots.R

# 2. Apply a filtering threshold
gatk ApplyVQSR \
    -R $GENOME_PATH/Homo_sapiens_assembly38.fasta \
    -V $DATA_PATH/tumor.raw.vcf.gz \
    -O $OUTPUT_VQSR/tumor.snp.vqsr.vcf.gz \
    --truth-sensitivity-filter-level 99.0 \
    --tranches-file $OUTPUT_VQSR/tumor.snp.tranches \
    --recal-file $OUTPUT_VQSR/tumor.snp.recal \
    -mode SNP 


# INDEL
# 1. Builf the models (without -an MQ)
gatk VariantRecalibrator \
    -R $GENOME_PATH/Homo_sapiens_assembly38.fasta \
    -V $DATA_PATH/tumor.raw.vcf.gz \
    -resource:mills,known=false,training=true,truth=true,prior=12.0 $KNOWN_VARIANTS_PATH/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
    -resource:axiomPoly,known=false,training=true,truth=false,prior=10.0 $KNOWN_VARIANTS_PATH/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz \
    -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $KNOWN_VARIANTS_PATH/Homo_sapiens_assembly38.dbsnp138.vcf \
    -an QD -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an DP \
    -mode INDEL \
    --max-gaussians 4 \
    -O $OUTPUT_VQSR/tumor.indel.recal \
    --tranches-file $OUTPUT_VQSR/tumor.indel.tranches \
    --rscript-file $OUTPUT_VQSR/tumor.indel.plots.R

# 2. Apply a filtering threshold
# Given the snp-filtered callset, this step results in the final filtered callset.
gatk ApplyVQSR \
    -R $GENOME_PATH/Homo_sapiens_assembly38.fasta \
    -V $OUTPUT_VQSR/tumor.snp.vqsr.vcf.gz \
    -O $OUTPUT_VQSR/tumor.indel.vqsr.vcf.gz \
    --truth-sensitivity-filter-level 99.0 \
    --tranches-file $OUTPUT_VQSR/tumor.indel.tranches \
    --recal-file $OUTPUT_VQSR/tumor.indel.recal \
    -mode INDEL

# The --max-gaussians parameter sets the expected number of clusters in modeling. 
# If a dataset gives fewer distinct clusters, e.g. as can happen for smaller data, then the tool will tell you there is insufficient data with a No data found error message. In this case, try decrementing the --max-gaussians value.

# 3. Select the variantsw which pass the filter of VQSR (FILTER==PASS)
bcftools view -f PASS -Oz $OUTPUT_VQSR/tumor.indel.vqsr.vcf.gz -o $DATA_PATH/tumor.final.vcf.gz
