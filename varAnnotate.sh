#!/bin/bash


# Author: Pin-Jou Wu
# Last update: 2021-01-26

# Annotate the vcf file by VEP (v101.0)
# The genome fasta file used here was downloaded from Ensembl: homo_sapiens.grch38.dna.primary_assembly.fa.gz
# If set --cache, cache files should be downloaded and decompressed in ~/.vep (default) or specific the path of the cache
# The version of cache files should match to the version of VEP
# Cache: GRCh38/hg38

# Reference
GENOME_PATH="PATH-to-GENOME-FASTA"
CACHE_PATH="~/.vep"

# Dataset
DATA_PATH="./hapCaller"
OUTPUT="./annotated"
mkdir annotated

vep -i $DATA_PATH/tumor.final.vcf.gz \
    --cache \
    --dir_cache $CACHE_PATH \
    --species homo_sapiens \
    -e \
    --fork 8 \
    --offline \
    --format vcf \
    --fasta  $GENOME_PATH/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz \
    -o $OUTPUT/tumor.variant_effect.vcf

# Extract the info in EXTRA column in the annotated vcf
# The scirpt CleanUPAnnoVariant.py is written by CYC Lab.
python CleanUPAnnoVariant.py -i $OUTPUT/tumor.variant_effect.vcf -o $OUTPUT/tumor.variant_effect.clean.vcf
