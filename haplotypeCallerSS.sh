#! /bin/bash

# This script is for single-sample data variant calling using GATK (HaplotypeCaller).

# Author: PJ Wu
# Last update: 2020-09-15

# Reference
GENOME_PATH="PATH-to-GENOME-FASTA"
DBSNP="PAHT-to-DBSNP"

# Dataset
DATA_PATH="PATH-to-ANAYSIS-READY-BAM"
INTERVALS="PATH-to-INTERVALS"

# Variant calling with snp reference, intervals and bamout to show realigned reads
# --dbsnp 1303/mgp.v3.snps.rsIDdbSNPv137.vcf.gz
gatk --java-options "-xmx4g" haplotypecaller \
    --dbsnp $DBSNP/mgp.v3.snps.rsiddbsnpv137.vcf.gz \
    --intervals $INTERVALS/mm_intervals.list \
    -r $GENOME_PATH/mus_musculus.grcm38.dna.primary_assembly.fa \
    -i $DATA_PATH/tumor_sorted_markdup_recal.bam \
    -o tumor.vcf.gz \
    -bamout tumor_Hapbamout.bam

# Variant calling with snp reference, intervals and bamout to show realigned reads
# --dbsnp 1807/mgp.v6.merged.norm.snp.indels.sfiltered.vcf.gz
# Failed because of "The provided VCF file is malformed"
# gatk --java-options "-Xmx4g" HaplotypeCaller \
    # --dbsnp /LVM_data/pinjouwu9325/ref/mm10/sanger/1807/mgp.v6.merged.norm.snp.indels.sfiltered.vcf.gz \
    # --intervals ./intervals.list \
    # -R /LVM_data/pinjouwu9325/ref/mm10/ensembl/Mus_musculus.GRCm38.dna.primary_assembly.fa \
    # -I ../mapped_reads/grcm38/MTCR1_sorted_markDup_recal.bam \
    # -O MTCQ1_wMergedSnpRef.vcf.gz \
  #   -bamout MTCQ1_haplo_bamout.bamout
