#!/bin/bash

# Annotate variants by VEP (version 100)
# If set --cache, cache files should be downloaded and decompressed in ~/.vep.
# The version of cache files should match to the version of VEP
vep -i tumor.PASS.vcf.gz --cache --species homo_sapiens -e --fork 20 -o tumor.PASS.variant_effect.txt
