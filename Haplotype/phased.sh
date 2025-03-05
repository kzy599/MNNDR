#!/bin/bash

wd="/home/kangziyi/poster/callSNP/"
Genome="GCF_963853765.1_xbMagGiga1.1_genomic.fna"
dicname="GCF_963853765.1_xbMagGiga1.1_genomic.dict"
genomedir="/home/kangziyi/poster/genome/"


bcftools +fixploidy "${wd}BasePoyster.snp.final.refiltered.recode.vcf" \
    -- -f 2 >"${wd}BasePoyster.snp.final.fixed.vcf"


beagle.jar \
    gt="${wd}BasePoyster.snp.final.fixed.vcf" \
    out="${wd}BasePoyster.snp.final.phased" \
    nthreads=40

gunzip BasePoyster.snp.final.phased.vcf.gz

