#!/bin/bash

wd="/home/kangziyi/poster/callSNP/"
Genome="GCF_963853765.1_xbMagGiga1.1_genomic.fna"
dicname="GCF_963853765.1_xbMagGiga1.1_genomic.dict"
genomedir="/home/kangziyi/poster/genome/"

gatk SelectVariants -V "${wd}BasePoyster.raw.vcf.gz" -select-type SNP -O "${wd}BasePoyster.raw.snp.vcf.gz"


gatk VariantFiltration \
  -V "${wd}BasePoyster.raw.snp.vcf.gz" \
  --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
  --filter-name "LowQualitySNP" \
  -O "${wd}raw.snp.filtered.vcf.gz"


picard CreateSequenceDictionary \
      R="$genomedir$Genome" \
      O="$genomedir$Genome"

gatk SelectVariants -R "$genomedir$Genome" -V "${wd}raw.snp.filtered.vcf.gz" --exclude-filtered -O "${wd}BasePoyster.raw.snp.pass.vcf.gz"


# Step 1: Apply Variant-Level Filters
vcftools --gzvcf "${wd}BasePoyster.raw.snp.pass.vcf.gz" \
    --max-missing 0.75 \
    --min-alleles 2 \
    --max-alleles 2 \
    --minQ 99 \
    --maf 0.05 \
    --recode \
    --recode-INFO-all \
    --out "${wd}BasePoyster.snp.filtered"

# Step 2: Calculate Per-Individual Missingness
vcftools --vcf "${wd}BasePoyster.snp.filtered.recode.vcf" \
    --missing-indv \
    --out "${wd}BasePoyster.snp.filtered"

# Step 3: Identify Individuals with >20% Missing Data
awk '$5 > 0.20 { print $1 }' "${wd}BasePoyster.snp.filtered.imiss" > "${wd}low_completeness.indv"

# Optional: Check if any individuals need to be removed
if [ -s "${wd}low_completeness.indv" ]; then
    echo "Removing individuals with less than 80% genotyping completeness:"
    cat "${wd}low_completeness.indv"
else
    echo "No individuals with less than 80% genotyping completeness found."
fi

# Step 4: Remove Identified Individuals
vcftools --vcf "${wd}BasePoyster.snp.filtered.recode.vcf" \
    --remove "${wd}low_completeness.indv" \
    --recode \
    --recode-INFO-all \
    --out "${wd}BasePoyster.snp.final"


# Step 5: Re-filter
vcftools --vcf "${wd}BasePoyster.snp.final.recode.vcf" \
    --max-missing 0.75 \
    --min-alleles 2 \
    --max-alleles 2 \
    --recode \
    --recode-INFO-all \
    --out "${wd}BasePoyster.snp.final.refiltered"

    

