#!/bin/bash

wd="/home/kangziyi/poster/callSNP/"
Genome="GCF_963853765.1_xbMagGiga1.1_genomic.fna"
genomedir="/home/kangziyi/poster/genome/"
sradir="/home/kangziyi/poster/seqdata/savedseq/"
outfile="BasePoyster.raw.vcf.gz"
mkdir -p "$wd"

find "$sradir" -name "SRR*" | cut -d '/' -f 7 | cut -d '.' -f 1 > "${wd}samplename.txt"

if sort "${wd}samplename.txt" | uniq -d | grep -q .; then
    # If there are duplicates, terminate the script with a message
    echo "Duplicates found! Terminating."
    exit 1
else
    # If no duplicates, continue executing
    echo "No duplicates found. Proceeding with the program."
    # Place any additional commands here
fi

mkdir -p "${wd}seqfq"
mkdir -p "${wd}seqvcf"
mkdir -p "${wd}seqbam"

find "$sradir" -name "SRR*" | xargs -I {} fasterq-dump {} --outdir "${wd}seqfq"

gtx index "$genomedir$Genome"

#easy code
#each line of samplename.txt contain an unique sample name

while read -r fn; do
    echo "Processing sample $fn"
    # Define FASTQ files
    #fq1="${wd}seqfq/${fn}_1.fastq"
    #fq2="${wd}seqfq/${fn}_2.fastq"
    fq_single=(${wd}seqfq/${fn}*.fastq)
    gtx wgs \
   -t 40 \
   -g \
   -o "${wd}seqvcf/${fn}.g.vcf.gz" \
   -b "${wd}seqbam/${fn}.bam" \
   -R "@RG\tID:${fn}\tSM:${fn}" \
   $genomedir$Genome \
   $fq_single \

done < "${wd}samplename.txt"

awk -v wd="$wd" '{print $0 "\t" wd "seqvcf/" $0 ".g.vcf.gz"}' "${wd}samplename.txt" > "${wd}sample_vcf_mapping.txt"


gtx joint --sample-name-map "${wd}sample_vcf_mapping.txt" --reference "$genomedir$Genome" -t 40 -o "$wd$outfile"
#
