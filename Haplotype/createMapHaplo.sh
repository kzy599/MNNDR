#!/bin/bash

CPM="ChPoMa.csv"
nQTL=1000
panelsize=20000
vcfpath="/home/kangziyi/poster/callSNP/BasePoyster.snp.final.phased.vcf"
iodir="/home/kangziyi/poster/BaseMapHaplo"
panelMap="panelMap.csv"
qtlMap="qtlMap.csv"
mergedMap="mergedMap.csv"
panelHaplo="panelHaplo.csv"
qtlHaplo="qtlHaplo.csv"
mergedHaplo="mergedHaplo.csv"


python3 CreateCPM.py -i $vcfpath -o $iodir -f $CPM -n $nQTL

python3 CreatePanel.py -io $iodir -c $CPM -q $qtlMap -p $panelMap -m $mergedMap -s $panelsize

python3 ExtractHaplo.py -i $vcfpath -io $iodir -q $qtlMap -p $panelMap -m $mergedMap -fq $qtlHaplo -fp $panelHaplo -fm $mergedHaplo