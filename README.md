# MT-DLnetwork
# Deep Learning-Based Genomic Prediction for Oyster Summer Mortality Resistance

This repository supports the unpublished research paper: **"Deep Learning and Multi-Trait Genomic Prediction Facilitate Selective Breeding for Summer Mortality Resistance in Pacific Oyster (Crassostrea gigas)"**.

## Prerequisites

Before running the program, ensure the following dependencies are installed:
- TensorFlow == 2.13.1
- AlphaSimR (latest version)
- Hiblup (latest version)
- Genotype processing tools:  
  `gtx.cat` (recommended) **or** GATK + samtools (alternative)

## Workflow Overview

### 1. Genotype Data Processing
Process raw genotype sequence data through:
1. `callSNP.sh` - SNP calling
2. `snpQC.sh` - Quality control
3. `phased.sh` - Haplotype phasing

### 2. SNP Panel Creation
```bash
createMapHaplo.sh
、、、

### 2. SNP Panel Creation
3. Dataset Generation
Rscript dataGenerator.r

Note: Configure these parameters in the script before execution:

Dominance degree magnitude
Epistatic-to-additive ratio
4. DL Network Optimization
main_tunning.py

Performs hyperparameter tuning for deep learning architectures.

5. Model Validation
Rscript runR.r

Important: Update hyperparameters in these files first:

main.py
utils.py
6. Result Analysis
Rscript analyzeResult.r

Generates performance comparisons and visualizations between DL networks and BLUP models.

Maintainers
For questions or support, please contact:

Ziyi Kang
Ocean University of China
kangziyi1998@163.com

Qi Li
Ocean University of China
qili66@ouc.edu.cn
