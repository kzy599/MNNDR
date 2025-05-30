# MNNDR

This repository supports the unpublished research paper: **"MNNDR: multi-Input neural network improves single- and multi-trait genomic prediction of disease resistance in aquaculture"**.

## Prerequisites

Before running the program, ensure the following dependencies are installed:
- TensorFlow == 2.13.1
- AlphaSimR (latest version)
- Hiblup (latest version)
- Genotype processing tools:  
  `gtx.cat` (recommended) **or** GATK + samtools (alternative)

## Workflow Overview
## Simulated datasets

### 1. Genotype Data Processing
Process raw genotype sequence data through:
1. `callSNP.sh` - SNP calling
2. `snpQC.sh` - Quality control
3. `phased.sh` - Haplotype phasing

### 2. SNP Panel Creation and Haplotype extraction for simulation
```bash
./createMapHaplo.sh
```


### 3. Dataset Generation
```R
Rscript dataGenerator.r
```
Note: Configure these parameters in the script before execution:

Dominance degree magnitude
Epistatic-to-additive ratio

### 4. DL Network Optimization
```python
python main_tunning.py --th
python mian_tunning.py --th --double
python main_tunning.py --th --double --linear
```
Performs hyperparameter tuning for DL.

### 5. Model Validation
```R
Rscript runR.r
```
Important: Update hyperparameters in these files first:

main.py
utils.py

### 6. Result Analysis
```R
Rscript analyzeResult.r
```
Generates performance comparisons and visualizations between DL networks and BLUP models.

## Real datasets

### 1. Unzip the `Crap1259.zip` or `Trout1935.zip`
### 2. Process the data format using `dataGenerator.r`
### 3. Perform 5-fold cross-validation using `runR.r`


### Maintainers
For questions or support, please contact:

Ziyi Kang
Ocean University of China
kangziyi1998@163.com

Qi Li
Ocean University of China
qili66@ouc.edu.cn
