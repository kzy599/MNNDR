# MT-DLnetwork
Deep Learning for Oyster Genomic Prediction
**Code for:**​ "Deep learning and multi-trait genomic prediction facilitate selective breeding for summer mortality resistance in Pacific oyster (Crassostrea gigas)"

🚀 Quick Start
Dependencies
Install these ​before running:

​TensorFlow 2.13.1:
bash
pip install tensorflow==2.13.1  
​AlphaSimR​ (R package for genomic simulation):
r
install.packages("AlphaSimR")  
​Hiblup: Follow Hiblup’s official guide.
​Genotype Tools: Use gtx.cat (or replace with ​GATK/samtools).
🔄 Pipeline Steps
1. ​Process Genotype Data​
Run in order:

​Call SNPs:
bash
./callSNP.sh  
​Quality Control:
bash
./snpQC.sh  
​Haplotype Phasing:
bash
./phased.sh  
2. ​Create SNP Panel​
Run createMapHaplo.sh to:
Build SNP panel.
Select QTLs and extract haplotype data.
3. ​Generate Dataset​
Adjust parameters in dataGenerator.r (dominance degree, epistatic ratio).
Run to create simulation data.
4. ​Tune DL Hyperparameters​
Optimize model settings:
bash
python main_tunning.py  
5. ​Cross-Validation​
Update main.py and utils.py with tuned parameters.
Run validation for DL and BLUP models:
bash
Rscript runR.r  
6. ​Analyze Results​
Plot metrics and compare models:
bash
Rscript analyzeResult.r  
📁 Directory Structure
├── genotype_processing/  # SNP calling, QC, phasing  
├── snp_panel/            # SNP/QTL selection  
├── data_generation/      # Simulate datasets  
├── dl_tuning/            # Hyperparameter optimization  
├── cross_validation/     # Model validation  
├── results_analysis/     # Visualizations  
└── utils/                # Helper functions  
💡 Notes
​Replace gtx.cat: Use GATK/samtools by editing callSNP.sh and phased.sh.
​GPU Support: Install CUDA/cuDNN for TensorFlow GPU acceleration.
​Citation: Cite this paper once published.
📬 Contact
For questions or collaborations:

​Ziyi Kang​
kangziyi1998@163.com
​Qi Li​
qili66@ouc.edu.cn
