#!/bin/bash

Rscript dataGenerator.r

source activate tf-gpu


python main_tunning.py --Th

python main_tunning.py --Th --doubleT

#python main_tunning.py --Th --doubleT --softm

python main_tunning.py --Th --doubleT --linear


#python main.py --Th

#python main.py --Th --doubleT

#python main.py --Th --doubleT --softm