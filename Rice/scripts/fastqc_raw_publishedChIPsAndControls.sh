#!/bin/bash
#SBATCH --mail-user=roule@upenn.edu
#SBATCH --mail-type=ALL

# create output directory
mkdir -p data/FASTQC/raw

input_list=("NFYC11_JY_2019_Rep1_1" "NFYC11_JY_2019_Rep1_2")
        
for x in "${input_list[@]}"; do
    ../../Software/fastqc_v0.11.9/FastQC/fastqc -o data/FASTQC/raw fastq/raw/${x}.fastq.gz
done
