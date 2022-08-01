#!/bin/bash
#SBATCH --mail-user=roule@upenn.edu
#SBATCH --mail-type=ALL





# hard-code directories
fastq_raw_dir="fastq/raw"
fastq_trim_dir="fastq/trimmed_4files"
fastqc_trim_dir="data/FASTQC/trimmed_4files"

# create output directories
mkdir -p ${fastq_trim_dir}
mkdir -p ${fastqc_trim_dir}

trimmomatic_install_dir="../../Software/"

# run Trimmomatic and remove TruSeq adapters
adapter_input_list=("NAC9_myc_JKK_2018_Rep1_" "NAC10_myc_JKK_2018_Rep1_")
for x in "${adapter_input_list[@]}"; do
    raw_fastq1="${fastq_raw_dir}/${x}1.fastq.gz"
    raw_fastq2="${fastq_raw_dir}/${x}2.fastq.gz"
    trim_fastq1="${fastq_trim_dir}/${x}1.trimmed_paired.fastq.gz"
    trim_fastq2="${fastq_trim_dir}/${x}1.trimmed_unpaired.fastq.gz"
    trim_fastq3="${fastq_trim_dir}/${x}2.trimmed_paired.fastq.gz"
    trim_fastq4="${fastq_trim_dir}/${x}2.trimmed_unpaired.fastq.gz"
    java -jar ${trimmomatic_install_dir}/Trimmomatic-0.39/Trimmomatic-0.39/trimmomatic-0.39.jar \
        PE -threads 3 -phred33 ${raw_fastq1} ${raw_fastq2} \
        ${trim_fastq1} ${trim_fastq2} ${trim_fastq3} ${trim_fastq4} \
        ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 TOPHRED33

    ../../Software/fastqc_v0.11.9/FastQC/fastqc -o ${fastqc_trim_dir} ${trim_fastq1}
    ../../Software/fastqc_v0.11.9/FastQC/fastqc -o ${fastqc_trim_dir} ${trim_fastq2}
    ../../Software/fastqc_v0.11.9/FastQC/fastqc -o ${fastqc_trim_dir} ${trim_fastq3}
    ../../Software/fastqc_v0.11.9/FastQC/fastqc -o ${fastqc_trim_dir} ${trim_fastq4}
done