#!/bin/bash
#SBATCH --mail-user=roule@upenn.edu
#SBATCH --mail-type=ALL

nthreads=3
PICARD_PATH="../../Software/picard/build/libs"
module load samtools/1.15
module load bowtie2/2.4.5


# create output directory
mkdir -p mapped/input

input_list=("U" "V" "X" "A")

for x in "${input_list[@]}"; do
    # run bowtie2
    bowtie2  --phred33 -q \
	-x meta/genome/bowtie2_genome_dir/IRGSP \
        -S mapped/input/input${x}.sam \
        -1 fastq/trimmed_4files/input${x}1.trimmed_paired.fastq.gz \
        -2 fastq/trimmed_4files/input${x}2.trimmed_paired.fastq.gz
    # sort the reads
    samtools sort -o mapped/input/input${x}.bam \
        mapped/input/input${x}.sam
    # index the bam file
    samtools index mapped/input/input${x}.bam
    # remove reads without MAPQ>=30
    samtools view -@ ${nthreads} -F 772 -q 30 \
        -b mapped/input/input${x}.bam \
	chr01 chr02 chr03 chr04 chr05 chr06 chr07 chr08 chr09 chr10 chr11 chr12 | \
        samtools sort - -o mapped/input/input${x}.filter.bam
    # index filtered reads
    samtools index mapped/input/input${x}.filter.bam
    # mark duplicates with picard
    java -jar ${PICARD_PATH}/picard.jar MarkDuplicates \
        -I mapped/input/input${x}.filter.bam \
        -O mapped/input/input${x}.dupmark.bam \
        -M mapped/input/input${x}.dup.qc \
        -VALIDATION_STRINGENCY LENIENT \
        -REMOVE_DUPLICATES false \
	-ASSUME_SORTED true
    # sort reads after marking the duplicates
    samtools sort -o mapped/input/input${x}.dupmark.sorted.bam \
        mapped/input/input${x}.dupmark.bam
    # index the sorted reads
    samtools index mapped/input/input${x}.dupmark.sorted.bam
done
