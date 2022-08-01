#!/bin/bash


module load samtools/1.15

genome="meta/genome/IRGSP-1.0_chr_count.txt"
indir="mapped/chip"
outdir="data/bigwigs/chip/individual_replicates"
normalizeby=10000000 # scaling factor

mkdir -p ${outdir}

# ChIP-Seq samples

while read line; do
    samp=$(echo $line | cut -d "," -f1)
    readSize=$(echo $line | cut -d "," -f2)
    fragSize=$(echo $line | cut -d "," -f3)
    extend=`awk -v f=${fragSize} -v r=${readSize} 'BEGIN{print f-r}'`
    orig_bam="${indir}/${samp}.dupmark.sorted.bam"
    if [[ ! -f "${outdir}/${samp}.bw" ]]; then
        # convert BAM to BED format
        bamToBed \
          -i ${orig_bam} \
          > ${outdir}/${samp}.bed

        # sort the BED file
        sort -k 1,1 \
          ${outdir}/${samp}.bed > \
          ${outdir}/${samp}.sorted.bed

        # extend the reads
        slopBed -i ${outdir}/${samp}.sorted.bed \
            -l 0 -r ${extend} -s -g ${genome} \
            > ${outdir}/${samp}.extend.bed

        # normalize signal and output BEDGRAPH
        totreads=`samtools view -c ${orig_bam}`
        scaling=`awk 'BEGIN{ print '"$normalizeby"' / '"$totreads"'}'`

        genomeCoverageBed -i \
           ${outdir}/${samp}.extend.bed \
           -g ${genome} -bg -scale $scaling | \
           awk 'BEGIN{OFS=FS="\t"} \
           {$4=sprintf("%.2f",$4)}{print}' \
           > ${outdir}/${samp}.bg

        # compress BEDGRAPH to BIGWIG FORMAT
        /home/roule/GreenScreen/Software/bedGraphToBigWig ${outdir}/${samp}.bg \
           ${genome} ${outdir}/${samp}.bw
    fi
done < meta/chip_readsize_fragsize.csv

# published ChIP-Seq controls


pubControl_list=("inputA_Rep1_" "inputC_Rep1_"
  "inputFGH_Rep1_" "inputFGH_Rep2_"
  "inputFGH_Rep3_" "NAC_myc_JKK_2018_Input_Rep1_")

for samp in "${pubControl_list[@]}"; do
    if [[ ! -f "${outdir}/${samp}.bw" ]]; then
        orig_bam="${indir}/${samp}.dupmark.sorted.bam"

        # normalize signal and output BEDGRAPH
        totreads=`samtools view -c ${orig_bam}`
        scaling=`awk 'BEGIN{ print '"$normalizeby"' / '"$totreads"'}'`

        genomeCoverageBed -ibam ${orig_bam} \
          -bg -scale $scaling | \
          awk 'BEGIN{OFS=FS="\t"} \
          {$4=sprintf("%.2f",$4)}{print}' \
          > ${outdir}/${samp}.bg

        # compress BEDGRAPH to BIGWIG FORMAT
        /home/roule/GreenScreen/Software/bedGraphToBigWig ${outdir}/${samp}.bg \
          ${genome} ${outdir}/${samp}.bw \
          && rm ${outdir}/${samp}.bg
    fi
done
