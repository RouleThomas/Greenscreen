#!/bin/bash


module load samtools/1.15

genome="meta/genome/IRGSP-1.0_chr_count.txt"
indir="mapped/chip"
outdir="data/bigwigs/chip/individual_replicates"
normalizeby=10000000 # scaling factor

mkdir -p ${outdir}




# published ChIP-Seq controls


pubControl_list=("NAC_input_JKK_2018_R1")

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
