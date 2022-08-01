#!/bin/bash

in_dir="mapped/chip"
out_dir="mapped/chip/downsample"
JVARKIT_PATH="/home/roule/GreenScreen/Software/jvarkit"
seed=42
mkdir -p ${out_dir}

module load samtools/1.15

# function to find the minimum
# value in an array
minIndex(){
   arr=("$@")
   min_val=${arr[0]}
   min_idx=0
   for i in ${!arr[@]}; do
        cur_val=${arr[${i}]}
        if [[ ${cur_val} -lt ${min_val} ]]; then
                min_val=${arr[$i]}
                min_idx=${i}
        fi
   done

}


# There is no need to down-sample anything
# with a single replicate
single_rep=("NAC_input_JKK_2018_Rep1_")
for samp in "${single_rep[@]}"; do
  cp ${in_dir}/${samp}.dupmark.sorted.bam \
    ${out_dir}/${samp}.dupmark.sorted.bam
  cp ${in_dir}/${samp}.dupmark.sorted.bam.bai \
    ${out_dir}/${samp}.dupmark.sorted.bam.bai
done

