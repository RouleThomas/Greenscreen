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
single_rep=("NAC5_myc_JKK_2018_Rep1_" 
"NAC6_myc_JKK_2018_Rep1_"
"NAC9_myc_JKK_2018_Rep1_"
"NAC10_myc_JKK_2018_Rep1_"
"NAC_myc_JKK_2018_Input_Rep1_"
"inputA_Rep1_"
"inputC_Rep1_")
for samp in "${single_rep[@]}"; do
  cp ${in_dir}/${samp}.dupmark.sorted.bam \
    ${out_dir}/${samp}.dupmark.sorted.bam
  cp ${in_dir}/${samp}.dupmark.sorted.bam.bai \
    ${out_dir}/${samp}.dupmark.sorted.bam.bai
done



# down-sample given three reps
pool_three=("inputFGH"
"NAC127_JY_2021"
"NAC129_JY_2021"
"SNAC1_LX_2019")
for samp in "${pool_three[@]}"; do
  depth1=`samtools view -c \
    ${in_dir}/${samp}_Rep1_.dupmark.sorted.bam`
  depth2=`samtools view -c \
    ${in_dir}/${samp}_Rep2_.dupmark.sorted.bam`
  depth3=`samtools view -c \
    ${in_dir}/${samp}_Rep3_.dupmark.sorted.bam`
  arrName=(${depth1} ${depth2} ${depth3})
  minIndex "${arrName[@]}"
  let "min_idx = $min_idx + 1"
  for (( rep=1; rep<=3; rep++ )); do
    if [[ ${min_idx} -eq ${rep} ]]; then
        # copy this replicate which
        # has the smallest read depth
        cp ${in_dir}/${samp}_Rep${rep}_.dupmark.sorted.bam \
          ${out_dir}/${samp}_Rep${rep}_.dupmark.sorted.bam
        cp ${in_dir}/${samp}_Rep${rep}_.dupmark.sorted.bam.bai \
          ${out_dir}/${samp}_Rep${rep}_.dupmark.sorted.bam.bai
    else
        # downsample these replicates
        java -jar ${JVARKIT_PATH}/dist/biostar145820.jar \
            --seed ${seed} -n ${min_val} \
            -o ${out_dir}/${samp}_Rep${rep}_.dupmark.bam \
            ${in_dir}/${samp}_Rep${rep}_.dupmark.sorted.bam
	samtools sort \
            -o ${out_dir}/${samp}_Rep${rep}_.dupmark.sorted.bam \
            ${out_dir}/${samp}_Rep${rep}_.dupmark.bam && \
            rm ${out_dir}/${samp}_Rep${rep}_.dupmark.bam
	samtools index ${out_dir}/${samp}_Rep${rep}_.dupmark.sorted.bam
    fi
  done
done

