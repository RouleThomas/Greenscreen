#!/bin/bash
#SBATCH --mail-user=roule@upenn.edu
#SBATCH --mail-type=ALL


raw_bam_dir="mapped/chip"
downsamp_bam_dir="${raw_bam_dir}/downsample"
bam_suffix="dupmark.sorted.bam"
macs2_out="data/macs2_out/chipPeaksPAIRED"
gs_regions="meta/gs_merge2000bp_call10_20inputs_corr.bed"
# average basepair q-value threshold (log5)
q=2



# make macs2 output directory
mkdir -p ${macs2_out}/noMask_qval_corr${q}
mkdir -p ${macs2_out}/gsMask_qval_corr${q}

while read line; do
        chip_name=$(echo $line | cut -d "," -f1)
        chip_nreps=$(echo $line | cut -d "," -f2)
        cntl_name=$(echo $line | cut -d "," -f3)
        cntl_nreps=$(echo $line | cut -d "," -f4)

        c_param="${downsamp_bam_dir}/${cntl_name}_R1.${bam_suffix}"
        for ((c=2; c<=${cntl_nreps}; c++ )); do
                c_param="${c_param} ${downsamp_bam_dir}/${cntl_name}_R${c}.${bam_suffix}"
        done

        # call peaks on individual replicates

        # run MACS2
        for ((t=1; t<=${chip_nreps}; t++ )); do
                t_raw_rep="${raw_bam_dir}/${chip_name}_R${t}.${bam_suffix}"
                if [[ ! -f "${macs2_out}/${chip_name}_R${t}_peaks.narrowPeak" ]];then
			macs2 callpeak -t ${t_raw_rep} \
        	                -c ${c_param} \
                	        -f BAMPE --keep-dup auto \
                        	--nomodel -g 373128865 \
	                        --outdir ${macs2_out} -n ${chip_name}_R${t}
		fi
                # remove all peaks that do not have an
                # average base pair q-value <=10^(-${q})
		awk -F"\t" -v q=${q} 'BEGIN{OFS="\t"} \
        		$9>=q && $1!="ChrC" && $1!="ChrM"{print}' \
                	${macs2_out}/${chip_name}_R${t}_peaks.narrowPeak > \
                       	${macs2_out}/noMask_qval_corr${q}/${chip_name}_R${t}_peaks.narrowPeak

                # remove all peaks that overlap greenscreen
                bedtools intersect -v -wa \
                        -a ${macs2_out}/noMask_qval_corr${q}/${chip_name}_R${t}_peaks.narrowPeak \
                        -b ${gs_regions} > \
                        ${macs2_out}/gsMask_qval_corr${q}/${chip_name}_R${t}_peaks.narrowPeak
        done


        # get pooled chip bams
        t_down_rep="${downsamp_bam_dir}/${chip_name}_R1.${bam_suffix}"
        t_pool_param="${t_down_rep}"
        for ((t=2; t<=${chip_nreps}; t++ )); do
                t_down_rep="${downsamp_bam_dir}/${chip_name}_R${t}.${bam_suffix}"
                t_pool_param="${t_pool_param} ${t_down_rep}"
        done
        echo "-t ${t_pool_param}"

        # call peaks on pooled ${chip_name}

        # run MACS2
        if [[ ! -f "${macs2_out}/${chip_name}_peaks.narrowPeak" ]];then
        	macs2 callpeak -t ${t_pool_param} \
                	-c ${c_param} \
	                -f BAMPE --keep-dup auto \
        	        --nomodel -g 373128865 \
                	--outdir ${macs2_out} -n ${chip_name}
	fi
        # remove all peaks that do not have an
        # average base pair q-value <=10^(-${q})
        awk -F"\t" -v q=${q} 'BEGIN{OFS="\t"} \
                $9>=q && $1!="ChrC" && $1!="ChrM"{print}' \
                ${macs2_out}/${chip_name}_peaks.narrowPeak > \
                ${macs2_out}/noMask_qval_corr${q}/${chip_name}_peaks.narrowPeak

        # remove all peaks that overlap greenscreen
        bedtools intersect -v -wa \
                -a ${macs2_out}/noMask_qval_corr${q}/${chip_name}_peaks.narrowPeak \
                -b ${gs_regions} > \
                ${macs2_out}/gsMask_qval_corr${q}/${chip_name}_peaks.narrowPeak

done < meta/chip_controls_fragsize_nreps.csv
