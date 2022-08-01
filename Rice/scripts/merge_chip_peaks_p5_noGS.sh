#!/bin/bash

chip_peaks_dir="data/macs2_out/chipPeaksPAIRED/noMask_qval_corrNACgood2"




cat \
	${chip_peaks_dir}/NAC9_myc_JKK_2018_R1_peaks.broadPeak \
	${chip_peaks_dir}/NAC10_myc_JKK_2018_R1_peaks.broadPeak \
	${chip_peaks_dir}/NAC127_JY_2021_R1_peaks.broadPeak \
	${chip_peaks_dir}/NAC127_JY_2021_R2_peaks.broadPeak \
	${chip_peaks_dir}/NAC127_JY_2021_R3_peaks.broadPeak \
	${chip_peaks_dir}/NAC129_JY_2021_R1_peaks.broadPeak \
	${chip_peaks_dir}/NAC129_JY_2021_R2_peaks.broadPeak \
	${chip_peaks_dir}/NAC129_JY_2021_R3_peaks.broadPeak \
	${chip_peaks_dir}/NAC5_myc_JKK_2018_R1_peaks.broadPeak \
	${chip_peaks_dir}/NAC6_myc_JKK_2018_R1_peaks.broadPeak \
	${chip_peaks_dir}/NFYC11_JY_2019_R1_peaks.broadPeak \
	${chip_peaks_dir}/SNAC1_LX_2019_R1_peaks.broadPeak \
	${chip_peaks_dir}/SNAC1_LX_2019_R2_peaks.broadPeak \
	${chip_peaks_dir}/SNAC1_LX_2019_R3_peaks.broadPeak | \
	sort -k1,1 -k2,2n | \
	bedtools merge -i - > \
	${chip_peaks_dir}/ChIPseq_Peaks.merged.bed
