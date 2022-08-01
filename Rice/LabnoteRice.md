# Generation of a Greenscreen for rice
Analyses of 20 inputs collected from random ChIP datasets which used *O.sativa* cv. Niponbare 
1. Import fastq of the 20 inputs
Modified and launch scripts:
```
sbatch scripts/import_raw_fasta_inputs.sh #Download Fastq
sbatch scripts/organize_raw_fasta_input.sh #Compress and rename (SRA name into inputA,B, etc...;) Fastq files
sbatch scripts/fastqc_raw_inputs.sh # #Quality control to check for adaptor presences
```
Some Fastq have presence of adaptors, some does not. For these which have adaptors remove them using either TruSeq3-SE.fa or TruSeq2-SE.fa fasta files: **I did both and pick the better one**
```
sbatch scripts/trimmomatic_inputs1.sh
```
2. Download the genome and gff files and import into cluster
Genome is found [here](https://rapdb.dna.affrc.go.jp/download/irgsp1.html)\
Genome indexed with bowtie2 
```
module load bowtie2/2.4.5
bowtie2-build meta/genome/IRGSP-1.0_genome.fasta meta/genome/bowtie2_genome_dir/IRGSP
```
3. Mapping
```
sbatch scripts/mapped_inputs1.sh
```
--> I launched 5 different scripts like this one, in paralell to save time.
$FAIL, mapping rate extremely low <3%. 
$Troubleshoots:
- Maybe genome was wrong, try with the MSUv7, same result
- Maybe presence of Spike-in DNA in the data, double check each paper: not the case
- Maybe contamination? Try extract some unmapped reads and blast it:
```
samtools view -f 4 file.bam > unmapped.sam
java -jar ../../Software/picard/build/libs/picard.jar SamToFastq I= mapped/input_test/inputX_unmapped.sam FASTQ= mapped/input_test/inputX_unmapped.fastq #Transform sam into fastq
```
Blast the first 3 unmapped into NCBI, Strong hit for rice... So that is not contamination!
- Maybe trimming issues? Yes and no! Most data was paired-end and not single-end!!!
$SOLUTION, some data was paired-end, and I treated them as single-end... \
I double check all in the SRA webtools\
That fail at the early start, in paired-end there is two fastq files...\
Need to use:\
```
fasterq-dump SRR8746746 -S #add -S flag to download both R1 and R2, and not only R1....
```
Launch scripts for download all:
```
sbatch scripts/import_raw_fasta_inputsA.sh 
```
$FAIL, Downloading the data did not work using a slurm job; asking  for vdb-config and lack resource issues... As time was running I did it the old-school way, one by one in bas... Then tydying data (rename and compress) and quality control :
```
scripts/organize_raw_fasta_input_PAIRA.sh 
sbatch scripts/fastqc_raw_inputs.sh
```
Some fastq have adaptor contents, other does not.\
Dependening data paired/single-end or adpator or not I used: TruSeq3-PE.fa or TruSeq3-SE.fa paired or single, respectively; and try TruSeq2 and TruSeq3 and picked the best

***NOTE:*** Paired-end data can be treated as single end. For that, use R1 mapped reads (or R2 mapped reads), and follow Greenscreen pipeline as it was single end. I prefered to treat paired-end as paired-end (just need to adapt all future tools; trimmomatic, MACS2...
Launch mapping as follow:
```
java -jar trimmomatic-0.35.jar PE -phred33 \
    input_forward.fq.gz input_reverse.fq.gz \
    output_forward_paired.fq.gz output_forward_unpaired.fq.gz \
    output_reverse_paired.fq.gz output_reverse_unpaired.fq.gz \
    ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
```
I made several scripts to launch in parralel:
```
mapped_input1PAIRED.sh --> mapped_input5PAIRED.sh
```
Mapping rate is good (>90%)

4. Check quality of mapping
Obtain genome and chr size:
```
sbatch scripts/measureContigLengthFromFasta.sh meta/genome/IRGSP-1.0_genome.fasta > meta/genome/IRGSP-1.0_chr_count.txt
```
$Troubleshoot, same issue as reported in Tutorial; output goes into the slurm file instead of the specified output; again, no time to have a look at it so just copy/paste into the output\
Input file modified: *meta/noMaskReads_Inputs_sampleSheet.csv*\
Launch ChIPQC script **in CondaGS environment**
```
conda activate CondaGS
/home/roule/R/R-4.2.0/bin/Rscript scripts/ChIPQC_forR42.R --indivReports -g IRGSP -c chr01 chr02 chr03 chr04 chr05 chr06 chr07 chr08 chr09 chr10 chr11 chr12 -a meta/genome/IRGSP-1.0_representative/transcripts_exon.gff -s meta/genome/IRGSP-1.0_chr_count.txt meta/noMaskReads_Inputs_sampleSheet.csv data/ChIPQCreport/20inputs_noMask
```
$FAIL, input B failed; 
> names' attribute [9] must be the same length as the vector [7]

$Troubleshoots; seems to be a common error that Sheng corrected [here](https://github.com/shengqh/ChIPQC). I installed the shengqh CHIPQC script instead of the one from the tutorial; and add "force = TRUE" to the installation, within the script.\

***NOTE:*** ChIPQC only work with Single-end file. It is usefull to predict the expected fragment size notably, needed for the Greenscreen pipeline. So for the Paired-end files, cannot CHIPQC but the macs2 can work with paired-end, it will ignore the estimated fragment length (via strand cross correlation) and use the fragment length of each mated pair.\
Run ChIPQC:
```
/home/roule/R/R-4.2.0/bin/Rscript scripts/ChIPQC_forR42.R --indivReports -g IRGSP -c chr01 chr02 chr03 chr04 chr05 chr06 chr07 chr08 chr09 chr10 chr11 chr12 -a meta/genome/IRGSP-1.0_representative/transcripts_exon.gff -s meta/genome/IRGSP-1.0_chr_count.txt meta/noMaskReads_Inputs_sampleSheet.csv data/ChIPQCreport/20inputs_noMask
```
5. Call peaks with MACS2
Estimate the mappable size of the genome with faCount:
```
rsync -aP rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/linux.x86_64/faCount ./
/home/roule/GreenScreen/Software/faCount  #Comand to use to run it
../rice/GreenscreenProject/meta/genome/IRGSP-1.0_genome.fasta
```
Total size= 373245519; N= 116654\
The effective genome size would then be the total number of base pairs minus the total number of ‘N’ = 373128865\
- For single-end:
```
    macs2 callpeak \
        -t mapped/input/input${x}.dupmark.sorted.bam \
        -f BAM --keep-dup auto --nomodel \
        --extsize ${readsize} --broad --nolambda \
        -g 101274395 -n input${x} \
        --outdir ${macs_out}

done < meta/input_readsizes.csv
```
- For paired-end:
```
while read line; do
    x=$(echo $line | cut -d "," -f1)

    # call peaks with MACS2 PAIRED data
 
    macs2 callpeak \
        -t mapped/input/input${x}.dupmark.sorted.bam \
        -f BAMPE --keep-dup auto --nomodel \
        --broad --nolambda \
        -g 373128865 -n input${x} \
        --outdir ${macs_out}

Done
```
**Can do some quality control on peak calling file; cross correlation plot, seq depth and NRF**\
- cross correlation plot (only work for single-end data) with macs2-predict
```
macs2 predictd -i mapped/input/inputW.dupmark.sorted.bam -g 373128865 --outdir mapped/input/
/home/roule/R/R-4.2.0/bin/Rscript mapped/input/predictd #After running previous comand you obtain an Rscript that you just have to run to obtain a plot
```
- Seq depth
```
samtools depth inputA.sam | awk '{c++;s+=$3}END{print s/c}'
```
- To calculate the Non-Redundant Fraction (NRF) – Number of distinct uniquely mapping reads (i.e. after removing duplicates) / Total number of reads.

6. Generate the greenscreen
```
sbatch scripts/generate_20input_greenscreenBed.sh 10 1000 10
```
Here we merge within 1000bp and keep peak found in 50% of inputs (i.e. 10 inputs).\
This has been modified, generated merge parameter from 500 to 50kb

7. Obtain some Greenscreen genome parameter/characteristics
- Percent of genome covered
Install bedops:
```
git clone https://github.com/bedops/bedops.git
cd bedops
make
make install
```
Then launch command as folow:
```
../../Software/bedops/bin/bedmap --echo --bases-uniq --delim '\t' data/macs2_out/inputControls/qval10/chr_size.bed data/macs2_out/inputControls/qval10/gs_merge1000bp_call10_20inputs.bed | awk 'BEGIN { genome_length = 0; masked_length = 0; } { genome_length += ($3 - $2); masked_length += $4; } END { print (masked_length / genome_length); }'
```
--> Multiply per 100 to get result in %.\
- Number of bp covered by Greenscreen mask
```
../../Software/bedops/bin/bedmap --echo --bases-uniq --delim '\t' data/macs2_out/inputControls/qval10/chr_size.bed data/macs2_out/inputControls/qval10/gs_merge10000bp_call10_20inputs.bed | awk 'BEGIN { genome_length = 0; masked_length = 0; } { genome_length += ($3 - $2); masked_length += $4; } END { print (masked_length); }'
```
- Number of genes covered
```
bedtools intersect -a data/macs2_out/inputControls/qval10/gs_merge1000bp_call10_20inputs.bed -b meta/genome/IRGSP-1.0_representative/locus.gff -c | awk '{sum += $11} END {print sum}'
```

8. Analyzed ChIPseq files and apply the generated Greenscreen to assess its efficiency\
8.1 Download the data one-by-one as script not working ($FAIL, no time to troubleshoot) and tidy/rename/compress it\
```
fasterq-dump XXX -S
sbatch scripts/fastqc_raw_publishedChIPsAndControls.sh
```
8.2 Check quality and trimmed accordingly (dependending on adaptor presence)\
If adaptors present I tested both TruSeq3 and TruSeq2 fasta sequences and keep the best one (for either PE and SE)\
Launch scripts:
```
trimmomatic_publishedChIPsAndControls1.sh #same until trimmomatic_publishedChIPsAndControls6.sh
```
8.3 Then mapped the trimmed reads:
```
mapped_IP1PAIRED.sh #Same until mapped_IP7PAIRED.sh
```
8.4 Downsamples when multiple replicates
```
downSampleBam_publishedChIPsAndControlsSINGLEREP.sh #No need to downsample; here just copy file to new "downsample" folder
downSampleBam_publishedChIPsAndControlsTRIPLEREP.sh
```
8.5 Call peaks with MACS2
The script were super buggy, so because I lack of time, I run it old-school in the comand line, as follow:
```
macs2 callpeak -t mapped/chip/NAC5_myc_JKK_2018_R1.dupmark.sorted.bam \
        	                -c mapped/chip/downsample/NAC_input_JKK_2018_R1.dupmark.sorted.bam \
                	        -f BAMPE --keep-dup auto \
                        	--nomodel -g 373128865 \
	                        --outdir data/macs2_out/chipPeaksPAIRED -n NAC5_myc_JKK_2018_R1

## very long
awk -F"\t" -v q=2 'BEGIN{OFS="\t"} \
                	data/macs2_out/chipPeaksPAIRED/NAC5_myc_JKK_2018_R1_peaks.narrowPeak > \
                       	data/macs2_out/chipPeaksPAIRED/noMask_qval_corrNACgood2/NAC5_myc_JKK_2018_R1_peaks.narrowPeak



bedtools intersect -v -wa \
                        -a data/macs2_out/chipPeaksPAIRED/noMask_qval_corr2/NAC5_myc_JKK_2018_R1_peaks.narrowPeak \
                        -b meta/gs_merge2000bp_call10_20inputs_corr.bed > \
                        data/macs2_out/chipPeaksPAIRED/gsMask_qval_corr2/NAC5_myc_JKK_2018_R1_peaks.narrowPeak

```
8.6 Generate the bigwig coverage files
```
sbatch scripts/bamToBigWig_publishedChIPsAndControls_replicates1.sh 
sbatch scripts/bamToBigWig_publishedChIPsAndControls_replicates1.sh
```
8.7 Called peaks taking into account their respective controls
```
scripts/macs2_callpeaks_publishedChIPsPAIREDpval5_corr.sh
```
***NOTE:*** For the paired-end data I need to indicate median fragment size  in meta/chip_controls_fragsize_nreps.csv; for that:
```
conda activate CondaGS
cd Software
pip install deeptools
bamPEFragmentSize -b mapped/chip/downsample/SNAC1_LX_2019_Rep1_.dupmark.sorted.bam #Output will give read and fragment median sizes
```
Then merge all the peaks into one bed file using:
```
sbatch merge_chip_peaks_p5_GS.sh
```
Generate a matrix taking all peak and coverage file into account :
```
conda activate CondaGS
pip install pandas # Need to be installed 
python3 scripts/readCorrelationPlot.py \
data/plotCorrelation/coverage_matrix_trueRep_peaks_merged.csv \
data/plotCorrelation/trueRep_peaks_merged_heatmap.png \
-lm ward --plot_numbers -k 2 -ri \
-sl meta/chip_trueReps_colorshapeLabels.csv \
-cf meta/chip_trueReps_expectedCluster.csv
```
Then, generate a heatmap from the matrix file:
```
pip install seaborn  # Need to be installed 
python3 scripts/readCorrelationPlot.py \
data/plotCorrelation/coverage_matrix_trueRep_peaks_merged.csv \
data/plotCorrelation/trueRep_peaks_merged_heatmap.png \
-lm ward --plot_numbers -k 2 -ri \
-sl meta/chip_trueReps_colorshapeLabels.csv \
-cf meta/chip_trueReps_expectedCluster.csv
```
Repeat the following scripts, changing on the top header qvalue parameters or Greenscreen to use
```
sbatch scripts/macs2_callpeaks_publishedChIPsPAIREDpval5_corr.sh
sbatch merge_chip_peaks_p5_GS.sh
python3 scripts/coverage_bed_matrix.py \
meta/chip_trueRep_bigwigs.csv \
data/macs2_out/chipPeaksPAIRED/gsMask_qval_corr5/ChIPseq_Peaks_gsqval5.merged.bed \
-o data/plotCorrelation \
-m coverage_matrix_trueRep_peaks_merged.csv

python3 scripts/readCorrelationPlot.py \
data/plotCorrelation/coverage_matrix_trueRep_peaks_merged.csv \
data/plotCorrelation/trueRep_peaks_merged_heatmap.png \
-lm ward --plot_numbers -k 2 -ri \
-sl meta/chip_trueReps_colorshapeLabels.csv \
-cf meta/chip_trueReps_expectedCluster.csv
```
