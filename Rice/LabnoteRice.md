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
Dependening data paired/single-end or adpator or not I used: TruSeq3-PE.fa or TruSeq3-SE.fa paired or single, respectively; and try TruSeq2 and TruSeq3 and pipck the best














