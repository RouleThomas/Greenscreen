# Greenscreen tutorial in Upenn HPC cluster
## Prerequisites
Klasfeld et al 2022 tutorial can be found [here](https://github.com/sklasfeld/GreenscreenProject.git).
Only 3 inputs will be analyzed to make the tutorial faster.\
**Frequently used command line:** 
```
srun --nodelist=node01 --mem=20g --pty bash -l #use node 1 or node3 otherwise few tools are not installed
squeue -u roule #check status of your jobs
sbatch FILENAME.sh #launch a job
```
**Additional tools to download/install**\
Create *Software* folder regrouping all aditional softwares
```
mkdir Software
cd Software
```
1. fastqc
Download FASTQC binary files, unzip and transfer to *software* folder using Mobaxterm
```
chmod +x PATH/FastQC/fastqc
fastqc_v0.11.9/FastQC/fastqc # comand to use to launch FASTQC
```
2. Picard
Clone Picard repo/
```
git clone https://github.com/broadinstitute/picard.git
cd picard/
./gradlew shadowJar
java -jar Software/picard/build/libs/picard.jar # Comand to use to launch Picard
```
## Get started
### Generation of the Greenscreen using 3 inputs
1. Download the 3 inputs data
Modification of the *import_raw_fasta_inputs.sh* script. 
2. Clean input reads (remove adapters and perform QC)
Modification of *trimmomatic_inputs.sh* script. $time, around 1h per sample of 1Go
3. Mapping
First, genome indexation with bowtie2
```
bowtie2-build meta/ArabidopsisGenome/TAIR10_Chr.all.fasta \
meta/ArabidopsisGenome/bowtie2_genome_dir/TAIR10
```
Then, start mapping with modification of *mapped_inputs.sh*/
Check mapped files do not show a strand bias (characterisit of input samples)/
Modify the meta file *noMaskReads_Inputs_sampleSheet.csv* to look only at the 3 mapped inputs/
Change gff format from unicode to Unix format
```
chmod +x scripts/translateUniCodeFile.py
python scripts/translateUniCodeFile.py meta/ArabidopsisGenome/Araport11_GFF3_genes_transposons.201606.gff meta/ArabidopsisGenome/Araport11_GFF3_genes_transposons.UPDATED.201606.gff
```
$FAIL, maybe because Python2 is on my system and need Python3?\
NOPE, $FAIL again; \
error message: \
>  from urllib.parse import unquote ImportError: No module named parse

$SOLUTION
Change the top lines\
```
import sys
reload(sys)
sys.setdefaultencoding('utf-8')
import argparse
from urlparse import unquote
from w3lib.html import replace_entities
```
Instead of:
```
import sys
import argparse
from urllib.parse import unquote
from w3lib.html import replace_entities
```
Run the script with:
```
python scripts/translateUniCodeFile.py meta/ArabidopsisGenome/Araport11_GFF3_genes_transposons.201606.gff meta/ArabidopsisGenome/Araport11_GFF3_genes_transposons.UPDATED.201606.gff
```
4. Obtain genome and Chr sizes
```
sbatch scripts/measureContigLengthFromFasta.sh meta/ArabidopsisGenome/TAIR10_Chr.all.fasta > meta/ArabidopsisGenome/TAIR10_chr_count.txt
```
$FAIL, The ouput goes into the slurm output, instead of the specified output "TAIR10_chr_count.txt"\
$SOLUTION, not clean... cp paste output fromm slurm into the specified output\
5. Use of ChIPQC to assess quality of ChIPseq data
$FAIL
The *ChIPQC.R* script crashed\
RSVG-related error message\
Maybe trouble installing package thourgh the script?\
$Troubleshoot, try installing package manually within R\
$FAIL, Try to change the shebang line from R script\
```
#!/usr/bin/Rscript –vanilla = ignore user-specific R settings
```
$FAIL, Try run everything from R\
$FAIL at installing rsvg package\
$Troubleshoot, try using other version of R (R3.2.3)\
$FAIL\
> In /usr/lib/pkgconfig/ I miss the librsvg-2.0.pc

$SOLUTION, need to install *librsvg2* on the cluster
```
sudo apt-get install -y librsvg2-dev
```
$FAIL, cannot as I need pasword and authorization,\
$SOLUTION, HPC Upenn team install it for me. But I could have installed it using a Conda Environment --> To use for the real analyses\
$FAIL, Bioconductor bug\
$Troubleshoot, try to force the installation
```
if( length(notInstalledBiocPackages) ) {
    for (biocPackage in listOfBiocPackages[ notInstalledBiocPackages ]) {
        BiocManager::install(biocPackage, force=TRUE)
    }
} 
```
$FAIL, try install package Rcurl\
$FAIL, try download binary files fopr Rcurl package and import it\
$SOLUTION, Install *libcurl-devel* on the system\

$FAIL with BiocManager, 
> ERROR: failed to lock directory ‘/home/roule/R/x86_64-pc-linux-gnu-library/4.0’ for modifying

$Troubleshoot, Try removing ‘/home/roule/R/x86_64-pc-linux-gnu-library/4.0/00LOCK-limma’
> ERROR: dependency ‘locfit’ is not available for package ‘DESeq2’ ERROR: dependencies ‘limma’, ‘locfit’ are not available for package ‘edgeR’removing ‘/home/roule/R/x86_64-pc-linux-gnu-library/4.0/edgeR’ ERROR: failed to lock directory ‘/home/roule/R/x86_64-pc-linux-gnu-library/4.0’ for modifying

$Troubleshoot, Try removing ‘/home/roule/R/x86_64-pc-linux-gnu-library/4.0/00LOCK-limma’
> ERROR: dependency ‘locfit’ is not available for package ‘DESeq2’ removing ‘/home/roule/R/x86_64-pc-linux-gnu-library/4.0/DESeq2’ ERROR: dependencies ‘limma’, ‘locfit’ are not available for package ‘edgeR’ removing ‘/home/roule/R/x86_64-pc-linux-gnu-library/4.0/edgeR’ ERROR: dependencies ‘limma’, ‘edgeR’, ‘DESeq2’ are not available for package ‘systemPipeR’ removing ‘/home/roule/R/x86_64-pc-linux-gnu-library/4.0/systemPipeR’ ERROR: dependencies ‘limma’, ‘locfit’, ‘systemPipeR’, ‘DESeq2’ are not available for package ‘DiffBind’ removing ‘/home/roule/R/x86_64-pc-linux-gnu-library/4.0/DiffBind’ ERROR: dependency ‘DiffBind’ is not available for package ‘ChIPQC’ removing ‘/home/roule/R/x86_64-pc-linux-gnu-library/4.0/ChIPQC’


$SOLUTION, was a version/compatible issues --> Need to use Biocmanager version 3.12 (and not 3.15) as we are in R4.0 (and not 4.4), see [here](https://www.bioconductor.org/about/release-announcements/)
```
BiocManager::install(version = "3.12", lib= "/home/roule/R/x86_64-pc-linux-gnu-library/4.0/")
```
$FAIL, at installing the ChIPQC package
```
BiocManager::install(c('ChIPQC'), lib="/home/roule/R/x86_64-pc-linux-gnu-library/4.0/")
```
> ERROR: failed to lock directory ‘/home/roule/R/x86_64-pc-linux-gnu-library/4.0’ for modifying

$Troubleshoot, goes a bit further but fail again
```
BiocManager::install(c('ChIPQC'), lib="/home/roule/R/x86_64-pc-linux-gnu-library/4.0/", INSTALL_opts = '--no-lock')
```
> ERROR: dependency ‘locfit’ is not available for package ‘DESeq2’

```
BiocManager::install(c('ChIPQC'), lib="/home/roule/R/x86_64-pc-linux-gnu-library/4.0/", INSTALL_opts = '--no-lock', dependencies=TRUE)
```
$FAIL, Try with another version of R (3.4.2) installing BioClite
```
source("https://bioconductor.org/biocLite.R")
BiocInstaller::biocLite(c("ChIPQC"))
```
$FAIL, try keeping dependencies
```
BiocInstaller::biocLite(c("ChIPQC"), lib="/home/roule/R/x86_64-pc-linux-gnu-library/3.4/", dependencies=TRUE)
```
$FAIL, try R 3.5.0
```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("ChIPQC")
```
$Troubleshoot, locfic require R version > 4.1 (see [here](http://cran.nexr.com/web/packages/locfit/index.html), but the ChIPQC.R works better with R4.0.4\
So let's try an older version of *locfit*
```
install.packages("devtools") 
install_version("locfit", version = "1.5-9.1", repos = "http://cran.us.r-project.org")
```
$SOLUTION, install older version of *locfit* with URL
```
packageurl <- "https://cran.r-project.org/src/contrib/Archive/locfit/locfit_1.5-9.tar.gz"
install.packages(packageurl, repos=NULL, type="source")
library(locfit)
BiocManager::install("ChIPQC")
```
$FAIL, bug at the end, cannot save PNG\
```
ChIPQC Rscript: Rscript scripts/ChIPQC.R --indivReports -g Araport11 -c Chr1 Chr2 Chr3 Chr4 Chr5 -a meta/ArabidopsisGenome/Araport11_GFF3_genes_transposons.201606.gff -s meta/ArabidopsisGenome/TAIR10_chr_count.txt meta/noMaskReads_Inputs_sampleSheet.csv data/ChIPQCreport/20inputs_noMask
```
> Saving 7 x 7 in image Error in .External2(C_X11, paste0("png::", filename), g$width, g$height,

$Troubleshoot, Check into R with\
```
capabilities()
```
Show that png not functional!
$Troubleshoot, Try to make it TRUE:
```
install.packages("png")
library(png)
```
$Troubleshoot, find a version of R where png is functional\
Version R/3.2.3\
Capabilities() show png TRUE!\
$FAIL, while installing rsvg package
$Troubleshoot, R3.5.2
- RSVG work
- Change R biocManager into V3.8
$FAIL at BiocManager::install("ChIPQC"), try version 3.5.2 but with v3.7 Biocmanager
$FAIL, probably because wrong Biocmanager installed: Change the script and use version 
$Troubleshoot:\
- Try R4.0.4 restart script and add: options(bitmapType='cairo')
- Try change the save .png into .pdf
- Try add a .Rprofile and add weird parameters:
```
curl -#LO https://www.rcac.purdue.edu/files/knowledge/run/examples/apps/r/Rprofile_example
mv -ib Rprofile_example ~/.Rprofile
.libPaths()
```
$FAIL at installing Bioconductor... Try install R package caro
```
png(tempfile(), type='cairo'); dev.off()
```
$Troubleshoots:
- Try disable X11
- Try R version 3.4.4 which capabilities OK for png (Fail at pkgs installation, installed it, but fail again)

$SOLUTION, Try download a version of R >4.2 with capabilities png OK, download it from R website
Download and transfer the .zip\
```
tar -zxvf path_to_file -C output_directory
chmod +x configure
./configure --with-pcre1 --with-readline=no
Make
```
Run with
```
/home/roule/R/R-4.2.0/bin/R
```
Now, launch the script using this version of R:\
```
/home/roule/R/R-4.2.0/bin/Rscript scripts/ChIPQC_forR42.R --indivReports -g Araport11 -c Chr1 Chr2 Chr3 Chr4 Chr5 -a meta/ArabidopsisGenome/Araport11_GFF3_genes_transposons.201606.gff -s meta/ArabidopsisGenome/TAIR10_chr_count.txt meta/noMaskReads_Inputs_sampleSheet.csv data/ChIPQCreport/20inputs_noMask
```

6. Vizualize mapped input sequences
**Download and install bedgraphtobigwig**\
```
URL=http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64
curl $URL/bedGraphToBigWig > /home/roule/GreenScreen/Software/bedGraphToBigWig
chmod +x /home/roule/GreenScreen/Software/bedGraphToBigWig
curl $URL/faToTwoBit > /home/roule/GreenScreen/Software/faToTwoBit
chmod +x /home/roule/GreenScreen/Software/faToTwoBit
```
To run command use:
```
/home/roule/GreenScreen/Software/bedGraphToBigWig
```
$FAIL, installing MACS2\
Need to install many stuff, so create a Conda environment as follow:
```
python3 -m venv /home/roule/GreenScreen/tutorial/GreenscreenProject/
source bin/activate #to activate and get into it
deactivate # to leave it
python3 -m pip install --upgrade pip #Install MACS2 within my custom conda environment
```
$FAIL, because MACS2 need Python2 (and not Python3)\
Install a custom version of Python2:
```
cd ~
mkdir tmp
cd tmp
wget https://www.python.org/ftp/python/2.7.15/Python-2.7.15.tgz
tar zxvf Python-2.7.15.tgz
cd Python-2.7.15 
./configure --prefix=$HOME/opt/python-2.7.15 --with-ensurepip=install
make
make install
```
Put this installed version as Python default:
```
Cd ~
export PATH=$HOME/opt/python-2.7.15/bin:$PATH
. ~/.bash_profile
which python #should indicate the path to 2.7.15 (~/opt/python-2.7.15/bin/python)
pip install virtualenv
virtualenv -p ~/opt/python-2.7.15/bin/python GSenv #To create Python2 environment
source GSenv/bin/activate #activate environment
python -m pip install --upgrade pip
pip install numpy
pip install MACS2
```
$FAIL, say Python must be 3.6 or over, so I installed a more recent verion:
```
wget https://www.python.org/ftp/python/3.8.4/Python-3.8.4.tgz
cd ~
mkdir tmp
cd tmp
wget https://www.python.org/ftp/python/3.8.4/Python-3.8.4.tgz
tar xvf Python-3.8.4.tgz
cd Python-3.8.4
./configure --enable-optimizations --with-ensurepip=install
make -j 8
/home/roule/tmp/Python-3.8.4/python # To run it
virtualenv --python=/home/roule/tmp/Python-3.8.4/python /home/roule/GreenScreen/tutorial/GreenscreenProject/ #To create virtual environment
source bin/activate # To activate
pip install git+https://github.com/ducksboard/libsaas.git#egg=libsaas
python3 -m pip install "MACS2"
```
> FAIL No module named '_ctypes'

$FAIL, need install *libfi*, I installed it as follow:
```
wget ftp://sourceware.org/pub/libffi/libffi-3.2.1.tar.gz
tar xvfz libffi-3.2.1.tar.gz
cd libffi-3.2.1
./configure --prefix=/usr/local/libffi/3_2_1
make
make install
```
Repeat installation but $FAIL again\
$SOLUTION, work within a CONDA environment
```
module load Anaconda/2019.10
conda create --name test1 python=3.6.8 r=3.6.0  #Create environment with desired version, use "conda search python" or "conda search R" to see available packages
conda create --name CondaGS python=3.6.8
conda init bash
conda activate CondaGS
python3 -m pip install "MACS2"
```
Modified script to call for peaks and launch it; then generate Greenscreen:
```
sbatch scripts/macs2_callpeaks_inputs.sh
sbatch scripts/generate_20input_greenscreenBed.sh 10 1000 2 #[1]=pvalue, [2]=merge parameter [3]=call for GS region if peak present in "2" inputs
```
The output is a bed format of the Arabidopsis Greenscreen region using 3 inputs
