# Microbial-Bioinformatics-101.2-A-Practical-Workshop-for-ShotGun-Metagenomics
The workshop and its computing environment have been specifically designed for Indonesian scientists working at The National Research and Innovation Agency (BRIN), Indonesia. This workshop offers BRIN staff comprehensive training in analyzing bacterial community sequence data using whole metagenome shotgun techniques.

#go to your working directory
cd /data/p283624

mkdir training

cd /data/p283624/training

## 00. Download sample metagenomic data

mkdir fastq

wget ftp.sra.ebi.ac.uk/vol1/fastq/ERR011/ERR011347/ERR011347_1.fastq.gz
wget ftp.sra.ebi.ac.uk/vol1/fastq/ERR011/ERR011347/ERR011347_2.fastq.gz

wget ftp.sra.ebi.ac.uk/vol1/fastq/ERR011/ERR011348/ERR011348_1.fastq.gz
wget ftp.sra.ebi.ac.uk/vol1/fastq/ERR011/ERR011348/ERR011348_2.fastq.gz

wget ftp.sra.ebi.ac.uk/vol1/fastq/ERR011/ERR011349/ERR011349_1.fastq.gz
wget ftp.sra.ebi.ac.uk/vol1/fastq/ERR011/ERR011349/ERR011349_2.fastq.gz

wget ftp.sra.ebi.ac.uk/vol1/fastq/ERR011/ERR011350/ERR011350_1.fastq.gz
wget ftp.sra.ebi.ac.uk/vol1/fastq/ERR011/ERR011350/ERR011350_2.fastq.gz

#how to copy and paste

cp *#files you want to copy #folder you are aiming for
cp *fastq.gz fastq

#how to cut and paste

mv *#files you want to cut/move #folder you are aiming for
mv *fastq.gz fastq

#how to copy and paste when the file you want to copy are not in your working directory

cp ~/#folder where your files exist/*.fastq.gz #folder you are aiming for
cp ~/data/p283624/training*.fastq.gz fastq

#how to move your folder when the file you want to copy are not in your working directory

mv ~/#folder where your files exist/*.fastq.gz #folder you are aiming for
mv ~/data/p283624/training*.fastq.gz fastq

###################### 01. Quality checking ############################

### Install fastQC
https://www.bioinformatics.babraham.ac.uk/projects/fastqc/

conda install -c bioconda fastqc

#Activate FastQC software
module load FastQC

#Checking if FastQC is loaded
fastqc -- help

#Checking the version of fastQC
fastqc -- version

### Run fastQC with all samples
mkdir QC #first we have to make a folder for the output
fastqc --outdir QC --threads 18 ./fastq/*.fastq.gz

fastqc --outdir #Name of output folder --threads #number of CPU ./#folder where your fastq files exist/*.fastq.gz

#threads is for the number of processors needed to compute such code, mine is 24 max (computer at RUG)
#to check the number of CPUS/core/processors : top, afterwards press 1

#However, peregrine in RUG does not allow direct computing for such code, we have to prepare ".sh" file and submit the job.

#so it will look like this:

nano fastqc.sh

#!/bin/bash

#SBATCH --job-name=fastqc
#SBATCH --time =24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=22
#SBATCH --mem= 100GB
#SBATCH -o fastqc-%j.out
#SBATCH -e fastqc-%j.error
#SBATCH -- mail-type=FAIL,END
#SBATCH -- mail-user=p.c.mawarda@rug.nl
#SBATCH --group = is_elements

module load FastQC

cd $pwd

fastqc --outdir QC --threads 22 ./fastq/*.fastq.gz

#Save the sh file by pressing Control X then Yes

#Check file again

cat fastqc.sh

#Run the code by typing

sbatch fastqc.sh

#After it's done, check your quality score: the good one is those in the green zone, usually >28
#Check your minimum length as well, the quality score and minimum length are parameters you need for Trimmomatic in slidingwindow and min length menu

################ 02. Quality filtering #####################

### Install Trimmomatic
https://github.com/usadellab/Trimmomatic
http://www.usadellab.org/cms/?page=trimmomatic

conda install -c bioconda trimmomatic

### Run the program
mkdir trim

module load Trimmomatic

#in github, you can find different adapter to trim, depending on the sequencing machine, what I used in this tutorial was TruSeq2 because the sequencing machine was GAII machine; HiSeq and MiSeq machine usually use TruSeq3 adapter

#Download the adapter from github: Click the adapter folder, search for TruSeq2, Click Raw, and copy the link above, wget, and paste

wget https://raw.githubusercontent.com/usadellab/Trimmomatic/main/adapters/TruSeq2-SE.fa

#Creating sh file to submit the job

nano trimmomatic.sh

#!/bin/bash

#SBATCH --job-name=trimmomatic
#SBATCH --time =24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=22
#SBATCH --mem= 100GB
#SBATCH -o trimmomatic-%j.out
#SBATCH -e trimmomatic-%j.error
#SBATCH -- mail-type=FAIL,END
#SBATCH -- mail-user=p.c.mawarda@rug.nl

module load Trimmomatic

cd $pwd

for i in fastq/*_1.fastq.gz; do base=$(echo "$i" | sed -e 's/fastq[/]\(.*\)_1.fastq.gz/\1/'); java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE -threads 22 $i fastq/${base}_2.fastq.gz trim/${base}_trim_1.fq trim/${base}_unpaired_1.fq trim/${base}_trim_2.fq trim/${base}_unpaired_2.fq ILLUMINACLIP:TruSeq2-PE.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:75; done

#Save the sh file by pressing Control X then Yes

#Check file again

cat trimmomatic.sh

#Run the code by typing

sbatch trimmomatic.sh

#From FastQC results, the sliding window is 15, Minimum Lenth is 75
#Interpretation Sliding Window, every 4 base if the quality score is below 15, then drop !

######################################### 03. Assembly #############################################

## There are 2 options: individual or co-assembly assembly:
## Mode 1: Individual assembly -> run the assembly on each sample separatelly. 

### 03.1 metaSPAdes
#### Install metaSPAdes assembler
https://github.com/ablab/spades
https://cab.spbu.ru/files/release3.15.5/manual.html

conda install -c bioconda spades

#or

wget https://github.com/ablab/spades/releases/download/v3.15.5/SPAdes-3.15.5-Linux.tar.gz
tar xzvf SPAdes-3.15.5-Linux.tar.gz

#### Run metaSPAdes Mode 1 #####

module load SPAdes

#creating sh file

nano assembly_mode1.sh

#!/bin/bash

#SBATCH --job-name=metaspades_mode1
#SBATCH --time =24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=22
#SBATCH --mem= 100GB
#SBATCH -o metaspades_mode1-%j.out
#SBATCH -e metaspades_mode1-%j.error
#SBATCH -- mail-type=FAIL,END
#SBATCH -- mail-user=p.c.mawarda@rug.nl

module load SPAdes
cd $pwd

for i in trim/*_trim_1.fq; do base=$(echo "$i" | sed -e 's/trim[/]\(.*\)_1.fq/\1/'); base2=$(echo "$i" | sed -e 's/trim[/]\(.*\)_trim_1.fq/\1/'); cat trim/${base2}_unpaired_1.fq trim/${base2}_unpaired_2.fq > trim/${base2}_unpaired.fq; metaspades.py -1 $i -2 trim/${base}_2.fq -s trim/${base2}_unpaired.fq -t 22 -m 250 -o ${base2}_metaspades; done

#Save the sh file by pressing Control X then Yes

#Check file again

cat assembly_mode1.sh

#Run the code by typing

sbatch assembly_mode1.sh

#### Run metaSPAdes Mode 2 #####

## Mode 2: Co-assembly -> Concatenate all samples into one file. However this needs higher RAM and heavy computation. Hence mode 1 is more preferred.

cat trim/ERR*_trim_1.fq > read_R1.fq
cat trim/ERR*_trim_2.fq > read_R2.fq
cat trim/ERR*_unpaired_1.fq > read_unpaired_R1.fq
cat trim/ERR*_unpaired_2.fq > read_unpaired_R2.fq
cat read_unpaired_*.fq > unpaired_reads.fq

module load SPAdes

#creating sh file

nano assembly_mode2.sh

#!/bin/bash

#SBATCH --job-name=metaspades_mode2
#SBATCH --time =24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=22
#SBATCH --mem= 100GB
#SBATCH -o metaspades_mode2-%j.out
#SBATCH -e metaspades_mode2-%j.error
#SBATCH -- mail-type=FAIL,END
#SBATCH -- mail-user=p.c.mawarda@rug.nl

module load SPAdes
cd $pwd

metaspades.py -1 read_R1.fq -2 read_R2.fq -s unpaired_reads.fq -t 22 -m 100 -o metaSPAdes

#Save the sh file by pressing Control X then Yes

#Check file again

cat assembly_mode2.sh

#Run the code by typing

sbatch assembly_mode2.sh

#In every folder, we will find different k mers, K21, K33, and K55


#K21: The sequence will be cut into pieces, with a length of 21 bases, and then reassembled.
#K33: After merging, the sequence will be cut again to a length of 33 bases, and then reassembled.
#K55: After merging, the sequence will be cut again to a length of 55 bases, and then reassembled. The presence of several Kmers is to optimize the final result.
#The circularity of our assembly can be determined by examining the scaffolds.gfa file using Bandage software.
#There are two important files from the assembly process: contigs.fasta and scaffolds.fasta. Scaffolds are contigs that we reassemble, making them longer.

### 03.2 MetaVelvet
#### Install metavelvet
http://metavelvet.dna.bio.keio.ac.jp/

conda install -c bioconda metavelvet

#### Helper scripts
wget https://raw.githubusercontent.com/dzerbino/velvet/master/contrib/observed-insert-length.pl/observed-insert-length.pl

## k-mers 57 #>> It can be adjusted according to your samples

## Mode 1 individual
for i in trim/*_trim_1.fq; do \
base=$(echo "$i" | sed -e 's/trim[/]\(.*\)_1.fq/\1/'); \
base2=$(echo "$i" | sed -e 's/trim[/]\(.*\)_trim_1.fq/\1/'); \
echo -e "velveth MetaVelveth_${base2}_k57 57 -fmtAuto -shortPaired -separate $i trim/${base}_2.fq"; \
echo -e "velvetg MetaVelveth_${base2}_k57 -ins_length auto -exp_cov auto -cov_cutoff auto"; done

perl observed-insert-length.pl MetaVelveth_ERR011347_k57 # get insert length

for i in MetaVelveth_*_k57; do \
meta-velvetg $i -ins_length XX; done

## Mode 2 co assembly
velveth MetaVelveth_k57 57 -fmtAuto -shortPaired -separate read_R1.fq read_R2.fq
velvetg MetaVelveth_k57 -ins_length auto -exp_cov auto -cov_cutoff auto
perl observed-insert-length.pl MetaVelveth_k57 # get insert length
meta-velvetg MetaVelveth -ins_length XX
