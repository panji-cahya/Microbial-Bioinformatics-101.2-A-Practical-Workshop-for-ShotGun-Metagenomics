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

### 01. Quality checking 

### 01.1 Install fastQC
https://www.bioinformatics.babraham.ac.uk/projects/fastqc/

conda install -c bioconda fastqc

#Activate FastQC software
module load FastQC

#Checking if FastQC is loaded
fastqc -- help

#Checking the version of fastQC
fastqc -- version

### 01.2 Run fastQC with all samples
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

## 02. Quality filtering 

### 02.1 Install Trimmomatic
https://github.com/usadellab/Trimmomatic
http://www.usadellab.org/cms/?page=trimmomatic

conda install -c bioconda trimmomatic

### 02.2 Run the program
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

### 03. Assembly 
There are 2 options: individual or co-assembly assembly:

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

#Before we begin, let's check how many conda environments we have:

Module load Anaconda3
conda env list

#Making conda environment

conda create -n bowtie2 #then click yes

conda activate bowtie2

conda install -c bioconda bowtie2 samtools #then click yes

#How to check if there's an installed program in the environment we are activating

conda list

#How to check if an installed program is active

bowtie2 --help

#The difference between ll and ls; ll lists all the files in the folder along with comprehensive information such as size, etc., while ls only lists the names of the files in that folder.

##### 01. Mapping reads to contigs with Bowtie2

https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml

#for short reads, we can use bowtie2 or Burrows Wheeler Aligner (BWA), for long reads we can use Minimap2

#Why do you have to map your contigs? maxbin2 and metabat2 needs information about abundance, and to know the abundance we have to map back the contigs.. how many reads in a contig we can map back..

###### Mode 1. Individual Assembly

mkdir bowtie2_ref

nano bowtie2_mode1.sh

#!/bin/bash

#SBATCH --job-name=bowtie2_mode1
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --partition=himem
#SBATCH --ntasks=22
#SBATCH --mem=100GB
#SBATCH -o bowtie2_mode1-%j.out
#SBATCH -e bowtie2_mode1-%j.error
#SBATCH --mail-type=FAIL,END

module load Anaconda3

conda activate bowtie2

for i in ERR*/contigs.fasta; do base=$(echo "$i" | cut -c1-9); bowtie2-build $i bowtie2_ref/$base; done

for i in ERR*/contigs.fasta; do base=$(echo "$i" | cut -c1-9); bowtie2 -x bowtie2_ref/$base -1 trim/${base}_trim_1.fq -2 trim/${base}_trim_2.fq -U trim/${base}_unpaired_1.fq,trim/${base}_unpaired_2.fq -p 22 | samtools view -bSu - | samtools sort -@ 22 -o ${base}_model1_sorted.bam; done

for x in *_model1_sorted.bam; do samtools index $x; done

#cut -c1-9 meaning remove 9 strings from behind

#Save the sh file by pressing Control X then Yes

#Check file again

cat bowtie2_mode1.sh

#Run the code by typing

sbatch bowtie2_mode1.sh


###### Mode 2. Co-assembly

nano bowtie2_mode2.sh

#!/bin/bash

#SBATCH --job-name=bowtie2_mode2
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --partition=himem
#SBATCH --ntasks=22
#SBATCH --mem=100GB
#SBATCH -o bowtie2_mode2-%j.out
#SBATCH -e bowtie2_mode2-%j.error
#SBATCH --mail-type=FAIL,END

module load Anaconda3

conda activate bowtie2

bowtie2-build metaSPAdes/contigs.fasta contigs_ref

list="ERR011347 ERR011348 ERR011349 ERR011350"

for i in $list; do bowtie2 -x contigs_ref -1 trim/${i}_trim_1.fq -2 trim/${i}_trim_2.fq -U trim/${i}_unpaired_1.fq,trim/${i}_unpaired_2.fq -p 22 | samtools view -bSu - | samtools sort -@ 22 -o ${i}_model2_sorted.bam; done

for x in *_model2_sorted.bam; do samtools index $x; done

#Save the sh file by pressing Control X then Yes

#Check file again

cat bowtie2_mode2.sh

#Run the code by typing

sbatch bowtie2_mode2.sh

#BAM visualization can be done using IGV. Download the BAM file, upload it to IGV, similar to .qzv files in QIIME2.


#### 02. Binning 


##### Metabat2 

https://bitbucket.org/berkeleylab/metabat/src/master/

#metabat is based on depth

conda create -n metabat2
conda activate metabat2
conda install -c bioconda metabat2

### Mode 1. Individual Binning Using Metabat 2

nano metabat2_mode1.sh

#!/bin/bash

#SBATCH --job-name=metabat2_mode1
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --partition=himem
#SBATCH --ntasks=22
#SBATCH --mem=100GB
#SBATCH -o metabat2_mode1-%j.out
#SBATCH -e metabat2_mode1-%j.error
#SBATCH --mail-type=FAIL,END

module load Anaconda3

conda activate metabat2

mkdir metabat2_bin

for i in ERR*/contigs.fasta; do base=$(echo "$i" | cut -c1-9); jgi_summarize_bam_contig_depths --outputDepth ${base}_model1_depth.txt ${base}_model1_sorted.bam; mkdir metabat2_bin/$base; metabat2 -i $i -a ${base}_model1_depth.txt -t 22 -o metabat2_bin/$base/$base; done

#Save the sh file by pressing Control X then Yes

#Check file again

cat metabat2_mode1.sh

#Run the code by typing

sbatch metabat2_mode1.sh

#### Mode 2. Co-assembly 

nano metabat2_mode2.sh

#!/bin/bash

#SBATCH --job-name=metabat2_mode2
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --partition=himem
#SBATCH --ntasks=22
#SBATCH --mem=100GB
#SBATCH -o metabat2_mode2-%j.out
#SBATCH -e metabat2_mode2-%j.error
#SBATCH --mail-type=FAIL,END

module load Anaconda3

conda activate metabat2

mkdir metabat2_co

jgi_summarize_bam_contig_depths --outputDepth co_depth.txt *_model2_sorted.bam
metabat2 -i metaSPAdes/contigs.fasta -a co_depth.txt -t 22 -o metabat2_co/metabat2_co

#Save the sh file by pressing Control X then Yes

#Check file again

cat metabat2_mode2.sh

#Run the code by typing

sbatch metabat2_mode2.sh

##### Maxbin2 

https://sourceforge.net/projects/maxbin2/

#maxbin is based on abundance

conda create -n maxbin2
conda activate maxbin2
conda install -c bioconda maxbin2

##### Mode 1. Individual Using Maxbin 2

nano maxbin2_mode1.sh

#!/bin/bash

#SBATCH --job-name=maxbin2_mode1
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --partition=himem
#SBATCH --ntasks=22
#SBATCH --mem=100GB
#SBATCH -o maxbin2_mode1-%j.out
#SBATCH -e maxbin2_mode1-%j.error
#SBATCH --mail-type=FAIL,END

module load Anaconda3

conda activate maxbin2

mkdir maxbin2_bin

for i in *_model1_sorted.bam; do base=$(echo "$i" | cut -c1-9); samtools idxstats $i > ${base}.idxstats.txt; sed -i '$ d' ${base}.idxstats.txt; cut -f1,3 ${base}.idxstats.txt > ${base}_model1.abund; mkdir maxbin2_bin/$base; run_MaxBin.pl -contig ${base}_metaspades/contigs.fasta -out maxbin2_bin/$base/$base -thread 22 -abund ${base}_model1.abund; done

#Save the sh file by pressing Control X then Yes

#Check file again

cat maxbin2_mode1.sh

#Run the code by typing

sbatch maxbin2_mode1.sh

#cut -f1,3 artinya select column no.1 dan 3

##### Mode 2. Co-assembly

nano list_abundant

ERR011347.abund
ERR011348.abund
ERR011349.abund
ERR011350.abund

nano maxbin2_mode2.sh

#!/bin/bash

#SBATCH --job-name=maxbin2_mode2
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --partition=himem
#SBATCH --ntasks=22
#SBATCH --mem=100GB
#SBATCH -o maxbin2_mode2-%j.out
#SBATCH -e maxbin2_mode2-%j.error
#SBATCH --mail-type=FAIL,END

module load Anaconda3
conda activate maxbin2

mkdir maxbin2_co

list="ERR011347 ERR011348 ERR011349 ERR011350"

for i in $list; do samtools idxstats ${i}_model2_sorted.bam > ${i}_model2.idxstats.txt; sed -i '$ d' ${i}_model2.idxstats.txt; cut -f1,3 ${i}_model2.idxstats.txt > ${i}.abund; done

run_MaxBin.pl -contig metaSPAdes/contigs.fasta -out maxbin2_co/maxbin2_co -thread 22 -abund_list list_abundant

#Save the sh file by pressing Control X then Yes

#Check file again

cat maxbin2_mode2.sh

#Run the code by typing

sbatch maxbin2_mode2.sh

###### 03. Check bins quality with checkM

https://github.com/Ecogenomics/CheckM
https://ecogenomics.github.io/CheckM/

conda create -n checkm
conda activate checkm
conda install -c bioconda checkm-genome

###### Metabat2

###### Mode 1. Individual

nano checkM_metabat2_mode1.sh

#!/bin/bash

#SBATCH --job-name=checkM_metabat2_mode1
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --partition=himem
#SBATCH --ntasks=22
#SBATCH --mem=100GB
#SBATCH -o checkM_metabat2_mode1-%j.out
#SBATCH -e checkM_metabat2_mode1-%j.error
#SBATCH --mail-type=FAIL,END

module load Anaconda3

conda activate checkm

cd metabat2_bin

for i in *; do checkm lineage_wf -r -t 22 --pplacer_threads 22 -x fa $i ${i}_checkM; checkm qa -t 22 ${i}_checkM/lineage.ms ${i}_checkM/ > ${i}.checkM; done

#-x fa (Adapting the file extension, some are .fa, and some are .fna.)

#Save the sh file by pressing Control X then Yes

#Check file again

cat checkM_metabat2_mode1.sh

#Run the code by typing

sbatch checkM_metabat2_mode1.sh

#Then we analyze .checkm to see the qualiity of MAGs. See the % completeness and contamination of the MAGS. Low if %completeness less than 50%, middle if it's 60-80%, high if it's >80%

###### Mode 2. Co-assembly

nano checkM_metabat2_mode2.sh

#!/bin/bash

#SBATCH --job-name=checkM_metabat2_mode2
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --partition=himem
#SBATCH --ntasks=22
#SBATCH --mem=100GB
#SBATCH -o checkM_metabat2_mode2-%j.out
#SBATCH -e checkM_metabat2_mode2-%j.error
#SBATCH --mail-type=FAIL,END

module load Anaconda3
conda activate checkm

checkm lineage_wf -r -t 22 --pplacer_threads 22 -x fa metabat2_co metabat2_co_checkM
checkm qa -t 22 metabat2_co_checkM/lineage.ms metabat2_co_checkM > metabat2_co.checkM

#-x fa (menyesuaikan ekstensi file nya, ada yang fa ada yang fna

#Save the sh file by pressing Control X then Yes

#Check file again

cat checkM_metabat2_mode2.sh

#Run the code by typing

sbatch checkM_metabat2_mode2.sh

###### Maxbin2

##### Mode 1. Individual Using Maxbin2

nano checkM_maxbin2_mode1.sh

#!/bin/bash

#SBATCH --job-name=checkM_maxbin2_mode1
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --partition=himem
#SBATCH --ntasks=22
#SBATCH --mem=100GB
#SBATCH -o checkM_maxbin2_mode1-%j.out
#SBATCH -e checkM_maxbin2_mode1-%j.error
#SBATCH --mail-type=FAIL,END

module load Anaconda3
conda activate checkm

cd maxbin2_bin

for i in *; do \
checkm lineage_wf -r -t 22 --pplacer_threads 22 -x fasta $i ${i}_checkM; \
checkm qa -t 22 ${i}_checkM/lineage.ms ${i}_checkM/ > ${i}.checkM; \
done

###### Mode 2. Co-assembly using MaxBin

nano checkM_maxbin2_mode2.sh

#!/bin/bash

#SBATCH --job-name=checkM_maxbin2_mode2
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --partition=himem
#SBATCH --ntasks=22
#SBATCH --mem=100GB
#SBATCH -o checkM_maxbin2_mode2-%j.out
#SBATCH -e checkM_maxbin2_mode2-%j.error
#SBATCH --mail-type=FAIL,END

module load Anaconda3
conda activate checkm

checkm lineage_wf -r -t 22 --pplacer_threads 22 -x fasta maxbin2_co maxbin2_co_checkM
checkm qa -t 22 maxbin2_co_checkM/lineage.ms maxbin2_co_checkM/ > maxbin2_co.checkM

### Taxonomy Annotation 

## First create diamond database with nr database

## Download taxonomy nr database

wget https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz
wget -c ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz 
tar -zxvf taxdump.tar.gz

nano create_diamond_db.sh

#!/bin/bash

#SBATCH --job-name=create_diamond_db
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --partition=himem
#SBATCH --ntasks=22
#SBATCH --mem=100GB
#SBATCH -o create_diamond_db-%j.out
#SBATCH -e create_diamond_db-%j.error
#SBATCH --mail-type=FAIL,END

module load Anaconda3
conda activate diamond

diamond makedb --in nr.gz --db nr.dmnd --taxonmap prot.accession2taxid.gz --taxonnodes nodes.dmp --taxonnames names.dmp --threads 22

sbatch create_diamond_db.sh

#this database is required at later steps after ORF prediction with prodigal

**********************************************************
**********************************************************
**********************************************************

## This tutorial uses mode 1: co-assembly as the input

**********************************************************
**********************************************************
**********************************************************


## ORF prediction with prodigal

######### if we do it from contigs ###########

nano prodigal.sh

#!/bin/bash

#SBATCH --job-name=prodigal
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --partition=himem
#SBATCH --ntasks=22
#SBATCH --mem=100GB
#SBATCH -o prodigal-%j.out
#SBATCH -e prodigal-%j.error
#SBATCH --mail-type=FAIL,END

module load Anaconda3
conda activate prodigal

mkdir cds
prodigal -q -m -p meta -i metaSPAdes/contigs.fasta -a cds/cds.faa -d cds/cds.fna -f gff -o cds/cds.gff

sbatch prodigal.sh

####### If we do it from MAGS

nano prodigal_mags.sh

#!/bin/bash

#SBATCH --job-name=prodigal_mags
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=22
#SBATCH --mem=50GB
#SBATCH -o prodigal_mags-%j.out
#SBATCH -e prodigal_mags-%j.error
#SBATCH --mail-type=FAIL,END

module load Anaconda3
conda activate prodigal

list_mags="metabat2_co.27 metabat2_co.16"

for i in $list_mags; do mkdir $i; prodigal -q -m -i metabat2_co/${i}.fa -a $i/${i}_cds_mags.faa -d $i/${i}_cds_mags.fna -f gff -o $i/${i}_cds_mags.gff; done

###### Convert gff to gtf ########


nano gff2gtf.sh

##################
#!/bin/bash

infile=$1

if [ "$infile" == "" ] ; then
    echo "Usage: gff2gtf.sh <gff file>"
    exit 0
fi

grep -v "#" $infile | grep "ID=" | cut -f1 -d ';' | sed 's/ID=//g' | cut -f1,4,5,7,9 |  awk -v OFS='\t' '{print $1,"Prodigal","CDS",$2,$3,".",$4,".","gene_id " $5}'
#####################

chmod 777 gff2gtf.sh

nano get_gtf.sh

#!/bin/bash

#SBATCH --job-name=get_gtf
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --partition=himem
#SBATCH --ntasks=22
#SBATCH --mem=100GB
#SBATCH -o get_gtf-%j.out
#SBATCH -e get_gtf-%j.error
#SBATCH --mail-type=FAIL,END

./gff2gtf.sh cds/cds.gff > cds/cds.gtf

sbatch get_gtf.sh


###### Taxonomy classification with nr

nano run_diamond_nr.sh

#!/bin/bash

#SBATCH --job-name=run_diamond_nr
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --partition=himem
#SBATCH --ntasks=22
#SBATCH --mem=500GB
#SBATCH -o run_diamond_nr-%j.out
#SBATCH -e run_diamond_nr-%j.error
#SBATCH --mail-type=FAIL,END

modul load Anaconda3
conda activate diamond

mkdir diamond_annotation
diamond blastp -q cdsprodigal/cds.faa -p 22 -d db/nr.dmnd -e 0.001 --id 60 -b 8 -k 10 -o diamond_annotation/diamond.out \
-f 6 qseqid sseqid pident length mismatch gapopen evalue bitscore staxids

sbatch run_diamond_nr.sh

####### How to know how many sequences before submitting to cog or kegg ##########

grep '^>' cds.faa | wc -l

## Get abundance information with subread

conda create -n subread
conda activate subread

conda install -c bioconda subread

nano count_abundance.sh

#!/bin/bash

#SBATCH --job-name=count_abundance
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --partition=himem
#SBATCH --ntasks=22
#SBATCH --mem=100GB
#SBATCH -o count_abundance-%j.out
#SBATCH -e count_abundance-%j.error
#SBATCH --mail-type=FAIL,END

modul load Anaconda3
conda activate subread 

featureCounts -T 22 -t CDS -g gene_id -a cdsprodigal/cds.gtf -o cdsprodigal/CDS_abundance.tsv ERR011347_sorted.bam ERR011348_sorted.bam ERR011349_sorted.bam ERR011350_sorted.bam

sbatch count_abundance.sh

##### Functional annotation with hmmer 


conda create -n hmmer
conda activate hmmer
conda install -c bioconda pfam_scan
conda install -c bioconda hmmer

wget https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
wget https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.dat.gz
wget https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/active_site.dat.gz

gunzip *.gz

hmmpress Pfam-A.hmm

wget http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.clans.tsv.gz
gunzip Pfam-A.clans.tsv.gz

nano run_hmmer.sh

#!/bin/bash

#SBATCH --job-name=run_hmmer
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --partition=himem
#SBATCH --ntasks=22
#SBATCH --mem=100GB
#SBATCH -o run_hmmer-%j.out
#SBATCH -e run_hmmer-%j.error
#SBATCH --mail-type=FAIL,END

modul load Anaconda3
conda activate hmmer

mkdir hmmer_annotation
cd hmmer_annotation

hmmsearch --domtblout pfamhmmer.out -E 1e-10 --cpu 22 Pfam-A.hmm data/p283624/panji_metagenomics/cdsprodigal/cds.faa

# OR

pfam_scan.pl -fasta ../cds/cds.faa -dir /home_sbi_cold/60110700001/training/hmmer_annotation -outfile pfamScan.out -as -cpu 18

sbatch run_hmmer.sh

