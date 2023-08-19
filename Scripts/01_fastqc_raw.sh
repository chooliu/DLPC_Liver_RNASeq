# ==============================================================================
# 01_fastqc_raw.sh
# ==============================================================================



# my project folder
cd /home/biostats_share/liucu/CIDAP20059SGhosh/Scripts



# make a script called "runfastqc_raw.sh"
cat > 01_runfastqc_raw.sh 
#!/bin/bash
FASTQC=/home/liucu/FastQC/fastqc
TMPFOLDER=/home/biostats_share/liucu/tmp
PROJECTDIR=/home/biostats_share/liucu/CIDAP20059SGhosh/
mkdir $PROJECTDIR/fastqc_raw

$FASTQC `ls $PROJECTDIR/raw_reads/*.fastq.gz` -o \
$PROJECTDIR/fastqc_raw/ --noextract -d $TMPFOLDER -t 10



# run fastqc on raw fastqs

sh 01_runfastqc_raw.sh > 01_runfastqc_raw.log



# aggregate with multiqc 

conda activate py38
cd $PROJECTDIR/fastqc_raw
multiqc .  -o multiqc_fastqc_raw



# software versions
# ------------------------------------------------------------------------------

# $FASTQC --version
#     FastQC v0.11.8
# conda activate py38; conda list;
#     multiqc, version 1.9 (pypi_0 conda install)




