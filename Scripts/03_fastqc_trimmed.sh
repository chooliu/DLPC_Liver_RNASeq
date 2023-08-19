# ==============================================================================
# 03_fastqc_trimmed.sh
# ==============================================================================



# my project folder
cd /home/biostats_share/liucu/CIDAP20059SGhosh/Scripts



# make a script called "runfastqc_trimmed.sh"
# ------------------------------------------------------------------------------
cat > 03_runfastqc_trimmed.sh 
#!/bin/bash
FASTQC=/home/liucu/FastQC/fastqc
TMPFOLDER=/home/biostats_share/liucu/tmp
PROJECTDIR=/home/biostats_share/liucu/CIDAP20059SGhosh/
mkdir $PROJECTDIR/fastqc_trimmed

$FASTQC `ls $PROJECTDIR/trimmed_reads/*.fastq.gz` \
  -o $PROJECTDIR/fastqc_trimmed/ --noextract -d $TMPFOLDER -t 12



# run fastqc on cutadapt-trimmed fastqs
# ------------------------------------------------------------------------------

sh 03_runfastqc_trimmed.sh > 03_runfastqc_trimmed.log




# aggregate with multiqc 
# ------------------------------------------------------------------------------

cd $PROJECTDIR/fastqc_trimmed
conda activate py38
multiqc . -o multiqc_fastqc_trimmed 



# software versions
# ------------------------------------------------------------------------------

# $FASTQC --version
#     FastQC v0.11.8
# conda activate py38; conda list;
#     multiqc, version 1.9 (pypi_0 conda install)
