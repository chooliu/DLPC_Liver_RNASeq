# ==============================================================================
# 02_cutadapt.sh
# ==============================================================================



# my project folder
cd /home/biostats_share/liucu/CIDAP20059SGhosh/Scripts



# helper script to run cutadapt
# ------------------------------------------------------------------------------
cat > 02_cutadapt_prep.sh 
#!/bin/bash
PROJECTDIR=/home/biostats_share/liucu/CIDAP20059SGhosh
CUTADAPT=~/.local/bin/cutadapt
TRIMMEDFOLDER=$PROJECTDIR/trimmed_reads
cd $PROJECTDIR/raw_reads

for READ1 in `ls *R1*.fastq.gz`
do
  
  READ2=$(echo | awk -v read1file=$READ1 '{ gsub("R1", "R2", read1file); print read1file }')
  OUT1=$(echo | awk -v read1file=$READ1 '{ print "trim_" read1file }')
  OUT2=$(echo | awk -v read2file=$READ2 '{ print "trim_" read2file }')
  
  echo $CUTADAPT -u 4 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
  --trim-n -j 24 -q 20 --nextseq-trim=20 -m 50 -o $TRIMMEDFOLDER/$OUT1 \
  -p $TRIMMEDFOLDER/$OUT2 $PROJECTDIR/raw_reads/$READ1 $PROJECTDIR/raw_reads/$READ2

done



# run cutadapt
# ------------------------------------------------------------------------------
sh 02_cutadapt_prep.sh &> 02_run_cutadapt.sh

conda activate py38
sh 02_run_cutadapt.sh > 02_run_cutadapt.log



# software versions
# ------------------------------------------------------------------------------

# $CUTADAPT --version
# 2.7
