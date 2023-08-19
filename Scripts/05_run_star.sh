# ==============================================================================
# 05_run_star.sh
# ==============================================================================


# proj folder
cd /home/biostats_share/liucu/CIDAP20059SGhosh/Scripts



# trimmed .fq --> align
# ------------------------------------------------------------------------------
cat > 05_prep_star_alignment_cmds.sh 
#!/bin/bash
STAR=/home/liucu/STAR/source/STAR
PROJECTDIR=/home/biostats_share/liucu/CIDAP20059SGhosh
TRIMMEDFOLDER=$PROJECTDIR/trimmed_reads
ALIGNEDFOLDER=$PROJECTDIR/aligned_reads

for R1trimmed in `ls $TRIMMEDFOLDER/*R1*.fastq.gz`
  do

    R2trimmed=$(echo | awk -v read1file=$READ1 '{ gsub("R1", "R2", read1file); print read1file }')
    OUTNAME=$(echo $R1trimmed | awk -F 'trim_|_L001_' '{print $2 }')
    
    echo $STAR --runThreadN 12 --outFilterScoreMinOverLread 0.33 \
      --outFilterMatchNminOverLread 0.33 --outSAMtype BAM Unsorted \
      --quantMode GeneCounts --genomeDir /home/biostats_share/liucu/GENCODE_Mouse_v26 \
      --readFilesCommand zcat --readFilesIn \
      $R1trimmed $R2trimmed --outFileNamePrefix $ALIGNEDFOLDER/${OUTNAME}_
done



# run STAR helpers
# ------------------------------------------------------------------------------

mkdir ../aligned_reads
sh 05_prep_star_alignment_cmds.sh > 05_run_star_alignment.sh
sh 05_run_star_alignment.sh > 05_run_star_alignment.log



# software versions
# ------------------------------------------------------------------------------

# $STAR --version
#     2.7.3a