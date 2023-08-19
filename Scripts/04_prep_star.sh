# ==============================================================================
# 04_prep_star.sh
# ==============================================================================



# install STAR
# ------------------------------------------------------------------------------
git clone https://github.com/alexdobin/STAR.git
cd STAR/source
make

# test working
/home/liucu/STAR/source/STAR --help 



# get GENCODE "PRI" for mouse
# ------------------------------------------------------------------------------
mkdir /home/biostats_share/liucu/GENCODE_Mouse
cd /home/biostats_share/liucu/GENCODE_Mouse
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M26/gencode.vM26.primary_assembly.annotation.gtf.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M26/GRCm39.primary_assembly.genome.fa.gz
gunzip *


  
# setup genome indices
# ------------------------------------------------------------------------------
cd /home/biostats_share/liucu/GENCODE_Mouse
STAR=/home/liucu/STAR/source/STAR
$STAR \
  --runThreadN 48 \
  --runMode genomeGenerate \
  --genomeDir /home/biostats_share/liucu/GENCODE_Mouse \
  --genomeFastaFiles /home/biostats_share/liucu/GENCODE_Mouse/GRCm39.primary_assembly.genome.fa \
  --sjdbGTFfile /home/biostats_share/liucu/GENCODE_Mouse/gencode.vM26.primary_assembly.annotation.gtf \
  --sjdbOverhang 150
  
  

# software versions
# ------------------------------------------------------------------------------

# $STAR --version
#     2.7.3a