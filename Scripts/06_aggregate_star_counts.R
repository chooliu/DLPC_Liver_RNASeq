# ==============================================================================
# 06_aggregate_star_counts.R
# ==============================================================================



# STAR counts per sample --> count matrix
# ------------------------------------------------------------------------------

R
library(tidyverse)
library(data.table)
library(rtracklayer)

setwd('/home/biostats_share/liucu/CIDAP20059SGhosh/')

star_filepaths <-
  list.files("aligned_reads", pattern = "*ReadsPerGene.out.tab",
             recursive = T, full.names = T)
sample_names <-
  star_filepaths %>%
  str_split(., pattern = "_|\\/") %>%
  map_chr(~ .[3]) %>%
  paste0("Sample", .)

gene_names <-
  read_tsv(star_filepaths[1], col_names = F) %>% select(1) %>% .$X1

counts <-
  map(star_filepaths, function(x) { fread(x) %>% .[ , V3] }) %>% # V3 = sense strand
  set_names(sample_names) %>%
  as_tibble %>%
  select(names(.) %>% sort) %>%
  bind_cols(Feature = gene_names, .) %>%
  slice(-(1:4)) %>%
  arrange(Feature)


dir.create("Output")
write_tsv(counts, "./Output/20210307_counts.txt")



# Ensembl --> gene symbol
# ------------------------------------------------------------------------------

gffannot <-
  readGFF("/home/biostats_share/liucu/GENCODE_Mouse_v26/gencode.vM26.annotation.gtf")
gene_annotations_short <-
  gffannot %>%
  transmute(ENSG = gene_id, symbol = gene_name, gene_type) %>%
  filter(!duplicated(paste0(ENSG, symbol)))

write_tsv(gene_annotations_short,
          "./Output/gencode_annotations.txt")




# software versions
# (different from main analysis sessionInfo() since run on biostats servers)
# ------------------------------------------------------------------------------

# sessionInfo()
# R version 4.0.3 (2020-10-10)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: CentOS Linux 7 (Core)
# 
# Matrix products: default
# BLAS/LAPACK: /usr/lib64/libopenblasp-r0.3.3.so
# 
# locale:
#   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
# [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
# [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
# [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
# [9] LC_ADDRESS=C               LC_TELEPHONE=C            
# [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
#   [1] parallel  stats4    stats     graphics  grDevices utils     datasets 
# [8] methods   base     
# 
# other attached packages:
#   [1] rtracklayer_1.50.0   GenomicRanges_1.42.0 GenomeInfoDb_1.26.7 
# [4] IRanges_2.24.1       S4Vectors_0.28.1     BiocGenerics_0.36.1 
# [7] data.table_1.14.2    forcats_0.5.1        stringr_1.4.0       
# [10] dplyr_1.0.7          purrr_0.3.4          readr_2.0.2         
# [13] tidyr_1.1.4          tibble_3.1.0         ggplot2_3.3.3       
# [16] tidyverse_1.3.1     
# 
# loaded via a namespace (and not attached):
#   [1] Rcpp_1.0.7                  lattice_0.20-41            
# [3] lubridate_1.7.10            Rsamtools_2.6.0            
# [5] Biostrings_2.58.0           assertthat_0.2.1           
# [7] utf8_1.1.4                  R6_2.5.0                   
# [9] cellranger_1.1.0            backports_1.2.1            
# [11] reprex_2.0.1                httr_1.4.2                 
# [13] pillar_1.6.3                zlibbioc_1.36.0            
# [15] rlang_0.4.11                readxl_1.3.1               
# [17] rstudioapi_0.13             Matrix_1.3-0               
# [19] BiocParallel_1.24.1         RCurl_1.98-1.3             
# [21] munsell_0.5.0               DelayedArray_0.16.3        
# [23] broom_0.7.9                 compiler_4.0.3             
# [25] modelr_0.1.8                pkgconfig_2.0.3            
# [27] SummarizedExperiment_1.20.0 tidyselect_1.1.1           
# [29] GenomeInfoDbData_1.2.4      matrixStats_0.61.0         
# [31] XML_3.99-0.6                fansi_0.4.1                
# [33] crayon_1.4.1                tzdb_0.1.2                 
# [35] dbplyr_2.1.1                withr_2.4.2                
# [37] GenomicAlignments_1.26.0    bitops_1.0-7               
# [39] grid_4.0.3                  jsonlite_1.7.2             
# [41] gtable_0.3.0                lifecycle_1.0.1            
# [43] DBI_1.1.1                   magrittr_2.0.1             
# [45] scales_1.1.1                cli_3.0.1                  
# [47] stringi_1.5.3               XVector_0.30.0             
# [49] fs_1.5.0                    xml2_1.3.2                 
# [51] ellipsis_0.3.2              generics_0.1.0             
# [53] vctrs_0.3.8                 tools_4.0.3                
# [55] Biobase_2.50.0              glue_1.6.2                 
# [57] MatrixGenerics_1.2.1        hms_1.1.1                  
# [59] colorspace_2.0-0            rvest_1.0.1                
# [61] haven_2.4.1                

