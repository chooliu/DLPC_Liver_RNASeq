# ==============================================================================
# 11_analysis_setup.R
# ==============================================================================



# libraries
# ------------------------------------------------------------------------------

library(tidyverse)
library(readxl)
library(lubridate)

library(DESeq2)
library(RUVSeq)
library(vegan)
library(WebGestaltR)

library(scales)
library(ggbeeswarm)
library(ggrepel)
library(ggthemes)
library(ComplexHeatmap)
library(viridis)
library(eulerr)
library(colorspace)
library(cowplot)

library(writexl)

# save version control history for scripts 11 - 30
# (scripts beginning with 0*_ run on linux server,
#  version info on end of each 0*_ script)
writeLines(capture.output(sessionInfo()), "R_sessionInfo.txt")



# palettes
# ------------------------------------------------------------------------------

palette_shape_group <-
  c(Chow = 1, Disease = 19, `Disease/Agonist` = 15,
    `DSS-PN` = 19, `DSS-PN/DLPC` = 15)

palette_color_group <-
  c(Chow = "black", Disease = "#f4a582", `Disease/Agonist` = "#b2182b",
    `DSS-PN` = "#f4a582", `DSS-PN/DLPC` = "#b2182b")

palette_pubumod <- 
  c('#8c510a','#bf812d','#dfc27d','#f6e8c3','#f5f5f5',
    '#c7eae5','#80cdc1','#35978f','#01665e') %>%
  .[-(4:6)] %>%
  set_names(metadata$Batch %>% unique)

palette_color_signif <-
  c(yes = "#2b8cbe", `TRUE` = "#2b8cbe",
    no = "grey", `FALSE` = "grey")

palette_color_signif_darker <-
  c(yes = "#0b628f", `TRUE` = "#0b628f",
    no = "#b3b3b3", `FALSE` = "#b3b3b3")

palette_shape_signif <-
  c(yes = 19, `TRUE` = 19,
    no = 1, `FALSE` = 1)

palette_color_M1M2 <-
  c(M1 = "#b2182b", M2 = "#e0e0e0")



