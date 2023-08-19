# ==============================================================================
# 18_webgestalt_gsea.R
# ==============================================================================



# for upload to WebGestalt server
dir.create("./Output/20210315_webgestalt")

# WebGestalt ORA 
# ------------------------------------------------------------------------------

writeLines(DEresults_RUVr_DvC %>% filter(!is.na(pvalue)) %>% .$symbol,
           "./Output/20210315_webgestalt/CvD_ref.txt")
writeLines(DEresults_RUVr_DAvC %>% filter(!is.na(pvalue)) %>% .$symbol,
           "./Output/20210315_webgestalt/CvDA_ref.txt")
writeLines(DEresults_RUVr_DAvD %>% filter(!is.na(pvalue)) %>% .$symbol,
           "./Output/20210315_webgestalt/DvDA_ref.txt")

writeLines(DEresults_RUVr_DvC %>% filter(padj < 0.05 & !is.na(pvalue)) %>% .$symbol,
           "./Output/20210315_webgestalt/CvD_RUVr_target.txt")
writeLines(DEresults_RUVr_DAvC %>% filter(padj < 0.05 & !is.na(pvalue)) %>% .$symbol,
           "./Output/20210315_webgestalt/CvDA_RUVr_target.txt")
writeLines(DEresults_RUVr_DAvD %>% filter(padj < 0.05 & !is.na(pvalue)) %>% .$symbol,
           "./Output/20210315_webgestalt/DvDA_RUVr_target.txt")




# WebGestalt GSEA
# ------------------------------------------------------------------------------

write_tsv(DEresults_RUVr_DvC %>% filter(!is.na(pvalue)) %>% select(symbol, log2FoldChange),
          "./Output/20210315_webgestalt/DvC_RUVr_gsea.rnk", col_names = F)
write_tsv(DEresults_RUVr_DAvC %>% filter(!is.na(pvalue)) %>% select(symbol, log2FoldChange),
          "./Output/20210315_webgestalt/DAvC_RUVr_gsea.rnk", col_names = F)
write_tsv(DEresults_RUVr_DAvD %>% filter(!is.na(pvalue)) %>% select(symbol, log2FoldChange),
          "./Output/20210315_webgestalt/DAvD_RUVr_gsea.rnk", col_names = F)




