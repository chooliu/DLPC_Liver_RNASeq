# ==============================================================================
# 30_export_analysis.R
# ==============================================================================



# final exports after rest of analysis is run
# exporting RUVr-normalized counts for visualization
# ------------------------------------------------------------------------------

symbol_in_order <-
  rownames(vst) %>%
  as_tibble() %>%
  left_join(gene_annotations_short %>% select(1:2), by = c("value" = "ENSG")) %>% 
  .$symbol
order(symbol_in_order)

assay(vst) %>%
  `colnames<-`(gsub("Sample", "Mouse", colnames(.))) %>%
  as_tibble(rownames = "ENSG_ID") %>%
  bind_cols(Symbol = symbol_in_order, .) %>%
  arrange(Symbol) %>%
  write_tsv("./Output/20210325_PNALI_Results_v2/rlog_expression_matrix.tsv")

ruv_results$normalizedCounts %>%
  `colnames<-`(gsub("Sample", "Mouse", colnames(.))) %>%
  as_tibble(rownames = "ENSG_ID") %>%
  bind_cols(Symbol = symbol_in_order, .) %>%
  arrange(Symbol) %>%
  write_tsv("./Output/20210325_PNALI_Results_v2/RUVcorrected_expression_matrix.tsv")




# helper fxns to export results
# DESeq2 output marked "v2" because of mislabeled sample metadata in "v1"
# ------------------------------------------------------------------------------

prep_DEresults <- function(resultsdat, higherin, altin) {
  resultsdat %>%
    transmute(ENSG = row, Symbol = symbol,
              `Gene Type` = gene_type,
              log2FoldChange,
              `Fold Change` = 2^log2FoldChange,
              HigherInThisGroup =
                if_else(log2FoldChange > 0, higherin, altin),
              `p-value` = pvalue, FDR = padj,
              Signif = if_else(FDR < 0.05, "*", ""))
}


prep_gsea <- function(gsea1, gsea2) {
    tmp <-
      bind_rows(
      read_tsv(gsea1),
      read_tsv(gsea2)) %>%
      arrange(pValue) %>%
      mutate(pValue = (pValue*1000 + 1)/1001)
    
    tmp$FDR <- p.adjust(tmp$pValue, "fdr")
    
    tmp %>%
      transmute(`Gene Set` = geneSet,
                Description = description,
                `Normalized Enrichment Score` = normalizedEnrichmentScore,
                `p-value` = pValue,
                FDR = FDR,
                `# Genes in Set` = size,
                `# Leading Edge Genes` = leadingEdgeNum,
                `Gene Symbols` = userId)
}


prep_ora <- function(filepath) {
  
  read_tsv(filepath) %>%
    transmute(`Gene Set` = geneSet,
              Description = description,
              `p-value` = pValue,
              FDR = FDR,
              `# Genes in Set` = size,
              `# Overlapping Signif Genes` = overlap,
              `Gene Symbols` = userId)
}



# export excel document
# ------------------------------------------------------------------------------

write_xlsx(
  list(DiffExp_DvC = DEresults_DvC %>% prep_DEresults(., "D", "C"),
       DiffExp_DAvC = DEresults_DAvC %>% prep_DEresults(., "D/A", "C"),
       DiffExp_DAvD = DEresults_DAvD %>% prep_DEresults(., "D/A", "D"),
       DiffExp_DvC_RUV = DEresults_RUVr_DvC %>% prep_DEresults(., "D", "C"),
       DiffExp_DAvC_RUV = DEresults_RUVr_DAvC %>% prep_DEresults(., "D/A", "C"),
       DiffExp_DAvD_RUV = DEresults_RUVr_DAvD %>% prep_DEresults(., "D/A", "D"),
                
       ORA_DvC = prep_ora("./Output/20210315_webgestalt/Project_CvD/enrichment_results_CvD.txt"),
       ORA_DAvD = prep_ora("./Output/20210315_webgestalt/Project_DvDA/enrichment_results_DvDA.txt"),
       ORA_DAvC = prep_ora("./Output/20210315_webgestalt/Project_CvDA/enrichment_results_CvDA.txt"),
       ORA_DvC_RUV = prep_ora("./Output/20210315_webgestalt/Project_CvD_RUVr/enrichment_results_CvD_RUVr.txt"),
       ORA_DAvD_RUV = prep_ora("./Output/20210315_webgestalt/Project_DvDA_RUVr/enrichment_results_DvDA_RUVr.txt"),
       ORA_DAvC_RUV = prep_ora("./Output/20210315_webgestalt/Project_CvDA_RUVr/enrichment_results_CvDA_RUVr.txt"),
       
       GSEA_DvC = prep_gsea("./Output/20210315_webgestalt/Project_wg_DvC_Reactome/enrichment_results_wg_result1615848429.txt",
                            "./Output/20210315_webgestalt/Project_wg_DvC_GO/enrichment_results_wg_result1615848443.txt"),
       GSEA_DAvC = prep_gsea("./Output/20210315_webgestalt/Project_wg_DAvC_Reactome/enrichment_results_wg_result1615849237.txt",
                             "./Output/20210315_webgestalt/Project_wg_DAvC_GO/enrichment_results_wg_result1615849227.txt"),
       GSEA_DAvD = prep_gsea("./Output/20210315_webgestalt/Project_wg_DAvD_Reactome/enrichment_results_wg_result1615850686.txt",
                             "./Output/20210315_webgestalt/Project_wg_DAvD_GO/enrichment_results_wg_result1615850694.txt"),
       GSEA_DvC_RUV = prep_gsea("./Output/20210315_webgestalt/Project_wg_DvC_RUVr_Reactome/enrichment_results_wg_result1616111565.txt",
                                "./Output/20210315_webgestalt/Project_wg_DvC_RUVr_GO/enrichment_results_wg_result1616111557.txt"),
       GSEA_DAvC_RUV = prep_gsea("./Output/20210315_webgestalt/Project_wg_DAvC_RUVr_Reactome/enrichment_results_wg_result1616111283.txt",
                                 "./Output/20210315_webgestalt/Project_wg_DAvC_RUVr_GO/enrichment_results_wg_result1616112609.txt"),
       GSEA_DAvD_RUV = prep_gsea("./Output/20210315_webgestalt/Project_wg_DAvD_RUVr_Reactome/enrichment_results_wg_result1616111531.txt",
                                 "./Output/20210315_webgestalt/Project_wg_DAvD_RUVr_GO/enrichment_results_wg_result1616111539.txt"),
                
       ChEA_DvC = export_chea3(chea3_DvC),
       ChEA_DAvC = export_chea3(chea3_DAvC),
       ChEA_DAvD = export_chea3(chea3_DAvD),
                
       ChEA_DvC_RUV = export_chea3(chea3_DvC_RUVr),
       ChEA_DAvC_RUV = export_chea3(chea3_DAvC_RUVr),
       ChEA_DAvD_RUV = export_chea3(chea3_DAvD_RUVr),
       MuSiC_Celltype = pvals_celltype) %>%
    set_names(., paste0(1:length(.), "_", names(.))),
       
  "Output/20210315_PNALI_Results_v2.xlsx")


system('open "Output/20210315_PNALI_Results_v2.xlsx"')


# cell composition estimates
# ------------------------------------------------------------------------------

# or alt. algorithms
# write_tsv(bisque_results$bulk.props %>% t %>% as_tibble(rownames = "Mouse"),

write_tsv(
  music_results$Est.prop.weighted %>% as_tibble(rownames = "Mouse"),
  "Output/20210325_PNALI_Results_v2/20210325_PNALI_celldeconv_MuSiC.tsv")



# minimal metadata for reproducing analysis 
# note that "Sample724" in count matrix has correct ID = "Sample725"
# ------------------------------------------------------------------------------

metadata %>%
  transmute(Sample = paste0("Sample", ID),
              Group) %>%
  write_tsv("Output/metadata.tsv")
