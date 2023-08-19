# ==============================================================================
# 16_heatmaps.R
# ==============================================================================



# anything signif in pairwise comparisons (union)
# ------------------------------------------------------------------------------

ensg_signif_in_any <-
  list(DEresults_RUVr_DAvC,
       DEresults_RUVr_DAvD,
       DEresults_RUVr_DvC) %>%
  map(~ filter(.x, padj < 0.05) %>%
        .$row) %>%
  unlist %>%
  unique

# could also add add "& log2FoldChange > 2" to filter()



# heatmap pre-RUVr correction
# 0-to-1 scaling
# ------------------------------------------------------------------------------

heatmap_mat_RUVr <-
  ruv_results$normalizedCounts %>% # alternativly can use vst-RUVr counts
  .[heatmap_list$row,
  ] %>%
  t %>%
  apply(., 2, function(x) { (x - min(x))/(max(x) - min(x))}) %>%
  `rownames<-`(NULL) %>% 
  `colnames<-`(NULL) %>%
  t

ruvr_heatmap <-
  Heatmap(heatmap_mat_RUVr,
          name = "Relative\nExpression",
          col = viridis(200, option = "D"),
          column_split = metadata$Group %>%
            fct_recode(Chow = "Chow", `DSS-PN` = "Disease", `DSS-PN/DLPC` = "Disease/Agonist"),
          show_column_dend = F, cluster_column_slices = F,
          column_title_gp = gpar(fontsize = 10),
          row_names_gp = gpar(fontsize = 8),
          column_names_gp = gpar(fontsize = 8),
          heatmap_legend_param = list(title_gp = gpar(fontsize = 9),
                                      labels_gp = gpar(fontsize = 8),
                                      at = c(0, 0.25, 0.5, 0.75, 1),
                                      labels = c("0 (low)", "0.25", "0.5", "0.75", "1 (high)"),
                                      legend_height = unit(2, "cm"))
  ) 


ggsave(filename = "Figures/Fig1A.png",
       plot = grid.grabExpr(draw(ruvr_heatmap)),
       width = 5, height = 4, dpi = 1200)

system('open "Figures/Fig1A.png"')







# M1/M2 polarization maker list from S.G.
# via email subject "The meeting sept 10 2021"; Sun, Oct 31, 2021, 3:38 PM
# ------------------------------------------------------------------------------

m1m2_list <-
  read_excel("./Data/list m1 m2 - CL - version2 - SG.xlsx") %>%
  mutate(ENSG_ID = if_else(!ENSG_ID %in% rownames(vst_RUVr),
                           "notpresent", ENSG_ID))
m1m2_list <-
  m1m2_list %>%
  filter(ENSG_ID %in% rownames(vst_RUVr)) %>%
  filter(!duplicated(ENSG_ID))

# consider using only M1/M2 DEGs in any set 
m1m2_list <-
  m1m2_list %>%
  filter(ENSG_ID %in% ensg_signif_in_any)

heatmap_mat_RUVr <-
  ruv_results$normalizedCounts %>%
  .[m1m2_list$ENSG_ID,
  ] %>%
  t %>%
  apply(., 2, function(x) { (x - min(x))/(max(x) - min(x))}) %>%
  `rownames<-`(NULL) %>% 
  `colnames<-`(m1m2_list$Manuscript_Name) %>%
  t

ruvr_heatmap <-
  Heatmap(heatmap_mat_RUVr,
          name = "Relative\nExpression",
          col = viridis(200, option = "D"),
          right_annotation =
            rowAnnotation(df = m1m2_list %>% transmute(`M1/M2` = M1M2) %>% as.data.frame(),
                          col = list(`M1/M2` = palette_color_M1M2),
                          annotation_legend_param =
                            list(title_gp = gpar(fontsize = 9),
                                 labels_gp = gpar(fontsize = 8)),
                          annotation_name_gp = gpar(fontsize = 9)),
          column_split = metadata$Group %>%
            fct_recode(Chow = "Chow", `DSS-PN` = "Disease", `DSS-PN/DLPC` = "Disease/Agonist"),
          show_column_dend = F, cluster_column_slices = F,
          column_title_gp = gpar(fontsize = 10),
          row_names_gp = gpar(fontsize = 8),
          column_names_gp = gpar(fontsize = 8),
          heatmap_legend_param = list(title_gp = gpar(fontsize = 9),
                                      labels_gp = gpar(fontsize = 8),
                                      at = c(0, 0.25, 0.5, 0.75, 1),
                                      labels = c("0 (low)", "0.25", "0.5", "0.75", "1 (high)"),
                                      legend_height = unit(2, "cm"))
  ) 

ggsave(filename = "Figures/Fig3A.png",
       plot = grid.grabExpr(draw(ruvr_heatmap)),
       width = 5, height = 7, dpi = 1200)

system('open "Figures/Fig3A.png"')