# ==============================================================================
# 14_deseq2_RUVr.R
# ==============================================================================



# vst-transform naive soln --> residuals for RUVr 
# then check k for RUVr
# ------------------------------------------------------------------------------

resid_vals <-
  assay(vst) %>%
  apply(., 1, function(y) {
    lm(y ~ Group, metadata) %>% resid
  })

ruv_out <-
  sapply(1:10,
       function(ktry) {
ruv_results <-
  RUVr(x = assay(vst),
     k = ktry,
     residuals = resid_vals %>% t,
     isLog = T)

dm_test <-
  assay(vst) %>%
  .[filter_genes_variance, ] %>%
  t %>% dist

  adonis(dm_test ~ ruv_results$W + metadata$Group) %>%
    .$aov.tab %>% .$R2 %>% .[1]
       }
)


plot(1 - ruv_out)

plot(ruv_out)
ggplot(metadata,
       aes(Batch, ruv_results$W[ , 1], color = Group)) +
  geom_boxplot() +
  geom_quasirandom(width = 0.2) +
  theme_few() +
  scale_color_manual(values = palette_color_group)



# final RUVr solution
# ------------------------------------------------------------------------------

ruv_results <-
  RUVr(x = assay(vst),
       k = 4,
       residuals = resid_vals %>% t,
       isLog = T)
adonis(dm_test ~ ruv_results$W + metadata$Group) %>% .$aov.tab %>% .$R2 %>% .[1]

pca_RUVr <-
  ruv_results$normalizedCounts %>%
  .[filter_genes_variance, ] %>%
  t %>%
  prcomp(scale = T)

percentexp_RUVr <-
  formatC(pca_RUVr$sdev^2 / sum(pca_RUVr$sdev^2) * 100, digits = 1, format = "f") %>%
  paste0("PC", 1:length(.), " (", ., "%)")
metadata_RUVr <-
  metadata %>% 
  bind_cols(PC1 = pca_RUVr$x[ , 1],
            PC2 = pca_RUVr$x[ , 2]) %>%
  bind_cols(RUV1 = ruv_results$W[ , 1],
            RUV2 = ruv_results$W[ , 2],
            RUV3 = ruv_results$W[ , 3],
            RUV4 = ruv_results$W[ , 4])



# check RUV-corrected PCA
# color by GROUP
# ------------------------------------------------------------------------------

pca_by_group <-
  ggplot(metadata_RUVr %>%
           mutate(Group = fct_recode(
             Group, `DSS-PN` = "Disease", `DSS-PN/DLPC` = "Disease/Agonist")),
         aes(PC1, PC2,
             color = Group, shape = Group)) +
  geom_point(size = 4, alpha = 0.75) +
  # geom_text_repel(aes(label = ID), show.legend = F, seed = 123) +
  theme_few() + theme(legend.position = "bottom") +
  xlab(percentexp_RUVr[1]) + ylab(percentexp_RUVr[2]) +
  geom_point(alpha = 0.25) +
  scale_color_manual("Group", values = palette_color_group %>% .[c(1, 4, 5)]) +
  scale_shape_manual("Group", values = palette_shape_group %>% .[c(1, 4, 5)]) +
  NULL

ggsave(pca_by_group,
       filename = "Figures/Fig1B.png",
       height = 5, width = 5, dpi = 1200)

# color by DATE
pca_by_date <-
  ggplot(metadata_RUVr,
         aes(PC1, PC2,
             color = Batch,
             shape = Group)) +
  geom_point(size = 4, alpha = 0.75) +
  geom_text_repel(aes(label = ID), show.legend = F, seed = 123) +
  theme_few() + theme(legend.position = "right") +
  xlab(percentexp[1]) + ylab(percentexp[2]) +
  geom_point(alpha = 0.25) +
  ggConvexHull::geom_convexhull(aes(fill = Batch), alpha = 0.5) +
  scale_fill_manual(values = palette_pubumod) +
  scale_color_manual(values = palette_pubumod) +
  # scale_color_manual("Group", values = palette_color_group) +
  scale_shape_manual("Group", values = palette_shape_group) +
  NULL

# 900 x 400
plot_grid(pca_by_group,
          pca_by_date,
          rel_widths = c(0.4, 0.6),
          labels = c("A", "B"))



# re-run DESeq2 for final results
# ------------------------------------------------------------------------------

deseqobj_RUVr <-
  DESeqDataSetFromMatrix(countData = filtered_counts_mat,
                         colData = metadata_RUVr,
                         design = ~ Group + RUV1 + RUV2 + RUV3 + RUV4) %>%
  DESeq()

vst_RUVr <-
  rlog(deseqobj_RUVr)


DEresults_RUVr_DvC <-
  results(deseqobj_RUVr, contrast = c("Group", "Disease", "Chow"),
          tidy = T, independentFiltering = F, cooksCutoff = F) %>%
  arrange(padj) %>%
  left_join(gene_annotations_short, by = c("row" = "ENSG"))

DEresults_RUVr_DAvC <-
  results(deseqobj_RUVr, contrast = c("Group", "Disease/Agonist", "Chow"),
          tidy = T, independentFiltering = F, cooksCutoff = F) %>%
  arrange(padj) %>%
  left_join(gene_annotations_short, by = c("row" = "ENSG"))

DEresults_RUVr_DAvD <-
  results(deseqobj_RUVr, contrast = c("Group", "Disease/Agonist", "Disease"),
          tidy = T, independentFiltering = F, cooksCutoff = F) %>%
  arrange(padj) %>%
  left_join(gene_annotations_short, by = c("row" = "ENSG"))



# viz DESeq2 results
# ------------------------------------------------------------------------------

DEresults_RUVr_DAvC$pvalue %>% hist
DEresults_RUVr_DvC$pvalue %>% hist
DEresults_RUVr_DAvD$pvalue %>% hist

# 500 x 350
list(`DAvC\n(with RUVr)` = DEresults_RUVr_DAvC %>% filter(padj < 0.05) %>% .$row,
     `DAvD\n(with RUVr)` = DEresults_RUVr_DAvD %>% filter(padj < 0.05) %>% .$row,
     `DvC (with RUVr)` = DEresults_RUVr_DvC %>% filter(padj < 0.05) %>% .$row ) %>%
  euler %>%
  plot(quantities = T)

# check +/- RUVr correction
list(`DvC\n(with RUVr)` = DEresults_RUVr_DvC %>% filter(padj < 0.05) %>% .$row,
     # DvC_RUVr = DEresults_RUVr_DvC %>% filter(padj < 0.05) %>% .$row,
     `DvC\n(no batch adj)` = DEresults_DvC %>% filter(padj < 0.05) %>% .$row#,
     # DvC = DEresults_DvC %>% filter(padj < 0.05) %>% .$row
     ) %>% # map_chr(length)
  euler %>%
  plot(quantities = T, fills = c("#fee8c8", "grey"))



# check effect size correlations
# ------------------------------------------------------------------------------

# compare DAvD
# 500 x 400

comparison_RUVr_plusminus <-
  left_join(DEresults_DAvD, DEresults_RUVr_DAvD, by = "row",
            suffix = c(".naive", ".ruvr"))

comparison_RUVr_tmp <-
  comparison_RUVr_plusminus %>%
  filter((log2FoldChange.naive * log2FoldChange.ruvr < 0) &
           (padj.ruvr < 0.05 | padj.naive < 0.05))

ggplot(comparison_RUVr_plusminus,
       aes(log2FoldChange.naive, log2FoldChange.ruvr,
           shape = padj.ruvr < 0.05,
           color = padj.ruvr < 0.05)) +
  geom_point(alpha = 0.5) +
  geom_text_repel(data = comparison_RUVr_tmp,
                  aes(label = symbol.ruvr),
                  color = "black", box.padding = 0.1) +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + 
  
  theme_few() +
  scale_color_manual(values = palette_color_signif) +
  scale_shape_manual(values = palette_shape_signif) +
  theme(legend.position = "none") +
  scale_x_continuous(expression(log[2]~"(FC), Naive"), breaks = pretty_breaks()) +
  scale_y_continuous(expression(log[2]~"(FC), RUVr-corrected"), breaks = pretty_breaks()) +
    ggtitle("DAvD Effects (Naive vs. RUVr)")


# compare DvC
comparison_RUVr_plusminus <-
  left_join(DEresults_DvC, DEresults_RUVr_DvC, by = "row",
            suffix = c(".naive", ".ruvr"))

comparison_RUVr_tmp <-
  comparison_RUVr_plusminus %>%
  filter((log2FoldChange.naive * log2FoldChange.ruvr < 0) &
           (padj.ruvr < 0.05 | padj.naive < 0.05))

ggplot(comparison_RUVr_plusminus,
       aes(log2FoldChange.naive, log2FoldChange.ruvr,
           shape = padj.ruvr < 0.05,
           color = padj.ruvr < 0.05)) +
  geom_point(alpha = 0.5) +
  geom_text_repel(data = comparison_RUVr_tmp,
                  aes(label = symbol.ruvr),
                  color = "black", box.padding = 0.1) +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + 
  theme_few() +
  scale_color_manual(values = palette_color_signif) +
  scale_shape_manual(values = palette_shape_signif) +
  theme(legend.position = "none") +
  scale_x_continuous(expression(log[2]~"(FC), Naive"), breaks = pretty_breaks()) +
  scale_y_continuous(expression(log[2]~"(FC), RUVr-corrected"), breaks = pretty_breaks()) +
  ggtitle("DvC Effects (Naive vs. RUVr)")



# compare DAvC
make_comparison_plot <-
  function(df1, df2, plottitle) {
  
  comparison_RUVr_plusminus <-
    left_join(df1, df2, by = "row",
              suffix = c(".naive", ".ruvr"))

comparison_RUVr_tmp <-
  comparison_RUVr_plusminus %>%
  # filter(!is.na(log2FoldChange.ruvr) & !is.na(log2FoldChange.ruvr)) %>%
  filter((log2FoldChange.naive * log2FoldChange.ruvr < 0) &
           (padj.ruvr < 0.05 | padj.naive < 0.05))

corval <-
  cor(comparison_RUVr_plusminus$log2FoldChange.naive,
      comparison_RUVr_plusminus$log2FoldChange.ruvr)

ggplot(comparison_RUVr_plusminus,
       aes(log2FoldChange.naive, log2FoldChange.ruvr,
           shape = padj.ruvr < 0.05,
           color = padj.ruvr < 0.05)) +
  geom_point(alpha = 0.5) +
  geom_text_repel(data = comparison_RUVr_tmp,
                  aes(label = symbol.ruvr),
                  color = "black", box.padding = 0.1) +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + 
  
  theme_few() +
  scale_color_manual(values = palette_color_signif) +
  scale_shape_manual(values = palette_shape_signif) +
  theme(legend.position = "none") +
  scale_x_continuous(expression(log[2]~"(FC), Naive"), breaks = pretty_breaks()) +
  scale_y_continuous(expression(log[2]~"(FC), RUVr-corrected"), breaks = pretty_breaks()) +
  ggtitle(plottitle) +
  annotate(geom = "text",
           label = paste0("Spearman r = ", formatC(corval, digits = 2), "\n"),
           x = Inf, y = -Inf, hjust = 1, vjust = 0)
}


# 1000 x 350
plot_grid(
  make_comparison_plot(DEresults_DvC, DEresults_RUVr_DvC,
                       "DvC Effects (Naive vs. RUVr)"),
  make_comparison_plot(DEresults_DAvC, DEresults_RUVr_DAvC,
                       "DAvC Effects (Naive vs. RUVr)"),
  make_comparison_plot(DEresults_DAvD, DEresults_RUVr_DAvD,
                       "DAvD Effects (Naive vs. RUVr)"),
  nrow = 1, labels = LETTERS[1:3])

