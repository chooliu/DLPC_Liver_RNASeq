# ==============================================================================
# 13_deseq2_naive.R
# ==============================================================================



# filter genes 
# ------------------------------------------------------------------------------

filter_genes_detected <-
  count_mat[ , -1] %>%
  apply(1, function(x) { sum(x == 0) }) %>% # hist(filter_genes_detected)  
  `<=`(., 10)

sample_libsize <-
  count_mat[ , -1] %>%
  colSums %>%
  `*`(1e-6)
filter_genes_cpm <-
  count_mat[ , -1] %>% 
  apply(., 1, function(x) { mean(x/sample_libsize) > 0.25 } )

filter_genes <-
  filter_genes_detected & filter_genes_cpm




# filter samples
# ------------------------------------------------------------------------------

# no filtering
# considered excluding metadata$ID != "724" before found metadata
filter_samples <-
  filter_samples <- rep(T, length(metadata$ID))




# create deseq2 obj
# ------------------------------------------------------------------------------

filtered_counts_mat <-
  count_mat[ , -1] %>%
  as.matrix %>%
  .[filter_genes, ] %>%
  apply(., 2, as.integer) %>%
  `rownames<-`(count_mat$Feature[filter_genes])


deseqobj <-
  DESeqDataSetFromMatrix(countData = filtered_counts_mat,
                         colData = metadata,
                         design = ~ Group) %>%
  DESeq()



# results tables
# ------------------------------------------------------------------------------

DEresults_DvC <-
  results(deseqobj, contrast = c("Group", "Disease", "Chow"),
          tidy = T, independentFiltering = F, cooksCutoff = F) %>%
  arrange(padj) %>%
  left_join(gene_annotations_short, by = c("row" = "ENSG"))

DEresults_DAvC <-
  results(deseqobj, contrast = c("Group", "Disease/Agonist", "Chow"),
          tidy = T, independentFiltering = F, cooksCutoff = F) %>%
  arrange(padj) %>%
  left_join(gene_annotations_short, by = c("row" = "ENSG"))

DEresults_DAvD <-
  results(deseqobj, contrast = c("Group", "Disease/Agonist", "Disease"),
          tidy = T, independentFiltering = F, cooksCutoff = F) %>%
  arrange(padj) %>%
  left_join(gene_annotations_short, by = c("row" = "ENSG"))



# check p-value distribution
DEresults_DAvC$pvalue %>% hist
DEresults_DvC$pvalue %>% hist
DEresults_DAvD$pvalue %>% hist


# check signif across pairwise comparisons
list(DAvC = DEresults_DAvC %>% filter(padj < 0.05) %>% .$row,
     DAvD = DEresults_DAvD %>% filter(padj < 0.05) %>% .$row,
     DvC = DEresults_DvC %>% filter(padj < 0.05) %>% .$row) %>% # map_chr(length)
  euler %>%
  plot(quantities = T)




# normalization for plotting
# ------------------------------------------------------------------------------

vst <-
  rlog(deseqobj)


filter_genes_variance <- 
  vst %>%
  assay %>%
  apply(., 1, mean) %>%
  order(decreasing = T) %>%
  .[1:2500]
pca <-
  assay(vst) %>%
  .[filter_genes_variance, ] %>%
  t %>%
  prcomp(scale = T)
percentexp <-
  formatC(pca$sdev^2 / sum(pca$sdev^2) * 100, digits = 1, format = "f") %>%
  paste0("PC", 1:length(.), " (", ., "%)")
metadata_for_pca <-
  metadata %>% 
  bind_cols(PC1 = pca$x[ , 1],
            PC2 = pca$x[ , 2])


# PCA, color by GROUP
# ------------------------------------------------------------------------------
pca_by_group <-
  ggplot(metadata_for_pca,
       aes(PC1, PC2,
           color = Group, shape = Group)) +
  geom_point(size = 4, alpha = 0.75) +
  geom_text_repel(aes(label = ID), show.legend = F, seed = 123) +
  theme_few() + theme(legend.position = "none") +
  xlab(percentexp[1]) + ylab(percentexp[2]) +
  geom_point(alpha = 0.25) +
  scale_color_manual("Group", values = palette_color_group) +
  scale_shape_manual("Group", values = palette_shape_group) +
  NULL



# PCA, color by DATE
# (main technical effect of suspicion)
# ------------------------------------------------------------------------------

pca_by_date <-
  ggplot(metadata_for_pca,
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
