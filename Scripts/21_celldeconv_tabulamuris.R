# ==============================================================================
# 21_celldeconv_tabulamuris.R
# ==============================================================================



# library versions
# (bulk --> cell type deconvolution not included in final manuscript
# due to high between deconv method variation)
# ------------------------------------------------------------------------------

library(tidyverse)

library(SCDC) # SCDC_0.0.0.9000
library(MuSiC) # MuSiC_0.1.1
library(xbioc) # xbioc_0.1.19
library(BisqueRNA) # BisqueRNA_1.0.5 

# & CIBERSORTx: access date March 15th 2021


# load reference single cell data
# ------------------------------------------------------------------------------

# https://github.com/czbiohub-sf/tabula-muris/tree/master
# Tabula muris download

# male mice marked with "_M" in mouse id,
# exclude female mice because they have hepatocytes only

tabula_muris_liver_counts <-
  fread("./Data/tabula_muris_5829687/Liver-counts.csv") %>%
  dplyr::select(X1, contains("_M"))
tabula_muris_liver_metadata <-
  tibble(cell = names(tabula_muris_liver_counts)[-1]) %>%
  left_join(
    read_csv("./Data/tabula_muris_5829687/annotations_facs.csv") %>%
      transmute(cell, celltype = cell_ontology_class, cluster.ids)
  ) %>%
  separate(
    col = "cell", sep = "\\.",
    into = c("well", "plate.barcode", "mouse.id"), remove = F
  ) %>%
  mutate(mouse.sex = if_else(grepl("_M", mouse.id), "M", "F"))

filter_no_cell_assign <-
  !is.na(tabula_muris_liver_metadata$celltype)

tabula_muris_liver_counts <-
  tabula_muris_liver_counts[, c(T, filter_no_cell_assign)]
tabula_muris_liver_metadata <-
  tabula_muris_liver_metadata %>%
  filter(!is.na(celltype))

tibble(sex = tiss$mouse.sex, cell = tiss$cell_ontology_class) %>%
  group_by(sex, cell) %>%
  tally()

refcelltypes <-
  tabula_muris_liver_metadata$celltype %>%
  unique() %>%
  sort()


# prep our normalized bulk RNA-seq in right format &&
# CIBERSORTx (upload to cibersortx.stanford.edu server)
# ------------------------------------------------------------------------------

symbol_in_order <-
  rownames(vst) %>%
  as_tibble() %>%
  left_join(gene_annotations_short %>% dplyr::select(1:2),
            by = c("value" = "ENSG")) %>% 
  .$symbol

bulk <-
  deseqobj %>%
  assay() %>%
  `colnames<-`(gsub("Sample", "Mouse", colnames(vst))) %>%
  `rownames<-`(rownames(vst)) %>%
  as_tibble(rownames = "ENSG") %>%
  bind_cols(Symbol = symbol_in_order, .) %>%
  arrange(Symbol) %>%
  dplyr::select(-ENSG) %>%
  filter(!is.na(Symbol)) %>%
  filter(!duplicated(Symbol))


bulkdata <-
  bulk[, -1] %>%
  as.matrix() %>%
  `rownames<-`(., bulk$Symbol) %>%
  getESET(., pdata = metadata, fdata = bulk$Symbol)

singlecell_set <-
  getESET(
    exprs = tabula_muris_liver_counts[, -1] %>%
      as.matrix() %>%
      `rownames<-`(tabula_muris_liver_counts$X1),
    pdata = data.frame(
      celltypes = tabula_muris_liver_metadata$celltype,
      subjects = tabula_muris_liver_metadata$mouse.id
    ),
    fdata = tabula_muris_liver_counts$X1
  )

exprs(singlecell_set) %>%
  as_tibble(rownames = "Symbol") %>%
  `colnames<-`(., c("Symbol", singlecell_set$celltypes)) %>%
  write_tsv("./Output/20210319_celldeconv/cibersort_scrna_ref.txt")

exprs(bulkdata) %>%
  as_tibble(rownames = "Symbol") %>%
  write_tsv("./Output/20210319_celldeconv/cibersort_realdat_mixture_COUNTS.txt")




# process CIBERSORT 
# ------------------------------------------------------------------------------

metadata_celldeconv_cibersortx <-
  metadata %>%
  left_join(read_csv("./Data/20210319_celldeconv/CIBERSORTx_Results.csv") %>%
              pivot_longer(cols = 2:ncol(.), names_to = "Celltype") %>%
              mutate(ID = gsub("Mouse", "", Mixture)),
            by = "ID"
  ) %>%
  filter(Celltype %in% refcelltypes)

cellcomp_plot_cibersortx <-
  ggplot(
    metadata_celldeconv_cibersortx %>% filter(!is.na(Celltype)),
    aes(Group, value * 100, fill = Group)
  ) +
  geom_quasirandom(width = 0.2, shape = 21, size = 2) +
  facet_grid(Celltype ~ .,
             scales = "free", switch = "y",
             labeller = label_wrap_gen()
  ) +
  scale_x_discrete(labels = c("C", "D", "D/A")) +
  scale_y_continuous("Percent", breaks = pretty_breaks(n = 5)) +
  scale_fill_manual(values = palette_color_group) +
  theme_few() +
  theme(legend.position = "none")



# SCDC
# ------------------------------------------------------------------------------

singlecellqc <-
  SCDC_qc(singlecell_set, ct.sub = refcelltypes, ct.varname = "celltypes")
singlecellqc$prop.qc %>% head

scdc_results <-
  SCDC_prop(
    sc.eset = singlecellqc$sc.eset.qc,
    bulk.eset = bulkdata,
    ct.sub = refcelltypes,
    ct.varname = "celltypes",
    Transform_bisque = F
  )

colnames(scdc_results$basis.mvw)
which.min(scdc_results$basis.mvw[, 2])

SCDC_out <-
  scdc_results$prop.est.mvw %>%
  as_tibble(rownames = "ID") %>%
  pivot_longer(cols = 2:ncol(.), names_to = "Celltype") %>%
  mutate(ID = gsub("Mouse", "", ID))

scdc_results$prop.est.mvw[, "hepatocyte"] %>% hist()

# 400 x 700
cellcomp_plot_scdc <-
  ggplot(
    metadata_celldeconv_SCDC_std,
    aes(Group, value * 100, fill = Group)
  ) +
  geom_quasirandom(width = 0.2, shape = 21, size = 2) +
  facet_grid(Celltype ~ .,
    scales = "free", switch = "y",
    labeller = label_wrap_gen()
  ) +
  scale_x_discrete(labels = c("C", "D", "D/A")) +
  scale_y_continuous("Percent", breaks = pretty_breaks(n = 5)) +
  scale_fill_manual(values = palette_color_group) +
  theme_few() +
  theme(legend.position = "none")




# MuSiC
# ------------------------------------------------------------------------------

music_results <-
  music_prop(
    bulk.eset = bulkdata,
    sc.eset = singlecell_set,
    # sc.eset = singlecellqc,
    clusters = "celltypes", samples = "subjects",
    select.ct = refcelltypes
  )

metadata_celldeconv_music <-
  metadata %>%
  left_join(music_results$Est.prop.weighted %>% as_tibble(rownames = "ID") %>%
    pivot_longer(cols = 2:ncol(.), names_to = "Celltype") %>%
    mutate(ID = gsub("Mouse", "", ID)),
  by = "ID"
  )

cellcomp_plot_music <-
  ggplot(
    metadata_celldeconv_music %>% filter(!is.na(Celltype)),
    aes(Group, value * 100, fill = Group)
  ) +
  geom_quasirandom(width = 0.2, shape = 21, size = 2) +
  facet_grid(Celltype ~ .,
    scales = "free", switch = "y",
    labeller = label_wrap_gen()
  ) +
  scale_y_continuous("Percent", breaks = pretty_breaks(n = 5)) +
  scale_fill_manual(values = palette_color_group) +
  scale_x_discrete(labels = c("C", "D", "D/A")) +
  theme_few() +
  theme(legend.position = "none")



# Bisque
# ------------------------------------------------------------------------------

bisque_results <-
  ReferenceBasedDecomposition(
    bulk.eset = bulkdata,
    sc.eset = singlecell_set,
    # sc.eset = singlecellqc,
    cell.types = "celltypes", subject.names = "subjects",
    use.overlap = F
  )

bisque_results$sc.props %>% head()

metadata_celldeconv_bisque <-
  metadata %>%
  left_join(bisque_results$bulk.props %>% t() %>% as_tibble(rownames = "ID") %>%
    pivot_longer(cols = 2:ncol(.), names_to = "Celltype") %>%
    mutate(ID = gsub("Mouse", "", ID)),
  by = "ID"
  )

cellcomp_plot_bisque <-
  ggplot(
    metadata_celldeconv_bisque %>% filter(!is.na(Celltype)),
    aes(Group, value * 100, fill = Group)
  ) +
  geom_quasirandom(width = 0.2, shape = 21, size = 2) +
  facet_grid(Celltype ~ .,
    scales = "free", switch = "y",
    labeller = label_wrap_gen()
  ) +
  scale_x_discrete(labels = c("C", "D", "D/A")) +
  scale_y_continuous("Percent", breaks = pretty_breaks(n = 5)) +
  scale_fill_manual(values = palette_color_group) +
  theme_few() +
  theme(legend.position = "none")





# across deconvolution method comparison
# calculate signif proportions across groups
# ------------------------------------------------------------------------------

# 800 x 600
plot_grid(cellcomp_plot_bisque + ggtitle("Bisque"),
  cellcomp_plot_cibersortx + ggtitle("CIBERSORTx"),
  cellcomp_plot_music + ggtitle("MuSiC"),
  cellcomp_plot_scdc + ggtitle("SCDC"),
  nrow = 1, labels = LETTERS[1:4]
)

generate_pvals_cellcomp <-
  function(cellcomp_df) {
    list(
      cellcomp_df,
      cellcomp_df %>% filter(Group %in% c("Chow", "Disease")),
      cellcomp_df %>% filter(Group %in% c("Chow", "Disease/Agonist")),
      cellcomp_df %>% filter(Group %in% c("Disease/Agonist", "Disease"))
    ) %>%
      map_dbl(~ kruskal.test(value ~ Group, .x)$p.value)
  }


pvals_celltype <-
  metadata_celldeconv_music %>% # use MuSiC based on simulation
  # metadata_celldeconv_bisque %>%
  # metadata_celldeconv_SCDC_std %>%
  # metadata_celldeconv_cibersortx %>%
  group_by(Celltype) %>%
  group_split() %>%
  map(~ generate_pvals_cellcomp(.x)) %>%
  set_names(., metadata_celldeconv_SCDC_std$Celltype %>% unique()) %>%
  as_tibble() %>%
  bind_cols(Test = c("Omnibus", "DvC", "DAvC", "DAvD"), .)

pvals_celltype_omnibus <-
  metadata_celldeconv_music %>%
  group_by(Celltype) %>%
  group_split() %>%
  map_dbl(~ aov(value ~ Group, .x) %>%
    summary() %>%
    unlist() %>%
    .["Pr(>F)1"])
