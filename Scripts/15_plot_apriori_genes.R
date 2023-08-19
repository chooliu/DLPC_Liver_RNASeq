# ==============================================================================
# 15_plot_apriori_genes.R
# ==============================================================================



# checking gene lists previously shown to be DEGs in model
# i.e., "validated" markers via qPCR / other methods
# ------------------------------------------------------------------------------

# list of genes via .xlsx 
# via e-mail subject: "RNA seq" Wed, Mar 10, 2021, 8:22 PM
apriori_genes <-
  read_excel("./Data/RNA seq gene list for now .xlsx") %>%
  .$`gene name` %>%
  .[!is.na(.)] %>% 
  sort

# via email subject: "RNA seq" Fri, Mar 12 2021, 4:59 PM
apriori_genes <-
  c("Beep", "Abcc2", "Abcg5", "abcg8", "il1b", "nr5a2",
    "cyp7a1", "cyp7b", "nr1h4", "nr1h4", "il1rn", "saa2",
    "il6", "il13", "chi3l3", "tgm2") %>%
  sort



# fxns to check if gene provided is (uniquely) in dataset
# ------------------------------------------------------------------------------

apriori_genes_matches <-
  sapply(apriori_genes,
         function(x) {
           DEresults_DvC %>%
             filter(grepl(paste0("^", x), symbol, ignore.case = T)) %>% .$row }) %>%
  .[map_lgl(.x = ., .f = function(x) { length(x) != 0 })]

apriori_genes_matches_symbolref <-
  sapply(apriori_genes,
         function(x) {
           DEresults_DvC %>%
             filter(grepl(paste0("^", x), symbol, ignore.case = T)) %>% .$symbol }) %>%
  .[map_lgl(.x = ., .f = function(x) { length(x) != 0 })]



# plot genes +/- RUVr correction
# (i.e., also checks how RUVr correction works on "validated" genes)
# ------------------------------------------------------------------------------

plot_apriori <-
  function(ensg_value, title, ruv = F, export = F) {
    toplot <-
      plotCounts(deseqobj,
               gene = ensg_value, "Group", returnData = T) %>%
      as_tibble(rownames = "ID")
    
    if (ruv) {
      toplot$count <- ruv_results$normalizedCounts[ensg_value, ]
    } else {
      toplot$count <- assay(vst)[ensg_value, ]
    }
    
    plot_out <-
      ggplot(data = toplot, aes(Group, count, color = Group)) +
      geom_quasirandom(size = 2.5, alpha = 0.6, width = 0.15) +
      scale_y_continuous("Normalized Gene Count", breaks = pretty_breaks(n = 5)) +
      scale_x_discrete(labels = c("C", "D", "D/A")) +
      geom_text_repel(aes(label = gsub("Sample", "", ID)),
                      size = 3) +
      scale_color_manual(values = palette_color_group) +
      theme_few() + theme(legend.position = "none") +
      ggtitle(title)
    
    if (export) {
      ggsave(plot_out, width = 4, height = 3,
             filename = paste0("./Output/20210316_apriori_plots/", title,
                               if_else(ruv, "_RUVr", ""), ".png"))
    }
    
    plot_out
  }




# test few specific examples
# ------------------------------------------------------------------------------

plot_grid(
  plot_apriori("ENSMUSG00000026542.7", "Apcs"),
  plot_apriori("ENSMUSG00000026542.7", "Apcs (RUVr-corrected)", T),
  plot_apriori("ENSMUSG00000021556.13", "Golm1"),
  plot_apriori("ENSMUSG00000021556.13", "Golm1 (RUVr-corrected)", T),
  nrow = 2,
  labels = c("A", "", "B", ""))

gene_annotations_short %>% filter(symbol == "Apcs")
gene_annotations_short %>% filter(symbol == "Hpx")
gene_annotations_short %>% filter(symbol == "Golm1")




# loop through all 
# ------------------------------------------------------------------------------

apriori_plots <-
  lapply(1:length(apriori_genes_matches),
       function(i) {
         plot_apriori(apriori_genes_matches[[i]],
                     str_to_sentence(names(apriori_genes_matches)[i]),
                     export = T) }
       )
  
  lapply(1:length(apriori_genes_matches),
         function(i) {
           plot_apriori(apriori_genes_matches[[i]],
                        str_to_sentence(names(apriori_genes_matches)[i]),
                        ruv = T, export = T) }
  )
  