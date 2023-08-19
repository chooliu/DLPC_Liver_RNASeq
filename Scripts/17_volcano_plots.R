# ==============================================================================
# 17_volcano_plots.R
# ==============================================================================



# helper fxn to make volcano plots
# ------------------------------------------------------------------------------

make_volcano_plot <-
  function(deresults_table,
         plottitle,
         righthandside) {
  
  tmp_row <-
    c(deresults_table %>%
        arrange(pvalue) %>%
        .[1:10, ] %>%
        .$row,
      deresults_table %>%
        arrange(-abs(log2FoldChange)) %>%
        .[1:10, ] %>%
        .$row)
  tmp <-
    deresults_table %>%
    filter(row %in% tmp_row) %>%
    filter(!is.na(pvalue)) %>%
    mutate(neglogten = -log10(pvalue),
           Signif = if_else(padj < 0.05, "yes", "no"))
    # optionally select genes to label
    # .[c(1:10, (nrow(.) - 10):nrow(.)), ]
    # & abs(log2FoldChange) > 4)
  
  # 600 x 500
  
  ggplot(deresults_table %>%
           filter(!is.na(pvalue)) %>%
           mutate(neglogten = -log10(pvalue),
                  Signif = if_else(padj < 0.05, "yes", "no")) %>%
           filter(!is.na(Signif)),
         aes(log2FoldChange, neglogten, color = Signif, shape = Signif)) +
    geom_vline(xintercept = 0) + geom_hline(yintercept = 0) +
    geom_point(alpha = 0.5) +
    geom_text_repel(data = tmp, aes(label = symbol),
                    box.padding = 0.05, segment.alpha = 0.4,
                    show.legend = F, color = "black", size = 3) +
    # annotate(geom = "text", x = -8, y = 10, label = "higher in D", hjust = 0, color = "black") +
    # annotate(geom = "text", x = 5, y = 10, label = "higher in D/A", hjust = 1, color = "black") +
    theme_few() + theme(legend.position = "none") +
    scale_color_manual(values = palette_color_signif) +
    scale_shape_manual(values = palette_shape_signif) +
    scale_y_continuous(expression(-log[10]~"(p-value)"), breaks = pretty_breaks(n = 5)) +
    xlab(expression(log[2]~"(Fold Change)")) +
    labs(title = plottitle,
         subtitle = paste0("right = higher in ", righthandside))
}



# make the volcano plots
# ------------------------------------------------------------------------------

volcano_plots <-
  plot_grid(
  make_volcano_plot(DEresults_DvC, "D versus C (no adjustment)", "D"),
  make_volcano_plot(DEresults_RUVr_DvC, "D versus C (RUVr-adjusted)", "D"),
  make_volcano_plot(DEresults_DAvC, "D/A versus C (no adjustment)", "D/A"),
  make_volcano_plot(DEresults_RUVr_DAvC, "D/A versus C (RUVr-adjusted)", "D/A"),
  make_volcano_plot(DEresults_DAvD, "D/A versus D (no adjustment)", "D/A"),
  make_volcano_plot(DEresults_RUVr_DAvD, "D/A versus D (RUVr-adjusted)", "D/A"),
  nrow = 3, ncol = 2,
  labels = LETTERS[1:6]
  )

ggsave(volcano_plots,
    file = "Output/20210318_volcanoplots.png",
    width = 8, height = 10)

system('open "Output/20210318_volcanoplots.png"')




# specific volcano plot requests (STAT6?)
# ------------------------------------------------------------------------------

volcano_requests <-
  read_excel("./Data/RNA seq -figure ENSG.xlsx", skip = 1) %>%
  .[1:16, ]

volcano_list_apriori <-
  volcano_requests$`...1`

make_volcano_plot_focus <-
  function(deresults_table,
           plottitle,
           righthandside) {
    
    ggplot(deresults_table %>%
             filter(row %in% volcano_list_apriori) %>%
             mutate(neglogten = -log10(pvalue),
                    Signif = if_else(padj < 0.05, "yes", "no"),
                    symbol = if_else(pvalue < 0.05, symbol, "")) %>% # can make false = "" for alt
             filter(!is.na(Signif)),
           aes(log2FoldChange, neglogten, color = Signif, shape = Signif)) +
      geom_vline(xintercept = 0, alpha = 0.5, lty = 1) +
      geom_hline(yintercept = 0, alpha = 0.5, lty = 1) +
      geom_hline(yintercept = -log10(0.05), alpha = 0.5, lty = 3) +
      geom_point(alpha = 0.5) +
      geom_text_repel(aes(label = symbol, color = Signif),
                      box.padding = 0.05, segment.alpha = 0.4,
                      show.legend = F, size = 3, max.overlaps = 15) +
      theme_few() + theme(legend.position = "none") +
      scale_color_manual(values = palette_color_signif_darker) +
      scale_shape_manual(values = palette_shape_signif) +
      scale_y_continuous(expression(-log[10]~"(p-value)"),
                         breaks = breaks_pretty(n = 5)) +
      scale_x_continuous(expression(log[2]~"(Fold Change)"),
                         breaks = breaks_pretty(n = 6),
                         expand = c(0.2, 0.2)) +
      labs(title = plottitle,
           subtitle = paste0("logFC > 0: higher expression in ", righthandside))
  }


volcano_plots <-
  plot_grid(
    make_volcano_plot_focus(DEresults_RUVr_DvC, "DSS-PN vs. Chow", "DSS-PN"),
    make_volcano_plot_focus(DEresults_RUVr_DAvC, "DSS-PN/DLPC vs. Chow", "DSS-PN/DLPC"),
    make_volcano_plot_focus(DEresults_RUVr_DAvD, "DSS-PN/DLPC vs. DSS-PN", "DSS-PN/DLPC"),
    nrow = 3, ncol = 1,
    labels = LETTERS[1:3]
  )

ggsave(volcano_plots,
       file = "Figures/Fig4.png",
       width = 4.5, height = 10, dpi = 1200)

system('open "Figures/Fig4.png"')



