# ==============================================================================
# 22_celldeconv_benchmark.R
# ==============================================================================



# since fairly different results in script 21,
# generate 1,000 mock bulk liver samples --> attempt to benchmark
# ------------------------------------------------------------------------------

# set.seed(1234)
# mockbulkdat <-
#   generateBulk_norep(singlecell_set,
#                      ct.varname = 'celltypes',
#                      sample = "subjects",
#                      ct.sub = allcelltypes,
#                      nbulk = 1000)
#
# save(mockbulkdat, file = "./Output/20210319_celldeconv/BENCHMARK.Rdata")

load("./Output/20210319_celldeconv//BENCHMARK.Rdata")



# normalized mock bulk data
# ------------------------------------------------------------------------------

mockbulkdat_cpm <-
  mockbulkdat

mock_sample_libsize <-
  exprs(mockbulkdat_cpm$pseudo_eset) %>%
  colSums()

exprs(mockbulkdat_cpm$pseudo_eset) <-
  exprs(mockbulkdat_cpm$pseudo_eset) %>%
  apply(., 1, function(x) {
    log(1 + x / mock_sample_libsize)
  }) %>%
  t()

singlecellqc <-
  SCDC_qc(singlecell_set, ct.sub = allcelltypes, ct.varname = "celltypes")
singlecellqc$prop.qc



# run each method on benchmark
# ------------------------------------------------------------------------------

scdc_results_benchmark <-
  SCDC_prop(
    sc.eset = singlecellqc$sc.eset.qc,
    bulk.eset = mockbulkdat$pseudo_eset,
    ct.sub = allcelltypes,
    ct.varname = "celltypes"
  )

music_results_benchmark <-
  music_prop(
    bulk.eset = mockbulkdat$pseudo_eset,
    sc.eset = singlecell_set,
    clusters = "celltypes", samples = "subjects",
    select.ct = allcelltypes
  )


bisque_results_benchmark <-
  ReferenceBasedDecomposition(
    bulk.eset = mockbulkdat$pseudo_eset,
    sc.eset = singlecell_set,
    cell.types = "celltypes", subject.names = "subjects",
    use.overlap = F
  )


# # CIBERSORTx export
# mockbulkdat$pseudo_bulk %>%
#   `colnames<-`(paste0("Mock", 1:1000)) %>%
#   as_tibble(rownames = "Symbol") %>%
#   write_tsv("./Output/20210319_celldeconv/cibersort_mixture_BENCHMARK.txt")

cibersortx_results_benchmark <-
  read_csv("./Data/20210325_celldeconv/CIBERSORTx_Benchmark_Results.csv") %>%
  select(Mixture, all_of(allcelltypes))



# compare to simulated ground truth
# ------------------------------------------------------------------------------

benchmark_fxn <- function(estimate, test) {
  list(
    Pearson =
      sapply(
        allcelltypes,
        function(i) {
          cor(estimate[, i], test[, i])
        }
      ),
    Spearman =
      sapply(
        allcelltypes,
        function(i) {
          cor(estimate[, i], test[, i], method = "spearman")
        }
      ),
    Bias =
      sapply(
        allcelltypes,
        function(i) {
          mean(estimate[, i] - test[, i])
        }
      ),
    MAD =
      sapply(
        allcelltypes,
        function(i) {
          mean(abs(estimate[, i] - test[, i]))
        }
      ),
    MAPD = sapply(
      allcelltypes,
      function(i) {
        mean(abs(estimate[, i] - test[, i]) / test[, i])
      }
    ),
    MSE = sapply(
      allcelltypes,
      function(i) {
        mean((estimate[, i] - test[, i])^2)
      }
    ),
    MAPE = sapply(
      allcelltypes,
      function(i) {
        mean((estimate[, i] - test[, i])^2 / estimate[, i])
      }
    )
  )
}

list_of_estimates <-
  list(
    SCDC = scdc_results_benchmark$prop.est.mvw,
    MUSIC = music_results_benchmark$Est.prop.weighted,
    BISQUE = bisque_results_benchmark$bulk.props %>% t(),
    CIBERSORTx = cibersortx_results_benchmark[, -1] %>% as.matrix()
  )

benchmark_results <-
  list(
    SCDC =
      benchmark_fxn(
        scdc_results_benchmark$prop.est.mvw,
        mockbulkdat$true_p
      ),
    MUSIC = benchmark_fxn(
      music_results_benchmark$Est.prop.weighted,
      mockbulkdat$true_p
    ),
    BISQUE = benchmark_fxn(
      bisque_results_benchmark$bulk.props %>% t(),
      mockbulkdat$true_p
    ),
    CIBERSORTx = benchmark_fxn(
      cibersortx_results_benchmark[, -1] %>% as.matrix(),
      mockbulkdat$true_p
    )
    
    # thought about doing an ensemble learning approach but
    # each individual version very different performance?
    # Ensemble = benchmark_fxn(Reduce("+", list_of_estimates) / 4,
    # mockbulkdat$true_p)
  )




benchmark_results %>%
  map(~ .[["MAPD"]]) %>%
  as_tibble() %>%
  bind_cols(Cell = allcelltypes, .)
benchmark_results %>%
  map(~ .[["MSE"]]) %>%
  as_tibble() %>%
  bind_cols(Cell = allcelltypes, .)
benchmark_results %>%
  map(~ .[["Pearson"]]) %>%
  as_tibble() %>%
  bind_cols(Cell = allcelltypes, .)
benchmark_results %>%
  map(~ .[["Spearman"]]) %>%
  as_tibble() %>%
  bind_cols(Cell = allcelltypes, .)


bisque_results_benchmark$bulk.props %>%
  `-`(t(mockbulkdat$true_p)) %>%
  .[3, ] %>%
  hist()


plot_est_cellcomp <-
  function(df_wide_fmt) {
    df_long_fmt <-
      bind_rows(
        mockbulkdat$true_p %>% as_tibble() %>% bind_cols(Method = "Simulated Data"),
        bisque_results_benchmark$bulk.props %>% t() %>% as_tibble() %>% bind_cols(Method = "Bisque"),
        cibersortx_results_benchmark[, -1] %>% bind_cols(Method = "CIBERSORTx"),
        music_results_benchmark$Est.prop.weighted %>% as_tibble() %>% bind_cols(Method = "MuSIC"),
        scdc_results_benchmark$prop.est.mvw %>% as_tibble() %>% bind_cols(Method = "SCDC"),
      ) %>%
      pivot_longer(cols = all_of(allcelltypes), names_to = "Celltype")

    df_long_fmt <-
      df_long_fmt %>%
      arrange(Method != "Simulated Data") %>%
      mutate(Method = fct_inorder(Method))

    ggplot(
      df_long_fmt,
      aes(Method, value * 100,
        color = Method,
        alpha = (Method == "Simulated Data")
      )
    ) +
      geom_quasirandom(width = 0.2, shape = 21) +
      facet_grid(Celltype ~ .,
        scales = "free", switch = "y",
        labeller = label_wrap_gen()
      ) +
      scale_y_continuous("Percent", breaks = pretty_breaks(n = 5)) +
      scale_color_manual(values = c("black", qualitative_hcl(palette = "Warm", n = 4))) +
      scale_alpha_manual(values = c(0.2, 1)) +
      theme_few() +
      theme(legend.position = "none")
  }

plot_est_cellcomp()





# some methods better w.r.t. correlation (relative cell proportions)
# but others better w.r.t. absolute values?
# ------------------------------------------------------------------------------

benchmark_pearson_plot <-
  benchmark_results %>%
  map(~ .[["Pearson"]]) %>%
  as_tibble() %>%
  bind_cols(Cell = allcelltypes, .) %>%
  pivot_longer(cols = 2:ncol(.)) %>%
  ggplot(data = ., aes(name, str_wrap(Cell, width = 20), fill = value)) +
  geom_tile(alpha = 0.5) +
  geom_text(aes(label = formatC(value, format = "f", digits = 2, ))) +
  scale_fill_continuous_diverging(name = "", "Tofino", limits = c(0, 1)) +
  theme_minimal() +
  xlab("Method") +
  ylab("Cell Type") +
  ggtitle("Pearson Correlation", "(larger values are better)")

benchmark_mad_plot <-
  benchmark_results %>%
  map(~ .[["MAPD"]]) %>%
  as_tibble() %>%
  bind_cols(Cell = allcelltypes, .) %>%
  pivot_longer(cols = 2:ncol(.)) %>%
  ggplot(data = ., aes(name, str_wrap(Cell, width = 20), fill = value)) +
  geom_tile(alpha = 0.5) +
  geom_text(aes(label = formatC(value, format = "f", digits = 2, ))) +
  scale_fill_continuous_diverging(name = "", "Tofino", rev = T) +
  theme_minimal() +
  xlab("Method") +
  ylab("Cell Type") +
  ggtitle("Mean Absolute Percent Diff.", "(smaller values are better)")

# 1000 x 350
plot_grid(benchmark_mad_plot,
  benchmark_pearson_plot,
  labels = c("A", "B")
)
