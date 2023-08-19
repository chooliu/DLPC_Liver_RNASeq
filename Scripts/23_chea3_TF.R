# ==============================================================================
# 23_chea3_TF.R
# ==============================================================================



# transcription factor binding analysis via ChEA3
# ultimately not included in final paper since had more actual wet lab data
# ------------------------------------------------------------------------------

library(httr) # httr_1.4.2
library(jsonlite) # jsonlite_1.7.2

url <- "https://maayanlab.cloud/chea3/api/enrich/"
encode <- "json"

run_chea3_api <- function(resultstable) {

  genes <-
    resultstable %>%
    filter(padj < 0.05) %>% .$symbol %>%
    unique

  if (length(genes) > 1000) {
  genes <-
    resultstable %>%
    filter(padj < 0.01) %>% .$symbol %>%
    unique
  }

  payload <- list(query_name = "myQuery", gene_set = genes)
  response <- POST(url = url, body = payload, encode = encode)
  json <- content(response, "text")
  results <- fromJSON(json, flatten = T)
  
  results

  }



# run ChEA3
# ------------------------------------------------------------------------------

chea3_DvC <-
  run_chea3_api(DEresults_DvC)
chea3_DAvC <-
  run_chea3_api(DEresults_DAvC)
chea3_DAvD <-
  run_chea3_api(DEresults_DAvD)


chea3_DvC_RUVr <-
  run_chea3_api(DEresults_RUVr_DvC)
chea3_DAvC_RUVr <-
  run_chea3_api(DEresults_RUVr_DAvC)
chea3_DAvD_RUVr <-
  run_chea3_api(DEresults_RUVr_DAvD)

save(chea3_DvC_RUVr,
     chea3_DAvC_RUVr,
     chea3_DAvD_RUVr,
     chea3_DvC,
     chea3_DAvC,
     chea3_DAvD,
     file = "./Output/20210315_webgestalt/chea3results.Rdata")

load("./Output/20210315_webgestalt/chea3results.Rdata")



# export ChEA3
# ------------------------------------------------------------------------------

prep_chea3 <- function(x) {
  tmp <- 
    x$`Integrated--meanRank` %>% select(3:6) %>%
    filter(!grepl("Enrichr", Library)) %>%
    set_names(c("Transcription Factor", "Score",
                "Reference Datasets", "Overlapping Gene Targets"))
  tmp[1:25, ]
}

