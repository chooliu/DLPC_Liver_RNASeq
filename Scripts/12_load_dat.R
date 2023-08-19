# ==============================================================================
# 12_load_dat.R
# ==============================================================================



# load annot, count matrix, metadata from S.G.
# ------------------------------------------------------------------------------

gene_annotations_short <- 
  read_tsv("./Data/gencode_annotations.txt")

count_mat <-
  read_tsv("./Output/20210307_counts.txt")

metadata <-
  read_excel("./Data/RNA seq sample ID.xlsx", skip = 4) %>%
  pivot_longer(1:3, names_to = "Group", values_to = "ID") %>%
  filter(!is.na(ID)) %>%
  mutate(Group = as.factor(Group) %>% fct_recode(`Disease/Agonist` = "Disease/agonist")) %>%
  select(ID, Group) %>%
  arrange(ID)


# sample labeled 724 in counts matrix is actually 725,
# and fix metadata typo
# ------------------------------------------------------------------------------
metadata$ID %>% .[duplicated(.)]
metadata$ID[metadata$ID == "726"][1] <- "725" # typo, two 726s one should be 725
metadata <- metadata %>% arrange(ID)

names(count_mat)[names(count_mat) == "Sample724"] <- "Sample725"

# check ordering
count_mat %>% names %>% gsub("Sample", "", .) %>% .[-1] %>% order
count_mat %>% names %>% gsub("Sample", "", .) %>% .[-1] %>%
  setdiff(metadata$ID %>% as.character())

count_mat %>% names %>% gsub("Sample", "", .) %>% .[-1] %>% identical(., metadata$ID)



# extra metadata about batch info (litter)
# ------------------------------------------------------------------------------
batch <-
  read_excel("./Data/batch no.xlsx") %>%
  transmute(ID = `Mouse No` %>% as.character,
            Date = `Date of Experiment` %>% as_date(),
            Month = month(Date),
            Year = year(Date),
            # Quarter = quarter(Date, with_year = T),
            Batch = paste0(Year, "_", str_pad(Month, width = 2, side = "left", pad = "0")),
            Notes = `...3`)

metadata <-
  metadata %>%
  left_join(., batch %>% select(ID, Batch))


