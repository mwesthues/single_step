# Goal: Add three new core sets (based on a previous admixture analysis) to the
# existing list of core samples.
# our predictions in increments of ten percentage points.
if (!require("pacman")) install.packages("pacman")
if (!require("devtools")) install.packages("devtools")
devtools::install_github("mwesthues/sspredr", update = TRUE)
options(java.parameters = "-Xmx25G")
pacman::p_load("tidyverse", "data.table", "dtplyr", "corehunter", "LEA", 
               "ggthemes")
pacman::p_load_gh("mwesthues/sspredr")

# This data frame contains the ancestry coefficients of the genotypes that are
# covered by both, the complete (SNP) and the incomplete (mRNA) predictors, 
# based on an admixture analysis assuming three ancestral populations.
ancestry_df <- "./data/derived/maizego/ancestry_k3.RDS" %>%
  readRDS()

# Genotypes with a share >= 50% of the second ancestral population (A2) comprise
# a far larger set (n = 109) than genotypes with either a share >= 50% of the
# first or the third ancestral population (n = 17, n = 19).
# To reduce the bias of the sample size on predictions using these three 
# ancestral populations, we will build a core set from the A2-genotypes of 
# similar size as the two other sets (n = 19).
a2 <- ancestry_df %>%
  filter(Ancestor == "A2") %>% 
  select(G) %>% 
  flatten_chr()

a2_core <- "./data/processed/maizego/imputed_snp_mat.RDS" %>%
  readRDS() %>%
  .[rownames(.) %in% a2, ] %>%
  sspredr::ensure_snp_quality(
    ., callfreq_check = FALSE, maf_check = TRUE, maf_threshold = 0.05,
    any_missing = FALSE, remove_duplicated = TRUE
  ) %>%
  genotypes(format = "biparental") %>%
  sampleCore(
    ., size = 19, obj = objective(type = "EN", measure = "MR", weight = 1)
  )

# Genotypes belonging largely to the first or third ancestral population will
# not be modified and simply be assigned to the list with the names of core set
# genotypes.
a1_a3 <- ancestry_df %>% 
  split(.$Ancestor) %>% 
  map("G") %>% 
  discard(names(.) == "A2") %>%
  set_names(paste0("selected_fraction_", names(.)))
a2_core <- a2_core %>%
  keep(names(.) == "sel") %>% 
  set_names(paste0("selected_fraction_A2"))
a_lst <- a1_a3 %>%
  prepend(a2_core, before = 2)

core_lst <- "./data/derived/maizego/scenario2_snp_core_list.RDS" %>% 
  readRDS() %>% 
  purrr::transpose() %>% 
  keep(names(.) == "sel") %>% 
  .[[1]] %>% 
  append(a_lst)
saveRDS(core_lst, "./data/derived/maizego/augmented_core_list.RDS")

