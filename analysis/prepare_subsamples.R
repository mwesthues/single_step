# Goal: Sample core sets of genotypes for the Yan-lab data and scenario two of
# our predictions in increments of ten percentage points.
if (!require("pacman")) install.packages("pacman")
if (!require("devtools")) install.packages("devtools")
#devtools::install_github("mwesthues/sspredr", update = FALSE)
pacman::p_load("tidyverse", "data.table", "dtplyr")
pacman::p_load_gh("mwesthues/sspredr")

# For the second scenario, keep only names of genotypes, which are covered by 
# all data types (i.e. phenotypic, genotypic and transcriptomic).
common_genotypes <- readRDS(
  "./data/derived/maizego/unique_snp-mrna-pheno_genotypes.RDS"
)
common_genotypes <- common_genotypes %>%
  reduce(intersect)

## --- Scenario 2 SNP preparation.
# Explore the kinship matrix.
genos <- "./data/processed/maizego/imputed_snp_mat.RDS" %>%
  readRDS() %>%
  rownames()

runs <- 10
fractions = seq(from = 0.1, to = 1, by = 0.1)



#### Functions
compute_totals_from_fractions <- function(nm, frac) {
  # nm: vector of names
  # frac: fractions
  totals <- nm %>%
    length() %>%
    map(function(i) i * frac) %>%
    flatten_dbl() %>%
    ceiling()
}


sample_fraction <- function(nm, frac) {
  # nm: names of genotypes
  # frac: vector of genotypes fractions
  totals <- compute_totals_from_fractions(
    nm = common_genotypes,
    frac = frac
  )
  totals %>%
    map(.f = ~sample(
      x = nm,
      size = .,
      replace = FALSE
    )) %>%
  set_names(., nm = as.character(frac)) %>%
  stack()
}


sample_loo <- function(geno, frac, iter) {
  train_geno <- map(geno, .f = ~ setdiff(geno, .))
  names(train_geno) <- geno
  train_lst <- rerun(.n = iter, {
    sub_train_lst <- map(train_geno, .f = ~sample_fraction(., frac = frac))
    sub_train_df <- bind_rows(sub_train_lst, .id = "TST_Geno")
    sub_train <- dplyr::rename(sub_train_df, TRN_Geno = values)
  })
  iter_seq <- seq_len(iter)
  names(train_lst) <-  as.character(iter_seq)
  train_df <- bind_rows(train_lst, .id = "Iter")
  as_data_frame(train_df)
}
###





set.seed(9434)
loocv_samples <- sample_loo(
  geno = genos,
  frac = 0.7,
  iter = runs
)

set.seed(9434)
existing_mrnas <- rerun(.n = runs, sample_fraction(
  nm = common_genotypes,
  frac = fractions
  )) %>%
  bind_rows(.id = "Rep") %>%
  as_data_frame()


