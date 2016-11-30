# Single-step prediction of hybrid performance based on Fernando et al.: A
# class of Bayesian methods to combine large numbers of genotyped and
# non-genotyped animals for whole-genome analyses. Genetics Selection Evolution
# 2014 46:50
#
# The approach will be adjusted in a sense that we do have complete genomic
# instead of complete pedigree data (i.e. A -> G) and that we have incomplete
# transcriptomic instead of incomplete genomic data (i.e. G -> T).
#
#
#
if (!require("pacman")) install.packages("pacman")
if (!require("devtools")) install.packages("devtools")
devtools::install_github("mwesthues/sspredr", update = TRUE)
pacman::p_load("BGLR","data.table", "parallel", "magrittr", "dplyr", "tidyr",
               "purrr", "testthat", "methods", "tibble", "ggplot2", "stringr",
               "stringi", "lubridate", "readr")
pacman::p_load_gh("mwesthues/sspredr")


## --------------------------------------------------------------------------
## COLLECT AND CHECK ARGUMENTS 
if (isTRUE(interactive())) {
  # Number of cores
  Sys.setenv("MOAB_PROCCOUNT" = "3")
  Sys.setenv("TRAIT" = "GTM")
  # Number of iterations in BGLR()
  Sys.setenv("ITER" = "30000")
  # Prediction model in BGLR()
  Sys.setenv("MODEL" = "BRR")
  # Algorithm to generate variance-covariance matrices.
  Sys.setenv("VCOV" = "RadenII")
  Sys.setenv("PI" = "0.5")
  Sys.setenv("PRIOR_PI_COUNT" = "10")
  # Main predictor. If 'Pred2' and 'Pred3' are empty, no imputation will take
  # place.
  Sys.setenv("PRED1" = "ped100")
  # If 'Pred3' is empty, 'Pred2' will be imputed via information from 'Pred1'.
  Sys.setenv("PRED2" = "snp77")
  # If not empty, this predictor will be imputed.
  Sys.setenv("PRED3" = "")
  # Number of genotypes to predict (only for testing!)
  Sys.setenv("RUNS" = "")
}

if (!interactive()) {
  job_id <- Sys.getenv("MOAB_JOBID")
} else job_id <- "interactive_00"

use_cores <- as.integer(Sys.getenv("MOAB_PROCCOUNT"))
init_traits <- as.character(Sys.getenv("TRAIT"))
init_iter <- as.integer(Sys.getenv("ITER"))
hypred_model <- as.character(Sys.getenv("MODEL"))
g_method <- as.character(Sys.getenv("VCOV"))
Pi <- as.numeric(Sys.getenv("PI"))
PriorPiCount <- as.integer(Sys.getenv("PRIOR_PI_COUNT"))
pred1 <- as.character(Sys.getenv("PRED1"))
pred2 <- as.character(Sys.getenv("PRED2"))
pred3 <- as.character(Sys.getenv("PRED3"))
runs <- as.character(Sys.getenv("RUNS"))


# Input tests
poss_traits <- c("GTM", "GTS", "FETT", "RPR", "STA", "XZ", "ADF")
if (isTRUE(nchar(init_traits) == 0)) {
  init_traits <- poss_traits
}
test_that("selected trait exists", {
  expect_true(all(init_traits %in% poss_traits))
})
#test_that("selected model is part of BGLR", {
#  expect_true(hypred_model %in% c("BRR", "BRR_Kernel", "BayesB", "BayesC"))
#})
test_that("kernel method is defined", {
  expect_true(g_method %in% c("RadenI", "RadenII", "Zhang", "none"))
})




## -- DATA SELECTION -----------------------------------------------------
# Load the predictor data.
pred_nms <- list(pred1 = pred1, pred2 = pred2, pred3 = pred3)
pred_lst <- pred_nms %>%
  keep(~nchar(.) != 0) %>%
  map(~paste0("./data/derived/predictor_subsets/", ., ".RDS")) %>%
  map(readRDS)

# Names of inbred lines for which at least one predictor type has records.
pred_geno_nms <- pred_lst %>%
  map(rownames) %>%
  reduce(union)

# Select hybrids for whose parent lines at least one predictor has records.
genos <- readRDS("./data/processed/common_genotypes.RDS")
hybrid <- genos %>%
  filter(Pool == "Hybrid") %>%
  mutate(G = stringr::str_replace(G, pattern = "DF_", replacement = "")) %>%
  separate(G, into = c("Dent", "Flint"), sep = "_") %>%
  filter(Dent %in% pred_geno_nms, Flint %in% pred_geno_nms) %>%
  unite(col = Hybrid, Dent, Flint, sep = "_", remove = FALSE) %>%
  mutate(Hybrid = paste0("DF_", Hybrid))

# Store the names of all Dent and Flint hybrid parents in a list in order to
# correctly augment the predictor matrices with their help.
hybrid_parent_nms <- hybrid %>%
  gather(key = "Group", value = "G", Dent, Flint) %>%
  split(.$Group) %>%
  map("G")


## -- PREDICTOR AND AGRONOMIC DATA PREPARATION ----------------------------
# Agronomic data
pheno <- readRDS("./data/processed/Pheno_stage2.RDS")
pheno <- pheno %>%
  as_data_frame() %>%
  filter(G %in% hybrid$Hybrid) %>%
  dplyr::select(G, EST, Trait) %>%
  spread(key = Trait, value = EST) %>%
  dplyr::select(-ADL) %>%
  as.data.frame %>%
  remove_rownames() %>%
  column_to_rownames(var = "G") %>%
  as.matrix

# Split the predictor matrices into Dent and Flint.
pred_lst <- pred_lst %>%
  map(as.data.frame) %>%
  map(rownames_to_column, var = "G") %>%
  map(as_tibble) %>%
  map(~left_join(x = .,
                 y = genos %>% 
                       select(G, Pool) %>%
                       distinct(G, .keep_all = TRUE),
                 by = "G")) %>%
  map(~split(., .$Pool)) %>%
  transpose() %>%
  at_depth(.depth = 2, ~select(., -Pool)) %>%
  at_depth(.depth = 2, .f = as.data.frame) %>%
  at_depth(.depth = 2, .f = column_to_rownames, var = "G") %>%
  at_depth(.depth = 2, as.matrix)
pred_lst[["Dent"]][["geno"]] <- hybrid_parent_nms$Dent
pred_lst[["Flint"]][["geno"]] <- hybrid_parent_nms$Flint



## -- GENERATE BGLR() INPUT MATRICES --------------------------------------
# Number of predictors used.
n_pred <- pred_nms %>%
  keep(~nchar(.) != 0) %>%
  as_vector() %>%
  length()

# If pedigree data constituted the first specified predictor, don't generate
# a kernel from them but use a simple Cholesky decomposition instead.
pedigree_first <- pred_nms %>%
  .[[1]] %>%
  str_detect(., pattern = "ped")

if (isTRUE(n_pred == 1)) {
  if (isTRUE(pedigree_first)) {
    eta_lst <- pred_lst %>%
      map(., ~complete_eta(x = .$pred1, geno = .$geno, 
                           as_kernel = FALSE,
                           is_pedigree = TRUE,
                           bglr_model = "BRR"))
  } else {
    eta_lst <- pred_lst %>%
      map(., ~complete_eta(x = .$pred1, geno = .$geno, 
                           as_kernel = TRUE,
                           bglr_model = "BRR"))
  }

} else if (isTRUE(n_pred == 2)) {
  if (isTRUE(pedigree_first)) {
    pred_lst %>%
      map(., ~impute_eta(x = .$pred1, y = .$pred2, geno = .$geno,
                         as_kernel = FALSE, bglr_model = "BRR")) %>%
      str()
  }
    eta_lst <- pred_lst %>%
      map(., ~impute_eta(x = .$pred1, y = .$pred2, geno = .$geno,
                         as_kernel = TRUE, bglr_model = "BRR"))

} else if (isTRUE(n_pred == 3)) {

  pred_lst %>%
    map(., ~impute2(ped = .$pred1, snp = .$pred2, mrna = .$pred3, 
                    geno = .$geno,
                    bglr_model = "BRR")) %>%
    at_depth(.depth = 2, .f = "X") %>%
    map(2) %>%
    map(~.[match(colnames(.), rownames(.)), ]) %>%
    map(~.[lower.tri(., diag = FALSE)]) %>%
    map(summary)
}


## --------------------------------------------------------------------------
#
## PREDICTION
# Determine how much time (hh:mm:ss format) has elapsed since script initiation.
get_elapsed_time <- function(start_time, tz = "CEST") {
  start_time <- as.POSIXct(start_time)
  dt <- difftime(Sys.time(), start_time, units = "secs")
  # Since you only want the H:M:S, we can ignore the date...
  # but you have to be careful about time-zone issues
  format(.POSIXct(dt, tz = tz), "%H:%M:%S")
}

# Define which runs shall be used, i.e., which genotypes shall be included as 
# test sets.
if (nchar(runs) != 0) {
  run_length <- runs %>%
    str_split(., pattern = "-") %>%
    as_vector() %>%
    as.integer() %>%
    invoke(.f = seq, .x = ., by = 1) 
} else {
  run_length <- seq_len(nrow(pheno))
}
param_df <- expand.grid(Trait = init_traits,
                        Iter = init_iter,
                        Run = run_length)
param_df$Trait <- as.character(param_df$Trait)

