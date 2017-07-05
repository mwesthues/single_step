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
if (isTRUE(interactive())) {
  .libPaths(c(.libPaths(), "~/R/x86_64-pc-linux-gnu-library/3.4/"))
}
if (!require("pacman")) install.packages("pacman")
if (!require("devtools")) install.packages("devtools")
pacman::p_load("BGLR","data.table", "parallel", "magrittr", "dplyr", "dtplyr",
               "tidyr", "purrr", "methods", "tibble", "lubridate")
#devtools::install_github("mwesthues/sspredr")
pacman::p_load_gh("mwesthues/sspredr")

# Load function to save session info on the software used for the analysis in
# this script.
source("./software/session_info.R")
# Prediction helper functions (outsourced from this script for improved
# readability).
source("./software/prediction_helper_functions.R")



## --------------------------------------------------------------------------
## COLLECT AND CHECK ARGUMENTS 
if (isTRUE(interactive())) {
  # Number of cores
  Sys.setenv("MOAB_PROCCOUNT" = "3")
  # Number of iterations in BGLR()
  Sys.setenv("ITER" = "3000") 
  # Prediction model in BGLR()
  Sys.setenv("MODEL" = "BRR")
  # Algorithm to generate variance-covariance matrices.
  Sys.setenv("VCOV" = "RadenII")
  Sys.setenv("PI" = "0.5")
  Sys.setenv("PRIOR_PI_COUNT" = "10")
  # Which interval of data shall be analyzed?
  Sys.setenv("INTERVAL" = "2")
  # Number of test runs.
  Sys.setenv("RUNS" = "1-3")
  # Output directory for temporary BGLR files
  Sys.setenv("TMP" = "./tmp")
}

if (!interactive()) {
  job_id <- Sys.getenv("MOAB_JOBID")
} else job_id <- "interactive_00"

# Save session info
write_session_info(directory = "./data/derived/session_info/",
                   job_id = job_id)

use_cores <- as.integer(Sys.getenv("MOAB_PROCCOUNT"))
material <- as.character(Sys.getenv("DATA_TYPE"))
traits <- as.character(Sys.getenv("TRAIT"))
iter <- as.integer(Sys.getenv("ITER"))
pred_model <- as.character(Sys.getenv("MODEL"))
g_method <- as.character(Sys.getenv("VCOV"))
Pi <- as.numeric(Sys.getenv("PI"))
PriorPiCount <- as.integer(Sys.getenv("PRIOR_PI_COUNT"))
temp <- paste0(as.character(Sys.getenv("TMP")), "/")
pred_interval <- as.character(Sys.getenv("INTERVAL"))
runs <- as.character(Sys.getenv("RUNS"))



# -- PREDICTOR INPUT SPECIFICATION -------------------------------------
temp_loc <- "./data/derived/prediction_runs/prediction_template.RDS"
original_prediction_template <- temp_loc %>%
  readRDS()

prediction_template <- original_prediction_template %>%
  filter(Interval == pred_interval) %>%
  mutate(Rnd_Level2 = as.character(Rnd_Level2))

if (isTRUE(nchar(runs) != 0)) {
  run_seq <- runs %>%
    strsplit(., split = "-") %>%
    flatten_chr() %>%
    as.integer() %>%
    max() %>%
    seq_len()
  prediction_template <- prediction_template %>%
    slice(run_seq)
}

# Get the individual levels of the level2 randomization procedure to load the
# correct subset of data.
rnd2_levels <- prediction_template %>%
  pull(Rnd_Level2) %>% 
  unique() %>%
  as.integer()

if (isTRUE(all(rnd2_levels %in% seq_len(25)))) {
  rnd_level2_df <- readRDS(
    "./data/derived/predictor_subsets/rnd_level2_1-25.RDS"
  )
} else if (isTRUE(all(rnd2_levels %in% seq(from = 26, to = 50, by = 1)))) {
  rnd_level2_df <- readRDS(
    "./data/derived/predictor_subsets/rnd_level2_26-50.RDS"
  )
}


# Get the log file for the ETA objects.
eta_log <- "./data/derived/log_prepare_subsamples.txt" %>%
  readr::read_tsv() %>%
  mutate(Pred2 = if_else(is.na(Pred2), true = "", false = Pred2))


## -- PHENOTYPIC DATA ----------------------------------------------------
inb_pheno <- readRDS("./data/derived/maizego/tst_pheno_tibble.RDS")
hyb_pheno <- "./data/derived/uhoh/agro_tibble.RDS" %>%
  readRDS() %>%
  separate(
    col = Genotype,
    into = c("Prefix", "Dent", "Flint"),
    sep = "_",
    remove = FALSE
  )


## -- PREDICTION -----------------------------------------------------------
# Define which runs shall be used, i.e., which genotypes shall be included as 
# test sets.
pred_seq <- prediction_template %>%
  nrow() %>%
  seq_len()

# Keep the time of the following computations.
start_time <- Sys.time()

# Analysis
mclapply(pred_seq, FUN = function(i) {

  # Collect the current set of parameters.
  template_i <- prediction_template %>%
    slice(i)

  # Collect the current ETA object for the prediction.
  eta_i <- eta_log %>%
    mutate(
      Rnd_Level1 = as.numeric(Rnd_Level1),
      Core_Fraction = as.character(Core_Fraction)
    ) %>%
    semi_join(
      y = template_i,
      by = c(
        "Extent", "Material", "Scenario", "Pred1", "Pred2", "Core_Fraction",
        "Rnd_Level1"
      )
    ) %>%
    pull(Data_Location) %>%
    readRDS()

  if (isTRUE(pull(template_i, Material) == "Hybrid")) {
    eta_order <- eta_i %>%
      map(.f = `[[`, "X") %>%
      map(rownames) %>%
      map(as_data_frame) %>%
      bind_cols() %>%
      rename(Dent = "value", Flint = "value1") %>%
      unite(col = Genotype, Dent, Flint, sep = "_") %>%
      mutate(Genotype = paste0("DF_", Genotype))
  } else if (isTRUE(pull(template_i, Material) == "Inbred")) {
    eta_order <- eta_i %>%
      map(.f = `[[`, "X") %>%
      map(rownames) %>%
      reduce(intersect) %>%
      as_data_frame() %>%
      rename(Genotype = "value")
  }

  # Select test and training set individuals for the current combination.
  rnd_level2_i <- rnd_level2_df %>%
    semi_join(
      y = template_i,
      by = c("Extent", "Material", "Scenario", "Rnd_Level2")
    ) %>%
    select(TST_Geno, TRN_Geno)

  # Vector with genotype names.
  tst_set_i <- rnd_level2_i %>%
    pull(TST_Geno) %>%
    unique()

  # Select the correct set of phenotypes.
  if (isTRUE(pull(template_i, Material) == "Hybrid")) {
    pheno_i <- hyb_pheno %>%
      filter(
        Trait == pull(template_i, Trait),
        Genotype %in% tst_set_i
      ) %>%
      select(-Prefix) %>%
      rename(Observed = "Value") %>%
      right_join(y = eta_order, by = "Genotype")
  } else if (isTRUE(pull(template_i, Material) == "Inbred")) {
    pheno_i <- inb_pheno %>%
      filter(
        Trait == pull(template_i, Trait),
        Genotype %in% tst_set_i
      ) %>%
      rename(Observed = "Value") %>%
      right_join(y = eta_order, by = "Genotype")
  }


  # Generate all leave-one-out cross-validation predictions for one parameter.
  eta_seq <- seq_along(tst_set_i)
  one_param_lst <- lapply(eta_seq, FUN = function(j) {

    tst_geno_j <- pheno_i %>%
      slice(j) %>%
      pull(Genotype)

    pheno_j <- pheno_i %>%
      mutate(Observed = if_else(
        Genotype == tst_geno_j,
        true = NA_real_,
        false = Observed
      )
    )

    mod_bglr <- BGLR::BGLR(
      y = pull(pheno_j, Observed),
      ETA = eta_i,
      nIter = iter,
      burnIn = iter / 2,
      saveAt = temp,
      verbose = FALSE
    )
    yhat <- mod_bglr$yHat
    sd_yhat <- mod_bglr$SD.yHat
    res <- slice(pheno_j, j)
    res$yhat <- yhat[j]
    res$var_yhat <- sd_yhat[j] ^ 2
    res$Observed <- pheno_i %>% slice(j) %>% pull(Observed)

    # Return the predictions.
    res
  })

  # Get the ID of the current scenario to reference it appropriately for later
  # retrieval.
  row_id <- prediction_template %>%
    slice(i) %>%
    pull(rowid)

  one_param_df <- one_param_lst %>%
    bind_rows() %>%
    select(Genotype, Observed, yhat, var_yhat) %>%
    bind_cols(slice(template_i, rep(1:n(), length(one_param_lst)))) %>%
    mutate(rowid = row_id)


  # Save output and link it with a log file.
  saveRDS(
    one_param_df,
    paste0("./data/derived/predictions/template_", row_id, ".RDS")
  )

}, mc.cores = use_cores)


# Write a corresponding log file for the predictions.
elapsed_time <- get_elapsed_time(start_time)
res_log <- data_frame(
  Start_Time = start_time,
  Finish_Time = Sys.time(),
  Elapsed_Time = elapsed_time,
  Interval = pred_interval,
  Runs = runs,
  JOB_ID = job_id
)

res_log_location <- "./data/derived/prediction_log.txt"
readr::write_tsv(
  res_log,
  path = res_log_location,
  append = if_else(file.exists(res_log_location), true = TRUE, false = FALSE)
)

