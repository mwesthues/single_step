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
  .libPaths("~/R/x86_64-pc-linux-gnu-library/3.4/")
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
  Sys.setenv("INTERVAL" = "FHN_1")
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
  dplyr::filter(Interval == pred_interval) %>%
  dplyr::mutate(Rnd_Level2 = as.character(Rnd_Level2))

if (isTRUE(nchar(runs) != 0)) {
  run_seq <- runs %>%
    strsplit(., split = "-") %>%
    purrr::flatten_chr() %>%
    as.integer() %>%
    max() %>%
    seq_len()
  prediction_template <- prediction_template %>%
    dplyr::slice(run_seq)
}

# Get the individual levels of the level2 randomization procedure to load the
# correct subset of data.
rnd2_levels <- prediction_template %>%
  dplyr::pull(Rnd_Level2) %>% 
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
  dplyr::mutate(Pred2 = if_else(is.na(Pred2), true = "", false = Pred2))


## -- PHENOTYPIC DATA ----------------------------------------------------
inb_pheno <- readRDS("./data/derived/maizego/tst_pheno_tibble.RDS")
hyb_pheno <- "./data/derived/uhoh/agro_tibble.RDS" %>%
  readRDS() %>%
  tidyr::separate(
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
    dplyr::slice(i)

  # Collect the current ETA object for the prediction.
  eta_i <- eta_log %>%
    dplyr::mutate(
      Rnd_Level1 = as.numeric(Rnd_Level1),
      Core_Fraction = as.character(Core_Fraction)
    ) %>%
    dplyr::semi_join(
      y = template_i,
      by = c(
        "Extent", "Material", "Scenario", "Pred1", "Pred2", "Core_Fraction",
        "Rnd_Level1"
      )
    ) %>%
    dplyr::pull(Data_Location) %>%
    readRDS()

  if (isTRUE(dplyr::pull(template_i, Material) == "Hybrid")) {
    # Remove duplicated entries. This is important for hybrids with multiple 
    # predictors.
    swap_logic <- purrr::compose(`!`, duplicated)
    duplicated_filter <- eta_i %>%
      purrr::map(.f = `[[`, "X") %>%
      purrr::map(rownames) %>%
      purrr::map(tibble::as_data_frame) %>%
      dplyr::bind_cols() %>%
      t() %>%
      swap_logic()

    eta_order <- eta_i %>%
      .[duplicated_filter] %>%
      purrr::map(.f = `[[`, "X") %>%
      purrr::map(rownames) %>%
      purrr::map(tibble::as_data_frame) %>%
      dplyr::bind_cols() %>%
      dplyr::rename(Dent = "value", Flint = "value1") %>%
      tidyr::unite(col = Genotype, Dent, Flint, sep = "_") %>%
      dplyr::mutate(Genotype = paste0("DF_", Genotype))
  } else if (isTRUE(pull(template_i, Material) == "Inbred")) {
    eta_order <- eta_i %>%
      purrr::map(.f = `[[`, "X") %>%
      purrr::map(rownames) %>%
      purrr::reduce(intersect) %>%
      tibble::as_data_frame() %>%
      dplyr::rename(Genotype = "value")
  }

  # Select test and training set individuals for the current combination.
  rnd_level2_i <- rnd_level2_df %>%
    dplyr::semi_join(
      y = template_i,
      by = c("Extent", "Material", "Scenario", "Rnd_Level2")
    ) %>%
    dplyr::select(TST_Geno, TRN_Geno)

  # Vector with genotype names.
  tst_set_i <- rnd_level2_i %>%
    dplyr::pull(TST_Geno) %>%
    base::unique()

  # Select the correct set of phenotypes.
  if (isTRUE(dplyr::pull(template_i, Material) == "Hybrid")) {
    pheno_i <- hyb_pheno %>%
      dplyr::filter(
        Trait == pull(template_i, Trait),
        Genotype %in% tst_set_i
      ) %>%
      dplyr::select(-Prefix) %>%
      dplyr::rename(Observed = "Value") %>%
      dplyr::right_join(y = eta_order, by = "Genotype")
  } else if (isTRUE(dplyr::pull(template_i, Material) == "Inbred")) {
    pheno_i <- inb_pheno %>%
      dplyr::filter(
        Trait == pull(template_i, Trait),
        Genotype %in% tst_set_i
      ) %>%
      dplyr::rename(Observed = "Value") %>%
      dplyr::right_join(y = eta_order, by = "Genotype")
  }


  # Generate all leave-one-out cross-validation predictions for one parameter.
  eta_seq <- seq_along(tst_set_i)
  one_param_lst <- lapply(eta_seq, FUN = function(j) {

    tst_geno_j <- pheno_i %>%
      dplyr::slice(j) %>%
      dplyr::pull(Genotype)

    # In the case of hybrid genotypes, filter T0 hybrids.
    if (isTRUE(dplyr::pull(template_i, Material) == "Hybrid")) {
      tst_dent <- tst_geno_j %>%
        strsplit(split = "_") %>%
        map_chr(`[[`, 2)

      tst_flint <- tst_geno_j %>%
        strsplit(split = "_") %>%
        map_chr(`[[`, 3)

      pheno_j <- pheno_i %>%
        dplyr::mutate(Observed = if_else(
          Genotype == tst_geno_j | Dent == tst_dent | Flint == tst_flint,
          true = NA_real_,
          false = Observed
        )
      )
    } else if (isTRUE(dplyr::pull(template_i, Material) == "Inbred")) {
      
      pheno_j <- pheno_i %>%
        dplyr::mutate(Observed = if_else(
          Genotype == tst_geno_j,
          true = NA_real_,
          false = Observed
        )
      )
    }



    mod_bglr <- BGLR::BGLR(
      y = dplyr::pull(pheno_j, Observed),
      ETA = eta_i,
      nIter = iter,
      burnIn = iter / 2,
      saveAt = temp,
      verbose = FALSE
    )
    yhat <- mod_bglr$yHat
    sd_yhat <- mod_bglr$SD.yHat
    res <- dplyr::slice(pheno_j, j)
    res$yhat <- yhat[j]
    res$var_yhat <- sd_yhat[j] ^ 2
    res$Observed <- pheno_i %>% dplyr::slice(j) %>% dplyr::pull(Observed)

    # Return the predictions.
    res
  })

  # Get the ID of the current scenario to reference it appropriately for later
  # retrieval.
  row_id <- prediction_template %>%
    dplyr::slice(i) %>%
    dplyr::pull(rowid)

  one_param_df <- one_param_lst %>%
    dplyr::bind_rows() %>%
    dplyr::select(Genotype, Observed, yhat, var_yhat) %>%
    dplyr::bind_cols(slice(template_i, rep(1:n(), length(one_param_lst)))) %>%
    dplyr::mutate(rowid = row_id)


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

