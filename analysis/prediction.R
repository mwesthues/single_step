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
  Sys.setenv("ITER" = "2000")
  # Algorithm to generate variance-covariance matrices.
  Sys.setenv("VCOV" = "RadenII")
  Sys.setenv("PI" = "0.5")
  Sys.setenv("PRIOR_PI_COUNT" = "10")
  # Which interval of data shall be analyzed?
  Sys.setenv("INTERVAL" = "CIC_3")
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
iter <- as.integer(Sys.getenv("ITER"))
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

# For scenario 'C' we have not created entries in the prediction template
# because it is exactly the same as that for scenario 'B'.
# In order to filter entries in the template anyway, we need to temporarily
# replace the letter 'C' with the letter 'B'.
if (grepl("CIC|FIC", x = pred_interval)) {
  original_pred_interval <- pred_interval
  substr(pred_interval, start = 3, stop = 3) <- "B"
}

prediction_template <- original_prediction_template %>%
  tibble::as_data_frame() %>%
  dplyr::filter(Interval == pred_interval) %>%
  dplyr::mutate(ETA_UUID = paste0("eta_", ETA_UUID))


if (grepl("CIC|FIC", x = original_pred_interval)) {

  # Function to replace the letter 'B' at position 3 with the letter 'C'.
  sub_b_with_c <- function(x) {
    substr(x, start = 3, stop = 3) <- "C"
    x
  }

  prediction_template %>%
    dplyr::mutate(
      Combi = sub_b_with_c(Combi),
      Interval = sub_b_with_c(Interval)
    )
  pred_interval <- sub_b_with_c(pred_interval)
}




# -- LOAD ADDITIONAL INFO -------------------------------------------------------
material <- pred_interval %>%
  substr(., start = 2, stop = 2)

core_fraction <- prediction_template %>%
  dplyr::distinct(Core_Fraction, .keep_all = FALSE) %>%
  dplyr::pull(Core_Fraction)


# List of genotypes.
full_geno_df <- "./data/derived/predictor_subsets/geno_df.RDS" %>%
  readRDS() %>%
  dplyr::mutate_at(
    .vars = vars(Extent, Material, Scenario),
    .funs = funs(substr(., start = 1, stop = 1)
  )) %>%
  tidyr::unite(
    col = Combi,
    c("Extent", "Material", "Scenario"),
    sep = "",
    remove = TRUE
  )



## -- ETA OBJECTS -------------------------------------------------------
# Load all ETA objects associated with the specified scenarios in the global
# environment.
# For the core set inbred lines we have plenty of ETA objects and loading them
# while processing predictions in parallel will likely cause a lot of
# unncecessary overhead.
eta_dir <- "./data/derived/predictor_subsets/eta/"
eta_uuid <- prediction_template %>%
  dplyr::distinct(ETA_UUID) %>%
  dplyr::pull(ETA_UUID)

eta_lst <- eta_uuid %>%
  purrr::map(~paste0(eta_dir, ., ".RDS")) %>%
  purrr::map(~readRDS(.)) %>%
  purrr::set_names(nm = eta_uuid)


if (grepl("CIC|FIC", x = pred_interval)) {

  # Load the data frame that specifies the main membership of genotypes to
  # to four putative ancestral populations based on a STRUCTURE analysis.
  # This is crucial for scenario C where we need to account for this structure
  # by adding a fixed effect term to each ETA object.
  structure_df <- "./data/derived/maizego/cluster_df_4pcs.RDS" %>%
    readRDS() %>%
    dplyr::select(G, main_cluster) %>%
    dplyr::mutate_at(vars(main_cluster), funs(as.factor))

  # Get the names of all genotypes that will be predicted in this particular
  # prediction run.
  # This is important to set up the fixed effect term for scenario C correctly.
  full_geno <- eta_lst %>%
    .[1] %>%
    purrr::pluck(1, function(x) rownames(x[[1]][[1]])) %>%
    tibble::as_data_frame() %>%
    dplyr::rename(G = "value") %>%
    dplyr::left_join(
      y = structure_df,
      by = "G"
    )


  # Create the fixed effect and add it to each ETA object in the list.
  fixeff <- list(fixed = list(
              ~factor(main_cluster),
              data = full_geno,
              model = "FIXED"
              )
            )

  eta_lst <- eta_lst %>%
    purrr::map(~purrr::prepend(., fixeff, before = 1))
}

eta_lst %>%
  list2env(., envir = globalenv())


## -- PHENOTYPIC DATA ----------------------------------------------------
if (material == "I") {

  pheno <- readRDS("./data/derived/maizego/tst_pheno_tibble.RDS")

} else if (material == "H") {

  pheno <- "./data/derived/uhoh/agro_tibble.RDS" %>%
    readRDS() %>%
    tidyr::separate(
      col = Genotype,
      into = c("Prefix", "Dent", "Flint"),
      sep = "_",
      remove = FALSE
    ) %>%
    dplyr::select(-Prefix)
}





# Sequence of all prediction scenarios that are scheduled in this particular
# run.
scenario_seq <- prediction_template %>%
  nrow() %>%
  seq_len()

if (nchar(runs) != 0){

  scenario_seq <- runs %>%
    base::strsplit(., split = "-") %>%
    purrr::flatten_chr() %>%
    as.integer() %>%
    max() %>%
    seq_len()
}


# Keep the time of the following computations.
start_time <- Sys.time()





scen_pred_lst <- parallel::mclapply(scenario_seq, FUN = function(i) {

  template_i <- prediction_template %>%
    dplyr::slice(i)

  if (material == "H") {

    # Load a table that connects the different variables to the data frame that
    # information on the corresponding TRN and TST genotypes.
    hybrid_scenario_df <- readRDS(
      "./data/derived/predictor_subsets/hybrid_scenario_df.RDS"
    ) %>%
    tibble::as_data_frame() %>%
    dplyr::mutate(ETA_UUID = paste0("eta_", ETA_UUID)) %>%
    dplyr::inner_join(
      y = template_i,
      by = c(
        "Core_Fraction",
        "Rnd_Level1",
        "Rnd_Level2",
        "Predictor",
        "ETA_UUID",
        "Combi",
        "TST_Geno",
        "TRN_Geno"
      )
    )

    # Load TRN and TST genotypes.
    geno_trn_df <- hybrid_scenario_df %>%
      dplyr::distinct(Hybrid_UUID, .keep_all = FALSE) %>%
      dplyr::pull(Hybrid_UUID) %>%
      purrr::map(
        ~paste0("./data/derived/predictor_subsets/hybrid_", ., ".RDS")
      ) %>%
      purrr::map(~readRDS(.)) %>%
      dplyr::bind_rows() %>%
      tibble::as_data_frame() %>%
      dplyr::mutate(ETA_UUID = paste0("eta_", ETA_UUID)) %>%
      dplyr::inner_join(
        y = hybrid_scenario_df %>% dplyr::select(-TST_Geno, -TRN_Geno),
        by = c(
          "Rnd_Level2",
          "Predictor",
          "ETA_UUID",
          "Combi",
          "Core_Fraction",
          "Rnd_Level1",
          "Hybrid_UUID"
        )
      )

    hyb_scen_geno_df <- full_geno_df %>%
      dplyr::right_join(
        y = template_i,
        by = "Combi"
      ) %>%
      dplyr::select(-TRN_Geno, -TST_Geno) %>%
      dplyr::rename(TRN_Geno = "G") %>%
      dplyr::distinct(TRN_Geno) %>%
      data.table::as.data.table() %>%
      data.table::setkey()

    # Separately for each combination of 'TST_Geno', 'Combi', 'Trait',
    # 'Core_Fraction', 'Rnd_Level1', 'Rnd_Level2', 'Predictor' and 'ETA_UUID'...
    # i) extract the training set
    # ii) define the hold-out set as its complement
    # iii) concatenate training set and hold-out set.
    geno_df <- geno_trn_df %>%
      data.table::data.table() %>%
      split(.,
        by = c(
          "TST_Geno",
          "Combi",
          "Trait",
          "Core_Fraction",
          "Rnd_Level1",
          "Rnd_Level2",
          "Predictor",
          "ETA_UUID"
        )
      ) %>%
      purrr::map(function(x) {
        x[, Set := "TRN", ]
        data.table::setkey(x)
        hold_out_trn <- data.table::copy(hyb_scen_geno_df[!x, on = "TRN_Geno"])
        tst_geno <- x[, (unique(TST_Geno)), ]
        unq_x_combi <- unique(
          x[, .(Rnd_Level2, Predictor, ETA_UUID, Combi, Core_Fraction, Rnd_Level1,
                Hybrid_UUID, Trait, Interval), ]
        )
        var_tbl <- unq_x_combi[rep(seq_len(.N), times = nrow(hold_out_trn)), ]
        hold_out <- cbind(hold_out_trn, var_tbl)
        hold_out[, `:=` (TST_Geno = tst_geno,
                         Set = "Hold_Out"), ]
        data.table::setkeyv(hold_out, cols = key(x))
        dplyr::bind_rows(x, hold_out)
      }) %>%
      data.table::rbindlist() %>%
      tibble::as_data_frame()

  } else if (material == "I" && all(core_fraction == "1")) {

    geno_trn_df <- "./data/derived/predictor_subsets/inbred_trn_df.RDS" %>%
      readRDS() %>%
      dplyr::mutate_at(
        .vars = vars(Extent, Material, Scenario),
        .funs = funs(substr(., start = 1, stop = 1)
      )) %>%
      tidyr::unite(
        col = Combi,
        c("Extent", "Material", "Scenario"),
        sep = "",
        remove = TRUE
      ) %>%
      tidyr::unite(Predictor, c("Pred1", "Pred2"), remove = TRUE) %>%
      dplyr::mutate_at(
        .vars = vars(Predictor),
        .funs = funs(gsub("_$", x = ., replacement = ""))
      ) %>%
      dplyr::rename(TRN_Geno = "G") %>%
      dplyr::mutate(
        TST_Geno = NA_character_,
        Rnd_Level2 = as.character(Rnd_Level2)
      ) %>%
      dplyr::inner_join(
        y = template_i %>% dplyr::select(-TST_Geno, -TRN_Geno),
        by = c(
          "Core_Fraction",
          "Rnd_Level1",
          "Rnd_Level2",
          "Predictor",
          "Combi"
        )
      ) %>%
      dplyr::mutate(Set = "TRN")

    # Generate hold-old genotypes, which are genotypes that are present but were
    # not selected for a given scenario.
    hold_out_df <- full_geno_df %>%
      dplyr::right_join(
        y = template_i,
        by = "Combi"
      ) %>%
      dplyr::anti_join(
        y = geno_trn_df,
        by = c(
          "G" = "TRN_Geno",
          "Combi",
          "Trait",
          "Core_Fraction",
          "Rnd_Level1",
          "Rnd_Level2",
          "Predictor",
          "ETA_UUID",
          "TST_Geno"
        )
      ) %>%
      dplyr::select(-TRN_Geno) %>%
      dplyr::rename(TRN_Geno = "G") %>%
      dplyr::mutate(Set = "Hold_Out")

    inbred_key <- c(
      "Combi",
      "Trait",
      "Core_Fraction",
      "Rnd_Level1",
      "Rnd_Level2",
      "Predictor",
      "ETA_UUID",
      "Interval"
    )

    trn_df <- dplyr::bind_rows(hold_out_df, geno_trn_df) %>%
      dplyr::select(-TST_Geno) %>%
      data.table::as.data.table() %>%
      data.table::setkeyv(., cols = inbred_key)

    tst_geno <- full_geno_df %>%
      dplyr::right_join(
        y = template_i,
        by = "Combi"
      ) %>%
      dplyr::select(-TST_Geno, -TRN_Geno) %>%
      dplyr::rename(TST_Geno = "G") %>%
      data.table::as.data.table() %>%
      data.table::setkeyv(., cols = inbred_key)

    geno_df <- merge(trn_df, tst_geno, all = TRUE, allow.cartesian = TRUE) %>%
      tibble::as_data_frame()

  } else if (material == "I" && any(core_fraction != "1")) {

    # For the inbred lines, generate leave-one-out cross-validation (LOOCV)
    # data frames with each genotypes occurring once as a TST genotype and
    # n times as a TRN genotype where n denotes the number of genotypes.
    # This will facilitate setting the correct genotypes to NA in the pertaining
    # data frame with phenotypic data and is crucial for predictions.
    geno_scen_df <- full_geno_df %>%
    dplyr::right_join(
      y = template_i,
      by = "Combi"
    )

    tst_geno <- geno_scen_df %>%
      dplyr::distinct(G) %>%
      dplyr::pull(G)

    # TRN_Geno, TST_Geno, Combi, Trait, Core_Fraction, Rnd_Level1, Rnd_Level2,
    # Predictor, ETA_UUID
    geno_df <- tst_geno %>%
      purrr::map(function(x) {
      trn_geno <- base::setdiff(tst_geno, x)
      tibble::data_frame(TRN_Geno = trn_geno, TST_Geno = x)
    }) %>%
    dplyr::bind_rows() %>%
    dplyr::full_join(
      y = geno_scen_df %>% dplyr::select(-TRN_Geno, -TST_Geno),
      by = c("TST_Geno" = "G")
    ) %>%
    dplyr::filter(complete.cases(.)) %>%
    dplyr::mutate(Set = "TRN")
  }


  # Ignoring training set genotypes, build a unique set of prediction runs to
  # loop over.
  unq_pred_template <- geno_df %>%
    dplyr::distinct(
      TST_Geno,
      Trait,
      Core_Fraction,
      Rnd_Level1,
      Rnd_Level2,
      Predictor,
      ETA_UUID,
      Interval
    )

  # Define a sequence with length equal to the number of unique test set
  # genotypes to build a loop from them where separate predictions are being
  # carried out for each fold of the loop.
  tst_seq <- unq_pred_template %>%
    nrow() %>%
    seq_len()





  ## -- BGLR ------------------------------------------------------------------
  tst_pred_lst <- lapply(tst_seq, FUN = function(j) {

    # Collect the current set of parameters.
    template_j <- unq_pred_template %>%
      dplyr::slice(j)

    # Collect the current ETA object for the prediction.
    eta_j <- template_j %>%
      dplyr::pull(ETA_UUID) %>%
      get(.)

    tst_j <- template_j %>%
      dplyr::select(TST_Geno, Trait) %>%
      dplyr::mutate(Set = "TST") %>%
      dplyr::rename(Genotype = "TST_Geno")

    # Combine the current test and the pre-built training and hold-out set.
    # Make sure to remove any genotype that was previously assigned to the
    # hold-out set when it is equal to the current test set.
    # Otherwise the data set `geno_i` would contain duplicate entries.
    geno_j <- geno_df %>%
      dplyr::inner_join(
        y = template_j,
        by = c(
          "TST_Geno",
          "Trait",
          "Core_Fraction",
          "Rnd_Level1",
          "Rnd_Level2",
          "Predictor",
          "ETA_UUID"
        )
      ) %>%
      dplyr::select(TRN_Geno, Trait, Set) %>%
      dplyr::rename(Genotype = "TRN_Geno") %>%
      dplyr::anti_join(
        y = tst_j %>% dplyr::mutate(Set = "Hold_Out"),
        by = c("Genotype", "Trait")
      ) %>%
      dplyr::bind_rows(tst_j)

    # Add missing genotypes as a third set (hold-out) and set their phenotypic
    # values to NA.
    pheno_j <- geno_j %>%
      dplyr::left_join(
        y = pheno,
        by = c("Genotype", "Trait")
      ) %>%
      dplyr::rename(Original_Observed = "Value") %>%
      dplyr::mutate(Observed_NA = dplyr::if_else(
        Set == "TRN",
        true = Original_Observed,
        false = NA_real_
        )
      )

    # Order phenotypic data according to the order of the ETA object.
    if (material == "I") {
      get_eta_rownames <- function(x) rownames(x[[1]])
      eta_geno_order <- eta_j %>%
        purrr::keep(names(.) != "fixed") %>%
        purrr::pluck(1, get_eta_rownames) %>%
        tibble::as_data_frame() %>%
        dplyr::rename(Genotype = "value")
      pheno_j <- eta_geno_order %>%
        dplyr::left_join(
          y = pheno_j,
          by = "Genotype"
        )
      stopifnot(
        identical(
          dplyr::pull(eta_geno_order, Genotype),
          dplyr::pull(pheno_j, Genotype)
        )
      )
    } else if (material == "H") {

      # In the case of multiple predictors we have "duplicates" of Dent and
      # Flint names, respectively. We need to filter a unique set of names
      # for each group in order to correctly align the order of genotypes in
      # ETA objects and in tables of phenotypic information.
      swap_logic <- purrr::compose(`!`, duplicated)
      duplicated_filter <- eta_j %>%
        purrr::map(.f = `[[`, "X") %>%
        purrr::map(~rownames(.)) %>%
        purrr::map(~tibble::as_data_frame(.)) %>%
        dplyr::bind_cols() %>%
        t() %>%
        swap_logic()

      eta_geno_order <- eta_j %>%
        .[duplicated_filter] %>%
        purrr::map(.f = `[[`, "X") %>%
        purrr::map(~rownames(.)) %>%
        purrr::map(~tibble::as_data_frame(.)) %>%
        dplyr::bind_cols() %>%
        dplyr::rename(Dent = "value", Flint = "value1") %>%
        tidyr::unite(col = Genotype, Dent, Flint, sep = "_") %>%
        dplyr::mutate(Genotype = paste0("DF_", Genotype))

      pheno_j <- eta_geno_order %>%
        dplyr::inner_join(y = pheno_j, by = "Genotype")

      stopifnot(
        identical(
          dplyr::pull(eta_geno_order, Genotype),
          dplyr::pull(pheno_j, Genotype)
        )
      )
    }


    # Run the BGLR model.
    mod_bglr <- BGLR::BGLR(
      y = dplyr::pull(pheno_j, Observed_NA),
      ETA = eta_j,
      nIter = iter,
      burnIn = iter / 2,
      saveAt = temp,
      verbose = FALSE
    )

    # Extract the index of the test set genotype from the phenotypic data in
    # order to extract the corresponding elements from the BGLR model ouptut.
    tst_idx <- pheno_j %>%
      dplyr::mutate(IDX = seq_len(n())) %>%
      dplyr::filter(Set == "TST") %>%
      dplyr::pull(IDX)

    # Save the results in a new data frame.
    yhat <- mod_bglr$yHat
    sd_yhat <- mod_bglr$SD.yHat
    res <- template_j
    res %>%
      dplyr::mutate(
        y = pheno_j %>%
          dplyr::slice(tst_idx) %>%
          dplyr::pull(Original_Observed),
        yhat = yhat[tst_idx],
        yhat_variance = sd_yhat[tst_idx] ^ 2
      )
  })

  tst_pred_df <- tst_pred_lst %>%
    dplyr::bind_rows()

}, mc.cores = use_cores, mc.preschedule = TRUE)



scen_pred_lst %>%
  dplyr::bind_rows() %>%
  saveRDS(
    object = .,
    file = paste0(
      "./data/derived/predictions/template_",
      job_id,
      ".RDS"
    )
  )

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

