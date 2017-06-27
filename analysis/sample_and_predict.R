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
               "tidyr", "purrr", "methods", "tibble", "lubridate",
               "readr")
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
  # Are you using inbred ("Inbred") data from the Yan-lab or hybrid ("Hybrid")
  # data from the UHOH-group?
  Sys.setenv("DATA_TYPE" = "Inbred")
  # Which agronomic trait do you want to evaluate? Leave blank if you want to
  # analyze all traits.
  Sys.setenv("TRAIT" = "100grainweight")
  # Number of iterations in BGLR()
  Sys.setenv("ITER" = "1000") 
  # Prediction model in BGLR()
  Sys.setenv("MODEL" = "BRR")
  # Algorithm to generate variance-covariance matrices.
  Sys.setenv("VCOV" = "RadenII")
  Sys.setenv("PI" = "0.5")
  Sys.setenv("PRIOR_PI_COUNT" = "10")
  # Main predictor. If 'Pred2' and 'Pred3' are empty, no imputation will take
  # place.
  Sys.setenv("PRED1" = "mrna")
  # If 'Pred3' is empty, 'Pred2' will be imputed via information from 'Pred1'.
  Sys.setenv("PRED2" = "")
  # Fraction of genotypes to be included in the core set.
  Sys.setenv("CORE_SET" = "1.0")
  # Which random core sample should be used (integer)
  Sys.setenv("RANDOM_SAMPLE" = "3")
  # Number of genotypes to predict (only for testing!)
  Sys.setenv("RUNS" = "1-2")
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
data_type <- as.character(Sys.getenv("DATA_TYPE"))
traits <- as.character(Sys.getenv("TRAIT"))
iter <- as.integer(Sys.getenv("ITER"))
pred_model <- as.character(Sys.getenv("MODEL"))
g_method <- as.character(Sys.getenv("VCOV"))
Pi <- as.numeric(Sys.getenv("PI"))
PriorPiCount <- as.integer(Sys.getenv("PRIOR_PI_COUNT"))
pred1 <- as.character(Sys.getenv("PRED1"))
pred2 <- as.character(Sys.getenv("PRED2"))
core_set <- as.character(Sys.getenv("CORE_SET"))
runs <- as.character(Sys.getenv("RUNS"))
temp <- paste0(as.character(Sys.getenv("TMP")), "/")


# Number of random samples per test set genotype.
core_nmb <- seq_len(10)

# -- PREDICTOR INPUT SPECIFICATION -------------------------------------
pred_combi <- list(pred1, pred2) %>%
  keep(nchar(.) != 0) %>%
  flatten_chr() %>%
  paste(., collapse = "_")

# Collect log files, which are indispensable for loading the data.
covered_genotypes <- "./data/derived/log_prepare_subsamples.txt" %>%
  read_tsv() %>%
  filter(
    Core_Set == core_set,
    Predictor == pred_combi,
    Random_Sample == core_nmb
  ) %>%
  pull(Data_Location) %>%
  gsub(x = ., pattern = "eta", replacement = "covered_genotypes") %>%
  readRDS()


pred_sets <- pred_combi %>%
  strsplit(split = "_") %>%
  flatten_chr()


## -- INPUT CHECKS AND UPDATES ------------------------------------------
if (isTRUE(core_set == "1.0" && length(pred_sets) > 1)) {
  stop("Imputation not possible if all predictors cover the same genotypes")
}

if (isTRUE(nchar(core_set) != 0 &&
           core_set != "1.0" &&
           data_type == "Hybrid")) {
  stop("Core sampling for hybrids is not yet supported")
}
  
if (isTRUE(data_type == "Hybrid")) {
  hybrid <- TRUE
} else hybrid <- FALSE


stopifnot(data_type %in% c("Inbred", "Hybrid"))
if (data_type == "Inbred") {
  if (any(grepl(pred_sets, pattern = "ped"))) {
    stop("No pedigree data for inbred lines")
  }
}

if (!all(pred_sets %in% c("ped", "snp", "mrna"))) {
  stop("Predictors are unknown")
}



## -- LOAD PREDICTOR DATA ------------------------------------------------
pheno <- readRDS("./data/derived/maizego/tst_pheno_tibble.RDS")
pheno <- pheno %>%
  filter(Genotype %in% covered_genotypes)





## -- AGRONOMIC DATA PREPARATION ------------------------------------------
# If only one trait was specified, use only that trait, otherwise select all
# traits.
if (nchar(traits) == 0 || length(traits) != 1) {
  traits <- pheno %>% select(Trait) %>% flatten_chr() %>% unique()
} else {
  traits <- traits 
}

# In the case of subsampling based on core samples, we may also have to subset
# the agronomic data afterwards.
# Therefore, we need to determine the largest set of genotypes covered by any
# selected predictor.
pheno_mat <- pheno %>%
  filter(Trait %in% traits,
         Genotype %in% covered_genotypes) %>%
  spread(key = Trait, value = Value) %>%
  as.data.frame() %>%
  remove_rownames() %>%
  column_to_rownames(var = "Genotype") %>%
  as.matrix()


## -- PREDICTION -----------------------------------------------------------
### Load pre-defined runs
pre_def_runs <- "./data/derived/predictor_subsets/loocv_samples.RDS" %>%
  readRDS()



# Define which runs shall be used, i.e., which genotypes shall be included as 
# test sets.
if (nchar(runs) != 0) {
  run_length <- runs %>%
    strsplit(., split = "-") %>%
    as_vector() %>%
    as.integer() %>%
    invoke(.f = seq, .x = ., by = 1) 
} else {
  run_length <- seq_len(nrow(pheno_mat))
}
param_df <- expand.grid(
  Tst_Geno = pre_def_runs %>% pull(TST_Geno) %>% unique() %>% .[run_length],
  Trait = traits,
  Iter = iter,
  Loocv_Run = pre_def_runs %>% pull(Iter) %>% unique(),
  Runs = length(run_length),
  Core_Number = core_nmb
)
param_df$Trait <- as.character(param_df$Trait)

if (isTRUE(data_type == "Hybrid")) {
  mother_idx <- 2
  father_idx <- 3
  split_char <- "_"
} else if (isTRUE(data_type == "Inbred")) {
  mother_idx <- NULL
  father_idx <- NULL
  split_char <- NULL
}


# Keep track of how long a job is running.
start_time <- Sys.time()
yhat_lst <- mclapply(seq_len(nrow(param_df)), FUN = function(i) {

  # Load pre-defined ETA objects.
  # Collect log files, which are indispensable for loading the data.
  eta <- "./data/derived/log_prepare_subsamples.txt" %>%
    read_tsv() %>%
    filter(
      Core_Set == core_set,
      Predictor == pred_combi,
      Random_Sample == core_nmb
    ) %>%
    pull(Data_Location) %>%
    readRDS()

  # Test genotype
  tst <- param_df %>%
    slice(i) %>%
    pull(Tst_Geno) %>%
    as.character()

  # Miscellaneous parameters
  loocv_run <- param_df %>%
    slice(i) %>%
    pull(Loocv_Run) %>%
    as.character()

  curr_core_nmb <- param_df %>%
    slice(i) %>%
    pull(Core_Number)

  # Training (sub)set
  trn <- pre_def_runs %>%
    filter(
      Iter == loocv_run,
      TST_Geno == tst
    ) %>%
  pull(TRN_Geno)

  # Test and training genotypes
  trn_tst <- c(tst, trn)

  # Operations required for input checks.
  eta <- lapply(eta, FUN = function(x) {
    sub_eta <- x$X
    x$X <- sub_eta[match(trn_tst, rownames(sub_eta)), , drop = FALSE]
    x
  })

  nrow_eta <- unique(vapply(eta, FUN = function(x) {
    nrow(x[["X"]])
  }, FUN.VALUE = integer(1)))

  pheno_mat <- pheno_mat[match(trn_tst, rownames(pheno_mat)), , drop = FALSE]

  # Ensure that the elements of ETA and the phenotyipic values are in the same
  # order.
  eta_rownms <- eta %>%
    map("X") %>%
    map(rownames) %>%
    reduce(intersect) 
  stopifnot(identical(rownames(pheno_mat), eta_rownms))

  # Trait
  trait <- param_df %>%
    slice(i) %>%
    pull(Trait) %>% 
    unique()

  # Set values in the test set as missing.
  y <- pheno_mat[, trait]
  y[tst] <- NA_real_
  stopifnot(sum(is.na(y)) == 1)

  # BGLR implementation
  res <- data.frame(
    Phenotype = NA_character_,
    Geno = NA_character_,
    Mother = NA_character_,
    Father = NA_character_,
    y = NA_complex_,
    yhat = NA_complex_,
    Loocv_run = loocv_run,
    Core_run = curr_core_nmb
  )

  # run the model (GBLUP)
  mod_BGLR <- BGLR::BGLR(
    y = y,
    ETA = eta,
    nIter = iter,
    burnIn = iter / 2,
    saveAt = temp,
    verbose = FALSE
  )

  idx <- which(names(mod_BGLR$yHat) == tst)
  yhat <- mod_BGLR$yHat
  sd_yhat <- mod_BGLR$SD.yHat

  # store results
  res$Phenotype <- trait
  res$Geno <- tst
  res$y <- pheno_mat[rownames(pheno_mat) == tst, trait]
  res$yhat <- yhat[idx]
  res$var_yhat <- sd_yhat[idx] ^ 2
  res
}, mc.cores = use_cores)
res <- rbindlist(yhat_lst)


# If not all three predictors were used, declare the missing ones as "none" and
# add them, which will facilitate further data analyses.
max_pred_number <- 2L
pred_nms <- c(pred_sets,
              rep("none", times = max_pred_number - length(pred_sets)))
res <- res %>%
  dplyr::rename(
    Trait = Phenotype,
    Dent = Mother,
    Flint = Father) %>%
  mutate(
    Data_Type = data_type,
    Pred1 = pred_nms[1],
    Pred2 = pred_nms[2],
    Core_Set = core_set,
    Job_ID = job_id,
    Runs = runs,
    Elapsed_Time = get_elapsed_time(start_time),
    Start_Time = as.character(start_time),
    Cores = use_cores
  ) %>%
  as.data.table()


log_file <- res %>%
  mutate(Iter = iter) %>%
  select(
    Job_ID,
    Pred1,
    Pred2,
    Data_Type,
    Core_Set,
    Elapsed_Time,
    Start_Time,
    Runs,
    Cores,
    Iter
  ) %>%
  unique()

# Reduce the size of the prediction object to the minimum possible size.
final <- res %>%
  select(
    Job_ID,
    Trait,
    Geno,
    Dent,
    Flint,
    y,
    yhat,
    Loocv_run,
    Core_run,
    var_yhat
  )


# Prediction results file
saveRDS(
  final, 
  file = paste0("./data/derived/predictions/", job_id, ".RDS"),
  compress = FALSE
)
                     
# Log file
log_location <- "./data/derived/uhoh_maizego_prediction_log.txt"
write.table(
  log_file,
  file = log_location,
  sep = "\t",
  row.names = FALSE,
  col.names = ifelse(
    file.exists(log_location),
    yes = FALSE,
    no = TRUE
  ),
  append = ifelse(file.exists(log_location), yes = TRUE, no = FALSE)
)

