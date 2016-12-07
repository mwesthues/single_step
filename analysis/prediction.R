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
  Sys.setenv("PRED3" = "mrna42")
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
pred_combi <- paste(c(pred1, pred2, pred3), collapse = "_")
if (isTRUE(transformation)) {
  pred_lst <- readRDS("./data/derived/transformed_pred_sub_list.RDS")
} else {
  pred_lst <- readRDS("./data/derived/pred_sub_list.RDS")
}

# Select the requested set of predictors.
pre_eta <- pred_lst %>%
  map(~keep(., names(.) == pred_combi))

# Select hybrids for whose parent lines at least one predictor has records.
hybrid <- pre_eta %>%
  at_depth(2, "geno") %>%
  flatten() %>%
  transpose() %>%
  map(paste, collapse = "_") %>%
  flatten_chr() %>%
  as_data_frame() %>%
  rename(G = value) %>%
  separate(G, into = c("Dent", "Flint"), sep = "_") %>%
  unite(col = Hybrid, Dent, Flint, sep = "_", remove = FALSE) %>%
  mutate(Hybrid = paste0("DF_", Hybrid))

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
      map(., ~complete_eta(x = .$pred1,
                           geno = .$geno,
                           as_kernel = TRUE,
                           is_pedigree = TRUE,
                           bglr_model = "BRR"))
  } else {
    eta_lst <- pred_lst %>%
      map(., ~complete_eta(x = .$pred1,
                           geno = .$geno,
                           as_kernel = TRUE,
                           is_pedigree = FALSE,
                           bglr_model = "BRR"))
  }

} else if (isTRUE(n_pred == 2)) {
  if (isTRUE(pedigree_first)) {
    eta_lst <- pred_lst %>%
      map(., ~impute_eta(x = .$pred1,
                         y = .$pred2,
                         geno = .$geno,
                         as_kernel = TRUE,
                         is_pedigree = TRUE,
                         bglr_model = "BRR"))
  } else {
    eta_lst <- pred_lst %>%
      map(., ~impute_eta(x = .$pred1,
                         y = .$pred2,
                         geno = .$geno,
                         as_kernel = TRUE,
                         is_pedigree = FALSE,
                         bglr_model = "BRR"))
  }

} else if (isTRUE(n_pred == 3)) {

  eta_lst <- pred_lst %>%
    map(., ~impute2(ped = .$pred1,
                    snp = .$pred2,
                    mrna = .$pred3,
                    geno = .$geno,
                    as_kernel = TRUE,
                    bglr_model = "BRR"))
}
eta <- flatten(eta_lst)

# Substitute hybrid names for the current names of the hybrid parents in the
# ETA objects. This is required for the cross-validation function to detect
# correct sorting of genotypes in the variance-covariance matrices with respect
# to the phenotypic traits.
hybrid_names <- eta_lst %>%
  map(~transpose(.)) %>%
  map("X") %>%
  at_depth(.depth = 2, rownames) %>%
  map(1) %>%
  reduce(.f = paste, sep = "_") %>%
  paste0("DF_", .)

eta <- lapply(eta, FUN = function(x) {
  rownames(x$X) <- hybrid_names
  x
})


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



start_time <- Sys.time()
# Keep track of how long a job is running.
yhat_lst <- mclapply(seq_len(nrow(param_df)), FUN = function(i) {
  run <- param_df[i, "Run"]
  trait <- param_df[i, "Trait"]
  iter <- as.integer(param_df[i, "Iter"])
  pred <- run_loocv(Pheno = pheno,
                    ETA = eta,
                    hybrid = TRUE,
                    mother_idx = 2,
                    father_idx = 3, 
                    split_char = "_",
                    trait = trait,
                    iter = iter,
                    speed_tst = FALSE,
                    run = run,
                    verbose = FALSE,
                    out_loc = "./tmp/")
  cbind(pred, Iter = iter)
}, mc.cores = use_cores)
res <- rbindlist(yhat_lst)
elapsed_time <- get_elapsed_time(start_time)
pred_nms <- pred_nms %>% 
  map_if(., .p = nchar(.) == 0, .f = function(x) {
    x <- "none"
    x
  }) %>%
  flatten_chr()
res <- res %>%
  rename(Trait = Phenotype,
         Dent = Mother,
         Flint = Father) %>%
  mutate(CV = "LOOCV",
         Pred1 = pred_nms[1],
         Pred2 = pred_nms[2],
         Pred3 = pred_nms[3],
         Job_ID = job_id,
         Runs = runs,
         Model = hypred_model,
         PI = Pi,
         PriorPiCount = PriorPiCount,
         Elapsed_Time = elapsed_time,
         Date = as.character(Sys.time()),
         Cores = use_cores) %>%
  as.data.table()

log_file <- unique(res[, .(Job_ID, Pred1, Pred2, Pred3, Elapsed_Time, Trait, 
                           Iter, CV, Model, PI, PriorPiCount, Date, Runs,
                           Cores), ])

# Reduce the size of the prediction object to the minimum possible size.
res[, `:=` (Iter = NULL, 
            Pred1 = NULL,
            Pred2 = NULL,
            Pred3 = NULL, 
            CV = NULL,
            Trait = NULL,
            Dent = NULL,
            Flint = NULL,
            Runs = NULL,
            Model = NULL,
            PI = NULL,
            PriorPiCount = NULL,
            Elapsed_Time = NULL,
            Date = NULL,
            Cores = NULL),
  ]
# Prediction results file
saveRDS(res, 
        file = paste0("./data/derived/predictions/", job_id, ".RDS"),
        compress = FALSE)
                     
# Log file
log_location <- "./data/derived/pred_log.txt"
write.table(log_file,
            file = log_location,
            sep = "\t", row.names = FALSE,
            col.names = ifelse(file.exists(log_location), yes = FALSE,
                               no = TRUE),
            append = ifelse(file.exists(log_location), yes = TRUE, no = FALSE))

