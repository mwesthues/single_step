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

# Load function to save session info on the software used for the analysis in
# this script.
source("./software/session_info.R")

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
  # Were the data transformed (TRUE) or not?
  Sys.setenv("TRANSFORMATION" = "TRUE")
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

# Save session info
write_session_info(directory = "./data/derived/session_info/",
                   job_id = job_id)

use_cores <- as.integer(Sys.getenv("MOAB_PROCCOUNT"))
init_traits <- as.character(Sys.getenv("TRAIT"))
init_iter <- as.integer(Sys.getenv("ITER"))
hypred_model <- as.character(Sys.getenv("MODEL"))
g_method <- as.character(Sys.getenv("VCOV"))
Pi <- as.numeric(Sys.getenv("PI"))
PriorPiCount <- as.integer(Sys.getenv("PRIOR_PI_COUNT"))
transformation <- as.logical(Sys.getenv("TRANSFORMATION"))
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
pred_combi <- list(pred1, pred2, pred3) %>%
  keep(nchar(.) != 0) %>%
  flatten_chr() %>%
  paste(., collapse = "_")

if (isTRUE(transformation)) {
  if (nchar(pred2) == 0 && nchar(pred3) == 0 && grepl("ped", x = pred1)) {
    stop("No transformation was applied to 'ped100'. Choose a different set!")
  }
}

if (isTRUE(transformation)) {
  pred_lst <- readRDS("./data/derived/transformed_pred_sub_list.RDS")
} else {
  pred_lst <- readRDS("./data/derived/pred_sub_list.RDS")
}

# Select the requested set of predictors.
pre_eta <- pred_lst %>%
  map(~keep(., names(.) == pred_combi)) %>%
  at_depth(.depth = 2, .f = function(x) {
    x$as_kernel <- TRUE
    x$bglr_model <- hypred_model
    x 
  }) %>% 
  flatten()
names(pre_eta) <- names(pred_lst)

# Select hybrids for whose parent lines at least one predictor has records.
hybrid <- pre_eta %>%
  map("geno") %>%
  transpose() %>%
  map(paste, collapse = "_") %>%
  flatten_chr() %>%
  as_data_frame() %>%
  rename(G = value) %>%
  separate(G, into = c("Dent", "Flint"), sep = "_") %>%
  unite(col = Hybrid, Dent, Flint, sep = "_", remove = FALSE) %>%
  mutate(Hybrid = paste0("DF_", Hybrid))

## -- AGRONOMIC DATA PREPARATION ----------------------------
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
# Determine the number of predictors to decide later on, which function to
# call for the set-up of the kernels.
raw_args <- pre_eta %>%
  map(names) %>%
  reduce(intersect)
pred_number <- raw_args %>%
  .[grep("ped|snp|mrna", x = .)] %>%
  length()

# If all predictors are being used, it is known that pedigree data are involved
# and no measures need to be taken.
# If only a subset of predigrees will be used, it needs to be determined
# whether pedigree data are involved or not, since subsequent functions need to
# 'know' if the predictor matrices need to be scaled or not.
if (pred_number != 3) {
  if (any(grepl("ped", x = raw_args))) {
    pre_eta <- pre_eta %>%
      map(.f = function(x) {
        x$is_pedigree <- TRUE
        x
      })
  } else {
    pre_eta <- pre_eta %>%
      map(.f = function(x) {
        x$is_pedigree <- FALSE
        x
      })
  }
}

# Depending on the number of predictors included in the set, choose the correct
# function for the set-up of the BGLR kernels.
eta_fun <- c(
  "complete_eta", "impute_eta", "impute2"
  ) %>%
  .[pred_number]

# Collect the formal arguments from the selected function and then rename the
# corresponding predictors in the 'pre_eta' object so that the call of the 
# kernel building function (e.g. "impute_eta") can be generalized.
eta_fun_args <- eta_fun %>%
  formals() %>%
  names()
eta_pred_nms <- eta_fun_args %>%
  .[!grepl("as_kernel|is_pedigree|geno|bglr_model", x = .)]
current_eta_names <- pre_eta %>%
  map(names) %>%
  reduce(intersect)
current_eta_names[seq_len(pred_number)] <- eta_pred_nms
pre_eta <- pre_eta %>%
  map(function(x) {
    names(x) <- current_eta_names
    x
  })
# Build the kernel for BGLR.
eta <- invoke_map(get(eta_fun), pre_eta) %>%
  flatten()

# Substitute hybrid names for the current names of the hybrid parents in the
# ETA objects. This is required for the cross-validation function to detect
# correct sorting of genotypes in the variance-covariance matrices with respect
# to the phenotypic traits.
eta <- lapply(eta, FUN = function(x) {
  rownames(x$X) <- hybrid$Hybrid
  x
})


## --------------------------------------------------------------------------
#
## PREDICTION
# Determine how much time (hh:mm:ss format) has elapsed since script initiation.
get_elapsed_time <- function(start_time, tz = "CEST") {
  start_time <- as.POSIXct(start_time)
  sec <- Sys.time() %>%
    difftime(time1 = ., time2 = start_time, units = "secs") %>%
    lubridate::seconds_to_period()
  paste0(sprintf("%02d", c(day(sec), hour(sec), minute(sec), second(sec))),
         collapse = ":")
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

# If not all three predictors were used, declare the missing ones as "none" and
# add them, which will facilitate further data analyses.
max_pred_number <- 3L
pred_nms <- pred_combi %>% 
  str_split(., pattern = "_") %>%
  flatten_chr()
pred_nms <- c(pred_nms,
              rep("none", times = max_pred_number - length(pred_nms)))
res <- res %>%
  rename(Trait = Phenotype,
         Dent = Mother,
         Flint = Father) %>%
  mutate(CV = "LOOCV",
         Pred1 = pred_nms[1],
         Pred2 = pred_nms[2],
         Pred3 = pred_nms[3],
         Transformation = transformation,
         Job_ID = job_id,
         Runs = runs,
         Model = hypred_model,
         PI = Pi,
         PriorPiCount = PriorPiCount,
         Elapsed_Time = get_elapsed_time(start_time),
         Start_Time = as.character(start_time),
         Cores = use_cores) %>%
  as.data.table()

log_file <- unique(res[, .(Job_ID, Pred1, Pred2, Pred3, Transformation, 
                           Elapsed_Time, Trait, Iter, CV, Model, PI,
                           PriorPiCount, Start_Time, Runs, Cores), ])

# Reduce the size of the prediction object to the minimum possible size.
res[, `:=` (Iter = NULL, 
            Pred1 = NULL,
            Pred2 = NULL,
            Pred3 = NULL, 
            Transformation = NULL,
            CV = NULL,
            Trait = NULL,
            Dent = NULL,
            Flint = NULL,
            Runs = NULL,
            Model = NULL,
            PI = NULL,
            PriorPiCount = NULL,
            Elapsed_Time = NULL,
            Start_Time = NULL,
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

