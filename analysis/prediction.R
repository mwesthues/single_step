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
pacman::p_load("BGLR", "caret", "corrplot", "Cubist", "data.table", "doMC", 
               "dplyr", "e1071", "earth", "elasticnet", "gbm", "kernlab", 
               "MASS", "Matrix", "matrixStats", "methods", "nnet", "parallel", 
               "party", "partykit", "pls", "randomForest", "rpart", "RWeka", 
               "testthat", "tibble", "tidyr")
pacman::p_load_gh("mwesthues/sspredr")


## --------------------------------------------------------------------------
## COLLECT AND CHECK ARGUMENTS 
if (isTRUE(interactive())) {
  # Number of cores
  Sys.setenv("MOAB_PROCCOUNT" = "4")
  # If TRUE, all agronomic data will be used, otherwise only agronomic data of
  # hybrids whose parents have both, mrna and snp data, will be used.
  Sys.setenv("ALL_PHENO" = "FALSE")
  Sys.setenv("TRAIT" = "GTM")
  # Number of iterations in BGLR()
  Sys.setenv("ITER" = "0")
  # Prediction model in BGLR()
  Sys.setenv("MODEL" = "pls")
  # Algorithm to generate variance-covariance matrices.
  Sys.setenv("VCOV" = "RadenII")
  Sys.setenv("PI" = "0")
  Sys.setenv("PRIOR_PI_COUNT" = "0")
  # 'mrna' or 'snp' 
  Sys.setenv("PREDICTOR" = "mrna")
  # Which fraction of Dent genotypes shall be set to 'NA' for cross-validation?
  Sys.setenv("DENT_NA_FRACTION" = "0.00")
  # Which fraction of Flint genotypes shall be set to 'NA' for cross-validation?
  Sys.setenv("FLINT_NA_FRACTION" = "0.00")
}

if (!interactive()) {
  job_id <- Sys.getenv("MOAB_JOBID")
} else job_id <- "interactive_00"

use_cores <- as.integer(Sys.getenv("MOAB_PROCCOUNT"))
all_pheno <- as.logical(Sys.getenv("ALL_PHENO"))
init_traits <- as.character(Sys.getenv("TRAIT"))
init_iter <- as.integer(Sys.getenv("ITER"))
hypred_model <- as.character(Sys.getenv("MODEL"))
g_method <- as.character(Sys.getenv("VCOV"))
Pi <- as.numeric(Sys.getenv("PI"))
PriorPiCount <- as.integer(Sys.getenv("PRIOR_PI_COUNT"))
predictor <- as.character(Sys.getenv("PREDICTOR"))
dent_na_frac <- as.numeric(Sys.getenv("DENT_NA_FRACTION" ))
flint_na_frac <- as.numeric(Sys.getenv("FLINT_NA_FRACTION" ))


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
if (any(c(dent_na_frac, flint_na_frac) != 0)) {
  stopifnot(isTRUE(predictor == "mrna"))
  if (isTRUE(all_pheno)) {
    stop(paste0("Missing data in parents not allowed when the full set of ", 
                "hybrids is used"))
  }
}


## -- DATA SELECTION -----------------------------------------------------
genos <- readRDS("./data/processed/common_genotypes.RDS")
if (isTRUE(all_pheno)) {
  dent <- genos$Dent$snp
  flint <- genos$Flint$snp
} else {
  if (length(setdiff(genos$Dent$mrna, genos$Dent$snp)) == 0) {
    dent <- genos$Dent$mrna
  } else stop("Not all genotypes have SNP data")
  if (length(setdiff(genos$Flint$mrna, genos$Flint$snp)) == 0) {
    flint <- genos$Flint$mrna
  } else stop("Not all genotypes have SNP data")
}

if (isTRUE(all_pheno)) {
  hybrid <- genos$Hybrid
} else {
  # Select hybrids for which both parents have genomic as well as transcriptomic
  # records.
  hybrid <- genos$Hybrid %>%
    as_data_frame() %>%
    separate(value, into = c("DF", "Dent", "Flint")) %>%
    mutate(Avail_Dent = ifelse(Dent %in% dent, yes = 1, no = 0),
           Avail_Flint = ifelse(Flint %in% flint, yes = 1, no = 0),
           Tested_Parents = Avail_Dent + Avail_Flint) %>%
    filter(Tested_Parents == 2) %>%
    unite(Hybrid, DF, Dent, Flint) %>%
    .$Hybrid
}


## -- PREDICTOR AND AGRONOMIC DATA PREPARATION ----------------------------
# Agronomic data
pheno <- readRDS("./data/processed/Pheno_stage2.RDS")
pheno <- pheno %>%
  as_data_frame() %>%
  filter(G %in% hybrid) %>%
  dplyr::select(G, EST, Trait) %>%
  spread(key = Trait, value = EST) %>%
  dplyr::select(-ADL) %>%
  as.data.frame %>%
  remove_rownames() %>%
  column_to_rownames(var = "G") %>%
  as.matrix

# Position of Dent in rownames(pheno): 2
# Position of Flint in rownames(pheno): 3
hetgrps <- c(2, 3)


# Endophenotypes
snp <- readRDS("./data/processed/snp_mat.RDS")
snp <- snp[rownames(snp) %in% c(dent, flint), ]
mrna <- readRDS("./data/processed/subset_mrna_blues.RDS")[["100%"]]
mrna <- t(mrna)
mrna <- mrna[grep("Exp|Sigma", x = rownames(mrna), invert = TRUE), ]
mrna <- mrna[rownames(mrna) %in% c(dent, flint), ]


# Set transcriptomic data to 'NA' on demand.
if (any(c(dent_na_frac, flint_na_frac) != 0)) {
  dent_na_nms <- sample(dent, size = ceiling(length(dent) * dent_na_frac), 
                        replace = FALSE)
  flint_na_nms <- sample(flint, size = ceiling(length(flint) * flint_na_frac), 
                         replace = FALSE)
  na_nms <- c(dent_na_nms, flint_na_nms)
  mrna <- mrna[!rownames(mrna) %in% na_nms, ]
} else {
  dent_na_nms <- vector(mode = "character", length = 0)
  flint_na_nms <- vector(mode = "character", length = 0)
}


## -- GENERATE BGLR() INPUT MATRICES --------------------------------------
eta_lst <- lapply(hetgrps, FUN = function(i) {
  grp_hyb <- vapply(strsplit(hybrid, split = "_"), FUN = "[[", i,
                    FUN.VALUE = character(1))
  # With imputation
  if (isTRUE(predictor == "mrna" && 
             (isTRUE(all_pheno) ||
              any(c(dent_na_frac, flint_na_frac) != 0)))) {
    grp_snp <- snp[rownames(snp) %in% grp_hyb, ]
    grp_mrna <- mrna[rownames(mrna) %in% grp_hyb, ]
    eta <- impute_eta(x = grp_snp, y = grp_mrna, geno = grp_hyb,
                      bglr_model = "BRR")
    # Without imputation 
  } else if (isTRUE(predictor == "mrna")) {
    grp_mrna <- snp[rownames(mrna) %in% grp_hyb, ]
    eta <- complete_eta(x = mrna, geno = grp_hyb, bglr_model = "BRR")
  } else if (isTRUE(predictor == "snp")) {
    grp_snp <- snp[rownames(snp) %in% grp_hyb, ]
    eta <- complete_eta(x = snp, geno = grp_hyb, bglr_model = "BRR")
  }
  eta
})
names(eta_lst) <- c("Dent", "Flint")
eta <- unlist(eta_lst, recursive = FALSE)
eta[] <- lapply(seq_along(eta), FUN = function(i) {
  dat <- eta[[i]]
  x <- dat[["X"]]
  rownames(x) <- hybrid
  dat[["X"]] <- x
  dat
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

param_df <- expand.grid(Trait = init_traits,
                        Iter = init_iter,
                        Run = seq_len(nrow(pheno)))
param_df$Trait <- as.character(param_df$Trait)
start_time <- Sys.time()

### CARET-PREDICTION
caret_mod_nms <- c("avNNet", "cubist", "earth", "enet", "gbm", "knn", "lm", 
                   "M5", "pls", "rf", "ridge", "rlm", "rpart2", "svmLinear", 
                   "svmPoly", "svmRadial", "treebag")

if (hypred_model %in% caret_mod_nms) {
  registerDoMC(use_cores)
  # Concatenate all individual predictor matrices to one big predictor matrix.
  pred_lst <- eta %>%
    unlist(recursive = FALSE) %>%
    .[grep(".model", x = names(.), invert = TRUE)] %>%
    .[grep("1", x = names(.), invert = TRUE)]
  # In order to avoid duplicated column names (which would otherwise be the 
  # case when concatenating Dent and Flint predictors), add the first letter
  # of the heterotic group to the predictor name.
  pred_lst[] <- lapply(seq_along(pred_lst), FUN = function(i) {
    x <- pred_lst[[i]]
    element_nm <- substr(names(pred_lst[i]), start = 1, stop = 1)
    colnames(x) <- paste0(element_nm, "_", colnames(x))
    x
  })
  pred_mat <- pred_lst %>%
    do.call(cbind, .)
  if (anyDuplicated(colnames(pred_mat))) stop("Duplicated predictor names")
  
  trait_vec <- param_df %>%
    dplyr::select(Trait) %>%
    unique.data.frame() %>%
    mutate(Trait = as.character(Trait)) %>%
    .$Trait
  
  
  res <- lapply(seq_len(nrow(param_df)), FUN = function(i) {
    trait <- trait_vec[i]
    y <- pheno[, trait]
    stopifnot(all.equal(names(y), rownames(pred_mat)))
    
    if (isTRUE(hypred_model == "pls")) {
      tuned_mod <- train(pred_mat, y,
                         method = "pls",
                         preProc = c("BoxCox", "center", "scale"),
                         tuneLength = 20,
                         trControl = trainControl(method = "LOOCV"))
    }
    if (isTRUE(hypred_model == "ridge")) {
      tune_grid <- data.frame(.lambda = seq(0, to = 1, length = 15))
      tuned_mod <- train(pred_mat, y,
                         method = "ridge",
                         preProc = c("BoxCox", "center", "scale"),
                         tuneGrid = tune_grid,
                         trControl = trainControl(method = "LOOCV"))
    }
    if (isTRUE(hypred_model == "enet")) {
      tune_grid <- data.frame(.lambda = seq(0, 0.01, 0.1),
                              .fraction = seq(0.05, to = 1, length = 20))
      tuned_mod <- train(pred_mat, y,
                         method = "enet",
                         preProc = c("BoxCox", "center", "scale"),
                         tuneGrid = tune_grid,
                         trControl = trainControl(method = "LOOCV"))
    }   
    if (isTRUE(hypred_model == "earth")) {
      tune_grid <- expand.grid(.degree = seq_len(2),
                               .nprune = 2:38)
      tuned_mod <- train(pred_mat, y,
                         method = "earth",
                         preProc = c("BoxCox"),
                         tuneGrid = tune_grid,
                         trControl = trainControl(method = "LOOCV"))
    }   
    if (isTRUE(hypred_model == "svmRadial")) {
     tuned_mod <- train(pred_mat, y,
                        method = "svmRadial",
                        preProc = c("center", "scale"),
                        tuneLength = 5,
                        trControl = trainControl(method = "LOOCV"))
    }
    if (isTRUE(hypred_model == "knn")) {
      tune_grid <- data.frame(.k = seq_len(20))
      tuned_mod <- train(pred_mat, y,
                         method = "knn",
                         preProc = c("BoxCox", "center", "scale", "nzv"),
                         tuneGrid = tune_grid,
                         trControl = trainControl(method = "LOOCV"))
    }   
    if (isTRUE(hypred_model == "M5")) {
      tuned_mod <- train(pred_mat, y,
                         method = "M5",
                         control = Weka_control(M = 10),
                         trControl = trainControl(method = "LOOCV"))
    }   
    if (isTRUE(hypred_model == "treebag")) {
      tuned_mod <- train(pred_mat, y,
                         method = "treebag",
                         nbagg = 50,
                         trControl = trainControl(method = "LOOCV"))
    }   
    if (isTRUE(hypred_model == "rf")) {
      tune_grid <- data.frame(mtry = floor(seq(10, to = ncol(training), 
                                               length = 10)))
      tuned_mod <- train(pred_mat, y,
                         method = "rf",
                         tuneGrid = tune_grid,
                         ntree = 1000,
                         importance = TRUE,
                         trControl = trainControl(method = "LOOCV"))
    }   
    if (isTRUE(hypred_model == "gbm")) {
      tune_grid <- expand.grid(interaction.depth = seq(1, to = 7, by = 2),
                               n.trees = seq(100, to = 1000, by = 50),
                               shrinkage = c(0.01, 0.1))
      tuned_mod <- train(pred_mat, y,
                         method = "gbm",
                         tuneGrid = tune_grid,
                         verbose = FALSE,
                         trControl = trainControl(method = "LOOCV"))
    }   
    if (isTRUE(hypred_model == "cubist")) {
      tune_grid <- expand.grid(committees = c(1:10, 20, 50, 75, 100),
                               neighbors = c(0, 1, 5, 9))
      tuned_mod <- train(pred_mat, y,
                         method = "cubist",
                         tuneGrid = tune_grid,
                         trControl = trainControl(method = "LOOCV"))
    }   
    tuned_mod
  })
  names(res) <- trait_vec
 
  elapsed_time <- get_elapsed_time(start_time)
  log_file <- expand.grid(Job_ID = job_id, Elapsed_Time = elapsed_time,
                          Trait = trait_vec, Iter = 0, CV = "LOOCV",
                          Model = hypred_model, PI = 0, PriorPiCount = 0,
                          Dent_NA_Fraction = 0, Flint_NA_Fraction = 0)
}



# Keep track of how long a job is running.
if (!hypred_model %in% caret_mod_nms) {
  pred_lst <- mclapply(seq_len(nrow(param_df)), FUN = function(i) {
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
  res <- rbindlist(pred_lst)
  elapsed_time <- get_elapsed_time(start_time)
  res <- res %>%
    rename(Trait = Phenotype,
           Dent = Mother,
           Flint = Father) %>%
    mutate(Dent_NA = ifelse(Dent %in% dent_na_nms, yes = "yes", no = "no"),
           Flint_NA = ifelse(Flint %in% flint_na_nms, yes = "yes", no = "no"),
           CV = "LOOCV",
           Job_ID = job_id,
           Model = hypred_model,
           PI = Pi,
           PriorPiCount = PriorPiCount,
           Predictor = predictor,
           Dent_NA_Fraction = dent_na_frac,
           Flint_NA_Fraction = flint_na_frac,
           Elapsed_Time = elapsed_time)
  
  log_file <- unique(res[, .(Job_ID, Elapsed_Time, Trait, Iter, CV, Model, 
                             PI, PriorPiCount, Dent_NA_Fraction,
                             Flint_NA_Fraction), ])
}


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

