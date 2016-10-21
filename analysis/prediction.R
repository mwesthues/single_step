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
pacman::p_load("BGLR","data.table", "parallel", "tidyverse")
pacman::p_load_gh("mwesthues/sspredr")


## --------------------------------------------------------------------------
## COLLECT AND CHECK ARGUMENTS 
if (isTRUE(interactive())) {
  # Number of cores
  Sys.setenv("MOAB_PROCCOUNT" = "4")
  Sys.setenv("TRAIT" = "GTM")
  # Number of iterations in BGLR()
  Sys.setenv("ITER" = "500")
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


## -- GENERATE BGLR() INPUT MATRICES --------------------------------------
eta_lst <- lapply(hetgrps, FUN = function(i) {
  grp_hyb <- vapply(strsplit(hybrid, split = "_"), FUN = "[[", i,
                    FUN.VALUE = character(1))
  # With imputation
  if (isTRUE(predictor == "mrna" && isTRUE(all_pheno))) {
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
    mutate(CV = "LOOCV",
           All_Pheno = all_pheno,
           Job_ID = job_id,
           Model = hypred_model,
           PI = Pi,
           PriorPiCount = PriorPiCount,
           Predictor = predictor,
           Elapsed_Time = elapsed_time) %>%
    as.data.table()
  
  log_file <- unique(res[, .(Job_ID, All_Pheno, Elapsed_Time, Trait, Iter, CV, 
                             Model, PI, PriorPiCount), ])
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

