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
pacman::p_load("data.table", "dplyr", "tidyr", "methods", "BGLR", "Matrix",
               "parallel", "matrixStats", "dtplyr", "tibble", "testthat")
pacman::p_load_gh("mwesthues/sspredr")


## --------------------------------------------------------------------------
## COLLECT AND CHECK ARGUMENTS 
if (isTRUE(interactive())) {
  # Number of cores
  Sys.setenv("MOAB_PROCCOUNT" = "4")
  Sys.setenv("TRAIT" = "")
  # Number of iterations in BGLR()
  Sys.setenv("ITER" = "10000")
  # Prediction model in BGLR()
  Sys.setenv("MODEL" = "BRR")
  # Algorithm to generate variance-covariance matrices.
  Sys.setenv("VCOV" = "RadenII")
  Sys.setenv("PI" = "0.5")
  Sys.setenv("PRIOR_PI_COUNT" = "10")
  # 'mrna' or 'snp' 
  Sys.setenv("PREDICTOR" = "mrna")
  # Which fraction of Dent genotypes shall be set to 'NA' for cross-validation?
  Sys.setenv("DENT_NA_FRACTION" = "0.05")
  # Which fraction of Flint genotypes shall be set to 'NA' for cross-validation?
  Sys.setenv("FLINT_NA_FRACTION" = "0.05")
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
predictor <- as.character(Sys.getenv("PREDICTOR"))
dent_na_frac <- as.numeric(Sys.getenv("DENT_NA_FRACTION" ))
flint_na_frac <- as.numeric(Sys.getenv("FLINT_NA_FRACTION" ))


# Input tests
poss_traits <- c("GTM", "GTS", "FETT", "RFA", "RPR", "STA", "XZ", "ADF")
if (isTRUE(nchar(init_traits) == 0)) {
  init_traits <- poss_traits
}
test_that("selected trait exists", {
  expect_true(all(init_traits %in% poss_traits))
})
test_that("selected model is part of BGLR", {
  expect_true(hypred_model %in% c("BRR", "BRR_Kernel", "BayesB", "BayesC"))
})
test_that("kernel method is defined", {
  expect_true(g_method %in% c("RadenI", "RadenII", "Zhang", "none"))
})
if (any(c(dent_na_frac, flint_na_frac) != 0)) {
  stopifnot(isTRUE(predictor == "mrna"))
}


## -- DATA SELECTION -----------------------------------------------------
genos <- readRDS("./data/processed/common_genotypes.RDS")
if (length(setdiff(genos$Dent$mrna, genos$Dent$snp)) == 0) {
  dent <- genos$Dent$mrna
} else stop("Not all genotypes have SNP data")
if (length(setdiff(genos$Flint$mrna, genos$Flint$snp)) == 0) {
  flint <- genos$Flint$mrna
} else stop("Not all genotypes have SNP data")

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


## -- PREDICTOR AND AGRONOMIC DATA PREPARATION ----------------------------
# Agronomic data
pheno <- readRDS("./data/processed/Pheno_stage2.RDS")
pheno <- pheno %>%
  filter(G %in% hybrid) %>%
  select(G, EST, Trait) %>%
  spread(key = Trait, value = EST) %>%
  select(-ADL) %>%
  as.data.frame %>%
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
}


## -- GENERATE BGLR() INPUT MATRICES --------------------------------------
eta_lst <- lapply(hetgrps, FUN = function(i) {
  grp_hyb <- vapply(strsplit(hybrid, split = "_"), FUN = "[[", i,
                    FUN.VALUE = character(1))
  if (isTRUE(predictor == "mrna" && 
             any(c(dent_na_frac, flint_na_frac) != 0))) {
    grp_snp <- snp[rownames(snp) %in% grp_hyb, ]
    grp_mrna <- mrna[rownames(mrna) %in% grp_hyb, ]
    eta <- impute_eta(x = grp_snp, y = grp_mrna, geno = grp_hyb,
                      bglr_model = hypred_model)
  } else if (isTRUE(predictor == "mrna")) {
    grp_mrna <- snp[rownames(mrna) %in% grp_hyb, ]
    eta <- complete_eta(x = mrna, geno = grp_hyb, bglr_model = hypred_model)
  } else if (isTRUE(predictor == "snp")) {
    grp_snp <- snp[rownames(snp) %in% grp_hyb, ]
    eta <- complete_eta(x = snp, geno = grp_hyb, bglr_model = hypred_model)
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
## PREDICTION
param_df <- expand.grid(Trait = init_traits,
                        Iter = init_iter,
                        Run = seq_len(nrow(pheno)))
param_df$Trait <- as.character(param_df$Trait)

# Keep track of how long a job is running.
start_time <- Sys.time()

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

# Determine how much time (hh:mm:ss format) has elapsed since script initiation.
get_elapsed_time <- function(start_time, tz = "CEST") {
  start_time <- as.POSIXct(start_time)
  dt <- difftime(Sys.time(), start_time, units = "secs")
  # Since you only want the H:M:S, we can ignore the date...
  # but you have to be careful about time-zone issues
  format(.POSIXct(dt, tz = tz), "%H:%M:%S")
}
elapsed_time <- get_elapsed_time(start_time)
 

res <- res %>%
  rename(Trait = Phenotype,
         Dent = Mother,
         Flint = Father) %>%
  mutate(Dent_NA = ifelse(Dent %in% dent_na_nms, yes = "yes", no = "no"),
         Flint_NA = ifelse(Flint %in% flint_na_nms, yes = "yes", no = "no"),
         CV = "LOOCV",
         Job_ID = job_id,
         BGLR_Model = hypred_model,
         PI = Pi,
         PriorPiCount = PriorPiCount,
         Predictor = predictor,
         Dent_NA_Fraction = dent_na_frac,
         Flint_NA_Fraction = flint_na_frac,
         Elapsed_Time = elapsed_time)

# Prediction results file
saveRDS(res, 
        file = paste0("./data/derived/predictions/", job_id, ".RDS"),
        compress = FALSE)
                     
# Log file
log_file <- unique(res[, .(Job_ID, Elapsed_Time, Trait, Iter, CV, BGLR_Model, 
                           PI, PriorPiCount, Dent_NA_Fraction,
                           Flint_NA_Fraction), ])
log_location <- "./data/derived/pred_log.txt"
write.table(log_file,
            file = log_location,
            sep = "\t", row.names = FALSE,
            col.names = ifelse(file.exists(log_location), yes = FALSE,
                               no = TRUE),
            append = ifelse(file.exists(log_location), yes = TRUE, no = FALSE))

