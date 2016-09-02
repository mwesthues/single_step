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
pacman::p_load("methods", "BGLR", "Matrix", "parallel", "data.table",
               "tidyr", "matrixStats", "testthat", "dtplyr", "dplyr", 
               "tibble")
devtools::install_github("mwesthues/sspredr", update = TRUE)
pacman::p_load_gh("mwesthues/sspredr")




## --------------------------------------------------------------------------
## COLLECT AND CHECK ARGUMENTS 
if (isTRUE(interactive())) {
  Sys.setenv("MOAB_PROCCOUNT" = "4")
  Sys.setenv("TRAIT" = "GTM")
  Sys.setenv("ITER" = "10000")
  Sys.setenv("MODEL" = "BRR")
  Sys.setenv("VCOV" = "RadenII")
  Sys.setenv("PI" = "0.5")
  Sys.setenv("PRIOR_PI_COUNT" = "10")
  Sys.setenv("IMPUTATION" = "TRUE")
  Sys.setenv("SPEED_TEST" = "FALSE")
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
imputation <- as.logical(Sys.getenv("IMPUTATION"))
speed_tst <- as.logical(Sys.getenv("SPEED_TEST"))


# Input tests
poss_traits <- c("GTM", "GTS", "ADL", "FETT", "RFA", "RPR", "STA", "XZ", "ADF")
test_that("selected trait exists", {
  expect_true(all(init_traits %in% poss_traits))
})
test_that("selected model is part of BGLR", {
  expect_true(hypred_model %in% c("BRR", "BRR_Kernel", "BayesB", "BayesC"))
})
test_that("kernel method is defined", {
  expect_true(g_method %in% c("RadenI", "RadenII", "Zhang", "none"))
})



## --------------------------------------------------------------------------
## LOAD CV SCHEME IF SPECIFIED
# If the cross-validation method is not LOOCV specify the CV scheme(s) and load
# it (them).
cv_name <- "cv1000_ps8081_trn=500_min_size=60_m=114_f=83.RDS"
cv <- readRDS(paste0("./data/processed/", cv_name))
runs <- sort(unique(cv$Run))
param_df <- expand.grid(Trait = init_traits,
                        Iter = init_iter,
                        Run = runs)
param_df$Trait <- as.character(param_df$Trait)



## --------------------------------------------------------------------------
## PHENOTYPIC DATA PREPARATION 
pheno <- readRDS("./data/processed/Pheno_stage2.RDS")
# Common genotypes
geno <- readRDS("./data/processed/common_genotypes.RDS")
pheno <- pheno %>%
  filter(G %in% geno$Hybrid) %>%
  select(G, EST, Trait) %>%
  spread(key = Trait, value = EST) %>%
  select(-ADL) %>%
  as.data.frame %>%
  column_to_rownames(var = "G") %>%
  as.matrix

# Position of Dent in rownames(pheno): 2
# Position of Flint in rownames(pheno): 3
hetgrps <- c(2, 3)


## --------------------------------------------------------------------------
## PREDICTOR AND AGRONOMIC DATA PREPARATION
# Load predictor data.
snp <- readRDS("./data/processed/snp_mat.RDS")

if (isTRUE(imputation)) {
  mrna <- readRDS("./data/processed/subset_mrna_blues.RDS")[["100%"]]
  mrna <- t(mrna)
  mrna <- mrna[grep("Exp|Sigma", x = rownames(mrna), invert = TRUE), ]
}

eta_lst <- lapply(hetgrps, FUN = function(i) {
  grp_hyb <- vapply(strsplit(geno$Hybrid, split = "_"), FUN = "[[", i,
                    FUN.VALUE = character(1))
  grp_snp <- snp[rownames(snp) %in% grp_hyb, ]

  if (isTRUE(imputation)) {
    grp_mrna <- mrna[rownames(mrna) %in% grp_hyb, ]
    eta <- impute_eta(x = grp_snp, y = grp_mrna, geno = grp_hyb,
                      bglr_model = hypred_model)
  } else {
     eta <- complete_eta(x = grp_snp, geno = grp_hyb, bglr_model = hypred_model)
  }
  eta
})
names(eta_lst) <- c("Dent", "Flint")
eta <- unlist(eta_lst, recursive = FALSE)
eta[] <- lapply(seq_along(eta), FUN = function(i) {
  dat <- eta[[i]]
  x <- dat[["X"]]
  rownames(x) <- geno$Hybrid
  dat[["X"]] <- x
  dat
})


## --------------------------------------------------------------------------
## CV1000
pred_lst <- mclapply(seq_len(nrow(param_df)), FUN = function(i) {
  run <- param_df[i, "Run"]
  trait <- param_df[i, "Trait"]
  iter <- param_df[i, "Iter"]
  pred <- run_cv(Pheno = pheno,
                 ETA = eta,
                 cv = cv,
                 mother_idx = 2,
                 father_idx = 3, 
                 split_char = "_",
                 trait = trait,
                 iter = iter,
                 speed_tst = FALSE,
                 run = run,
                 verbose = FALSE,
                 out_loc = "./tmp/")
  cbind(pred, Trait = trait, Iter = iter, CV = cv_name)
}, mc.cores = use_cores)
res <- rbindlist(pred_lst)
res[, `:=`(Job_ID = job_id,
           BGLR_Model = hypred_model,
           PI = Pi,
           PriorPiCount = PriorPiCount,
           Imputation = imputation), ]

# Prediction results file
saveRDS(res, 
        file = paste0("./data/processed/predictions/", job_id, ".RDS"),
        compress = FALSE)
                      
# Log file
log_file <- unique(res[, .(Job_ID, Trait, Iter, CV, BGLR_Model, PI, 
                           PriorPiCount, Imputation), ])
log_location <- "./data/processed/prediction_log/log_list.txt"
write.table(log_file,
            file = log_location,
            sep = "\t", row.names = FALSE,
            col.names = ifelse(file.exists(log_location), yes = FALSE,
                               no = TRUE),
            append = ifelse(file.exists(log_location), yes = TRUE, no = FALSE))

