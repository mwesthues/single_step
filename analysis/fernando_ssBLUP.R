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
               "tidyr", "matrixStats", "testthat", "dtplyr", "dplyr")
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
  Sys.setenv("ONLY_PROFILED" = "TRUE")
}

use_cores <- as.integer(Sys.getenv("MOAB_PROCCOUNT"))
init_traits <- as.character(Sys.getenv("TRAIT"))
init_iter <- as.integer(Sys.getenv("ITER"))
hypred_model <- as.character(Sys.getenv("MODEL"))
g_method <- as.character(Sys.getenv("VCOV"))
Pi <- as.numeric(Sys.getenv("PI"))
PriorPiCount <- as.integer(Sys.getenv("PRIOR_PI_COUNT"))
imputation <- as.logical(Sys.getenv("IMPUTATION"))
speed_tst <- as.logical(Sys.getenv("SPEED_TEST"))
only_profiled <- as.logical(Sys.getenv("ONLY_PROFILED"))


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

# Print start time to a log file
start_time <- Sys.time()
print(start_time)



## --------------------------------------------------------------------------
## LOAD CV SCHEME IF SPECIFIED
# If the cross-validation method is not LOOCV specify the CV scheme(s) and load
# it (them).
cv <- readRDS(paste0("./data/processed/cv1000_ps8081_trn=500_min_", 
                     "size=60_m=114_f=83.RDS"))
param_df <- expand.grid(Phenotype = init_traits, 
                        Iter = init_iter)
for (i in seq_len(ncol(param_df))) param_df[, i] <- as.character(param_df[, i])
numParam <- nrow(param_df)



## --------------------------------------------------------------------------
## PHENOTYPIC DATA PREPARATION 
pheno <- readRDS("./data/processed/Pheno_stage2.RDS")
# Common genotypes
geno <- readRDS("./data/processed/common_genotypes.RDS")
pheno <- pheno[G %in% geno$Hybrid, .(G, EST, Trait), ]
pheno <- pheno %>% 
  spread(Trait, value = EST) %>%
  select(-ADL)
pheno <- data.frame(pheno)
rownames(pheno) <- pheno$G
pheno$G <- NULL
pheno <- as.matrix(pheno)


## --------------------------------------------------------------------------
## PREDICTOR AND AGRONOMIC DATA PREPARATION
# Load predictor data.
snp <- readRDS("./data/processed/snp_mat.RDS")
mrna <- readRDS("./data/processed/subset_mrna_blues.RDS")[["100%"]]
mrna <- t(mrna)
mrna <- mrna[grep("Exp|Sigma", x = rownames(mrna), invert = TRUE), ]
hetgrps <- c(2, 3)
eta_lst <- lapply(hetgrps, FUN = function(i) {
  grp_hyb <- vapply(strsplit(geno$Hybrid, split = "_"), FUN = "[[", i,
                    FUN.VALUE = character(1))
  grp_snp <- snp[rownames(snp) %in% grp_hyb, ]
  grp_mrna <- mrna[rownames(mrna) %in% grp_hyb, ]
  eta <- impute_eta(x = grp_snp, y = grp_mrna, geno = grp_hyb,
                    bglr_model = hypred_model)
  eta
  #unlist(eta, recursive = FALSE)
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
str(eta)


# Dent
hist(eta[[2]][["X"]][upper.tri(eta[[2]][["X"]], diag = FALSE)])
hist(eta[[3]][["X"]][upper.tri(eta[[3]][["X"]], diag = FALSE)])

hist(eta[[5]][["X"]][upper.tri(eta[[5]][["X"]], diag = FALSE)])
hist(eta[[6]][["X"]][upper.tri(eta[[6]][["X"]], diag = FALSE)])


## --------------------------------------------------------------------------
## CV1000
runs <- unique(cv$Run)
runs <- seq_len(16)
pred_lst <- mclapply(seq_along(runs), FUN = function(i) {
  run <- runs[i]
  pred <- run_cv(Pheno = pheno,
                 ETA = eta,
                 cv = cv,
                 mother_idx = 2,
                 father_idx = 3, 
                 split_char = "_",
                 trait = "GTM",
                 iter = 2500L,
                 speed_tst = FALSE,
                 run = run,
                 verbose = FALSE,
                 out_loc = "./tmp/")
  pred
}, mc.cores = use_cores)
rbindlist(pred_lst)













# If specified, keep only the intersect of parents that have both, genomic
# as well as transcriptomic records.
if (isTRUE(only_profiled)) {
  comgeno <- Reduce("intersect", lapply(endo_lst, FUN = rownames))
  comgeno <- intersect(comgeno, c(dent, flint))
  endo_lst[] <- lapply(seq_along(endo_lst), FUN = function(i) {
    dat <- endo_lst[[i]]
    dat[match(comgeno, rownames(dat)), ]
  })
}

# If specified, use only SNP data for prediction.
if (!isTRUE(imputation)) {
  endo_lst <- endo_lst[["A"]]
}

# Build a kernel from the genomic data
hetgrps <- c("Dent", "Flint")
if (isTRUE(imputation)) {
  ETA <- lapply(hetgrps, FUN = function(hetgrp) {
    impute_ETA(snp = snp, mrna = mrna, grp_geno = grp_geno, hetgrp = hetgrp,
               bglr_model = hypred_model)
  })
} else {
  ETA <- lapply(hetgrps, FUN = function(hetgrp) {
    create_snp_ETA(snp = snp, grp_geno = grp_geno, hetgrp = hetgrp, 
                   bglr_model = hypred_model)
  })
}
ETA <- unlist(ETA, recursive = FALSE)




## ---------------------------------------------------------------------------
## SPEED TEST
if (isTRUE(speed_tst)) {
  trait <- as.character(param_df[1, "Phenotype"])
  iter <- as.integer(param_df[1, "Iter"])
  burnin <- iter / 2
  runs <- 1L
  # Operations required for input checks.
  nrow_eta <- unique(vapply(ETA, FUN = function(x) {
    nrow(x[["X"]])
  }, FUN.VALUE = integer(1)))

  if (!exists("cv_method")) {
    speed <- run_bglr_loocv(Pheno = y_mat, ETA = ETA, trait = trait, 
                            iter = iter, ncores = 1, speed_tst = TRUE,
                            out = "./tmp/")
    speed_df <- data.frame(Elapsed = speed[["elapsed"]],
                           Trait = trait,
                           Iter = iter,
                           CV_Method = cv_method,
                           Imputation = imputation,
                           Only_Profiled = only_profiled,
                           VCOV = g_method,
                           Pi = Pi,
                           PriorPiCount = PriorPiCount,
                           Model = hypred_model)
    if (exists("./data/derived/speed_tests.txt")) {
      write.table(speed_df, file = "./data/derived/speed_tests.txt",
                  append = TRUE, sep = "\t", col.names = FALSE)
    } else {
      write.table(speed_df, file = "./data/derived/speed_tests.txt",
                  append = FALSE, sep = "\t", col.names = FALSE)
    }
  } else if (isTRUE(cv_method == "CV800")) {
    cur_cv_nm <- as.character(param_df[1, "CV_Scheme"])
    cv <- cv_lst[[cur_cv_nm]]
    # Input testing.
    test_that("BGLR input is correct", {
      expect_is(cv, "data.frame")
      expect_named(cv, c("Sample_ID", "Run", "Set"))
      expect_type(trait, "character")
      expect_length(trait, 1)
      expect_type(iter, "integer")
      expect_is(y_mat, "matrix")
      expect_gt(nrow(y_mat), ncol(y_mat))
      expect_is(ETA, "list")
      expect_match(rownames(y_mat), regexp = "DF_[0-9]*_[0-9]*")
      expect_length(nrow_eta, 1)
      expect_identical(nrow_eta, expected = nrow(y_mat))
    })
    speed <- run_bglr_cv800(Pheno = y_mat, ETA = ETA, trait = trait,
                            iter = iter, burnin = burnin, 
                            cv = cv, run = runs, speed_tst = TRUE,
                            out = "./tmp/")
    speed_df <- data.frame(Elapsed = speed[["elapsed"]],
                           Trait = trait,
                           Iter = iter,
                           CV_Method = cv_method,
                           Imputation = imputation,
                           Only_Profiled = only_profiled,
                           VCOV = g_method,
                           Pi = Pi,
                           PriorPiCount = PriorPiCount,
                           Model = hypred_model)
    if (file.exists("./data/derived/speed_tests.txt")) {
      write.table(speed_df, file = "./data/derived/speed_tests.txt",
                  append = TRUE, sep = "\t", col.names = FALSE, 
                  row.names = FALSE)
    } else {
      write.table(speed_df, file = "./data/derived/speed_tests.txt",
                  append = FALSE, sep = "\t", col.names = TRUE,
                  row.names = FALSE)
    }
  }
}



use_cores <- as.integer(Sys.getenv("MOAB_PROCCOUNT"))
init_traits <- as.character(Sys.getenv("TRAIT"))
init_iter <- as.integer(Sys.getenv("ITER"))
hypred_model <- as.character(Sys.getenv("MODEL"))
g_method <- as.character(Sys.getenv("VCOV"))
Pi <- as.numeric(Sys.getenv("PI"))
PriorPiCount <- as.integer(Sys.getenv("PRIOR_PI_COUNT"))
imputation <- as.logical(Sys.getenv("IMPUTATION"))
cv_method <- as.character(Sys.getenv("CV_METHOD"))
speed_tst <- as.logical(Sys.getenv("SPEED_TEST"))
if (isTRUE(cv_method == "CV800")) {
  cv_scheme <- as.character(Sys.getenv("CV_SCHEME"))
}
only_profiled <- as.logical(Sys.getenv("ONLY_PROFILED"))








## ---------------------------------------------------------------------------
## COMPUTATION OF PREDICTIVE ABILITY
if (!isTRUE(speed_tst)) {
  if (!exists("cv_method")) {
    full_lst <- vector(mode = "list", length = nrow(param_df))
    for (i in seq_len(nrow(param_df))) {
      trait <- as.character(param_df[i, "Phenotype"])
      iter <- as.integer(param_df[i, "Iter"])
      full_lst[[i]] <- mclapply(seq_len(nrow(y_mat)), FUN = function(run) {
        run_bglr_loocv(Pheno = y_mat, ETA = ETA, trait = trait, iter = iter,
                       ncores = use_cores, speed_tst = FALSE, run = run,
                       out = "./tmp/")
      })
    }
    DT <- rbindlist(full_lst)
  } else if (isTRUE(cv_method == "CV800")) {
  
    #### Run CV800
    # Operations required for input checks.
    nrow_eta <- unique(vapply(ETA, FUN = function(x) {
      nrow(x[["X"]])
    }, FUN.VALUE = integer(1)))
  
    full_lst <- vector(mode = "list", length = nrow(param_df))
    for (i in seq_len(nrow(param_df))) {
      trait <- as.character(param_df[i, "Phenotype"])
      iter <- as.integer(param_df[i, "Iter"])
      cur_cv_nm <- as.character(param_df[i, "CV_Scheme"])
      burnin <- iter / 2
      cv <- cv_lst[[cur_cv_nm]]
      runs <- length(unique(cv$Run))
  
      # Input testing.
      test_that("BGLR input is correct", {
        expect_equal(runs, 800)
        expect_is(cv, "data.frame")
        expect_named(cv, c("Sample_ID", "Run", "Set"))
        expect_type(trait, "character")
        expect_length(trait, 1)
        expect_type(iter, "integer")
        expect_is(y_mat, "matrix")
        expect_gt(nrow(y_mat), ncol(y_mat))
        expect_is(ETA, "list")
        expect_match(rownames(y_mat), regexp = "DF_[0-9]*_[0-9]*")
        expect_length(nrow_eta, 1)
        expect_identical(nrow_eta, expected = nrow(y_mat))
      })
      bglr_out <- mclapply(seq_len(runs), FUN = function(run) {
      run_bglr_cv800(Pheno = y_mat, ETA = ETA, trait = trait,
                     iter = iter, burnin = burnin, 
                     cv = cv, run = run, speed_tst = FALSE, out = "./tmp/")
      }, mc.cores = use_cores)
      scen_DT <- rbindlist(bglr_out)
      scen_DT[, `:=`(CV_Scheme = cur_cv_nm,
                     Trait = trait,
                     Iter = iter,
                     Model = hypred_model,
                     VCOV = g_method,
                     Pi = Pi, 
                     PriorPiCount = PriorPiCount,
                     Imputed = imputation), ]
      full_lst[[i]] <- scen_DT
    }
    DT <- rbindlist(full_lst)
  }
  # Save the final output to a file.
  if (cv_method == "LOOCV") {
    cv_scheme <- "loocv"
  }
  out_nm <- paste0("Trait=", init_traits, "_Iter=", init_iter,
                   "_Model=", hypred_model, "_VCOV=", g_method, 
                   "_Pi=", Pi, "_PriorPiCount=", PriorPiCount,
                   "_Imputed=", imputation, "_CV_Scheme=", cv_scheme)
  saveRDS(DT, file = paste0("./data/derived/ssBLUP/", out_nm, ".RDS"))
  
  # When did the job finish?
  print(Sys.time())
}




