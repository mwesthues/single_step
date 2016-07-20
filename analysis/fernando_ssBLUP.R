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
pacman::p_load("methods", "BGLR", "Matrix", "parallel", "data.table", "plyr",
               "reshape2", "matrixStats", "testthat")

# Load functions (including Zhang, RadenI and RadenII as part of 'gmat'.
source("./analysis/snp_functions.R")
# Load functions for the creation of BGLR-kernels.
source("./analysis/fernando_ssBLUP_functions.R")




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
  Sys.setenv("CV_METHOD" = "CV800")
  if (isTRUE(Sys.getenv("CV_METHOD") == "CV800")) {
    Sys.setenv("CV_SCHEME" = "custom")
  }
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
cv_method <- as.character(Sys.getenv("CV_METHOD"))
speed_tst <- as.logical(Sys.getenv("SPEED_TEST"))
if (isTRUE(cv_method == "CV800")) {
  cv_scheme <- as.character(Sys.getenv("CV_SCHEME"))
}
only_profiled <- as.logical(Sys.getenv("ONLY_PROFILED"))
test_that("classes of global parameters are correct", {
  expect_is(use_cores, "integer")
  expect_is(init_traits, "character")
  expect_is(init_iter, "integer")
  expect_is(hypred_model, "character")
  expect_is(g_method, "character")
  expect_is(Pi, "numeric")
  expect_is(PriorPiCount, "integer")
  expect_is(imputation, "logical")
  expect_is(cv_method, "character")
  expect_is(speed_tst, "logical")
  expect_is(only_profiled, "logical")
})


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
if (exists("cv_scheme")) {
  if (cv_scheme == "custom") {
    cv_scheme <- readRDS("./data/input/cust_cv.RDS")
  }
  test_that("CV scheme is specified", {
    expect_gte(length(cv_scheme), expected = 1)
    expect_true(all(cv_scheme %in% list.files("./data/processed/")))
  })
  }  
  cat(cv_scheme, sep = "\n")

  # Load the cross-validation scheme
  cv_lst <- lapply(seq_along(cv_scheme), FUN = function(i) {
    readRDS(paste0("./data/processed/", cv_scheme[i]))
  })
  names(cv_lst) <- cv_scheme

  # Combine all factor levels and store the values in a data frame so that 
  # every possible computation has its own row.
  param_df <- expand.grid(Phenotype = init_traits, 
                          Iter = init_iter, 
                          CV_Scheme = cv_scheme)
} else if (!exists("cv_scheme")) {
  param_df <- expand.grid(Phenotype = init_traits, 
                          Iter = init_iter)
}

for (i in 1:ncol(param_df)) param_df[, i] <- as.character(param_df[, i])
numParam <- nrow(param_df)



## --------------------------------------------------------------------------
## PHENOTYPIC DATA PREPARATION 
pheno_fls <- list.files("./data/processed/", 
                        pattern = "pheno_BLUE")
y_lst <- lapply(seq_along(pheno_fls), FUN = function(i) {
  Pheno <- read.table(paste0("./data/processed/", pheno_fls[i]),
                      header = TRUE)
  trait_pos <- regexpr("(?<=stage2_)[A-Z]+(?=.txt)", text = pheno_fls[i], 
                       perl = TRUE)
  trait_nm <- substring(pheno_fls[i], first = trait_pos,
                        last = trait_pos + attr(trait_pos, "match.length") - 1)
  Pheno$trait <- trait_nm 
  Pheno <- Pheno[Pheno$check == 0 & Pheno$dent.GTP != 0 &
                 Pheno$flint.GTP != 0, ]
  Pheno <- droplevels(Pheno)
  Pheno <- Pheno[grep("DF_", as.character(Pheno$G)), ]
  Pheno$G <- as.character(Pheno$G)
  Pheno$dent.GTP <- as.character(Pheno$dent.GTP)
  Pheno$flint.GTP <- as.character(Pheno$flint.GTP)
  Pheno
})
Pheno <- rbindlist(y_lst)



## --------------------------------------------------------------------------
## PREDICTOR AND AGRONOMIC DATA PREPARATION
# Load predictor data.
snp <- readRDS("./data/processed/snp_mat.RDS")
mrna <- readRDS("./data/processed/subset_mrna_blues.RDS")[["100%"]]
mrna <- t(mrna)
mrna <- mrna[grep("Exp|Sigma", x = rownames(mrna), invert = TRUE), ]
mrna <- mrna[rownames(mrna) %in% rownames(snp), ]
endo_lst <- list(A = snp, M = mrna)

comhybrid <- Pheno[grepl(paste0(rownames(snp), collapse = "|"), 
                         x = Pheno[, G, ]), unique(G), ]
ptlhybrid <- Pheno[dent.GTP %in% rownames(snp) & flint.GTP %in% 
                   rownames(snp), ][, G, ]
comhybrid <- intersect(comhybrid, ptlhybrid)
Pheno <- droplevels(Pheno[Pheno$G %in% comhybrid, , ])
# Store agronomic phenotypes in a matrix.
Y <- dcast.data.table(Pheno, formula = G ~ trait, value.var = "EST")
y_mat <- as.matrix(data.frame(Y[, .(ADF, FETT, GTM, GTS, RPR, STA, XZ), ]))
rownames(y_mat) <- Y[, G, ]
y_mat <- y_mat[complete.cases(y_mat), ]
phenotypes <- colnames(y_mat)
# Extract a vector of parental dent and flint lines, respectively, which have 
# the same length as the number of hybrids. These vectors will be used 
# subsequently for the augmentation of the predictor matrices.
hybrid <- rownames(y_mat)
dent <- sapply(strsplit(rownames(y_mat), split = "_"), FUN = "[[", 2)
flint <- sapply(strsplit(rownames(y_mat), split = "_"), FUN = "[[", 3)
grp_geno <- list(Dent = dent, Flint = flint)

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




